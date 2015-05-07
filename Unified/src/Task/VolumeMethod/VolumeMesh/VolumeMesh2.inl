template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::
  LoadGeom(Vector* vertexPositions, IndexType* cellIndices, IndexType verticesCount, IndexType cellsCount,
           EdgePairIndices* contactEdges, IndexType* contactEdgesCount, IndexType contactTypesCount,
           BoundaryEdge*    boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount,
           MediumParameters* mediumParameters,
           IndexType *internalContactTypes)
{
  printf("Building geom mesh topology \n");

  GeomMesh<Space>::LoadGeom(vertexPositions, cellIndices, verticesCount, cellsCount);
  GeomMesh<Space>::BuildTopologyInfos();
  GeomMesh<Space>::BuildAdditionalTopology(
    contactEdges,  contactEdgesCount,  contactTypesCount,
    boundaryEdges, boundaryEdgesCount, boundaryTypesCount,
    internalContactTypes);

  Initialize();
  printf("Building volume method additional matrices \n");
  BuildMatrices();

  // medium parameters setting
  if (mediumParameters) 
  {
    cellMediumParameters.resize(cellsCount);
    std::copy(mediumParameters, mediumParameters + cellsCount, cellMediumParameters.begin());
  }

  printf("Loading done \n");
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::BuildMatrices()
{
  #pragma omp parallel for
  for(int functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
  {
    for(int functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
    {
      printf(".");
      cellVolumeIntegrals[functionIndex0 * functionsCount + functionIndex1] =
        functionSpace.ComputeCellVolumeIntegral       (functionIndex1, functionIndex0);

      Vector derivativeIntegral = functionSpace.ComputeDerivativeVolumeIntegral (functionIndex0, functionIndex1);
      xDerivativeVolumeIntegrals [functionIndex0 * functionsCount + functionIndex1] = derivativeIntegral.x;
      yDerivativeVolumeIntegrals [functionIndex0 * functionsCount + functionIndex1] = derivativeIntegral.y;

      for(IndexType srcEdgeNumber = 0; srcEdgeNumber < 3; srcEdgeNumber++)
      {
        outgoingFlux.srcEdges[srcEdgeNumber].surfaceIntegral[functionIndex0 * functionsCount + functionIndex1] =
          functionSpace.ComputeOutgoingFlux(srcEdgeNumber, functionIndex0, functionIndex1);
        for(IndexType dstEdgeNumber = 0; dstEdgeNumber < 3; dstEdgeNumber++)
        {
          incomingFlux.srcEdges[srcEdgeNumber].dstEdges[dstEdgeNumber].surfaceIntegral[functionIndex0 * functionsCount + functionIndex1] =
            functionSpace.ComputeIncomingFlux(srcEdgeNumber, dstEdgeNumber, functionIndex0, functionIndex1);
        }
      }
    }
  }
  MatrixInverse(cellVolumeIntegrals, cellVolumeIntegralsInv, functionsCount);

  Scalar flux[functionsCount * functionsCount];
  Scalar derivatives[functionsCount * functionsCount];

  MatrixCopy(xDerivativeVolumeIntegrals, derivatives, functionsCount, functionsCount);
  MatrixMulMatrix(derivatives, cellVolumeIntegralsInv, xDerivativeVolumeIntegrals, functionsCount, functionsCount, functionsCount);

  MatrixCopy(yDerivativeVolumeIntegrals, derivatives, functionsCount, functionsCount);
  MatrixMulMatrix(derivatives, cellVolumeIntegralsInv, yDerivativeVolumeIntegrals, functionsCount, functionsCount, functionsCount);

  for(IndexType srcEdgeNumber = 0; srcEdgeNumber < 3; srcEdgeNumber++)
  {
    MatrixCopy(outgoingFlux.srcEdges[srcEdgeNumber].surfaceIntegral, flux, functionsCount, functionsCount);
    MatrixMulMatrix(flux, cellVolumeIntegralsInv, outgoingFlux.srcEdges[srcEdgeNumber].surfaceIntegral, functionsCount, functionsCount, functionsCount);

    for(IndexType dstEdgeNumber = 0; dstEdgeNumber < 3; dstEdgeNumber++)
    {
      MatrixCopy(incomingFlux.srcEdges[srcEdgeNumber].dstEdges[dstEdgeNumber].surfaceIntegral, flux, functionsCount, functionsCount);
      MatrixMulMatrix(flux, cellVolumeIntegralsInv, incomingFlux.srcEdges[srcEdgeNumber].dstEdges[dstEdgeNumber].surfaceIntegral, functionsCount, functionsCount, functionsCount);
    }
  }

  Scalar testMatrix[functionsCount * functionsCount];
  Scalar err = 0;

  MatrixMulMatrix(cellVolumeIntegrals, cellVolumeIntegralsInv, testMatrix, functionsCount, functionsCount, functionsCount);
  for(IndexType i = 0; i < functionsCount; i++)
  {
    for(IndexType j = 0; j < functionsCount; j++)
    {
      if(i == j)
        err += fabs(testMatrix[i + j * functionsCount] - Scalar(1.0));
      else
        err += fabs(testMatrix[i + j * functionsCount]);
    }
  }
  printf("\nPrecomputations complete, volume integral matrix error is : %f\n", err);

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    cellVolumeAverageIntegrals[functionIndex] = functionSpace.ComputeCellVolumeIntegral(functionIndex);

    for(IndexType edgeNumber = 0; edgeNumber < 3; edgeNumber++)
    {
      edgeAverages[edgeNumber].surfaceIntegral[functionIndex] =
        functionSpace.ComputeEdgeFlux(edgeNumber, functionIndex);
    }
  }
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::
  GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState)
{
  #pragma omp parallel 
  {
    int threadIndex  = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    Scalar tmpMatrix0[dimsCount * dimsCount];
    Scalar tmpMatrix1[dimsCount * dimsCount];

    Scalar tmpPhaseVelocityMatrix0[dimsCount * functionsCount];
    Scalar tmpPhaseVelocityMatrix1[dimsCount * functionsCount];
    Scalar currCellValues         [dimsCount * functionsCount];
    Scalar correspondingCellValues[dimsCount * functionsCount];
    Scalar timeDerivatives        [dimsCount * functionsCount];

    Scalar leftPositiveMultMatrix[dimsCount * dimsCount];
    Scalar leftNegativeMultMatrix[dimsCount * dimsCount];

    Scalar edgeTransformMatrix    [dimsCount * dimsCount];
    Scalar edgeTransformMatrixInv [dimsCount * dimsCount];

    Scalar xMixedMatrix           [dimsCount * dimsCount];
    Scalar yMixedMatrix           [dimsCount * dimsCount];

    Scalar boundaryMatrix         [dimsCount * dimsCount];

    Scalar leftContactMatrix      [dimsCount * dimsCount];
    Scalar rightContactMatrix     [dimsCount * dimsCount];

    Scalar xMatrix                [dimsCount * dimsCount];
    Scalar yMatrix                [dimsCount * dimsCount];

    Scalar boundaryInfoValues     [dimsCount * functionsCount];
    Scalar ghostValues            [dimsCount * functionsCount];
    Scalar sourceValues           [dimsCount * functionsCount];
    Scalar sourcePointValues      [dimsCount * functionsCount];

    Vector cellVertices[Space::NodesPerCell];

    for(int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        continue;
      }

      bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
      for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          timeDerivatives[valueIndex * functionsCount + functionIndex] = Scalar(0.0);
          currCellValues[valueIndex * functionsCount + functionIndex] = 
            useHalfStepSolution ?
            halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
            cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];
        }
      }

      GetCellVertices(cellIndex, cellVertices);
      Scalar invJacobian = Scalar(1.0) / fabs(GetCellDeformJacobian(cellVertices));

      IndexType cellIncidentNodes[Space::NodesPerCell];
      GetFixedCellIndices(cellIndex, cellIncidentNodes);

      /*
        There are regular cells, where we compute solution,
        and domain boundary cells which are taken from neighbouring domains and should not be computed here.
      */
      bool regularCell = IsCellRegular(cellIndex);

      if (regularCell)
      {
        for(IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
        {
          IndexType edgeNodeIndices[Space::NodesPerEdge];
          GetCellEdgeNodes(cellIncidentNodes, edgeNumber, edgeNodeIndices);

          Vector edgeGlobalVertices[Space::NodesPerEdge];
          for(IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; nodeNumber++)
          {
            edgeGlobalVertices[nodeNumber] = nodes[edgeNodeIndices[nodeNumber]].pos;
          }

          system.BuildEdgeTransformMatrix   (edgeGlobalVertices, edgeTransformMatrix   );
          system.BuildEdgeTransformMatrixInv(edgeGlobalVertices, edgeTransformMatrixInv);

          Scalar edgeLen = (edgeGlobalVertices[1] - edgeGlobalVertices[0]).Len(); //GetFaceSquare(faceGlobalVertices);
          Vector edgeNormal = GetEdgeExternalNormal(cellIndex, edgeNumber);

          IndexType correspondingCellIndex  = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex;
          IndexType correspondingEdgeNumber = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber;
          IndexType interactionType         = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;

          if (interactionType == IndexType(-1))
          {
            continue;
          }

          // interior matrix
          Scalar xInteriorMatrix[dimsCount * dimsCount];
          system.BuildXnInteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            edgeNormal, xInteriorMatrix);

          // exterior matrix
          Scalar xExteriorMatrix[dimsCount * dimsCount];
          system.BuildXnExteriorMatrix(
            cellMediumParameters[cellIndex], 
            //cellMediumParameters[cellIndex], 
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            edgeNormal, xExteriorMatrix);

          // outgoing flux
          MatrixMulMatrix(edgeTransformMatrix, xInteriorMatrix,    tmpMatrix0,     dimsCount, dimsCount, dimsCount);
          MatrixMulMatrix(tmpMatrix0, edgeTransformMatrixInv, leftPositiveMultMatrix, dimsCount, dimsCount, dimsCount);

          MatrixMulMatrix(leftPositiveMultMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
          MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcEdges[edgeNumber].surfaceIntegral, tmpPhaseVelocityMatrix1, 
            dimsCount, functionsCount, functionsCount);

          MatrixMulScalar(tmpPhaseVelocityMatrix1, edgeLen * invJacobian, dimsCount, functionsCount);
          MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

          if(correspondingCellIndex == IndexType(-1))
          {
            // bounary condition
            system.BuildBoundaryMatrix(interactionType, boundaryMatrix);
            IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);

            if(!allowDynamicCollisions || dynamicContactType == IndexType(-1)) //set 1 for regular boundary, 0 for dynamic collisions
            {
              MatrixMulMatrix(edgeTransformMatrix, xExteriorMatrix,    tmpMatrix0,     dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix0, boundaryMatrix, tmpMatrix1, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix1, edgeTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

              MatrixMulMatrix(leftNegativeMultMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcEdges[edgeNumber].surfaceIntegral, 
                              tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, edgeLen * invJacobian, dimsCount, functionsCount);
              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

              // external force/velocity 
              BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);

              typedef BoundaryFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > FunctorWrapper;
              FunctorWrapper wrapper(functor, time, this, cellIndex, edgeNumber);

              functionSpace.template Decompose< FunctorWrapper, dimsCount >(wrapper, boundaryInfoValues);

              MatrixMulMatrix(edgeTransformMatrix, xExteriorMatrix,    tmpMatrix1,     dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix1, edgeTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

              MatrixMulMatrix(leftNegativeMultMatrix, boundaryInfoValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcEdges[edgeNumber].surfaceIntegral, 
                              tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, Scalar(2.0) * edgeLen * invJacobian, dimsCount, functionsCount);

              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
            } else
            {
              system.BuildContactMatrices(dynamicContactType, leftContactMatrix, rightContactMatrix);

              BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);

              Scalar boundaryEdgeMatrix[dimsCount * dimsCount];
              Scalar leftContactEdgeMatrix[dimsCount * dimsCount];
              Scalar rightContactEdgeMatrix[dimsCount * dimsCount];

              MatrixMulMatrix(boundaryMatrix,     edgeTransformMatrixInv, boundaryEdgeMatrix,     dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(leftContactMatrix,  edgeTransformMatrixInv, leftContactEdgeMatrix,  dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(rightContactMatrix, edgeTransformMatrixInv, rightContactEdgeMatrix, dimsCount, dimsCount, dimsCount);

              GhostCellFunctionGetter<VolumeMeshT> functionGetter(
                this, 
                cellIndex,
                edgeGlobalVertices[0], edgeNormal,
                edgeTransformMatrixInv,
                boundaryEdgeMatrix,
                functor,
                leftContactEdgeMatrix,
                rightContactEdgeMatrix,
                time, dynamicContactType);

              functionSpace.template Decompose<GhostCellFunctionGetter<VolumeMeshT>, dimsCount>(functionGetter, ghostValues);

              MatrixMulMatrix(edgeTransformMatrix, xExteriorMatrix,    tmpMatrix0,     dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix0, ghostValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);

              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcEdges[edgeNumber].surfaceIntegral, 
                              tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, edgeLen * invJacobian, dimsCount, functionsCount);
              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
            }
          } else
          {
            // contact condition
            bool useHalfStepSolutionForCorrespondingCell = timeHierarchyLevelsManager.UseHalfStepSolutionForNeighbour(
              cellIndex, solverState, auxCell, correspondingCellIndex);
            
            for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
            {
              for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
              {
                correspondingCellValues[valueIndex * functionsCount + functionIndex] = 
                  useHalfStepSolutionForCorrespondingCell ?
                  halfStepCellSolutions[correspondingCellIndex].basisVectors[functionIndex].values[valueIndex] :
                  cellSolutions[correspondingCellIndex].basisVectors[functionIndex].values[valueIndex];
              }
            }

            system.BuildContactMatrices(interactionType, leftContactMatrix, rightContactMatrix);

            // for glue contact left matrix equals 0
            bool zeroLeftContactMatrix = true;
            for (IndexType i = 0; i < dimsCount * dimsCount; ++i)
            {
              if (fabs(leftContactMatrix[i]) > std::numeric_limits<Scalar>::epsilon()) zeroLeftContactMatrix = false;
            }

            if (!zeroLeftContactMatrix)
            {
              // interior side contribution
              MatrixMulMatrix(edgeTransformMatrix, xExteriorMatrix,    tmpMatrix0,     dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix0, leftContactMatrix, tmpMatrix1, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix1, edgeTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

              MatrixMulMatrix(leftNegativeMultMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcEdges[edgeNumber].surfaceIntegral, 
                              tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, edgeLen * invJacobian, dimsCount, functionsCount);
              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
            }

            // exterior side contribution
            MatrixMulMatrix(edgeTransformMatrix, xExteriorMatrix,    tmpMatrix0,     dimsCount, dimsCount, dimsCount);
            MatrixMulMatrix(tmpMatrix0, rightContactMatrix, tmpMatrix1, dimsCount, dimsCount, dimsCount);
            MatrixMulMatrix(tmpMatrix1, edgeTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

            MatrixMulMatrix(leftNegativeMultMatrix, correspondingCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
            MatrixMulMatrix(tmpPhaseVelocityMatrix0, 
              incomingFlux.srcEdges[edgeNumber].dstEdges[correspondingEdgeNumber].surfaceIntegral,
              tmpPhaseVelocityMatrix1,
              dimsCount, functionsCount, functionsCount);

            MatrixMulScalar(tmpPhaseVelocityMatrix1, edgeLen * invJacobian, dimsCount, functionsCount);
            MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
          }
        }
      }

      system.BuildXMatrix(cellMediumParameters[cellIndex], xMatrix);
      system.BuildYMatrix(cellMediumParameters[cellIndex], yMatrix);

      Vector refXDerivatives = (GetRefXDerivatives(cellVertices)) * invJacobian;
      Vector refYDerivatives = (GetRefYDerivatives(cellVertices)) * invJacobian;

      for (IndexType i = 0; i < dimsCount; i++)
      {
        for (IndexType j = 0; j < dimsCount; j++)
        {
          xMixedMatrix[i * dimsCount + j] = xMatrix[i * dimsCount + j] * refXDerivatives.x +
            yMatrix[i * dimsCount + j] * refXDerivatives.y;

          yMixedMatrix[i * dimsCount + j] = xMatrix[i * dimsCount + j] * refYDerivatives.x +
            yMatrix[i * dimsCount + j] * refYDerivatives.y;
        }
      }

      MatrixMulMatrix(xMixedMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
      MatrixMulMatrix(tmpPhaseVelocityMatrix0, xDerivativeVolumeIntegrals, tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);
      MatrixSubMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

      MatrixMulMatrix(yMixedMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
      MatrixMulMatrix(tmpPhaseVelocityMatrix0, yDerivativeVolumeIntegrals, tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);
      MatrixSubMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

      typename SystemT::SourceFunctorT* sourceFunctor = system.GetSourceFunctor();
      if (sourceFunctor)
      {
        typedef SourceFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > SourceFunctorWrapper;
        SourceFunctorWrapper wrapper(sourceFunctor, time, this, cellIndex);
        functionSpace.template Decompose< SourceFunctorWrapper, dimsCount >(wrapper, sourceValues);
        MatrixSubMatrix(timeDerivatives, sourceValues, dimsCount, functionsCount);
      }

      // point sources
      for (IndexType sourceIndex = 0; sourceIndex < system.pointSources.size(); ++sourceIndex)
      {
        Vector point = system.pointSources[sourceIndex]->GetPoint();
        if (PointInCell<Scalar>(cellVertices, point))
        {
          Vector refPoint = GlobalToRefVolumeCoords(point, cellVertices);
          Scalar values[dimsCount];
          (*system.pointSources[sourceIndex])(time, values);
          for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
          {
            for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
            {
              sourcePointValues[valueIndex * functionsCount + functionIndex] =
                values[valueIndex] * functionSpace.GetBasisFunctionValue(refPoint, functionIndex);
            }
          }
          MatrixMulMatrix(sourcePointValues, cellVolumeIntegralsInv, tmpPhaseVelocityMatrix0, dimsCount, functionsCount, functionsCount);
          MatrixSubMatrix(timeDerivatives, tmpPhaseVelocityMatrix0, dimsCount, functionsCount);
        }
      }

      for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          derivatives[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex]
            = -timeDerivatives[valueIndex * functionsCount + functionIndex];
        }
      }
      targetCellIndex++;
    }
  }
}

template<typename FunctionSpace, typename System>
inline typename Space2::Scalar VolumeMesh<Space2, FunctionSpace, System>::
  GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const
{
  return (cellVertices[1].x - cellVertices[0].x) * (cellVertices[2].y - cellVertices[0].y) -
         (cellVertices[2].x - cellVertices[0].x) * (cellVertices[1].y - cellVertices[0].y);
}

template<typename FunctionSpace, typename System>
inline typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GetRefXDerivatives(Vector cellVertices[Space::NodesPerCell]) const //(dξ/dx, dξ/dy, dξ/dz) * J
{
    return Vector(
      cellVertices[2].y - cellVertices[0].y,
      cellVertices[0].x - cellVertices[2].x);
}

template<typename FunctionSpace, typename System>
inline typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GetRefYDerivatives(Vector cellVertices[Space::NodesPerCell]) const //(dη/dx, dη/dy, dη/dz) * J
{
    return Vector(
      cellVertices[0].y - cellVertices[1].y,
      cellVertices[1].x - cellVertices[0].x);
}

template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const //x -> ?
{
  Scalar invJacobian = Scalar(1.0) / GetCellDeformJacobian(cellVertices);
  Vector refXDerivatives = GetRefXDerivatives(cellVertices);
  Vector refYDerivatives = GetRefYDerivatives(cellVertices);

  return invJacobian * Vector(
    cellVertices[2].x * cellVertices[0].y - cellVertices[0].x * cellVertices[2].y +
    globalCoords.x * refXDerivatives.x + globalCoords.y * refXDerivatives.y,

    cellVertices[0].x * cellVertices[1].y - cellVertices[1].x * cellVertices[0].y +
    globalCoords.x * refYDerivatives.x + globalCoords.y * refYDerivatives.y);
}

template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const // ? -> x
{
  return Vector(
    cellVertices[0].x + (cellVertices[1].x - cellVertices[0].x) * refCoords.x +
                        (cellVertices[2].x - cellVertices[0].x) * refCoords.y,
    cellVertices[0].y + (cellVertices[1].y - cellVertices[0].y) * refCoords.x +
                        (cellVertices[2].y - cellVertices[0].y) * refCoords.y
    );
}

template<typename FunctionSpace, typename System>
typename System::ValueType VolumeMesh<Space2, FunctionSpace, System>::
  GetEdgeAverageSolution(IndexType cellIndex, IndexType edgeNumber) const
{
  typename System::ValueType result;
  result.SetZeroValues();

  for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
  {
    for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++) 
    {
      Scalar basisFunctionCoefficient = 
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];

      result.values[valueIndex] += basisFunctionCoefficient * edgeAverages[edgeNumber].surfaceIntegral[functionIndex]; 
    }
  }
  return result;
}

template<typename FunctionSpace, typename System>
bool VolumeMesh<Space2, FunctionSpace, System>::IsCellRegular(IndexType cellIndex) const
{
  bool regularCell = true;
  for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
  {
    IndexType interactionType = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;
    if (interactionType == IndexType(-1)) regularCell = false;
  }
  return regularCell;
}
