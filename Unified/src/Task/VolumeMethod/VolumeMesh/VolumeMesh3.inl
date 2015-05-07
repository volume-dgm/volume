template<typename FunctionSpace, typename System>
void VolumeMesh<Space3, FunctionSpace, System>::
  LoadGeom(Vector* vertexPositions, IndexType* cellIndices, IndexType verticesCount, IndexType cellsCount,
           FacePairIndices* contactFaces,  IndexType* contactFacesCount,  IndexType contactTypesCount,
           BoundaryFace*    boundaryFaces, IndexType* boundaryFacesCount, IndexType boundaryTypesCount,
           MediumParameters* mediumParameters,
           IndexType *internalContactTypes)
{
  printf("Building geom mesh topology \n");
  GeomMesh<Space>::LoadGeom(vertexPositions, cellIndices, verticesCount, cellsCount);
  GeomMesh<Space>::BuildTopologyInfos();
  GeomMesh<Space>::BuildAdditionalTopology(
    contactFaces, contactFacesCount, contactTypesCount,
    boundaryFaces, boundaryFacesCount, boundaryTypesCount, 
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
void VolumeMesh<Space3, FunctionSpace, System>::BuildMatrices()
{
  #pragma omp parallel for
  for(int functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
  {
    for(int functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
    {
      printf(".");
      fflush (stdout);
      cellVolumeIntegrals[functionIndex0 * functionsCount + functionIndex1] =
        functionSpace.ComputeCellVolumeIntegral       (functionIndex1, functionIndex0);

      Vector3 derivativeIntegral = functionSpace.ComputeDerivativeVolumeIntegral (functionIndex0, functionIndex1);
      xDerivativeVolumeIntegrals [functionIndex0 * functionsCount + functionIndex1] = derivativeIntegral.x;
      yDerivativeVolumeIntegrals [functionIndex0 * functionsCount + functionIndex1] = derivativeIntegral.y;
      zDerivativeVolumeIntegrals [functionIndex0 * functionsCount + functionIndex1] = derivativeIntegral.z;

      for(IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
      {
        outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral[functionIndex0 * functionsCount + functionIndex1] =
          functionSpace.ComputeOutgoingFlux(srcFaceNumber, functionIndex0, functionIndex1);
        for(IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
        {
          for(IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
          {
            incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral[functionIndex0 * functionsCount + functionIndex1] =
              functionSpace.ComputeIncomingFlux(srcFaceNumber, dstFaceNumber, orientationNumber, functionIndex0, functionIndex1);
          }
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

  MatrixCopy(zDerivativeVolumeIntegrals, derivatives, functionsCount, functionsCount);
  MatrixMulMatrix(derivatives, cellVolumeIntegralsInv, zDerivativeVolumeIntegrals, functionsCount, functionsCount, functionsCount);

  for(IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
  {
    MatrixCopy(outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral, flux, functionsCount, functionsCount);
    MatrixMulMatrix(flux, cellVolumeIntegralsInv, outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral, functionsCount, functionsCount, functionsCount);

    for(IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
    {
      for(IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
      {
        MatrixCopy(incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral, flux, functionsCount, functionsCount);
        MatrixMulMatrix(flux, cellVolumeIntegralsInv, incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral, functionsCount, functionsCount, functionsCount);
      }
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
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space3, FunctionSpace, System>::
  GetCurrDerivatives(Scalar *derivatives, const SolverState& solverState)
{
  #pragma omp parallel 
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    Scalar tmpMatrix0[dimsCount * dimsCount];
    Scalar tmpMatrix1[dimsCount * dimsCount];

    Scalar tmpPhaseVelocityMatrix0[dimsCount * functionsCount];
    Scalar tmpPhaseVelocityMatrix1[dimsCount * functionsCount];
    Scalar currCellValues[dimsCount * functionsCount];
    Scalar correspondingCellValues[dimsCount * functionsCount];
    Scalar timeDerivatives[dimsCount * functionsCount];

    Scalar leftPositiveMultMatrix[dimsCount * dimsCount];
    Scalar leftNegativeMultMatrix[dimsCount * dimsCount];

    Scalar faceTransformMatrix[dimsCount * dimsCount];
    Scalar faceTransformMatrixInv[dimsCount * dimsCount];

    Scalar xMixedMatrix[dimsCount * dimsCount];
    Scalar yMixedMatrix[dimsCount * dimsCount];
    Scalar zMixedMatrix[dimsCount * dimsCount];

    Scalar boundaryMatrix[dimsCount * dimsCount];

    Scalar leftContactMatrix[dimsCount * dimsCount];
    Scalar rightContactMatrix[dimsCount * dimsCount];

    Scalar xMatrix[dimsCount * dimsCount];
    Scalar yMatrix[dimsCount * dimsCount];
    Scalar zMatrix[dimsCount * dimsCount];

    Scalar boundaryInfoValues[dimsCount * functionsCount];
    Scalar ghostValues[dimsCount * functionsCount];
    Scalar sourceValues[dimsCount * functionsCount];
    Scalar sourcePointValues[dimsCount * functionsCount];

    Vector cellVertices[Space::NodesPerCell];

    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        continue;
      }
      GetCellVertices(cellIndex, cellVertices);

      Scalar invJacobian = Scalar(1.0) / fabs(GetCellDeformJacobian(cellVertices));

      system.BuildXMatrix(cellMediumParameters[cellIndex], xMatrix);
      system.BuildYMatrix(cellMediumParameters[cellIndex], yMatrix);
      system.BuildZMatrix(cellMediumParameters[cellIndex], zMatrix);

      Vector refXDerivatives = GetRefXDerivatives(cellVertices) * invJacobian;
      Vector refYDerivatives = GetRefYDerivatives(cellVertices) * invJacobian;
      Vector refZDerivatives = GetRefZDerivatives(cellVertices) * invJacobian;

      for (IndexType i = 0; i < dimsCount; i++)
      {
        for (IndexType j = 0; j < dimsCount; j++)
        {
          xMixedMatrix[i * dimsCount + j] = xMatrix[i * dimsCount + j] * refXDerivatives.x +
                                            yMatrix[i * dimsCount + j] * refXDerivatives.y +
                                            zMatrix[i * dimsCount + j] * refXDerivatives.z;

          yMixedMatrix[i * dimsCount + j] = xMatrix[i * dimsCount + j] * refYDerivatives.x +
                                            yMatrix[i * dimsCount + j] * refYDerivatives.y +
                                            zMatrix[i * dimsCount + j] * refYDerivatives.z;

          zMixedMatrix[i * dimsCount + j] = xMatrix[i * dimsCount + j] * refZDerivatives.x +
                                            yMatrix[i * dimsCount + j] * refZDerivatives.y +
                                            zMatrix[i * dimsCount + j] * refZDerivatives.z;
        }
      }

      bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          timeDerivatives[valueIndex * functionsCount + functionIndex] = Scalar(0.0);
          currCellValues[valueIndex * functionsCount + functionIndex] =
            useHalfStepSolution ?
            halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
            cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];
        }
      }

      IndexType cellIncidentNodes[Space::NodesPerCell];
      GetFixedCellIndices(cellIndex, cellIncidentNodes);

      /*
        There are regular cells, where we compute solution,
        and domain boundary cells which are taken from neighbouring domains and should not be computed here.
      */
      bool regularCell = IsCellRegular(cellIndex);

      if (regularCell)
      {
        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
        {
          IndexType faceNodeIndices[Space::NodesPerFace];
          GetCellFaceNodes(cellIncidentNodes, faceNumber, faceNodeIndices);

          Vector faceGlobalVertices[Space::NodesPerFace];
          for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; nodeNumber++)
          {
            faceGlobalVertices[nodeNumber] = nodes[faceNodeIndices[nodeNumber]].pos;
          }

          system.BuildFaceTransformMatrix(faceGlobalVertices, faceTransformMatrix);
          system.BuildFaceTransformMatrixInv(faceGlobalVertices, faceTransformMatrixInv);

          Scalar faceDeformJacobian = Scalar(2.0) * GetFaceSquare(faceGlobalVertices);
          Vector faceNormal = GetFaceExternalNormal(cellIndex, faceNumber);

          IndexType correspondingCellIndex       = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
          IndexType correspondingFaceNumber      = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber;
          IndexType correspondingFaceOrientation = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].orientation;
          IndexType interactionType              = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;

          if (interactionType == IndexType(-1)) continue;

          // interior matrix
          Scalar xInteriorMatrix[dimsCount * dimsCount];
          system.BuildXnInteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xInteriorMatrix);

          // exterior matrix
          Scalar xExteriorMatrix[dimsCount * dimsCount];
          system.BuildXnExteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex], 
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xExteriorMatrix);

          // outgoing flux
          MatrixMulMatrix(faceTransformMatrix, xInteriorMatrix, tmpMatrix0, dimsCount, dimsCount, dimsCount);
          MatrixMulMatrix(tmpMatrix0, faceTransformMatrixInv, leftPositiveMultMatrix, dimsCount, dimsCount, dimsCount);

          MatrixMulMatrix(leftPositiveMultMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
          MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcFaces[faceNumber].surfaceIntegral, tmpPhaseVelocityMatrix1,
            dimsCount, functionsCount, functionsCount);

          MatrixMulScalar(tmpPhaseVelocityMatrix1, faceDeformJacobian * invJacobian, dimsCount, functionsCount);
          MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

          if (correspondingCellIndex == IndexType(-1))
          {
            // bounary condition
            system.BuildBoundaryMatrix(interactionType, boundaryMatrix);
            IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);

            if (!allowDynamicCollisions || dynamicContactType == IndexType(-1))
            {
              // regular boundary
              MatrixMulMatrix(faceTransformMatrix, xExteriorMatrix, tmpMatrix0, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix0, boundaryMatrix, tmpMatrix1, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix1, faceTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

              MatrixMulMatrix(leftNegativeMultMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcFaces[faceNumber].surfaceIntegral,
                tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, faceDeformJacobian * invJacobian, dimsCount, functionsCount);
              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

              // external force/velocity 
              BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);

              typedef BoundaryFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > FunctorWrapper;
              FunctorWrapper wrapper(functor, time, this, cellIndex, faceNumber);

              functionSpace.template Decompose< FunctorWrapper, dimsCount >(wrapper, boundaryInfoValues);

              MatrixMulMatrix(faceTransformMatrix, xExteriorMatrix, tmpMatrix0, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix0, faceTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

              MatrixMulMatrix(leftNegativeMultMatrix, boundaryInfoValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcFaces[faceNumber].surfaceIntegral,
                tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, Scalar(2.0) * faceDeformJacobian * invJacobian, dimsCount, functionsCount);

              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
            } else
            {
              // dynamic collisions
              system.BuildContactMatrices(dynamicContactType, leftContactMatrix, rightContactMatrix);

              BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);

              Scalar boundaryFaceMatrix[dimsCount * dimsCount];
              Scalar leftContactFaceMatrix[dimsCount * dimsCount];
              Scalar rightContactFaceMatrix[dimsCount * dimsCount];

              MatrixMulMatrix(boundaryMatrix, faceTransformMatrixInv, boundaryFaceMatrix, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(leftContactMatrix, faceTransformMatrixInv, leftContactFaceMatrix, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(rightContactMatrix, faceTransformMatrixInv, rightContactFaceMatrix, dimsCount, dimsCount, dimsCount);

              GhostCellFunctionGetter<VolumeMeshT> functionGetter(
                this,
                cellIndex,
                faceGlobalVertices[0], faceNormal,
                faceTransformMatrixInv,
                boundaryFaceMatrix,
                functor,
                leftContactFaceMatrix,
                rightContactFaceMatrix,
                time, dynamicContactType);

              functionSpace.template Decompose<GhostCellFunctionGetter<VolumeMeshT>, dimsCount>(functionGetter, ghostValues);

              MatrixMulMatrix(faceTransformMatrix, xExteriorMatrix, tmpMatrix0, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix0, ghostValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);

              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcFaces[faceNumber].surfaceIntegral,
                tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, faceDeformJacobian * invJacobian, dimsCount, functionsCount);
              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
            }
          } else
          {
            // contact condition
            bool useHalfStepSolutionForCorrespondingCell = timeHierarchyLevelsManager.UseHalfStepSolutionForNeighbour(
              cellIndex, solverState, auxCell, correspondingCellIndex);

            for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
            {
              for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
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
              MatrixMulMatrix(faceTransformMatrix, xExteriorMatrix, tmpMatrix0, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix0, leftContactMatrix, tmpMatrix1, dimsCount, dimsCount, dimsCount);
              MatrixMulMatrix(tmpMatrix1, faceTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

              MatrixMulMatrix(leftNegativeMultMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
              MatrixMulMatrix(tmpPhaseVelocityMatrix0, outgoingFlux.srcFaces[faceNumber].surfaceIntegral,
                tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);

              MatrixMulScalar(tmpPhaseVelocityMatrix1, faceDeformJacobian * invJacobian, dimsCount, functionsCount);
              MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
            }

            // exterior side contribution
            MatrixMulMatrix(faceTransformMatrix, xExteriorMatrix, tmpMatrix0, dimsCount, dimsCount, dimsCount);
            MatrixMulMatrix(tmpMatrix0, rightContactMatrix, tmpMatrix1, dimsCount, dimsCount, dimsCount);
            MatrixMulMatrix(tmpMatrix1, faceTransformMatrixInv, leftNegativeMultMatrix, dimsCount, dimsCount, dimsCount);

            MatrixMulMatrix(leftNegativeMultMatrix, correspondingCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
            MatrixMulMatrix(tmpPhaseVelocityMatrix0,
              incomingFlux.srcFaces[faceNumber].dstFaces[correspondingFaceNumber].orientations[correspondingFaceOrientation].surfaceIntegral,
              tmpPhaseVelocityMatrix1,
              dimsCount, functionsCount, functionsCount);

            MatrixMulScalar(tmpPhaseVelocityMatrix1, faceDeformJacobian * invJacobian, dimsCount, functionsCount);
            MatrixAddMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);
          }
        }
      }
      MatrixMulMatrix(xMixedMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
      MatrixMulMatrix(tmpPhaseVelocityMatrix0, xDerivativeVolumeIntegrals, tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);
      MatrixSubMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

      MatrixMulMatrix(yMixedMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
      MatrixMulMatrix(tmpPhaseVelocityMatrix0, yDerivativeVolumeIntegrals, tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);
      MatrixSubMatrix(timeDerivatives, tmpPhaseVelocityMatrix1, dimsCount, functionsCount);

      MatrixMulMatrix(zMixedMatrix, currCellValues, tmpPhaseVelocityMatrix0, dimsCount, dimsCount, functionsCount);
      MatrixMulMatrix(tmpPhaseVelocityMatrix0, zDerivativeVolumeIntegrals, tmpPhaseVelocityMatrix1, dimsCount, functionsCount, functionsCount);
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

      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          derivatives[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex]
            = -timeDerivatives[valueIndex * functionsCount + functionIndex];
        }
      }
      targetCellIndex++;
    }
  }
}

template< typename FunctionSpace, typename System>
typename Space3::Scalar VolumeMesh<Space3, FunctionSpace, System>::
  GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const
{
  // deform jacobian equals 6 * tetrahedron volume
  return 
    cellVertices[0].x * ( cellVertices[1].y * (cellVertices[3].z - cellVertices[2].z) +
                          cellVertices[2].y * (cellVertices[1].z - cellVertices[3].z) +
                          cellVertices[3].y * (cellVertices[2].z - cellVertices[1].z)) +

    cellVertices[1].x * ( cellVertices[0].y * (cellVertices[2].z - cellVertices[3].z) +
                          cellVertices[2].y * (cellVertices[3].z - cellVertices[0].z) +
                          cellVertices[3].y * (cellVertices[0].z - cellVertices[2].z)) +

    cellVertices[2].x * ( cellVertices[0].y * (cellVertices[3].z - cellVertices[1].z) +
                          cellVertices[1].y * (cellVertices[0].z - cellVertices[3].z) +
                          cellVertices[3].y * (cellVertices[1].z - cellVertices[0].z)) +

    cellVertices[3].x * ( cellVertices[0].y * (cellVertices[1].z - cellVertices[2].z) +
                          cellVertices[1].y * (cellVertices[2].z - cellVertices[0].z) +
                          cellVertices[2].y * (cellVertices[0].z - cellVertices[1].z));
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const //x -> ?
{
  Scalar invJacobian = Scalar(1.0) / GetCellDeformJacobian(cellVertices);
  Vector refXDerivatives = GetRefXDerivatives(cellVertices);
  Vector refYDerivatives = GetRefYDerivatives(cellVertices);
  Vector refZDerivatives = GetRefZDerivatives(cellVertices);
  
  return invJacobian * Vector(
    cellVertices[0].x * (cellVertices[3].y * cellVertices[2].z - cellVertices[2].y * cellVertices[3].z) +
    cellVertices[2].x * (cellVertices[0].y * cellVertices[3].z - cellVertices[3].y * cellVertices[0].z) +
    cellVertices[3].x * (cellVertices[2].y * cellVertices[0].z - cellVertices[0].y * cellVertices[2].z) +

    globalCoords.x * refXDerivatives.x + globalCoords.y * refXDerivatives.y + globalCoords.z * refXDerivatives.z,


    cellVertices[0].y * (cellVertices[3].x * cellVertices[1].z - cellVertices[1].x * cellVertices[3].z) +
    cellVertices[1].y * (cellVertices[0].x * cellVertices[3].z - cellVertices[3].x * cellVertices[0].z) +
    cellVertices[3].y * (cellVertices[1].x * cellVertices[0].z - cellVertices[0].x * cellVertices[1].z) +

    globalCoords.x * refYDerivatives.x + globalCoords.y * refYDerivatives.y + globalCoords.z * refYDerivatives.z,
    

    cellVertices[0].z * (cellVertices[2].x * cellVertices[1].y - cellVertices[1].x * cellVertices[2].y) +
    cellVertices[1].z * (cellVertices[0].x * cellVertices[2].y - cellVertices[2].x * cellVertices[0].y) +
    cellVertices[2].z * (cellVertices[1].x * cellVertices[0].y - cellVertices[0].x * cellVertices[1].y) +

    globalCoords.x * refZDerivatives.x + globalCoords.y * refZDerivatives.y + globalCoords.z * refZDerivatives.z
    );
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const // ? -> x
{
  return Vector(
    cellVertices[0].x + (cellVertices[1].x - cellVertices[0].x) * refCoords.x + 
                        (cellVertices[2].x - cellVertices[0].x) * refCoords.y + 
                        (cellVertices[3].x - cellVertices[0].x) * refCoords.z,
    cellVertices[0].y + (cellVertices[1].y - cellVertices[0].y) * refCoords.x + 
                        (cellVertices[2].y - cellVertices[0].y) * refCoords.y + 
                        (cellVertices[3].y - cellVertices[0].y) * refCoords.z,
    cellVertices[0].z + (cellVertices[1].z - cellVertices[0].z) * refCoords.x + 
                        (cellVertices[2].z - cellVertices[0].z) * refCoords.y + 
                        (cellVertices[3].z - cellVertices[0].z) * refCoords.z
    );
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GetRefXDerivatives(Vector cellVertices[Space::NodesPerCell]) const //(dξ/dx, dξ/dy, dξ/dz) * J
{
    return Vector(
      cellVertices[0].y * (cellVertices[2].z - cellVertices[3].z) + 
      cellVertices[2].y * (cellVertices[3].z - cellVertices[0].z) + 
      cellVertices[3].y * (cellVertices[0].z - cellVertices[2].z),
      cellVertices[0].x * (cellVertices[3].z - cellVertices[2].z) + 
      cellVertices[2].x * (cellVertices[0].z - cellVertices[3].z) + 
      cellVertices[3].x * (cellVertices[2].z - cellVertices[0].z),
      cellVertices[0].x * (cellVertices[2].y - cellVertices[3].y) + 
      cellVertices[2].x * (cellVertices[3].y - cellVertices[0].y) + 
      cellVertices[3].x * (cellVertices[0].y - cellVertices[2].y));
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GetRefYDerivatives(Vector cellVertices[Space::NodesPerCell]) const //(dη/dx, dη/dy, dη/dz) * J
{
    return Vector(
      cellVertices[0].y * (cellVertices[3].z - cellVertices[1].z) + 
      cellVertices[1].y * (cellVertices[0].z - cellVertices[3].z) + 
      cellVertices[3].y * (cellVertices[1].z - cellVertices[0].z),
      cellVertices[0].x * (cellVertices[1].z - cellVertices[3].z) + 
      cellVertices[1].x * (cellVertices[3].z - cellVertices[0].z) + 
      cellVertices[3].x * (cellVertices[0].z - cellVertices[1].z),
      cellVertices[0].x * (cellVertices[3].y - cellVertices[1].y) + 
      cellVertices[1].x * (cellVertices[0].y - cellVertices[3].y) + 
      cellVertices[3].x * (cellVertices[1].y - cellVertices[0].y));
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GetRefZDerivatives(Vector cellVertices[Space::NodesPerCell]) const //(dζ/dx, dζ/dy, dζ/dz) * J
{
    return Vector(
      cellVertices[0].y * (cellVertices[1].z - cellVertices[2].z) + 
      cellVertices[1].y * (cellVertices[2].z - cellVertices[0].z) + 
      cellVertices[2].y * (cellVertices[0].z - cellVertices[1].z),
      cellVertices[0].x * (cellVertices[2].z - cellVertices[1].z) + 
      cellVertices[1].x * (cellVertices[0].z - cellVertices[2].z) + 
      cellVertices[2].x * (cellVertices[1].z - cellVertices[0].z),
      cellVertices[0].x * (cellVertices[1].y - cellVertices[2].y) + 
      cellVertices[1].x * (cellVertices[2].y - cellVertices[0].y) + 
      cellVertices[2].x * (cellVertices[0].y - cellVertices[1].y));
}


template<typename FunctionSpace, typename System>
bool VolumeMesh<Space3, FunctionSpace, System>::IsCellRegular(IndexType cellIndex) const
{
  bool regularCell = true;
  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
  {
    IndexType interactionType = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;
    if (interactionType == IndexType(-1)) regularCell = false;
  }
  return regularCell;
}
