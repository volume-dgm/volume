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
      cellVolumeIntegrals(functionIndex0, functionIndex1) =
        functionSpace->ComputeCellVolumeIntegral(functionIndex1, functionIndex0);

      Vector derivativeIntegral = functionSpace->ComputeDerivativeVolumeIntegral(functionIndex0, functionIndex1);
      xDerivativeVolumeIntegrals(functionIndex0, functionIndex1) = derivativeIntegral.x;
      yDerivativeVolumeIntegrals(functionIndex0, functionIndex1) = derivativeIntegral.y;
      zDerivativeVolumeIntegrals(functionIndex0, functionIndex1) = derivativeIntegral.z;

      for(IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
      {
        outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral(functionIndex0, functionIndex1) =
          functionSpace->ComputeOutgoingFlux(srcFaceNumber, functionIndex0, functionIndex1);
        for(IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
        {
          for(IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
          {
            incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral(functionIndex0, functionIndex1) =
              functionSpace->ComputeIncomingFlux(srcFaceNumber, dstFaceNumber, orientationNumber, functionIndex0, functionIndex1);
          }
        }
      }
    }
  }
  cellVolumeIntegralsInv = cellVolumeIntegrals.inverse();

  xDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;
  yDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;
  zDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;

  for(IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
  {
    outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral *= cellVolumeIntegralsInv;
    for(IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
    {
      for(IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
      {
        incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral *= cellVolumeIntegralsInv;
      }
    }
  }

  Eigen::Matrix<Scalar, functionsCount, functionsCount> testMatrix;
  Scalar err = 0;

  testMatrix = cellVolumeIntegrals * cellVolumeIntegralsInv;
  for(IndexType i = 0; i < functionsCount; i++)
  {
    for(IndexType j = 0; j < functionsCount; j++)
    {
      if(i == j)
        err += fabs(testMatrix(i, j) - Scalar(1.0));
      else
        err += fabs(testMatrix(i, j));
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

    Eigen::Matrix<Scalar, dimsCount, functionsCount> currCellValues;
    Eigen::Matrix<Scalar, dimsCount, functionsCount> correspondingCellValues;
    Eigen::Matrix<Scalar, dimsCount, functionsCount> timeDerivatives;

    Eigen::Matrix<Scalar, dimsCount, dimsCount> faceTransformMatrix;
    Eigen::Matrix<Scalar, dimsCount, dimsCount> faceTransformMatrixInv;

    Eigen::Matrix<Scalar, dimsCount, dimsCount> xInteriorMatrix;
    Eigen::Matrix<Scalar, dimsCount, dimsCount> xExteriorMatrix;

    Eigen::Matrix<Scalar, 1, dimsCount> boundaryMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> leftContactMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> rightContactMatrix;

    Eigen::Matrix<Scalar, dimsCount, dimsCount> xMatrix;
    Eigen::Matrix<Scalar, dimsCount, dimsCount> yMatrix;
    Eigen::Matrix<Scalar, dimsCount, dimsCount> zMatrix;

    Eigen::Matrix<Scalar, dimsCount, dimsCount> xMixedMatrix;
    Eigen::Matrix<Scalar, dimsCount, dimsCount> yMixedMatrix;
    Eigen::Matrix<Scalar, dimsCount, dimsCount> zMixedMatrix;

    Eigen::Matrix<Scalar, dimsCount, functionsCount> boundaryInfoValues;
    Eigen::Matrix<Scalar, dimsCount, functionsCount> ghostValues;
    Eigen::Matrix<Scalar, dimsCount, functionsCount> sourceValues;
    Eigen::Matrix<Scalar, dimsCount, functionsCount> sourcePointValues;

    Vector cellVertices[Space::NodesPerCell];

    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        continue;
      }

      timeDerivatives.setZero();

      bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          currCellValues(valueIndex, functionIndex) =
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
          system.BuildXnInteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xInteriorMatrix);

          // exterior matrix
          system.BuildXnExteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex], 
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xExteriorMatrix);

          // outgoing flux
          timeDerivatives.noalias() += (faceDeformJacobian * invJacobian) * (faceTransformMatrix * xInteriorMatrix * faceTransformMatrixInv *
            currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral);

          if (correspondingCellIndex == IndexType(-1))
          {
            // bounary condition
            system.BuildBoundaryMatrix(interactionType, boundaryMatrix);
            IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);

            if (!allowDynamicCollisions || dynamicContactType == IndexType(-1))
            {
              // regular boundary          
              timeDerivatives.noalias() += (faceDeformJacobian * invJacobian) * (faceTransformMatrix * xExteriorMatrix * 
                boundaryMatrix.asDiagonal() * faceTransformMatrixInv * currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral);

              // external force/velocity 
              BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);
              if (functor)
              {
                typedef BoundaryFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > FunctorWrapper;
                FunctorWrapper wrapper(functor, time, this, cellIndex, faceNumber);

                functionSpace->template Decompose< FunctorWrapper, dimsCount >(wrapper, boundaryInfoValues.data());

                timeDerivatives.noalias() += (Scalar(2.0) * faceDeformJacobian * invJacobian) * 
                  (faceTransformMatrix * xExteriorMatrix * faceTransformMatrixInv * boundaryInfoValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral);
              }
            } else
            {
              // dynamic collisions
              /*
              system.BuildContactMatrices(dynamicContactType, leftContactMatrix, rightContactMatrix);

              BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);

              Eigen::Matrix<Scalar, dimsCount, dimsCount> boundaryFaceMatrix     = boundaryMatrix.asDiagonal()     * faceTransformMatrixInv;
              Eigen::Matrix<Scalar, dimsCount, dimsCount> leftContactFaceMatrix  = leftContactMatrix.asDiagonal()  * faceTransformMatrixInv;
              Eigen::Matrix<Scalar, dimsCount, dimsCount> rightContactFaceMatrix = rightContactMatrix.asDiagonal() * faceTransformMatrixInv;

              GhostCellFunctionGetter<VolumeMeshT> functionGetter(
                this,
                cellIndex,
                faceGlobalVertices[0], faceNormal,
                faceTransformMatrixInv.data(),
                boundaryFaceMatrix.data(),
                functor,
                leftContactFaceMatrix.data(),
                rightContactFaceMatrix.data(),
                time, dynamicContactType);

              functionSpace->template Decompose<GhostCellFunctionGetter<VolumeMeshT>, dimsCount>(functionGetter, ghostValues.data());

              timeDerivatives.noalias() += (faceDeformJacobian * invJacobian) * (faceTransformMatrix * xExteriorMatrix * 
                ghostValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral); */
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
                correspondingCellValues(valueIndex, functionIndex) =
                  useHalfStepSolutionForCorrespondingCell ?
                  halfStepCellSolutions[correspondingCellIndex].basisVectors[functionIndex].values[valueIndex] :
                  cellSolutions[correspondingCellIndex].basisVectors[functionIndex].values[valueIndex];
              }
            }

            system.BuildContactMatrices(interactionType, leftContactMatrix, rightContactMatrix);
            // for glue contact left matrix equals 0
            if (!leftContactMatrix.isZero(std::numeric_limits<Scalar>::epsilon()))
            {
              // interior side contribution
              timeDerivatives.noalias() += (faceDeformJacobian * invJacobian) * (faceTransformMatrix * xExteriorMatrix * leftContactMatrix.asDiagonal() * 
                faceTransformMatrixInv * currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral);
            }

            // exterior side contribution
            timeDerivatives.noalias() += (faceDeformJacobian * invJacobian) * (faceTransformMatrix * xExteriorMatrix * rightContactMatrix.asDiagonal() * 
              faceTransformMatrixInv * correspondingCellValues * 
              incomingFlux.srcFaces[faceNumber].dstFaces[correspondingFaceNumber].orientations[correspondingFaceOrientation].surfaceIntegral);
          }
        }
      }

      system.BuildXMatrix(cellMediumParameters[cellIndex], xMatrix);
      system.BuildYMatrix(cellMediumParameters[cellIndex], yMatrix);
      system.BuildZMatrix(cellMediumParameters[cellIndex], zMatrix);

      Vector refXDerivatives = GetRefXDerivatives(cellVertices) * invJacobian;
      Vector refYDerivatives = GetRefYDerivatives(cellVertices) * invJacobian;
      Vector refZDerivatives = GetRefZDerivatives(cellVertices) * invJacobian;

      xMixedMatrix.noalias() = xMatrix * refXDerivatives.x + yMatrix * refXDerivatives.y + zMatrix * refXDerivatives.z;
      yMixedMatrix.noalias() = xMatrix * refYDerivatives.x + yMatrix * refYDerivatives.y + zMatrix * refYDerivatives.z;
      zMixedMatrix.noalias() = xMatrix * refZDerivatives.x + yMatrix * refZDerivatives.y + zMatrix * refZDerivatives.z;

      timeDerivatives.noalias() -=
        xMixedMatrix * currCellValues * xDerivativeVolumeIntegrals +
        yMixedMatrix * currCellValues * yDerivativeVolumeIntegrals + 
        zMixedMatrix * currCellValues * zDerivativeVolumeIntegrals;

      typename SystemT::SourceFunctorT* sourceFunctor = system.GetSourceFunctor();
      if (sourceFunctor)
      {
        typedef SourceFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > SourceFunctorWrapper;
        SourceFunctorWrapper wrapper(sourceFunctor, time, this, cellIndex);
        functionSpace->template Decompose< SourceFunctorWrapper, dimsCount >(wrapper, sourceValues.data());
        timeDerivatives.noalias() -= sourceValues;
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
              sourcePointValues(valueIndex, functionIndex) = 
                values[valueIndex] * functionSpace->GetBasisFunctionValue(refPoint, functionIndex);
            }
          }
          timeDerivatives.noalias() -= sourcePointValues * cellVolumeIntegralsInv;
        } 
      }

      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          derivatives[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex]
            = -timeDerivatives(valueIndex, functionIndex);
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
