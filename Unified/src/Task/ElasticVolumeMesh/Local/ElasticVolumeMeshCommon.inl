template<typename Space, typename FunctionSpace>
ElasticVolumeMeshCommon<Space, FunctionSpace>::
  ElasticVolumeMeshCommon(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  DifferentialSystem<Scalar>(solver->GetPhasesCount(), hierarchyLevelsCount), volumeMesh(solver->GetPhasesCount(), hierarchyLevelsCount)
{
  this->solver = solver;
  this->tolerance = tolerance;
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::LoadState(const std::vector< IniStateMaker<ElasticSpaceType>* >& stateMakers, Scalar mult)
{
  #pragma omp parallel for
  for(int cellIndex = 0; cellIndex < int(volumeMesh.cells.size()); cellIndex++)
  {
    for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
    {
      for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++) 
      {
        volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] = Scalar(0.0);
      }
    }

    Scalar totalCellValues[functionsCount * dimsCount];
    std::fill(totalCellValues, totalCellValues + functionsCount * dimsCount, Scalar(0.0));

    typedef FunctionGetter<ElasticVolumeMeshCommon<Space, FunctionSpace>, IniStateMaker<ElasticSpaceType> > ElasticGetter;
    for (IndexType stateMakerIndex = 0; stateMakerIndex < stateMakers.size(); ++stateMakerIndex)
    {
      ElasticGetter getter(this, stateMakers[stateMakerIndex], cellIndex);
      Scalar cellValues[functionsCount * dimsCount];
      volumeMesh.functionSpace->template Decompose<ElasticGetter, dimsCount>(getter, cellValues);
      std::transform(cellValues, cellValues + functionsCount * dimsCount,
                     totalCellValues, totalCellValues, std::plus<Scalar>());
    }

    for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
    {
      for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++) 
      {
        volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] = 
          totalCellValues[valueIndex * functionsCount + functionIndex];
      }
    }
  }
}

template<typename Space, typename FunctionSpace>
typename ElasticVolumeMeshCommon<Space, FunctionSpace>::Elastic 
ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElasticRef(IndexType cellIndex, Vector refCoords) const
{
  return volumeMesh.GetRefCellSolution(cellIndex, refCoords);
}

template<typename Space, typename FunctionSpace>
typename ElasticVolumeMeshCommon<Space, FunctionSpace>::Elastic 
ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElastic(IndexType cellIndex, Vector globalPoint, bool halfStepSolution) const
{
  return volumeMesh.GetCellSolution(cellIndex, globalPoint, halfStepSolution);
}

template<typename Space, typename FunctionSpace>
typename ElasticVolumeMeshCommon<Space, FunctionSpace>::Elastic 
ElasticVolumeMeshCommon<Space, FunctionSpace>::GetAverageCellElastic(IndexType cellIndex) const
{
  return volumeMesh.GetCellAverageSolution(cellIndex);
}

template<typename Space, typename FunctionSpace>
int ElasticVolumeMeshCommon<Space, FunctionSpace>::GetDimentionsCount(const SolverState& solverState) const
{
  /* (velocity_x, velocity_y, sigma_xx, sigma_yy, sigma_xy) * cellsCount + 
     (position_x, position_y) * nodesCount */

  if(allowMovement && solverState.IsLastStep())
  {
    return volumeMesh.GetDimentionsCount(solverState) + volumeMesh.nodes.size() * Space::Dimension;
  } else
  {
    return volumeMesh.GetDimentionsCount(solverState);
  }
}

template<typename Space, typename FunctionSpace>
int ElasticVolumeMeshCommon<Space, FunctionSpace>::GetMaxDimentionsCount() const
{
  if(allowMovement)
  {
    return volumeMesh.GetMaxDimentionsCount() + volumeMesh.nodes.size() * Space::Dimension;
  } else
  {
    return volumeMesh.GetMaxDimentionsCount(); 
  }
}


template<typename Space, typename FunctionSpace>
void  ElasticVolumeMeshCommon<Space, FunctionSpace>::GetCurrCoords(Scalar& time, 
  Scalar* currCoords, Scalar* oldCoords, const SolverState& solverState)
{
  volumeMesh.GetCurrCoords(time, currCoords, oldCoords, solverState);

  if(allowMovement && solverState.IsLastStep())
  {    
    Vector* nodePositions = (Vector*)(currCoords + volumeMesh.GetDimentionsCount(solverState));
    for(IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      nodePositions[nodeIndex] = volumeMesh.nodes[nodeIndex].pos;
    }
  }
}

template<typename Space, typename FunctionSpace>
void  ElasticVolumeMeshCommon<Space, FunctionSpace>::GetCurrCoords(Scalar& time, Scalar* currCoords) const
{
  volumeMesh.GetCurrCoords(time, currCoords);

  if(allowMovement)
  {
    Vector* nodePositions = (Vector*)(currCoords + volumeMesh.GetMaxDimentionsCount());
    for(IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      nodePositions[nodeIndex] = volumeMesh.nodes[nodeIndex].pos;
    }
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState)
{
  volumeMesh.SetCurrCoords(time, newCoords, solverState);

  if(allowMovement && solverState.IsLastStep())
  {
    Vector* nodePositions = (Vector*)(newCoords + volumeMesh.GetDimentionsCount(solverState));
    for(IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      volumeMesh.nodes[nodeIndex].pos = nodePositions[nodeIndex];
    }

    UnfoldMesh(minGridHeight, unfoldIterationsCount);
    volumeMesh.UpdateAABBTree();
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetCurrCoords(Scalar time, 
  const Scalar* newCoords, const Scalar* oldCoords, const SolverState& solverState)
{
  volumeMesh.SetCurrCoords(time, newCoords, oldCoords, solverState);

  if(allowMovement && solverState.IsLastStep())
  {
    Vector* nodePositions = (Vector*)(newCoords + volumeMesh.GetDimentionsCount(solverState));
    for(IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      volumeMesh.nodes[nodeIndex].pos = nodePositions[nodeIndex];
    }

    UnfoldMesh(minGridHeight, unfoldIterationsCount);
    volumeMesh.UpdateAABBTree();
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetCurrCoords(Scalar time, const Scalar* oldCoords)
{
  volumeMesh.SetCurrCoords(time, oldCoords);

  if(allowMovement)
  {
    Vector* nodePositions = (Vector*)(oldCoords + volumeMesh.GetMaxDimentionsCount());
    for(IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      volumeMesh.nodes[nodeIndex].pos = nodePositions[nodeIndex];
    }
    volumeMesh.UpdateAABBTree();
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState)
{
  volumeMesh.GetCurrDerivatives(derivatives, solverState);

  if(allowMovement && solverState.IsLastStep())
  {
    Vector* nodeVelocities = (Vector*)(derivatives + volumeMesh.GetDimentionsCount(solverState));
    for(IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      nodeVelocities[nodeIndex] = Vector::zeroVector();
    }
    for(int cellIndex = 0; cellIndex < int(volumeMesh.cells.size()); cellIndex++)
    {
      IndexType cellIndices[Space::NodesPerCell];
      volumeMesh.GetFixedCellIndices(cellIndex, cellIndices);

      Vector cellVertices[Space::NodesPerCell];
      for(IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; nodeNumber++)
      {
        nodeVelocities[cellIndices[nodeNumber]] += InterpolateElastic(cellIndex, volumeMesh.nodes[cellIndices[nodeNumber]].pos).GetVelocity();
      }
    }

    for (IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      nodeVelocities[nodeIndex] /= Scalar(volumeMesh.GetIncidentCellsCount(nodeIndex));
    }
    
    if (allowDestruction)
    {
      std::fill(isNodeVelocityFound.begin(), isNodeVelocityFound.end(), false);

      std::vector<IndexType> nodeGroupPool;
      nodeGroupPool.reserve(16);

      for (IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); ++nodeIndex)
      {
        if (!isNodeVelocityFound[nodeIndex])
        {
          IndexType nodeGroupSize = volumeMesh.GetNodeGroup(nodeIndex, nodeGroupPool);

          Vector groupMeanVelocity = Vector::zero();
          for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSize; ++nodeNumber)
          {
            groupMeanVelocity += nodeVelocities[nodeGroupPool[nodeNumber]];
          }
          groupMeanVelocity /= Scalar(nodeGroupSize);
          for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSize; ++nodeNumber)
          {
            isNodeVelocityFound[nodeGroupPool[nodeNumber]] = true;
            nodeVelocities[nodeGroupPool[nodeNumber]] = groupMeanVelocity;
          }
        }
      }
    }
    
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetDetectors(const std::vector<Vector>& globalDetectorsPositions)
{
  detectorsPositions.resize(globalDetectorsPositions.size(), Vector::zeroVector());
  detectorsCells.resize(globalDetectorsPositions.size(), IndexType(-1));

  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    Vector cellVertices[Space::NodesPerCell];
    volumeMesh.GetCellVertices(cellIndex, cellVertices);
          
    for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
    {
      if (detectorsCells[detectorIndex] == IndexType(-1) && 
          PointInCell<Scalar>(cellVertices, globalDetectorsPositions[detectorIndex]))
      {
        detectorsCells[detectorIndex] = cellIndex;
        // mesh is movable -> detection position is not constant
        Vector refCoords = volumeMesh.GlobalToRefVolumeCoords(globalDetectorsPositions[detectorIndex], cellVertices);
        detectorsPositions[detectorIndex] = refCoords;
      }
    }
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::GetDetectorsData(IndexType detectorIndex, Scalar* data)
{
  std::fill(data, data + volumeMesh.dimsCount, Scalar(0.0));
  if (detectorsCells[detectorIndex] != IndexType(-1))
  {      
    IndexType cellIndex = detectorsCells[detectorIndex];
    Vector    refCoords = detectorsPositions[detectorIndex];
    Elastic     elastic = InterpolateElasticRef(cellIndex, refCoords);
    // elastic.MakeDimention(tensionDimensionlessMult, velocityDimensionlessMult);

    for (IndexType valueIndex = 0; valueIndex < volumeMesh.dimsCount; ++ valueIndex)
    {
      data[valueIndex] = elastic.values[valueIndex];
    }
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCellMaterial(IndexType cellIndex, Scalar alpha)
{
  if (!volumeMesh.cellMediumParameters[cellIndex].destroyed)
  {
    volumeMesh.cellMediumParameters[cellIndex].destroyed = true;
    /*
      from Chelnokov PhD work...

      mju` = alpha * mju 
      K = lambda + 2/3 * mju = const -> lambda`
    */
    volumeMesh.cellMediumParameters[cellIndex].lambda -= 2 * (1 - alpha) / 3 * volumeMesh.cellMediumParameters[cellIndex].mju;
    volumeMesh.cellMediumParameters[cellIndex].mju *= alpha;
  }
}

template<typename Space, typename FunctionSpace>
typename ElasticVolumeMeshCommon<Space, FunctionSpace>::ElasticSystemType* ElasticVolumeMeshCommon<Space, FunctionSpace>::GetSystem()
{
  return &(volumeMesh.system);
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::MoveSceneToSnapshotRegion()
{
  Vector massCenter = volumeMesh.GetMassCenter();
  #pragma omp parallel for
  for (int nodeIndex = 0; nodeIndex < int(volumeMesh.nodes.size()); ++nodeIndex)
  {
    volumeMesh.nodes[nodeIndex].pos -= massCenter;
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::Initialize(bool allowMovement,
  Scalar collisionWidth,
  IndexType unfoldIterationsCount, Scalar minGridHeight,
  Scalar tensionErrorMult,
  Scalar velocityErrorMult,
  Scalar positionErrorMult,
  Scalar tensionDimensionlessMult,
  Scalar velocityDimensionlessMult)
{
  volumeMesh.collisionWidth = collisionWidth;
  volumeMesh.allowDynamicCollisions = allowMovement;

  if (volumeMesh.allowDynamicCollisions)
  {
    volumeMesh.BuildAABBTree();
  }

  if (allowDestruction)
  {
    initialAdditionalCellInfos = volumeMesh.additionalCellInfos;
    isNodeVelocityFound.resize(volumeMesh.nodes.size());
  }

  this->allowMovement = allowMovement;

  this->unfoldIterationsCount = unfoldIterationsCount;
  this->minGridHeight = minGridHeight;

  this->tensionErrorMult = tensionErrorMult;
  this->velocityErrorMult = velocityErrorMult;
  this->positionErrorMult = positionErrorMult;

  this->tensionDimensionlessMult = tensionDimensionlessMult;
  this->velocityDimensionlessMult = velocityDimensionlessMult;

  ComputeElasticMults();
  solver->SetSystem(this);
}

template<typename Space, typename FunctionSpace>
typename Space::Scalar
  ElasticVolumeMeshCommon<Space, FunctionSpace>::GetErrorValue(
    Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState& solverState, const Scalar*)
{
  Scalar staticError = volumeMesh.GetErrorValue(time, coords0, coords1, solverState, elasticMults);

  if (allowMovement)
  {
    Scalar dynamicError = 0;
    Vector* node0Positions = (Vector*)(coords0 + volumeMesh.GetDimentionsCount(solverState));
    Vector* node1Positions = (Vector*)(coords1 + volumeMesh.GetDimentionsCount(solverState));

    #pragma omp parallel for
    for (int nodeIndex = 0; nodeIndex < int(volumeMesh.nodes.size()); nodeIndex++)
    {
      Scalar nodeError = (node0Positions[nodeIndex] - node1Positions[nodeIndex]).Len() * positionErrorMult;
      #pragma omp flush(dynamicError)
      if (nodeError > dynamicError)
      {
        #pragma omp critical
        {
          if (nodeError > dynamicError) dynamicError = nodeError;
        }
      }
    }
    return dynamicError + staticError;
  } else
  {
    return staticError;
  }
}

template<typename MeshType>
struct DampingCorrector
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  DampingCorrector(MeshType* const mesh, IndexType cellIndex): mesh(mesh)
  {
    this->cellIndex = cellIndex;
  }

  void operator()(const Vector& refPoint, Scalar* values) const
  {
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, refPoint);
    Vector  velocity = elastic.GetVelocity();
    elastic.SetVelocity(velocity * mesh->GetDamping());
    std::copy(elastic.values, elastic.values + mesh->dimsCount, values);
  }
  MeshType* const mesh;
  IndexType cellIndex;
};

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::HandleDamping()
{
  typedef ElasticVolumeMeshCommon<Space, FunctionSpace> ElasticVolumeMeshCommonType;
  typedef DampingCorrector< ElasticVolumeMeshCommonType > DampingCorrectorType;
  ApplyCorrector<ElasticVolumeMeshCommonType, DampingCorrectorType>(this);
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetDamping(Scalar damping)
{
  this->damping = damping;
}

template<typename Space, typename FunctionSpace>
typename Space::Scalar ElasticVolumeMeshCommon<Space, FunctionSpace>::GetDamping() const
{
  return damping;
}
