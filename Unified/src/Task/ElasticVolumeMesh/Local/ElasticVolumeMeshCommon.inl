template<typename Space, typename FunctionSpace>
ElasticVolumeMeshCommon<Space, FunctionSpace>::
  ElasticVolumeMeshCommon(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  DifferentialSystem<Scalar>(solver->GetPhasesCount(), hierarchyLevelsCount), volumeMesh(solver->GetPhasesCount(), hierarchyLevelsCount),
  nodeGroupManager(99) // TODO
{
  this->solver = solver;
  this->tolerance = tolerance;
  this->minHeightInMesh = 0;

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    precomputedBasisFunctions.basisFunctions[functionIndex] = 
      volumeMesh.functionSpace->GetBasisPolynomial(functionIndex);
    for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
    {
      IndexVector derivatives = IndexVector::zeroVector();
      derivatives[dimIndex] = 1;

      precomputedBasisFunctions.basisDerivativeFunctions[functionIndex][dimIndex] = 
        volumeMesh.functionSpace->GetBasisFunctionDerivative(derivatives, functionIndex);
    }
  }
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

    if (volumeMesh.cellMediumParameters[cellIndex].fixed) continue;

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
ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElasticRef(IndexType cellIndex, IndexType pointIndex, typename VolumeMeshTypeCommon::PointType pointType) const
{
  return volumeMesh.GetRefCellSolution(cellIndex, pointIndex, pointType);
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

    #pragma omp parallel for
    for(int nodeIndex = 0; nodeIndex < int(volumeMesh.nodes.size()); nodeIndex++)
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

    #pragma omp parallel for
    for(int nodeIndex = 0; nodeIndex < int(volumeMesh.nodes.size()); nodeIndex++)
    {
      nodePositions[nodeIndex] = volumeMesh.nodes[nodeIndex].pos;
    }
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::UpdateAABBTree(IndexType globalStepIndex)
{
  if (globalStepIndex % updateCollisionInfoPeriod == 0)
  {
    volumeMesh.UpdateAABBTree();
    volumeMesh.UpdateCollisionsInfo();
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState)
{
  volumeMesh.SetCurrCoords(time, newCoords, solverState);

  if(allowMovement && solverState.IsLastStep())
  {
    Vector* nodePositions = (Vector*)(newCoords + volumeMesh.GetDimentionsCount(solverState));

    #pragma omp parallel for
    for(int nodeIndex = 0; nodeIndex < int(volumeMesh.nodes.size()); nodeIndex++)
    {
      volumeMesh.nodes[nodeIndex].pos = nodePositions[nodeIndex];
    }
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

    #pragma omp parallel for
    for(int nodeIndex = 0; nodeIndex < int(volumeMesh.nodes.size()); nodeIndex++)
    {
      volumeMesh.nodes[nodeIndex].pos = nodePositions[nodeIndex];
    }

    UnfoldMesh(minGridHeight, unfoldIterationsCount);
    UpdateAABBTree(globalStepIndex);
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetCurrCoords(Scalar time, const Scalar* oldCoords)
{
  volumeMesh.SetCurrCoords(time, oldCoords);

  if(allowMovement)
  {
    Vector* nodePositions = (Vector*)(oldCoords + volumeMesh.GetMaxDimentionsCount());
    #pragma omp parallel for
    for(int nodeIndex = 0; nodeIndex < int(volumeMesh.nodes.size()); nodeIndex++)
    {
      volumeMesh.nodes[nodeIndex].pos = nodePositions[nodeIndex];
    }

    UnfoldMesh(minGridHeight, unfoldIterationsCount);
    UpdateAABBTree(globalStepIndex);
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState)
{
  volumeMesh.GetCurrDerivatives(derivatives, solverState);

  if(allowMovement && solverState.IsLastStep())
  {
    Vector* nodeVelocities = (Vector*)(derivatives + volumeMesh.GetDimentionsCount(solverState));

    std::fill_n(nodeVelocities, volumeMesh.nodes.size(), Vector::zeroVector());

    for(int cellIndex = 0; cellIndex < int(volumeMesh.cells.size()); cellIndex++)
    {
      if (!volumeMesh.cellMediumParameters[cellIndex].fixed)
      {
        IndexType cellIndices[Space::NodesPerCell];
        volumeMesh.GetFixedCellIndices(cellIndex, cellIndices);

        for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; nodeNumber++)
        {
          Vector v = InterpolateElasticRef(cellIndex, nodeNumber, VolumeMeshTypeCommon::CellNode).GetVelocity();
          // Vector v = InterpolateElastic(cellIndex, volumeMesh.nodes[cellIndices[nodeNumber]].pos).GetVelocity();

          nodeVelocities[cellIndices[nodeNumber]] += v;
        }
      }
    }

    for (IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); nodeIndex++)
    {
      nodeVelocities[nodeIndex] /= Scalar(volumeMesh.GetIncidentCellsCount(nodeIndex));
    }

    if (allowDiscreteDestruction || allowContinuousDestruction)
    {
      for (IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); ++nodeIndex)
      {
        IndexType nodeGroupSize = nodeGroupManager.GetGroupSize(nodeIndex);

        if (nodeGroupSize > 0)
        {
          Vector groupMeanVelocity = Vector::zero();
          const IndexType* nodeGroupPool = nodeGroupManager.GetGroup(nodeIndex);

          for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSize; ++nodeNumber)
          {
            groupMeanVelocity += nodeVelocities[nodeGroupPool[nodeNumber]];
          }
          groupMeanVelocity /= Scalar(nodeGroupSize);
          for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSize; ++nodeNumber)
          {
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

  if (allowContinuousDestruction || allowDiscreteDestruction)
  {
    isCellBroken.resize(volumeMesh.cells.size());
    nodeGroupManager.SetVolumeMesh(&volumeMesh);
    nodeGroupManager.BuildNodeGroups();
  }

  if (allowPlasticity || allowContinuousDestruction)
  {
    plasticDeforms.resize(volumeMesh.cells.size());
  }

  minHeightInMesh = volumeMesh.GetMinHeight();

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

  void operator()(IndexType basisPointIndex, Scalar* values) const
  {
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, basisPointIndex, MeshType::VolumeMeshTypeCommon::Basis);
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

template<typename Space, typename FunctionSpace>
bool ElasticVolumeMeshCommon<Space, FunctionSpace>::ProcessPlasticity(const Scalar k, Elastic& elastic, bool updateElastic)
{
  const Scalar pressure = elastic.GetPressure();

  Tensor tension = elastic.GetTension();
  tension += pressure;
  Scalar ss = DoubleConvolution(tension, tension);

  bool plasticMode = ss > 2 * Sqr(k);

  if (plasticMode && updateElastic)
  {
    tension *= sqrt(2 * Sqr(k) / ss);
    tension += -pressure;
    elastic.SetTension(tension);
  }
  return plasticMode;
}

/*
s:s < 2k^2 -- von Mises criterion;
k = k0 - a * pressure, a > 0.
*/
template<typename MeshType>
struct PlasticityCorrector
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::Tensor    Tensor;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  PlasticityCorrector(MeshType* const mesh, IndexType cellIndex):
    mesh(mesh), cellIndex(cellIndex)
  {}

  void operator()(IndexType basisPointIndex, Scalar* values) const
  {
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, basisPointIndex, MeshType::VolumeMeshTypeCommon::Basis);
    const bool brittle = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
    if (!brittle)
    {
      const Scalar k0 = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.k0;
      const Scalar  a = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.a;
      const Scalar pressure = elastic.GetPressure();
      const Scalar k = k0 + a * pressure;
    
      mesh->ProcessPlasticity(k, elastic, true);
    }
    std::copy(elastic.values, elastic.values + mesh->dimsCount, values);
  }

  MeshType* const mesh;
  IndexType cellIndex;
};

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::HandlePlasticity(Scalar dt)
{
  if (allowContinuousDestruction)
  {
    #pragma omp parallel for
    for (int cellIndex = 0; cellIndex < int(volumeMesh.cells.size()); ++cellIndex)
    {
      Scalar plasticDeformRate = 0;

      Vector cellVertices[Space::NodesPerCell];
      volumeMesh.GetCellVertices(cellIndex, cellVertices);
      Scalar jacobian = volumeMesh.GetCellDeformJacobian(cellVertices);

      const Scalar k0 = volumeMesh.cellMediumParameters[cellIndex].plasticity.k0;
      const Scalar  a = volumeMesh.cellMediumParameters[cellIndex].plasticity.a;

      /*for (IndexType pointIndex = 0; pointIndex < volumeMesh.quadraturePoints.size(); ++pointIndex)
      {
        const bool brittle = volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
        if (!brittle)
        {
          Elastic elastic = volumeMesh.GetRefCellSolution(cellIndex, pointIndex, VolumeMeshTypeCommon::Quadrature);
          // Elastic elastic = volumeMesh.GetRefCellSolution(cellIndex, volumeMesh.quadraturePoints[pointIndex]);
          const Scalar pressure = elastic.GetPressure();
          const Scalar k = k0 + a * pressure;

          if (ProcessPlasticity(k, elastic, false))
          {
            plasticDeformRate += GetDeformRate(cellIndex, volumeMesh.quadraturePoints[pointIndex]) *
              volumeMesh.quadratureWeights[pointIndex] *
              jacobian;
          }
        }
      }
      plasticDeforms[cellIndex] += plasticDeformRate * dt / volumeMesh.GetVolume(cellIndex);
      */      
      /*auto deformFunc = [&](Vector point, Scalar *pointDeformRate) -> void
      {
        *pointDeformRate = GetDeformRate(cellIndex, point);
      };
      volumeMesh.functionSpace->GetCellAverage<decltype(deformFunc), 1>(func, &plasticDeformRate);*/ //works the same
      plasticDeformRate = volumeMesh.functionSpace->GetCellAverage(
        [&](Vector refPoint, Scalar *pointDeformRate) -> void
        {
          *pointDeformRate = 0;
          const bool brittle = volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
          if (!brittle)
          {
            Elastic elastic = volumeMesh.GetRefCellSolution(cellIndex, refPoint);
            //Elastic elastic = volumeMesh.GetRefCellSolution(cellIndex, pointIndex, VolumeMeshTypeCommon::Quadrature);
            const Scalar pressure = elastic.GetPressure();
            const Scalar k = k0 + a * pressure;

            if (ProcessPlasticity(k, elastic, false))
            {
              *pointDeformRate = GetDeformRate(cellIndex, refPoint);
              //*pointDeformRate = GetDeformRate(cellIndex, pointIndex);
            }
          }
        });

      plasticDeforms[cellIndex] += plasticDeformRate * dt; //* jacobian / volumeMesh.GetVolume(cellIndex) -- possibly needed, but i'm not sure

      if (plasticDeforms[cellIndex] > volumeMesh.cellMediumParameters[cellIndex].plasticity.maxPlasticDeform)
      {
        DestroyCellMaterial(cellIndex, volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
      }
    }
  }
  typedef PlasticityCorrector< ElasticVolumeMeshCommon<Space, FunctionSpace> > PlasticityCorrectorType;
  ApplyCorrector< ElasticVolumeMeshCommon<Space, FunctionSpace>, PlasticityCorrectorType>(this);
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::MakeRhoCorrection(const Scalar minTimeStep)
{
  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    if (volumeMesh.isCellAvailable[cellIndex] && !volumeMesh.cellMediumParameters[cellIndex].fixed)
    {
      Scalar minHeight = volumeMesh.GetMinHeight(cellIndex);
      Scalar currTimeStep = minHeight / volumeMesh.cellMediumParameters[cellIndex].GetPSpeed();

      if (currTimeStep < minTimeStep)
      {
        volumeMesh.cellMediumParameters[cellIndex].invRho *= sqrt(currTimeStep / minTimeStep);
      }
    }
    /*
    Scalar currentVolume = volumeMesh.GetVolume(cellIndex);
    Scalar initialMass = volumeMesh.cellMediumParameters[cellIndex].initialVolume * volumeMesh.cellMediumParameters[cellIndex].initialRho;
    volumeMesh.cellMediumParameters[cellIndex].invRho = std::max(currentVolume / initialMass, Scalar(1.0) / volumeMesh.cellMediumParameters[cellIndex].initialRho);
    */
  }
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::HandleMaterialErosion()
{
  // MakeRhoCorrection();

  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    bool lowDensity =  volumeMesh.GetVolume(cellIndex) > erosion.rhoReduction * volumeMesh.cellMediumParameters[cellIndex].initialVolume;
    if (volumeMesh.isCellAvailable[cellIndex] && 
      !volumeMesh.cellMediumParameters[cellIndex].fixed &&
      (volumeMesh.GetAspectRatio(cellIndex) > erosion.cellAspectRatio || 
       volumeMesh.GetMinHeight(cellIndex) < erosion.minHeightRatio * minHeightInMesh ||
      (allowPlasticity && plasticDeforms[cellIndex] > erosion.maxPlasticDeformation) ||
      lowDensity ||
      !erosion.boxOfInterest.Includes(volumeMesh.GetCellAABB(cellIndex))))
    {
      volumeMesh.cellSolutions[cellIndex].SetToZero();

      /*
      const IndexType MaxContactCellsCount = 1024;
      IndexType contactCells[MaxContactCellsCount];
      IndexType contactCellsCount = 0;
      
      Elastic elastic = GetAverageCellElastic(cellIndex);
      if (elastic.GetPressure() > 0 && !lowDensity)
      {
        ContactFinder< VolumeMeshType > contactFinder(&volumeMesh, cellIndex);
        contactFinder.Find(contactCells, &contactCellsCount);

        for (IndexType contactCellNumber = 0; contactCellNumber < contactCellsCount; ++contactCellNumber)
        {
          IndexType contactCellIndex = contactCells[contactCellNumber];

          if (!volumeMesh.cellMediumParameters[contactCellIndex].fixed)
          {
            volumeMesh.cellMediumParameters[contactCellIndex].invRho =
              Scalar(1.0) / (1 / volumeMesh.cellMediumParameters[contactCellIndex].invRho +
                1 / volumeMesh.cellMediumParameters[cellIndex].invRho / Scalar(contactCellsCount) *
                volumeMesh.GetVolume(cellIndex) / volumeMesh.GetVolume(contactCellIndex));
            // TODO: also splash potential energy of this cell between neighbours
          }
        }
      } */
      
      volumeMesh.RemoveCellFromAABBTree(cellIndex);
      volumeMesh.isCellAvailable[cellIndex] = false;

      for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
      {
        IndexType contactTypeIndex = volumeMesh.GetInteractionType(cellIndex, faceNumber);
        IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);

        if (contactTypeIndex == IndexType(-1) || correspondingCellIndex == IndexType(-1)) continue;
        IndexType dynamicBoundaryType = volumeMesh.system.GetContactDynamicBoundaryType(contactTypeIndex);
        if (dynamicBoundaryType != IndexType(-1))
        {
          if (volumeMesh.isCellAvailable[correspondingCellIndex])
          {
            //  volumeMesh.AddToAABBTree(correspondingCellIndex);
          }

          DestroyFace(cellIndex, faceNumber, dynamicBoundaryType);
        }
      }

      nodeGroupManager.RemoveCell(cellIndex);
    }
  }
}

template<typename Space, typename FunctionSpace>
typename Space::Scalar ElasticVolumeMeshCommon<Space, FunctionSpace>::GetDeformRate(IndexType cellIndex, Vector refCoords) const
{
  Scalar strainRate[Space::Dimension * Space::Dimension] = { 0 };
  
  Vector cellVertices[Space::NodesPerCell];
  volumeMesh.GetCellVertices(cellIndex, cellVertices);

  Vector refDerivatives[Space::Dimension];
  volumeMesh.GetRefDerivatives(cellVertices, refDerivatives);

  /*for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
  {
    for (IndexType functionIndex = 0; functionIndex < FunctionSpace::functionsCount; ++functionIndex)
    {
      Scalar basisFunctionDerivativeValue = precomputedBasisFunctions.basisDerivativeFunctions[functionIndex][dimIndex].GetValue(&(refCoords.x));
      for (IndexType i = 0; i < Space::Dimension; ++i)
      {
        for (IndexType j = 0; j < Space::Dimension; ++j)
        {
          strainRate[i * Space::Dimension + j] += Scalar(0.5) * basisFunctionDerivativeValue * (
            volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[dimsCount - Space::Dimension + i] * refDerivatives[dimIndex][j] +
            volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[dimsCount - Space::Dimension + j] * refDerivatives[dimIndex][i] );
        }
      }
    }
  }*/
  Vector basisFunctionGradients[functionsCount];
  for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
  {
    for (IndexType functionIndex = 0; functionIndex < FunctionSpace::functionsCount; ++functionIndex)
    {
      basisFunctionGradients[functionIndex][dimIndex] = precomputedBasisFunctions.basisDerivativeFunctions[functionIndex][dimIndex].GetValue(&(refCoords.x));
    }
  } 

  Space::AsymmetricTensor refVelocityGradient;

  for (IndexType velocityComponentIndex = 0; velocityComponentIndex < Space::Dimension; ++velocityComponentIndex)
  {
    refVelocityGradient[velocityComponentIndex] = Vector::zero();
    for (IndexType functionIndex = 0; functionIndex < FunctionSpace::functionsCount; ++functionIndex)
    {
      refVelocityGradient[velocityComponentIndex] += 
        basisFunctionGradients[functionIndex] * volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].GetVelocity()[velocityComponentIndex];
    }
  }

  Space::AsymmetricTensor velocityGradient;
  for (IndexType velocityComponentIndex = 0; velocityComponentIndex < Space::Dimension; ++velocityComponentIndex)
  {
    velocityGradient[velocityComponentIndex] = Vector::zero();
    for (IndexType dimIndex = 0; dimIndex < Space::Dimension; dimIndex++)
    {
      velocityGradient[velocityComponentIndex] +=
        refDerivatives[dimIndex] * refVelocityGradient[velocityComponentIndex][dimIndex];
    }
  }

  Space::AsymmetricTensor strainRateTensor;
  for (int i = 0; i < Space::Dimension; ++i)
  {
    for (int j = i; j < Space::Dimension; ++j)
    {
      strainRateTensor[i][j] =
        (velocityGradient[i][j] + velocityGradient[j][i]) * Scalar(0.5);
    }
  }
  for (int i = 0; i < Space::Dimension; ++i)
  {
    strainRateTensor[i][i] -= Scalar(1.0) / Scalar(Space::Dimension) * velocityGradient[i][i];
  }

  /*Scalar effectiveStrain = 0;

  for (int i = 0; i < Space::Dimension; ++i)
  {
    for (int j = i; j < Space::Dimension; ++j)
    {
      effectiveStrain += Sqr(strainRate[i * Space::Dimension + j]);
    }
  }

  return pow(Scalar(2.0 / 3) * effectiveStrain, Scalar(0.5));*/
  Scalar effectiveStrain = 0;

  for (int i = 0; i < Space::Dimension; ++i)
  {
    for (int j = i; j < Space::Dimension; ++j)
    {
      effectiveStrain += Sqr(strainRateTensor[i][j]);
    }
  }

  return sqrt(effectiveStrain);
}

template<typename Space, typename FunctionSpace>
typename Space::Scalar ElasticVolumeMeshCommon<Space, FunctionSpace>::GetTimeStepPrediction()
{
  return volumeMesh.GetTimeStepPrediction();
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::SetGlobalStepIndex(IndexType globalStepIndex)
{
  this->globalStepIndex = globalStepIndex;
}
