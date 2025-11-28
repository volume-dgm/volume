/* Simplified Task::Run() that (potentially (at some point)) 
    tests all steps and does not write snapshots */
template<typename Space, unsigned int order>
void Task<Space, order>::TestRun() 
{
  network->init(nullptr, nullptr);
  SetThreadsCount();
  network->setReceiveListener(this);

  ASSERT_TRUE(settingsFile);
  settings.ParseDirect(settingsFile);

  LoadNodesSchedule();

  for(IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber) 
  {
    solvers.push_back(SolversFactory<Scalar, IndexType>::Build(settings.solver.integrator));
  }

  distributedElasticMeshes.resize(GetCurrentNodeDomainsCount());
  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    distributedElasticMeshes[domainNumber] = new DistributedElasticMesh(
      solvers[domainNumber],
      settings.solver.tolerance,
      nodesSchedule[GetNodeId()].domainsIndices[domainNumber],
      GetDomainsCount(), settings.solver.hierarchyLevelsCount, 
      settings.solver.allowMovement,
      settings.solver.allowDiscreteDestruction || settings.solver.allowContinuousDestruction);
    distributedElasticMeshes[domainNumber]->allowPlasticity  = settings.solver.allowPlasticity;
    distributedElasticMeshes[domainNumber]->allowContinuousDestruction = settings.solver.allowContinuousDestruction;
    distributedElasticMeshes[domainNumber]->allowDiscreteDestruction   = settings.solver.allowDiscreteDestruction;

    distributedElasticMeshes[domainNumber]->updateCollisionInfoPeriod = settings.solver.updateCollisionInfoPeriod;

    distributedElasticMeshes[domainNumber]->erosion = settings.solver.erosion;
    distributedElasticMeshes[domainNumber]->dynamicContactBox = settings.solver.dynamicContactBox;
  }
  detectorsData.resize(GetCurrentNodeDomainsCount());

  LoadMeshes();
  BuildSourceTerms();
  BuildPointSources();
  BuildIniStateMakers();
  LoadInitialState();


  Scalar timeStep = settings.solver.maxTimeStep * settings.solver.maxScale;
  bool forceStep  = false;
  Scalar currTime = Scalar(0.0);

  // all meshes have the same phases count
  solverState.phasesCount = distributedElasticMeshes[0]->solver->GetPhasesCount();
  solverState.hierarchyLevelsCount = settings.solver.hierarchyLevelsCount;

  int lastInitialGlobalStep = -1;
  Scalar lastInitialCurrTime = Scalar(-1);
  const IndexType MaxHierarchyLevel = 1 << (settings.solver.hierarchyLevelsCount - 1);
  bool phaseAdvanced = true;

  printf("Entering main cycle\n");
  PrintSolverState(solverState, currTime, timeStep);

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    distributedElasticMeshes[domainNumber]->RebuildTimeHierarchyLevels(solverState.globalStepIndex, settings.solver.allowMovement);
    distributedElasticMeshes[domainNumber]->solver->InitStep(timeStep, settings.solver.tolerance, true);
    distributedElasticMeshes[domainNumber]->solver->InitStep(solverState);
    distributedElasticMeshes[domainNumber]->SetDamping(exp(-settings.solver.damping * timeStep));
  }

  while (currTime < settings.task.destinationTime || settings.task.destinationTime < 0)
  {
    if (phaseAdvanced && solverState.AllHierarchyLevelsCompleted())
    {
      if (lastInitialCurrTime < currTime && lastInitialGlobalStep < solverState.globalStepIndex)
      {
        lastInitialCurrTime = currTime;
        lastInitialGlobalStep = solverState.globalStepIndex;
      }
    }

    bool allPhasesCompleted;
    SolverState nextState = solverState.GetNext(&allPhasesCompleted);
    if (!solverState.IsUseful())
    {
      solverState = nextState;
      continue;
    }

    phaseAdvanced = true;
    for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber) 
    {
      phaseAdvanced = phaseAdvanced && distributedElasticMeshes[domainNumber]->solver->AdvancePhase(solverState);
    }

    if (!allPhasesCompleted || !phaseAdvanced)
    {
      SynchronizeMeshes();
      if (phaseAdvanced) 
      {
        solverState = nextState;
      } 
    } else 
    {
      bool globalStepSuccessful = true;

      for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
      {
        Scalar localError = distributedElasticMeshes[domainNumber]->solver->GetLastStepError();
        bool stepSuccessful = forceStep || (localError < Scalar(2.0));

        globalStepSuccessful = network->template Negotiate<bool, LogicSumComparator>(stepSuccessful, logicSumComparator);
      }

      if (globalStepSuccessful)
      {
        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
        {
          distributedElasticMeshes[domainNumber]->SetGlobalStepIndex(solverState.globalStepIndex);
          distributedElasticMeshes[domainNumber]->solver->AdvanceStep(solverState);
        }
        currTime = distributedElasticMeshes[0]->solver->GetCurrTime();

        if (solverState.IsPreInitial())
        {
          for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber) 
          {
            if (settings.solver.allowContinuousDestruction || settings.solver.allowDiscreteDestruction) FindDestructions(domainNumber);
            if (settings.mesh.moveMassCenter)    distributedElasticMeshes[domainNumber]->MoveSceneToSnapshotRegion();
            if (settings.solver.allowPlasticity) distributedElasticMeshes[domainNumber]->HandlePlasticity(
                                                 distributedElasticMeshes[domainNumber]->solver->GetCurrStep());
            if(settings.solver.damping > 0)      distributedElasticMeshes[domainNumber]->HandleDamping();
          }
        }

        SynchronizeMeshes();
        solverState = nextState;
      } else 
      {
        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
        {
          distributedElasticMeshes[domainNumber]->solver->RevertStep(lastInitialCurrTime);
        }
        solverState.SetInitialState(lastInitialGlobalStep);
        currTime = lastInitialCurrTime;
      }

      if (!globalStepSuccessful || solverState.AllHierarchyLevelsCompleted())
      {
        Scalar globalDesiredStep = std::numeric_limits<Scalar>::max() * Scalar(0.5);

        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber) 
        {
          Scalar localDesiredStep = distributedElasticMeshes[domainNumber]->solver->GetTimeStepPrediction();
          localDesiredStep = minValueComparator(globalDesiredStep, localDesiredStep);
          globalDesiredStep = network->template Negotiate<Scalar, MinValueComparator>(localDesiredStep, minValueComparator);
        }

        timeStep = globalDesiredStep;
        forceStep = false;

        if (timeStep > settings.solver.maxTimeStep) timeStep = settings.solver.maxTimeStep;
        if (timeStep < settings.solver.maxTimeStep * settings.solver.maxScale)
        {
          timeStep = settings.solver.maxTimeStep * settings.solver.maxScale;
          forceStep = true;
        }

        for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
        {
          if (settings.solver.allowContinuousDestruction || settings.solver.allowDiscreteDestruction || 
              settings.solver.allowMovement || settings.solver.allowPlasticity)
          {
            distributedElasticMeshes[domainNumber]->RebuildTimeHierarchyLevels(solverState.globalStepIndex, settings.solver.allowMovement);
          }
          distributedElasticMeshes[domainNumber]->solver->InitStep(timeStep, settings.solver.tolerance, globalStepSuccessful);
          distributedElasticMeshes[domainNumber]->SetDamping(exp(-settings.solver.damping * timeStep));
        }
      }

      for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
      {
        distributedElasticMeshes[domainNumber]->solver->InitStep(solverState);
      }
      PrintSolverState(solverState, currTime, timeStep);
      (void) GetPointSolution(Vector(Scalar(0.0), Scalar(0.0)));
    }
  }
  network->finalize();
  printf("reached end\n");
}

template<typename Space, unsigned int order>
void Task<Space, order>::CompareTestSolutionSnapshot(Scalar currTime, IndexType stepIndex, std::string solutionFile)
{
  int resolutionX = 100, resolutionY = 100;

  Scalar eps = Scalar(1e-3);
  Vector boxPoint1 = Vector(Scalar(-2.0) + eps, Scalar(-2.0) + eps);
  Vector boxPoint2 = Vector(Scalar(2.0) - eps, Scalar(2.0) - eps);

  std::vector<Elastic> sampleData;
  sampleData.resize(resolutionX*resolutionY);

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber) 
  {
    distributedElasticMeshes[domainNumber]->MakeSnapshot(sampleData.data(),
      resolutionX, resolutionY,
      boxPoint1, boxPoint2, false);
    
    Vector stepSize = Vector((boxPoint2.x - boxPoint1.x)/Scalar(resolutionX - 1), (boxPoint2.y - boxPoint1.y)/Scalar(resolutionY - 1));
    Scalar diff = 0;

    for (int y = 0; y < resolutionY; ++y)
    {
      for (int x = 0; x < resolutionX; ++x)
      {

      }
    }
  }

}

template<typename Space, unsigned int order>
typename ElasticSystem<Space>::ValueType Task<Space, order>::GetPointSolution(Vector poi, bool halfStepSolution)
{
  IndexType poiDomain = 0, poiCellIndex = 0;
  Vector poiVertices[Space::NodesPerCell];
  bool cellFound = false;

  for (IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    for (IndexType cellIndex = 0; cellIndex < distributedElasticMeshes[domainNumber]->volumeMesh.cells.size(); ++cellIndex) 
    {
      Vector points[Space::NodesPerCell];
      distributedElasticMeshes[domainNumber]->volumeMesh.GetCellVertices(cellIndex, points);
      if (PointInCell(points, poi)) 
      {
        distributedElasticMeshes[domainNumber]->volumeMesh.GetCellVertices(cellIndex, poiVertices);
        poiCellIndex = cellIndex;
        poiDomain = domainNumber;
        cellFound = true;
      }
    }
  }

  EXPECT_TRUE(cellFound);

  printf(" point found: domain %zu cellIndex %zu\n", poiDomain, poiCellIndex);

  Vector refCelCoords = distributedElasticMeshes[poiDomain]->volumeMesh.GlobalToRefVolumeCoords(poi, poiVertices);

  typename ElasticSystemType::ValueType result(Scalar(0.0));

  for(IndexType functionIndex = 0; functionIndex < distributedElasticMeshes[poiDomain]->volumeMesh.functionsCount; ++functionIndex)
  {
    Scalar basisFuncValue = distributedElasticMeshes[poiDomain]->volumeMesh.functionSpace->GetBasisFunctionValue(refCelCoords, functionIndex);
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex) 
    {
      Scalar basisFuncCoeff = 
        (halfStepSolution ? distributedElasticMeshes[poiDomain]->volumeMesh.halfStepCellSolutions[poiCellIndex].basisVectors[functionIndex].values[valueIndex] :
         distributedElasticMeshes[poiDomain]->volumeMesh.cellSolutions[poiCellIndex].basisVectors[functionIndex].values[valueIndex]);
        
      result.values[valueIndex] += basisFuncCoeff*basisFuncValue;
    }
  }

  printf("  solution stress XX at poi: %g\n", result.GetXX());

  return result;
}

