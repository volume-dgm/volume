template<typename Space, unsigned int order>
void Task<Space, order>::TestRunDecomposed()
{

  /* network initialization */
  network->init(nullptr, nullptr);
  SetThreadsCount();
  network->setReceiveListener(this);

  ASSERT_TRUE(settingsFile);
  settings.ParseDirect(settingsFile);

  SetupNodes();

}

template<typename Space, unsigned int order>
void Task<Space, order>::SetupNodes()
{
  LoadNodesSchedule();

  for(IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    solvers.push_back(SolversFactory<Scalar, IndexType>::Build(settings.solver.integrator));
    EXPECT_TRUE(solvers[domainNumber] != nullptr);
  }

  distributedElasticMeshes.resize(GetCurrentNodeDomainsCount());
  for(IndexType domainNumber = 0; domainNumber < GetCurrentNodeDomainsCount(); ++domainNumber)
  {
    distributedElasticMeshes[domainNumber] = new DistributedElasticMesh
    (
      solvers[domainNumber],
      settings.solver.tolerance,
      nodesSchedule[GetNodeId()].domainIndices[domainNumber],
      GetDomainsCount(), settings.solver.hierarchyLevelsCount,
      settings.solver.allowMovement,
      settings.solver.allowDiscreteDestruction || settings.solver.allowContinuousDestruction
    );
    distributedElasticMeshes[domainNumber]->allowPlasticity  = settings.solver.allowPlasticity;
    distributedElasticMeshes[domainNumber]->allowContinuousDestruction = settings.solver.allowContinuousDestruction;
    distributedElasticMeshes[domainNumber]->allowDiscreteDestruction   = settings.solver.allowDiscreteDestruction;
    distributedElasticMeshes[domainNumber]->updateCollisionInfoPeriod = settings.solver.updateCollisionInfoPeriod;
    distributedElasticMeshes[domainNumber]->erosion = settings.solver.erosion;
    distributedElasticMeshes[domainNumber]->dynamicContactBox = settings.solver.dynamicContactBox;
  }
  detectorsData.resize(GetCurrentNodeDomainsCount());
}

template<typename Space, unsigned int order>
void Task<Space, order>::IniState()
{
  LoadMeshes();
  BuildSourceTerms();
  BuildPointSources();
  BuildIniStateMakers();
  LoadInitialState();

  Scalar timeStep = settings.solver.maxTimeStep * settings.solver.maxScale;
  bool forceStep = false;
  Scalar currTime = Scalar(0.0);

  solverState.phasesCount = distributedElasticMeshes[0]->solver->GetPhasesCount();
  solverState.hierarchyLevelsCount = settings.solver.hierarchyLevelsCount;

  int lastInitialGlobalStep = -1;
  Scalar lastInitialCurrTime = Scalar(-1);
  const IndexType MaxHierarchyLevel = 1 << (settings.solver.hierarchyLevelsCount - 1);
  bool phaseAdvanced = true;

}