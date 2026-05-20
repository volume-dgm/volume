#pragma once

#include <limits>
#include <algorithm>
#include "../../ElasticSystem/ElasticSystem.h"
#include "../../ElasticSystem/IniStates.h"
#include "../../../IO/Vtk/BasicVtkWriter.h"
#include "../../../Maths/Spaces.h"
#include "../../VolumeMethod/VolumeMesh/VolumeMesh.h"
#include "../../Task/SettingsParser/SolverSettingsParser.h"
#include "../../GeomMesh/NodeGroupManager.h"

template<typename Space, typename FunctionSpace>
struct ElasticVolumeMeshCommon: public DifferentialSystem<typename Space::Scalar>
{
  SPACE_TYPEDEFS

  typedef          ElasticSystem<Space>                                ElasticSystemType;
  typedef typename ElasticSystemType::ElasticSpace                     ElasticSpaceType;
  typedef typename ElasticSystemType::MediumParameters                 MediumParameters;
  typedef          VolumeMesh<Space, FunctionSpace, ElasticSystemType> VolumeMeshType;
  typedef          VolumeMeshCommon<Space, FunctionSpace, ElasticSystemType> VolumeMeshTypeCommon;
  typedef typename ElasticSpaceType::Elastic                           Elastic;
  typedef typename VolumeMeshType::Node                                Node;
  typedef typename VolumeMeshType::Cell                                Cell;
  typedef typename VolumeMeshType::CellSolution                        CellSolution;
  typedef          Space                                               SpaceT;

  const static int functionsCount = FunctionSpace::functionsCount;
  const static int dimsCount      = VolumeMeshType::dimsCount;

  ElasticVolumeMeshCommon(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  std::vector<Vector> detectorsPositions;
  std::vector<IndexType> detectorsCells;

  VolumeMeshType volumeMesh;
  DifferentialSolver<Scalar>* solver;
  Scalar tolerance;

  IndexType unfoldIterationsCount;
  Scalar minGridHeight;

  Scalar tensionErrorMult;
  Scalar velocityErrorMult;
  Scalar positionErrorMult;

  Scalar tensionDimensionlessMult;
  Scalar velocityDimensionlessMult;

  bool allowMovement;
  bool allowPlasticity;
  bool allowContinuousDestruction;
  bool allowDiscreteDestruction;

  IndexType updateCollisionInfoPeriod;

  void Initialize(bool allowMovement,
    Scalar collisionWidth,
    IndexType unfoldIterationsCount, Scalar minGridHeight,
    Scalar tensionErrorMult,
    Scalar velocityErrorMult,
    Scalar positionErrorMult,
    Scalar tensionDimensionlessMult,
    Scalar velocityDimensionlessMult);

  void LoadState(const std::vector< IniStateMaker<ElasticSpaceType>* >& stateMakers, Scalar mult);

  Elastic InterpolateElasticRef(IndexType cellIndex, Vector refCoords) const;
  Elastic InterpolateElasticRef(IndexType cellIndex, IndexType pointIndex, typename VolumeMeshTypeCommon::PointType pointType) const;
  Elastic InterpolateElastic(IndexType cellIndex, Vector globalPoint, bool halfStepSolution = false) const;

  Elastic GetAverageCellElastic(IndexType cellIndex) const;
  
  int GetDimentionsCount(const SolverState& solverState) const override;
  int GetMaxDimentionsCount() const override;

  void GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState) override;
  void GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState& solverState) override;
  void GetCurrCoords(Scalar& time, Scalar* currCoords) const override;

  void SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState) override;
  void SetCurrCoords(Scalar time, const Scalar* newCoords, const Scalar* oldCoords, const SolverState& solverState) override;
  void SetCurrCoords(Scalar time, const Scalar* oldCoords) override;
  
  void SetDetectors(const std::vector<Vector>& globalDetectorsPositions);
  void GetDetectorsData(IndexType detectorIndex, Scalar* data);

  Scalar GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState&, const Scalar* mults) override;

  void FindDestructions(std::vector<bool>* isCellBroken);
  virtual void UnfoldMesh(Scalar minHeight, IndexType iterationsCount) = 0;
  virtual void HandleContinuousDestruction() = 0;

  void HandleDamping();
  void SetDamping(Scalar damping);
  Scalar GetDamping() const;

  void DestroyCell (IndexType cellIndex, Scalar powderShearMult);
  ElasticSystemType* GetSystem();
  void MoveSceneToSnapshotRegion();

  /* 
    sqrt(2/3 * strain_{ij} : strain_{ij})
  */
  Scalar GetDeformRate(IndexType cellIndex, Vector refCoords) const;

  // stress correction due to plasticity in compliance with  Prandtl-Ruess model
  void HandlePlasticity(Scalar dt);

  typename SolverSettings<Space>::Erosion erosion;
  typename SolverSettings<Space>::DynamicContactBox dynamicContactBox;

  void   GetCellEnergy(IndexType cellIndex, Scalar& kineticEnergy, Scalar& potentialEnergy) const;
  Vector GetTotalImpulse() const;
  Scalar GetTotalEnergy(Scalar& kineticEnergy, Scalar& potentialEnergy) const;
  Scalar GetTotalMass() const;

  void HandleMaterialErosion();

  virtual Scalar GetTimeStepPrediction() override;

  virtual void DestroyFace(IndexType cellIndex, IndexType faceNumber, IndexType dynamicBoundaryType) = 0;

  // total work of the plasticity for each cell
  std::vector<Scalar> plasticDeforms;

  std::vector<bool> isCellBroken;

  bool ProcessPlasticity(const Scalar k, Elastic& elastic, bool updateElastic = false);

  void MakeRhoCorrection(Scalar minTimeStep = 2e-5);

  void SetGlobalStepIndex(IndexType globalStepIndex);

protected:
  void ComputeElasticMults()
  {
    ComputeElasticMults(Overload<Space>());
  }

  IndexType globalStepIndex;

  void UpdateAABBTree(IndexType globalStepIndex);

  NodeGroupManager<Space, VolumeMeshType> nodeGroupManager;

  // for error computation
  Scalar elasticMults[dimsCount];
  Scalar damping;
  Scalar minHeightInMesh;

  struct PrecomputedBasisFunctions
  {
    Polynomial<Scalar, IndexType, Space::Dimension> basisFunctions[functionsCount];
    Polynomial<Scalar, IndexType, Space::Dimension> basisDerivativeFunctions[functionsCount][Space::Dimension];
  } precomputedBasisFunctions;

private:
  void ComputeElasticMults(Overload<Space2>)
  {
    for (int valueIndex = 0; valueIndex < volumeMesh.dimsCount; ++valueIndex)
    {
      Scalar mult = Scalar(1.0);
      switch (valueIndex)
      {
        case 0: case 1: case 2:
        {
          mult = tensionErrorMult;
        } break;
        case 3: case 4:
        {
          mult = velocityErrorMult;
        } break;
      }
      elasticMults[valueIndex] = mult;
    }
  }

  void ComputeElasticMults(Overload<Space3>)
  {
    for (int valueIndex = 0; valueIndex < volumeMesh.dimsCount; ++valueIndex)
    {
      Scalar mult = Scalar(1.0);
      switch (valueIndex)
      {
        case 0: case 1: case 2: case 3: case 4: case 5:
        {
          mult = tensionErrorMult;
        } break;
        case 6: case 7: case 8:
        {
          mult = velocityErrorMult;
        } break;
      }
      elasticMults[valueIndex] = mult;
    }
  }
protected:
  // for test only
  bool CheckNodeGroupInfo(IndexType nodeIndex)
  {
    IndexType nodeGroupSize = nodeGroupManager.GetGroupSize(nodeIndex);

    if (nodeGroupSize == 0) return true;

    std::vector<IndexType> pool;
    pool.reserve(60);

    int test[99] = { 0 };
    IndexType size = volumeMesh.GetNodeGroup(nodeIndex, pool);
    
    if (size != nodeGroupSize)
    {
      std::cout << nodeIndex << "\n";
      return false;
    }

    const IndexType* nodeGroupP = nodeGroupManager.GetGroup(nodeIndex);
    std::copy(nodeGroupP, nodeGroupP + size, test);

    std::sort(pool.begin(), pool.end());
    std::sort(test, test + size);

    for (int i = 0; i < size; ++i)
    {
      if (test[i] != pool[i])
      {
        return false;
      }
    }
    return true;
  }

  bool CheckFace(IndexType cellIndex, IndexType faceNumber)
  {
    IndexType faceIndices[Space::NodesPerFace];
    volumeMesh.GetCellFaceNodes(cellIndex, faceNumber, faceIndices);

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerFace; ++nodeNumber)
    {
      if (!CheckNodeGroupInfo(faceIndices[nodeNumber])) return false;
    }
    return true;
  }

  bool CheckCell(IndexType cellIndex)
  {
    IndexType cellIndices[Space::NodesPerCell];
    volumeMesh.GetFixedCellIndices(cellIndex, cellIndices);

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (!CheckNodeGroupInfo(cellIndices[nodeNumber])) return false;
    }
    return true;
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename Space, typename FunctionSpace>
struct ElasticVolumeMesh;

template<typename FunctionSpace>
struct ElasticVolumeMesh<Space2, FunctionSpace>: public ElasticVolumeMeshCommon<Space2, FunctionSpace>
{
  SPACE2_TYPEDEFS
  typedef          Space2                                                  Space;
  typedef          ElasticSystem<Space>                                    ElasticSystemType;
  using ElasticSpaceType = ElasticSystemType::ElasticSpace;
  using Elastic = ElasticSpaceType::Elastic;
  typedef          VolumeMesh<Space, FunctionSpace, ElasticSystemType>     VolumeMeshType;
  typedef typename VolumeMeshType::EdgePairIndices                         EdgePairIndices;
  typedef typename VolumeMeshType::BoundaryEdge                            BoundaryEdge;
  typedef typename VolumeMeshType::EdgeLocation                            EdgeLocation;
  typedef typename VolumeMeshType::EdgeLocationPair                        EdgeLocationPair;
  typedef typename VolumeMeshType::EdgeIndices                             EdgeIndices;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::tensionDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::velocityDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetAverageCellElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetSystem;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCell;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::erosion;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::minHeightInMesh;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::HandleMaterialErosion;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowContinuousDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowDiscreteDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::FindDestructions;

  ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  virtual void MakeSnapshot(Elastic* destData,
    IndexType width, IndexType height, 
    Vector boxPoint1, Vector boxPoint2, bool halfStepSolution);

  virtual void MakeSnapshot(Elastic* destData,
    const Vector& origin, const Vector& spacing, 
    const Vector& boxPoint1, const Vector& boxPoint2,
    bool halfStepSolution);

  Elastic GetAverageEdgeElastic (IndexType cellIndex, IndexType edgeNumber) const;

  void HandleContinuousDestruction() override;

  void DestroyFace(IndexType cellIndex, IndexType edgeNumber, IndexType dynamicBoundaryType) override
  {
    IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, edgeNumber);
    IndexType correspondingEdgeNumber = volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber;

    volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType = dynamicBoundaryType;

    if (correspondingCellIndex != IndexType(-1) && correspondingEdgeNumber != IndexType(-1))
    {
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringEdges[correspondingEdgeNumber].correspondingCellIndex = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringEdges[correspondingEdgeNumber].correspondingEdgeNumber = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringEdges[correspondingEdgeNumber].interactionType = dynamicBoundaryType;
    }
    this->nodeGroupManager.UpdateFace(cellIndex, edgeNumber);
  }

private:
  void UnfoldMesh(Scalar minHeight, IndexType iterationsCount) override;
};

template <typename FunctionSpace>
struct ElasticVolumeMesh<Space3, FunctionSpace>: public ElasticVolumeMeshCommon<Space3, FunctionSpace>
{
  SPACE3_TYPEDEFS
  typedef          Space3                                       Space;
  typedef          ElasticSystem<Space3>                        ElasticSystemType;
  using ElasticSpaceType = ElasticSystemType::ElasticSpace;
  using MediumParameters = ElasticSystemType::MediumParameters;

  typedef          VolumeMesh<Space3, FunctionSpace, ElasticSystemType>    VolumeMeshType;
  typedef typename VolumeMeshType::FacePairIndices                         FacePairIndices;
  typedef typename VolumeMeshType::BoundaryFace                            BoundaryFace;
  typedef typename VolumeMeshType::Node                                    Node;
  typedef typename VolumeMeshType::Cell                                    Cell;
  typedef typename VolumeMeshType::CellSolution                            CellSolution;
  using Elastic = ElasticSpaceType::Elastic;
  typedef typename VolumeMeshType::FaceLocation                            FaceLocation;
  typedef typename VolumeMeshType::FaceLocationPair                        FaceLocationPair;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::tensionDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::velocityDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetAverageCellElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetSystem;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCell;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::HandlePlasticity;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::HandleMaterialErosion;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowContinuousDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::allowDiscreteDestruction;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::FindDestructions;

  ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  void UnfoldMesh(Scalar minHeight, IndexType iterationsCount) override;

  virtual void MakeSnapshot(Elastic* destData,
    const Vector& origin, const Vector& spacing, 
    const Vector& boxPoint1, const Vector& boxPoint2,
    bool halfStepSolution);

  void HandleContinuousDestruction() override;

  void DestroyFace(IndexType cellIndex, IndexType faceNumber, IndexType dynamicBoundaryType) override
  {
    IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);
    IndexType correspondingFaceNumber = volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber;

    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].orientation = IndexType(-1);
    volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType = dynamicBoundaryType;

    if (correspondingCellIndex != IndexType(-1) && correspondingFaceNumber != IndexType(-1))
    {
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].correspondingCellIndex = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].correspondingFaceNumber = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].orientation = IndexType(-1);
      volumeMesh.additionalCellInfos[correspondingCellIndex].neighbouringFaces[correspondingFaceNumber].interactionType = dynamicBoundaryType;
    }
    this->nodeGroupManager.UpdateFace(cellIndex, faceNumber);
  }

};

template<typename MeshType, typename StateMaker>
struct FunctionGetter
{
  typedef typename MeshType::SpaceT Space;
  SPACE_TYPEDEFS

  FunctionGetter(MeshType* mesh, StateMaker* stateMaker, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->stateMaker = stateMaker;
    this->cellIndex = cellIndex;
    mesh->volumeMesh.GetCellVertices(cellIndex, cellVertices);
  }

  void operator()(const Vector& point, Scalar* values)
  {
    Vector globalPoint = mesh->volumeMesh.RefToGlobalVolumeCoords(point, cellVertices);

    Scalar lambda = mesh->volumeMesh.cellMediumParameters[cellIndex].lambda;
    Scalar mju    = mesh->volumeMesh.cellMediumParameters[cellIndex].mju;
    Scalar invRho = mesh->volumeMesh.cellMediumParameters[cellIndex].invRho;
    IndexType submeshIndex = mesh->volumeMesh.cellMediumParameters[cellIndex].submeshIndex;

    typename MeshType::Elastic elastic = stateMaker->GetValue(globalPoint, lambda, mju, invRho, submeshIndex);
    for (IndexType valueIndex = 0; valueIndex < mesh->dimsCount; ++valueIndex)
    {
      values[valueIndex] = elastic.values[valueIndex];
    }
  }

  StateMaker* stateMaker;
  MeshType*   mesh;
  IndexType   cellIndex;
  Vector      cellVertices[Space::NodesPerCell];
};

template <typename MeshType, typename Corrector>
void ApplyCorrector(MeshType* const mesh)
{
  using Scalar = typename MeshType::Scalar;
  using IndexType = typename MeshType::IndexType;
  using Vector = typename MeshType::Vector;

  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(mesh->volumeMesh.cells.size()); ++cellIndex)
  {
    Corrector corrector(mesh, cellIndex);
    Scalar cellValues[mesh->functionsCount * mesh->dimsCount];
    
    mesh->volumeMesh.functionSpace->template DecomposePrecomputed< Corrector, MeshType::dimsCount >(corrector, cellValues);
    //mesh->volumeMesh.functionSpace->template Decompose< Corrector, MeshType::dimsCount >(corrector, cellValues);


    for (IndexType functionIndex = 0; functionIndex < mesh->functionsCount; ++functionIndex)
    {
      for (IndexType valueIndex = 0; valueIndex < mesh->dimsCount; ++valueIndex)
      {
        mesh->volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
          cellValues[valueIndex * mesh->functionsCount + functionIndex];
      }
    }
  }
}

/****************************************************
 *                                                  *
 *       ElasticVolumeMeshCommon.inl                *
 *                                                  *
 ****************************************************/

template<typename Space, typename FunctionSpace>
ElasticVolumeMeshCommon<Space, FunctionSpace>::
  ElasticVolumeMeshCommon(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  DifferentialSystem<Scalar>(solver->GetPhasesCount(), hierarchyLevelsCount), volumeMesh(solver->GetPhasesCount(), hierarchyLevelsCount), 
  allowMovement(false), 
  allowPlasticity(false), 
  allowContinuousDestruction(false), 
  allowDiscreteDestruction(false),
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
    for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
    {
      for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex) 
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

    for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
    {
      for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex) 
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
    Vector* nodePositions = reinterpret_cast<Vector*>(currCoords + volumeMesh.GetDimentionsCount(solverState));

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
    Vector* nodePositions = reinterpret_cast<Vector*>(currCoords + volumeMesh.GetMaxDimentionsCount());

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
    const Vector* nodePositions = reinterpret_cast<const Vector*>(newCoords + volumeMesh.GetDimentionsCount(solverState));

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
    const Vector* nodePositions = reinterpret_cast<const Vector*>(newCoords + volumeMesh.GetDimentionsCount(solverState));

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
    const Vector* nodePositions = reinterpret_cast<const Vector*>(oldCoords + volumeMesh.GetMaxDimentionsCount());
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
    Vector* nodeVelocities = reinterpret_cast<Vector*>(derivatives + volumeMesh.GetDimentionsCount(solverState));

    std::fill_n(nodeVelocities, volumeMesh.nodes.size(), Vector::zeroVector());

    for(int cellIndex = 0; cellIndex < int(volumeMesh.cells.size()); cellIndex++)
    {
      if (!volumeMesh.cellMediumParameters[cellIndex].fixed)
      {
        IndexType cellIndices[Space::NodesPerCell];
        volumeMesh.GetFixedCellIndices(cellIndex, cellIndices);

        for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
        {
          Vector v = InterpolateElasticRef(cellIndex, nodeNumber, VolumeMeshTypeCommon::CellNode).GetVelocity();
          // Vector v = InterpolateElastic(cellIndex, volumeMesh.nodes[cellIndices[nodeNumber]].pos).GetVelocity();

          nodeVelocities[cellIndices[nodeNumber]] += v;
        }
      }
    }

    for (IndexType nodeIndex = 0; nodeIndex < volumeMesh.nodes.size(); ++nodeIndex)
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

/*template<typename MeshType>
struct ShearNullifierFunctionGetter
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  ShearNullifierFunctionGetter(MeshType* mesh, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->cellIndex = cellIndex;
  }

  void operator()(IndexType basisPointIndex, Scalar* values) const
  {
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, basisPointIndex, MeshType::VolumeMeshTypeCommon::Basis);
    Vector  velocity = elastic.GetVelocity();
    elastic.SetVelocity(velocity * mesh->GetDamping());
    std::copy(elastic.values, elastic.values + mesh->dimsCount, values);
  }
private:
  MeshType*   mesh;
  IndexType   cellIndex;

  Vector      cellVertices[Space::NodesPerCell];
};*/

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCell(IndexType cellIndex, Scalar powderShearMult)
{
  if (!volumeMesh.cellMediumParameters[cellIndex].destroyed)
  {
    //destroying material
    volumeMesh.cellMediumParameters[cellIndex].destroyed = true;
    /*
      from Chelnokov PhD work...

      mju` = alpha * mju 
      K = lambda + 2/3 * mju = const -> lambda`
    */
    /*volumeMesh.cellMediumParameters[cellIndex].lambda -= 2 * (1 - powderShearMult) / 3 * volumeMesh.cellMediumParameters[cellIndex].mju;
    volumeMesh.cellMediumParameters[cellIndex].mju *= powderShearMult;

    //destroying elastic

    Scalar cellValues[functionsCount * dimsCount];
    auto shearNullifier = [&](IndexType basisPointIndex, Scalar* values)
    {
      Elastic elastic = this->InterpolateElasticRef(cellIndex, basisPointIndex, VolumeMeshTypeCommon::Basis);
      //Scalar pressure = elastic.GetPressure();
      //Tensor tension = Tensor(-pressure);
      //elastic.SetTension(tension);
      Scalar eigenValues[Space::Dimension];
      elastic.GetTension().GetEigenValues(eigenValues);
      Scalar averagePressure = 0;
      for (int coordIndex = 0; coordIndex < Space::Dimension; coordIndex++)
      {
        averagePressure += eigenValues[coordIndex];
      }
      averagePressure /= Scalar(Space::Dimension);
      elastic.SetTension(Tensor(averagePressure));
      std::copy(elastic.values, elastic.values + dimsCount, values);
    };

    volumeMesh.functionSpace->template DecomposePrecomputed< decltype(shearNullifier), dimsCount >
      (shearNullifier, cellValues);

    for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
    {
      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
          cellValues[valueIndex * functionsCount + functionIndex];
      }
    }*/
    volumeMesh.cellMediumParameters[cellIndex].plasticity.k0 *= powderShearMult;
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
    const Vector* node0Positions = reinterpret_cast<const Vector*>(coords0 + volumeMesh.GetDimentionsCount(solverState));
    const Vector* node1Positions = reinterpret_cast<const Vector*>(coords1 + volumeMesh.GetDimentionsCount(solverState));

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

      auto func = [&](Vector refPoint)
      {
        Scalar pointDeformRate = 0;
        const bool brittle = volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
        if (!brittle)
        {
          Elastic elastic = volumeMesh.GetRefCellSolution(cellIndex, refPoint);
          const Scalar pressure = elastic.GetPressure();
          const Scalar k = k0 + a * pressure;

          if (ProcessPlasticity(k, elastic, false))
          {
            pointDeformRate = GetDeformRate(cellIndex, refPoint);
          }
        }
        return pointDeformRate;
      };

      plasticDeformRate = volumeMesh.functionSpace->template GetCellAverage<decltype(func), Scalar>(func);

      plasticDeforms[cellIndex] += plasticDeformRate * dt; //* jacobian / volumeMesh.GetVolume(cellIndex) -- possibly needed, but i'm not sure

      if (plasticDeforms[cellIndex] > volumeMesh.cellMediumParameters[cellIndex].plasticity.maxPlasticDeform)
      {
        DestroyCell(cellIndex, volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
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

  AsymmetricTensor refVelocityGradient;

  for (IndexType velocityComponentIndex = 0; velocityComponentIndex < Space::Dimension; ++velocityComponentIndex)
  {
    refVelocityGradient[velocityComponentIndex] = Vector::zero();
    for (IndexType functionIndex = 0; functionIndex < FunctionSpace::functionsCount; ++functionIndex)
    {
      refVelocityGradient[velocityComponentIndex] += 
        basisFunctionGradients[functionIndex] * volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].GetVelocity()[velocityComponentIndex];
    }
  }

  AsymmetricTensor velocityGradient;
  for (IndexType velocityComponentIndex = 0; velocityComponentIndex < Space::Dimension; ++velocityComponentIndex)
  {
    velocityGradient[velocityComponentIndex] = Vector::zero();
    for (IndexType dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
    {
      velocityGradient[velocityComponentIndex] +=
        refDerivatives[dimIndex] * refVelocityGradient[velocityComponentIndex][dimIndex];
    }
  }

  AsymmetricTensor strainRateTensor;
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


template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::FindDestructions(std::vector<bool>* isCellBroken)
{
  if (allowContinuousDestruction)
  {
    HandleContinuousDestruction();
    for (IndexType cellIndex = 0; cellIndex < volumeMesh.cellMediumParameters.size(); ++cellIndex)
    {
      (*isCellBroken)[cellIndex] = volumeMesh.cellMediumParameters[cellIndex].destroyed;
    }
  }

  if (allowDiscreteDestruction)
  {
    #pragma omp parallel
    {
      int threadIndex = omp_get_thread_num();
      // todo
      IndexType stateIndex = threadIndex * volumeMesh.GetMaxHierarchyLevel() * volumeMesh.GetHierarchyLevelsCount() * 2 + 1;
      int segmentBegin = volumeMesh.threadSegmentBegins[stateIndex];
      int segmentEnd = volumeMesh.threadSegmentEnds[stateIndex];

      for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
      {
        Elastic elastic = volumeMesh.GetCellAverageSolution(cellIndex);

        Scalar principalStresses[Space::Dimension];
        elastic.GetTension().GetEigenValues(principalStresses);

        Vector principalNormals[Space::Dimension];
        elastic.GetTension().GetEigenValues(principalStresses, principalNormals);

        IndexType maxStressIndex = IndexType(-1);

        for (IndexType stressIndex = 0; stressIndex < Space::Dimension; ++stressIndex)
        {
          if (maxStressIndex == IndexType(-1) || principalStresses[stressIndex] > principalStresses[maxStressIndex])
          {
            maxStressIndex = stressIndex;
          }
        }

        bool normalsUsed[Space::FacesPerCell];

        Scalar maxCosOfAngle = Scalar(-1.0);
        IndexType bestFaceNumber = IndexType(-1);
        IndexType dynamicBoundaryTypes[Space::FacesPerCell];

        // if this cell is deleted due to erosion
        if (!volumeMesh.isCellAvailable[cellIndex] || volumeMesh.cellMediumParameters[cellIndex].fixed) continue;
        if (volumeMesh.cellMediumParameters[cellIndex].destroyed) continue;

        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
        {
          normalsUsed[faceNumber] = false;
          Scalar maxLongitudinalStress;
          Scalar maxShearStress;

          IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);
          IndexType correspondingFaceNumber = volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber);
          IndexType interactionType = volumeMesh.GetInteractionType(cellIndex, faceNumber);

          if (interactionType == IndexType(-1)) continue;

          dynamicBoundaryTypes[faceNumber] = volumeMesh.system.GetContactDynamicBoundaryType(interactionType);

          if (dynamicBoundaryTypes[faceNumber] != IndexType(-1))
          {
            volumeMesh.system.GetContactCriticalInfo(interactionType, maxShearStress, maxLongitudinalStress);
          }
          else
          {
            continue;
          }

          if (principalStresses[maxStressIndex] < fabs(maxLongitudinalStress))
          {
            continue;
          }

          const Vector& mainNormal = principalNormals[maxStressIndex];
          Vector faceNormal = volumeMesh.GetFaceExternalNormal(cellIndex, faceNumber).GetNorm();

          if (maxCosOfAngle < fabs(mainNormal * faceNormal))
          {
            maxCosOfAngle = fabs(mainNormal * faceNormal);
            bestFaceNumber = faceNumber;
          }
        }

        if (bestFaceNumber != IndexType(-1))
        {
          IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, bestFaceNumber);
          IndexType correspondingFaceNumber = volumeMesh.GetCorrespondingFaceNumber(cellIndex, bestFaceNumber);
          if (correspondingCellIndex == IndexType(-1) ||
            correspondingFaceNumber == IndexType(-1) ||
            volumeMesh.cellMediumParameters[correspondingCellIndex].destroyed) continue;

          // check if it is last non-destructed edge for this cell
          IndexType incidentCellsCount = 0;
          for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
          {
            if (volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber) != IndexType(-1) &&
              volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber) != IndexType(-1))
            {
              ++incidentCellsCount;
            }
          }

          if (incidentCellsCount == 1)
          {
            DestroyCell(cellIndex, volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
            (*isCellBroken)[cellIndex] = true;
            continue;
          }

          if (volumeMesh.isCellAvailable[cellIndex])
          {
            //  volumeMesh.AddToAABBTree(cellIndex);
          }

          if (volumeMesh.isCellAvailable[correspondingCellIndex])
          {
            //  volumeMesh.AddToAABBTree(correspondingCellIndex);
          }

          #pragma omp critical
          {
            DestroyFace(cellIndex, bestFaceNumber, dynamicBoundaryTypes[bestFaceNumber]);
          }
        }
      }
    }
  }

  HandleMaterialErosion();
}

template<typename Space, typename FunctionSpace>
typename Space::Vector ElasticVolumeMeshCommon<Space, FunctionSpace>::GetTotalImpulse() const
{
  Vector totalImpulse = Vector::zero();

  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    auto velocityFunc = [&](IndexType basisPointIndex)
    {
      return volumeMesh.GetRefCellSolution(cellIndex, basisPointIndex, VolumeMeshTypeCommon::Basis).GetVelocity();
    };

    Vector intVel = volumeMesh.functionSpace->template GetCellAveragePrecomputed<decltype(velocityFunc), Vector>(velocityFunc);

    totalImpulse += volumeMesh.GetVolume(cellIndex) * intVel / volumeMesh.cellMediumParameters[cellIndex].invRho;
  }
  return totalImpulse;
}

template<typename Space, typename FunctionSpace>
void ElasticVolumeMeshCommon<Space, FunctionSpace>::GetCellEnergy(IndexType cellIndex, Scalar& kineticEnergy, Scalar& potentialEnergy) const
{
  const Scalar lambda = volumeMesh.cellMediumParameters[cellIndex].lambda;
  const Scalar mu = volumeMesh.cellMediumParameters[cellIndex].mju;
  const Scalar rho = 1 / volumeMesh.cellMediumParameters[cellIndex].invRho;
  const Tensor identityTensor = Tensor(1);

  Scalar cellVolume = volumeMesh.GetVolume(cellIndex);

  auto kineticEnergyFunc = [&](IndexType basisPointIndex)
  {
    auto cellSolution = volumeMesh.GetRefCellSolution(cellIndex, basisPointIndex, VolumeMeshTypeCommon::Basis);
    return Scalar(0.5) * cellSolution.GetVelocity().SquareLen() * rho;
  };

  auto potentialEnergyFunc = [&](IndexType basisPointIndex)
  {
    auto cellSolution = volumeMesh.GetRefCellSolution(cellIndex, basisPointIndex, VolumeMeshTypeCommon::Basis);
    Tensor t = cellSolution.GetTension();
    Scalar s = lambda / (Space::Dimension * lambda + 2 * mu);

    // from Chelnokov phd thesis
    return Scalar(0.25) / mu * (DoubleConvolution(t, t) - s * Sqr(DoubleConvolution(t, identityTensor)));
  };

  kineticEnergy = volumeMesh.functionSpace->template GetCellAveragePrecomputed<decltype(kineticEnergyFunc), Scalar>(kineticEnergyFunc)   * cellVolume;
  potentialEnergy = volumeMesh.functionSpace->template GetCellAveragePrecomputed<decltype(potentialEnergyFunc), Scalar>(potentialEnergyFunc) * cellVolume;
}

template<typename Space, typename FunctionSpace>
typename Space::Scalar ElasticVolumeMeshCommon<Space, FunctionSpace>::GetTotalEnergy(Scalar& kineticEnergy, Scalar& potentialEnergy) const
{
  kineticEnergy = potentialEnergy = 0;

  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    Scalar kineticCellEnergy;
    Scalar potentialCellEnergy;
    GetCellEnergy(cellIndex, kineticCellEnergy, potentialCellEnergy);
    kineticEnergy += kineticCellEnergy;
    potentialEnergy += potentialCellEnergy;
  }
  Scalar totalEnergy = kineticEnergy + potentialEnergy;
  return totalEnergy;
}

template<typename Space, typename FunctionSpace>
typename Space::Scalar ElasticVolumeMeshCommon<Space, FunctionSpace>::GetTotalMass() const
{
  Scalar totalMass = 0;
  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    if (volumeMesh.IsCellInAABBTree(cellIndex))
    {
      totalMass += volumeMesh.GetVolume(cellIndex) / volumeMesh.cellMediumParameters[cellIndex].invRho;
    }
  }
  return totalMass;
}

/****************************************************
 *                                                  *
 *       ElasticVolumeMesh2.inl                     *
 *                                                  *
 ****************************************************/

template<typename FunctionSpace>
ElasticVolumeMesh<Space2, FunctionSpace>::ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  ElasticVolumeMeshCommon<Space2, FunctionSpace>(solver, tolerance, hierarchyLevelsCount)
{
}

template<typename FunctionSpace>
inline typename ElasticVolumeMesh<Space2, FunctionSpace>::Elastic 
ElasticVolumeMesh<Space2, FunctionSpace>::GetAverageEdgeElastic(IndexType cellIndex, IndexType edgeNumber) const
{
  return volumeMesh.GetEdgeAverageSolution(cellIndex, edgeNumber);
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::MakeSnapshot(Elastic* destData, 
  IndexType width, IndexType height, 
  Vector boxPoint1, Vector boxPoint2, bool halfStepSolution)
{
  for(IndexType i = 0; i < width * height; i++)
  {
    for(IndexType j = 0; j < volumeMesh.dimsCount; j++)
    {
      destData[i].values[j] = Scalar(0.0);
    }
  }

  for(IndexType i = 0; i < volumeMesh.cells.size(); i++)
  {
    Vector points[3];
    volumeMesh.GetCellVertices(i, points);

    Scalar 
      xmin = points[0].x, xmax = points[0].x,
      ymin = points[0].y, ymax = points[0].y;
    for(int j = 0; j < 3; j++)
    {
      if(points[j].x < xmin) xmin = points[j].x;
      if(points[j].x > xmax) xmax = points[j].x;
      if(points[j].y < ymin) ymin = points[j].y;
      if(points[j].y > ymax) ymax = points[j].y;
    }

    int ixmin = int(((xmin - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * Scalar(width  - 1));
    int ixmax = int(((xmax - boxPoint1.x) / (boxPoint2.x - boxPoint1.x)) * Scalar(width  - 1)) + 1;
    int iymin = int(((ymin - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * Scalar(height - 1));
    int iymax = int(((ymax - boxPoint1.y) / (boxPoint2.y - boxPoint1.y)) * Scalar(height - 1)) + 1;

    if(ixmin < 0) ixmin = 0;
    if(ixmin > int(width - 1)) ixmin = int(width - 1);
    if(ixmax < 0) ixmax = 0;
    if(ixmax > int(width - 1)) ixmax = int(width - 1);

    if(iymin < 0) iymin = 0;
    if(iymin > int(height - 1)) iymin = int(height - 1);
    if(iymax < 0) iymax = 0;
    if(iymax > int(height - 1)) iymax = int(height - 1);

    for(int x = ixmin; x <= ixmax; x++)
    {
      for(int y = iymin; y <= iymax; y++)
      {
        Vector point(boxPoint1.x + (boxPoint2.x - boxPoint1.x) * (Scalar(x) / Scalar(width - 1)),
                      boxPoint1.y + (boxPoint2.y - boxPoint1.y) * (Scalar(y) / Scalar(height - 1)));

        if (PointInCell<Scalar>(points, point))
        {
          Elastic e = InterpolateElastic(i, point, halfStepSolution);
          e.MakeDimension(tensionDimensionlessMult, velocityDimensionlessMult);
          destData[y * width + x] = e;
        }
      }
    }
  }
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::MakeSnapshot(Elastic* destData, 
  const Vector& origin, const Vector& spacing,
  const Vector& boxPoint1, const Vector& boxPoint2,
  bool halfStepSolution)
{
  IntAABB globalSnapshotArea;
  AABB snapshotArea;
  snapshotArea.Set(boxPoint1, boxPoint2);

  IndexVector snapshotSize = ::ComputeSnapshotAABBSize<Space2>(origin, spacing, snapshotArea, 
    &globalSnapshotArea);

  int globalIxmin = globalSnapshotArea.boxPoint1.x;
  int globalIymin = globalSnapshotArea.boxPoint1.y;

  for(IndexType elasticIndex = 0; elasticIndex < snapshotSize.GetVolume(); ++elasticIndex)
  {
    std::fill_n(destData[elasticIndex].values, volumeMesh.dimsCount, Scalar(0.0));
  }

  for(IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    Vector points[Space::NodesPerCell];
    volumeMesh.GetCellVertices(cellIndex, points);

    Scalar 
      xmin = points[0].x, xmax = points[0].x,
      ymin = points[0].y, ymax = points[0].y;

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if(points[nodeNumber].x < xmin) xmin = points[nodeNumber].x;
      if(points[nodeNumber].x > xmax) xmax = points[nodeNumber].x;
      if(points[nodeNumber].y < ymin) ymin = points[nodeNumber].y;
      if(points[nodeNumber].y > ymax) ymax = points[nodeNumber].y;
    }

    int ixmin = int((xmin - origin.x) / spacing.x) - 1;
    int ixmax = int((xmax - origin.x) / spacing.x) + 1;

    int iymin = int((ymin - origin.y) / spacing.y) - 1;
    int iymax = int((ymax - origin.y) / spacing.y) + 1;

    for(int x = ixmin; x <= ixmax; x++)
    {
      for(int y = iymin; y <= iymax; y++)
      {
        if(y - globalIymin >= 0 && y - globalIymin < int(snapshotSize.y) &&
           x - globalIxmin >= 0 && x - globalIxmin < int(snapshotSize.x))
        {
          Vector point(origin.x + spacing.x * Scalar(x),
                       origin.y + spacing.y * Scalar(y));

          if(PointInCell<Scalar>(points, point) && snapshotArea.Includes(point))
          {
            Elastic e = InterpolateElastic(cellIndex, point, halfStepSolution);
            e.MakeDimension(tensionDimensionlessMult, velocityDimensionlessMult);
            destData[(y - globalIymin) * snapshotSize.x + x - globalIxmin] = e;
          }
        }
      }
    }
  }
}

/*
  if a cell is destructed then only compress and shear tensions exist

  s:s < 2k^2 -- von Mises criterion;
  k = k0 - a * pressure, a > 0.
*/
template<typename MeshType>
struct ContinuousDestructionCorrector
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::Tensor    Tensor;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  ContinuousDestructionCorrector(MeshType* mesh, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->cellIndex = cellIndex;
  }

  void operator()(IndexType basisPointIndex, Scalar* values) const
  {
    bool brittle = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, basisPointIndex, MeshType::VolumeMeshTypeCommon::Basis);

    if (brittle)
    {
      const Scalar k0 = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.k0;
      const Scalar  a = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.a;
      const Scalar pressure = elastic.GetPressure();
      const Scalar k = k0 + a * pressure;

      if (mesh->ProcessPlasticity(k, elastic, false))
      {
        mesh->DestroyCell(cellIndex, mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
      }
    }

    if (mesh->volumeMesh.cellMediumParameters[cellIndex].destroyed)
    {
      Scalar mainStresses[2];
      Scalar maxTangentStress = Scalar(0.5) * (sqrt(Sqr(elastic.GetXX() - elastic.GetYY()) + 4 * Sqr(elastic.GetXY())));

      mainStresses[0] = Scalar(0.5) * (elastic.GetXX() + elastic.GetYY()) + maxTangentStress;
      mainStresses[1] = Scalar(0.5) * (elastic.GetXX() + elastic.GetYY()) - maxTangentStress;

      Vector normal = Vector(mainStresses[0] - elastic.GetYY(), elastic.GetXY()).GetNorm();

      for (IndexType stressIndex = 0; stressIndex < 2; ++stressIndex)
      {
        mainStresses[stressIndex] = std::min<Scalar>(0, mainStresses[stressIndex]);
      }

      // rotate from local main tensions axes to global coordinates
      elastic.SetTension(mainStresses[0] * Sqr(normal.x) + mainStresses[1] * Sqr(normal.y), 
        mainStresses[0] * Sqr(normal.y) + mainStresses[1] * Sqr(normal.x),
        (mainStresses[0] - mainStresses[1]) * normal.x * normal.y);
    }
    
    std::copy(elastic.values, elastic.values + mesh->volumeMesh.dimsCount, values);
  }

  MeshType* mesh;
  IndexType cellIndex;
};

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::HandleContinuousDestruction()
{
  typedef ContinuousDestructionCorrector< ElasticVolumeMesh<Space2, FunctionSpace> > DestructionCorrectorType;
  ApplyCorrector< ElasticVolumeMesh<Space2, FunctionSpace>, DestructionCorrectorType>(this);
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space2, FunctionSpace>::UnfoldMesh(Scalar minHeight, IndexType iterationsCount)
{
  bool found = true;
  for(IndexType iterationIndex = 0; (iterationIndex < iterationsCount) && found; iterationIndex++)
  {
    found = false;
    for(IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); cellIndex++)
    {
      IndexType cellIndices[3];
      volumeMesh.GetFixedCellIndices(cellIndex, cellIndices);

      for(IndexType nodeNumber = 0; nodeNumber < 3; nodeNumber++)
      {
        Vector &currVertex = volumeMesh.nodes [cellIndices[nodeNumber]].pos;
        Vector &edgeVertex0 = volumeMesh.nodes[cellIndices[(nodeNumber + 1) % 3]].pos;
        Vector &edgeVertex1 = volumeMesh.nodes[cellIndices[(nodeNumber + 2) % 3]].pos;
        // Vector oppositeEdgeVertices[2];

        Scalar side;

        Vector projectionPoint;
        if(ProjectPointAgainstLine(edgeVertex0, edgeVertex1, currVertex, projectionPoint, side) &&
           (projectionPoint - currVertex).SquareLen() < Sqr(minHeight))
        {
          found = true;
          Scalar len0 = (projectionPoint - edgeVertex0).Len();
          Scalar len1 = (projectionPoint - edgeVertex1).Len();
          Scalar ratio = len0 / (len0 + len1);

          Scalar pushDepth = minHeight - (projectionPoint - currVertex).Len();
          Vector pushVector = (projectionPoint - currVertex).GetNorm() * pushDepth;

          Scalar coef = 0.5;
          currVertex  -= pushVector * 0.5 * coef;
          edgeVertex0 += pushVector * 0.5 * (Scalar(1.0) - ratio) * coef;
          edgeVertex1 += pushVector * 0.5 * (              ratio) * coef;
        }
      }
    }
  }
}

/****************************************************
 *                                                  *
 *       ElasticVolumeMesh3.inl                     *
 *                                                  *
 ****************************************************/

template<typename FunctionSpace>
ElasticVolumeMesh<Space3, FunctionSpace>::ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  ElasticVolumeMeshCommon<Space3, FunctionSpace>(solver, tolerance, hierarchyLevelsCount)
{
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::UnfoldMesh(Scalar minHeight, IndexType iterationsCount)
{
}

template<typename MeshType>
struct ContinuousDestructionCorrector3
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::Tensor    Tensor;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  ContinuousDestructionCorrector3(MeshType* mesh, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->cellIndex = cellIndex;
  }

  void operator()(IndexType basisPointIndex, Scalar* values) const
  {
    bool brittle = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, basisPointIndex, MeshType::VolumeMeshType::Basis);

    if (brittle)
    {
      const Scalar k0 = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.k0;
      const Scalar  a = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.a;
      const Scalar pressure = elastic.GetPressure();
      const Scalar k = k0 + a * pressure;

      if (mesh->ProcessPlasticity(k, elastic, false))
      {
        mesh->DestroyCell(cellIndex, mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
      }
    }

    if (mesh->volumeMesh.cellMediumParameters[cellIndex].destroyed)
    {
      Scalar principalStresses[Space3::Dimension];

      Scalar test[Space3::Dimension];
      elastic.GetTension().GetEigenValues(test);

      bool needToCorrect = false;
      for (IndexType stressIndex = 0; stressIndex < Space3::Dimension; ++stressIndex)
      {
        if (principalStresses[stressIndex] > 0)
        {
          needToCorrect = true;
          break;
        }
      }

      if (needToCorrect)
      {
        Vector principalNormals[Space3::Dimension];
        elastic.GetTension().GetEigenValues(principalStresses, principalNormals);

        for (IndexType stressIndex = 0; stressIndex < Space3::Dimension; ++stressIndex)
        {
          principalStresses[stressIndex] = std::min<Scalar>(0, principalStresses[stressIndex]);
        }

        Tensor t;
        t.xx = principalStresses[0];
        t.yy = principalStresses[1];
        t.zz = principalStresses[2];

        std::swap(principalNormals[0].y, principalNormals[1].x);
        std::swap(principalNormals[0].z, principalNormals[2].x);
        std::swap(principalNormals[1].z, principalNormals[2].y);

        t = t.RotateAxes(principalNormals[0], principalNormals[1], principalNormals[2]);
        // rotate from local principal tensions axes to global coordinates
        elastic.SetTension(t);
      }
    }

    std::copy(elastic.values, elastic.values + mesh->volumeMesh.dimsCount, values);
  }

  MeshType* mesh;
  IndexType cellIndex;
};

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::HandleContinuousDestruction()
{
  typedef ContinuousDestructionCorrector3< ElasticVolumeMesh<Space3, FunctionSpace> > DestructionCorrectorType;
  ApplyCorrector< ElasticVolumeMesh<Space3, FunctionSpace>, DestructionCorrectorType>(this);
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::MakeSnapshot(Elastic* destData,
  const Vector& origin, const Vector& spacing,
  const Vector& boxPoint1, const Vector& boxPoint2,
  bool halfStepSolution)
{
  IntAABB globalSnapshotArea;
  AABB snapshotArea;
  snapshotArea.Set(boxPoint1, boxPoint2);

  IndexVector snapshotSize = ::ComputeSnapshotAABBSize<Space3>(origin, spacing, snapshotArea,
    &globalSnapshotArea);

  int globalIxmin = globalSnapshotArea.boxPoint1.x;
  int globalIymin = globalSnapshotArea.boxPoint1.y;
  int globalIzmin = globalSnapshotArea.boxPoint1.z;

  for (IndexType elasticIndex = 0; elasticIndex < snapshotSize.GetVolume(); ++elasticIndex)
  {
    std::fill_n(destData[elasticIndex].values, volumeMesh.dimsCount, Scalar(0.0));
  }

  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    Vector points[Space::NodesPerCell];
    volumeMesh.GetCellVertices(cellIndex, points);

    Scalar
      xmin = points[0].x, xmax = points[0].x,
      ymin = points[0].y, ymax = points[0].y,
      zmin = points[0].z, zmax = points[0].z;

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (points[nodeNumber].x < xmin) xmin = points[nodeNumber].x;
      if (points[nodeNumber].x > xmax) xmax = points[nodeNumber].x;
      if (points[nodeNumber].y < ymin) ymin = points[nodeNumber].y;
      if (points[nodeNumber].y > ymax) ymax = points[nodeNumber].y;
      if (points[nodeNumber].z < zmin) zmin = points[nodeNumber].z;
      if (points[nodeNumber].z > zmax) zmax = points[nodeNumber].z;
    }

    int ixmin = int((xmin - origin.x) / spacing.x) - 1;
    int ixmax = int((xmax - origin.x) / spacing.x) + 1;

    int iymin = int((ymin - origin.y) / spacing.y) - 1;
    int iymax = int((ymax - origin.y) / spacing.y) + 1;

    int izmin = int((zmin - origin.z) / spacing.z) - 1;
    int izmax = int((zmax - origin.z) / spacing.z) + 1;

    for (int x = ixmin; x <= ixmax; x++)
    {
      for (int y = iymin; y <= iymax; y++)
      {
        for (int z = izmin; z <= izmax; z++)
        {
          if (x - globalIxmin >= 0 && x - globalIxmin < int(snapshotSize.x) &&
              y - globalIymin >= 0 && y - globalIymin < int(snapshotSize.y) &&              
              z - globalIzmin >= 0 && z - globalIzmin < int(snapshotSize.z))
          {
            Vector point(origin.x + spacing.x * Scalar(x),
                         origin.y + spacing.y * Scalar(y),
                         origin.z + spacing.z * Scalar(z));

            if (PointInCell<Scalar>(points, point) && snapshotArea.Includes(point))
            {
              Elastic e = InterpolateElastic(cellIndex, point, halfStepSolution);
              e.MakeDimension(tensionDimensionlessMult, velocityDimensionlessMult);
              destData[(z - globalIzmin) * snapshotSize.x * snapshotSize.y +
                       (y - globalIymin) * snapshotSize.x + 
                       (x - globalIxmin)] = e;
            }
          }
        }
      }
    }
  }
}
