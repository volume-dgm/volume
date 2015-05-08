#pragma once

#include <limits>
#include <algorithm>
#include "../../ElasticSystem/ElasticSystem.h"
#include "../../ElasticSystem/IniStates.h"
#include "../../../Vtk/BasicVtkWriter.h"
#include "../../../Maths/Spaces.h"
#include "../../VolumeMethod/VolumeMesh/VolumeMesh.h"

template<typename Space, typename FunctionSpace>
struct ElasticVolumeMeshCommon: public DifferentialSystem<typename Space::Scalar>
{
  SPACE_TYPEDEFS

  typedef          ElasticSystem<Space>                                ElasticSystemType;
  typedef typename ElasticSystemType::ElasticSpace                     ElasticSpaceType;
  typedef typename ElasticSystemType::MediumParameters                 MediumParameters;
  typedef          VolumeMesh<Space, FunctionSpace, ElasticSystemType> VolumeMeshType;
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
  bool allowDestruction;

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
  Elastic InterpolateElastic(IndexType cellIndex, Vector globalPoint, bool halfStepSolution = false) const;

  Elastic GetAverageCellElastic(IndexType cellIndex) const;
  
  int GetDimentionsCount(const SolverState& solverState) const;
  int GetMaxDimentionsCount() const;

  void GetCurrDerivatives(Scalar* derivatives, const SolverState& solverState);
  void GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState& solverState);
  void GetCurrCoords(Scalar& time, Scalar* currCoords) const;

  void SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState);
  void SetCurrCoords(Scalar time, const Scalar* newCoords, const Scalar* oldCoords, const SolverState& solverState);
  void SetCurrCoords(Scalar time, const Scalar* oldCoords);
  
  void SetDetectors(const std::vector<Vector>& globalDetectorsPositions);
  void GetDetectorsData(IndexType detectorIndex, Scalar* data);

  Scalar GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState&, const Scalar* mults);

  virtual void HandlePlasticity() = 0;
  virtual void UnfoldMesh(Scalar minHeight, IndexType iterationsCount) = 0;

  void HandleDamping();
  void SetDamping(Scalar damping);
  Scalar GetDamping() const;

  void DestroyCellMaterial(IndexType cellIndex, Scalar alpha);
  ElasticSystemType* GetSystem();
  void MoveSceneToSnapshotRegion();

protected:
  // before any destructions
  std::vector< AdditionalCellInfo<Space> > initialAdditionalCellInfos;

  void ComputeElasticMults()
  {
    ComputeElasticMults(Overload<Space>());
  }
  // for error computation
  Scalar elasticMults[dimsCount];
  Scalar damping;

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
  typedef typename ElasticSystemType::ElasticSpace                         ElasticSpaceType;
  typedef typename ElasticSpaceType::Elastic                               Elastic;
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
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCellMaterial;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::initialAdditionalCellInfos;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;

  ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  virtual void MakeSnapshot(Elastic* destData,
    IndexType width, IndexType height, 
    Vector boxPoint1, Vector boxPoint2, bool halfStepSolution);

  virtual void MakeSnapshot(Elastic* destData,
    const Vector& origin, const Vector& spacing, 
    const Vector& boxPoint1, const Vector& boxPoint2,
    bool halfStepSolution);

  Elastic GetAverageEdgeElastic (IndexType cellIndex, IndexType edgeNumber) const;

  // if stress components on the edge are over threshold then contact will be substituted for free boundary
  void FindDestructions(EdgePairIndices* contactEdges, IndexType* contactEdgesCount, IndexType contactTypesCount,
                        std::vector<bool>* isContactBroken, std::vector<bool>* isCellBroken);

  void HandleContinuousDestruction();

  // stress correction due to plasticity in compliance with  Prandtl-Ruess model
  void HandlePlasticity();

  // replaces set of discrete separate cells to continuous destroyed material
  void MergeSeparateCells();
  void DisjoinMergedCells(EdgeLocationPair contactEdgesLocation, IndexType dynamicBoundaryType);
  bool IsProperContact(EdgeLocationPair contactEdgesLocation);

private:
  void UnfoldMesh(Scalar minHeight, IndexType iterationsCount);
};

template <typename FunctionSpace>
struct ElasticVolumeMesh<Space3, FunctionSpace>: public ElasticVolumeMeshCommon<Space3, FunctionSpace>
{
  SPACE3_TYPEDEFS
  typedef          Space3                                       Space;
  typedef          ElasticSystem<Space3>                        ElasticSystemType;
  typedef typename ElasticSystemType::ElasticSpace              ElasticSpaceType;
  typedef typename ElasticSystemType::MediumParameters          MediumParameters;

  typedef          VolumeMesh<Space3, FunctionSpace, ElasticSystemType>    VolumeMeshType;
  typedef typename VolumeMeshType::FacePairIndices                         FacePairIndices;
  typedef typename VolumeMeshType::BoundaryFace                            BoundaryFace;
  typedef typename VolumeMeshType::Node                                    Node;
  typedef typename VolumeMeshType::Cell                                    Cell;
  typedef typename VolumeMeshType::CellSolution                            CellSolution;
  typedef typename ElasticSpaceType::Elastic                               Elastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::InterpolateElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::tensionDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::velocityDimensionlessMult;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetAverageCellElastic;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::GetSystem;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::DestroyCellMaterial;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::initialAdditionalCellInfos;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using ElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;

  ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount);

  void HandlePlasticity();

  void UnfoldMesh(Scalar minHeight, IndexType iterationsCount);

  virtual void MakeSnapshot(Elastic* destData,
    const Vector& origin, const Vector& spacing, 
    const Vector& boxPoint1, const Vector& boxPoint2,
    bool halfStepSolution);

  void FindDestructions(FacePairIndices* contactEdges, IndexType* contactFacesCount, IndexType contactTypesCount,
    std::vector<bool>* isContactBroken, std::vector<bool>* isCellBroken);
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

    typename MeshType::Elastic elastic = stateMaker->GetValue(globalPoint, lambda, mju, invRho);
    for (IndexType valueIndex = 0; valueIndex < mesh->dimsCount; valueIndex++)
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
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Vector    Vector;

  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < (int)mesh->volumeMesh.cells.size(); ++cellIndex)
  {
    Corrector corrector(mesh, cellIndex);
    Scalar cellValues[mesh->functionsCount * mesh->dimsCount];
    
    mesh->volumeMesh.functionSpace.template Decompose< Corrector, MeshType::dimsCount >(corrector, cellValues);

    for (IndexType functionIndex = 0; functionIndex < mesh->functionsCount; functionIndex++)
    {
      for (IndexType valueIndex = 0; valueIndex < mesh->dimsCount; valueIndex++)
      {
        mesh->volumeMesh.cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
          cellValues[valueIndex * mesh->functionsCount + functionIndex];
      }
    }
  }
}

#include "ElasticVolumeMeshCommon.inl"
#include "ElasticVolumeMesh2.inl"
#include "ElasticVolumeMesh3.inl"
