#pragma once
#include "../../../Maths/MatrixMaths.h"

#include "../../../DifferentialSolvers/DifferentialSystem.h"
#include "../../../DifferentialSolvers/DifferentialSolver.h"
#include "../../GeomMesh/TimeHierarchyLevelsManager.h"

#include "../FunctionGetters/FunctionGetters.h"
#include "../../GeomMesh/GeomMesh/GeomMesh.h"
#include <algorithm>
#include <omp.h>

#include "Eigen/Dense"

template <typename Space, typename FunctionSpace, typename System>
class VolumeMeshCommon: public DifferentialSystem<typename Space::Scalar>, public GeomMesh<Space>
{
public:
  SPACE_TYPEDEFS

  const static int dimsCount = System::dimsCount;
  const static int functionsCount = FunctionSpace::functionsCount;
  typedef typename System::MediumParameters          MediumParameters;
  typedef typename GeomMesh<Space>::Node             Node;
  typedef          GeomMesh<Space>                   GeomMeshT;

  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<int>    IntTypeCellInfo;
  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<Scalar> ScalarTypeCellInfo;

  VolumeMeshCommon(int solverPhasesCount, int hierarchyLevelsCount):
    DifferentialSystem<Scalar>(solverPhasesCount, hierarchyLevelsCount), GeomMesh<Space>()
  {
    functionSpace = new FunctionSpace;
  }

  virtual ~VolumeMeshCommon()
  {
    delete functionSpace;
  }

  struct CellSolution
  {
  public:
    typename System::ValueType basisVectors[functionsCount];
    CellSolution(const CellSolution& other)
    {
      Copy(other);
    }

    CellSolution() {}

    CellSolution& operator=(const CellSolution& other)
    {
      Copy(other);
      return *this;
    }
  private:
    void Copy(const CellSolution& other)
    {
      for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
      {
        std::copy(other.basisVectors[functionIndex].values, 
          other.basisVectors[functionIndex].values + dimsCount, basisVectors[functionIndex].values);
      }
    }
  };

  using DifferentialSystem<Scalar>::GetHierarchyLevelsCount;
  using DifferentialSystem<Scalar>::GetMaxHierarchyLevel;
  using DifferentialSystem<Scalar>::GetSolverPhasesCount;
  using GeomMesh<Space>::UpdateAABBTree;

  using GeomMesh<Space>::nodes;
  using GeomMesh<Space>::cells;
  using GeomMesh<Space>::aabbTree;

  using GeomMesh<Space>::additionalCellInfos;
  using GeomMesh<Space>::GetMassCenter;
  using GeomMesh<Space>::GetCellVertices;
  using GeomMesh<Space>::GetFixedCellIndices;

  typename System::ValueType GetRefCellSolution(IndexType cellIndex, Vector refCoords, bool halfStepCellSolution = false) const;
  typename System::ValueType GetCellSolution(IndexType cellIndex, Vector globalPoint, bool halfStepCellSolution = false) const;

  typename System::ValueType GetRefCellSolution(Scalar* coeffs, Vector refCoords) const;
  typename System::ValueType GetCellSolution(Scalar* coeffs, Vector globalPoint) const;

  typename System::MediumParameters GetRefCellParams(Scalar* coeffs, Vector refCoords) const;

  int  GetDimentionsCount(const SolverState&) const;
  int  GetMaxDimentionsCount() const;

  void GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState&);
  void GetCurrCoords(Scalar& time, Scalar* currCoords) const;

  void SetCurrCoords(Scalar time, const Scalar* newCoords, const Scalar* oldCoords, const SolverState&);
  void SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState&);
  void SetCurrCoords(Scalar time, const Scalar* oldCoords);

  Scalar GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState&, const Scalar* mults);

  virtual Vector GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const = 0; // x -> ξ
  typename System::ValueType GetCellAverageSolution(IndexType cellIndex) const;

  // associatedPermutation = 0 1 2 -> 1 2 0
  void TransformCellSolution(IndexType cellIndex, IndexType associatedPermutation[Space::NodesPerCell], CellSolution* cellSolution);

  FunctionSpace* functionSpace;
  System system;

  Scalar time;
  Scalar collisionWidth;
  bool allowDynamicCollisions;

  std::vector<CellSolution>       cellSolutions;
  std::vector<MediumParameters>   cellMediumParameters;

  TimeHierarchyLevelsManager<Space> timeHierarchyLevelsManager;
  void RebuildTimeHierarchyLevels(IndexType globalStepIndex, bool allowCollisions, bool externalInitialization = false);

protected:
  std::vector<Scalar> cellMaxWaveSpeeds;

  std::vector<CellSolution>       halfStepCellSolutions;
  std::vector<CellSolution>       bufferCellSolutions;
  std::vector<char>               inBuffer;

  Eigen::Matrix<Scalar, functionsCount, functionsCount> cellVolumeIntegrals;
  Eigen::Matrix<Scalar, functionsCount, functionsCount> cellVolumeIntegralsInv;
  Eigen::Matrix<Scalar, functionsCount, functionsCount> xDerivativeVolumeIntegrals;
  Eigen::Matrix<Scalar, functionsCount, functionsCount> yDerivativeVolumeIntegrals;

  Scalar cellVolumeAverageIntegrals[functionsCount]; // for computing average values

  void Initialize();

  std::vector<IndexType> hierarchyDimentionsCount;
  std::vector<IndexType> threadCellsCount;
  std::vector<IndexType> threadCellOffsets;
  std::vector<IndexType> threadSegmentBegins;
  std::vector<IndexType> threadSegmentEnds;

  /*
    There are regular cells, where we compute solution,
    and domain boundary cells which are taken from neighbouring domains and should not be computed here.
  */
  virtual bool IsCellRegular(IndexType cellIndex) const = 0;

  bool debugMode;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#include "VolumeMeshCommon.inl"

template <typename Space, typename FunctionSpace, typename System>
class VolumeMesh;

template <typename FunctionSpace, typename System>
class VolumeMesh<Space2, FunctionSpace, System>: public VolumeMeshCommon<Space2, FunctionSpace, System>
{
public:
  SPACE2_TYPEDEFS
  typedef Space2                                     Space;
  typedef System                                     SystemT;
  typedef VolumeMesh<Space, FunctionSpace, SystemT>  VolumeMeshT;
  typedef typename System::MediumParameters          MediumParameters;
  typedef typename VolumeMeshCommon<Space, FunctionSpace, SystemT>::GeomMeshT GeomMeshT;

  typedef typename GeomMesh<Space2>::EdgeLocationPair EdgeLocationPair;
  typedef typename GeomMesh<Space2>::EdgeLocation     EdgeLocation;
  typedef typename GeomMesh<Space2>::EdgePairIndices  EdgePairIndices;
  typedef typename GeomMesh<Space2>::BoundaryEdge     BoundaryEdge;
  typedef typename GeomMesh<Space2>::EdgeIndices      EdgeIndices;

  using VolumeMeshCommon<Space, FunctionSpace, System>::additionalCellInfos;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetMassCenter;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetCellVertices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetFixedCellIndices;

  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetCellEdgeNodes;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetEdgeExternalNormal;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetCellEdgeMiddle;

  using VolumeMeshCommon<Space, FunctionSpace, System>::GetHierarchyLevelsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetMaxHierarchyLevel;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetSolverPhasesCount;

  using VolumeMeshCommon<Space, FunctionSpace, System>::nodes;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cells;
  using VolumeMeshCommon<Space, FunctionSpace, System>::aabbTree;

  using VolumeMeshCommon<Space, FunctionSpace, System>::dimsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::functionsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::functionSpace;
  using VolumeMeshCommon<Space, FunctionSpace, System>::system;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellMediumParameters;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeIntegralsInv;
  using VolumeMeshCommon<Space, FunctionSpace, System>::xDerivativeVolumeIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::yDerivativeVolumeIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeAverageIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::Initialize;

  using VolumeMeshCommon<Space, FunctionSpace, System>::hierarchyDimentionsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadCellsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadCellOffsets;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadSegmentBegins;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadSegmentEnds;

  using VolumeMeshCommon<Space, FunctionSpace, System>::timeHierarchyLevelsManager;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellSolutions;
  using VolumeMeshCommon<Space, FunctionSpace, System>::halfStepCellSolutions;
  using VolumeMeshCommon<Space, FunctionSpace, System>::allowDynamicCollisions;
  using VolumeMeshCommon<Space, FunctionSpace, System>::time;

public:
  VolumeMesh(int solverPhasesCount, int hierarchyLevelsCount):
    VolumeMeshCommon<Space2, FunctionSpace, System>(solverPhasesCount, hierarchyLevelsCount)
  {}

  void GetCurrDerivatives(Scalar* derivatives, const SolverState&);

  void LoadGeom(Vector* vertexPositions, IndexType* cellIndices, IndexType verticesCount, IndexType cellsCount,
    EdgePairIndices*  contactEdges, IndexType* contactEdgesCount, IndexType contactTypesCount,
    BoundaryEdge*     boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount,
    MediumParameters* cellMediumParameters,
    IndexType *internalContactTypes);

  typename System::ValueType GetEdgeAverageSolution(IndexType cellIndex, IndexType edgeNumber) const;

  Vector GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const; // x -> ξ
  Vector RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const; // ξ -> x

private:
  void BuildMatrices();

  Scalar GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const;
  Vector    GetRefXDerivatives(Vector cellVertices[Space::NodesPerCell]) const; //(dξ/dx, dξ/dy) * J
  Vector    GetRefYDerivatives(Vector cellVertices[Space::NodesPerCell]) const; //(dη/dx, dη/dy) * J

  bool IsCellRegular(IndexType cellIndex) const;

  struct OutgoingFlux
  {
    struct SrcEdgeFlux
    {
      Eigen::Matrix<Scalar,
        VolumeMeshCommon<Space2, FunctionSpace, System>::functionsCount, 
        VolumeMeshCommon<Space2, FunctionSpace, System>::functionsCount> surfaceIntegral;
    };
    SrcEdgeFlux srcEdges[Space::EdgesPerCell];
  } outgoingFlux;

  struct IncomingFlux
  {
    struct SrcEdgeFlux
    {
      struct DstEdgeFlux
      {
        Eigen::Matrix<Scalar,
          VolumeMeshCommon<Space2, FunctionSpace, System>::functionsCount,
          VolumeMeshCommon<Space2, FunctionSpace, System>::functionsCount> surfaceIntegral;
      };
      DstEdgeFlux dstEdges[Space::EdgesPerCell];
    };
    SrcEdgeFlux srcEdges[Space::EdgesPerCell];
  } incomingFlux;

  struct EdgeAverage
  {
    Scalar surfaceIntegral[VolumeMeshCommon<Space2, FunctionSpace, System>::functionsCount];
  } edgeAverages[Space::EdgesPerCell];

  // for quadtature integration
  void InitializeQuadrature(); 
  std::vector<Scalar> quadratureWeights;
  std::vector<Scalar> quadraturePoints;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include "VolumeMesh2.inl"

template<typename FunctionSpace, typename System>
class VolumeMesh<Space3, FunctionSpace, System>: public VolumeMeshCommon<Space3, FunctionSpace, System>
{
public:
  SPACE3_TYPEDEFS
  typedef Space3 Space;
  typedef System SystemT;

  typedef typename GeomMesh<Space3>::FacePairIndices FacePairIndices;
  typedef typename GeomMesh<Space3>::BoundaryFace    BoundaryFace;
  typedef VolumeMesh<Space3, FunctionSpace, SystemT> VolumeMeshT;
  typedef typename System::MediumParameters          MediumParameters;
  typedef typename VolumeMeshCommon<Space, FunctionSpace, SystemT>::GeomMeshT GeomMeshT;

  using VolumeMeshCommon<Space, FunctionSpace, System>::additionalCellInfos;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetMassCenter;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetCellVertices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetFixedCellIndices;

  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetFaceExternalNormal;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetCellFaceNodes;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetFaceSquare;

  using VolumeMeshCommon<Space, FunctionSpace, System>::nodes;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cells;
  using VolumeMeshCommon<Space, FunctionSpace, System>::aabbTree;

  using VolumeMeshCommon<Space, FunctionSpace, System>::dimsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::functionsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::functionSpace;
  using VolumeMeshCommon<Space, FunctionSpace, System>::system;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellMediumParameters;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeIntegralsInv;
  using VolumeMeshCommon<Space, FunctionSpace, System>::xDerivativeVolumeIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::yDerivativeVolumeIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeAverageIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::Initialize;

  using VolumeMeshCommon<Space, FunctionSpace, System>::GetHierarchyLevelsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetMaxHierarchyLevel;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetSolverPhasesCount;

  using VolumeMeshCommon<Space, FunctionSpace, System>::hierarchyDimentionsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadCellsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadCellOffsets;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadSegmentBegins;
  using VolumeMeshCommon<Space, FunctionSpace, System>::threadSegmentEnds;

  using VolumeMeshCommon<Space, FunctionSpace, System>::timeHierarchyLevelsManager;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellSolutions;
  using VolumeMeshCommon<Space, FunctionSpace, System>::halfStepCellSolutions;
  using VolumeMeshCommon<Space, FunctionSpace, System>::allowDynamicCollisions;
  using VolumeMeshCommon<Space, FunctionSpace, System>::time;


  VolumeMesh(int solverPhasesCount, int hierarchyLevelsCount):
    VolumeMeshCommon<Space3, FunctionSpace, System>(solverPhasesCount, hierarchyLevelsCount)
  {}

  void LoadGeom(Vector *vertexPositions, IndexType *cellIndices, IndexType verticesCount, IndexType cellsCount,
    FacePairIndices *contactFaces, IndexType *contactFacesCount, IndexType contactTypesCount,
    BoundaryFace    *boundaryFaces, IndexType *boundaryFacesCount, IndexType boundaryTypesCount,
    MediumParameters* cellMediumParameters,
    IndexType *internalContactTypes);

  void GetCurrDerivatives(Scalar* derivatives, const SolverState&);

  Vector GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const; // x -> ξ
  Vector RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const; // ξ -> x
  typename System::ValueType GetFaceAverageSolution(IndexType cellIndex, IndexType faceNumber) const;

private:
  void BuildMatrices();
  bool IsCellRegular(IndexType cellIndex) const;

  Scalar GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const;
  Vector    GetRefXDerivatives(Vector cellVertices[Space::NodesPerCell]) const; //(dξ/dx, dξ/dy, dξ/dz) * J
  Vector    GetRefYDerivatives(Vector cellVertices[Space::NodesPerCell]) const; //(dη/dx, dη/dy, dη/dz) * J
  Vector    GetRefZDerivatives(Vector cellVertices[Space::NodesPerCell]) const; //(dζ/dx, dζ/dy, dζ/dz) * J

  Eigen::Matrix<Scalar, functionsCount, functionsCount> zDerivativeVolumeIntegrals;

  struct OutgoingFlux
  {
    struct SrcFaceFlux
    {
      Eigen::Matrix<Scalar,
        VolumeMeshCommon<Space3, FunctionSpace, System>::functionsCount,
        VolumeMeshCommon<Space3, FunctionSpace, System>::functionsCount> surfaceIntegral;
    };
    SrcFaceFlux srcFaces[Space::FacesPerCell];
  } outgoingFlux;

  struct IncomingFlux
  {
    struct SrcFaceFlux
    {
      struct DstFaceFlux
      {
        struct OrientationFlux
        {
          Eigen::Matrix<Scalar,
            VolumeMeshCommon<Space3, FunctionSpace, System>::functionsCount,
            VolumeMeshCommon<Space3, FunctionSpace, System>::functionsCount> surfaceIntegral;
        };
        OrientationFlux orientations[3];
      };
      DstFaceFlux dstFaces[Space::FacesPerCell];
    };
    SrcFaceFlux srcFaces[Space::FacesPerCell];
  } incomingFlux;

  // for quadtature integration
  void InitializeQuadrature();
  std::vector<Scalar> quadratureWeights;
  std::vector<Vector2> quadraturePoints;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include "VolumeMesh3.inl"
