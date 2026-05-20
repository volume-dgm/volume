#pragma once
#include "../../../Maths/MatrixMaths.h"

#include "../../../DifferentialSolvers/DifferentialSystem.h"
#include "../../../DifferentialSolvers/DifferentialSolver.h"
#include "../../GeomMesh/TimeHierarchyLevelsManager.h"

#include "../FunctionGetters/FunctionGetters.h"
#include "../../VolumeMethod/FunctionGetters/ContactFinder.h"
#include "../../GeomMesh/GeomMesh/GeomMesh.h"
#include <algorithm>
#include <mpi.h>
#include <omp.h>

#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "../../../Maths/QuadraturePrecomputer.h"
#include "../Cell.h"

// #define USE_SPARSE_MATRIX_FOR_DERIVATIVES

template <typename Space, typename FunctionSpace, typename System>
class VolumeMeshCommon: public DifferentialSystem<typename Space::Scalar>, public GeomMesh<Space>
{
public:
  SPACE_TYPEDEFS

  typedef System                                     SystemT;
  const static int dimsCount = System::dimsCount;
  const static int functionsCount = FunctionSpace::functionsCount;
  typedef typename System::MediumParameters          MediumParameters;
  typedef typename GeomMesh<Space>::Node             Node;
  typedef          GeomMesh<Space>                   GeomMeshT;

  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<int>    IntTypeCellInfo;
  typedef typename AdditionalCellInfo<Space>:: template AuxInfo<Scalar> ScalarTypeCellInfo;

  #ifdef USE_DYNAMIC_MATRICIES
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXDimFunc;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXFunc;
  #else
    using MatrixXDimFunc = Eigen::Matrix<Scalar, dimsCount, functionsCount>;
    using MatrixXFunc    = Eigen::Matrix<Scalar, functionsCount, functionsCount>;
  #endif

  VolumeMeshCommon(int solverPhasesCount, int hierarchyLevelsCount) :
    DifferentialSystem<Scalar>(solverPhasesCount, hierarchyLevelsCount), GeomMesh<Space>(),
    xDerivativeVolumeIntegralsSparse(functionsCount, functionsCount),
    yDerivativeVolumeIntegralsSparse(functionsCount, functionsCount),
    debugMode(false)
  {
    functionSpace = new FunctionSpace;

    // *2 because of we compute integrals from product of two functions of order N-th degree
    QuadraturePrecomputer::BuildQuadrature<typename Space::BorderSpace>(2 * FunctionSpace::order,
      quadratureWeightsForBorder, quadraturePointsForBorder);

    for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
    {
      std::vector<Vector> basisPoints = functionSpace->GetBasisPoints();
      for (IndexType pointIndex = 0; pointIndex < basisPoints.size(); ++pointIndex)
      {
        basisPointFunctionValues[functionIndex].push_back(functionSpace->GetBasisFunctionValue(basisPoints[pointIndex], functionIndex));
      }

      for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
      {
        cellNodeBasisFunctionValues[functionIndex].push_back(functionSpace->GetBasisFunctionValue(::Cell<Space>::GetNode(nodeNumber), functionIndex));
      }
    }

    allowDynamicCollisions = false;
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

    void SetToZero()
    {
      for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
      {
        std::fill(basisVectors[functionIndex].values, basisVectors[functionIndex].values + dimsCount, 0);
      }
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
  using GeomMesh<Space>::RemoveCellFromAABBTree;
  using GeomMesh<Space>::IsCellInAABBTree;
  using GeomMesh<Space>::treeNodeCellIndices;

  using GeomMesh<Space>::additionalCellInfos;
  using GeomMesh<Space>::GetMassCenter;
  using GeomMesh<Space>::GetMinHeight;
  using GeomMesh<Space>::GetCellVertices;
  using GeomMesh<Space>::GetFixedCellIndices;
  using GeomMesh<Space>::GetGhostCellVertices;
  using GeomMesh<Space>::GetAspectRatio;
  using GeomMesh<Space>::GetCorrespondingCellIndex;
  using GeomMesh<Space>::GetCorrespondingFaceNumber;
  using GeomMesh<Space>::GetInteractionType;
  using GeomMesh<Space>::GetCellFaceNodes;
  using GeomMesh<Space>::AddToAABBTree;
  using GeomMesh<Space>::GetCellAABB;
  using GeomMesh<Space>::GetVolume;

  enum PointType {Basis, CellNode};

  typename System::ValueType GetRefCellSolution(IndexType cellIndex, Vector refCoords, bool halfStepCellSolution = false) const;
  typename System::ValueType GetRefCellSolution(IndexType cellIndex, IndexType pointIndex, PointType pointType, bool halfStepCellSolution = false) const;
  typename System::ValueType GetCellSolution(IndexType cellIndex, Vector globalPoint, bool halfStepCellSolution = false) const;

  typename System::ValueType GetRefCellSolution(Scalar* coeffs, Vector refCoords) const;
  typename System::ValueType GetCellSolution(IndexType cellIndex, Scalar* coeffs, Vector globalPoint) const;

  typename System::MediumParameters GetRefCellParams(Scalar* coeffs, Vector refCoords) const;

  int  GetDimentionsCount(const SolverState&) const override;
  int  GetMaxDimentionsCount() const override;

  void GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState&) override;
  void GetCurrCoords(Scalar& time, Scalar* currCoords) const override;

  void SetCurrCoords(Scalar time, const Scalar* newCoords, const Scalar* oldCoords, const SolverState&) override;
  void SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState&) override;
  void SetCurrCoords(Scalar time, const Scalar* oldCoords) override;

  Scalar GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState&, const Scalar* mults) override;

  virtual Vector GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const = 0; // x -> ξ
  typename System::ValueType GetCellAverageSolution(IndexType cellIndex) const;

  // associatedPermutation = 0 1 2 -> 1 2 0
  void TransformCellSolution(IndexType cellIndex, IndexType associatedPermutation[Space::NodesPerCell], CellSolution* cellSolution);

  Scalar GetTimeStepPrediction() override;

  FunctionSpace* functionSpace;
  System system;

  Scalar time;
  Scalar collisionWidth;
  bool allowDynamicCollisions;

  // if cell isn`t available then all it derivatives is not calculated and merely is set to zeros. for example it is used for erosion. 
  std::vector<bool>               isCellAvailable;

  std::vector<CellSolution>       cellSolutions;
  std::vector<MediumParameters>   cellMediumParameters;

  TimeHierarchyLevelsManager<Space> timeHierarchyLevelsManager;
  void RebuildTimeHierarchyLevels(IndexType globalStepIndex, bool allowCollisions, bool externalInitialization = false);

  virtual Scalar GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const = 0;

  void BuildAABBTree(const Vector& boxPoint1, const Vector& boxPoint2);

  // for quadtature integration
  std::vector<Scalar> quadratureWeightsForBorder;
  std::vector<typename BorderSpace::Vector> quadraturePointsForBorder;

  struct CollisionsInfo
  {
    void Initialize(IndexType cellsCount)
    {
      CollisionNode emptyElem;
      collisionNodes.resize(cellsCount, emptyElem);
      poolSize = 0;
    }

    void Clear()
    {
      CollisionNode emptyElem;
      std::fill(collisionNodes.begin(), collisionNodes.end(), emptyElem);
      poolSize = 0;
    }

    struct CollisionNode
    {
      CollisionNode(): offset(0), count(0)
      {
      }
      IndexType offset;
      IndexType count;
    };

    std::vector<IndexType> pool;
    std::vector<CollisionNode> collisionNodes;
    IndexType poolSize;
  } collisionsInfo;

  void UpdateCollisionsInfo();

protected:
  
  std::vector<Scalar> cellMaxWaveSpeeds;

  std::vector<CellSolution>       halfStepCellSolutions;
  std::vector<CellSolution>       bufferCellSolutions;
  std::vector<char>               inBuffer;

  Eigen::Matrix<Scalar, functionsCount, functionsCount> cellVolumeIntegrals;
  Eigen::Matrix<Scalar, functionsCount, functionsCount> cellVolumeIntegralsInv;

  MatrixXFunc xDerivativeVolumeIntegrals;
  MatrixXFunc yDerivativeVolumeIntegrals;

  Eigen::SparseMatrix<Scalar> xDerivativeVolumeIntegralsSparse;
  Eigen::SparseMatrix<Scalar> yDerivativeVolumeIntegralsSparse;

  Scalar cellVolumeAverageIntegrals[functionsCount]; // for computing average values

  void Initialize();

  std::vector<IndexType> hierarchyDimentionsCount;
  std::vector<IndexType> threadCellsCount;
  std::vector<IndexType> threadWorkloads;
  std::vector<IndexType> threadCellOffsets;
  std::vector<IndexType> threadSegmentBegins;
  std::vector<IndexType> threadSegmentEnds;

  /* precomputed basis values for each basis decomposer`s point */
  std::vector<Scalar> basisPointFunctionValues[functionsCount];
  std::vector<Scalar> cellNodeBasisFunctionValues[functionsCount];

  /*
    There are regular cells, where we compute solution,
    and domain boundary cells which are taken from neighbouring domains and should not be computed here.
  */
  virtual bool IsCellRegular(IndexType cellIndex) const = 0;

  bool debugMode;

  // only this cells may collide with each other.
  bool IsReadyForCollisionCell(IndexType cellIndex) const;

private:
  // has non -1 dynamic boundary type
  bool IsCellHasDynamicBoundary(IndexType cellIndex) const;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/****************************************************
 *                                                  *
 *              VolumeMeshCommon.inl                *
 *                                                  *
 ****************************************************/

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellSolution(IndexType cellIndex, Vector refCoords, bool halfStepSolution) const
{ 
  typename System::ValueType result(Scalar(0.0));

  for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refCoords, functionIndex);
    for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        (halfStepSolution ? halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
          cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
      result.values[valueIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellSolution(IndexType cellIndex, IndexType pointIndex, PointType pointType, bool halfStepSolution) const
{
  typename System::ValueType result(Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = 0;

    switch (pointType)
    {
      case Basis: 
        basisFunctionValue = basisPointFunctionValues[functionIndex][pointIndex]; 
       break;
      case CellNode:
        basisFunctionValue = cellNodeBasisFunctionValues[functionIndex][pointIndex];
      break;
      default:
        assert(0);
    }

    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        (halfStepSolution ? halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
          cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
      result.values[valueIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetCellSolution(IndexType cellIndex, Vector globalPoint, bool halfStepSolution) const
{
  Vector points[Space::NodesPerCell];
  GetCellVertices(cellIndex, points);
  Vector refCoords = GlobalToRefVolumeCoords(globalPoint, points);
  return GetRefCellSolution(cellIndex, refCoords, halfStepSolution);
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellSolution(Scalar* coeffs, Vector refCoords) const
{
  typename System::ValueType result;
  result.SetZeroValues();

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refCoords, functionIndex);
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient = coeffs[valueIndex * functionsCount + functionIndex];
      result.values[valueIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::MediumParameters VolumeMeshCommon<Space, FunctionSpace, System>::
  GetRefCellParams(Scalar* coeffs, Vector refCoords) const
{
  typename System::MediumParameters result;
  std::fill(result.params, result.params + MediumParameters::ParamsCount, Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refCoords, functionIndex);
    for (IndexType paramIndex = 0; paramIndex < MediumParameters::ParamsCount; ++paramIndex)
    {
      Scalar basisFunctionCoefficient = coeffs[paramIndex * functionsCount + functionIndex];
      result.params[paramIndex] += basisFunctionCoefficient * basisFunctionValue;
    }
  }
  return result;
}

template <typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetCellSolution(IndexType cellIndex, Scalar* coeffs, Vector globalPoint) const
{
  Vector points[Space::NodesPerCell];
  GetCellVertices(cellIndex, points);
  Vector refCoords = GlobalToRefVolumeCoords(globalPoint, points);
  return GetRefCellSolution(coeffs, refCoords);
}

template <typename Space, typename FunctionSpace, typename System>
int VolumeMeshCommon<Space, FunctionSpace, System>::
  GetDimentionsCount(const SolverState& solverState) const
{
  return hierarchyDimentionsCount[solverState.Index()];
}

template <typename Space, typename FunctionSpace, typename System>
int VolumeMeshCommon<Space, FunctionSpace, System>::
  GetMaxDimentionsCount() const
{
  return dimsCount * functionsCount * cells.size();
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::GetCurrCoords(Scalar& time, Scalar* currCoords) const
{
  time = this->time;
  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
  {
    for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
      {
        currCoords[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] = 
          cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];
      }
    }
  } 
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::GetCurrCoords(Scalar& time, Scalar* currCoords, Scalar* oldCoords, const SolverState& solverState)
{
  time = this->time;
  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];
    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    bool auxCell;
    bool useHalfStepSolution;

    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            currCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] = 
              (useHalfStepSolution ?
               halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]:
               cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
          }
        }

        useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            oldCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] = 
              (useHalfStepSolution ?
               halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]:
               cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]);
          }
        }

        ++targetCellIndex;
      }
    }
  }

  const bool writeNeedToUpdate = debugMode;
  if (writeNeedToUpdate)
  {
    std::vector<int> needToUpdate(cells.size(), 0);
   
    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, true);
        needToUpdate[cellIndex] = useHalfStepSolution ? 1 : -1;
      }
    }


    MeshVtkWriter< Space, IntTypeCellInfo > meshWriter;

    std::ostringstream needToUpdateFileName;
    needToUpdateFileName << "out/GetData[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(needToUpdateFileName.str(), nodes, cells, needToUpdate);
  }
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::SetCurrCoords(Scalar time, const Scalar* oldCoords)
{
  this->time = time;

  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
  {
    for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
      {
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] = 
          oldCoords[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
      }
    }
  }
  std::fill(inBuffer.begin(), inBuffer.end(), 0);
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::SetCurrCoords(Scalar time, const Scalar* newCoords, const SolverState& solverState)
{
  const bool writeNeedToUpdate = debugMode;
  if (writeNeedToUpdate)
  {
    std::vector<int> needToUpdate(cells.size(), 0);
    std::vector<int> auxCells(cells.size(), 0);

    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        needToUpdate[cellIndex] = useHalfStepSolution ? 1 : -1;
        auxCells[cellIndex] = auxCell ? 1 : -1;
      }
    }

    MeshVtkWriter< Space, IntTypeCellInfo > meshWriter;

    std::ostringstream needToUpdateFileName;
    needToUpdateFileName << "out/NeedToUpdate[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(needToUpdateFileName.str(), nodes, cells, needToUpdate);

    std::ostringstream auxCellsFileName;
    auxCellsFileName << "out/AuxCell[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(auxCellsFileName.str(), nodes, cells, auxCells);
  }

  this->time = time;

  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;
    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            (useHalfStepSolution ? 
              halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
              cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]) = 
              newCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
          }
        }
        ++targetCellIndex;
      }
    }
  }
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::SetCurrCoords(Scalar time, 
  const Scalar* newCoords, const Scalar* oldCoords, const SolverState& solverState)
{
  const bool writeNeedToUpdate = debugMode;
  if (writeNeedToUpdate)
  {
    std::vector<int> needToUpdate(cells.size(), 0);
    std::vector<int> auxCells(cells.size(), 0);

    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);
        needToUpdate[cellIndex] = useHalfStepSolution ? 1 : -1;
        auxCells[cellIndex] = auxCell ? 1 : -1;
      }
    }

    MeshVtkWriter< Space, IntTypeCellInfo > meshWriter;

    std::ostringstream needToUpdateFileName;
    needToUpdateFileName << "out/NeedToUpdate[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(needToUpdateFileName.str(), nodes, cells, needToUpdate);

    std::ostringstream auxCellsFileName;
    auxCellsFileName << "out/AuxCell[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(auxCellsFileName.str(), nodes, cells, auxCells);
  }

  this->time = time;

  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;
    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        const bool useHalfStepSolution = timeHierarchyLevelsManager.UseHalfStepSolution(cellIndex, solverState, auxCell, false);

        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            if (auxCell || solverState.hierarchyPhase == 1)
            {
              (useHalfStepSolution ? 
                halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] :
                cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex]) =
                  oldCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
            } else 
            {
              // !aux && hierarchyPhase == 0
              halfStepCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
                newCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
            }

            if (!auxCell && GetHierarchyLevelsCount() > 1 && solverState.hierarchyPhase == 1)
            {
              bufferCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
                newCoords[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex];
              inBuffer[cellIndex] = 1;
            }
          }
        }
        ++targetCellIndex;
      }
    }
  }

  // load buffer cell solutions
  if (GetHierarchyLevelsCount() > 1 && solverState.IsPreInitial())
  {
    #pragma omp parallel for
    for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
    {
      if (inBuffer[cellIndex])
      {
        for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
        {
          for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
          {
            cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex] =
              bufferCellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];
          }
        }
        inBuffer[cellIndex] = 0;
      }
    }
  }
}

template <typename Space, typename FunctionSpace, typename System>
typename Space::Scalar VolumeMeshCommon<Space, FunctionSpace, System>::
  GetErrorValue(Scalar time, const Scalar* coords0, const Scalar* coords1, const SolverState& solverState, const Scalar* mults)
{
  Scalar maxError = Scalar(0.0);
  const bool writeErrors = debugMode;
  std::vector<Scalar> errors;
  if (writeErrors)
  {
    errors.resize(cells.size(), 0);
  }

  #pragma omp parallel
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd   = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    for(int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      bool auxCell;
      if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell))
      {
        if (!auxCell && IsCellRegular(cellIndex))
        {
          Scalar cellError = Scalar(0.0);
          for(IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
          {
            Scalar mult = (mults) ? mults[valueIndex] : Scalar(1.0);
            for(IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
            {
              cellError += mult * fabs(
                (coords0[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] -
                 coords1[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex])
                * sqrt(cellVolumeIntegrals(functionIndex, functionIndex)));
            }
          }
          if (writeErrors) errors[cellIndex] = cellError;
          #pragma omp flush(maxError)
          if (cellError > maxError)
          {
            #pragma omp critical
            {
              if(cellError > maxError) maxError = cellError;
            }
          }
        }
        ++targetCellIndex;
      }
    }
  }

  /*
    for (int cellIndex = 0; cellIndex < cells.size(); cellIndex++)
    {
      Scalar cellError = Scalar(0.0);
      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        Scalar mult = (mults) ? mults[valueIndex] : Scalar(1.0);
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          cellError += mult * fabs(
            (coords0[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex] -
            coords1[cellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex])
            * sqrt(cellVolumeIntegrals[functionsCount * functionIndex + functionIndex]));
        }
      }
      if (cellError > maxError) maxError = cellError;
    } */

  if (writeErrors)
  {
    MeshVtkWriter< Space, ScalarTypeCellInfo > meshWriter;

    std::ostringstream errorsFileName;
    errorsFileName << "out/Errors[" 
             << "globalStepIndex_" << solverState.globalStepIndex << " "
             << "hierarchyPhase_"  << solverState.hierarchyPhase  << " "
             << "hierarchyLevel_"  << solverState.hierarchyLevel  << " "
             << "phaseIndex_"      << solverState.phaseIndex
             << "].vtk";
    meshWriter.Write(errorsFileName.str(), nodes, cells, errors);
  }

  return maxError;
}

template <typename Space, typename FunctionSpace, typename System>
typename Space::Scalar VolumeMeshCommon<Space, FunctionSpace, System>::GetTimeStepPrediction()
{
  #pragma omp parallel for
  for (int cellIndex = 0; cellIndex < int(cells.size()); ++cellIndex)
  {
    cellMaxWaveSpeeds[cellIndex] = system.GetMaxWaveSpeed(cellMediumParameters[cellIndex]);
  }

  Scalar minTimeStep = std::numeric_limits<Scalar>::max();
  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
    {
      Scalar cellTimeStep = GetMinHeight(cellIndex) / cellMaxWaveSpeeds[cellIndex];
      minTimeStep = std::min(minTimeStep, cellTimeStep);
    }
  }
  return minTimeStep / Scalar(FunctionSpace::order + 0.5);
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::RebuildTimeHierarchyLevels(IndexType globalStepIndex, 
  bool allowCollisions, bool externalInitialization)
{
  // hierarchy phase count equals 2
  int threadsCount = 1;
  #pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
    {
      threadsCount = omp_get_num_threads();
    }
  }

  std::fill(threadWorkloads.begin(),     threadWorkloads.end(),     0);
  std::fill(threadCellsCount.begin(),    threadCellsCount.end(),    0);
  std::fill(threadCellOffsets.begin(),   threadCellOffsets.end(),   0);
  std::fill(threadSegmentBegins.begin(), threadSegmentBegins.end(), 0);
  std::fill(threadSegmentEnds.begin(),   threadSegmentEnds.end(),   0);

  // build time hierachy levels
  timeHierarchyLevelsManager.Initialize(cells.size(), GetHierarchyLevelsCount(), GetSolverPhasesCount(), externalInitialization);


  Scalar minTimeStep = GetTimeStepPrediction();

  if (globalStepIndex % 100 == 0)
    std::cout << "  Min timestep: " << minTimeStep << std::endl;

  if (GetHierarchyLevelsCount() > 1)
  {
    timeHierarchyLevelsManager.BuildTimeHierarchyLevels(this, cellMaxWaveSpeeds, allowCollisions, minTimeStep);

    const bool writeTimeHierarchyLevels = true;
    if (writeTimeHierarchyLevels)
    {
      timeHierarchyLevelsManager.SaveToVtk(this, globalStepIndex);
    }

    std::cout << "Computational cost of " << GetHierarchyLevelsCount() << "-level time hierarchy equals " << 
      timeHierarchyLevelsManager.GetComputationalTotalCost() << std::endl;
  }

  // precomputing
  SolverState solverState(GetHierarchyLevelsCount(), GetSolverPhasesCount());
  for (int globalStepIndex = 0; globalStepIndex < GetMaxHierarchyLevel(); ++globalStepIndex)
  {
    for (int hierarchyPhase = 0; hierarchyPhase < 2; ++hierarchyPhase)
    {
      for (int hierarchyLevel = 0; hierarchyLevel < GetHierarchyLevelsCount(); ++hierarchyLevel)
      {
        solverState.SetState(globalStepIndex, hierarchyPhase, hierarchyLevel);
        IndexType totalCellCount = 0;

        #pragma omp parallel for reduction(+:totalCellCount)
        for (int threadIndex = 0; threadIndex < threadsCount; ++threadIndex)
        {
          IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
          IndexType segmentBegin = threadSegmentBegins[stateIndex] = (threadIndex + 0) * cells.size() / threadsCount;
          IndexType segmentEnd   = threadSegmentEnds[stateIndex]   = (threadIndex + 1) * cells.size() / threadsCount;
          threadCellsCount[stateIndex] = 0;
          for (IndexType cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
          {
            // ***
            if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState))
            {
              ++threadCellsCount[stateIndex];
              if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
              {
                ++threadWorkloads[stateIndex];
              }
            }
          }
          totalCellCount += threadCellsCount[stateIndex];
        }
        hierarchyDimentionsCount[solverState.Index()] = dimsCount * functionsCount * totalCellCount;
        
        // cells per thread balancing
        bool balanced;
        do {
          balanced = false;
          for (int threadIndex = 0; threadIndex + 1 < threadsCount; ++threadIndex)
          {
            IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
            IndexType nextThreadStateIndex = (threadIndex + 1) * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
            while (threadWorkloads[stateIndex] > threadWorkloads[nextThreadStateIndex] + 1)
            {
              --threadSegmentEnds[stateIndex];
              --threadSegmentBegins[nextThreadStateIndex];

              IndexType cellIndex = threadSegmentEnds[stateIndex];
              if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState))
              {
                --threadCellsCount[stateIndex];
                ++threadCellsCount[nextThreadStateIndex];

                if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
                {
                  --threadWorkloads[stateIndex];
                  ++threadWorkloads[nextThreadStateIndex];
                  balanced = true;
                }
              }
            }
            while (threadWorkloads[stateIndex] + 1 < threadWorkloads[nextThreadStateIndex])
            {
              IndexType cellIndex = threadSegmentEnds[stateIndex];
              if (timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState))
              {
                ++threadCellsCount[stateIndex];
                --threadCellsCount[nextThreadStateIndex];

                if (!cellMediumParameters[cellIndex].fixed && isCellAvailable[cellIndex])
                {
                  ++threadWorkloads[stateIndex];
                  --threadWorkloads[nextThreadStateIndex];
                  balanced = true;
                }
              }
              ++threadSegmentEnds[stateIndex];
              ++threadSegmentBegins[nextThreadStateIndex];
            }
          }
        } while (balanced); 
      }
    }
  } 

  for (int globalStepIndex = 0; globalStepIndex < GetMaxHierarchyLevel(); ++globalStepIndex)
  {
    for (int hierarchyPhase = 0; hierarchyPhase < 2; ++hierarchyPhase)
    {
      for (int hierarchyLevel = 0; hierarchyLevel < GetHierarchyLevelsCount(); ++hierarchyLevel)
      {
        solverState.SetState(globalStepIndex, hierarchyPhase, hierarchyLevel);
        for (int threadIndex = 0; threadIndex + 1 < threadsCount; ++threadIndex)
        {
          IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
          IndexType nextThreadStateIndex = (threadIndex + 1) * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
          threadCellOffsets[nextThreadStateIndex] = threadCellOffsets[stateIndex] + threadCellsCount[stateIndex];
        }
      }
    }
  }
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::
  TransformCellSolution(IndexType cellIndex,
    IndexType associatedPermutation[Space::NodesPerCell],
    CellSolution* cellSolution)
{
  bool ordered = true;
  for (IndexType nodeNumber = 0; nodeNumber + 1 < Space::NodesPerCell; ++nodeNumber)
  {
    if (associatedPermutation[nodeNumber] > associatedPermutation[nodeNumber + 1]) ordered = false;
  }
  assert(ordered);
  // TODO for network communication
}

template <typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::Initialize()
{
  time = Scalar(0.0);
  collisionWidth = Scalar(0.0);
  debugMode = false;

  cellSolutions.resize(cells.size());
  if (GetMaxHierarchyLevel() > 1)
  {
    halfStepCellSolutions.resize(cells.size());
    bufferCellSolutions.resize(cells.size());
  }

  hierarchyDimentionsCount.resize(GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);

  cellMaxWaveSpeeds.resize(cells.size());

  IndexType threadsCount = 1;
  #pragma omp parallel
  {
    threadsCount = omp_get_num_threads();
  }

  isCellAvailable.resize(cells.size(), true);
  inBuffer.resize(cells.size(), 0);
  threadWorkloads.resize(threadsCount     * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadCellsCount.resize(threadsCount    * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadCellOffsets.resize(threadsCount   * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadSegmentBegins.resize(threadsCount * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
  threadSegmentEnds.resize(threadsCount   * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2);
}

template<typename Space, typename FunctionSpace, typename System>
typename System::ValueType VolumeMeshCommon<Space, FunctionSpace, System>::
  GetCellAverageSolution(IndexType cellIndex) const
{
  typename System::ValueType result;
  std::fill_n(result.values, dimsCount, Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];

      result.values[valueIndex] += basisFunctionCoefficient * cellVolumeAverageIntegrals[functionIndex] / ::Cell<Space>::GetRefCellVolume();
    }
  }
  return result;
}

template<typename Space, typename FunctionSpace, typename System>
bool VolumeMeshCommon<Space, FunctionSpace, System>::IsCellHasDynamicBoundary(IndexType cellIndex) const
{
  bool isCellHasDynamicBoundary = false;

  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
  {
    IndexType correspondingCellIndex = GetCorrespondingCellIndex(cellIndex, faceNumber);
    IndexType correspondingFaceNumber = GetCorrespondingFaceNumber(cellIndex, faceNumber);
    IndexType interactionType = GetInteractionType(cellIndex, faceNumber);

    if (correspondingCellIndex == IndexType(-1) && correspondingFaceNumber == IndexType(-1) && interactionType != IndexType(-1))
    {
      IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);
      if (dynamicContactType != IndexType(-1))
      {
        isCellHasDynamicBoundary = true;
        break;
      }
    }
  }
  return isCellHasDynamicBoundary;
}

template<typename Space, typename FunctionSpace, typename System>
bool VolumeMeshCommon<Space, FunctionSpace, System>::IsReadyForCollisionCell(IndexType cellIndex) const
{
  bool isReadyForCollisionCell = IsCellHasDynamicBoundary(cellIndex);
  
  for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
  {
    IndexType correspondingCellIndex = GetCorrespondingCellIndex(cellIndex, faceNumber);

    if (correspondingCellIndex != IndexType(-1) && IsCellHasDynamicBoundary(correspondingCellIndex))
    {
      isReadyForCollisionCell = true;
      break;
    }
  }
  return isReadyForCollisionCell;
}

template<typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::BuildAABBTree(const Vector& boxPoint1, const Vector& boxPoint2)
{
  printf("Building AABB tree for boundary cells\n");
  treeNodeCellIndices.resize(cells.size(), IndexType(-1));

  AABB dynamicContactBox(boxPoint1, boxPoint2);

  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    // if (IsReadyForCollisionCell(cellIndex))
    Vector cellVertices[Space::NodesPerCell];
    GetCellVertices(cellIndex, cellVertices);
    bool needToAdd = false;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (dynamicContactBox.Includes(cellVertices[nodeNumber]))
      {
        needToAdd = true;
        break;
      }
    }

    if (needToAdd)
    {
      treeNodeCellIndices[cellIndex] = aabbTree.InsertNode(GetCellAABB(cellIndex));
      aabbTree.SetUserData(treeNodeCellIndices[cellIndex], cellIndex);
    }
  }

  if (allowDynamicCollisions)
  {
    collisionsInfo.Initialize(cells.size());
  }

  UpdateCollisionsInfo();
}

template<typename Space, typename FunctionSpace, typename System>
void VolumeMeshCommon<Space, FunctionSpace, System>::UpdateCollisionsInfo()
{
  if (!allowDynamicCollisions) return;

  collisionsInfo.Clear();

  for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
  {
    if (!isCellAvailable[cellIndex] || cellMediumParameters[cellIndex].fixed) continue;

    bool isCellContacted = false;
    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
    {
      if (GetCorrespondingCellIndex(cellIndex, faceNumber) == IndexType(-1) &&
        GetCorrespondingFaceNumber(cellIndex, faceNumber) == IndexType(-1))
      {
        IndexType interactionType = GetInteractionType(cellIndex, faceNumber);
        if (interactionType != IndexType(-1))
        {
          IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);
          if (dynamicContactType != IndexType(-1))
          {
            isCellContacted = true;
            break;
          }
        }
      }
    }

    if (isCellContacted)
    {
      ContactFinder< VolumeMeshCommon<Space, FunctionSpace, System> > contactFinder(this, cellIndex);

      const IndexType MaxCollidedCellsCount = 999;
      IndexType collidedCellsCount = 0;
      IndexType collidedCellsPool[MaxCollidedCellsCount];
      contactFinder.Find(collidedCellsPool, &collidedCellsCount, MaxCollidedCellsCount); //second arguments has a pointer

      assert(collidedCellsCount <= MaxCollidedCellsCount);

      if (collidedCellsCount + collisionsInfo.poolSize > collisionsInfo.pool.size())
      {
        collisionsInfo.pool.resize(collidedCellsCount + collisionsInfo.poolSize);
      }

      collisionsInfo.collisionNodes[cellIndex].count = collidedCellsCount;
      collisionsInfo.collisionNodes[cellIndex].offset = collisionsInfo.poolSize;
      collisionsInfo.poolSize += collidedCellsCount;

      if (collidedCellsCount > 0)
      {
        std::copy(collidedCellsPool,
          collidedCellsPool + collidedCellsCount,
          collisionsInfo.pool.data() + collisionsInfo.collisionNodes[cellIndex].offset);
      }
    }
  }
}

/****************************************************
 *                                                  *
 *              Class VolumeMesh2                   *
 *                                                  *
 ****************************************************/

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

  typedef typename System::MatrixXDim MatrixXDim;
  typedef typename VolumeMeshCommon<Space, FunctionSpace, SystemT>::MatrixXDimFunc MatrixXDimFunc;
  typedef typename VolumeMeshCommon<Space, FunctionSpace, SystemT>::MatrixXFunc MatrixXFunc;

  using EdgeLocationPair = GeomMesh<Space2>::EdgeLocationPair;
  using EdgeLocation     = GeomMesh<Space2>::EdgeLocation;
  using EdgePairIndices  = GeomMesh<Space2>::EdgePairIndices;
  using BoundaryEdge     = GeomMesh<Space2>::BoundaryEdge;
  using EdgeIndices      = GeomMesh<Space2>::EdgeIndices;

  using VolumeMeshCommon<Space, FunctionSpace, System>::additionalCellInfos;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetMassCenter;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetCellVertices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetFixedCellIndices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetRefCellSolution;

  using VolumeMeshCommon<Space, FunctionSpace, System>::collisionsInfo;

  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetCellEdgeNodes;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetEdgeExternalNormal;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetCellEdgeMiddle;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetGhostCellVertices;

  using VolumeMeshCommon<Space, FunctionSpace, System>::GetHierarchyLevelsCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetMaxHierarchyLevel;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetSolverPhasesCount;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetCellDeformJacobian;

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
  using VolumeMeshCommon<Space, FunctionSpace, System>::xDerivativeVolumeIntegralsSparse;
  using VolumeMeshCommon<Space, FunctionSpace, System>::yDerivativeVolumeIntegralsSparse;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeAverageIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::Initialize;
  using VolumeMeshCommon<Space, FunctionSpace, System>::collisionWidth;

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
  using VolumeMeshCommon<Space, FunctionSpace, System>::quadratureWeightsForBorder;
  using VolumeMeshCommon<Space, FunctionSpace, System>::quadraturePointsForBorder;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetAspectRatio;
  using VolumeMeshCommon<Space, FunctionSpace, System>::isCellAvailable;
  using VolumeMeshCommon<Space, FunctionSpace, System>::IsReadyForCollisionCell;
  using VolumeMeshCommon<Space, FunctionSpace, System>::AddToAABBTree;

public:
  VolumeMesh(int solverPhasesCount, int hierarchyLevelsCount):
    VolumeMeshCommon<Space2, FunctionSpace, System>(solverPhasesCount, hierarchyLevelsCount)
  {}

  void GetCurrDerivatives(Scalar* derivatives, const SolverState&) override;

  void LoadGeom(Vector* vertexPositions, IndexType* cellIndices, IndexType verticesCount, IndexType cellsCount,
    EdgePairIndices*  contactEdges, IndexType* contactEdgesCount, IndexType contactTypesCount,
    BoundaryEdge*     boundaryEdges, IndexType* boundaryEdgesCount, IndexType boundaryTypesCount,
    MediumParameters* cellMediumParameters,
    IndexType *internalContactTypes);

  typename System::ValueType GetEdgeAverageSolution(IndexType cellIndex, IndexType edgeNumber) const;
  typename System::ValueType GetFaceAverageSolution(IndexType cellIndex, IndexType edgeNumber) const;

  Vector GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const override; // x -> ξ
  Vector RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const; // ξ -> x

  Scalar GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const override;
  Vector    GetRefXDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const; //(dξ/dx, dξ/dy) * J
  Vector    GetRefYDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const; //(dη/dx, dη/dy) * J

  void GetRefDerivatives(Vector cellVertices[Space::NodesPerCell], Vector* refDerivatives) const;

private:
  void BuildMatrices();
  bool IsCellRegular(IndexType cellIndex) const override;

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

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/****************************************************
 *                                                  *
 *              VolumeMesh2.inl                     *
 *                                                  *
 ****************************************************/
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
  #pragma omp parallel
  for(int functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
  {
    #pragma omp for
    for(int functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
    {
      printf(".");
      cellVolumeIntegrals(functionIndex0, functionIndex1) =
        functionSpace->ComputeCellVolumeIntegral       (functionIndex1, functionIndex0);

      Vector derivativeIntegral = functionSpace->ComputeDerivativeVolumeIntegral (functionIndex0, functionIndex1);
      xDerivativeVolumeIntegrals (functionIndex0, functionIndex1) = derivativeIntegral.x;
      yDerivativeVolumeIntegrals (functionIndex0, functionIndex1) = derivativeIntegral.y;

      for(IndexType srcEdgeNumber = 0; srcEdgeNumber < 3; srcEdgeNumber++)
      {
        outgoingFlux.srcEdges[srcEdgeNumber].surfaceIntegral(functionIndex0, functionIndex1) =
          functionSpace->ComputeOutgoingFlux(srcEdgeNumber, functionIndex0, functionIndex1);
        for(IndexType dstEdgeNumber = 0; dstEdgeNumber < 3; dstEdgeNumber++)
        {
          incomingFlux.srcEdges[srcEdgeNumber].dstEdges[dstEdgeNumber].surfaceIntegral(functionIndex0, functionIndex1) =
            functionSpace->ComputeIncomingFlux(srcEdgeNumber, dstEdgeNumber, functionIndex0, functionIndex1);
        }
      }
    }
  }

  cellVolumeIntegralsInv = cellVolumeIntegrals.inverse();

  xDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;
  yDerivativeVolumeIntegrals *= cellVolumeIntegralsInv;

  #ifdef USE_SPARSE_MATRIX_FOR_DERIVATIVES
    xDerivativeVolumeIntegralsSparse.reserve(Eigen::VectorXi::Constant(functionsCount, functionsCount));
    yDerivativeVolumeIntegralsSparse.reserve(Eigen::VectorXi::Constant(functionsCount, functionsCount));

    for (int functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
    {
      for (int functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
      {
        if (fabs(xDerivativeVolumeIntegrals(functionIndex0, functionIndex1)) > std::numeric_limits<Scalar>::epsilon())
          xDerivativeVolumeIntegralsSparse.insert(functionIndex0, functionIndex1) = xDerivativeVolumeIntegrals(functionIndex0, functionIndex1);
        if (fabs(yDerivativeVolumeIntegrals(functionIndex0, functionIndex1)) > std::numeric_limits<Scalar>::epsilon())
          yDerivativeVolumeIntegralsSparse.insert(functionIndex0, functionIndex1) = yDerivativeVolumeIntegrals(functionIndex0, functionIndex1);
      }
    }

    xDerivativeVolumeIntegralsSparse.makeCompressed();
    yDerivativeVolumeIntegralsSparse.makeCompressed();
  #endif

  for(IndexType srcEdgeNumber = 0; srcEdgeNumber < 3; srcEdgeNumber++)
  {
    outgoingFlux.srcEdges[srcEdgeNumber].surfaceIntegral *= cellVolumeIntegralsInv;

    for(IndexType dstEdgeNumber = 0; dstEdgeNumber < 3; dstEdgeNumber++)
    {
      incomingFlux.srcEdges[srcEdgeNumber].dstEdges[dstEdgeNumber].surfaceIntegral *= cellVolumeIntegralsInv;
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

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    cellVolumeAverageIntegrals[functionIndex] = functionSpace->ComputeCellVolumeIntegral(functionIndex);

    for(IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
    {
      edgeAverages[edgeNumber].surfaceIntegral[functionIndex] =
        functionSpace->ComputeEdgeFlux(edgeNumber, functionIndex);
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

    MatrixXDimFunc currCellValues;
    MatrixXDimFunc correspondingCellValues;
    MatrixXDimFunc timeDerivatives;
    MatrixXDimFunc edgeFlux;

    MatrixXDim edgeTransformMatrix;
    MatrixXDim edgeTransformMatrixInv;

    MatrixXDim xnAuxMatrix;
    MatrixXDim xInteriorMatrix;
    MatrixXDim xExteriorMatrix;

    MatrixXDim xMatrix;
    MatrixXDim yMatrix;

    MatrixXDim xMixedMatrix;
    MatrixXDim yMixedMatrix;

    Eigen::Matrix<Scalar, 1, dimsCount> boundaryMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> leftContactMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> rightContactMatrix;

    MatrixXDimFunc boundaryInfoValues;
    MatrixXDimFunc sourceValues;
    MatrixXDimFunc sourcePointValues;

    MatrixXDimFunc flux;

    Vector cellVertices[Space::NodesPerCell];
    Scalar tmp[dimsCount];

    for(int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      /*
        There are regular cells, where we compute solution,
        and domain boundary cells which are taken from neighbouring domains and should not be computed here.
      */

      bool auxCell;
      if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell)) continue;

      timeDerivatives.setZero();

      if (IsCellRegular(cellIndex) && isCellAvailable[cellIndex] && !cellMediumParameters[cellIndex].fixed)
      {
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

        system.BuildXMatrix(cellMediumParameters[cellIndex], xMatrix);
        system.BuildYMatrix(cellMediumParameters[cellIndex], yMatrix);

        for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; edgeNumber++)
        {
          edgeFlux.setZero();

          IndexType edgeNodeIndices[Space::NodesPerEdge];
          GetCellEdgeNodes(cellIncidentNodes, edgeNumber, edgeNodeIndices);

          Vector edgeGlobalVertices[Space::NodesPerEdge];
          for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerEdge; nodeNumber++)
          {
            edgeGlobalVertices[nodeNumber] = nodes[edgeNodeIndices[nodeNumber]].pos;
          }

          system.BuildEdgeTransformMatrix(edgeGlobalVertices, edgeTransformMatrix);
          system.BuildEdgeTransformMatrixInv(edgeGlobalVertices, edgeTransformMatrixInv);

          Scalar edgeLen = (edgeGlobalVertices[1] - edgeGlobalVertices[0]).Len(); //GetFaceSquare(faceGlobalVertices);
          Vector edgeNormal = GetEdgeExternalNormal(cellIndex, edgeNumber);

          IndexType correspondingCellIndex = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingCellIndex;
          IndexType correspondingEdgeNumber = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].correspondingEdgeNumber;
          IndexType interactionType = additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;

          if (interactionType == IndexType(-1)) continue;

          system.BuildXnAuxMatrix(cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            xnAuxMatrix);

          // interior matrix
          system.BuildXnInteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            edgeNormal, xnAuxMatrix, xInteriorMatrix);

          // exterior matrix
          system.BuildXnExteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex], 
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            edgeNormal, xnAuxMatrix, xExteriorMatrix);

          if (correspondingCellIndex == IndexType(-1))
          {
            // bounary condition
            IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);

            // external force/velocity 
            BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);
            if (functor)
            {
              typedef BoundaryFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > FunctorWrapper;
              FunctorWrapper wrapper(functor, time, this, cellIndex, edgeNumber);

              functionSpace->template Decompose< FunctorWrapper, dimsCount >(wrapper, boundaryInfoValues.data());

              edgeFlux.noalias() += Scalar(2.0) * (xExteriorMatrix * edgeTransformMatrixInv * boundaryInfoValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral);
            }

            if (!allowDynamicCollisions || dynamicContactType == IndexType(-1)) //set 1 for regular boundary, 0 for dynamic collisions
            {
              system.BuildBoundaryMatrix(interactionType, boundaryMatrix);
              edgeFlux.noalias() +=  (xInteriorMatrix + xExteriorMatrix * boundaryMatrix.asDiagonal()) * 
                edgeTransformMatrixInv * currCellValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral;
            }
            else
            {
              GhostCellFunctionGetter<VolumeMeshT> functionGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::Solution);

              GhostCellFunctionGetter<VolumeMeshT> paramsGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::MediumParams);

              flux.setZero();

              // quadrature integration of numerical flux
              for (IndexType pointIndex = 0; pointIndex < quadraturePointsForBorder.size(); ++pointIndex)
              {
                Vector globalPoint = (edgeGlobalVertices[1] - edgeGlobalVertices[0]) * quadraturePointsForBorder[pointIndex] + edgeGlobalVertices[0];

                Vector refPoint = GlobalToRefVolumeCoords(globalPoint, cellVertices);
                typename System::ValueType interiorSolution = GetRefCellSolution(cellIndex, refPoint);
                
                MatrixMulVector(edgeTransformMatrixInv.data(), interiorSolution.values, tmp, dimsCount, dimsCount);
                std::copy(tmp, tmp + dimsCount, interiorSolution.values);

                MediumParameters exteriorParams;
                std::fill(exteriorParams.params, exteriorParams.params + MediumParameters::ParamsCount, 0);
                IndexType collidedCellIndex = IndexType(-1);

                if (collisionWidth > std::numeric_limits<Scalar>::epsilon())
                  collidedCellIndex = paramsGetter(globalPoint, exteriorParams.params);
                else
                {
                  if (collisionsInfo.collisionNodes[cellIndex].count > 0)
                  {
                    collidedCellIndex = paramsGetter(globalPoint,
                      collisionsInfo.pool.data() + collisionsInfo.collisionNodes[cellIndex].offset,
                      collisionsInfo.collisionNodes[cellIndex].count, exteriorParams.params);
                  }
                }

                typename System::ValueType exteriorSolution;

                if (collidedCellIndex != IndexType(-1))
                {
                  if (!functionGetter.TryGhostCell(globalPoint + edgeNormal * collisionWidth, collidedCellIndex, exteriorSolution.values))
                  {
                    assert(0);
                  }

                  MatrixMulVector(edgeTransformMatrixInv.data(), exteriorSolution.values, tmp, dimsCount, dimsCount);
                  std::copy(tmp, tmp + dimsCount, exteriorSolution.values);
                }

                typename System::ValueType riemannSolution =
                  system.GetRiemannSolution(interiorSolution, exteriorSolution,
                    cellMediumParameters[cellIndex], exteriorParams,
                    interactionType, // boundary type
                    dynamicContactType);

                for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
                {
                  Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refPoint, functionIndex);
                  for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
                  {
                    flux(valueIndex, functionIndex) +=
                      quadratureWeightsForBorder[pointIndex] *
                      riemannSolution.values[valueIndex] * basisFunctionValue;
                  }
                }
              }

              edgeFlux.noalias() += xMatrix * flux * cellVolumeIntegralsInv;
            }
          }
          else
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
              edgeFlux.noalias() += xExteriorMatrix * leftContactMatrix.asDiagonal() * edgeTransformMatrixInv *
                currCellValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral;
            }

            // exterior side contribution
            edgeFlux.noalias() += xExteriorMatrix * rightContactMatrix.asDiagonal() * edgeTransformMatrixInv *
              correspondingCellValues * incomingFlux.srcEdges[edgeNumber].dstEdges[correspondingEdgeNumber].surfaceIntegral;

            // outgoing flux
            edgeFlux.noalias() += xInteriorMatrix * edgeTransformMatrixInv *
              currCellValues * outgoingFlux.srcEdges[edgeNumber].surfaceIntegral;
          }
          timeDerivatives.noalias() += (edgeLen * invJacobian) * (edgeTransformMatrix * edgeFlux);
        }

        Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
        Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;

        xMixedMatrix.noalias() = xMatrix * refXDerivatives.x + yMatrix * refXDerivatives.y;
        yMixedMatrix.noalias() = xMatrix * refYDerivatives.x + yMatrix * refYDerivatives.y;

        #ifdef USE_SPARSE_MATRIX_FOR_DERIVATIVES
        timeDerivatives.noalias() -=
          xMixedMatrix * currCellValues * xDerivativeVolumeIntegralsSparse +
          yMixedMatrix * currCellValues * yDerivativeVolumeIntegralsSparse;
        #else
        timeDerivatives.noalias() -=
          xMixedMatrix * currCellValues * xDerivativeVolumeIntegrals +
          yMixedMatrix * currCellValues * yDerivativeVolumeIntegrals;
        #endif

        typename SystemT::SourceFunctorT* sourceFunctor = system.GetSourceFunctor();
        if (sourceFunctor)
        {
          typedef SourceFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > SourceFunctorWrapper;
          SourceFunctorWrapper wrapper(sourceFunctor, time, this, cellVertices);
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
      }

      for(IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          derivatives[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex]
            = -timeDerivatives(valueIndex, functionIndex);
        }
      }
      ++targetCellIndex;
    }
  }
}

template<typename FunctionSpace, typename System>
typename Space2::Scalar VolumeMesh<Space2, FunctionSpace, System>::
  GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const
{
  return (cellVertices[1].x - cellVertices[0].x) * (cellVertices[2].y - cellVertices[0].y) -
         (cellVertices[2].x - cellVertices[0].x) * (cellVertices[1].y - cellVertices[0].y);
}

template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GetRefXDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dξ/dx, dξ/dy, dξ/dz) * J
{
    return Vector(
      cellVertices[2].y - cellVertices[0].y,
      cellVertices[0].x - cellVertices[2].x);
}

template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GetRefYDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dη/dx, dη/dy, dη/dz) * J
{
    return Vector(
      cellVertices[0].y - cellVertices[1].y,
      cellVertices[1].x - cellVertices[0].x);
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space2, FunctionSpace, System>::GetRefDerivatives(
  Vector cellVertices[Space::NodesPerCell],
  Vector* refDerivatives) const
{
  Scalar invJacobian = Scalar(1.0) / fabs(GetCellDeformJacobian(cellVertices));
  refDerivatives[0] = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
  refDerivatives[1] = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;
}


template<typename FunctionSpace, typename System>
typename Space2::Vector VolumeMesh<Space2, FunctionSpace, System>::
  GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const //x -> ?
{
  Scalar invJacobian = Scalar(1.0) / GetCellDeformJacobian(cellVertices);
  Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices);
  Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices);

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
  typename System::ValueType result(Scalar(0.0));

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
typename System::ValueType VolumeMesh<Space2, FunctionSpace, System>::
GetFaceAverageSolution(IndexType cellIndex, IndexType faceNumber) const
{
  return GetEdgeAverageSolution(cellIndex, faceNumber);
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

/****************************************************
 *                                                  *
 *              Class VolumeMesh3                   *
 *                                                  *
 ****************************************************/

template<typename FunctionSpace, typename System>
class VolumeMesh<Space3, FunctionSpace, System>: public VolumeMeshCommon<Space3, FunctionSpace, System>
{
public:
  SPACE3_TYPEDEFS
  typedef Space3 Space;
  typedef System SystemT;

  using FacePairIndices = GeomMesh<Space3>::FacePairIndices;
  using BoundaryFace = GeomMesh<Space3>::BoundaryFace;
  typedef VolumeMesh<Space3, FunctionSpace, SystemT> VolumeMeshT;
  typedef typename System::MediumParameters          MediumParameters;
  typedef typename VolumeMeshCommon<Space, FunctionSpace, SystemT>::GeomMeshT GeomMeshT;

  typedef typename System::MatrixXDim MatrixXDim;
  typedef typename VolumeMeshCommon<Space, FunctionSpace, SystemT>::MatrixXDimFunc MatrixXDimFunc;
  typedef typename VolumeMeshCommon<Space, FunctionSpace, SystemT>::MatrixXFunc MatrixXFunc;

  using VolumeMeshCommon<Space, FunctionSpace, System>::collisionsInfo;
  using VolumeMeshCommon<Space, FunctionSpace, System>::additionalCellInfos;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetMassCenter;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetCellVertices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetFixedCellIndices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetRefCellSolution;

  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetFaceExternalNormal;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetCellFaceNodes;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetFaceSquare;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetGhostCellVertices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GeomMeshT::GetCellFaceVertices;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetCellDeformJacobian;

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
  using VolumeMeshCommon<Space, FunctionSpace, System>::xDerivativeVolumeIntegralsSparse;
  using VolumeMeshCommon<Space, FunctionSpace, System>::yDerivativeVolumeIntegralsSparse;
  using VolumeMeshCommon<Space, FunctionSpace, System>::cellVolumeAverageIntegrals;
  using VolumeMeshCommon<Space, FunctionSpace, System>::Initialize;
  using VolumeMeshCommon<Space, FunctionSpace, System>::collisionWidth;

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
  using VolumeMeshCommon<Space, FunctionSpace, System>::quadratureWeightsForBorder;
  using VolumeMeshCommon<Space, FunctionSpace, System>::quadraturePointsForBorder;
  using VolumeMeshCommon<Space, FunctionSpace, System>::GetAspectRatio;
  using VolumeMeshCommon<Space, FunctionSpace, System>::isCellAvailable;
  using VolumeMeshCommon<Space, FunctionSpace, System>::IsReadyForCollisionCell;
  using VolumeMeshCommon<Space, FunctionSpace, System>::AddToAABBTree;

  VolumeMesh(int solverPhasesCount, int hierarchyLevelsCount):
    VolumeMeshCommon<Space3, FunctionSpace, System>(solverPhasesCount, hierarchyLevelsCount)
  {}

  void LoadGeom(Vector *vertexPositions, IndexType *cellIndices, IndexType verticesCount, IndexType cellsCount,
    FacePairIndices *contactFaces, IndexType *contactFacesCount, IndexType contactTypesCount,
    BoundaryFace    *boundaryFaces, IndexType *boundaryFacesCount, IndexType boundaryTypesCount,
    MediumParameters* cellMediumParameters,
    IndexType *internalContactTypes);

  void GetCurrDerivatives(Scalar* derivatives, const SolverState&) override;

  Vector GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const override; // x -> ξ
  Vector RefToGlobalVolumeCoords(Vector refCoords, Vector cellVertices[Space::NodesPerCell]) const; // ξ -> x
  typename System::ValueType GetFaceAverageSolution(IndexType cellIndex, IndexType faceNumber) const;

  Scalar GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const override;
  Vector    GetRefXDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const; //(dξ/dx, dξ/dy, dξ/dz) * J
  Vector    GetRefYDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const; //(dη/dx, dη/dy, dη/dz) * J
  Vector    GetRefZDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const; //(dζ/dx, dζ/dy, dζ/dz) * J
  void      GetRefDerivatives(Vector cellVertices[Space::NodesPerCell], Vector* refDerivatives) const;


private:
  void BuildMatrices();
  bool IsCellRegular(IndexType cellIndex) const override;

  MatrixXFunc zDerivativeVolumeIntegrals;

  struct OutgoingFlux
  {
    struct SrcFaceFlux
    {
      MatrixXFunc surfaceIntegral;
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
          MatrixXFunc surfaceIntegral;
        };
        OrientationFlux orientations[3];
      };
      DstFaceFlux dstFaces[Space::FacesPerCell];
    };
    SrcFaceFlux srcFaces[Space::FacesPerCell];
  } incomingFlux;

  struct FaceAverage
  {
    Scalar surfaceIntegral[VolumeMeshCommon<Space3, FunctionSpace, System>::functionsCount];
  } faceAverages[Space::FacesPerCell];

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/****************************************************
 *                                                  *
 *              VolumeMesh3.inl                     *
 *                                                  *
 ****************************************************/

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
  #ifdef USE_DYNAMIC_MATRICIES
  xDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  yDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);

  for (IndexType srcFaceNumber = 0; srcFaceNumber < 4; srcFaceNumber++)
  {
    outgoingFlux.srcFaces[srcFaceNumber].surfaceIntegral.resize(functionsCount, functionsCount);
    for (IndexType dstFaceNumber = 0; dstFaceNumber < 4; dstFaceNumber++)
    {
      for (IndexType orientationNumber = 0; orientationNumber < 3; orientationNumber++)
      {
        incomingFlux.srcFaces[srcFaceNumber].dstFaces[dstFaceNumber].orientations[orientationNumber].surfaceIntegral.resize(functionsCount, functionsCount);
      }
    }
  }

  xDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  yDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  zDerivativeVolumeIntegrals.resize(functionsCount, functionsCount);
  #endif

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


  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    cellVolumeAverageIntegrals[functionIndex] = functionSpace->ComputeCellVolumeIntegral(functionIndex);
    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
    {
      faceAverages[faceNumber].surfaceIntegral[functionIndex] =
        functionSpace->ComputeFaceFlux(faceNumber, functionIndex);
    }
  }
}

template<typename FunctionSpace, typename System>
void VolumeMesh<Space3, FunctionSpace, System>::
GetCurrDerivatives(Scalar *derivatives, const SolverState& solverState)
{
  Scalar res[32] = { 0 };

  #pragma omp parallel 
  {
    int threadIndex = omp_get_thread_num();
    IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
    int segmentBegin = threadSegmentBegins[stateIndex];
    int segmentEnd = threadSegmentEnds[stateIndex];

    IndexType offset = threadCellOffsets[stateIndex];
    IndexType targetCellIndex = offset + 0;

    MatrixXDimFunc currCellValues;
    MatrixXDimFunc correspondingCellValues;
    MatrixXDimFunc timeDerivatives;
    MatrixXDimFunc faceFlux;

    MatrixXDim faceTransformMatrix;
    MatrixXDim faceTransformMatrixInv;

    MatrixXDim xnAuxMatrix;
    MatrixXDim xInteriorMatrix;
    MatrixXDim xExteriorMatrix;

    Eigen::Matrix<Scalar, 1, dimsCount> boundaryMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> leftContactMatrix;
    Eigen::Matrix<Scalar, 1, dimsCount> rightContactMatrix;

    MatrixXDim xMatrix;
    MatrixXDim yMatrix;
    MatrixXDim zMatrix;

    MatrixXDim xMixedMatrix;
    MatrixXDim yMixedMatrix;
    MatrixXDim zMixedMatrix;

    MatrixXDimFunc boundaryInfoValues;
    MatrixXDimFunc sourceValues;
    MatrixXDimFunc sourcePointValues;

    Scalar tmp[dimsCount];
    MatrixXDimFunc flux;

    #ifdef USE_DYNAMIC_MATRICIES
    currCellValues.resize(dimsCount, functionsCount);
    correspondingCellValues.resize(dimsCount, functionsCount);
    timeDerivatives.resize(dimsCount, functionsCount);

    faceTransformMatrix.resize(dimsCount, dimsCount);
    faceTransformMatrixInv.resize(dimsCount, dimsCount);

    xnAuxMatrix.resize(dimsCount, dimsCount);
    xInteriorMatrix.resize(dimsCount, dimsCount);
    xExteriorMatrix.resize(dimsCount, dimsCount);

    xMatrix.resize(dimsCount, dimsCount);
    yMatrix.resize(dimsCount, dimsCount);
    zMatrix.resize(dimsCount, dimsCount);

    xMixedMatrix.resize(dimsCount, dimsCount);
    yMixedMatrix.resize(dimsCount, dimsCount);
    zMixedMatrix.resize(dimsCount, dimsCount);

    boundaryInfoValues.resize(dimsCount, functionsCount);
    ghostValues.resize(dimsCount, functionsCount);
    sourceValues.resize(dimsCount, functionsCount);
    sourcePointValues.resize(dimsCount, functionsCount);

    flux.resize(dimsCount, functionsCount);
    #endif

    Vector cellVertices[Space::NodesPerCell];

    double beginPhase = MPI_Wtime();
    for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
    {
      /*
        There are regular cells, where we compute solution,
        and domain boundary cells which are taken from neighbouring domains and should not be computed here.
      */
      bool auxCell;
      if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell)) continue;

      timeDerivatives.setZero();

      if (IsCellRegular(cellIndex) && isCellAvailable[cellIndex] && !cellMediumParameters[cellIndex].fixed)
      {
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

        system.BuildXMatrix(cellMediumParameters[cellIndex], xMatrix);
        system.BuildYMatrix(cellMediumParameters[cellIndex], yMatrix);
        system.BuildZMatrix(cellMediumParameters[cellIndex], zMatrix);

        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
        {
          Vector faceGlobalVertices[Space::NodesPerFace];
          GetCellFaceVertices(cellIndex, faceNumber, faceGlobalVertices);

          system.BuildFaceTransformMatrix(faceGlobalVertices, faceTransformMatrix);
          system.BuildFaceTransformMatrixInv(faceGlobalVertices, faceTransformMatrixInv);

          Scalar faceDeformJacobian = Scalar(2.0) * GetFaceSquare(faceGlobalVertices);

          Vector faceNormal = GetFaceExternalNormal(faceGlobalVertices);

          IndexType correspondingCellIndex = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
          IndexType correspondingFaceNumber = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingFaceNumber;
          IndexType correspondingFaceOrientation = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].orientation;
          IndexType interactionType = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;

          if (interactionType == IndexType(-1)) continue;

          faceFlux.setZero();

          system.BuildXnAuxMatrix(cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            xnAuxMatrix);

          // interior matrix
          system.BuildXnInteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex],
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xnAuxMatrix, xInteriorMatrix);

          // exterior matrix
          system.BuildXnExteriorMatrix(
            cellMediumParameters[cellIndex],
            //cellMediumParameters[cellIndex], 
            cellMediumParameters[(correspondingCellIndex == IndexType(-1)) ? cellIndex : correspondingCellIndex],
            faceNormal, xnAuxMatrix, xExteriorMatrix);

          if (correspondingCellIndex == IndexType(-1))
          {
            // bounary condition
            IndexType dynamicContactType = system.GetBoundaryDynamicContactType(interactionType);

            // external force/velocity 
            BoundaryInfoFunctor<Space>* functor = system.GetBoundaryInfoFunctor(interactionType);
            if (functor)
            {
              typedef BoundaryFunctionGetter< VolumeMesh<Space, FunctionSpace, System> > FunctorWrapper;
              FunctorWrapper wrapper(functor, time, this, cellIndex, faceNumber);

              // TODO: replace it with integration over face
              functionSpace->template Decompose< FunctorWrapper, dimsCount >(wrapper, boundaryInfoValues.data());
              faceFlux.noalias() += Scalar(2.0) * xExteriorMatrix * faceTransformMatrixInv * boundaryInfoValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
            }

            if (!allowDynamicCollisions || dynamicContactType == IndexType(-1))
            {
              system.BuildBoundaryMatrix(interactionType, boundaryMatrix);
              // regular boundary          
              faceFlux.noalias() += (xInteriorMatrix + xExteriorMatrix * boundaryMatrix.asDiagonal()) *
                faceTransformMatrixInv * currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
            }
            else
            {
              // dynamic collision
              GhostCellFunctionGetter<VolumeMeshT> functionGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::Solution);

              // functionSpace->template Decompose< GhostCellFunctionGetter<VolumeMeshT>, dimsCount >(functionGetter, ghostValues.data());

              GhostCellFunctionGetter<VolumeMeshT> paramsGetter(this, cellIndex, time,
                GhostCellFunctionGetter<VolumeMeshT>::MediumParams);

              flux.setZero();

              // quadrature integration of numerical flux
              for (IndexType pointIndex = 0; pointIndex < quadraturePointsForBorder.size(); ++pointIndex)
              {
                Vector globalPoint = (faceGlobalVertices[1] - faceGlobalVertices[0]) * quadraturePointsForBorder[pointIndex].x +
                  (faceGlobalVertices[2] - faceGlobalVertices[0]) * quadraturePointsForBorder[pointIndex].y +
                  faceGlobalVertices[0];

                Vector refPoint = GlobalToRefVolumeCoords(globalPoint, cellVertices);
                typename System::ValueType interiorSolution = GetRefCellSolution(cellIndex, refPoint);

                MatrixMulVector(faceTransformMatrixInv.data(), interiorSolution.values, tmp, dimsCount, dimsCount);
                std::copy(tmp, tmp + dimsCount, interiorSolution.values);

                MediumParameters exteriorParams;
                std::fill(exteriorParams.params, exteriorParams.params + MediumParameters::ParamsCount, 0);
                IndexType collidedCellIndex = IndexType(-1);

                if (collisionWidth > std::numeric_limits<Scalar>::epsilon())
                  collidedCellIndex = paramsGetter(globalPoint + faceNormal * collisionWidth, exteriorParams.params);
                else
                {
                  if (collisionsInfo.collisionNodes[cellIndex].count > 0)
                  {
                    collidedCellIndex = paramsGetter(globalPoint + faceNormal * collisionWidth,
                      collisionsInfo.pool.data() + collisionsInfo.collisionNodes[cellIndex].offset,
                      collisionsInfo.collisionNodes[cellIndex].count, exteriorParams.params);
                  }
                }

                typename System::ValueType exteriorSolution; //GetRefCellSolution(ghostValues.data(), ghostRefPoint);
                if (collidedCellIndex != IndexType(-1))
                {
                  if (!functionGetter.TryGhostCell(globalPoint + faceNormal * collisionWidth, collidedCellIndex, exteriorSolution.values))
                  {
                    assert(0);
                  }

                  MatrixMulVector(faceTransformMatrixInv.data(), exteriorSolution.values, tmp, dimsCount, dimsCount);
                  std::copy(tmp, tmp + dimsCount, exteriorSolution.values);
                }

                typename System::ValueType riemannSolution =
                  system.GetRiemannSolution(interiorSolution, exteriorSolution,
                    cellMediumParameters[cellIndex], exteriorParams,
                    interactionType, // boundary type
                    dynamicContactType);

                for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
                {
                  Scalar basisFunctionValue = functionSpace->GetBasisFunctionValue(refPoint, functionIndex);
                  for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
                  {
                    flux(valueIndex, functionIndex) +=
                      quadratureWeightsForBorder[pointIndex] *
                      riemannSolution.values[valueIndex] * basisFunctionValue;
                  }
                }
              }

              faceFlux.noalias() += xMatrix * flux * cellVolumeIntegralsInv;
            }
          }
          else
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

            /*
            // README: uncomment when glide internal contact is needed
            system.BuildContactMatrices(interactionType, leftContactMatrix, rightContactMatrix);
            // for glue contact left matrix equals 0
            if (!leftContactMatrix.isZero(std::numeric_limits<Scalar>::epsilon()))
            {
            // interior side contribution
            faceFlux.noalias() += xExteriorMatrix * leftContactMatrix.asDiagonal() *
            faceTransformMatrixInv * currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
            }
            */

            // exterior side contribution
            faceFlux.noalias() += xExteriorMatrix * // rightContactMatrix.asDiagonal() *
              faceTransformMatrixInv * correspondingCellValues *
              incomingFlux.srcFaces[faceNumber].dstFaces[correspondingFaceNumber].orientations[correspondingFaceOrientation].surfaceIntegral;

            // outgoing flux
            faceFlux.noalias() += xInteriorMatrix * faceTransformMatrixInv *
              currCellValues * outgoingFlux.srcFaces[faceNumber].surfaceIntegral;
          }

          timeDerivatives.noalias() += (faceDeformJacobian * invJacobian) * (faceTransformMatrix * faceFlux);
        }

        Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
        Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;
        Vector refZDerivatives = GetRefZDerivativesMulJacobian(cellVertices) * invJacobian;

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
          SourceFunctorWrapper wrapper(sourceFunctor, time, this, cellVertices);
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
      }

      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
        {
          derivatives[targetCellIndex * functionsCount * dimsCount + valueIndex * functionsCount + functionIndex]
            = -timeDerivatives(valueIndex, functionIndex);
        }
      }
      ++targetCellIndex;
    }
    double endPhase = MPI_Wtime();
    res[threadIndex] = endPhase - beginPhase;
  }

  /*
  Scalar minTime = res[0];
  Scalar maxTime = res[0];
  int minIndex = 0;
  const int threadsCount = 8;
  for (int i = 0; i < threadsCount; ++i)
  {
    if (res[i] < minTime)
    {
      minTime = res[i];
      minIndex = i;
    }
    maxTime = std::max(res[i], maxTime);

  }
  std::cout << "\n";

  if (solverState.globalStepIndex % 20 == 0)
  {
    std::cout << "VolumeMesh: GetCurrDerivatives " << (maxTime - minTime) / minTime << " " << minIndex << "\n";

    for (IndexType threadIndex = 0; threadIndex < threadsCount; ++threadIndex)
    {
      IndexType stateIndex = threadIndex * GetMaxHierarchyLevel() * GetHierarchyLevelsCount() * 2 + solverState.Index();
      int segmentBegin = threadSegmentBegins[stateIndex];
      int segmentEnd = threadSegmentEnds[stateIndex];

      int count = 0;

      for (int cellIndex = segmentBegin; cellIndex < segmentEnd; ++cellIndex)
      {
        bool auxCell;
        if (!timeHierarchyLevelsManager.NeedToUpdate(cellIndex, solverState, &auxCell)) continue;

        if (IsCellRegular(cellIndex) && isCellAvailable[cellIndex] && !cellMediumParameters[cellIndex].fixed)
        {
          //count += 6 + 4;
          count++;
        }
        continue;

        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; faceNumber++)
        {
          IndexType correspondingCellIndex = additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].correspondingCellIndex;
          if (correspondingCellIndex != IndexType(-1))
          {
            count += 6;
          }
          else
          {
            count += 3;
          }
        }

      }
      std::cout << count << " ";
    }
    std::cout << "\n";

    for (IndexType threadIndex = 0; threadIndex < 8; ++threadIndex)
    {
      std::cout << res[threadIndex] << " ";
    }
    std::cout << "\n";
  }
  */
}

template< typename FunctionSpace, typename System>
typename Space3::Scalar VolumeMesh<Space3, FunctionSpace, System>::
GetCellDeformJacobian(Vector cellVertices[Space::NodesPerCell]) const
{
  // deform jacobian equals 6 * tetrahedron volume
  return
    cellVertices[0].x * (cellVertices[1].y * (cellVertices[3].z - cellVertices[2].z) +
    cellVertices[2].y * (cellVertices[1].z - cellVertices[3].z) +
    cellVertices[3].y * (cellVertices[2].z - cellVertices[1].z)) +

    cellVertices[1].x * (cellVertices[0].y * (cellVertices[2].z - cellVertices[3].z) +
    cellVertices[2].y * (cellVertices[3].z - cellVertices[0].z) +
    cellVertices[3].y * (cellVertices[0].z - cellVertices[2].z)) +

    cellVertices[2].x * (cellVertices[0].y * (cellVertices[3].z - cellVertices[1].z) +
    cellVertices[1].y * (cellVertices[0].z - cellVertices[3].z) +
    cellVertices[3].y * (cellVertices[1].z - cellVertices[0].z)) +

    cellVertices[3].x * (cellVertices[0].y * (cellVertices[1].z - cellVertices[2].z) +
    cellVertices[1].y * (cellVertices[2].z - cellVertices[0].z) +
    cellVertices[2].y * (cellVertices[0].z - cellVertices[1].z));
}

template<typename FunctionSpace, typename System>
typename Space3::Vector VolumeMesh<Space3, FunctionSpace, System>::
  GlobalToRefVolumeCoords(Vector globalCoords, Vector cellVertices[Space::NodesPerCell]) const //x -> ?
{
  Scalar invJacobian = Scalar(1.0) / GetCellDeformJacobian(cellVertices);
  Vector refXDerivatives = GetRefXDerivativesMulJacobian(cellVertices);
  Vector refYDerivatives = GetRefYDerivativesMulJacobian(cellVertices);
  Vector refZDerivatives = GetRefZDerivativesMulJacobian(cellVertices);
  
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
  GetRefXDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dξ/dx, dξ/dy, dξ/dz) * J
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
  GetRefYDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dη/dx, dη/dy, dη/dz) * J
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
  GetRefZDerivativesMulJacobian(Vector cellVertices[Space::NodesPerCell]) const //(dζ/dx, dζ/dy, dζ/dz) * J
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
void  VolumeMesh<Space3, FunctionSpace, System>::GetRefDerivatives(
  Vector cellVertices[Space::NodesPerCell],
  Vector* refDerivatives) const
{
  Scalar invJacobian = Scalar(1.0) / fabs(GetCellDeformJacobian(cellVertices));
  refDerivatives[0] = GetRefXDerivativesMulJacobian(cellVertices) * invJacobian;
  refDerivatives[1] = GetRefYDerivativesMulJacobian(cellVertices) * invJacobian;
  refDerivatives[2] = GetRefZDerivativesMulJacobian(cellVertices) * invJacobian;
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

template<typename FunctionSpace, typename System>
typename System::ValueType VolumeMesh<Space3, FunctionSpace, System>::GetFaceAverageSolution(IndexType cellIndex, IndexType faceNumber) const
{
  typename System::ValueType result(Scalar(0.0));

  for (IndexType functionIndex = 0; functionIndex < functionsCount; ++functionIndex)
  {
    for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
    {
      Scalar basisFunctionCoefficient =
        cellSolutions[cellIndex].basisVectors[functionIndex].values[valueIndex];

      result.values[valueIndex] += basisFunctionCoefficient * faceAverages[faceNumber].surfaceIntegral[functionIndex] * Scalar(2.0);
    }
  }

  return result;
}
