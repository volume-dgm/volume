#pragma once

template <typename GhostCellFunctionGetter>
struct CollisionProcessor
{
  typedef typename GhostCellFunctionGetter::MeshTypeT MeshType;
  typedef typename MeshType::Space   Space;

  typedef typename MeshType::SystemT     System;
  typedef typename System::ValueType     ValueType;
  SPACE_TYPEDEFS

  const static int dimsCount = MeshType::dimsCount;

  CollisionProcessor(GhostCellFunctionGetter* functionGetter,
    const Vector& reflectedPoint, const ValueType& innerSolution, Scalar* values):
    functionGetter(functionGetter), reflectedPoint(reflectedPoint), innerSolution(innerSolution), values(values), found(false)
  {}
  void operator()(int nodeIndex)
  {
    if (!found)
    {
      IndexType neighbourCellIndex = functionGetter->mesh->aabbTree.GetUserData(nodeIndex);
      if (functionGetter->TryGhostCell(reflectedPoint, innerSolution, neighbourCellIndex, values)) found = true;
    }
  }

  GhostCellFunctionGetter* functionGetter;
  Vector reflectedPoint;
  ValueType innerSolution;
  Scalar* values;
  bool found;
};

template<typename MeshType>
struct GhostCellFunctionGetter
{
  typedef typename MeshType::Space   Space;

  typedef typename MeshType::SystemT     System;
  typedef MeshType MeshTypeT;
  typedef typename System::ValueType     ValueType;
  SPACE_TYPEDEFS

  typedef CollisionProcessor< GhostCellFunctionGetter<MeshType> > CollisionProcessorT;

  const static int dimsCount = MeshType::dimsCount;

  GhostCellFunctionGetter(
    MeshType *mesh,
    IndexType cellIndex,
    Vector edgePoint, Vector edgeNormal, // face for 3d
    Scalar *edgeTransformMatrixInv,
    Scalar *boundaryEdgeMatrix,
    BoundaryInfoFunctor<Space>* functor,
    Scalar *leftContactEdgeMatrix,
    Scalar *rightContactEdgeMatrix,
    Scalar currTime,
    IndexType dynamicInteractionType):
    edgeNormal(edgeNormal)
  {
    this->mesh = mesh;
    this->cellIndex = cellIndex;
    this->edgePoint = edgePoint;
    this->edgeTransformMatrixInv = edgeTransformMatrixInv;
    this->boundaryEdgeMatrix = boundaryEdgeMatrix;
    this->functor = functor;
    this->leftContactEdgeMatrix = leftContactEdgeMatrix;
    this->rightContactEdgeMatrix = rightContactEdgeMatrix;
    this->currTime = currTime;
    this->dynamicInteractionType = dynamicInteractionType;

    mesh->GetCellVertices(cellIndex, this->cellVertices);
    mesh->GetFixedCellIndices(cellIndex, this->cellIndices);
  }

  void operator()(const Vector& point, Scalar *values)
  {
    Vector globalPoint = mesh->RefToGlobalVolumeCoords(point, cellVertices);
    Vector reflectedPoint = globalPoint - edgeNormal * (edgeNormal * (globalPoint - edgePoint)) * Scalar(2.0);

    // only points located on corresponding faces make contribution for incoming flux 
    if ((globalPoint - reflectedPoint).Len() > mesh->collisionWidth * Scalar(1e-3)) 
    {
      std::fill(values, values + dimsCount, 0);
      return;
    }

    reflectedPoint += edgeNormal * mesh->collisionWidth;
    ValueType innerSolution = mesh->GetCellSolution(cellIndex, globalPoint);

    CollisionProcessorT collisionProcessor(this, reflectedPoint, innerSolution, values);

    mesh->aabbTree.template FindCollisions<CollisionProcessorT>(AABB(reflectedPoint, reflectedPoint), collisionProcessor);

    if (collisionProcessor.found) return;

    // point is outside of all cells. counting as a ghost point
    MatrixMulVector(boundaryEdgeMatrix, innerSolution.values, values, dimsCount, dimsCount);

    if (functor)
    {
      ValueType externalInfo;
      functor->operator()(globalPoint, edgeNormal, /*innerSolution,*/ currTime, externalInfo.values);
      Scalar externalValues[dimsCount];
      MatrixMulVector(edgeTransformMatrixInv, externalInfo.values, externalValues, dimsCount, dimsCount);

      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        values[valueIndex] += externalValues[valueIndex];
      }
    }
  }

  Vector cellVertices[Space::NodesPerCell];
  IndexType cellIndices[Space::NodesPerCell];
  IndexType cellIndex;

  MeshType *mesh;
  Vector edgePoint;
  const Vector edgeNormal;
  Scalar *edgeTransformMatrixInv;
  Scalar *boundaryEdgeMatrix;
  BoundaryInfoFunctor<Space>* functor;
  Scalar *leftContactEdgeMatrix;
  Scalar *rightContactEdgeMatrix;
  Scalar currTime;
  IndexType dynamicInteractionType;

  bool TryGhostCell(const Vector& reflectedPoint, const ValueType& innerSolution, int neighbourCellIndex,
    Scalar* values)
  {
    IndexType neighbourCellIndices[Space::NodesPerCell];
    mesh->GetFixedCellIndices(neighbourCellIndex, neighbourCellIndices);

    bool hasSharedVertex = false;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      for (IndexType neighbourNodeNumber = 0; neighbourNodeNumber < Space::NodesPerCell; ++neighbourNodeNumber)
      {
        if (cellIndices[nodeNumber] == neighbourCellIndices[neighbourNodeNumber]) hasSharedVertex = true;
      }
    }
    if (hasSharedVertex) return false;

    Vector neighbourCellVertices[Space::NodesPerCell];
    mesh->GetCellVertices(neighbourCellIndex, neighbourCellVertices);

    //point is inside a colliding cell, counting it as a contact
    if (PointInCell<Scalar>(neighbourCellVertices, reflectedPoint))
    {
      ValueType neighbourSolution = mesh->GetCellSolution(neighbourCellIndex, reflectedPoint);
      if (!mesh->system.IsProperContact(innerSolution, neighbourSolution, edgeNormal)) return false;
      Scalar innerContactValues[dimsCount];
      Scalar outerContactValues[dimsCount];
      MatrixMulVector(rightContactEdgeMatrix, neighbourSolution.values, outerContactValues, dimsCount, dimsCount);
      MatrixMulVector(leftContactEdgeMatrix, (Scalar*)innerSolution.values, innerContactValues, dimsCount, dimsCount);
      for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
      {
        values[valueIndex] = innerContactValues[valueIndex] + outerContactValues[valueIndex];
      }
      mesh->system.CorrectContact(values, edgeNormal, dynamicInteractionType);
      return true;
    }
    return false;
  }
};
