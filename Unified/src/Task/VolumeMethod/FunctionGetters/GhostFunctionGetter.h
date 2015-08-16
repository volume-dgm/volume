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
    const Vector& testPoint, Scalar* values):
    functionGetter(functionGetter), testPoint(testPoint), values(values), found(false)
  {}
  void operator()(int nodeIndex)
  {
    if (!found)
    {
      IndexType neighbourCellIndex = functionGetter->mesh->aabbTree.GetUserData(nodeIndex);
      if (functionGetter->TryGhostCell(testPoint, neighbourCellIndex, values)) found = true;
    }
  }

  GhostCellFunctionGetter* functionGetter;
  Vector testPoint;
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
  typedef typename System::MediumParameters MediumParameters;
  SPACE_TYPEDEFS

  typedef CollisionProcessor< GhostCellFunctionGetter<MeshType> > CollisionProcessorT;

  const static int dimsCount = MeshType::dimsCount;

  enum GetterType
  {
    Solution, MediumParams
  };

  GhostCellFunctionGetter(
    MeshType* mesh,
    IndexType cellIndex,
    const Vector& edgeNormal,
    Vector* ghostCellVertices,
    Scalar currTime,
    GetterType getterType):
      cellIndex(cellIndex),
      ghostCellVertices(ghostCellVertices),
      edgeNormal(edgeNormal),
      mesh(mesh),
      currTime(currTime),
      getterType(getterType)
  {
    mesh->GetFixedCellIndices(cellIndex, cellIndices);
  }

  IndexType cellIndex;
  IndexType cellIndices[Space::NodesPerCell];
  Vector* ghostCellVertices;
  Vector edgeNormal;
  MeshType* mesh;
  Scalar currTime;
  GetterType getterType;

  void operator()(const Vector& point, Scalar* values)
  {
    switch (getterType)
    {
      case Solution: std::fill(values, values + dimsCount, 0); break;
      case MediumParams: std::fill(values, values + MediumParameters::ParamsCount, 0); break; // lambda, mju, invRho 
    }

    Vector globalPoint = mesh->RefToGlobalVolumeCoords(point, ghostCellVertices) + edgeNormal * mesh->collisionWidth;    
    CollisionProcessorT collisionProcessor(this, globalPoint, values);
    mesh->aabbTree.template FindCollisions<CollisionProcessorT>(AABB(globalPoint, globalPoint), collisionProcessor);

    if (collisionProcessor.found) return;
  }

  bool GetValue(const Vector& globalPoint, Scalar* values)
  {
    switch (getterType)
    {
      case Solution: std::fill(values, values + dimsCount, 0); break;
      case MediumParams: std::fill(values, values + MediumParameters::ParamsCount, 0); break; // lambda, mju, invRho 
    }

    CollisionProcessorT collisionProcessor(this, globalPoint, values);
    mesh->aabbTree.template FindCollisions<CollisionProcessorT>(AABB(globalPoint, globalPoint), collisionProcessor);

    return collisionProcessor.found;
  }

  bool TryGhostCell(const Vector& testPoint, int neighbourCellIndex,
    Scalar* values)
  {
    if (cellIndex == neighbourCellIndex) return false;

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
    if (PointInCell<Scalar>(neighbourCellVertices, testPoint))
    {
      ValueType neighbourSolution = mesh->GetCellSolution(neighbourCellIndex, testPoint);
      switch (getterType)
      {
        case Solution:
          for (IndexType valueIndex = 0; valueIndex < dimsCount; valueIndex++)
          {
            values[valueIndex] = neighbourSolution.values[valueIndex];
          }
        break;
        case MediumParams:
          for (IndexType paramIndex = 0; paramIndex < MediumParameters::ParamsCount; ++paramIndex)
          {
            values[paramIndex] = mesh->cellMediumParameters[neighbourCellIndex].params[paramIndex];
          }
        break;
      }
      return true;
    }
    return false;
  }
};
