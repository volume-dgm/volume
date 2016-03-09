#pragma once

#include "../../../Maths/Collisions.h"

template <typename ContactFinder>
struct ContactProcessor
{
  typedef typename ContactFinder::MeshTypeT MeshType;
  typedef typename MeshType::Space          Space;
  typedef typename MeshType::SystemT        System;
  typedef typename System::ValueType        ValueType;
  SPACE_TYPEDEFS

  const static int dimsCount = MeshType::dimsCount;

  ContactProcessor(ContactFinder* contactFinder,
    IndexType cellIndex, IndexType* contactedCells, IndexType* contactedCellsCount):
    contactFinder(contactFinder), cellIndex(cellIndex), contactedCells(contactedCells), contactedCellsCount(contactedCellsCount)
  {
    (*contactedCellsCount) = 0;
  }

  void operator()(int nodeIndex)
  {
    IndexType neighbourCellIndex = contactFinder->mesh->aabbTree.GetUserData(nodeIndex);
    if (contactFinder->TryCell(cellIndex, neighbourCellIndex))
    {
      if (contactedCellsCount)
      {
        if (contactedCells)
          contactedCells[*contactedCellsCount] = neighbourCellIndex;

        (*contactedCellsCount)++;
      }
    }
  }

  ContactFinder* contactFinder;
  IndexType cellIndex;
  IndexType* contactedCells; // pointer to array of colliding cells
  IndexType* contactedCellsCount;
};


template<typename MeshType>
struct ContactFinder
{
  typedef typename MeshType::Space          Space;
  typedef typename MeshType::SystemT        System;
  typedef MeshType MeshTypeT;
  typedef typename System::ValueType        ValueType;
  typedef typename System::MediumParameters MediumParameters;
  SPACE_TYPEDEFS

  typedef ContactProcessor< ContactFinder<MeshType> > ContactProcessorT;

  const static int dimsCount = MeshType::dimsCount;

  ContactFinder(
    MeshType* mesh,
    IndexType cellIndex):
    mesh(mesh),
    cellIndex(cellIndex)
  {
    mesh->GetFixedCellIndices(cellIndex, cellIndices);
    mesh->GetCellVertices(cellIndices, cellVertices);
  }

  MeshType* mesh;
  IndexType cellIndex;
  IndexType cellIndices[Space::NodesPerCell];
  Vector cellVertices[Space::NodesPerCell];

  void Find(IndexType* collidedCells, IndexType* collidedCellsCount)
  {
    AABB aabb = mesh->GetCellAABB(cellIndex);
    ContactProcessorT contactProcessor(this, cellIndex, collidedCells, collidedCellsCount);
    mesh->aabbTree.template FindCollisions<ContactProcessorT>(aabb, contactProcessor);
  }

  bool TryCell(IndexType cellIndex, IndexType neighbourCellIndex)
  {
    if (cellIndex == neighbourCellIndex) return false;

    Vector neighbourCellVertices[Space::NodesPerCell];
    mesh->GetCellVertices(neighbourCellIndex, neighbourCellVertices);

    return CellsCollide(cellVertices, neighbourCellVertices, 0);
    // return CellsCollide(cellVertices, neighbourCellVertices); // for 2d only
  }
};

