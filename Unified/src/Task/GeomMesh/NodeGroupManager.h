#pragma once

#include <vector>

template <typename Space, typename VolumeMeshType>
class NodeGroupManager
{
public:
  SPACE_TYPEDEFS

  NodeGroupManager(IndexType maxNodeGroupSize): 
    maxNodeGroupSize(maxNodeGroupSize), volumeMesh(0)
  {
  }

  void SetVolumeMesh(VolumeMeshType* volumeMesh)
  {
    this->volumeMesh = volumeMesh;

    nodeGroupPool.resize(volumeMesh->nodes.size() * maxNodeGroupSize, 0);
    nodeGroupSizes.resize(volumeMesh->nodes.size(), 0);
    nodeGroupIndices.resize(volumeMesh->nodes.size(), IndexType(-1));
  }

  void BuildNodeGroups()
  {
    IndexType groupsCount = 0;

    for (IndexType nodeIndex = 0; nodeIndex < volumeMesh->nodes.size(); ++nodeIndex)
    {
      if (nodeGroupIndices[nodeIndex] == IndexType(-1))
      {
        std::vector<IndexType> pool;
        pool.reserve(maxNodeGroupSize);

        nodeGroupSizes[groupsCount] = volumeMesh->GetNodeGroup(nodeIndex, pool);

        for (IndexType nodeNumber = 0; nodeNumber < nodeGroupSizes[groupsCount]; ++nodeNumber)
        {
          nodeGroupIndices[pool[nodeNumber]] = groupsCount;
          nodeGroupPool[groupsCount * maxNodeGroupSize + nodeNumber] = pool[nodeNumber];
        }

        groupsCount++;
      }
    }
  }

  void RemoveNode(IndexType nodeIndex)
  {
    IndexType groupIndex = nodeGroupIndices[nodeIndex];
    assert(groupIndex != IndexType(-1));

    IndexType nodeNumber = 0;
    while (nodeNumber < nodeGroupSizes[groupIndex] && nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber] != nodeIndex)
      nodeNumber++;

    assert(nodeNumber < nodeGroupSizes[groupIndex]);

    std::swap(nodeGroupPool[groupIndex * maxNodeGroupSize + nodeNumber], nodeGroupPool[groupIndex * maxNodeGroupSize + nodeGroupSizes[groupIndex] - 1]);
    nodeGroupPool[groupIndex * maxNodeGroupSize + nodeGroupSizes[groupIndex] - 1] = 0;
    nodeGroupSizes[groupIndex]--;
    nodeGroupIndices[nodeIndex] = IndexType(-1);
  }

  void RemoveCell(IndexType cellIndex)
  {
    IndexType nodeIndices[Space::NodesPerCell];
    volumeMesh->GetFixedCellIndices(cellIndex, nodeIndices);

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      RemoveNode(nodeIndices[nodeNumber]);
    }
  }

  IndexType GetGroupSize(IndexType nodeIndex) const
  {
    IndexType groupIndex = nodeGroupIndices[nodeIndex];
    return groupIndex == IndexType(-1) ? 1 : nodeGroupSizes[nodeIndex];
  }

  const IndexType* GetGroup(IndexType nodeIndex) const
  {
    return nodeGroupPool.data() + nodeIndex * maxNodeGroupSize;
  }

private:
  IndexType maxNodeGroupSize;
  VolumeMeshType* volumeMesh;
  std::vector<IndexType> nodeGroupPool;
  std::vector<IndexType> nodeGroupIndices;
  std::vector<IndexType> nodeGroupSizes;
};
