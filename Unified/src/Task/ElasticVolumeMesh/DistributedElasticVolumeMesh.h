#pragma once

#include "DistributedElasticVolumeMeshCommon.h"

template<typename Space, typename FunctionSpace>
class DistributedElasticVolumeMesh;

template<typename FunctionSpace>
class DistributedElasticVolumeMesh<Space2, FunctionSpace>: public DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>
{
public:
  SPACE2_TYPEDEFS
  typedef Space2 Space;

  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::EdgeLocationPair EdgeLocationPair;
  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::EdgeLocation     EdgeLocation;
  typedef typename ElasticVolumeMesh<Space, FunctionSpace>::EdgeIndices      EdgeIndices;
  typedef typename DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::Node Node;
  typedef typename DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::Cell Cell;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::detectorsPositions;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainIndex;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::syncDataSizes;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::sendNodesInfo;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::MakeSnapshot;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::transitionInfos;

  DistributedElasticVolumeMesh(DifferentialSolver<Scalar>* solver, 
    Scalar tolerance,
    IndexType domainIndex,
    IndexType domainsCount,
    int hierarchyLevelsCount,
    bool sendNodesInfo = false,
    bool sendEdgesInfo = false):
  DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>(solver, 
    tolerance,
    domainIndex,
    domainsCount,
    hierarchyLevelsCount,
    sendNodesInfo),
    sendEdgesInfo(sendEdgesInfo)
  {
  }

  struct EdgeSyncData
  {
    typename GeomMesh<Space>::Edge edge;
    IndexType interactionType;
  };

  void UpdateDomainData(const char* const data);
  void BuildSyncData(IndexType dstDomainIndex, char* const data);

private:
  bool sendEdgesInfo;

  void BuildEdgesSyncData(IndexType dstDomainIndex, char* const data);
  void UpdateEdgesData(IndexType edgesCount, const EdgeSyncData* const edgesData);

  void ComputeSyncDataSizes() override;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/****************************************************
 *                                                  *
 *    DistributedElasticVolumeMesh2.inl             *
 *                                                  *
 ****************************************************/

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::UpdateDomainData(const char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>::UpdateDomainData(data);
  // edges
  IndexType edgesCount = *((IndexType*)(data + offset));
  offset += sizeof(IndexType);
  UpdateEdgesData(edgesCount, (EdgeSyncData*)(data + offset));
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::BuildSyncData(IndexType dstDomainIndex, char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>::BuildSyncData(dstDomainIndex, data);
  if (sendEdgesInfo)
  {
    BuildEdgesSyncData(dstDomainIndex, data + offset);
  } else
  {
    ((IndexType*)(data + offset))[0] = 0;
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::BuildEdgesSyncData(IndexType dstDomainIndex, char* const data)
{
  *((IndexType*)data) = transitionInfos[dstDomainIndex].cells.size() * Space::EdgesPerCell;
  EdgeSyncData* edgesData = (EdgeSyncData*)(data + sizeof(IndexType));

  for (IndexType cellNumber = 0; cellNumber < transitionInfos[dstDomainIndex].cells.size(); ++cellNumber)
  {
    IndexType cellIndex = volumeMesh.GetCellIndex(transitionInfos[dstDomainIndex].cells[cellNumber].incidentNodes);
    assert(cellIndex != IndexType(-1));

    for (IndexType edgeNumber = 0; edgeNumber < Space::EdgesPerCell; ++edgeNumber)
    {
      EdgeSyncData& edgeSyncData = edgesData[cellNumber * Space::EdgesPerCell + edgeNumber];
      edgeSyncData.interactionType = volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[edgeNumber].interactionType;
      volumeMesh.GetCellEdgeNodes(cellIndex, edgeNumber, edgeSyncData.edge.incidentNodes);
    }
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::UpdateEdgesData(IndexType edgesCount, const EdgeSyncData* const edgesData)
{
  for (IndexType edgeNumber = 0; edgeNumber < edgesCount; ++edgeNumber)
  {
    IndexType interactionType = edgesData[edgeNumber].interactionType;
    // TODO
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space2, FunctionSpace>::ComputeSyncDataSizes()
{
  DistributedElasticVolumeMeshCommon<Space2, FunctionSpace>::ComputeSyncDataSizes();
  for (IndexType dstDomainIndex = 0; dstDomainIndex < domainsCount; ++dstDomainIndex)
  {
    syncDataSizes[dstDomainIndex] += 
      sizeof(IndexType) + // edges count
      (sendEdgesInfo ? transitionInfos[dstDomainIndex].cells.size() * Space2::EdgesPerCell * sizeof(EdgeSyncData) : 0);
  }
}


template <typename FunctionSpace>
class DistributedElasticVolumeMesh<Space3, FunctionSpace>: public DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>
{
public:
  SPACE3_TYPEDEFS
  typedef Space3 Space;

  using ElasticVolumeMesh<Space, FunctionSpace>::FindDestructions;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::dimsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::functionsCount;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::volumeMesh;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::detectorsPositions;

  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainIndex;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::domainsCount;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::syncDataSizes;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::sendNodesInfo;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::MakeSnapshot;
  using DistributedElasticVolumeMeshCommon<Space, FunctionSpace>::transitionInfos;

  struct FaceSyncData
  {
    typename GeomMesh<Space>::Face face;
    IndexType interactionType;
  };

  DistributedElasticVolumeMesh(DifferentialSolver<Scalar>* solver, 
    Scalar tolerance,
    IndexType domainIndex,
    IndexType domainsCount,
    int hierarchyLevelsCount,
    bool sendNodesInfo = false,
    bool sendFacesInfo = false):
  DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>(solver, 
    tolerance,
    domainIndex,
    domainsCount,
    hierarchyLevelsCount,
    sendNodesInfo),
    sendFacesInfo(sendFacesInfo)
  {
  }

  void UpdateDomainData(const char* const data);
  void BuildSyncData(IndexType dstDomainIndex, char* const data);

private:
  bool sendFacesInfo;

  void BuildFacesSyncData(IndexType dstDomainIndex, char* const data);
  void UpdateFacesData(IndexType facesCount, const FaceSyncData* const facesData);

  void ComputeSyncDataSizes() override;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


/****************************************************
 *                                                  *
 *    DistributedElasticVolumeMesh3.inl             *
 *                                                  *
 ****************************************************/

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::UpdateDomainData(const char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>::UpdateDomainData(data);
  // faces
  IndexType facesCount = *((IndexType*)(data + offset));
  offset += sizeof(IndexType);
  UpdateFacesData(facesCount, (FaceSyncData*)(data + offset));
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::BuildSyncData(IndexType dstDomainIndex, char* const data)
{
  IndexType offset = DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>::BuildSyncData(dstDomainIndex, data);
  if (sendFacesInfo)
  {
    BuildFacesSyncData(dstDomainIndex, data + offset);
  } else
  {
    ((IndexType*)(data + offset))[0] = 0;
  }
}


template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::BuildFacesSyncData(IndexType dstDomainIndex, char* const data)
{
  *((IndexType*)data) = transitionInfos[dstDomainIndex].cells.size() * Space::FacesPerCell;
  FaceSyncData* facesData = (FaceSyncData*)(data + sizeof(IndexType));

  for (IndexType cellNumber = 0; cellNumber < transitionInfos[dstDomainIndex].cells.size(); ++cellNumber)
  {
    IndexType cellIndex = volumeMesh.GetCellIndex(transitionInfos[dstDomainIndex].cells[cellNumber].incidentNodes);
    assert(cellIndex != IndexType(-1));

    for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
    {
      FaceSyncData& faceSyncData = facesData[cellNumber * Space::FacesPerCell + faceNumber];
      faceSyncData.interactionType = volumeMesh.additionalCellInfos[cellIndex].neighbouringFaces[faceNumber].interactionType;
      volumeMesh.GetCellFaceNodes(cellIndex, faceNumber, faceSyncData.face.incidentNodes);
    }
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::UpdateFacesData(IndexType facesCount, const FaceSyncData* const facesData)
{
  //assert(0); // TODO
  for (IndexType faceNumber = 0; faceNumber < facesCount; ++faceNumber)
  {
    IndexType interactionType = facesData[faceNumber].interactionType;
    // TODO
  }
}

template<typename FunctionSpace>
void DistributedElasticVolumeMesh<Space3, FunctionSpace>::ComputeSyncDataSizes()
{
  DistributedElasticVolumeMeshCommon<Space3, FunctionSpace>::ComputeSyncDataSizes();
  for (IndexType dstDomainIndex = 0; dstDomainIndex < domainsCount; ++dstDomainIndex)
  {
    syncDataSizes[dstDomainIndex] += 
      sizeof(IndexType) + // faces count
      (sendFacesInfo ? transitionInfos[dstDomainIndex].cells.size() * Space3::FacesPerCell * sizeof(FaceSyncData) : 0);
  }
}
