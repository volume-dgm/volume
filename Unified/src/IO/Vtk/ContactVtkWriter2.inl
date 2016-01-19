template <typename FunctionSpace>
typename ContactVtkWriter<Space2, FunctionSpace>::OutputData 
  ContactVtkWriter<Space2, FunctionSpace>::ConstructOutputData(const std::vector<Node>& nodes, 
  const std::vector<Cell>& cells, const std::vector<bool>& isCellBroken,
  const std::vector<Scalar>& plasticDeforms,
  Scalar velocityDimensionlessMult,
  const std::vector<EdgePairIndices>&  contactEdges, 
  const std::vector<IndexType>&  contactEdgesCount, IndexType contactTypesCount, 
  const std::vector<bool>& isContactBroken,
  const std::vector<BoundaryEdge>&    boundaryEdges, 
  const std::vector<IndexType>& boundaryEdgesCount, IndexType boundaryTypesCount,
  const ElasticSystem<Space2>& system,
  bool drawContacts,
  bool drawRegularGlueContacts,
  bool drawBrokenContacts,
  // for debuging
  bool drawContactBindings)
{
  OutputData outputData;

  // points
  IndexType additionalNodesCount = 0;
  for (IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; ++contactTypeIndex)
  {
    if (drawContacts && drawContactBindings && (contactTypeIndex > 0 || drawRegularGlueContacts))
    {
      additionalNodesCount += contactEdgesCount[contactTypeIndex] * 2;
    }
  }

  for (IndexType nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex) 
  {
    Node node(nodes[nodeIndex].pos * velocityDimensionlessMult);
    outputData.nodes.push_back(node);
  }

  IndexType offset = 0;
  for (IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; ++contactTypeIndex)
  {
    if (drawContacts && drawContactBindings && (contactTypeIndex > 0 || drawRegularGlueContacts))
    {
      for (IndexType contactEdgeNumber = 0; contactEdgeNumber < contactEdgesCount[contactTypeIndex]; ++contactEdgeNumber)
      {
        IndexType contactEdgeIndex = contactEdgeNumber + offset;
        for (IndexType edgeNumber = 0; edgeNumber < 2; ++edgeNumber) 
        {
          IndexType nodeIndex0 = contactEdges[contactEdgeIndex].edges[edgeNumber].nodeIndices[0];
          IndexType nodeIndex1 = contactEdges[contactEdgeIndex].edges[edgeNumber].nodeIndices[1];
          assert(nodeIndex0 != IndexType(-1) && nodeIndex1 != IndexType(-1));

          Node edgeMiddleNode((nodes[nodeIndex0].pos + nodes[nodeIndex1].pos) * Scalar(0.5) * velocityDimensionlessMult);
          outputData.nodes.push_back(edgeMiddleNode);
        }
      }
    }
    offset += contactEdgesCount[contactTypeIndex];
  }

  // if the edge between v1-v2 is broken then v1 and v2 marks as broken
  outputData.nodeData.resize(nodes.size() + additionalNodesCount, 0);

  // contacts
  if (drawContacts)
  {
    offset = 0;
    for (IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; ++contactTypeIndex)
    {
      if (contactTypeIndex > 0 || drawRegularGlueContacts)
      {
        for (IndexType contactEdgeNumber = 0; contactEdgeNumber < contactEdgesCount[contactTypeIndex]; ++contactEdgeNumber)
        {
          IndexType contactEdgeIndex = contactEdgeNumber + offset;
          bool isEdgeBroken = isContactBroken[contactEdgeIndex];
          for (IndexType edgeNumber = 0; edgeNumber < 2; ++edgeNumber) 
          {
            IndexType nodeIndex0 = contactEdges[contactEdgeIndex].edges[edgeNumber].nodeIndices[0];
            IndexType nodeIndex1 = contactEdges[contactEdgeIndex].edges[edgeNumber].nodeIndices[1];
            assert(nodeIndex0 != IndexType(-1) && nodeIndex1 != IndexType(-1));

            if (isEdgeBroken)
            {
              outputData.nodeData[nodeIndex0] = outputData.nodeData[nodeIndex1] = 1;
            }

            if (drawBrokenContacts && isEdgeBroken)
            {
              Cell cell;
              cell.incidentNodes[0] = nodeIndex0;
              cell.incidentNodes[1] = cell.incidentNodes[2] = nodeIndex1;
              outputData.cells.push_back(cell);
              outputData.cellData.push_back(typename OutputData::CellData(1));
            }
          }
          if (drawContactBindings)
          {
            Cell cell;
            cell.incidentNodes[0] = nodes.size() + 2 * contactEdgeIndex;
            cell.incidentNodes[1] = cell.incidentNodes[2] = nodes.size() + 2 * contactEdgeIndex + 1;
            outputData.cells.push_back(cell);
            outputData.cellData.push_back(typename OutputData::CellData(0));
          }
        }
      }
      offset += contactEdgesCount[contactTypeIndex];
    }
  }

  // boundaries
  offset = 0;
  for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryTypesCount; ++boundaryTypeIndex)
  {
    for (IndexType boundaryEdgeIndex = 0; boundaryEdgeIndex < boundaryEdgesCount[boundaryTypeIndex]; ++boundaryEdgeIndex)
    {
      IndexType nodeIndex0 = boundaryEdges[boundaryEdgeIndex + offset].nodeIndices[0];
      IndexType nodeIndex1 = boundaryEdges[boundaryEdgeIndex + offset].nodeIndices[1];

      if (nodeIndex0 != IndexType(-1) && nodeIndex1 != IndexType(-1) &&
          system.GetBoundaryType(boundaryTypeIndex) != BoundaryConditions::Absorb)
      {
        Cell cell;
        cell.incidentNodes[0] = nodeIndex0;
        cell.incidentNodes[1] = cell.incidentNodes[2] = nodeIndex1;
        outputData.cells.push_back(cell);
        outputData.cellData.push_back(typename OutputData::CellData(0));
      }
    }
    offset += boundaryEdgesCount[boundaryTypeIndex];
  }

  if (drawBrokenContacts)
  {
    for (IndexType cellIndex = 0; cellIndex < cells.size(); ++cellIndex)
    {
      if (isCellBroken[cellIndex])
      {
        typename OutputData::CellData data;
        data.isCellBroken = 1;
        data.plasticDeforms = plasticDeforms[cellIndex];
        outputData.cells.push_back(cells[cellIndex]);
        outputData.cellData.push_back(data);
      }
    }
  }

  return outputData;
}

template<typename FunctionSpace>
void ContactVtkWriter<Space2, FunctionSpace>::Write(const std::string& fileName,
                                                  const std::vector<Node>& nodes, const std::vector<Cell>& cells,
                                                  const std::vector<bool>& isCellBroken,
                                                  const std::vector<Scalar>& plasticDeforms,
                                                  Scalar velocityDimensionlessMult,
                                                  const std::vector<EdgePairIndices>& contactEdges, 
                                                  const std::vector<IndexType>& contactEdgesCount, 
                                                  IndexType contactTypesCount,
                                                  const std::vector<bool>& isContactBroken,
                                                  const std::vector<BoundaryEdge>&    boundaryEdges, 
                                                  const std::vector<IndexType>& boundaryEdgesCount, 
                                                  IndexType boundaryTypesCount,
                                                  const ElasticSystem<Space2>& system,
                                                  bool drawContacts, bool drawRegularGlueContacts,
                                                  bool drawBrokenContacts, bool drawContactBindings)
{
  OutputData outputData = ConstructOutputData(nodes, cells, isCellBroken, plasticDeforms, velocityDimensionlessMult,
    contactEdges, contactEdgesCount, contactTypesCount, isContactBroken,
    boundaryEdges, boundaryEdgesCount, boundaryTypesCount,
    system, drawContacts, drawRegularGlueContacts, drawBrokenContacts, drawContactBindings);

  SaveToFile(fileName, outputData);
}

template<typename FunctionSpace>
typename ContactVtkWriter<Space2, FunctionSpace>::OutputData
  ContactVtkWriter<Space2, FunctionSpace>::ConstructOutputData(
  ElasticVolumeMesh<Space2, FunctionSpace>* mesh,
  const ElasticSystem<Space2>& system)
{
  OutputData outputData;

  for (IndexType nodeIndex = 0; nodeIndex < mesh->volumeMesh.nodes.size(); ++nodeIndex)
  {
    Node node(mesh->volumeMesh.nodes[nodeIndex].pos * mesh->velocityDimensionlessMult);
    outputData.nodes.push_back(node);
  }

  // boundaries
  for (IndexType cellIndex = 0; cellIndex < mesh->volumeMesh.cells.size(); ++cellIndex)
  {
    if (mesh->volumeMesh.isCellAvailable[cellIndex])
    {
      for (IndexType faceNumber = 0; faceNumber < Space2::FacesPerCell; ++faceNumber)
      {
        IndexType correspondingCellIndex = mesh->volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[faceNumber].correspondingCellIndex;
        IndexType correspondingEdgeNumber = mesh->volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[faceNumber].correspondingEdgeNumber;
        IndexType interactionType = mesh->volumeMesh.additionalCellInfos[cellIndex].neighbouringEdges[faceNumber].interactionType;

        if (correspondingCellIndex == IndexType(-1) && correspondingEdgeNumber == IndexType(-1) &&
          system.GetBoundaryType(interactionType) != BoundaryConditions::Absorb)
        {
          IndexType edgeNodeNumbers[Space2::NodesPerFace];
          mesh->volumeMesh.GetCellEdgeNodes(cellIndex, faceNumber, edgeNodeNumbers);

          Cell cell;
          cell.incidentNodes[0] = edgeNodeNumbers[0];
          cell.incidentNodes[1] = cell.incidentNodes[2] = edgeNodeNumbers[1];
          outputData.cells.push_back(cell);
          outputData.cellData.push_back(typename OutputData::CellData(0));
        }
      }
    }
  }

  if (mesh->isCellBroken.size() == mesh->volumeMesh.cells.size())
  {
    for (IndexType cellIndex = 0; cellIndex < mesh->volumeMesh.cells.size(); ++cellIndex)
    {
      if (mesh->isCellBroken[cellIndex] && mesh->volumeMesh.isCellAvailable[cellIndex])
      {
        typename OutputData::CellData data;
        data.isCellBroken = 1;
        data.plasticDeforms = mesh->plasticDeforms[cellIndex];
        outputData.cells.push_back(mesh->volumeMesh.cells[cellIndex]);
        outputData.cellData.push_back(data);
      }
    }
  }

  return outputData;
}

template <typename FunctionSpace>
void ContactVtkWriter<Space2, FunctionSpace>::Write(const std::string& fileName,
  ElasticVolumeMesh<Space2, FunctionSpace>* mesh,
  const ElasticSystem<Space2>& system)
{
  OutputData outputData = ConstructOutputData(mesh, system);
  SaveToFile(fileName, outputData);
}

template <typename FunctionSpace>
void ContactVtkWriter<Space2, FunctionSpace>::SaveToFile(const std::string& fileName, const OutputData& outputData) const
{
  fstream file(fileName.c_str(), fstream::out | fstream::binary);

  if (file.fail()) {
    return;
  }

  // header lines
  file << "# vtk DataFile Version 4.0\n";
  file << "vtk output\n";
  file << "ASCII\n";

  file << "DATASET POLYDATA\n";

  file << "POINTS " << outputData.nodes.size()
    << " " << TypeName() << "\n";
  for (IndexType nodeIndex = 0; nodeIndex < outputData.nodes.size(); ++nodeIndex)
  {
    file << outputData.nodes[nodeIndex].pos.x << " "
      << outputData.nodes[nodeIndex].pos.y << " "
      << 0 << " ";
  }
  file << "\n";

  file << "POLYGONS " << outputData.cells.size() << " "
    << outputData.cells.size() * (3 /*vertex count in a cell*/
      + 1 /*count of connectivity indicies*/) << "\n";

  for (IndexType cellIndex = 0; cellIndex < outputData.cells.size(); ++cellIndex)
  {
    file << 3 << " ";
    for (IndexType nodeNumber = 0; nodeNumber < 3; ++nodeNumber)
    {
      file << outputData.cells[cellIndex].incidentNodes[nodeNumber] << " ";
    }
  }
  file << "\n";

  file << "CELL_DATA " << outputData.cellData.size() << "\n";
  file << "SCALARS ISBROKEN INT 1\n";
  file << "LOOKUP_TABLE default\n";
  for (IndexType cellIndex = 0; cellIndex < outputData.cellData.size(); ++cellIndex)
  {
    file << outputData.cellData[cellIndex].isCellBroken << std::endl;
  }

  file << "SCALARS PLASTIC_DEFORMS float 1\n";
  file << "LOOKUP_TABLE default\n";
  for (IndexType cellIndex = 0; cellIndex < outputData.cellData.size(); ++cellIndex)
  {
    file << outputData.cellData[cellIndex].plasticDeforms << std::endl;
  }
}
