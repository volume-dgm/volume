template <typename Space, typename FunctionSpace>
void ContactVtkWriterBase<Space, FunctionSpace>::Write(const std::string& fileName,
  ElasticVolumeMesh<Space, FunctionSpace>* mesh,
  const ElasticSystem<Space>& system)
{
  OutputData outputData = ConstructOutputData(mesh, system);
  SaveToFile(fileName, outputData);
}

template<typename Space, typename FunctionSpace>
typename ContactVtkWriterBase<Space, FunctionSpace>::OutputData
ContactVtkWriterBase<Space, FunctionSpace>::ConstructOutputData(
  ElasticVolumeMesh<Space, FunctionSpace>* mesh,
  const ElasticSystem<Space>& system)
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
      for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
      {
        IndexType correspondingCellIndex = mesh->volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);
        IndexType correspondingFaceNumber = mesh->volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber);
        IndexType interactionType = mesh->volumeMesh.GetInteractionType(cellIndex, faceNumber);

        if (correspondingCellIndex == IndexType(-1) && correspondingFaceNumber == IndexType(-1) &&
          system.GetBoundaryType(interactionType) != BoundaryConditions::Absorb)
        {
          IndexType faceNodeNumbers[Space::NodesPerFace];
          mesh->volumeMesh.GetCellFaceNodes(cellIndex, faceNumber, faceNodeNumbers);

          Cell cell;
          for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
          {
            cell.incidentNodes[nodeNumber] = faceNodeNumbers[nodeNumber % Space::NodesPerFace];
          }

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
      if (mesh->volumeMesh.isCellAvailable[cellIndex])
      {
        typename OutputData::CellData data;
        data.isCellBroken = mesh->isCellBroken[cellIndex] ? 1 : 0;
        data.plasticDeforms = mesh->plasticDeforms.size() ? mesh->plasticDeforms[cellIndex] : 0;
        outputData.cells.push_back(mesh->volumeMesh.cells[cellIndex]);
        outputData.cellData.push_back(data);
      }
    }
  }

  return outputData;
}

template <typename Space, typename FunctionSpace>
void ContactVtkWriterBase<Space, FunctionSpace>::SaveToFile(const std::string& fileName, const OutputData& outputData) const
{
  std::fstream file(fileName.c_str(), std::fstream::out | std::fstream::binary);

  if (file.fail())
  {
    return;
  }

  // header lines
  file << "# vtk DataFile Version 4.0\n";
  file << "vtk output\n";
  file << "ASCII\n";

  file << "DATASET UNSTRUCTURED_GRID\n";

  // points
  file << "POINTS " << outputData.nodes.size() << " " << TypeName() << "\n";
  for (IndexType index = 0; index < outputData.nodes.size(); ++index)
  {
    for (IndexType coordIndex = 0; coordIndex < 3; ++coordIndex)
    {
      file << outputData.nodes[index].pos.Get(coordIndex) << " ";
    }
  }
  file << "\n";

  // cells
  file << "CELLS " << outputData.cells.size() << " "
    << outputData.cells.size() * (Space::NodesPerCell /*vertex count in a cell*/
      + 1 /*count of connectivity indices*/) << "\n";

  for (IndexType cellIndex = 0; cellIndex < outputData.cells.size(); ++cellIndex)
  {
    file << Space::NodesPerCell;
    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      file << " " << outputData.cells[cellIndex].incidentNodes[nodeNumber];
    }
    file << "\n";
  }

  file << "CELL_TYPES " << outputData.cells.size() << "\n";
  for (IndexType cellIndex = 0; cellIndex < outputData.cells.size(); ++cellIndex)
  {
    if (Space::NodesPerCell == 4)
    {
      file << 10; // tetrahedron
    } else
    if (Space::NodesPerCell == 3)
    {
      file << 5; // triangle
    }
    file << std::endl;
  }

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

  file.close();
}
