#pragma once

#include "MeshIOBase.h"
#include "../../../../IO/Vtk/MeshVtkWriter.h"
#include "../../../../Utils/Utils.h"

template <typename Space>
struct MeshIO: public MeshIOBase<Space>
{
  SPACE_TYPEDEFS
  MeshIO();
  virtual ~MeshIO() = default;
  virtual void Save(const std::string& fileName, IO::FileType fileType = IO::Binary);
  virtual void Load(const std::string& fileName, IO::FileType fileType = IO::Binary);

  // file should be opened
  void Save(std::fstream& file, IO::FileType fileType);
  void Load(std::fstream& file, IO::FileType fileType);

  std::vector<Vector>    vertices;
  std::vector<IndexType> indices;

  IndexType contactTypesCount;
  IndexType boundaryTypesCount;

  std::vector<Vector> detectorsPositions;

  IndexType GetCellsCount() const;
  Vector GetCellCenter(IndexType cellIndex) const;
  
  bool operator==(const MeshIO<Space>& other) const;

  void SaveContacts(const std::string& vtkFileName);
  void SaveBoundaries(const std::string& vtkFileName);
  void SaveDetectors(const std::string& vtkFileName);

private:
  void SaveContacts(std::fstream& file, IO::FileType fileType);
  void SaveBoundaries(std::fstream& file, IO::FileType fileType);
  void LoadContacts(std::fstream& file, IO::FileType fileType);
  void LoadBoundaries(std::fstream& file, IO::FileType fileType);
};

/****************************************************
 *                                                  *
 *                  MeshIO.inl                      *
 *                                                  *
 ****************************************************/

template <typename Space>
MeshIO<Space>::MeshIO(): 
  contactTypesCount(0),
  boundaryTypesCount(0)
{
}

template <typename Space>
void MeshIO<Space>::Save(const std::string& fileName, IO::FileType fileType)
{
  std::fstream file;

  switch (fileType)
  {
    case IO::Ascii: file.open(fileName.c_str(), std::fstream::out); break;
    case IO::Binary: file.open(fileName.c_str(), std::fstream::out | std::fstream::binary); break;
  }
  
  if (file.fail())
  {
    std::cerr << "Can`t open file " << fileName << " for saving" << std::endl;
    return;
  }

  Save(file, fileType);
  file.close();
}

template <typename Space>
void MeshIO<Space>::Load(const std::string& fileName, IO::FileType fileType)
{
  std::fstream file;

  switch (fileType)
  {
    case IO::Ascii: file.open(fileName.c_str(), std::fstream::in); break;
    case IO::Binary: file.open(fileName.c_str(), std::fstream::in | std::fstream::binary); break;
  }
  
  if (file.fail())
  {
    std::cerr << "Can`t open file " << fileName << " for loading" << std::endl;
    throw;
  }

  Load(file, fileType);
  file.close();
}

template <typename Space>
typename Space::IndexType MeshIO<Space>::GetCellsCount() const
{
  return indices.size() / Space::NodesPerCell;
}

template <typename Space>
typename Space::Vector MeshIO<Space>::GetCellCenter(IndexType cellIndex) const
{
  Vector center = Vector::zero();
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    center += vertices[indices[cellIndex * Space::NodesPerCell + nodeNumber]];
  }
  center /= Scalar(Space::NodesPerCell);
  return center;
}

template <typename Space>
bool MeshIO<Space>::operator==(const MeshIO<Space>& other) const
{
  if (vertices.size() != other.vertices.size()) return false;

  const Scalar eps = 1e-3;
  for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
  {
    Scalar meanLen = (vertices[vertexIndex] + other.vertices[vertexIndex]).Len() * Scalar(0.5);
    if ((vertices[vertexIndex] - other.vertices[vertexIndex]).Len() > meanLen * eps) return false;
  }

  if (detectorsPositions.size() != other.detectorsPositions.size()) return false;
  for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
  {
    Scalar meanLen = (detectorsPositions[detectorIndex] + other.detectorsPositions[detectorIndex]).Len() * Scalar(0.5);
    if ((detectorsPositions[detectorIndex] - other.detectorsPositions[detectorIndex]).Len() > meanLen * eps) return false;
  }

  return MeshIOBase<Space>::operator==(other) && 
    indices == other.indices &&
    contactTypesCount == other.contactTypesCount &&
    boundaryTypesCount == other.boundaryTypesCount;
}

// file should be opened
template <typename Space>
void MeshIO<Space>::Save(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      // vertices
      file << vertices.size() << std::endl;
      for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file << vertices[vertexIndex][componentNumber] << " ";
        }
        file << std::endl;
      }
      file << std::endl;

      // indices
      file << indices.size() << std::endl;
      for (IndexType index = 0; index < indices.size(); ++index)
      {
        file << indices[index] << " ";
      }
      file << std::endl;

      // contacts
      SaveContacts(file, IO::Ascii);
      std::cout << std::endl;

      // boundaries
      SaveBoundaries(file, IO::Ascii);
      std::cout << std::endl;

      // detectors
      file << detectorsPositions.size() << std::endl;
      for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file << detectorsPositions[detectorIndex][componentNumber] << " ";
        }
        file << std::endl;
      }
      file << std::endl;
    break;
    case IO::Binary:
      // vertices
      IO::Write(file, vertices.size());
      IO::WriteVector(file, vertices);

      IO::Write(file, indices.size());
      IO::WriteVector(file, indices);

      // contacts
      SaveContacts(file, IO::Binary);

      // boundaries
      SaveBoundaries(file, IO::Binary);

      // detectors
      IO::Write(file, detectorsPositions.size());
      IO::WriteVector(file, detectorsPositions);
    break;
  }
}

template <typename Space>
void MeshIO<Space>::Load(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
    {
      // vertices
      IndexType verticesCount;
      file >> verticesCount;
      vertices.resize(verticesCount);

      for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file >> vertices[vertexIndex][componentNumber];
        }
      }

      // indices
      IndexType indicesCount;
      file >> indicesCount;
      indices.resize(indicesCount);
      for (IndexType index = 0; index < indices.size(); ++index)
      {
        file >> indices[index];
      }

      // contacts
      LoadContacts(file, IO::Ascii);

      // boundaries
      LoadBoundaries(file, IO::Ascii);

      // detectors
      IndexType detectorsCount;
      file >> detectorsCount;
      detectorsPositions.resize(detectorsCount);
      for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file >> detectorsPositions[detectorIndex][componentNumber];
        }
      }
    } break;

    case IO::Binary:
    {
      // vertices
      IndexType verticesCount;
      IO::Read(file, verticesCount);
      vertices.resize(verticesCount);
      IO::Read(file, vertices.data(), verticesCount);

      // indices
      IndexType indicesCount;
      IO::Read(file, indicesCount);
      indices.resize(indicesCount);
      IO::Read(file, indices.data(), indicesCount);

      // contacts
      LoadContacts(file, IO::Binary);

      // boundaries
      LoadBoundaries(file, IO::Binary);

      // detectors
      IndexType detectorsCount;
      IO::Read(file, detectorsCount);
      detectorsPositions.resize(detectorsCount);
      IO::Read(file, detectorsPositions.data(), detectorsCount);
    } break;
  }
}

template <typename Space>
void MeshIO<Space>::SaveDetectors(const std::string& vtkFileName)
{
  Scalar minY = vertices[0].y;
  Scalar maxY = vertices[0].y;
  for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
  {
    minY = std::min(minY, vertices[vertexIndex].y);
    maxY = std::max(maxY, vertices[vertexIndex].y);
  }

  std::vector<Vector> vertices;
  std::vector<IndexType> indices;

  const Scalar Height = (maxY - minY) / 50;

  for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
  {
    vertices.push_back(detectorsPositions[detectorIndex]);
    vertices.push_back(detectorsPositions[detectorIndex] + Vector::yAxis() * Height);
    indices.push_back(2 * detectorIndex + 0);
    indices.push_back(2 * detectorIndex + 1);
  }

  MeshVtkWriter<Space> writer;
  writer.WriteFaces(vtkFileName, vertices, indices);
}

/****************************************************
 *                                                  *
 *                  MeshIO2.inl                     *
 *                                                  *
 ****************************************************/

template <>
void MeshIO<Space2>::SaveContacts(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      using std::endl;
      file << contactEdges.size() << endl;
      for (IndexType contactEdgesIndex = 0; contactEdgesIndex < contactEdges.size(); ++contactEdgesIndex)
      {
        for (IndexType edgePairNumber = 0; edgePairNumber < 2; ++edgePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
          {
            file << contactEdges[contactEdgesIndex].edges[edgePairNumber].nodeIndices[nodeNumber] << " ";
          }
        }
      }
      file << endl;

      assert(contactEdgesCount.size() == contactTypesCount);

      file << contactEdgesCount.size() << endl;
      for (IndexType index = 0; index < contactEdgesCount.size(); ++index)
      {
        file << contactEdgesCount[index] << " ";
      } 
      file << endl;
    break;
    case IO::Binary:
      IO::Write(file, contactEdges.size());
      IO::WriteVector(file, contactEdges);

      assert(contactEdgesCount.size() == contactTypesCount);
      IO::Write(file, contactEdgesCount.size());
      IO::WriteVector(file, contactEdgesCount);
    break;
  }
}

template <>
void MeshIO<Space2>::SaveBoundaries(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      using std::endl;
      file << boundaryEdges.size() << endl;
      for (IndexType boundaryEdgesIndex = 0; boundaryEdgesIndex < boundaryEdges.size(); ++boundaryEdgesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
        {
          file << boundaryEdges[boundaryEdgesIndex].nodeIndices[nodeNumber] << " ";
        }
      }
      file << endl;

      assert(boundaryEdgesCount.size() == boundaryTypesCount);
      file << boundaryEdgesCount.size() << endl;
      for (IndexType index = 0; index < boundaryEdgesCount.size(); ++index)
      {
        file << boundaryEdgesCount[index] << " ";
      } 
      file << endl;
    break;
    case IO::Binary:
      IO::Write(file, boundaryEdges.size());
      IO::WriteVector(file, boundaryEdges);

      assert(boundaryEdgesCount.size() == boundaryTypesCount);
      IO::Write(file, boundaryEdgesCount.size());
      IO::WriteVector(file, boundaryEdgesCount);
    break;
  }
}

template <>
void MeshIO<Space2>::LoadContacts(std::fstream& file, IO::FileType fileType)
{
  IndexType contactEdgesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> contactEdgesSize;
      contactEdges.resize(contactEdgesSize);

      for (IndexType contactEdgesIndex = 0; contactEdgesIndex < contactEdges.size(); ++contactEdgesIndex)
      {
        for (IndexType edgePairNumber = 0; edgePairNumber < 2; ++edgePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
          {
            file >> contactEdges[contactEdgesIndex].edges[edgePairNumber].nodeIndices[nodeNumber];
          }
        }
      }

      file >> contactTypesCount;
      contactEdgesCount.resize(contactTypesCount);
      for (IndexType index = 0; index < contactEdgesCount.size(); ++index)
      {
        file >> contactEdgesCount[index];
      }
    break;
    case IO::Binary:
      IO::Read(file, contactEdgesSize);
      contactEdges.resize(contactEdgesSize);
      IO::Read(file, contactEdges.data(), contactEdgesSize);

      IO::Read(file, contactTypesCount);
      contactEdgesCount.resize(contactTypesCount);
      IO::Read(file, contactEdgesCount.data(), contactEdgesCount.size());
    break;
  }
}

template <>
void MeshIO<Space2>::LoadBoundaries(std::fstream& file, IO::FileType fileType)
{
  IndexType boundaryEdgesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> boundaryEdgesSize;
      boundaryEdges.resize(boundaryEdgesSize);

      for (IndexType boundaryEdgesIndex = 0; boundaryEdgesIndex < boundaryEdges.size(); ++boundaryEdgesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
        {
          file >> boundaryEdges[boundaryEdgesIndex].nodeIndices[nodeNumber];
        }
      }

      file >> boundaryTypesCount;
      boundaryEdgesCount.resize(boundaryTypesCount);
      for (IndexType index = 0; index < boundaryEdgesCount.size(); ++index)
      {
        file >> boundaryEdgesCount[index];
      }
    break;
    case IO::Binary:
      IO::Read(file, boundaryEdgesSize);
      boundaryEdges.resize(boundaryEdgesSize);
      IO::Read(file, boundaryEdges.data(), boundaryEdgesSize);

      IO::Read(file, boundaryTypesCount);
      boundaryEdgesCount.resize(boundaryTypesCount);
      IO::Read(file, boundaryEdgesCount.data(), boundaryTypesCount);
    break;
  }
}

template <>
void MeshIO<Space2>::SaveContacts(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space2>::AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> contactTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType contactType = 0; contactType < contactTypesCount; ++contactType)
  {
    for (IndexType contactNumber = 0; contactNumber < contactEdgesCount[contactType]; ++contactNumber)
    {
      EdgePairIndices contactEdge = contactEdges[offset + contactNumber];
      for (IndexType edgeNumber = 0; edgeNumber < 2; ++edgeNumber)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
        {
          indices.push_back(contactEdge.edges[edgeNumber].nodeIndices[nodeNumber]);
        }
        contactTypes.push_back(contactType);
      }
    }
    offset += contactEdgesCount[contactType];
  }
  cellInfos.resize(contactTypes.size());
  for (IndexType contactIndex = 0; contactIndex < cellInfos.size(); ++contactIndex)
  {
    cellInfos[contactIndex].count = 1;
    cellInfos[contactIndex].data = &(contactTypes[contactIndex]);
  }
  MeshVtkWriter<Space2, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}

template <>
void MeshIO<Space2>::SaveBoundaries(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space2>::AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> boundaryTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType boundaryType = 0; boundaryType < boundaryTypesCount; ++boundaryType)
  {
    for (IndexType boundaryNumber = 0; boundaryNumber < boundaryEdgesCount[boundaryType]; ++boundaryNumber)
    {
      BoundaryEdge boundaryEdge = boundaryEdges[offset + boundaryNumber];
      for (IndexType nodeNumber = 0; nodeNumber < Space2::NodesPerEdge; ++nodeNumber)
      {
        indices.push_back(boundaryEdge.nodeIndices[nodeNumber]);
      }
      boundaryTypes.push_back(boundaryType);
    }
    offset += boundaryEdgesCount[boundaryType];
  }
  cellInfos.resize(boundaryTypes.size());
  for (IndexType boundaryIndex = 0; boundaryIndex < cellInfos.size(); ++boundaryIndex)
  {
    cellInfos[boundaryIndex].count = 1;
    cellInfos[boundaryIndex].data = &(boundaryTypes[boundaryIndex]);
  }
  MeshVtkWriter<Space2, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}

/****************************************************
 *                                                  *
 *                  MeshIO3.inl                     *
 *                                                  *
 ****************************************************/

template <>
void MeshIO<Space3>::SaveContacts(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      file << contactFaces.size() << std::endl;
      for (IndexType contactFacesIndex = 0; contactFacesIndex < contactFaces.size(); ++contactFacesIndex)
      {
        for (IndexType facePairNumber = 0; facePairNumber < 2; ++facePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
          {
            file << contactFaces[contactFacesIndex].faces[facePairNumber].nodeIndices[nodeNumber] << " ";
          }
        }
      }
      file << std::endl;

      assert(contactFacesCount.size() == contactTypesCount);

      file << contactFacesCount.size() << std::endl;
      for (IndexType contactTypeIndex = 0; contactTypeIndex < contactFacesCount.size(); ++contactTypeIndex)
      {
        file << contactFacesCount[contactTypeIndex] << " ";
      } 
      file << std::endl;
    break;
    case IO::Binary:
      IO::Write(file, contactFaces.size());
      IO::WriteVector(file, contactFaces);
      
      assert(contactTypesCount == contactFacesCount.size());
      IO::Write(file, contactTypesCount);
      IO::WriteVector(file, contactFacesCount);
    break;
  }
}

template <>
void MeshIO<Space3>::SaveBoundaries(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
      file << boundaryFaces.size() << std::endl;
      for (IndexType boundaryFacesIndex = 0; boundaryFacesIndex < boundaryFaces.size(); ++boundaryFacesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
        {
          file << boundaryFaces[boundaryFacesIndex].nodeIndices[nodeNumber] << " ";
        }
      }
      file << std::endl;

      assert(boundaryFacesCount.size() == boundaryTypesCount);
      file << boundaryFacesCount.size() << std::endl;
      for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryFacesCount.size(); ++boundaryTypeIndex)
      {
        file << boundaryFacesCount[boundaryTypeIndex] << " ";
      } 
      file << std::endl;
    break;
    case IO::Binary:
      IO::Write(file, boundaryFaces.size());
      IO::WriteVector(file, boundaryFaces);

      assert(boundaryFacesCount.size() == boundaryTypesCount);
      IO::Write(file, boundaryFacesCount.size());
      IO::WriteVector(file, boundaryFacesCount);
    break;
  }
}

template <>
void MeshIO<Space3>::LoadContacts(std::fstream& file, IO::FileType fileType)
{
  IndexType contactFacesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> contactFacesSize;
      contactFaces.resize(contactFacesSize);

      for (IndexType contactFacesIndex = 0; contactFacesIndex < contactFaces.size(); ++contactFacesIndex)
      {
        for (IndexType facePairNumber = 0; facePairNumber < 2; ++facePairNumber)
        {
          for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
          {
            file >> contactFaces[contactFacesIndex].faces[facePairNumber].nodeIndices[nodeNumber];
          }
        }
      }

      file >> contactTypesCount;
      contactFacesCount.resize(contactTypesCount);
      for (IndexType contactTypeIndex = 0; contactTypeIndex < contactFacesCount.size(); ++contactTypeIndex)
      {
        file >> contactFacesCount[contactTypeIndex];
      }
    break;
    case IO::Binary:
      IO::Read(file, contactFacesSize);
      contactFaces.resize(contactFacesSize);
      IO::Read(file, contactFaces.data(), contactFacesSize);

      IO::Read(file, contactTypesCount);
      contactFacesCount.resize(contactTypesCount);
      IO::Read(file, contactFacesCount.data(), contactFacesCount.size());
    break;
  }
}

template <>
void MeshIO<Space3>::LoadBoundaries(std::fstream& file, IO::FileType fileType)
{
  IndexType boundaryFacesSize;
  switch (fileType)
  {
    case IO::Ascii:
      file >> boundaryFacesSize;
      boundaryFaces.resize(boundaryFacesSize);

      for (IndexType boundaryFacesIndex = 0; boundaryFacesIndex < boundaryFaces.size(); ++boundaryFacesIndex)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
        {
          file >> boundaryFaces[boundaryFacesIndex].nodeIndices[nodeNumber];
        }
      }

      file >> boundaryTypesCount;
      boundaryFacesCount.resize(boundaryTypesCount);
      for (IndexType boundaryTypeIndex = 0; boundaryTypeIndex < boundaryFacesCount.size(); ++boundaryTypeIndex)
      {
        file >> boundaryFacesCount[boundaryTypeIndex];
      }
    break;
    case IO::Binary:
      IO::Read(file, boundaryFacesSize);
      boundaryFaces.resize(boundaryFacesSize);
      IO::Read(file, boundaryFaces.data(), boundaryFacesSize);

      IO::Read(file, boundaryTypesCount);
      boundaryFacesCount.resize(boundaryTypesCount);
      IO::Read(file, boundaryFacesCount.data(), boundaryTypesCount);
    break;
  }
}

template <>
void MeshIO<Space3>::SaveContacts(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space3>::AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> contactTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType contactType = 0; contactType < contactTypesCount; ++contactType)
  {
    for (IndexType contactNumber = 0; contactNumber < contactFacesCount[contactType]; ++contactNumber)
    {
      FacePairIndices contactFace = contactFaces[offset + contactNumber];
      for (IndexType faceNumber = 0; faceNumber < 2; ++faceNumber)
      {
        for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
        {
          indices.push_back(contactFace.faces[faceNumber].nodeIndices[nodeNumber]);
        }
        contactTypes.push_back(contactType);
      }
    }
    offset += contactFacesCount[contactType];
  }
  cellInfos.resize(contactTypes.size());
  for (IndexType contactIndex = 0; contactIndex < cellInfos.size(); ++contactIndex)
  {
    cellInfos[contactIndex].count = 1;
    cellInfos[contactIndex].data = &(contactTypes[contactIndex]);
  }
  MeshVtkWriter<Space3, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}

template <>
void MeshIO<Space3>::SaveBoundaries(const std::string& vtkFileName)
{
  typedef AdditionalCellInfo<Space3>:: AuxInfo<IndexType> CellInfo;
  std::vector<CellInfo>  cellInfos;
  std::vector<IndexType> boundaryTypes;
  std::vector<IndexType> indices;
  IndexType offset = 0;
  for (IndexType boundaryType = 0; boundaryType < boundaryTypesCount; ++boundaryType)
  {
    for (IndexType boundaryNumber = 0; boundaryNumber < boundaryFacesCount[boundaryType]; ++boundaryNumber)
    {
      BoundaryFace boundaryFace = boundaryFaces[offset + boundaryNumber];
      for (IndexType nodeNumber = 0; nodeNumber < Space3::NodesPerFace; ++nodeNumber)
      {
        indices.push_back(boundaryFace.nodeIndices[nodeNumber]);
      }
      boundaryTypes.push_back(boundaryType);
    }
    offset += boundaryFacesCount[boundaryType];
  }
  cellInfos.resize(boundaryTypes.size());
  for (IndexType boundaryIndex = 0; boundaryIndex < cellInfos.size(); ++boundaryIndex)
  {
    cellInfos[boundaryIndex].count = 1;
    cellInfos[boundaryIndex].data = &(boundaryTypes[boundaryIndex]);
  }
  MeshVtkWriter<Space3, CellInfo> writer;
  writer.WriteFaces(vtkFileName, vertices, indices, cellInfos);
}
