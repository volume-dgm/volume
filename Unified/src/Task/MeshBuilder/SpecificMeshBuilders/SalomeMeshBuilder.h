#pragma once

#include "UnstructuredMeshBuilder.h"
#include <vector>
#include <fstream>

template <typename Space>
class SalomeMeshBuilder: public UnstructuredMeshBuilder<Space>
{
public:
  SPACE_TYPEDEFS

  SalomeMeshBuilder(const std::string& meshFileName, bool duplicateNodes = false, Scalar splitNodeBiasRatio = 0.0): 
    UnstructuredMeshBuilder<Space>(duplicateNodes, splitNodeBiasRatio),
    meshFileName(meshFileName)
  {
  }

  virtual void LoadGeom(MeshIO<Space>* const mesh)
  {
    mesh->Load(std::string(AddExtensionToFileName(meshFileName, ".mesh")).c_str(), IO::Ascii);
  }

  virtual void BuildMediumParams(MeshIO<Space>* const mesh, std::vector<char>* mediumParams)
  {
    std::fstream paramsFile(std::string(AddExtensionToFileName(meshFileName, ".params")).c_str(), std::fstream::in | std::fstream::binary);
    
    mediumParams->resize(mesh->GetCellsCount());
    if (paramsFile.fail())
    {
      std::cerr << "Can`t read params file" << std::endl;
      std::fill(mediumParams->begin(), mediumParams->end(), 0);
      return;
    }
    paramsFile.read(mediumParams->data(), mediumParams->size());
  }

  void BuildDetectors(MeshIO<Space>* const mesh)
  {/*
    const Scalar Distance = 50;
    const int DetectorsCount = 131;
    for (int detectorIndex = -DetectorsCount / 2; detectorIndex <= DetectorsCount / 2; ++detectorIndex)
    {
      mesh->detectorsPositions.push_back(Vector(Distance * detectorIndex, Scalar(-1)));
    }*/
    /*
    // Wu
    mesh->detectorsPositions.push_back(Vector(10.2, 12));
    mesh->detectorsPositions.push_back(Vector(6, 12));
    */
    // mesh->detectorsPositions.push_back(Vector(2, 0, 0));
    // mesh->detectorsPositions.push_back(Vector(2, 0, -1));

    //what the fuck?! detectors are supported by salome plugin
    //mesh->detectorsPositions.push_back(Vector(8, 0));
  }

private:
  std::string meshFileName;
};
