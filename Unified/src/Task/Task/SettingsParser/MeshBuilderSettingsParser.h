#pragma once

#include "ParserUtil.h"
#include "../../../Maths/Spaces.h"

template <typename Space>
struct MeshBuilderSettings
{
  SPACE_TYPEDEFS

  MeshBuilderSettings(): 
    partitionAlgorithmName("METIS_PartMeshDual"),  // METIS_PartMeshNodal
    salomeFileName("none")
  {}

  void Parse(TiXmlElement *meshBuilderElement)
  {
    TiXmlElement* partitionInfoElement = meshBuilderElement->FirstChildElement("Partition");
    assert(partitionInfoElement);

    ParseString(partitionInfoElement, "algorithm", &partitionAlgorithmName);
  
    TiXmlElement* salomeFileElement = meshBuilderElement->FirstChildElement("SalomeFile");

    assert(salomeFileElement);
    ParseString(salomeFileElement, "fileName", &salomeFileName);
  }

  std::string partitionAlgorithmName;
  std::string salomeFileName;
};