#pragma once

#include "ParserUtil.h"
#include "../../../Maths/Spaces.h"

template <typename Space>
struct ResultCombinerSettings
{
  SPACE_TYPEDEFS

  ResultCombinerSettings():
    snapshotsCount(1)
  {}

  void Parse(TiXmlElement *resultCombinerElement)
  {
    snapshotsCount = 0;
    TiXmlElement* snapshotsElement = resultCombinerElement->FirstChildElement("Snapshots");
    if(snapshotsElement)
    {
      ParseUnsigned(snapshotsElement, "count", &snapshotsCount);
    }

    TiXmlElement* seismogramsElement = resultCombinerElement->FirstChildElement("Seismograms");
    if(seismogramsElement)
    {
      ParseString(seismogramsElement, "detectorsFileName", &detectorsFileName);
      ParseString(seismogramsElement, "refDetectorsFileName", &refDetectorsFileName);

      ParseString(seismogramsElement, "detectorsLocationsFile", &detectorsLocationsFile);

      ParseString(seismogramsElement, "velocityScvName", &velocityScvName);
      ParseString(seismogramsElement, "pressureScvName", &pressureScvName);

      ParseString(seismogramsElement, "velocityCoordSegyName", &velocityCoordSegyName);
    }
  }

  int snapshotsCount;

  std::string detectorsFileName;
  std::string refDetectorsFileName;

  std::string detectorsLocationsFile;

  std::string velocityScvName;
  std::string pressureScvName;

  std::string velocityCoordSegyName;
};