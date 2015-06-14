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
    TiXmlElement* snapshotsElement = resultCombinerElement->FirstChildElement("Snapshots");
    assert(snapshotsElement);

    ParseUnsigned(snapshotsElement, "count", &snapshotsCount);
  }

  int snapshotsCount;
};