#pragma once

#include "../../../../3rdparty/tinyxml/tinyxml.h"
#include "../../../../3rdparty/tinyxml/tinystr.h"
#include "ParserUtil.h"

#include <string>
#include <limits>

template<typename Space>
struct SolverSettings
{
  SPACE_TYPEDEFS

  SolverSettings(): tolerance(0.0), maxScale(1.0), maxTimeStep(1.0),
    tensionErrorMult(1.0), velocityErrorMult(1.0), positionErrorMult(1.0),
    velocityDimensionlessMult(1.0), tensionDimensionlessMult(1.0),
    allowMovement(false), allowPlasticity(false), 
    allowDiscreteDestruction(false), allowContinuousDestruction(false),
    integrator("Euler"), precision("double"), polynomialsOrder(1), hierarchyLevelsCount(1), damping(0.0)
  {
  }

  Scalar tolerance;
  Scalar maxScale;
  Scalar maxTimeStep;

  Scalar tensionErrorMult;
  Scalar velocityErrorMult;
  Scalar positionErrorMult;

  Scalar velocityDimensionlessMult;
  Scalar tensionDimensionlessMult;

  bool allowMovement;
  bool allowPlasticity;
  bool allowDiscreteDestruction;
  bool allowContinuousDestruction;

  std::string integrator;
  std::string precision;
  IndexType polynomialsOrder;
  IndexType hierarchyLevelsCount;
  Scalar damping;

  struct Erosion
  {
    Erosion() : cellAspectRatio(std::numeric_limits<Scalar>::max() / 2),
      minHeightRatio(0),
      maxPlasticDeformation(1.0),
      rhoReduction(std::numeric_limits<Scalar>::max() / 2)
    {
    }
    Scalar cellAspectRatio;
    Scalar minHeightRatio;
    Scalar maxPlasticDeformation;
    Scalar rhoReduction;
  } erosion;

  void Parse(TiXmlElement* solverElement);
};

template<typename Space>
void SolverSettings<Space>::Parse(TiXmlElement *solverElement)
{
  if (ParseString(solverElement, "integrator", &integrator) != TIXML_SUCCESS)
  {
    std::cerr << "There is no integrator attribute";
  }
  
  if (ParseUnsigned(solverElement, "hierarchyLevelsCount", &hierarchyLevelsCount) != TIXML_SUCCESS)
  {
    std::cerr << "There is no hierarchyLevelsCount attribute";
  }

  if (ParseString(solverElement, "precision", &precision) != TIXML_SUCCESS)
  {
    std::cerr << "There is no precision attribute";
  }

  if (ParseUnsigned(solverElement, "polynomialsOrder", &polynomialsOrder) != TIXML_SUCCESS)
  {
    std::cerr << "There is no polynomialsOrder attribute";
  }

  if (ParseScalar(solverElement, "tolerance", &tolerance) != TIXML_SUCCESS)
  {
    std::cerr << "There is no tolerance attribute";
  }

  if (ParseScalar(solverElement, "maxScale", &maxScale) != TIXML_SUCCESS)
  {
    std::cerr << "There is no maxScale attribute";
  }

  if (ParseScalar(solverElement, "maxTimeStep", &maxTimeStep) != TIXML_SUCCESS)
  {
    std::cerr << "There is no maxTimeStep attribute";
  }

  if (ParseBool(solverElement, "allowMovement", &allowMovement) != TIXML_SUCCESS)
  {
    std::cerr << "There is no allowMovement attribute";
  }

  if (ParseBool(solverElement, "allowPlasticity", &allowPlasticity) != TIXML_SUCCESS)
  {
    std::cerr << "There is no allowPlasticity attribute";
  }
      
  if (ParseBool(solverElement, "allowContinuousDestruction", &allowContinuousDestruction) != TIXML_SUCCESS)
  {
    std::cerr << "There is no allowContinuousDestruction attribute";
  }

  if (ParseBool(solverElement, "allowDiscreteDestruction", &allowDiscreteDestruction) != TIXML_SUCCESS)
  {
    std::cerr << "There is no allowDiscreteDestruction attribute";
  }

  if (ParseScalar(solverElement, "tensionErrorMult", &tensionErrorMult) != TIXML_SUCCESS)
  {
    std::cerr << "There is no tensionErrorMult attribute";
  }

  if (ParseScalar(solverElement, "velocityErrorMult", &velocityErrorMult) != TIXML_SUCCESS)
  {
    std::cerr << "There is no velocityErrorMult attribute";
  }

  if (ParseScalar(solverElement, "positionErrorMult", &positionErrorMult) != TIXML_SUCCESS && allowMovement)
  {
    std::cerr << "allowMovement is true, however positionErrorMult is not specified";
  }

  ParseScalar(solverElement, "velocityDimensionlessMult", &velocityDimensionlessMult);  
  ParseScalar(solverElement, "tensionDimensionlessMult", &tensionDimensionlessMult);
  ParseScalar(solverElement, "damping", &damping);

  TiXmlElement* erosionElement = solverElement->FirstChildElement("Erosion");
  if (erosionElement)
  {
    ParseScalar(erosionElement, "cellAspectRatio", &erosion.cellAspectRatio);
    ParseScalar(erosionElement, "minHeightRatio",  &erosion.minHeightRatio);
    ParseScalar(erosionElement, "maxPlasticDeformation", &erosion.maxPlasticDeformation);
    ParseScalar(erosionElement, "rhoReduction", &erosion.rhoReduction);
  }
}
