#include <string>

template<typename Space>
struct SolverSettings
{
  SPACE_TYPEDEFS

  SolverSettings(): tolerance(0.0), maxScale(1.0), maxTimeStep(1.0),
    tensionErrorMult(1.0), velocityErrorMult(1.0), positionErrorMult(1.0),
    velocityDimensionlessMult(1.0), tensionDimensionlessMult(1.0),
    allowMovement(false), allowPlasticity(false), allowDestruction(false),
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
  bool allowDestruction;

  std::string integrator;
  std::string precision;
  IndexType polynomialsOrder;
  IndexType hierarchyLevelsCount;
  Scalar damping;

  void Parse(TiXmlElement* solverElement);
};

template<typename GeomSpace>
void SolverSettings<GeomSpace>::Parse(TiXmlElement *solverElement)
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
      
  if (ParseBool(solverElement, "allowDestruction", &allowDestruction) != TIXML_SUCCESS)
  {
    std::cerr << "There is no allowDestructions attribute";
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
}
