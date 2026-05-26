#pragma once
#include "ParserUtil.h"

#include "LambdaParser.hpp"
#include "tinyxml.h"
#include <algorithm>
#include <optional>
#include <string>

template<typename Space>
struct LambdaSettings
{
  SPACE_TYPEDEFS

  using FuncType = std::function<GenericDataType(std::shared_ptr<Thunk>)>;
  using FuncTypeAppl = std::pair<std::shared_ptr<Thunk>, std::shared_ptr<Thunk>>;

  XMLLambdaParserInstance parser;
  LambdaSettings() 
  {
    parser = XMLLambdaParserInstance();
    //parser.debug_log_eval = true;
    parser.AssignGlobal("lambda global . 'Vector2'", "lambda global x y . insert (insert table 'x' x) 'y' y");
    parser.AssignGlobal("lambda global . 'Vector3'", "lambda global x y z . insert (insert (insert table 'x' x) 'y' y) 'z' z");
  }

  std::shared_ptr<Thunk> vectorX = std::make_shared<Thunk>(parser.Expression("lambda global vec . get vec 'x'"));
  std::shared_ptr<Thunk> vectorY = std::make_shared<Thunk>(parser.Expression("lambda global vec . get vec 'y'"));
  std::shared_ptr<Thunk> vectorZ = std::make_shared<Thunk>(parser.Expression("lambda global vec . get vec 'z'"));

  template<typename First, typename... Vars>
  GenericDataType ApplyVars(GenericDataType* lambdaExpr, First var, Vars... vars);
  GenericDataType ApplyVars(GenericDataType* lambdaExpr)
  {
    if (lambdaExpr->eval().myType == GenericDataType::FUNCTION)
    {
      std::cerr << "error: not enough data passed to lambda (or some lambda has way too many variables)\n";
      throw;
    }
    return lambdaExpr->eval();
  }

  Vector LambdaToVector(GenericDataType* vectorExpr, const Vector& meshPoint=Vector::zero());
  Scalar LambdaToScalar(GenericDataType* scalarExpr, const Vector& meshPoint=Vector::zero());

  void LoadGlobals(const std::string& fileName);
  void Parse(TiXmlElement* lambdaParserElement);
  int ParseToGenericDT(TiXmlElement* element, const std::string& name, GenericDataType* value);
  int ParseToGenericDT(TiXmlElement* element, const std::string& name, std::optional<GenericDataType>* value);

  void toCorrectNumberString(std::string& num);
};

template<typename Space>
template<typename First, typename... Vars>
GenericDataType LambdaSettings<Space>::ApplyVars(GenericDataType* lambdaExpr, First var, Vars... vars)
{
  if (lambdaExpr->eval().myType == GenericDataType::NUMBER) return lambdaExpr->eval();

  GenericDataType applyFirst = GenericDataType();
  applyFirst.myType = GenericDataType::FUNCTIONAPPL;
  applyFirst.funcAppl_val = std::make_shared<FuncTypeAppl>(std::make_shared<Thunk>(*lambdaExpr), std::make_shared<Thunk>(GenericDataType(var)));
  return LambdaSettings<Space>::ApplyVars(&applyFirst, vars...);
}

template<typename Space>
typename Space::Scalar LambdaSettings<Space>::LambdaToScalar(GenericDataType* scalarExpr, const typename Space::Vector& meshPoint)
{
  GenericDataType applyCoords = LambdaSettings<Space>::ApplyVars(scalarExpr, meshPoint.Get(0), meshPoint.Get(1), meshPoint.Get(2));
  return typename Space::Scalar(applyCoords.num_val);
}

template<>
typename Space2::Vector LambdaSettings<Space2>::LambdaToVector(GenericDataType* vectorExpr, const Space2::Vector& meshPoint)
{
  GenericDataType applyCoords = LambdaSettings<Space2>::ApplyVars(vectorExpr, meshPoint.Get(0), meshPoint.Get(1));

  GenericDataType getX = GenericDataType();
  getX.myType = GenericDataType::FUNCTIONAPPL;
  getX.funcAppl_val = std::make_shared<FuncTypeAppl>(vectorX, std::make_shared<Thunk>(applyCoords));
  GenericDataType getY = GenericDataType();
  getY.myType = GenericDataType::FUNCTIONAPPL;
  getY.funcAppl_val = std::make_shared<FuncTypeAppl>(vectorY, std::make_shared<Thunk>(applyCoords));
  return typename Space2::Vector(getX.eval().num_val, getY.eval().num_val);
}

template<>
typename Space3::Vector LambdaSettings<Space3>::LambdaToVector(GenericDataType* vectorExpr, const Space3::Vector& meshPoint)
{
  GenericDataType applyCoords = LambdaSettings<Space3>::ApplyVars(vectorExpr, meshPoint.Get(0), meshPoint.Get(1), meshPoint.Get(2));

  GenericDataType getX = GenericDataType();
  getX.myType = GenericDataType::FUNCTIONAPPL;
  getX.funcAppl_val = std::make_shared<FuncTypeAppl>(vectorX, std::make_shared<Thunk>(applyCoords));
  GenericDataType getY = GenericDataType();
  getY.myType = GenericDataType::FUNCTIONAPPL;
  getY.funcAppl_val = std::make_shared<FuncTypeAppl>(vectorY, std::make_shared<Thunk>(applyCoords));
  GenericDataType getZ = GenericDataType();
  getZ.myType = GenericDataType::FUNCTIONAPPL;
  getZ.funcAppl_val = std::make_shared<FuncTypeAppl>(vectorZ, std::make_shared<Thunk>(applyCoords));
  return typename Space3::Vector(getX.eval().num_val, getY.eval().num_val, getZ.eval().num_val);
}

template<typename Space>
void LambdaSettings<Space>::LoadGlobals(const std::string& fileName)
{
  TiXmlDocument globalsFile;
  bool fileLoadStatus = globalsFile.LoadFile(fileName);  

  if (!fileLoadStatus)
  {
    std::cerr << "Loading " << fileName << " fails: " << std::string(globalsFile.ErrorDesc()) << 
                 " in row " << globalsFile.ErrorRow() << std::endl;
    throw;
  }

  TiXmlElement* lambdasElement = globalsFile.FirstChildElement("LambdaParser");
  if (!lambdasElement)
  {
    std::cerr << "There is no LambdaParser element in " << fileName << std::endl;
    throw;
  }

  Parse(lambdasElement);
}

template<typename Space>
void LambdaSettings<Space>::Parse(TiXmlElement* lambdaParserElement)
{
  std::string tmpName;
  if (ParseString(lambdaParserElement, "LoadFromFile", &tmpName) == TIXML_SUCCESS)
  {
    LoadGlobals(tmpName);
  }

  TiXmlElement* assignGlobalElement = lambdaParserElement->FirstChildElement("AssignGlobal");
  while(assignGlobalElement)
  {
    std::string globalName, globalExpr;
    ParseString(assignGlobalElement, "name", &globalName); 
    std::string formatName = "lambda global . '" + globalName + "'";
    ParseString(assignGlobalElement, "value", &globalExpr);
    if (globalExpr.find("lambda") != std::string::npos) //treat as lambda expression and assign as is
    {
      parser.AssignGlobal(formatName, globalExpr);
    } else 
    {
      toCorrectNumberString(globalExpr);
      parser.AssignGlobal(formatName, "lambda global . " + globalExpr);
    }

    assignGlobalElement = assignGlobalElement->NextSiblingElement("AssignGlobal");
  }

  TiXmlElement* globalsFileElement = lambdaParserElement->FirstChildElement("LoadGlobals");
  while(globalsFileElement)
  {
    ParseString(globalsFileElement, "fileName", &tmpName);
    LoadGlobals(tmpName);

    globalsFileElement = globalsFileElement->NextSiblingElement("LoadGlobals");
  }
}

template<typename Space>
int LambdaSettings<Space>::ParseToGenericDT(TiXmlElement* element, const std::string& name, GenericDataType* value)
{
  std::string parsedString;
  int returnCode = ParseString(element, name, &parsedString);
  if (returnCode == TIXML_NO_ATTRIBUTE) return returnCode;
  if (parsedString.find("lambda") != std::string::npos)
  {
    *value = parser.Expression(parsedString);
  } else 
  {
    toCorrectNumberString(parsedString);
    *value = Parser::GenericDataTypeFromString(parsedString, parser.funcs);
  }
  return TIXML_SUCCESS;
}

template<typename Space>
int LambdaSettings<Space>::ParseToGenericDT(TiXmlElement* element, const std::string& name, std::optional<GenericDataType>* value)
{
  std::string parsedString;
  int returnCode = ParseString(element, name, &parsedString);
  if (returnCode == TIXML_NO_ATTRIBUTE) return returnCode;
  if (parsedString.find("lambda") != std::string::npos)
  {
    *value = std::make_optional<GenericDataType>(parser.Expression(parsedString));
  } else 
  {
    toCorrectNumberString(parsedString);
    *value = std::make_optional<GenericDataType>(Parser::GenericDataTypeFromString(parsedString, parser.funcs));
  }
  return TIXML_SUCCESS;
}

template<typename Space>
void LambdaSettings<Space>::toCorrectNumberString(std::string& num)
{
  if (num.find(',') != std::string::npos) 
  { 
    std::replace(num.begin(), num.end(), ',', '.');
  }
  num = std::to_string(std::stod(num));
  std::replace(num.begin(), num.end(), '.', ',');
}