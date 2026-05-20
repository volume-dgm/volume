#pragma once

#include "../../../../3rdparty/tinyxml/tinyxml.h"
#include "../../../../3rdparty/tinyxml/tinystr.h"
#include <cstddef>
#include <limits>
#include <sstream>
#include <filesystem>
#include "ParserUtil.h"
#include "MeshSettingsGenericDTParser.h"
#include "ScheduleSettingsParser.h"
#include "SnapshotSettingsParser.h"
#include "SolverSettingsParser.h"
#include "TaskGenericDTSettingsParser.h"
#include "MeshBuilderSettingsParser.h"
#include "ResultCombinerSettingsParser.h"
#include "LambdaSettingsParser.hpp"

namespace fs = std::filesystem;

struct BasicSettings
{
  virtual ~BasicSettings() = default;
  std::string taskName;
  fs::path configFilePath;
  fs::path outDirPath;
  int configDimsCount; 
  int configPolyOrder;

  virtual void Parse(const std::string& fileName);

protected:
  void SetupPaths(const std::string& inputString) 
  {
    if (inputString.find(".xml") == std::string::npos)
    {
      this->taskName = inputString;
      this->configFilePath = "config/" + inputString + ".xml";
    } else 
    {
      this->taskName = fs::path(inputString).stem().string();
      this->configFilePath = inputString;
    }
    this->outDirPath = "out/" + this->taskName;
  }

  static TiXmlElement* GetSettingXmlElement(TiXmlDocument& taskFile, const std::string& fileName)
  {
    bool taskLoadOkay = taskFile.LoadFile(fileName);

    if (!taskLoadOkay)
    {
      std::cerr << "Loading " << fileName << " fails with following error: " <<
        std::string(taskFile.ErrorDesc()) <<
        " in row " << taskFile.ErrorRow() << std::endl;
      throw;
    }

    TiXmlElement* settingsElement = taskFile.FirstChildElement("Settings");
    if (!settingsElement)
    {
      std::cerr << "There is no Settings element";
      throw;
    }
    return settingsElement;
  }
};


template <typename Space>
struct Settings : public BasicSettings
{
  SPACE_TYPEDEFS

  LambdaSettings          <Space>   lambdaParser;
  MeshSettings            <Space>   mesh{&lambdaParser};
  ScheduleSettings        <Space>   schedule;
  TaskSettings            <Space>   task{&lambdaParser};
  SolverSettings          <Space>   solver;
  DetectorsSettings       <Space>   detectors;
  MeshBuilderSettings     <Space>   meshBuilder;
  ResultCombinerSettings  <Space>   resultCombiner;
  std::vector< SnapshotSettings<Space> > snapshots;

  void Parse(const std::string& fileName) override;

private:
  void ParseSettingsFile();
};

void BasicSettings::Parse(const std::string& fileName)
{
  TiXmlDocument taskFile;
  TiXmlElement* settingsElement = GetSettingXmlElement(taskFile, fileName);
  ParseString(settingsElement, "fileName", &taskName);
  ParseUnsigned(settingsElement, "dimsCount", &configDimsCount);
  ParseUnsigned(settingsElement, "polynomialsOrder", &configPolyOrder);

  if (taskName.empty())
  {
    std::cerr << "You should choose settings file to compute\n";
    throw;
  }
}

template<typename Space>
void Settings<Space>::Parse(const std::string& fileName)
{
  SetupPaths(fileName);
  Settings<Space>::ParseSettingsFile();
}

template<typename Space>
void Settings<Space>::ParseSettingsFile()
{
  TiXmlDocument settingsFile;
  TiXmlElement* settingsElement = GetSettingXmlElement(settingsFile, configFilePath);

  if(settingsElement->QueryIntAttribute("dimsCount", &configDimsCount) != TIXML_SUCCESS)
  {
    std::cerr << "dimsCount not provided in Settings element in config" << std::endl;
    std::cerr << "exiting just in case (don't know and don't wan to check whether it will break anything)" << std::endl;
    throw;
  }

  TiXmlElement* lambdaParserElement = settingsElement->FirstChildElement("LambdaParser");
  if (lambdaParserElement)
  {
    lambdaParser.Parse(lambdaParserElement);
  } 

  TiXmlElement* meshInfoElement = settingsElement->FirstChildElement("Mesh");
  if (meshInfoElement)
  {
    mesh.Parse(meshInfoElement);
  } else
  {
    std::cout << "There is no Mesh section\n";
  }

  TiXmlElement* snapshotElement;
  for (IndexType snapshotIndex = 0; ; ++snapshotIndex)
  {
    if (snapshotIndex == 0)
    {
      snapshotElement = settingsElement->FirstChildElement("Snapshot");
    } else
    {
      snapshotElement = snapshotElement->NextSiblingElement("Snapshot");
    }

    SnapshotSettings<Space> snapshot;
    if (!snapshotElement)
    {
      snapshot.data      .used = false;
      snapshot.mesh      .used = false;
      snapshot.contacts  .used = false;
      if (snapshotIndex == 0)
      {
        std::cout << "There is no Snapshot section\n";
      }
      break;
    }else
    {
      snapshot.Parse(snapshotElement);
    }
    snapshots.push_back(snapshot);
  }

  TiXmlElement* scheduleElement = settingsElement->FirstChildElement("Schedule");
  if (scheduleElement)
  {
    schedule.Parse(scheduleElement);
  } else
  {
    std::cout << "There is no Schedule element\n";
  }
    
  TiXmlElement* taskElement = settingsElement->FirstChildElement("Task");
  if (taskElement)
  {
    task.Parse(taskElement);
  } else
  {
    std::cout << "There is no Task element\n";
  }

  TiXmlElement* solverElement = settingsElement->FirstChildElement("Solver");
  if (solverElement)
  {
    solver.Parse(solverElement);
  } else
  {
    std::cout << "There is no Solver section\n";
  }

  TiXmlElement* meshBuilderElement = settingsElement->FirstChildElement("MeshBuilder");
  if (meshBuilderElement)
  {
    meshBuilder.Parse(meshBuilderElement);
  } else
  {
    std::cout << "There is no MeshBuilder section\n";
  }

  TiXmlElement* resultCombinerElement = settingsElement->FirstChildElement("ResultCombiner");
  if (resultCombinerElement)
  {
    resultCombiner.Parse(resultCombinerElement);
  } else
  {
    std::cout << "There is no ResultCombiner section\n";
  }

  detectors.Parse(settingsElement);
}
