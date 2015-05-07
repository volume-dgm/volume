#include "../../Utils/Utils.h"
#include "../../Vtk/BasicVtkWriter.h"
#include "../../Vtk/SnapshotVtkWriter.h"
#include "../../Vtk/VtkReader.h"
#include "../Task/SettingsParser/SettingsParser.h"
#include "../ElasticSystem/ElasticSystem.h"
#include "SampleIO.h"

#include <fstream>
#include <vector>
#include <stdio.h>
#include <algorithm>

typedef Space2 DefaultSpace;

 #define SPACE_FROM_SETTINGS

template <typename Space>
class ResultCombinerTask
{
public:
  SPACE_TYPEDEFS

  struct Location
  {
    IndexType domainIndex;
    IndexType localIndex;
  };

  typedef ElasticSystem<Space>            ElasticSystemType;
  typedef typename ElasticSystemType::ElasticSpace ElasticSpace;
  typedef typename ElasticSpace::Elastic           Elastic;

  const static int dimsCount = ElasticSystemType::dimsCount;

  void Run()
  {
    resultCombinerSettings.Parse("resultcombiner.xml");
    settings.Parse("../Task/task.xml");

    IndexType domainsCount = settings.schedule.domainsCount;

    std::vector<Elastic> elastic, result;
    SnapshotVtkWriter<ElasticSpace> snapshotWriter;
    VtkReader<ElasticSpace> snapshotReader;

    for (IndexType snapshotIndex = 0; snapshotIndex < settings.snapshots.size(); ++snapshotIndex)
    {
      IndexType snapshotSize = ::ComputeSnapshotSize<Space>(settings.snapshots[snapshotIndex].data.origin,
        settings.snapshots[snapshotIndex].data.spacing,
        AABB(settings.snapshots[snapshotIndex].data.boxPoint1,
        settings.snapshots[snapshotIndex].data.boxPoint2));
      if (result.size() < snapshotSize) result.resize(snapshotSize);

      for (int stepIndex = 0; stepIndex < resultCombinerSettings.maxStepsCount; ++stepIndex)
      {
        #pragma omp parallel for
        for (int elasticIndex = 0; elasticIndex < int(snapshotSize); ++elasticIndex)
        {
          result[elasticIndex].SetZeroValues();
        }

        char stepString[256];
        sprintf(stepString, "%.6lu", stepIndex);
        bool isOk = true;

        for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
        {
          char domainString[256];
          sprintf(domainString, "%.3lu", domainIndex);

          if (settings.snapshots[snapshotIndex].data.used)
          {
            std::string dataFilename = settings.snapshots[snapshotIndex].data.filename;
            ReplaceSubstring(dataFilename, "<step>", std::string(stepString));
            ReplaceSubstring(dataFilename, "<domain>", std::string(domainString));
            std::string pathName = "../Task/" + dataFilename;
            bool res = false;

            if (EndWith(dataFilename, ".vti"))
            {
              res = snapshotReader.Read(pathName, &elastic);
            }
            if (!res)
            {
              isOk = false;
              break;
            }

            #pragma omp parallel for
            for (int elasticIndex = 0; elasticIndex < int(snapshotSize); ++elasticIndex)
            {
              for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
              {
                if (fabs(elastic[elasticIndex].values[valueIndex]) > fabs(result[elasticIndex].values[valueIndex]))
                {
                  result[elasticIndex].values[valueIndex] = elastic[elasticIndex].values[valueIndex];
                }
              }
            }
            if (remove(pathName.c_str()) != 0)
            {
              std::cout << "Can`t remove " << pathName << std::endl;
            }
          }
        }

        if (isOk)
        {
          std::string dataFilename = settings.snapshots[snapshotIndex].data.filename;
          ReplaceSubstring(dataFilename, "<step>", std::string(stepString));
          ReplaceSubstring(dataFilename, "<domain>", "");
          ReplaceSubstring(dataFilename, ".vtk", ".vti");
          std::cout << dataFilename << std::endl;

          Vector boxPoint1 = settings.snapshots[snapshotIndex].data.boxPoint1;
          Vector boxPoint2 = settings.snapshots[snapshotIndex].data.boxPoint2;
          Vector origin = settings.snapshots[snapshotIndex].data.origin;
          Vector spacing = settings.snapshots[snapshotIndex].data.spacing;

          snapshotWriter.Write(("../Task/" + dataFilename).c_str(),
            &result[0],
            origin, spacing,
            AABB(boxPoint1, boxPoint2),
            AABB(boxPoint1, boxPoint2),
            settings.snapshots[snapshotIndex].data.writeVelocity,
            settings.snapshots[snapshotIndex].data.writeTension,
            settings.snapshots[snapshotIndex].data.writeWaves);
        }
      }
    }
    CombineDetectors();
  }
private:
  Settings<Space> settings;
  ResultCombinerSettings<Space> resultCombinerSettings;

  void CombineDetectors()
  {
    SPACE_TYPEDEFS
      if (!settings.detectors.used) return;

    std::cout << "Combine detectors\n";

    std::string detectorsFileNamePattern = std::string("../Task/") +
      settings.detectors.filename;

    std::fstream detectorsLocationsFile("detectorsLocations.txt", std::ios_base::in);
    assert(detectorsLocationsFile.fail() == false);

    IndexType detectorsCount = 0;
    detectorsLocationsFile >> detectorsCount;

    std::vector<Location> detectorsLocations(detectorsCount, Location());

    for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
    {
      detectorsLocationsFile >>
        detectorsLocations[detectorIndex].domainIndex >>
        detectorsLocations[detectorIndex].localIndex;
    }

    std::vector<FILE*> files(detectorsCount, NULL);
    for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
    {
      char detectorString[9];
      sprintf(detectorString, "%.3lu", detectorsLocations[detectorIndex].localIndex);
      char domainString[9];
      sprintf(domainString, "%.3lu", detectorsLocations[detectorIndex].domainIndex);

      std::string fileName = detectorsFileNamePattern;
      ReplaceSubstring(fileName, "<detector>", std::string(detectorString));
      ReplaceSubstring(fileName, "<domain>", std::string(domainString));
      if ((files[detectorIndex] = fopen(fileName.c_str(), "r")) == NULL)
      {
        std::cerr << "Can`t open file: " << fileName << "\n";
      }
      assert(files[detectorIndex] != NULL);
    }

    std::fstream outputVelocity("../Task/out/velocity_detector.csv", std::ios_base::out);
    std::fstream outputPressure("../Task/out/pressure_detector.csv", std::ios_base::out);

    assert(outputVelocity.fail() == false);
    assert(outputPressure.fail() == false);

    outputVelocity << "Time;";
    outputPressure << "Time;";

    for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
    {
      WriteVelocityHandle<Space>(outputVelocity);
      outputPressure << "Vx;";
    }
    outputVelocity << "\n";
    outputPressure << "\n";

    bool eof;
    do
    {
      eof = false;
      for (IndexType detectorIndex = 0; detectorIndex < detectorsCount; ++detectorIndex)
      {
        Vector v = Vector::zero();
        Scalar time = 0;
        Tensor sigma;

        ReadSample<Space>(files[detectorIndex], time, v, sigma);

        if (feof(files[detectorIndex]))
        {
          eof = true;
          continue;
        }

        if (detectorIndex == 0)
        {
          outputVelocity << time << ";";
          outputPressure << time << ";";
        }
        WriteSample<Space>(outputVelocity, outputPressure, time, v, sigma);
      }
      outputVelocity << "\n";
      outputPressure << "\n";
    } while (!eof);

    outputVelocity.close();
    outputPressure.close();
  }
};

int main()
{
  BasicSettings basicSettings; //Space2 settings reader should read Space3 settings file just fine. probably.
  basicSettings.Parse("../Task/task.xml");

  #ifdef SPACE_FROM_SETTINGS
  {
    switch (basicSettings.configDimsCount)
    {
      case 2:
      {
        ResultCombinerTask<Space2> resultCombinerTask;
        resultCombinerTask.Run();
      } break;
      case 3:
      {
        ResultCombinerTask<Space3> resultCombinerTask;
        resultCombinerTask.Run();
      } break;
      default: 
        std::cerr << "Wrong dims count in config"; 
      break;
    }
  }
  #else
    ResultCombinerTask<DefaultSpace> resultCombinerTask;
    resultCombinerTask.Run();
  #endif

  return 0;
}
