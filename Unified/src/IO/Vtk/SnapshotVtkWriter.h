#pragma once

#include "BasicVtkWriter.h"
#include <vector>
#include "../../Utils/Base64.h"
#include <zlib.h>
#include <zconf.h>

using std::string;
using std::fstream;
using std::vector;

template <typename ElasticSpace>
class SnapshotVtkWriter: public BasicVtkWriter<ElasticSpace>
{
  public:
    typedef typename ElasticSpace::SpaceType Space;
    SPACE_TYPEDEFS

    typedef typename ElasticSpace::Elastic   Elastic;

    void Write(const string& fileName,
               Elastic* data,
               const Vector& origin, const Vector& spacing,
               const AABB& wholeBox, const AABB& box,
               bool writeVelocity, bool writeTension, bool writeWaves, bool asScalars = false);

    void WriteSnapsotCollectionHeader(std::fstream& file)
    {
      file << "<?xml version=\"1.0\"?>\n";
      file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      file << "<Collection>\n";
    }

    void WriteSnapsotCollectionLastTags(std::fstream& file)
    {
      file << "</Collection>";
      file << "</VTKFile>";
    }

  private:
    void WriteCompressedData(fstream& file, unsigned char* src, size_t srcLen);

    string FloatTypeName() const
    {
      if (sizeof(Scalar) == 4)
      {
        return "\"Float32\"";
      } else
      if (sizeof(Scalar) == 8)
      {
        return "\"Float64\"";
      } else
      {
        throw 42;
      }
    }

    std::vector<Scalar> rawData;
    std::vector<unsigned char> compressedData;
};


template <typename ElasticSpace>
void SnapshotVtkWriter<ElasticSpace>::Write(const string& fileName,
                                            Elastic* data,    
                                            const Vector& origin, const Vector& spacing,
                                            const AABB& wholeBox, const AABB& box,
                                            bool writeVelocity, bool writeTension, bool writeWaves, bool asScalars)
{
  fstream file(fileName.c_str(), fstream::out | fstream::binary);

  if (file.fail()) return;

  IntAABB iaabb;
  ::ComputeSnapshotSize<Space>(origin, spacing, wholeBox, &iaabb);

  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
  
  file << "<ImageData WholeExtent=\"";
  // it`s impossible to create 2d snapshot directly. so we have to set third component in both 2d/3d cases.  
  for (IndexType dimIndex = 0; dimIndex < 3; ++dimIndex)
  {
    file << iaabb.boxPoint1.Get(dimIndex) << " " << iaabb.boxPoint2.Get(dimIndex) << " ";
  }
  file << "\"";

  file << " Origin=\"";
  for (IndexType dimIndex = 0; dimIndex < 3; ++dimIndex)
  {
    file << origin.Get(dimIndex)  << " ";
  }
  file << "\"";

  file << " Spacing=\"";
  for (IndexType dimIndex = 0; dimIndex < 3; ++dimIndex)
  {
    file << spacing.Get(dimIndex) << " ";
  }
  file << "\">\n";

  IndexType snapshotSize = ::ComputeSnapshotSize<Space>(origin, spacing, box, &iaabb);

  IndexType tensorElementsCount = Space::Dimension * (Space::Dimension + 1) / 2;
  IndexType vectorElementsCount = 3; //we have to write 3 components regardless of space used

  IndexType maxElementsCount = 0;

  if(writeVelocity)
  {
    maxElementsCount = std::max(maxElementsCount, vectorElementsCount);
  }
  if(writeTension)
  {
    maxElementsCount = std::max(maxElementsCount, tensorElementsCount);
  }
  if(writeWaves)
  {
    maxElementsCount = std::max(maxElementsCount, IndexType(1));
  }
  assert(maxElementsCount != 0);

  if (rawData.size() < snapshotSize * maxElementsCount) rawData.resize(snapshotSize * maxElementsCount);

  file << "<Piece Extent=\"";
  for (IndexType dimIndex = 0; dimIndex < 3; ++dimIndex)
  {
    file << iaabb.boxPoint1.Get(dimIndex) << " " << iaabb.boxPoint2.Get(dimIndex) << " ";
  }
  file << "\">\n";

  std::string dataTypeString = "<PointData";

  //this section is "active" components, not all components present. actually both "Scalars=..." and "Vectors=..." can be skipped safely
  if(writeTension)
  {
    dataTypeString += " Scalars=\"Stress\"";
  }else
  {
    if(asScalars && writeVelocity)
    {
      dataTypeString += " Scalars=\"Velocity\"";
    }else
    {
      if(writeWaves)
      {
        dataTypeString += " Scalars=\"abs(pWave)\"";
      }
    }
  }
  if(writeVelocity && !asScalars)
  {
    dataTypeString += " Vectors=\"Velocity\"";
  }
  dataTypeString += ">\n";
  file << dataTypeString;

  if(writeVelocity)
  {
    file << "<DataArray NumberOfComponents=\"3\" type=" << FloatTypeName() << " Name=\"Velocity\" format=\"binary\">";
    for (IndexType elasticIndex = 0; elasticIndex < snapshotSize; ++elasticIndex)
    {
      for (IndexType dimIndex = 0; dimIndex < vectorElementsCount; ++dimIndex)
      {
        rawData[3 * elasticIndex + dimIndex] = data[elasticIndex].GetVelocity().Get(dimIndex);
      }
    }
    WriteCompressedData(file, (unsigned char*)(&rawData[0]), snapshotSize * vectorElementsCount * sizeof(Scalar));
    file << "</DataArray>\n";
  }

  if(writeTension)
  {
    file << "<DataArray NumberOfComponents=\"" << tensorElementsCount
         << "\" type=" << FloatTypeName() << " Name=\"Stress\" format=\"binary\">";
    for (IndexType elasticIndex = 0; elasticIndex < snapshotSize; ++elasticIndex)
    {
      const Tensor& tension = data[elasticIndex].GetTension();
      for (IndexType tensorElementIndex = 0; tensorElementIndex < tensorElementsCount; ++tensorElementIndex)
      {
        rawData[tensorElementsCount * elasticIndex + tensorElementIndex] = tension[tensorElementIndex];
      }
    }
    WriteCompressedData(file, (unsigned char*)(&rawData[0]), snapshotSize * tensorElementsCount * sizeof(Scalar));
    file << "</DataArray>\n"; 
  }

  if(writeWaves)
  {
    //pWave
    file << "<DataArray NumberOfComponents=\"1\" type=" << FloatTypeName() << " Name=\"abs(pWave)\" format=\"binary\">";
    for (IndexType elasticIndex = 0; elasticIndex < snapshotSize; ++elasticIndex)
    {
      const Tensor& tension = data[elasticIndex].GetTension();
      const Vector& velocity = data[elasticIndex].GetVelocity();
      Vector energyFlow = tension(velocity);
      Vector energyGradient = energyFlow.GetNorm();

      Scalar pWave = energyGradient * velocity; //p-wave intensity is v * gradient
      rawData[elasticIndex] = fabs(pWave);
    }
    WriteCompressedData(file, (unsigned char*)(&rawData[0]), snapshotSize * sizeof(Scalar));
    file << "</DataArray>\n";

    //sWave
    file << "<DataArray NumberOfComponents=\"1\" type=" << FloatTypeName() << " Name=\"sWave\" format=\"binary\">";
    for (IndexType elasticIndex = 0; elasticIndex < snapshotSize; ++elasticIndex)
    {
      const Tensor& tension = data[elasticIndex].GetTension();
      const Vector& velocity = data[elasticIndex].GetVelocity();
      Vector energyFlow = tension(velocity);
      Vector energyGradient = energyFlow.GetNorm();

      Scalar sWave = sqrt((energyGradient ^ velocity) * (energyGradient ^ velocity)); //s-wave intensity is v ^ gradient

      rawData[elasticIndex] = sWave;
    }
    WriteCompressedData(file, (unsigned char*)(&rawData[0]), snapshotSize * sizeof(Scalar));
    file << "</DataArray>\n";

    //abs(pWave)
    file << "<DataArray NumberOfComponents=\"1\" type=" << FloatTypeName() << " Name=\"pWave\" format=\"binary\">";
    for (IndexType elasticIndex = 0; elasticIndex < snapshotSize; ++elasticIndex)
    {
      const Tensor& tension = data[elasticIndex].GetTension();
      const Vector& velocity = data[elasticIndex].GetVelocity();
      Vector energyFlow = tension(velocity);
      Vector energyGradient = energyFlow.GetNorm();

      Scalar pWave = energyGradient * velocity; //p-wave intensity is v * gradient
      rawData[elasticIndex] = pWave;
    }
    WriteCompressedData(file, (unsigned char*)(&rawData[0]), snapshotSize * sizeof(Scalar));
    file << "</DataArray>\n";
  }

  file << "</PointData>\n";

  file << "</Piece>\n";
  file << "</ImageData>\n";
  file << "</VTKFile>\n";
}

template <typename ElasticSpace>
void SnapshotVtkWriter<ElasticSpace>::WriteCompressedData(fstream& file, unsigned char* source, size_t sourceLen)
{
  uLongf compressedDataLen = ::compressBound(sourceLen);
  if (compressedData.size() < compressedDataLen) compressedData.resize(compressedDataLen);

  ::compress(&compressedData[0], &compressedDataLen, source, sourceLen);

  uint32_t prefix[4];
  prefix[0] = 1;
  prefix[1] = (uint32_t)sourceLen;
  prefix[2] = (uint32_t)sourceLen;
  prefix[3] = (uint32_t)compressedDataLen;

  size_t outputSize;
  char* output = ::base64_encode((const unsigned char*)prefix, 4 * sizeof(uint32_t), &outputSize);
  if(output == NULL)
  {
    printf("Fatal error: base64_encode returned NULL trying to encode %lu bytes", 4 * sizeof(uint32_t));
    return;
  }
  file.write(output, outputSize);

  output = ::base64_encode((const unsigned char*)(&compressedData[0]), compressedDataLen, &outputSize);
  if(output == NULL)
  {
    printf("Fatal error: base64_encode returned NULL trying to encode %lu bytes", compressedDataLen);
    return;
  }
  file.write(output, outputSize);
}
