#pragma once

#include "BasicVtkWriter.h"
#include "../../Maths/Spaces.h"
#include "../../Task/GeomMesh/GeomMesh/GeomMesh.h"
#include "../../Task/ElasticSystem/ElasticSystem.h"

using std::string;
using std::fstream;

template<typename Space, typename FunctionSpace>
class ContactVtkWriterBase: public BasicVtkWriter<Space>
{
public:
  SPACE_TYPEDEFS

  typedef typename GeomMesh<Space>::Node              Node;
  typedef typename GeomMesh<Space>::Cell              Cell;

  using BasicVtkWriter<Space>::TypeName;

  struct OutputData
  {
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<IndexType> nodeData;
    struct CellData
    {
      CellData(IndexType isCellBroken, Scalar plasticDeforms = 0) :
        isCellBroken(isCellBroken), plasticDeforms(plasticDeforms)
      {}

      CellData()
      {}

      IndexType isCellBroken;
      Scalar plasticDeforms;
    };
    std::vector<CellData> cellData;
  };

  void Write(const std::string& fileName,
    ElasticVolumeMesh<Space, FunctionSpace>* mesh,
    const ElasticSystem<Space>& system);

protected:
  void SaveToFile(const std::string& fileName, const OutputData& outputData) const;

private:
  OutputData ConstructOutputData(ElasticVolumeMesh<Space, FunctionSpace>* mesh,
    const ElasticSystem<Space>& system);
};

template<typename Space, typename FunctionSpace>
class ContactVtkWriter;

template<typename FunctionSpace>
class ContactVtkWriter<Space2, FunctionSpace>: public ContactVtkWriterBase<Space2, FunctionSpace>
{
public:
  SPACE2_TYPEDEFS

  typedef typename ContactVtkWriterBase<Space2, FunctionSpace>::Node  Node;
  typedef typename ContactVtkWriterBase<Space2, FunctionSpace>::Cell  Cell;
  typedef typename GeomMesh<Space2>::EdgePairIndices   EdgePairIndices;
  typedef typename GeomMesh<Space2>::BoundaryEdge      BoundaryEdge;

  using ContactVtkWriterBase<Space2, FunctionSpace>::Write;
  using ContactVtkWriterBase<Space2, FunctionSpace>::SaveToFile;

  void Write(const std::string& fileName,
              const std::vector<Node>& nodes, const std::vector<Cell>& cells, const std::vector<bool>& isCellBroken,
              const std::vector<Scalar>& plasticDeforms,
              Scalar velocityDimensionlessMult,
              const std::vector<EdgePairIndices>&  contactEdges, 
              const std::vector<IndexType>&  contactEdgesCount, IndexType contactTypesCount, 
              const std::vector<bool>& isContactBroken,
              const std::vector<BoundaryEdge>&    boundaryEdges, 
              const std::vector<IndexType>& boundaryEdgesCount, IndexType boundaryTypesCount,
              const ElasticSystem<Space2>& system,
              bool drawContacts,
              bool drawRegularGlueContacts,
              bool drawBrokenContacts,
              // for debuging
              bool drawContactBindings = false);

private:

  OutputData ConstructOutputData(const std::vector<Node>& nodes, 
  const std::vector<Cell>& cells, const std::vector<bool>& isCellBroken,
  const std::vector<Scalar>& plasticDeforms,
  Scalar velocityDimensionlessMult,
  const std::vector<EdgePairIndices>&  contactEdges, 
  const std::vector<IndexType>&  contactEdgesCount, IndexType contactTypesCount, 
  const std::vector<bool>& isContactBroken,
  const std::vector<BoundaryEdge>&    boundaryEdges, 
  const std::vector<IndexType>& boundaryEdgesCount, IndexType boundaryTypesCount,
  const ElasticSystem<Space2>& system,
  bool drawContacts,
  bool drawRegularGlueContacts,
  bool drawBrokenContacts,
  // for debuging
  bool drawContactBindings = false);
};

template <typename FunctionSpace>
class ContactVtkWriter<Space3, FunctionSpace>: public ContactVtkWriterBase<Space3, FunctionSpace>
{
public:
  SPACE3_TYPEDEFS

  typedef typename ContactVtkWriterBase<Space3, FunctionSpace>::Node              Node;
  typedef typename ContactVtkWriterBase<Space3, FunctionSpace>::Cell              Cell;
  typedef typename GeomMesh<Space3>::FacePairIndices   FacePairIndices;
  typedef typename GeomMesh<Space3>::BoundaryFace      BoundaryFace;

  using ContactVtkWriterBase<Space3, FunctionSpace>::Write;
};

#include "ContactVtkWriterBase.inl"
#include "ContactVtkWriter2.inl"
#include "ContactVtkWriter3.inl"
