#pragma once

#include "BasicVtkWriter.h"
#include "../../Maths/Spaces.h"
#include "../../Task/GeomMesh/GeomMesh/GeomMesh.h"
#include "../../Task/ElasticSystem/ElasticSystem.h"

using std::string;
using std::fstream;

template<typename Space, typename FunctionSpace>
class ContactVtkWriter;

template<typename FunctionSpace>
class ContactVtkWriter<Space2, FunctionSpace>: public BasicVtkWriter<Space2>
{
public:
  SPACE2_TYPEDEFS

  typedef typename GeomMesh<Space2>::Node              Node;
  typedef typename GeomMesh<Space2>::Cell              Cell;
  typedef typename GeomMesh<Space2>::EdgePairIndices   EdgePairIndices;
  typedef typename GeomMesh<Space2>::BoundaryEdge      BoundaryEdge;

  using BasicVtkWriter<Space2>::TypeName;

  void Write(const std::string& fileName,
              const std::vector<Node>& nodes, const std::vector<Cell>& cells, const std::vector<bool>& isCellBroken,
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

  struct OutputData
  {
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<IndexType> nodeData;
    std::vector<IndexType> cellData;
  };

  OutputData ConstructOutputData(const std::vector<Node>& nodes, 
  const std::vector<Cell>& cells, const std::vector<bool>& isCellBroken,
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

#include "ContactVtkWriter2.inl"

template <typename FunctionSpace>
class ContactVtkWriter<Space3, FunctionSpace>: public BasicVtkWriter<Space2>
{
public:
  SPACE3_TYPEDEFS
};

#include "ContactVtkWriter3.inl"
