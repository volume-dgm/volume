template <typename Space>
DistributedMeshIOCommon<Space>::DistributedMeshIOCommon(IndexType domainsCount): 
  MeshIO<Space>(), domainsCount(domainsCount), transitionInfos(domainsCount)
{
}

template <typename Space>
bool DistributedMeshIOCommon<Space>::operator==(const DistributedMeshIOCommon& other) const
{
  return MeshIO<Space>::operator==(other) &&
    index == other.index &&
    domainsCount == other.domainsCount;
}

template <typename Space>
void DistributedMeshIOCommon<Space>::Load(const std::string& fileName, IO::FileType fileType)
{
  std::fstream file;

  switch (fileType)
  {
    case IO::Ascii:  file.open(fileName.c_str(), std::fstream::in);                        break;
    case IO::Binary: file.open(fileName.c_str(), std::fstream::in | std::fstream::binary); break;
  }
  
  if (file.fail()) return;

  MeshIO<Space>::Load(file, fileType);

  switch (fileType)
  {
    case IO::Ascii:  file >> domainsCount;         break;
    case IO::Binary: IO::Read(file, domainsCount); break;
  }
  LoadTransitionInfo(file, fileType, &transitionInfos);
  file.close();
}

template <typename Space>
void DistributedMeshIOCommon<Space>::Save(const std::string& fileName, IO::FileType fileType)
{
  std::fstream file;

  switch (fileType)
  {
    case IO::Ascii:  file.open(fileName.c_str(), std::fstream::out);                        break;
    case IO::Binary: file.open(fileName.c_str(), std::fstream::out | std::fstream::binary); break;
  }
  
  if (file.fail()) return;

  MeshIO<Space>::Save(file, fileType);

  switch (fileType)
  {
    case IO::Ascii:  file << domainsCount << std::endl;     break;
    case IO::Binary: IO::Write(file, domainsCount);         break;
  }
  SaveTransitionInfo(file, fileType, transitionInfos);
  file.close();
}

template <typename Space>
void DistributedMeshIOCommon<Space>::SaveTransitionInfo(std::fstream& file, IO::FileType fileType, 
  const std::vector< TransitionInfo<Space> >& infos)
{  
  for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
  {
    IO:: template Save<IndexType>(     file, fileType, infos[domainIndex].nodesIndices);
    IO:: template Save<Cell>(          file, fileType, infos[domainIndex].cells);
    IO:: template Save<TransitionNode>(file, fileType, infos[domainIndex].transitionNodes);
  }
}

template <typename Space>
void DistributedMeshIOCommon<Space>::LoadTransitionInfo(std::fstream& file, IO::FileType fileType, 
  std::vector< TransitionInfo<Space> >* const infos)
{
  for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
  {
    IO:: template Load<IndexType>(     file, fileType, &((*infos)[domainIndex].nodesIndices));
    IO:: template Load<Cell>(          file, fileType, &((*infos)[domainIndex].cells));
    IO:: template Load<TransitionNode>(file, fileType, &((*infos)[domainIndex].transitionNodes));
  }
}
