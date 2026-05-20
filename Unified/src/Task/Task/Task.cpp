// #define _HAS_EXCEPTIONS 0
// #define _ITERATOR_DEBUG_LEVEL 0 
// #include <xmmintrin.h>
// #include <vld.h>

#include "Task.h"
#include <cstdlib>
#include <string>

int main(int argc, char* argv[])
{
  int dimCount, polyOrder;
  std::string taskName;
  
  if(argc == 1) 
  {
    BasicSettings settings; 
    settings.Parse("task.xml");

    dimCount = settings.configDimsCount;
    polyOrder = settings.configPolyOrder;
    taskName = settings.taskName;
  } else 
  {
    dimCount = std::strtol(argv[1], nullptr, 10);
    polyOrder = std::strtol(argv[2], nullptr, 10);
    taskName = argv[3];
  }

  switch(dimCount)
  {
    case 2:
    {
      switch (polyOrder)
      {
        case 1: {Task<Space2, 1> task(&taskName); task.Run();} break; 
        case 2: {Task<Space2, 2> task(&taskName); task.Run();} break; 
        case 3: {Task<Space2, 3> task(&taskName); task.Run();} break; 
        case 4: {Task<Space2, 4> task(&taskName); task.Run();} break; 
        case 5: {Task<Space2, 5> task(&taskName); task.Run();} break; 
        default: std::cerr << "Unknown polynomial order"; break;
      }
    }break;
    case 3:
    {
      switch (polyOrder)
      {
        case 1: {Task<Space3, 1> task(&taskName); task.Run();} break; 
        case 2: {Task<Space3, 2> task(&taskName); task.Run();} break; 
        case 3: {Task<Space3, 3> task(&taskName); task.Run();} break; 
        case 4: {Task<Space3, 4> task(&taskName); task.Run();} break; 
        case 5: {Task<Space3, 5> task(&taskName); task.Run();} break; 
        default: std::cerr << "Unknown polynomial order"; break;
      }
    }break;
    default: std::cerr << "Unknown dims count"; break;
  }

  return 0;
}
