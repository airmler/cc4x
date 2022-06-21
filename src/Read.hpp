#include "Util.hpp"

#ifndef READ_DEFINED
#define READ_DEFINED


namespace Read {
  struct input {
    std::string fileName;
  };
  struct output {
    tensor<Complex> **T;
  };
  void run(input const& in, output &out);
  void getAmplitudesType(std::string fileName);

  struct yamlData {
    std::string fileName;
    std::string scalarType;
    std::string fileType;
    int order = 0;
    std::vector<int64_t> lens;
    // The following is metaData which may or may not be given in the tensor
    int No = 0;
    int Nv = 0;
    int halfGrid = 0;
    ivec kMesh = {1,1,1}; // if no mesh is provided we work with {1,1,1} mesh
  }; 



}


#endif
