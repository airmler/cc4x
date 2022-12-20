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
    size_t order = 0;
    std::vector<size_t> lens;
    // The following is metaData which may or may not be given in the tensor
    int No = 0;
    int Nv = 0;
    int halfGrid = 0;
    std::vector<size_t> kMesh = {1,1,1};
  }; 



}


#endif
