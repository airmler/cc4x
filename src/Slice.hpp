#include "Util.hpp"

#ifndef SLICE_DEFINED
#define SLICE_DEFINED


namespace Slice {
  struct input {
    tensor<Complex> *I;
    std::vector<size_t> partitionPoint;
  };
  struct output {
    tensor<Complex> **A;
    tensor<Complex> **B;
    tensor<Complex> **C;
    tensor<Complex> **D;
  };
  void run(input const& in, output &out);

  struct sliceDim {
    std::vector< std::vector<size_t> > srcBegin;
    std::vector< std::vector<size_t> > srcEnd;
    std::vector< std::vector<size_t> > dstBegin;
    std::vector< std::vector<size_t> > dstEnd;
  };
}


#endif
