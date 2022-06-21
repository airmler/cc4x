#include "Util.hpp"

#ifndef SLICE_DEFINED
#define SLICE_DEFINED


namespace Slice {
  struct input {
    tensor<Complex> *I;
    std::vector<int64_t> partitionPoint;
  };
  struct output {
    tensor<Complex> **A;
    tensor<Complex> **B;
    tensor<Complex> **C;
    tensor<Complex> **D;
  };
  void run(input const& in, output &out);

  struct sliceDim {
    std::vector<i64vec> srcBegin;
    std::vector<i64vec> srcEnd;
    std::vector<i64vec> dstBegin;
    std::vector<i64vec> dstEnd;
  };
}


#endif
