#include "Util.hpp"

#ifndef SLICE_DEFINED
#define SLICE_DEFINED


namespace Slice {
  struct input {
    CTF::bsTensor<Complex> *I;
    std::vector<int64_t> partitionPoint;
  };
  struct output {
    CTF::bsTensor<Complex> **A;
    CTF::bsTensor<Complex> **B;
    CTF::bsTensor<Complex> **C;
    CTF::bsTensor<Complex> **D;
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
