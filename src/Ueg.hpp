#include "Util.hpp"

#ifndef UEG_DEFINED
#define UEG_DEFINED

using iarr  = array<int64_t,3>;
using darr  = array<double,4>;


namespace Ueg {
  struct input {
    int64_t No, Nv, NF;
    double rs;
  };
  struct output {
    tensor<Complex> **coulombVertex;
    tensor<Complex> **eps;
  };
  void run(input const& in, output &out);

}


#endif
