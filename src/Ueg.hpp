#include "Util.hpp"

#ifndef UEG_DEFINED
#define UEG_DEFINED


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
