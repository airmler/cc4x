#include "Util.hpp"

#ifndef DRCCD_DEFINED
#define DRCCD_DEFINED


namespace Drccd {
  struct input {
    tensor<Complex> *Vpphh;
    tensor<Complex> *Vphhp;
    tensor<Complex> *Vhhpp;
    tensor<Complex> *epsi;
    tensor<Complex> *epsa;
  };
  struct output {
  };
  void run(input const& in, output &out);

}


#endif
