#include "Util.hpp"

#ifndef INTEGRALS_DEFINED
#define INTEGRALS_DEFINED


namespace Integrals {
  struct input {
    tensor<Complex> *hh;
    tensor<Complex> *ph;
    tensor<Complex> *hp;
    tensor<Complex> *pp;
  };
  struct output {
    tensor<Complex> **Vhhhh;
    tensor<Complex> **Vhhhp;
    tensor<Complex> **Vhhph;
    tensor<Complex> **Vhhpp;
    tensor<Complex> **Vhphp;
    tensor<Complex> **Vhppp;
    tensor<Complex> **Vphhh;
    tensor<Complex> **Vphhp;
    tensor<Complex> **Vphph;
    tensor<Complex> **Vphpp;
    tensor<Complex> **Vpphh;
    tensor<Complex> **Vpphp;
    tensor<Complex> **Vppph;
    tensor<Complex> **Vpppp;
  };
  void run(input const& in, output &out);

}

#endif
