#include "Util.hpp"

#ifndef CCSD_DEFINED
#define CCSD_DEFINED


namespace Ccsd {
  struct input {
    tensor<Complex> *Vhhhh;
    tensor<Complex> *Vhhhp;
    tensor<Complex> *Vhhph;
    tensor<Complex> *Vhhpp;
    tensor<Complex> *Vphhh;
    tensor<Complex> *Vphhp;
    tensor<Complex> *Vphph;
    tensor<Complex> *Vpphh;
    tensor<Complex> *epsi;
    tensor<Complex> *epsa;
    tensor<Complex> *hhVertex;
    tensor<Complex> *phVertex;
    tensor<Complex> *hpVertex;
    tensor<Complex> *ppVertex;
  };
  struct output {
  };
  void run(input const& in, output &out);

}


#endif
