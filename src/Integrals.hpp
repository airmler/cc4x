#include "Util.hpp"

#ifndef INTEGRALS_DEFINED
#define INTEGRALS_DEFINED


namespace Integrals {
  struct input {
    CTF::bsTensor<Complex> *pp;
    CTF::bsTensor<Complex> *ph;
    CTF::bsTensor<Complex> *hp;
    CTF::bsTensor<Complex> *hh;
  };
  struct output {
    CTF::bsTensor<Complex> **Vpphh;
    CTF::bsTensor<Complex> **Vphhp;
    CTF::bsTensor<Complex> **Vhhpp;
  };
  void run(input const& in, output &out);

}


#endif
