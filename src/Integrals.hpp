#include "Util.hpp"

#ifndef INTEGRALS_DEFINED
#define INTEGRALS_DEFINED


namespace Integrals {
  struct input {
    CTF::bsTensor<Complex> *hh;
    CTF::bsTensor<Complex> *ph;
    CTF::bsTensor<Complex> *hp;
    CTF::bsTensor<Complex> *pp;
  };
  struct output {
    CTF::bsTensor<Complex> **Vhhhh;
    CTF::bsTensor<Complex> **Vhhhp;
    CTF::bsTensor<Complex> **Vhhph;
    CTF::bsTensor<Complex> **Vhhpp;
    CTF::bsTensor<Complex> **Vhphp;
    CTF::bsTensor<Complex> **Vhppp;
    CTF::bsTensor<Complex> **Vphhh;
    CTF::bsTensor<Complex> **Vphhp;
    CTF::bsTensor<Complex> **Vphph;
    CTF::bsTensor<Complex> **Vphpp;
    CTF::bsTensor<Complex> **Vpphh;
    CTF::bsTensor<Complex> **Vpphp;
    CTF::bsTensor<Complex> **Vppph;
    CTF::bsTensor<Complex> **Vpppp;
  };
  void run(input const& in, output &out);

}

#endif
