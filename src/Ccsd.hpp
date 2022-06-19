#include "Util.hpp"

#ifndef CCSD_DEFINED
#define CCSD_DEFINED


namespace Ccsd {
  struct input {
    CTF::bsTensor<Complex> *Vhhhh;
    CTF::bsTensor<Complex> *Vhhhp;
    CTF::bsTensor<Complex> *Vhhph;
    CTF::bsTensor<Complex> *Vhhpp;
    CTF::bsTensor<Complex> *Vhphp;
    CTF::bsTensor<Complex> *Vhppp;
    CTF::bsTensor<Complex> *Vphhh;
    CTF::bsTensor<Complex> *Vphhp;
    CTF::bsTensor<Complex> *Vphph;
    CTF::bsTensor<Complex> *Vphpp;
    CTF::bsTensor<Complex> *Vpphh;
    CTF::bsTensor<Complex> *Vpphp;
    CTF::bsTensor<Complex> *Vppph;
    CTF::bsTensor<Complex> *Vpppp;
    CTF::bsTensor<Complex> *epsi;
    CTF::bsTensor<Complex> *epsa;
  };
  struct output {
  };
  void run(input const& in, output &out);

}


#endif
