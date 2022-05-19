#include "Util.hpp"

#ifndef DRCCD_DEFINED
#define DRCCD_DEFINED


namespace Drccd {
  struct input {
    CTF::bsTensor<Complex> *Vpphh;
    CTF::bsTensor<Complex> *Vphhp;
    CTF::bsTensor<Complex> *Vhhpp;
    CTF::bsTensor<Complex> *epsi;
    CTF::bsTensor<Complex> *epsa;
  };
  struct output {
  };
  void run(input const& in, output &out);

}


#endif
