#include "Util.hpp"

#ifndef SCF_DEFINED
#define SCF_DEFINED


namespace Scf {

  void evalEnergy(
    CTF::bsTensor<Complex> &T, CTF::bsTensor<Complex> &V, std::string what = ""
  );

  void getDabij( CTF::bsTensor<Complex> &Dabij
               , CTF::bsTensor<Complex> &epsi
               , CTF::bsTensor<Complex> &epsa);

  void getDai( CTF::bsTensor<Complex> &Dai
             , CTF::bsTensor<Complex> &epsi
             , CTF::bsTensor<Complex> &epsa);



}


#endif
