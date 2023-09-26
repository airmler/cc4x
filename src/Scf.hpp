#include "Util.hpp"

#ifndef SCF_DEFINED
#  define SCF_DEFINED

namespace Scf {

void evalEnergy(tensor<Complex> &T, tensor<Complex> &V, std::string what = "");

void getDabij(tensor<Complex> &Dabij,
              tensor<Complex> &epsi,
              tensor<Complex> &epsa);

void getDai(tensor<Complex> &Dai, tensor<Complex> &epsi, tensor<Complex> &epsa);

} // namespace Scf

#endif
