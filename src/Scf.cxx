#include "Util.hpp"
#include "Scf.hpp"
#include "cc4x.hpp"

namespace Scf{

  Complex   inv(Complex x){ return 1.0/x; }

  void evalEnergy(
    tensor<Complex> &T, tensor<Complex> &V, std::string what
  ){

    tensor<Complex> energy(0, {}, cc4x::kmesh->getNZC(0), cc4x::dw, "e");
    std::complex<double> direct, exchange;
    energy.contract(2.0, T, "abij", V, "ijab", 0.0, "");
    //std::memcpy(&direct, energy.tensors[0]->data, sizeof(Complex));
    energy.read_all(&direct);
    energy.contract(-1.0, T, "abij", V, "jiab", 0.0, "");
    energy.read_all(&exchange);
    int Nk(cc4x::kmesh->Nk);
    int N3(Nk*Nk*Nk);
    if (what.size()) LOG() << what << "\n";
    LOG() << "Energy (d/x): " << real(direct+exchange)/N3;
    LOG() << " ( " << real(direct)/N3;
    LOG() << " / " << real(exchange)/N3;
    LOG() <<  " )" << std::endl;
    if (std::abs(imag(direct)) > 1e-5){
      LOG() << "\nWarning. direct part is imaginary " << imag(direct)/N3 << '\n';
    }
    if (std::abs(imag(exchange)) > 1e-5){
      LOG() << "\nWarning. exchange part is imaginary " << imag(exchange)/N3 << '\n';
    }
  }

  void getDabij( tensor<Complex> &Dabij
               , tensor<Complex> &epsi
               , tensor<Complex> &epsa){

    std::function<Complex(const Complex)> fInv(&inv);
    Dabij.sum(1.0, epsi, "i", 0.0, "abij");
    Dabij.sum(1.0, epsi, "j", 1.0, "abij");
    Dabij.sum(-1., epsa, "a", 1.0, "abij");
    Dabij.sum(-1., epsa, "b", 1.0, "abij");
    Dabij.sum(1.0, Dabij, "abij", 0.0, "abij", fInv);
  }

  void getDai( tensor<Complex> &Dai
             , tensor<Complex> &epsi
             , tensor<Complex> &epsa){

    std::function<Complex(const Complex)> fInv(&inv);
    Dai.sum(1.0, epsi, "i", 0.0, "ai");
    Dai.sum(-1., epsa, "a", 1.0, "ai");
    Dai.sum(1.0, Dai, "ai", 0.0, "ai", fInv);
  }

}
