#include "Drccd.hpp"
#include "Util.hpp"
#include <cc4x.hpp>

namespace Drccd{

  Complex conjo(Complex x){ return std::conj(x); }
  Complex   inv(Complex x){ return 1.0/x; }

  void run(input const& in, output& out){
    // sanity checks
    if ( in.Vpphh == NULL || in.Vpphh == (CTF::bsTensor<Complex>*)0xfafa
      || in.Vphhp == NULL || in.Vphhp == (CTF::bsTensor<Complex>*)0xfafa
      || in.Vhhpp == NULL || in.Vhhpp == (CTF::bsTensor<Complex>*)0xfafa
      || in.epsi == NULL  || in.epsi  == (CTF::bsTensor<Complex>*)0xfafa
      || in.epsa == NULL  || in.epsa  == (CTF::bsTensor<Complex>*)0xfafa
       ) {
      THROW("Input of Drccd not valid");
    }

    int64_t h(cc4x::No), p(cc4x::Nv);
    auto nzc4(cc4x::kmesh->getNZC(4));
    CTF::bsTensor<Complex> Dabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Dabij");
    CTF::bsTensor<Complex> Tabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Tabij");
    CTF::bsTensor<Complex> energy(0, {}, cc4x::kmesh->getNZC(0), cc4x::dw, "energy");

    CTF::Univar_Function<Complex> fConj(&conjo);
    CTF::Univar_Function<Complex> fInv(&inv);


    Dabij.sum(1.0, *in.epsi, "i", 0.0, "abij");
    Dabij.sum(1.0, *in.epsi, "j", 1.0, "abij");
    Dabij.sum(-1., *in.epsa, "a", 1.0, "abij");
    Dabij.sum(-1., *in.epsa, "b", 1.0, "abij");
    Dabij.sum(1.0, Dabij, "abij", 0.0, "abij", fInv);

    Tabij.contract(1.0, Dabij, "abij", *in.Vpphh, "abij", 0.0, "abij");
    //Tabij.contract(1.0, Dabij, "abij", *in.Vhhpp, "ijab", 0.0, "abij");
    //Tabij.sum(1.0, Tabij, "abij", 0.0, "abij", fConj);

    std::complex<double> direct, exchange;
    energy.contract(2.0, Tabij, "abij", *in.Vhhpp, "ijab", 0.0, "");
    energy.read_all(&direct);
    energy.contract(-1.0, Tabij, "abij", *in.Vhhpp, "jiab", 0.0, "");
    energy.read_all(&exchange);
    int Nk(cc4x::kmesh->Nk);
    int N3(Nk*Nk*Nk);
    std::cout << "direct: " << real(direct);
    if (std::abs(imag(direct)) > 1e-5)
      std::cout << " + imag! " << imag(direct);
    std::cout << std::endl;
    std::cout << "exchange: " << real(exchange);
    if (std::abs(imag(exchange)) > 1e-5)
      std::cout << " + imag! " << imag(exchange);
    std:cout << std::endl;

    return;
  }

}
