#include "Drccd.hpp"
#include "Util.hpp"
#include <cc4x.hpp>

namespace Drccd{

  Complex conjo(Complex x){ return std::conj(x); }
  Complex   inv(Complex x){ return 1.0/x; }

  void evalEnergy(
    CTF::bsTensor<Complex> &T, CTF::bsTensor<Complex> &V, std::string what
  ){

    CTF::bsTensor<Complex> energy(0, {}, cc4x::kmesh->getNZC(0), cc4x::dw, "e");
    std::complex<double> direct, exchange;
    energy.contract(2.0, T, "abij", V, "ijab", 0.0, "");
    //std::memcpy(&direct, energy.tensors[0]->data, sizeof(Complex));
    energy.read_all(&direct);
    energy.contract(-1.0, T, "abij", V, "jiab", 0.0, "");
    energy.read_all(&exchange);
    int Nk(cc4x::kmesh->Nk);
    int N3(Nk*Nk*Nk);
    LOG() << what << "\n";
    LOG() << "Total: " << real(direct+exchange)/N3;
    LOG() << " ; direct: " << real(direct)/N3;
    if (std::abs(imag(direct)) > 1e-5){
      LOG() << " + imag! " << imag(direct)/N3;
    }
    LOG() << "exchange: " << real(exchange)/N3;
    if (std::abs(imag(exchange)) > 1e-5){
      LOG() << " + imag! " << imag(exchange)/N3;
    }
    std::cout << std::endl;
  }


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
    CTF::bsTensor<Complex> Rabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Rabij");
    CTF::bsTensor<Complex> Tabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Tabij");
    CTF::bsTensor<Complex> Calid(4, {p,h,h,p}, nzc4, cc4x::dw, "Calid");
    CTF::bsTensor<Complex> Xijab(4, {h,h,p,p}, nzc4, cc4x::dw, "Xijab");
    CTF::bsTensor<Complex> energy(0, {}, cc4x::kmesh->getNZC(0), cc4x::dw, "energy");

    CTF::Univar_Function<Complex> fConj(&conjo);
    CTF::Univar_Function<Complex> fInv(&inv);


    Dabij.sum(1.0, *in.epsi, "i", 0.0, "abij");
    Dabij.sum(1.0, *in.epsi, "j", 1.0, "abij");
    Dabij.sum(-1., *in.epsa, "a", 1.0, "abij");
    Dabij.sum(-1., *in.epsa, "b", 1.0, "abij");
    Dabij.sum(1.0, Dabij, "abij", 0.0, "abij", fInv);

    Tabij.contract(1.0, Dabij, "abij", *in.Vpphh, "abij", 0.0, "abij");

    evalEnergy(Tabij, *in.Vhhpp, "MP2");

    Rabij.sum(1.0, *in.Vpphh, "abij", 0.0, "abij");
    Rabij.contract(2.0, *in.Vphhp, "akic", Tabij, "cbkj", 0.0, "abij");
    Rabij.contract(2.0, *in.Vphhp, "bkjc", Tabij, "acik", 1.0, "abij");


//    Xijab.sum(2.0, *in.Vhhpp, "ijab", 0.0, "ijab");
//    Calid.contract(2.0, Xijab, "klcd", Tabij, "acik", 1.0, "alid");
//    Rabij.contract(1.0, Calid, "alid", Tabij, "dblj", 1.0, "abij");


    Tabij.contract(1.0, Dabij, "abij", Rabij, "abij", 0.0, "abij");
    //Tabij.contract(1.0, Dabij, "abij", *in.Vhhpp, "ijab", 0.0, "abij");

    evalEnergy(Tabij, *in.Vhhpp, "Drccd");


    return;
  }

}
