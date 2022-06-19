#include "Drccd.hpp"
#include "Util.hpp"
#include "Scf.hpp"
#include "cc4x.hpp"

namespace Drccd{

  void run(input const& in, output& out){
    bool linearized(false);
    size_t maxIter(20);
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

    Scf::getDabij(Dabij, *in.epsi, *in.epsa);

    Tabij.contract(1.0, Dabij, "abij", *in.Vpphh, "abij", 0.0, "abij");

    Scf::evalEnergy(Tabij, *in.Vhhpp, "MP2");

    for (size_t i(0); i < maxIter; i++){
      Rabij.sum(1.0, *in.Vpphh, "abij", 0.0, "abij");
      Rabij.contract(2.0, *in.Vphhp, "akic", Tabij, "cbkj", 1.0, "abij");
      Rabij.contract(2.0, *in.Vphhp, "bkjc", Tabij, "acik", 1.0, "abij");
      if (!linearized) {
        Xijab.sum(2.0, *in.Vhhpp, "ijab", 0.0, "ijab");
        Calid.contract(2.0, Xijab, "klcd", Tabij, "acik", 0.0, "alid");
        Rabij.contract(1.0, Calid, "alid", Tabij, "dblj", 1.0, "abij");
      }
      Tabij.contract(1.0, Dabij, "abij", Rabij, "abij", 0.0, "abij");

      Scf::evalEnergy(Tabij, *in.Vhhpp);
    }

    Scf::evalEnergy(Tabij, *in.Vhhpp, "Drccd");

    return;
  }

}
