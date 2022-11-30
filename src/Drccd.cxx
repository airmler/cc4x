#include "Drccd.hpp"
#include "Util.hpp"
#include "Scf.hpp"
#include "cc4x.hpp"

namespace Drccd{

  void run(input const& in, output& out){
    Timings chrono;
    chrono["drccd"].start();

    bool linearized(false);
    // sanity checks
    if ( in.Vpphh == NULL || in.Vpphh == (tensor<Complex>*)0xfafa
      || in.Vphhp == NULL || in.Vphhp == (tensor<Complex>*)0xfafa
      || in.Vhhpp == NULL || in.Vhhpp == (tensor<Complex>*)0xfafa
      || in.epsi == NULL  || in.epsi  == (tensor<Complex>*)0xfafa
      || in.epsa == NULL  || in.epsa  == (tensor<Complex>*)0xfafa
       ) {
      THROW("Input of Drccd not valid");
    }

    int64_t h(cc4x::No), p(cc4x::Nv);
    auto nzc4(cc4x::kmesh->getNZC(4));
    tensor<Complex> Dabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Dabij");
    tensor<Complex> Rabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Rabij");
    tensor<Complex> Tabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Tabij");
    tensor<Complex> Calid(4, {p,h,h,p}, nzc4, cc4x::dw, "Calid");
    tensor<Complex> Xijab(4, {h,h,p,p}, nzc4, cc4x::dw, "Xijab");

    tensor<Complex> Kabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Kabij");
    Scf::getDabij(Dabij, *in.epsi, *in.epsa);

    Tabij.contract(1.0, Dabij, "abij", *in.Vpphh, "abij", 0.0, "abij");

    Scf::evalEnergy(Tabij, *in.Vhhpp, "MP2");

    for (size_t i(0); i < cc4x::iterations; i++){
      Rabij.contract(2.0, *in.Vphhp, "akic", Tabij, "cbkj", 0.0, "abij");
//      Rabij.contract(2.0, *in.Vphhp, "bkjc", Tabij, "acik", 1.0, "abij");
      Kabij.sum(1.0, Rabij, "baji", 0.0, "abij");
      Rabij.sum(1.0, Kabij, "abij", 1.0, "abij");
      Rabij.sum(1.0, *in.Vpphh, "abij", 1.0, "abij");
      if (!linearized) {
        Xijab.sum(2.0, *in.Vhhpp, "ijab", 0.0, "ijab");
        Calid.contract(2.0, Xijab, "klcd", Tabij, "acik", 0.0, "alid");
        Rabij.contract(1.0, Calid, "alid", Tabij, "dblj", 1.0, "abij");
      }
      Tabij.contract(1.0, Dabij, "abij", Rabij, "abij", 0.0, "abij");

      Scf::evalEnergy(Tabij, *in.Vhhpp);
    }

    Scf::evalEnergy(Tabij, *in.Vhhpp, "Drccd");

    chrono["drccd"].stop();
    for (auto t: chrono)
      LOG() << t.first << " " << t.second.count() << "\n";

    return;
  }

}
