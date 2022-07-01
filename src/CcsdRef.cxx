#include "CcsdRef.hpp"
#include "Util.hpp"
#include "Scf.hpp"
#include "cc4x.hpp"

namespace CcsdRef{

  void run(input const& in, output& out){
    Timings chrono;
    chrono["ccsd"].start();
    // sanity checks
    if ( in.Vhhhh == NULL    || in.Vhhhh == (tensor<Complex>*)0xfafa
      || in.Vhhhp == NULL    || in.Vhhhp == (tensor<Complex>*)0xfafa
      || in.Vhhph == NULL    || in.Vhhph == (tensor<Complex>*)0xfafa
      || in.Vhhpp == NULL    || in.Vhhpp == (tensor<Complex>*)0xfafa
      || in.Vhppp == NULL    || in.Vhppp == (tensor<Complex>*)0xfafa
      || in.Vphhp == NULL    || in.Vphhp == (tensor<Complex>*)0xfafa
      || in.Vphph == NULL    || in.Vphph == (tensor<Complex>*)0xfafa
      || in.Vphpp == NULL    || in.Vphpp == (tensor<Complex>*)0xfafa
      || in.Vpphh == NULL    || in.Vpphh == (tensor<Complex>*)0xfafa
      || in.Vppph == NULL    || in.Vppph == (tensor<Complex>*)0xfafa
      || in.Vpppp == NULL    || in.Vpppp == (tensor<Complex>*)0xfafa
      || in.epsi ==  NULL    || in.epsi  == (tensor<Complex>*)0xfafa
      || in.epsa ==  NULL    || in.epsa  == (tensor<Complex>*)0xfafa
       ) {
      THROW("Input of CcsdRef not valid");
    }

		int64_t h(cc4x::No), p(cc4x::Nv);
    auto nzc4(cc4x::kmesh->getNZC(4));
    auto nzc3(cc4x::kmesh->getNZC(3));
    auto nzc2(cc4x::kmesh->getNZC(2));
    auto nzc1(cc4x::kmesh->getNZC(1));


    tensor<Complex> Dai(2, {p,h}, nzc2, cc4x::dw, "Dai");
    tensor<Complex> Dabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Dabij");
    tensor<Complex> Tai(2, {p,h}, nzc2, cc4x::dw, "Tai");
    tensor<Complex> Tabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Tabij");

    tensor<Complex> Kac(2, {p,p}, nzc2, cc4x::dw, "Kac");
    tensor<Complex> Kki(2, {h,h}, nzc2, cc4x::dw, "Kki");
    tensor<Complex> Lac(2, {p,p}, nzc2, cc4x::dw, "Lac");
    tensor<Complex> Lki(2, {h,h}, nzc2, cc4x::dw, "Lki");
    tensor<Complex> Kck(2, {p,h}, nzc2, cc4x::dw, "Kck");
    tensor<Complex> Xakic(4, {p,h,h,p}, nzc4, cc4x::dw, "Xakic");
    tensor<Complex> Xakci(4, {p,h,p,h}, nzc4, cc4x::dw, "Xakci");
    tensor<Complex> Xklij(4, {h,h,h,h}, nzc4, cc4x::dw, "Xklij");
    tensor<Complex> Xabcd(4, {p,p,p,p}, nzc4, cc4x::dw, "Xabij");
    tensor<Complex> Xakij(4, {p,h,h,h}, nzc4, cc4x::dw, "Xakij");
    tensor<Complex> Xabic(4, {p,p,h,p}, nzc4, cc4x::dw, "Xabic");

    tensor<Complex> Xabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Xabij");
    tensor<Complex> Yabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Yabij");
    tensor<Complex> Rai(2, {p,h}, nzc2, cc4x::dw, "Rai");
    tensor<Complex> Rabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Rabij");
    tensor<Complex> Kabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Kabij");


    Scf::getDai(Dai, *in.epsi, *in.epsa);
    Scf::getDabij(Dabij, *in.epsi, *in.epsa);

    Tabij.contract(1.0, Dabij, "abij", *in.Vpphh, "abij", 0.0, "abij");

    Scf::evalEnergy(Tabij, *in.Vhhpp, "MP2");

    for (size_t i(0); i < cc4x::iterations; i++){
      Xabij.sum(1.0, Tabij, "abij", 0.0, "abij");
      Xabij.contract(1.0, Tai, "ai", Tai, "bj", 1.0, "abij");
      Yabij.sum(1.0, Tabij, "abij", 0.0, "abij");
      Yabij.contract(2.0, Tai, "ai", Tai, "bj", 1.0, "abij");

      // add kappa and lambda contributions
      chrono["ccsd - laka"].start();
      Kac.contract(-2.0, *in.Vhhpp, "klcd", Xabij, "adkl", 0.0, "ac");
      Kac.contract( 1.0, *in.Vhhpp, "kldc", Xabij, "adkl", 1.0, "ac");
      Lac.sum(1.0, Kac, "ac", 0.0, "ac");

      Lac.contract( 2.0, *in.Vphpp, "akcd", Tai, "dk", 1.0, "ac");
      Lac.contract(-1.0, *in.Vphpp, "akdc", Tai, "dk", 1.0, "ac");
      Kki.contract( 2.0, *in.Vhhpp, "klcd", Xabij, "cdil", 0.0, "ki");
      Kki.contract(-1.0, *in.Vhhpp, "kldc", Xabij, "cdil", 1.0, "ki");
      Lki.sum(1.0, Kki, "ki", 0.0, "ki");
      Lki.contract( 2.0, *in.Vhhhp, "klic", Tai,"cl", 1.0, "ki");
      Lki.contract(-1.0, *in.Vhhph, "klci", Tai,"cl", 1.0, "ki");

      Rabij.contract( 1.0, Lac, "ac", Tabij, "cbij", 0.0, "abij");
      Rabij.contract(-1.0, Lki, "ki", Tabij, "abkj", 1.0, "abij");
      chrono["ccsd - laka"].stop();

      // add T1 -> R2 contributions
      chrono["ccsd - t1r2"].start();
      Xakij.sum(1.0, *in.Vphhh, "akij", 0.0, "akij");
      Xakij.contract(1.0, *in.Vphhp, "akic", Tai, "cj", 1.0, "akij");
      Rabij.contract(-1.0, Xakij, "akij", Tai, "bk", 1.0, "abij");
      Xabic.sum(1.0, *in.Vpphp, "abic", 0.0, "abic");
      Xabic.contract(-1.0, *in.Vhphp, "kbic", Tai, "ak", 1.0, "abic");
      Rabij.contract(1.0, Xabic, "abic", Tai, "cj", 1.0, "abij");
      chrono["ccsd - t1r2"].stop();
      // Xakic
      chrono["ccsd - akic"].start();
      Xakic.sum(1.0, *in.Vphhp, "akic", 0.0, "akic");
      Xakic.contract(-1.0, *in.Vhhhp, "lkic", Tai, "al", 1.0, "akic");
      Xakic.contract( 1.0, *in.Vphpp, "akdc", Tai, "di", 1.0, "akic");
      Xakic.contract(-0.5, *in.Vhhpp, "lkdc", Yabij, "dail", 1.0, "akic");
      Xakic.contract( 1.0, *in.Vhhpp, "lkdc", Tabij, "adil", 1.0, "akic");
      Xakic.contract(-0.5, *in.Vhhpp, "lkcd", Tabij, "adil", 1.0, "akic");
      Rabij.contract( 2.0, Xakic, "akic", Tabij, "cbkj", 1.0, "abij");
      Rabij.contract(-1.0, Xakic, "akic", Tabij, "bckj", 1.0, "abij");
      chrono["ccsd - akic"].stop();

      // Xakci
      chrono["ccsd - akci"].start();
      Xakci.sum(1.0, *in.Vphph, "akci", 0.0, "akci");
      Xakci.contract(-1.0, *in.Vhhph, "lkci", Tai, "al", 1.0, "akci");
      Xakci.contract( 1.0, *in.Vphpp, "akcd", Tai, "di", 1.0, "akci");
      Xakci.contract(-0.5, *in.Vhhpp, "lkcd", Yabij, "dail", 1.0, "akci");
      Rabij.contract(-1.0, Xakci, "akci", Tabij, "cbkj", 1.0, "abij");
      Rabij.contract(-1.0, Xakci, "bkci", Tabij, "ackj", 1.0, "abij");
      chrono["ccsd - akci"].stop();

      // Permutation and add V
      Kabij.sum(1.0, Rabij, "abij", 0.0, "baji");
      Rabij.sum(1.0, Kabij, "abij", 1.0, "abij");
      Rabij.sum(1.0, *in.Vpphh, "abij", 1.0, "abij");

      // Xabcd
      chrono["ccsd - pp-ladder"].start();
      Xabcd.sum(1.0, *in.Vpppp, "abcd", 0.0, "abcd");
      Xabcd.contract(-1.0, *in.Vphpp, "akcd", Tai, "bk", 1.0, "abcd");
      Xabcd.contract(-1.0, *in.Vhppp, "kbcd", Tai, "ak", 1.0, "abcd");
      Rabij.contract(1.0, Xabcd, "abcd", Xabij, "cdij", 1.0, "abij");
      chrono["ccsd - pp-ladder"].stop();


      // Xijkl
      chrono["ccsd - hh-ladder"].start();
      Xklij.sum(1.0, *in.Vhhhh, "klij", 0.0, "klij");
      Xklij.contract(1.0, *in.Vhhhp, "klic", Tai, "cj", 1.0, "klij");
      Xklij.contract(1.0, *in.Vhhph, "klcj", Tai, "ci", 1.0, "klij");
      Xklij.contract(1.0, *in.Vhhpp, "klcd", Xabij, "cdij", 1.0, "klij");
      Rabij.contract(1.0, Xklij, "klij", Xabij, "abkl", 1.0, "abij");
      chrono["ccsd - hh-ladder"].stop();

      // Singles contribution
      chrono["ccsd - singles"].start();
      Rai.contract( 1.0, Kac, "ac", Tai, "ci", 0.0, "ai");
      Rai.contract(-1.0, Kki, "ki", Tai, "ak", 1.0, "ai");

      Kck.contract( 2.0, *in.Vhhpp, "klcd", Tai, "dl", 0.0, "ck");
      Kck.contract(-1.0, *in.Vhhpp, "kldc", Tai, "dl", 1.0, "ck");

      Rai.contract( 2.0, Kck, "ck", Tabij, "caki", 1.0, "ai");
      Rai.contract(-1.0, Kck, "ck", Tabij, "caik", 1.0, "ai");
      Yabij.contract(1.0, Tai, "ci", Tai, "ak", 0.0, "caik");
      Rai.contract(1.0, Kck, "ck", Yabij, "caik", 1.0, "ai");

      Rai.contract( 2.0, *in.Vphhp, "akic", Tai, "ck", 1.0, "ai");
      Rai.contract(-1.0, *in.Vphph, "akci", Tai, "ck", 1.0, "ai");

      Rai.contract( 2.0, *in.Vphpp, "akcd", Xabij, "cdik", 1.0, "ai");
      Rai.contract(-1.0, *in.Vphpp, "akdc", Xabij, "cdik", 1.0, "ai");

      Rai.contract(-2.0, *in.Vhhhp, "klic", Xabij, "ackl", 1.0, "ai");
      Rai.contract( 1.0, *in.Vhhhp, "lkic", Xabij, "ackl", 1.0, "ai");
      chrono["ccsd - singles"].stop();

      // add energy denominator
      Tabij.contract(1.0, Dabij, "abij", Rabij, "abij", 0.0, "abij");
      Tai.contract(1.0, Dai, "ai", Rai, "ai", 0.0, "ai");
      // As long as we dont want to use a fancy mixer we can write like that
      Xabij.sum(1.0, Tabij, "abij", 0.0, "abij");
      Xabij.contract(1.0, Tai, "ai", Tai, "bj", 1.0, "abij");
      Scf::evalEnergy(Xabij, *in.Vhhpp);
    }

    Scf::evalEnergy(Xabij, *in.Vhhpp, "Ccsd");

    chrono["ccsd"].stop();
    for (auto t: chrono)
      LOG() << t.first << " " << t.second.count() << "\n";

    return;
  }

}
