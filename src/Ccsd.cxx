#include "Ccsd.hpp"
#include "Util.hpp"
#include "Scf.hpp"
#include "cc4x.hpp"

namespace Ccsd{
  Complex conjo(Complex x){ return std::conj(x); }
  std::function<int(const ivec &, const ivec &)>
  compare(const ivec p)
  {
    return [p] (const ivec &a, const ivec &b) -> int {
      size_t n(a.size());
      ivec c(n), d(n);
      for (size_t i(0); i < n - 1; i++){
        c[i] = a[p[i]]; d[i] = b[p[i]];
      }
      return c < d;
    };
  }



  void run(input const& in, output& out){
    Timings chrono;
    chrono["ccsd"].start();
    // sanity checks
    if ( in.Vhhhh == NULL    || in.Vhhhh == (tensor<Complex>*)0xfafa
      || in.Vhhhp == NULL    || in.Vhhhp == (tensor<Complex>*)0xfafa
      || in.Vhhph == NULL    || in.Vhhph == (tensor<Complex>*)0xfafa
      || in.Vhhpp == NULL    || in.Vhhpp == (tensor<Complex>*)0xfafa
      || in.Vphhh == NULL    || in.Vphhh == (tensor<Complex>*)0xfafa
      || in.Vphhp == NULL    || in.Vphhp == (tensor<Complex>*)0xfafa
      || in.Vphph == NULL    || in.Vphph == (tensor<Complex>*)0xfafa
      || in.Vpphh == NULL    || in.Vpphh == (tensor<Complex>*)0xfafa
      || in.Vphpp == NULL    || in.Vphpp == (tensor<Complex>*)0xfafa
      || in.epsi ==  NULL    || in.epsi  == (tensor<Complex>*)0xfafa
      || in.epsa ==  NULL    || in.epsa  == (tensor<Complex>*)0xfafa
      || in.hhVertex == NULL || in.hhVertex == (tensor<Complex>*)0xfafa
      || in.phVertex == NULL || in.phVertex == (tensor<Complex>*)0xfafa
      || in.hpVertex == NULL || in.hpVertex == (tensor<Complex>*)0xfafa
      || in.ppVertex == NULL || in.ppVertex == (tensor<Complex>*)0xfafa
       ) {
      THROW("Input of Ccsd not valid");
    }

		int64_t h(cc4x::No), p(cc4x::Nv), g(in.ppVertex->lens[0]);
    auto nzc4(cc4x::kmesh->getNZC(4));
    auto nzc3(cc4x::kmesh->getNZC(3));
    auto nzc2(cc4x::kmesh->getNZC(2));
    auto nzc1(cc4x::kmesh->getNZC(1));


    tensor<Complex> cThhVertex(3, {g,h,h}, nzc3, cc4x::dw, "cTGhh");
    tensor<Complex> cTphVertex(3, {g,p,h}, nzc3, cc4x::dw, "cTGph");
    tensor<Complex> cThpVertex(3, {g,h,p}, nzc3, cc4x::dw, "cTGhp");
    tensor<Complex> cTppVertex(3, {g,p,p}, nzc3, cc4x::dw, "cTGpp");

    tensor<Complex> Dai(2, {p,h}, nzc2, cc4x::dw, "Dai");
    tensor<Complex> Dabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Dabij");
    tensor<Complex> Rai(2, {p,h}, nzc2, cc4x::dw, "Rai");
    tensor<Complex> Rabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Rabij");
    tensor<Complex> Tai(2, {p,h}, nzc2, cc4x::dw, "Tai");
    tensor<Complex> Tabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Tabij");
    tensor<Complex> Kabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Kabij");


    tensor<Complex> Kac(2, {p,p}, nzc2, cc4x::dw, "Kac");
    tensor<Complex> Kki(2, {h,h}, nzc2, cc4x::dw, "Kki");
    tensor<Complex> Lac(2, {p,p}, nzc2, cc4x::dw, "Lac");
    tensor<Complex> Lki(2, {h,h}, nzc2, cc4x::dw, "Lki");
    tensor<Complex> Kck(2, {p,h}, nzc2, cc4x::dw, "Kck");
    tensor<Complex> Xabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Xabij");
    tensor<Complex> Yabij(4, {p,p,h,h}, nzc4, cc4x::dw, "Yabij");
    tensor<Complex> Xakic(4, {p,h,h,p}, nzc4, cc4x::dw, "Xakic");
    tensor<Complex> Xakci(4, {p,h,p,h}, nzc4, cc4x::dw, "Xakci");
    tensor<Complex> Xklij(4, {h,h,h,h}, nzc4, cc4x::dw, "Xklij");
    tensor<Complex> Xabcd(4, {p,p,p,p}, nzc4, cc4x::dw, "Xabij");
    tensor<Complex> Xakij(4, {p,h,h,h}, nzc4, cc4x::dw, "Xakij");
    tensor<Complex> Xabic(4, {p,p,h,p}, nzc4, cc4x::dw, "Xabic");

    tensor<Complex> Ghh(3, {g,h,h}, nzc3, cc4x::dw, "Ghh");
    tensor<Complex> Gph(3, {g,p,h}, nzc3, cc4x::dw, "Gph");
    tensor<Complex> Ghp(3, {g,h,p}, nzc3, cc4x::dw, "Ghp");
    tensor<Complex> Gpp(3, {g,p,p}, nzc3, cc4x::dw, "Gpp");
    tensor<Complex> G(1, {g}, nzc1, cc4x::dw, "G");

    std::function<Complex(const Complex)> fConj(&conjo);
    auto flip(in.hpVertex->nonZeroCondition);
    // sort in a way that the second column in the fastest, third second fastest
    // sorting like that introduces a flip from ia->ai of the nonZeroConditions
    std::sort(flip.begin(), flip.end(), compare({1,2,0}));

    cThhVertex.sum(1.0, *in.hhVertex, "gij", 1.0, "gji", in.hhVertex->nonZeroCondition, flip, fConj);
    cTphVertex.sum(1.0, *in.hpVertex, "gij", 1.0, "gji", in.hpVertex->nonZeroCondition, flip, fConj);
    cThpVertex.sum(1.0, *in.phVertex, "gij", 1.0, "gji", in.phVertex->nonZeroCondition, flip, fConj);
    cTppVertex.sum(1.0, *in.ppVertex, "gij", 1.0, "gji", in.ppVertex->nonZeroCondition, flip, fConj);

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


//TODO: The following five lines have to be replaced by the CCSDREF implementation
//    (because of the flawed logic in ctf-bs)
//      G.contract(1.0,   *in.hpVertex, "Gkd", Tai, "dk", 0.0, "G");
//      Lac.contract(2.0, cTppVertex, "Gac", G, "G", 1.0, "ac");
//      Gpp.contract(1.0, *in.hpVertex, "Gkc", Tai, "dk", 0.0, "Gcd");
//      Lac.contract(-1., cTppVertex, "Gad", Gpp, "Gcd", 1.0, "ac");
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

//TODO:The following five lines have to be replaced by the CCSDREF implementation
//    (because of the flawed logic in ctf-bs)
      //Gph.sum(1.0, cTphVertex, "Gai", 0.0, "Gai");
      //Gph.contract(-1.0, cThhVertex, "Gki", Tai, "ak", 1.0, "Gai");
      //Ghp.contract(1.0, *in.ppVertex, "Gbc", Tai, "cj", 0.0, "Gjb");
      //Rabij.contract(0.0, Gph, "Gai", Ghp, "Gjb", 1.0, "abij");

// Note: We use the coulomb integral available and conjugate it
      Xabic.sum(1.0, *in.Vphpp, "ciba", 0.0, "abic");
      Xabic.sum(1.0, Xabic, "abic", 0.0, "abic", fConj);
      Xabic.contract(-1.0, *in.Vphph, "bkci", Tai, "ak", 1.0, "abic");
      Rabij.contract(1.0, Xabic, "abic", Tai, "cj", 1.0, "abij");

      chrono["ccsd - t1r2"].stop();
      // Xakic
      chrono["ccsd - akic"].start();

//TODO:The following five lines have to be replaced by the CCSDREF implementation
//    (because of the flawed logic in ctf-bs)
//      Gph.contract( 1.0, cTppVertex, "Gad", Tai, "di", 1.0, "Gai");
//      Gph.contract(-0.5, cThpVertex, "Gld", Yabij, "dail", 1.0, "Gai");
//      Xakic.contract(1.0, Gph, "Gai", *in.hpVertex, "Gkc", 0.0, "akic");
//      Xakic.contract(1.0, *in.Vhhpp, "lkdc", Tabij, "adil", 1.0, "akic");
//      Xakic.contract(-0.5,*in.Vhhpp, "lkcd", Tabij, "adil", 1.0, "akic");
      Xakic.sum(1.0, *in.Vphhp, "akic", 0.0, "akic");
      Xakic.contract(-1.0, *in.Vhhhp, "lkic", Tai, "al", 1.0, "akic");
      Xakic.contract( 1.0, *in.Vphpp, "akdc", Tai, "di", 1.0, "akic");
      Xakic.contract(-0.5, *in.Vhhpp, "lkdc", Yabij, "dail", 1.0, "akic");
      Xakic.contract( 1.0, *in.Vhhpp, "lkdc", Tabij, "adil", 1.0, "akic");
      Xakic.contract(-0.5, *in.Vhhpp, "lkcd", Tabij, "adil", 1.0, "akic");

      Rabij.contract(2.0, Xakic, "akic", Tabij, "cbkj", 1.0, "abij");
      Rabij.contract(-1.0, Xakic, "akic", Tabij, "bckj", 1.0, "abij");
      chrono["ccsd - akic"].stop();

      // Xakci
      chrono["ccsd - akci"].start();
      Gpp.sum(1.0, cTppVertex, "Gab", 0.0, "Gab");
      Gpp.contract(-1.0, cThpVertex, "Glc", Tai, "al", 1.0, "Gac");
      Ghh.sum(1.0, *in.hhVertex, "Gij", 0.0, "Gij");
      Ghh.contract(1.0, *in.hpVertex, "Gkd", Tai, "di", 1.0, "Gki");
      Xakci.contract( 1.0, Gpp, "Gac", Ghh, "Gki", 0.0, "akci");
      Xakci.contract(-0.5, *in.Vhhpp, "lkcd", Tabij, "dail", 1.0, "akci");
      Rabij.contract(-1.0, Xakci, "akci", Tabij, "cbkj", 1.0, "abij");
      Rabij.contract(-1.0, Xakci, "bkci", Tabij, "ackj", 1.0, "abij");
      chrono["ccsd - akci"].stop();

      // Permutation and add V
      Kabij.sum(1.0, Rabij, "abij", 0.0, "baji");
      Rabij.sum(1.0, Kabij, "abij", 1.0, "abij");

      Rabij.sum(1.0, *in.Vpphh, "abij", 1.0, "abij");

      // Xabcd
      chrono["ccsd - pp-ladder"].start();
      //number of slices
      int nS(int(std::ceil(1.0*p/cc4x::Nx)));
      //prepare the dressed vertices in a vector
      std::vector<tensor<Complex> *> slcTppVertex(nS);
      std::vector<tensor<Complex> *> slppVertex(nS);

      for (int i(0); i < nS; i++){
        int iStart = i*cc4x::Nx, iEnd = std::min((i+1)*cc4x::Nx, (int) p);
        int y = iEnd - iStart;
        slcTppVertex[i] = new tensor<Complex>(3,{g,cc4x::Nx,p}, nzc3, cc4x::dw, "slicecT");

        slcTppVertex[i]->slice(
          {0,0,0}, {g,y,p}, 0.0, Gpp, {0,iStart,0}, {g,iEnd,p}, 1.0
        );
      }
      Gpp.sum(1.0, *in.ppVertex, "Gab", 0.0, "Gab");
      Gpp.contract(-1.0, *in.hpVertex, "Gkb", Tai, "ak", 1.0, "Gab");
      for (int i(0); i < nS; i++){
        int iStart = i*cc4x::Nx, iEnd = std::min((i+1)*cc4x::Nx, (int) p);
        int y = iEnd - iStart;
        slppVertex[i] = new tensor<Complex>(3,{g,cc4x::Nx,p}, nzc3, cc4x::dw, "slicecT");
        slppVertex[i]->slice(
          {0,0,0}, {g,y,p}, 0.0, Gpp, {0,iStart,0}, {g,iEnd,p}, 1.0
        );
      }
      for (int m(0); m < nS; m++)
      for (int n(m); n < nS; n++){
        //adapt the size of Vxycd and co
        int a(n*cc4x::Nx), b(m*cc4x::Nx);
        int y(slppVertex[m]->lens[1]), x(slcTppVertex[n]->lens[1]);
        tensor<Complex> Vxycd(4, {x,y,p,p}, nzc4, cc4x::dw, "Vxycd");
        tensor<Complex> Rxyij(4, {x,y,h,h}, nzc4, cc4x::dw, "Rxyij");
        tensor<Complex> Ryxji(4, {y,x,h,h}, nzc4, cc4x::dw, "Rxyji");
        Vxycd.contract(
          1.0, *slcTppVertex[n], "Gxc", *slppVertex[m], "Gyd", 0.0, "xycd"
        );
        Rxyij.contract(1.0, Vxycd, "xycd", Xabij, "cdij", 0.0, "xyij");
        Rabij.slice({a,b,0,0}, {a+x,b+y,h,h}, 1.0, Rxyij, {0,0,0,0}, {x,y,h,h}, 1.0);
        if (a > b){
          Ryxji.sum(1.0, Rxyij, "xyij", 0.0, "yxji");
          Rabij.slice({b,a,0,0}, {b+y,a+x,h,h}, 1.0, Ryxji, {0,0,0,0}, {x,y,h,h}, 1.0);
        }

      }

      chrono["ccsd - pp-ladder"].stop();


      // Xijkl
      chrono["ccsd - hh-ladder"].start();
      Xklij.sum(1.0, *in.Vhhhh, "klij", 0.0, "klij");
      Xklij.contract(1.0, *in.Vhhhp, "klic", Tai, "cj", 1.0, "klij");
      Xklij.contract(1.0, *in.Vhhhp, "lkjc", Tai, "ci", 1.0, "klij");
      Rabij.contract(1.0, Xklij, "klij", Xabij, "abkl", 1.0, "abij");
      Xklij.contract(1.0, *in.Vhhpp, "klcd", Xabij, "cdij", 0.0, "klij");
      Rabij.contract(1.0, Xklij, "klij", Tabij, "abkl", 1.0, "abij");

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

//TODO:The following four lines have to be replaced by the CCSDREF implementation
//    (because of the flawed logic in ctf-bs)
//      Gph.contract( 1.0, *in.hpVertex, "Gkd", Xabij, "cdik", 0.0, "Gci",true);
//      Rai.contract( 2.0, cTppVertex, "Gca", Gph, "Gci", 1.0, "ai",true);
//      Gph.contract( 1.0, *in.hpVertex, "Gkc", Xabij, "cdik", 0.0, "Gdi");
//      Rai.contract(-1.0, cTppVertex, "Gad", Gph, "Gdi", 1.0, "ai");
      Rai.contract( 2.0, *in.Vphpp, "akcd", Xabij, "cdik", 1.0, "ai");
      Rai.contract(-1.0, *in.Vphpp, "akdc", Xabij, "cdik", 1.0, "ai");

      Rai.contract(-2.0, *in.Vhhhp, "klic", Xabij, "ackl", 1.0, "ai");
      Rai.contract( 1.0, *in.Vhhhp, "lkic", Xabij, "ackl", 1.0, "ai");
      chrono["ccsd - singles"].stop();

      // add energy denominator
      Tai.contract(1.0, Dai, "ai", Rai, "ai", 0.0, "ai");
      Tabij.contract(1.0, Dabij, "abij", Rabij, "abij", 0.0, "abij");
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
