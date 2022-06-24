#include "Integrals.hpp"
#include "Util.hpp"
#include <cc4x.hpp>

namespace Integrals{

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
    // sanity checks
    if ( in.pp == NULL || in.pp == (tensor<Complex>*)0xfafa
      || in.ph == NULL || in.ph == (tensor<Complex>*)0xfafa
      || in.hp == NULL || in.hp == (tensor<Complex>*)0xfafa
      || in.hh == NULL || in.hh == (tensor<Complex>*)0xfafa
       ) {
      THROW("Input of Slice not valid");
    }
    auto gph(in.ph->lens);
    auto ghp(in.hp->lens);
    auto ghh(in.hh->lens);
    auto gpp(in.pp->lens);
    auto nzc3(cc4x::kmesh->getNZC(3));
    tensor<Complex> conjGph(3, gph, nzc3, cc4x::dw, "cGph");
    tensor<Complex> conjGhp(3, ghp, nzc3, cc4x::dw, "cGhp");
    tensor<Complex> conjGhh(3, ghh, nzc3, cc4x::dw, "cGhh");
    tensor<Complex> conjGpp(3, gpp, nzc3, cc4x::dw, "cGpp");

    std::function<Complex(const Complex)> fConj(&conjo);
    auto flip(in.hp->nonZeroCondition);

    // sort in a way that the second column in the fastest, third second fastest
    // sorting like that introduces a flip from ia->ai of the nonZeroConditions
    std::sort(flip.begin(), flip.end(), compare({1,2,0}));

    chrono["Integrals - pre"].start();
    conjGhp.sum(1.0, *in.ph, "Gai", 0.0, "Gia", in.ph->nonZeroCondition, flip, fConj);
    conjGph.sum(1.0, *in.hp, "Gia", 0.0, "Gai", in.hp->nonZeroCondition, flip, fConj);
    conjGhh.sum(1.0, *in.hh, "Gia", 0.0, "Gai", in.hh->nonZeroCondition, flip, fConj);
    conjGpp.sum(1.0, *in.pp, "Gia", 0.0, "Gai", in.pp->nonZeroCondition, flip, fConj);
    chrono["Integrals - pre"].stop();
    int64_t h(cc4x::No), p(cc4x::Nv);
    auto nzc4(cc4x::kmesh->getNZC(4));
    auto Vpphh = new tensor<Complex>(4, {p,p,h,h}, nzc4, cc4x::dw, "Vpphh");
    auto Vphhp = new tensor<Complex>(4, {p,h,h,p}, nzc4, cc4x::dw, "Vphhp");
    auto Vhhpp = new tensor<Complex>(4, {h,h,p,p}, nzc4, cc4x::dw, "Vhhpp");

    chrono["Integrals - Ccsd"].start();
    chrono["Integrals - Drccd"].start();
    Vpphh->contract(1.0, conjGph, "Gai", *in.ph, "Gbj", 0.0, "abij");
    Vphhp->contract(1.0, conjGph, "Gaj", *in.hp, "Gib", 0.0, "aijb");
    Vhhpp->contract(1.0, conjGhp, "Gia", *in.hp, "Gjb", 0.0, "ijab");
    chrono["Integrals - Drccd"].stop();

    *out.Vpphh = Vpphh;
    *out.Vphhp = Vphhp;
    *out.Vhhpp = Vhhpp;

    auto Vhhhh = new tensor<Complex>(4, {h,h,h,h}, nzc4, cc4x::dw, "Vhhhh");
    auto Vhhhp = new tensor<Complex>(4, {h,h,h,p}, nzc4, cc4x::dw, "Vhhhp");
    auto Vhhph = new tensor<Complex>(4, {h,h,p,h}, nzc4, cc4x::dw, "Vhhph");
    auto Vhphp = new tensor<Complex>(4, {h,p,h,p}, nzc4, cc4x::dw, "Vhphp");
    auto Vhppp = new tensor<Complex>(4, {h,p,p,p}, nzc4, cc4x::dw, "Vhppp");
    auto Vphhh = new tensor<Complex>(4, {p,h,h,h}, nzc4, cc4x::dw, "Vphhh");
    auto Vphph = new tensor<Complex>(4, {p,h,p,h}, nzc4, cc4x::dw, "Vphph");
    auto Vphpp = new tensor<Complex>(4, {p,h,p,p}, nzc4, cc4x::dw, "Vphpp");
    auto Vppph = new tensor<Complex>(4, {p,p,p,h}, nzc4, cc4x::dw, "Vppph");
    auto Vpphp = new tensor<Complex>(4, {p,p,h,p}, nzc4, cc4x::dw, "Vpphp");
    auto Vpppp = new tensor<Complex>(4, {p,p,p,p}, nzc4, cc4x::dw, "Vpppp");

    Vhhhh->contract(1.0, conjGhh, "Gai", *in.hh, "Gbj", 0.0, "abij");
    Vhhhp->contract(1.0, conjGhh, "Gaj", *in.hp, "Gib", 0.0, "aijb");
    Vhhph->contract(1.0, conjGhp, "Gia", *in.hh, "Gjb", 0.0, "ijab");
    Vhphp->contract(1.0, conjGhh, "Gai", *in.pp, "Gbj", 0.0, "abij");
    Vhppp->contract(1.0, conjGhp, "Gai", *in.pp, "Gbj", 0.0, "abij");
    Vphhh->contract(1.0, conjGph, "Gaj", *in.hh, "Gib", 0.0, "aijb");
    Vphph->contract(1.0, conjGpp, "Gaj", *in.hh, "Gib", 0.0, "aijb");
    Vphpp->contract(1.0, conjGpp, "Gia", *in.hp, "Gjb", 0.0, "ijab");
    Vpphp->contract(1.0, conjGph, "Gai", *in.pp, "Gbj", 0.0, "abij");
    Vppph->contract(1.0, conjGpp, "Gai", *in.ph, "Gbj", 0.0, "abij");
    chrono["Integrals - pppp"].start();
    Vpppp->contract(1.0, conjGpp, "Gaj", *in.pp, "Gib", 0.0, "aijb");
    chrono["Integrals - Ccsd"].stop();
    chrono["Integrals - pppp"].stop();


    for (auto t: chrono)
      LOG() << t.first << " " << t.second.count() << "\n";

    *out.Vhhhh = Vhhhh;
    *out.Vhhhp = Vhhhp;
    *out.Vhhph = Vhhph;
    *out.Vhphp = Vhphp;
    *out.Vhppp = Vhppp;
    *out.Vphhh = Vphhh;
    *out.Vphph = Vphph;
    *out.Vphpp = Vphpp;
    *out.Vpphp = Vpphp;
    *out.Vppph = Vppph;
    *out.Vpppp = Vpppp;



    return;
  }

}
