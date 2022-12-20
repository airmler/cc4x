#include "Integrals.hpp"
#include "Util.hpp"
#include <cc4x.hpp>

namespace Integrals{
  Complex conjo(Complex x){ return std::conj(x); }

  std::function<bool(const std::vector<size_t> &, const std::vector<size_t> &)>
  _compare(const std::vector<size_t> p)
  {
    return [p] (const std::vector<size_t> &a, const std::vector<size_t> &b) -> bool {
      size_t n(p.size());
      std::vector<size_t> c(n), d(n);
      for (size_t i(0); i < n; i++){
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
    auto flip(in.hp->nzc);
    std::vector<size_t> remap(in.hp->nzc.size());

    // sort in a way that the second column in the fastest, third second fastest
    // sorting like that introduces a flip from ia->ai of the nonZeroConditions
    std::sort(flip.begin(), flip.end(), _compare({1,2,0}));
    for (size_t i(0); i < nzc3.size(); i++) remap[i] = flip[i].back();


    chrono["Integrals - pre"].start();
    conjGhp.sum(1.0, *in.ph, "Gai", 0.0, "Gia", remap, fConj);
    conjGph.sum(1.0, *in.hp, "Gia", 0.0, "Gai", remap, fConj);
    conjGhh.sum(1.0, *in.hh, "Gia", 0.0, "Gai", remap, fConj);
    conjGpp.sum(1.0, *in.pp, "Gia", 0.0, "Gai", remap, fConj);

    auto minusGhh(conjGhp.nzc);
    for (auto &e: minusGhh) e[0] = cc4x::kmesh->getMinusIdx(e[0]);
    conjGhh.relabelBlocks(minusGhh);
    conjGhp.relabelBlocks(minusGhh);
    conjGph.relabelBlocks(minusGhh);
    conjGpp.relabelBlocks(minusGhh);

    chrono["Integrals - pre"].stop();
    size_t h(cc4x::No), p(cc4x::Nv);
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

    if (cc4x::ccsd) {
      auto Vhhhh = new tensor<Complex>(4, {h,h,h,h}, nzc4, cc4x::dw, "Vhhhh");
      auto Vhhhp = new tensor<Complex>(4, {h,h,h,p}, nzc4, cc4x::dw, "Vhhhp");
      auto Vhhph = new tensor<Complex>(4, {h,h,p,h}, nzc4, cc4x::dw, "Vhhph");
      auto Vhphp = new tensor<Complex>(4, {h,p,h,p}, nzc4, cc4x::dw, "Vhphp");
      auto Vphhh = new tensor<Complex>(4, {p,h,h,h}, nzc4, cc4x::dw, "Vphhh");
      auto Vphph = new tensor<Complex>(4, {p,h,p,h}, nzc4, cc4x::dw, "Vphph");

      Vhhhh->contract(1.0, conjGhh, "Gai", *in.hh, "Gbj", 0.0, "abij");
      Vhhhp->contract(1.0, conjGhh, "Gaj", *in.hp, "Gib", 0.0, "aijb");
      Vhhph->contract(1.0, conjGhp, "Gia", *in.hh, "Gjb", 0.0, "ijab");
      Vhphp->contract(1.0, conjGhh, "Gai", *in.pp, "Gbj", 0.0, "abij");
      Vphhh->contract(1.0, conjGph, "Gaj", *in.hh, "Gib", 0.0, "aijb");
      Vphph->contract(1.0, conjGpp, "Gaj", *in.hh, "Gib", 0.0, "aijb");

      *out.Vhhhh = Vhhhh;
      *out.Vhhhp = Vhhhp;
      *out.Vhhph = Vhhph;
      *out.Vhphp = Vhphp;
      *out.Vphhh = Vphhh;
      *out.Vphph = Vphph;
    }
    if (cc4x::ref) {
      auto Vhppp = new tensor<Complex>(4, {h,p,p,p}, nzc4, cc4x::dw, "Vhppp");
      auto Vpphp = new tensor<Complex>(4, {p,p,h,p}, nzc4, cc4x::dw, "Vpphp");
      auto Vppph = new tensor<Complex>(4, {p,p,p,h}, nzc4, cc4x::dw, "Vppph");
      auto Vpppp = new tensor<Complex>(4, {p,p,p,p}, nzc4, cc4x::dw, "Vpppp");
      auto Vphpp = new tensor<Complex>(4, {p,h,p,p}, nzc4, cc4x::dw, "Vphpp");

      Vhppp->contract(1.0, conjGhp, "Gai", *in.pp, "Gbj", 0.0, "abij");
      Vpphp->contract(1.0, conjGph, "Gai", *in.pp, "Gbj", 0.0, "abij");
      Vppph->contract(1.0, conjGpp, "Gai", *in.ph, "Gbj", 0.0, "abij");
      Vphpp->contract(1.0, conjGpp, "Gac", *in.hp, "Gkd", 0.0, "akcd");
      Vpppp->contract(1.0, conjGpp, "Gaj", *in.pp, "Gib", 0.0, "aijb");

      *out.Vhppp = Vhppp;
      *out.Vpphp = Vpphp;
      *out.Vppph = Vppph;
      *out.Vphpp = Vphpp;
      *out.Vpppp = Vpppp;
    }
    chrono["Integrals - Ccsd"].stop();

    for (auto t: chrono)
      LOG() << t.first << " " << t.second.count() << "\n";


    return;
  }

}
