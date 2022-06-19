#include "Integrals.hpp"
#include "Util.hpp"
#include <cc4x.hpp>

namespace Integrals{

  Complex conjo(Complex x){ return std::conj(x); }
  Complex   inv(Complex x){ return 1.0/x; }

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
    // sanity checks
    if ( in.pp == NULL || in.pp == (CTF::bsTensor<Complex>*)0xfafa
      || in.ph == NULL || in.ph == (CTF::bsTensor<Complex>*)0xfafa
      || in.hp == NULL || in.hp == (CTF::bsTensor<Complex>*)0xfafa
      || in.hh == NULL || in.hh == (CTF::bsTensor<Complex>*)0xfafa
       ) {
      THROW("Input of Slice not valid");
    }
    auto gph(in.ph->lens);
    auto ghp(in.hp->lens);
    auto ghh(in.hh->lens);
    auto gpp(in.pp->lens);
    auto nzc3(cc4x::kmesh->getNZC(3));
    CTF::bsTensor<Complex> conjGph(3, gph, nzc3, cc4x::dw, "cGph");
    CTF::bsTensor<Complex> conjGhp(3, ghp, nzc3, cc4x::dw, "cGhp");
    CTF::bsTensor<Complex> conjGhh(3, ghh, nzc3, cc4x::dw, "cGhh");
    CTF::bsTensor<Complex> conjGpp(3, gpp, nzc3, cc4x::dw, "cGpp");


    CTF::Univar_Function<Complex> fConj(&conjo);
    auto flip(in.hp->nonZeroCondition);

    // sort in a way that the second column in the fastest, third second fastest
    // sorting like that introduces a flip from ia->ai of the nonZeroConditions
    std::sort(flip.begin(), flip.end(), compare({1,2,0}));
    conjGhp.sum(1.0, *in.ph, "Gai", 0.0, "Gia", in.ph->nonZeroCondition, flip, fConj);
    conjGph.sum(1.0, *in.hp, "Gia", 0.0, "Gai", in.hp->nonZeroCondition, flip, fConj);
    conjGhh.sum(1.0, *in.hh, "Gia", 0.0, "Gai", in.hh->nonZeroCondition, flip, fConj);
    conjGpp.sum(1.0, *in.pp, "Gia", 0.0, "Gai", in.pp->nonZeroCondition, flip, fConj);
    int64_t h(cc4x::No), p(cc4x::Nv);
    auto nzc4(cc4x::kmesh->getNZC(4));
    auto Vpphh = new CTF::bsTensor<Complex>(4, {p,p,h,h}, nzc4, cc4x::dw, "Vpphh");
    auto Vphhp = new CTF::bsTensor<Complex>(4, {p,h,h,p}, nzc4, cc4x::dw, "Vphhp");
		auto Vhhpp = new CTF::bsTensor<Complex>(4, {h,h,p,p}, nzc4, cc4x::dw, "Vhhpp");

    Vpphh->contract(1.0, conjGph, "Gai", *in.ph, "Gbj", 0.0, "abij");
    Vphhp->contract(1.0, conjGph, "Gaj", *in.hp, "Gib", 0.0, "aijb");
    Vhhpp->contract(1.0, conjGhp, "Gia", *in.hp, "Gjb", 0.0, "ijab");

    *out.Vpphh = Vpphh;
    *out.Vphhp = Vphhp;
    *out.Vhhpp = Vhhpp;

    auto Vhhhh = new CTF::bsTensor<Complex>(4, {h,h,h,h}, nzc4, cc4x::dw, "Vhhhh");
    auto Vhhhp = new CTF::bsTensor<Complex>(4, {h,h,h,p}, nzc4, cc4x::dw, "Vhhhp");
		auto Vhhph = new CTF::bsTensor<Complex>(4, {h,h,p,h}, nzc4, cc4x::dw, "Vhhph");
    auto Vhphp = new CTF::bsTensor<Complex>(4, {h,p,h,p}, nzc4, cc4x::dw, "Vhphp");
    auto Vhppp = new CTF::bsTensor<Complex>(4, {h,p,p,p}, nzc4, cc4x::dw, "Vhppp");
    auto Vphhh = new CTF::bsTensor<Complex>(4, {p,h,h,h}, nzc4, cc4x::dw, "Vphhh");
    auto Vphph = new CTF::bsTensor<Complex>(4, {p,h,p,h}, nzc4, cc4x::dw, "Vphph");
		auto Vphpp = new CTF::bsTensor<Complex>(4, {p,h,p,p}, nzc4, cc4x::dw, "Vphpp");
		auto Vppph = new CTF::bsTensor<Complex>(4, {p,p,p,h}, nzc4, cc4x::dw, "Vppph");
		auto Vpphp = new CTF::bsTensor<Complex>(4, {p,p,h,p}, nzc4, cc4x::dw, "Vpphp");
		auto Vpppp = new CTF::bsTensor<Complex>(4, {p,p,p,p}, nzc4, cc4x::dw, "Vpppp");

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
    Vpppp->contract(1.0, conjGpp, "Gaj", *in.pp, "Gib", 0.0, "aijb");

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
