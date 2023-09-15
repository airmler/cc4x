#include "Ueg.hpp"
#include <math.h>
#include "Util.hpp"
#include <cc4x.hpp>



namespace Ueg{

  using iarr  = array<int64_t,3>;
  using darr  = array<double,4>;



  double evalMadelung(const double v){
    double kappa = pow(v,-1.0/3.0);
    double term2 = M_PI / (kappa*kappa*v);
    double term4 = 2 * kappa/sqrt(M_PI);
    double boxLength = 1.0/kappa;
    double recipsum = 0.0;
    double realsum = 0.0;
    for (int64_t l1=-6; l1 <= 6; ++l1)
    for (int64_t l2=-6; l2 <= 6; ++l2)
    for (int64_t l3=-6; l3 <= 6; ++l3){
      int64_t n2 = l1*l1 + l2*l2 + l3*l3;
      double modr = boxLength * sqrt((double)n2);
      double k2 = kappa*kappa*n2;
      if (n2 > 0){
       recipsum -= 1.0/(M_PI*k2)*exp(-M_PI*M_PI*k2/kappa/kappa)/v;
       realsum -= erfc(kappa*modr)/modr;
      }
    }
    return realsum + term2 + term4 + recipsum;
  }


  int64_t sL(const iarr a) {  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
  double sL(const darr a){    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
  double Vijji(const darr a, const darr b, const double v){
    darr q({a[0]-b[0], a[1]-b[1], a[2]-b[2]});
    if ( sL(q) < 1e-8 ) return evalMadelung(v);
    return 4.0*M_PI/v/sL(q);
  }

  void run(input const& in, output& out){
    auto No(in.No), Nv(in.Nv), Np(No+Nv);
    auto rs(in.rs); auto NF(in.NF);
    cc4x::No = No;
    cc4x::Nv = Nv;
    cc4x::complexT = true;
    int64_t maxG = std::pow(5.0*Np,1.0/3.0);
    std::vector<iarr> iGrid;
    for (int64_t g1(-maxG); g1 <= maxG; g1++)
    for (int64_t g2(-maxG); g2 <= maxG; g2++)
    for (int64_t g3(-maxG); g3 <= maxG; g3++)
      iGrid.push_back({g1, g2, g3});

    sort(iGrid.begin(), iGrid.end(), [](iarr a, iarr b){ return sL(a) < sL(b); });
    if (iGrid.size() < Np ) LOG() <<"BUG related to Np & maxG\n";
    if (sL(iGrid[No]) == sL(iGrid[No-1]))
      LOG() << "WARNING: occupied orbitals form not a closed shell\n";
    if (sL(iGrid[Np]) == sL(iGrid[Np-1]))
      LOG() << "WARNING: virtual orbitals form not a closed shell\n";
    iGrid.resize(Np);

    // define volume, lattice Constant, and reciprocal lattice constant
    double v(rs*rs*rs/3.0*4.0*M_PI*No*2);
    double a(pow(v,1./3.));
    double b(2.0*M_PI/a);

    std::vector<darr> dGrid;
    // here we can introduce a possible shift of the mesh
    for (auto i: iGrid)
      dGrid.push_back( { b*i[0], b*i[1], b*i[2], 0.0} );

    // now we can write the hartree fock energy in the 4th entry
    for (auto &d: dGrid){
      d[3] = 0.5*sL(d); // add the kinetic energy
      double exchE(0.0);
      for (int64_t o(0); o < No; o++)
        exchE += Vijji(d, dGrid[o], v);
      d[3] -= exchE;
    }
    double refE(0.0);
    for (int64_t o(0); o < No; o++) {
      refE += dGrid[o][3];
      refE += 0.5*sL(dGrid[o]);
    }

    LOG() << "Hartree Fock energy: " << refE/No/2 << " Ha per electron\n";
    LOG() << "HOMO: " << dGrid[No-1][3] << " , LUMO: " << dGrid[No][3] << '\n';
    // work on the eigen energies
    // We have to perform a hack here and make the energy vector larger
    auto _Nk(cc4x::kmesh->Nk);
    std::vector<Complex> energies(Np*_Nk);
    auto eps = new tensor<Complex>(1, {Np}, cc4x::kmesh->getNZC(1), cc4x::world, "eps");

    for (int64_t k(0); k < _Nk; k++)
    for (int64_t d(0); d < Np; d++)
      energies[d + k*Np] = {dGrid[d][3], 0.0};

    std::vector<int64_t> idx(Np);
    if (!cc4x::world->rank) std::iota(idx.begin(), idx.end(), 0);
    if (cc4x::world->rank)  eps->write(0,  idx, energies);
    else                 eps->write(Np, idx, energies);

    *out.eps = eps;

    // work on the coulomb vertex
    iarr maxMom({0,0,0});
    for (int64_t p(0); p < Np; p++)
    for (int64_t q(0); q < Np; q++){
      iarr d = { iGrid[p][0] - iGrid[q][0]
               , iGrid[p][1] - iGrid[q][1]
               , iGrid[p][2] - iGrid[q][2]
               };
      maxMom = std::max(maxMom, d, [](iarr a, iarr b) { return sL(a) < sL(b);});
    }

    int64_t maxR = sL(maxMom);
    maxG = max( {maxMom[0], maxMom[1], maxMom[2]}
              , [](int64_t a, int64_t b){ return std::abs(a) < std::abs(b);});
    maxG = std::abs(maxG);
    std::map<iarr,int64_t> momMap;
    int64_t index(0);
    for (int64_t g1(-maxG); g1 <= maxG; g1++)
    for (int64_t g2(-maxG); g2 <= maxG; g2++)
    for (int64_t g3(-maxG); g3 <= maxG; g3++){
      iarr t({g1,g2,g3});
      if ( sL(t) > maxR ) continue;
      momMap[t] = index++;
    }
    if (NF == 0) {
      NF = momMap.size();
      cc4x::NF = NF;
    }

    auto cV = new tensor<Complex>(3, {NF,Np,Np}, cc4x::kmesh->getNZC(3), cc4x::world, "cVertex");


    // Writing CoulombVertex to buffer
    // We have to do it mpi-able...otherwise we will
    // not be able to write it to a ctf tensor
    double fac(4.0*M_PI/v);

    int64_t np = cc4x::world->np;
    int64_t rank = cc4x::world->rank;
    // We slice the number of states for all the mpi processes
    int64_t slices(Np/np);
    std::vector<size_t> slicePerRank(np);
    for (size_t r(0); r < np; r++){
      size_t lslice(slices);
      for (size_t i(0); i < Np - slices*np; i++) if (r == i){
        lslice++;
      }
      slicePerRank[r] = lslice;
    }
    slices = slicePerRank[rank];
    //allocate only a buffer of needed size
    std::vector<Complex> vData(NF*Np*slices,{0,0});
    // determine begin and end of the rank's slices
    auto sbegin( std::accumulate( slicePerRank.begin()
                                , slicePerRank.begin() + rank
                                , 0UL
                                , std::plus<int64_t>()
                                )
               );

    for (int64_t s(0); s < slices; s++)
    for (int64_t q(0); q < Np; q++){
      auto p(s+sbegin);
      iarr d = { iGrid[q][0] - iGrid[p][0]
               , iGrid[q][1] - iGrid[p][1]
               , iGrid[q][2] - iGrid[p][2]
               };
      // This is a hack!
      // If NF is chosen by the user we will not have an overflow
      int64_t ii = momMap[d] % NF;
      double res;
      (sL(d)) ? res = fac/( b*b*sL(d) ) : res = evalMadelung(v);
      vData[ii+q*NF+s*NF*Np] = { sqrt(res), 0.0};
    }

    idx.resize(vData.size());
    std::iota(idx.begin(), idx.end(), sbegin*Np*NF);
    for (auto k(0); k < _Nk; k++)
      cV->write(idx.size(), idx, vData, k);


    *out.coulombVertex = cV;


    return;
  }

}
