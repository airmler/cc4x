#include <string>
#include <iostream>

#include <ctf.hpp>

#include <cc4x.hpp>
#include <Read.hpp>
#include <Kmesh.hpp>
#include <Slice.hpp>
#include <Integrals.hpp>
#include <Drccd.hpp>

// this is a insane hack.once in the lifetime of the universe, we will fail
// ever tried. ever failed. no matter. try again. fail again. fail better.
#define NULL_TENSOR ((CTF::bsTensor<Complex>*)0xfafa)


int cc4x::No;
int cc4x::Nv;
CTF::World * cc4x::dw = NULL;
kMesh * cc4x::kmesh = NULL;

int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
 

  //cc4x::dw = new CTF::World("normal", 72);
  cc4x::dw = new CTF::World();

  CTF::bsTensor<Complex> *eps = NULL_TENSOR;
  CTF::bsTensor<Complex> *coulombVertex = NULL_TENSOR;

  CTF::bsTensor<Complex> *epsi = NULL_TENSOR;
  CTF::bsTensor<Complex> *epsa = NULL_TENSOR;

  CTF::bsTensor<Complex> *hhVertex = NULL_TENSOR;
  CTF::bsTensor<Complex> *phVertex = NULL_TENSOR;
  CTF::bsTensor<Complex> *hpVertex = NULL_TENSOR;
  CTF::bsTensor<Complex> *ppVertex = NULL_TENSOR;

  CTF::bsTensor<Complex> *Vpphh = NULL_TENSOR;
  CTF::bsTensor<Complex> *Vphhp = NULL_TENSOR;
  CTF::bsTensor<Complex> *Vhhpp = NULL_TENSOR;

  try {
    {
      std::cout << "read eigen" << std::endl;
      Read::input in({"EigenEnergies.yaml"});
      Read::output out({&eps});
      Read::run(in, out);
    }
    {
      std::cout << "read coulomb" << std::endl;
      Read::input in({"CoulombVertex.yaml"});
      Read::output out({&coulombVertex});
      Read::run(in, out);
    }
    {
      std::cout << "slice eps" << std::endl;
      Slice::input in({eps, {cc4x::No}});
      Slice::output out({&epsi, &epsa});
      Slice::run(in, out);
    }
    {
      std::cout << "slice ck" << std::endl;
      Slice::input in({coulombVertex, {0, cc4x::No, cc4x::No}});
      Slice::output out({&hhVertex, &phVertex, &hpVertex, &ppVertex});
      Slice::run(in, out);
    }
    {
      std::cout << "eval Integrals" << std::endl;
      Integrals::input in({hhVertex, phVertex, hpVertex, ppVertex});
      Integrals::output out({&Vpphh, &Vphhp, &Vhhpp});
      Integrals::run(in, out);
    }
    {
      std::cout << "drccd" << std::endl;
      Drccd::input in({Vpphh, Vphhp, Vhhpp, epsi, epsa});
      Drccd::output out({});
      Drccd::run(in, out);
    }
  } catch (...) {
    std::cout << "WHAT THE FUCK" << std::endl;
  }



  MPI_Finalize();
  return 0;
}
