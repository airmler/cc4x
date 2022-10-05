#include <ctf.hpp>
#include <Kmesh.hpp>

#ifndef CC4X_DEFINED
#define CC4X_DEFINED

class cc4x {
  public: 
   static bool verbose;
   static bool complexT;
   static bool drccd;
   static bool ccsd;
   static bool ref;
   static int No;
   static int Nv;
   static int Nx;
   static int NF;
   static int iterations;
   static CTF::World *dw;
   static kMesh *kmesh;
};



#endif
