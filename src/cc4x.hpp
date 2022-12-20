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
   static size_t No;
   static size_t Nv;
   static size_t Nx;
   static size_t NF;
   static size_t iterations;
   static CTF::World *dw;
   static kMesh *kmesh;
};



#endif
