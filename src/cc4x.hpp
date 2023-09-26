#include <Kmesh.hpp>

#ifndef CC4X_DEFINED
#  define CC4X_DEFINED

class cc4x {
public:
  static bool verbose;
  static bool complexT;
  static bool drccd;
  static bool ccsd;
  static bool ref;
  static int64_t No;
  static int64_t Nv;
  static int64_t Nx;
  static int64_t NF;
  static int64_t iterations;
  static World *world;
  static kMesh *kmesh;
};

#endif
