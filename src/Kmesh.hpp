#include "Util.hpp"

#ifndef MESH_DEFINED
#define MESH_DEFINED



class kMesh {
  public:
   int Nk;
   ivec mesh;


   kMesh();
   kMesh(ivec _mesh);

   int backfold(const int i, const int N);
   ivec backfold(const ivec v);
   int kToIdx(const ivec kPoint);
   ivec idxToK(const int i);
   int idxMinusIdx(const int i, const int j);
   int getForthIdx(const int k, const int i, const int j);
   std::vector<ivec> getNZC(const int d);
   void print();
};


#endif
