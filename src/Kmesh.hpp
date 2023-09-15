#include "Util.hpp"

#ifndef MESH_DEFINED
#define MESH_DEFINED



class kMesh {
  public:
   int64_t Nk;
   std::vector<int64_t> mesh;


   kMesh();
   kMesh(std::vector<int64_t> _mesh);

   int64_t backfold(const int64_t i, const int64_t N);
   std::vector<int64_t> backfold(const std::vector<int64_t> v);
   int64_t kToIdx(const std::vector<int64_t> kPoint);
   std::vector<int64_t> idxToK(const int64_t i);
   int64_t idxMinusIdx(const int64_t i, const int64_t j);
   int64_t getForthIdx(const int64_t k, const int64_t i, const int64_t j);
   int64_t getMinusIdx(const int64_t i);
   std::vector< std::vector<int64_t> > getNZC(const int64_t d);
   void print();
};


#endif
