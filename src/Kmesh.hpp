#include "Util.hpp"

#ifndef MESH_DEFINED
#define MESH_DEFINED



class kMesh {
  public:
   size_t Nk;
   std::vector<size_t> mesh;


   kMesh();
   kMesh(std::vector<size_t> _mesh);

   size_t backfold(const size_t i, const size_t N);
   std::vector<size_t> backfold(const std::vector<size_t> v);
   size_t kToIdx(const std::vector<size_t> kPoint);
   std::vector<size_t> idxToK(const size_t i);
   size_t idxMinusIdx(const size_t i, const size_t j);
   size_t getForthIdx(const size_t k, const size_t i, const size_t j);
   size_t getMinusIdx(const size_t i);
   std::vector< std::vector<size_t> > getNZC(const size_t d);
   void print();
};


#endif
