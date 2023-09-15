#include <Kmesh.hpp>
#include <iostream>
#include <cassert>
#include <numeric>

  kMesh::kMesh(){};

  kMesh::kMesh(std::vector<int64_t> _mesh){
    mesh.resize(_mesh.size());
    mesh = _mesh;
    Nk = std::accumulate(mesh.begin(), mesh.end(), 1L, std::multiplies<int64_t>());
  }

  int64_t kMesh::backfold(const int64_t i, const int64_t N) { return (i+4*N)%N; }
  std::vector<int64_t> kMesh::backfold(const std::vector<int64_t> in){
    return { backfold(in[0], mesh[0])
           , backfold(in[1], mesh[1])
           , backfold(in[2], mesh[2])
           };
  }
  int64_t kMesh::kToIdx(const std::vector<int64_t> kPoint){
    return backfold(kPoint[0], mesh[0])
         + backfold(kPoint[1], mesh[1])*mesh[0]
         + backfold(kPoint[2], mesh[2])*mesh[0]*mesh[1];
  }
  std::vector<int64_t> kMesh::idxToK(const int64_t i){
    return {  i%mesh[0]
           , (i/mesh[0])%mesh[1]
           ,  i/mesh[0]/ mesh[1]
           };
  }
  int64_t kMesh::idxMinusIdx(const int64_t i, const int64_t j){
    std::vector<int64_t> ki(idxToK(i));
    std::vector<int64_t> kj(idxToK(j));
    std::vector<int64_t> d(
      {  ( kj[0] - ki[0] + mesh[0])
      ,  ( kj[1] - ki[1] + mesh[1])
      ,  ( kj[2] - ki[2] + mesh[2])
      }                  );
    return kToIdx(backfold(d));
  }
  int64_t kMesh::getForthIdx(const int64_t k, const int64_t i, const int64_t j){
    // kk + ko = ki + kj ==>  ko = ki + kj - kk
    std::vector<int64_t> ki(idxToK(i));
    std::vector<int64_t> kj(idxToK(j));
    std::vector<int64_t> kk(idxToK(k));
    std::vector<int64_t> d(
      { ( ki[0] + kj[0] - kk[0] + mesh[0] )
      , ( ki[1] + kj[1] - kk[1] + mesh[1] )
      , ( ki[2] + kj[2] - kk[2] + mesh[2] )
      }                  );
    return kToIdx(backfold(d));
  }
  int64_t kMesh::getMinusIdx(const int64_t i){
    std::vector<int64_t> k(idxToK(i));
    std::vector<int64_t> d( { -k[0] + mesh[0]
                           ,  -k[1] + mesh[1]
                           ,  -k[2] + mesh[2]
                           });
    return kToIdx(backfold(d));
  }
  std::vector< std::vector<int64_t> > kMesh::getNZC(const int64_t dim){
    if ( dim == 0) {
      std::vector< std::vector<int64_t> > out(1, std::vector<int64_t>(1));
      out[0] = {0};
      return out;
    }
    if ( dim == 1) {
      std::vector< std::vector<int64_t> > out(Nk, std::vector<int64_t>(2));
      for (int64_t o(0); o < Nk; o++) out[o] = {o, o};
      return out;
    }
    if ( dim == 2) {
      std::vector< std::vector<int64_t> > out(Nk, std::vector<int64_t>(3));
      for (int64_t o(0); o < Nk; o++) out[o] = {o, o, o};
      return out;
    }
    if ( dim == 3) {
      std::vector< std::vector<int64_t> > out(Nk*Nk, std::vector<int64_t>(4));
      int64_t it(0);
      for (int64_t m(0); m < Nk; m++)
      for (int64_t n(0); n < Nk; n++)
        out[n + m * Nk] = {idxMinusIdx(n,m), n, m, it++};
        //out[n + m * Nk] = {-1, n, m};
      return out;
    }
    if ( dim == 4) {
      std::vector< std::vector<int64_t> > out(Nk*Nk*Nk, std::vector<int64_t>(5));
      int64_t it(0);
      for (int64_t m(0); m < Nk; m++)
      for (int64_t n(0); n < Nk; n++)
      for (int64_t o(0); o < Nk; o++)
        out[o + n * Nk + m * Nk * Nk] = {m, n, o, getForthIdx(o,n,m), it++};
      return out;
    }
    assert(0);
    return {};
  }
  void kMesh::print(){
    std::cout << "Print the k-mesh indices:" << std::endl;
    int cnt(0);
    for (int64_t i(0); i < mesh[0]; i++)
    for (int64_t j(0); j < mesh[1]; j++)
    for (int64_t k(0); k < mesh[2]; k++)
      std::cout << cnt++ << ": " << i << " " << j << " " << k << std::endl;
    std::cout << "=====================\n";
  }


