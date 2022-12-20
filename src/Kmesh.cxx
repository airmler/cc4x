#include <Kmesh.hpp>
#include <iostream>
#include <cassert>
#include <numeric>

  kMesh::kMesh(){};

  kMesh::kMesh(std::vector<size_t> _mesh){
    mesh.resize(_mesh.size());
    mesh = _mesh;
    Nk = std::accumulate(mesh.begin(), mesh.end(), 1UL, std::multiplies<size_t>());
  }

  size_t kMesh::backfold(const size_t i, const size_t N) { return (i+4*N)%N; }
  std::vector<size_t> kMesh::backfold(const std::vector<size_t> in){
    return { backfold(in[0], mesh[0])
           , backfold(in[1], mesh[1])
           , backfold(in[2], mesh[2])
           };
  }
  size_t kMesh::kToIdx(const std::vector<size_t> kPoint){
    return backfold(kPoint[0], mesh[0])
         + backfold(kPoint[1], mesh[1])*mesh[0]
         + backfold(kPoint[2], mesh[2])*mesh[0]*mesh[1];
  }
  std::vector<size_t> kMesh::idxToK(const size_t i){
    return {  i%mesh[0]
           , (i/mesh[0])%mesh[1]
           ,  i/mesh[0]/ mesh[1]
           };
  }
  size_t kMesh::idxMinusIdx(const size_t i, const size_t j){
    std::vector<size_t> ki(idxToK(i));
    std::vector<size_t> kj(idxToK(j));
    std::vector<size_t> d(
      { (size_t) ((int64_t) kj[0] - (int64_t) ki[0] + (int64_t) mesh[0])
      , (size_t) ((int64_t) kj[1] - (int64_t) ki[1] + (int64_t) mesh[1])
      , (size_t) ((int64_t) kj[2] - (int64_t) ki[2] + (int64_t) mesh[2])
      }                  );
    return kToIdx(backfold(d));
  }
  size_t kMesh::getForthIdx(const size_t k, const size_t i, const size_t j){
    // kk + ko = ki + kj ==>  ko = ki + kj - kk
    std::vector<size_t> ki(idxToK(i));
    std::vector<size_t> kj(idxToK(j));
    std::vector<size_t> kk(idxToK(k));
    std::vector<size_t> d(
      { (size_t)( (int64_t) ki[0] + (int64_t) kj[0] - (int64_t) kk[0] + (int64_t) mesh[0] )
      , (size_t)( (int64_t) ki[1] + (int64_t) kj[1] - (int64_t) kk[1] + (int64_t) mesh[1] )
      , (size_t)( (int64_t) ki[2] + (int64_t) kj[2] - (int64_t) kk[2] + (int64_t) mesh[2] )
      }                  );
    return kToIdx(backfold(d));
  }
  size_t kMesh::getMinusIdx(const size_t i){
    std::vector<size_t> k(idxToK(i));
    std::vector<size_t> d( { (size_t) (- (int64_t) k[0] + (int64_t) mesh[0])
                           , (size_t) (- (int64_t) k[1] + (int64_t) mesh[1])
                           , (size_t) (- (int64_t) k[2] + (int64_t) mesh[2])
                           });
    return kToIdx(backfold(d));
  }
  std::vector< std::vector<size_t> > kMesh::getNZC(const size_t dim){
    if ( dim == 0) {
      std::vector< std::vector<size_t> > out(1, std::vector<size_t>(1));
      out[0] = {0};
      return out;
    }
    if ( dim == 1) {
      std::vector< std::vector<size_t> > out(Nk, std::vector<size_t>(2));
      for (size_t o(0); o < Nk; o++) out[o] = {o, o};
      return out;
    }
    if ( dim == 2) {
      std::vector< std::vector<size_t> > out(Nk, std::vector<size_t>(3));
      for (size_t o(0); o < Nk; o++) out[o] = {o, o, o};
      return out;
    }
    if ( dim == 3) {
      std::vector< std::vector<size_t> > out(Nk*Nk, std::vector<size_t>(4));
      size_t it(0);
      for (size_t m(0); m < Nk; m++)
      for (size_t n(0); n < Nk; n++)
        out[n + m * Nk] = {idxMinusIdx(n,m), n, m, it++};
        //out[n + m * Nk] = {-1, n, m};
      return out;
    }
    if ( dim == 4) {
      std::vector< std::vector<size_t> > out(Nk*Nk*Nk, std::vector<size_t>(5));
      size_t it(0);
      for (size_t m(0); m < Nk; m++)
      for (size_t n(0); n < Nk; n++)
      for (size_t o(0); o < Nk; o++)
        out[o + n * Nk + m * Nk * Nk] = {m, n, o, getForthIdx(o,n,m), it++};
      return out;
    }
    assert(0);
    return {};
  }
  void kMesh::print(){
    std::cout << "Print the k-mesh indices:" << std::endl;
    int cnt(0);
    for (size_t i(0); i < mesh[0]; i++)
    for (size_t j(0); j < mesh[1]; j++)
    for (size_t k(0); k < mesh[2]; k++)
      std::cout << cnt++ << ": " << i << " " << j << " " << k << std::endl;
    std::cout << "=====================\n";
  }


