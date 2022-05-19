#include <Kmesh.hpp>
#include <iostream>
#include <cassert>
#include <numeric>

  kMesh::kMesh(){};

  kMesh::kMesh(ivec _mesh){
    mesh.resize(_mesh.size());
    mesh = _mesh;
    Nk = std::accumulate(mesh.begin(), mesh.end(), 1, std::multiplies<int>());
  }

  int kMesh::backfold(const int i, const int N) { return (i+4*N)%N; }
  ivec kMesh::backfold(const ivec in){
    return { backfold(in[0], mesh[0])
           , backfold(in[1], mesh[1])
           , backfold(in[2], mesh[2])
           };
  }
  int kMesh::kToIdx(const ivec kPoint){
    return backfold(kPoint[0], mesh[0])
         + backfold(kPoint[1], mesh[1])*mesh[0]
         + backfold(kPoint[2], mesh[2])*mesh[0]*mesh[1];
  }
  ivec kMesh::idxToK(const int i){
    return {  i%mesh[0]
           , (i/mesh[0])%mesh[1]
           ,  i/mesh[0]/ mesh[1]
           };
  }
  int kMesh::idxMinusIdx(const int i, const int j){
    ivec ki(idxToK(i));
    ivec kj(idxToK(j));
    ivec d( { kj[0] - ki[0]
            , kj[1] - ki[1]
            , kj[2] - ki[2]
            } );
    return kToIdx(backfold(d));
  }
  int kMesh::getForthIdx(const int k, const int i, const int j){
    // kk + ko = ki + kj ==>  ko = ki + kj - kk
    ivec ki(idxToK(i));
    ivec kj(idxToK(j));
    ivec kk(idxToK(k));
    ivec d( { ki[0] + kj[0] - kk[0]
            , ki[1] + kj[1] - kk[1]
            , ki[2] + kj[2] - kk[1]
            } );
    return kToIdx(backfold(d));
  }
  std::vector<ivec> kMesh::getNZC(const int dim){
    if ( dim == 0) return {{}};
    if ( dim == 1) {
      std::vector<ivec> out(Nk, ivec(1));
      for (int o(0); o < Nk; o++) out[o] = {o};
      return out;
    }
    if ( dim == 2) {
      std::vector<ivec> out(Nk, ivec(2));
      for (int o(0); o < Nk; o++) out[o] = {o, o};
      return out;
    }
    if ( dim == 3) {
      std::vector<ivec> out(Nk*Nk, ivec(3));
      for (int m(0); m < Nk; m++)
      for (int n(0); n < Nk; n++)
        out[n + m * Nk] = {idxMinusIdx(n,m), n, m};
      return out;
    }
    if ( dim == 4) {
      std::vector<ivec> out(Nk*Nk*Nk, ivec(4));
      for (int m(0); m < Nk; m++)
      for (int n(0); n < Nk; n++)
      for (int o(0); o < Nk; o++)
        out[o + n * Nk + m * Nk * Nk] = {getForthIdx(o,n,m), o, n, m};
      return out;
    }
    assert(0);
    return {};
  }
  void kMesh::print(){
    std::cout << "Print the k-mesh indices:" << std::endl;
    int cnt(0);
    for (int i(0); i < mesh[0]; i++)
    for (int j(0); j < mesh[1]; j++)
    for (int k(0); k < mesh[2]; k++)
      std::cout << cnt++ << ": " << i << " " << j << " " << k << std::endl;
    std::cout << "=====================\n";
  }


