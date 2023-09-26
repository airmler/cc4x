#pragma once
#include <ibst.h>

template <typename F>
class tensor;

class World {
  ibst::World *wrld;

public:
  World() {
    wrld = new ibst::World();
    np = wrld->np;
    rank = wrld->rank;
  }
  MPI_Comm comm() { return wrld->comm; }
  int np, rank;

  template <typename F>
  friend class tensor;
};

template <typename F>
class tensor {
public:
  ibst::Tensor<F> *machine_tensor;

  tensor(int64_t const _order,
         std::vector<int64_t> const _lens,
         std::vector<std::vector<int64_t>> const _nzc,
         World *_wrld,
         std::string const _name) {

    nzc = _nzc;
    lens = _lens;
    machine_tensor = new ibst::Tensor<F>(_order, lens, nzc, _wrld->wrld, _name);
  }

  ~tensor() { delete machine_tensor; };

  void contract(F alpha,
                const tensor &A,
                std::string const idxA,
                const tensor &B,
                std::string const idxB,
                F beta,
                std::string const idxC,
                bool verbose = false) {

    machine_tensor->contract(alpha,
                             *A.machine_tensor,
                             idxA,
                             *B.machine_tensor,
                             idxB,
                             beta,
                             idxC,
                             verbose);
  }

  void sum(F alpha,
           const tensor &A,
           std::string const idxA,
           F beta,
           std::string const idxB,
           std::vector<int64_t> remap,
           const std::function<F(const F)> &fseq,
           bool verbose = false) {
    machine_tensor
        ->sum(alpha, *A.machine_tensor, idxA, beta, idxB, remap, fseq, verbose);
  }

  void sum(F alpha,
           const tensor &A,
           std::string const idxA,
           F beta,
           std::string const idxB,
           std::vector<int64_t> remap,
           bool verbose = false) {
    machine_tensor
        ->sum(alpha, *A.machine_tensor, idxA, beta, idxB, remap, verbose);
  }

  void sum(F alpha,
           const tensor &A,
           std::string const idxA,
           F beta,
           std::string const idxB,
           const std::function<F(const F)> &fseq,
           bool verbose = false) {
    machine_tensor
        ->sum(alpha, *A.machine_tensor, idxA, beta, idxB, fseq, verbose);
  }

  void sum(F alpha,
           const tensor &A,
           std::string const idxA,
           F beta,
           std::string const idxB,
           bool verbose = false) {
    machine_tensor->sum(alpha, *A.machine_tensor, idxA, beta, idxB, verbose);
  }

  void read_dense_from_file(MPI_File &file) {
    machine_tensor->read_dense_from_file(file);
  }

  void write(int64_t npair,
             std::vector<int64_t> const global_idx,
             std::vector<F> const data,
             int64_t const block = -1) {
    machine_tensor->write(npair, global_idx, data, block);
  }

  void read(int64_t npair,
            std::vector<int64_t> const global_idx,
            std::vector<F> data,
            int64_t const block = -1) {
    machine_tensor->read(npair, global_idx, data, block);
  }

  void read(F &data) { machine_tensor->read(data); }

  void slice(std::vector<int64_t> const offsets,
             std::vector<int64_t> const ends,
             F beta,
             tensor const &A,
             std::vector<int64_t> const offsets_A,
             std::vector<int64_t> const ends_A,
             F alpha) {
    machine_tensor->slice(offsets,
                          ends,
                          beta,
                          *A.machine_tensor,
                          offsets_A,
                          ends_A,
                          alpha);
  }

  void relabelBlocks(std::vector<std::vector<int64_t>> nzcOut,
                     bool verbose = false) {
    machine_tensor->relabelBlocks(nzcOut, verbose);
  }

  std::vector<int64_t> lens;
  std::vector<std::vector<int64_t>> nzc;
};
