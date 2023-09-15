#pragma once

#include <ctf.hpp>

  template <typename F>
  class tensor;

  // TODO: maybe replace everywhere int or int64_t by this
  using TensorIdx = int64_t;

  class World {
    CTF::World * wrld;
    public:
    World() {
      wrld = new CTF::World();
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
      CTF::Tensor<F> *machine_tensor;

      tensor(int                                       order,
             std::vector<int64_t> const                len,
             std::vector< std::vector<int64_t> > const nonZero,
             World                                   * _wrld,
             std::string const                         name) {
        std::vector<int> sym(order);
        machine_tensor =
          new CTF::Tensor<F>(order, len.data(), sym.data(), *(_wrld->wrld), name.c_str());
        lens = len;
        int64_t i(0);
        for (auto nz: nonZero){
          std::vector<int64_t> r(nz);
          r.push_back(i++);
          nzc.push_back(r);
        }
      }

      ~tensor() { delete machine_tensor; }

      void
      contract(F          alpha,
               tensor     &A,
               std::string const idx_A,
               tensor     &B,
               std::string const idx_B,
               F          beta,
               std::string const idx_C,
               bool       verbose = false) {

        machine_tensor->contract(alpha,
                                 *A.machine_tensor,
                                 idx_A.c_str(),
                                 *B.machine_tensor,
                                 idx_B.c_str(),
                                 beta,
                                 idx_C.c_str()
                                );

      }

      void sum(F                                   alpha,
               tensor                            & A,
               std::string const                   idx_A,
               F                                   beta,
               std::string const                   idx_B,
               std::vector<int64_t>                remap,
               std::function<F(const F)>         & fseq,
               bool                                verbose=false) {

        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A.c_str(),
                            beta,
                            idx_B.c_str(),
                            CTF::Univar_Function<F>(fseq)
                           );

      }

      void sum(F                                   alpha,
               tensor                            & A,
               std::string const                   idx_A,
               F                                   beta,
               std::string const                   idx_B,
               std::vector<int64_t>                remap,
               bool                                verbose=false) {

        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A.c_str(),
                            beta,
                            idx_B.c_str()
                           );

      }

      void sum(F                           alpha,
               tensor                    & A,
               std::string const           idx_A,
               F                           beta,
               std::string const           idx_B,
               std::function<F(const F)> & fseq,
               bool                        verbose=false) {

        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A.c_str(),
                            beta,
                            idx_B.c_str(),
                            CTF::Univar_Function<F>(fseq)
                           );

      }

      void sum(F                  alpha,
               tensor           & A,
               std::string const  idx_A,
               F                  beta,
               std::string const  idx_B,
               bool               verbose=false) {

        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A.c_str(),
                            beta,
                            idx_B.c_str()
                           );

      }

      void read_dense_from_file(MPI_File &file) {
        machine_tensor->read_dense_from_file(file);
      }

      void write(int64_t                      npair,
                 std::vector<int64_t> const & global_idx,
                 std::vector<F> const       & data,
                 int64_t const                block = 0) {
        machine_tensor->write(npair, global_idx.data(), data.data());
      }

      //TODO: provide the second read function
      void read(F *data) {
        machine_tensor->read_all(data);
      }

      void slice(std::vector<int64_t> const offsets,
                 std::vector<int64_t> const ends,
                 F                          beta,
                 tensor               const &A,
                 std::vector<int64_t> const offsets_A,
                 std::vector<int64_t> const ends_A,
                 F                          alpha ) {

        machine_tensor->slice( offsets.data()
                             , ends.data()
                             , beta
                             , *A.machine_tensor
                             , offsets_A.data()
                             , ends_A.data()
                             , alpha);
      }

      void relabelBlocks( std::vector< std::vector<int64_t> > nonZeroOut
                        , bool verbose = false
                        ){
        /* .... */
      }


    std::vector<int64_t> lens;
    std::vector< std::vector<int64_t> > nzc;

  };


