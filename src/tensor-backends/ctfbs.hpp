#pragma once

#include <ctf.hpp>



  template <typename F>
  class tensor {
    public:
      CTF::bsTensor<F> *machine_tensor;

      tensor(int                                      order,
             std::vector<int64_t> const               len,
             std::vector< std::vector<int> > const    nonZero,
             CTF::World *                             wrld,
             char const *                             name) {

        machine_tensor = new CTF::bsTensor<F>(order, len, nonZero, wrld, name);
        lens = len;
        int i(0);
        for (auto nz: nonZero){
          std::vector<int> r(nz);
          r.push_back(i++);
          nonZeroCondition.push_back(r);
        }
      }
      ~tensor() { delete machine_tensor; };
      void
      contract(F          alpha,
               tensor     &A,
               char const *idx_A,
               tensor     &B,
               char const *idx_B,
               F          beta,
               char const *idx_C,
               bool       verbose=false) {

        machine_tensor->contract(alpha,
                                 *A.machine_tensor,
                                 idx_A,
                                 *B.machine_tensor,
                                 idx_B,
                                 beta,
                                 idx_C,
                                 verbose);

      }

      void sum(F                               alpha,
               tensor                          & A,
               char const *                    idx_A,
               F                               beta,
               char const *                    idx_B,
               std::vector< std::vector<int> > nonZeroA,
               std::vector< std::vector<int> > nonZeroB,
               std::function<F(const F)>       &fseq,
               bool                            verbose=false){
        CTF::Univar_Function<F> f(fseq);
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A,
                            beta,
                            idx_B,
                            nonZeroA,
                            nonZeroB,
                            f,
                            verbose);
      }

      void sum(F                               alpha,
               tensor                          & A,
               char const *                    idx_A,
               F                               beta,
               char const *                    idx_B,
               std::vector< std::vector<int> > nonZeroA,
               std::vector< std::vector<int> > nonZeroB,
               bool                            verbose=false){
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A,
                            beta,
                            idx_B,
                            nonZeroA,
                            nonZeroB,
                            verbose);
      }

      void sum(F                               alpha,
               tensor                          & A,
               char const *                    idx_A,
               F                               beta,
               char const *                    idx_B,
               std::function<F(const F)>       &fseq,
               bool                            verbose=false){

        CTF::Univar_Function<F> f(fseq);
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A,
                            beta,
                            idx_B,
                            f,
                            verbose);
      }

      void sum(F                               alpha,
               tensor                          & A,
               char const *                    idx_A,
               F                               beta,
               char const *                    idx_B,
               bool                            verbose=false){
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idx_A,
                            beta,
                            idx_B,
                            verbose);
      }

      void read_dense_from_file(MPI_File &file){
        machine_tensor->read_dense_from_file(file);
      }

      void write(int64_t         npair,
                 int64_t const * global_idx,
                 F const       * data){
        machine_tensor->write(npair, global_idx, data);
      }

      void read_all(F *data){
        machine_tensor->read_all(data);
      }

      void slice(std::vector<int64_t> const offsets,
                 std::vector<int64_t> const ends,
                 F                          beta,
                 tensor               const &A,
                 std::vector<int64_t> const offsets_A,
                 std::vector<int64_t> const ends_A,
                 F                          alpha ){
        machine_tensor->slice(offsets, ends, beta, *A.machine_tensor, offsets_A, ends_A, alpha);
      }

    std::vector<int64_t> lens;
    std::vector< std::vector<int> > nonZeroCondition;

  };


