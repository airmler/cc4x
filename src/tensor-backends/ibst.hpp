#pragma once
#include <ibst.h>
#include <ctf.hpp>

  class World {
   public:
    ibst::iWorld * ibst_world;
    static ibst::iWorld wrld(bool dryRun);

  };

  template <typename F>
  class tensor {
    public:
      ibst::Tensor<F> *machine_tensor;

      tensor(size_t const                       _order,
             std::vector<size_t>  const         _lens,
             std::vector< std::vector<size_t> > const _nzc,
             World       &                      _wrld,
             std::string const                  _name) {

        nzc = _nzc;
        lens = _lens;
        machine_tensor =
//          new ibst::Tensor<F>(_order, lens, nzc, _name);
          new ibst::Tensor<F>(_order, lens, nzc, _wrld.wrld, _name);
      }

      ~tensor() { delete machine_tensor; };

      void
      contract(F                 alpha,
               const tensor     &A,
               std::string const idxA,
               const tensor     &B,
               std::string const idxB,
               F                 beta,
               std::string const idxC,
               bool              verbose=false) {

        machine_tensor->contract(alpha,
                                 *A.machine_tensor,
                                 idxA,
                                 *B.machine_tensor,
                                 idxB,
                                 beta,
                                 idxC,
                                 verbose);

      }

      void sum(F                          alpha,
               const tensor              &A,
               std::string const          idxA,
               F                          beta,
               std::string const          idxB,
               std::vector<size_t>        remap,
               const std::function<F(const F)> &fseq,
               bool                       verbose=false){
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idxA,
                            beta,
                            idxB,
                            remap,
                            fseq,
                            verbose);
      }

      void sum(F                   alpha,
               const tensor          &   A,
               std::string const   idxA,
               F                   beta,
               std::string const   idxB,
               std::vector<size_t> remap,
               bool                verbose=false){
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idxA,
                            beta,
                            idxB,
                            remap,
                            verbose);
      }

      void sum(F                            alpha,
               const tensor               & A,
               std::string const            idxA,
               F                            beta,
               std::string const            idxB,
               const std::function<F(const F)>   &fseq,
               bool                         verbose=false){
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idxA,
                            beta,
                            idxB,
                            fseq,
                            verbose);
      }

      void sum(F                  alpha,
               const tensor     & A,
               std::string const  idxA,
               F                  beta,
               std::string const  idxB,
               bool               verbose=false){
        machine_tensor->sum(alpha,
                            *A.machine_tensor,
                            idxA,
                            beta,
                            idxB,
                            verbose);
      }

      void read_dense_from_file(MPI_File &file){
        machine_tensor->read_dense_from_file(file);
      }

      void write(size_t                    npair,
                 std::vector<size_t> const global_idx,
                 std::vector<F>      const data){
        machine_tensor->write(npair, global_idx, data);
      }

      void write(size_t                    npair,
                 std::vector<size_t> const global_idx,
                 std::vector<F>      const data,
                 int64_t             const block){
        machine_tensor->write(npair, global_idx, data, block);
      }

      void read(size_t                    npair,
                std::vector<size_t> const global_idx,
                std::vector<F>            data){
        machine_tensor->read(npair, global_idx, data);
      }

      void read(F &data) {
        machine_tensor->read(data);
      }

      void slice(std::vector<size_t> const offsets,
                 std::vector<size_t> const ends,
                 F                         beta,
                 tensor const             &A,
                 std::vector<size_t> const offsets_A,
                 std::vector<size_t> const  ends_A,
                 F                         alpha ){
        machine_tensor->slice(
          offsets, ends, beta, *A.machine_tensor, offsets_A, ends_A, alpha
        );
      }

      void relabelBlocks( std::vector< std::vector<size_t> > nzcOut
                        , bool verbose = false
                        ){
        machine_tensor->relabelBlocks(nzcOut, verbose);
      }

    std::vector<size_t> lens;
    std::vector< std::vector<size_t> > nzc;

  };


