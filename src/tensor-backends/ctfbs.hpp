#pragma once

#include <ctf.hpp>

namespace cc4x {


  template <typename F>
  class Tensor {
    public:
      CTF::bsTensor<F> *machine_tensor;

      Tensor(int                                      order,
             std::vector<int64_t> const               len,
             std::vector< std::vector<int> > const    nonZero) {

        // TODO
        // build lens
        // build sym
        // build world
        machine_tensor = new CTF::bsTensor<F>(oder,....);
      }

    void
    contract(F          alpha,
             Tensor     &A,
             char const *idx_A,
             Tensor     &B,
             char const *idx_B,
             F          beta,
             char const *idx_C) {

      machine_tensor->contract(alpha,
                               *A.machine_tensor,
                               idx_A,
                               *B.machine_tensor,
                               idx_B,
                               beta,
                               idx_C);

    }

  };


}


