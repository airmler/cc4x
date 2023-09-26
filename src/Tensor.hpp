#pragma once

#ifdef CTF_BACKEND
#  include "tensor-backends/ctf.hpp"
#elif CTFBS_BACKEND
#  include "tensor-backends/ctfbs.hpp"
#elif IBST_BACKEND
#  include "tensor-backends/ibst.hpp"
#else
#  error "backend not done"
#endif
