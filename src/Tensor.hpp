#pragma once

#ifdef CTF_BACKEND
#  include "tensor-backends/ctf.hpp"
#ifdef CTFBS_BACKEND
#  include "tensor-backends/ctfbs.hpp"
#else
#  error "backend not done"
#endif
