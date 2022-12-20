#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <ctf.hpp>
#include <Tensor.hpp>
#include <iomanip>
#include <yaml-cpp/yaml.h>
#include <chrono>
#include <map>

#ifndef UTIL_DEFINED
#define UTIL_DEFINED


using Complex = std::complex<double>;

#define THROW(msg) std::cout << msg << std::endl; throw msg;

#define LOG() if (!cc4x::dw->rank) std::cout << std::setprecision(8)

// taken from Atrip code of AGallo
struct Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Event = std::chrono::time_point<Clock>;
  std::chrono::duration<double> duration;
  Event _start;
  inline void start() noexcept { _start = Clock::now(); }
  inline void stop() noexcept { duration += Clock::now() - _start; }
  inline void clear() noexcept { duration *= 0; }
  inline double count() const noexcept { return duration.count(); }
};
using Timings = std::map<std::string, Timer>;


#endif
