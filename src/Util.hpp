#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <ctf.hpp>
#include <yaml-cpp/yaml.h>

using ivec = std::vector<int>;
using i64vec = std::vector<int64_t>;
using Complex = std::complex<double>;

#define THROW(msg) std::cout << msg << std::endl; throw msg;

#define LOG() if (!cc4x::dw->rank) std::cout

