/*
template <typename F>
struct InputTensor {
  bool initialize = false;
  CTF::Tensor<F> *data;
};

template <typename F>
InputTensor<F>* make_def_tensor() {
  auto a = new InputTensor<F>();
  a->initialize = true;
  return a;
}
*/
#define THROW(msg) std::cout << msg << std::endl; throw msg;

  //auto epsi = make_def_tensor<Complex>();

