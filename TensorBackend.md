# Description of the tensor API

In the following the type `F` is the tensor template type
of `tensor`. For the current version of the code the tensor framework needs to support the following functionalities:
- a constructor
- a destructor (currently missing)
- a ```contract``` function for ```C["ij"] = A["ik"] B["kj"]```
- four ```sum``` functions; ```B["ij"] = A["ji"]``` without block sparsity there are only two different versions of sum.
- ```read_dense_from_file``` to read a tensor from a binary file
- ```write``` to write data from C++ array to tensor
- ```read_all``` to read data from the tensor to a C++ pointer (Note: the whole (distributed) tensor is written on every rank).
- ```slice``` cut out a slice of a tensor


## Constructors

A tensor should be a class providing
the following constructors

```cpp
Tensor(int                                      order,
       std::vector<int64_t> const               len,
       std::vector< std::vector<int> > const    nonZero
       CTF::World *                             wlrd
       char const *                             name);
```

- `order` number of dimensions of tensor
- `len` edge lengths of tensor
- `nonZero` non-zero conditions for possible block sparsity
- `wrld` currently still the CTF world
- `name` name of the tensor

For a simple matrix the constructor would read ```M(2, {40,50}, {0, 0}, wrld_, "M");``` 
Have a look at ```tensor-backends/ctf.hpp:21``` to see how to handle non-zero conditions if the underlying backend does not support block sparsity.

## Member functions


### The contract function

```cpp
void
contract(F          alpha,
         Tensor     &A,
         char const *idx_A,
         Tensor     &B,
         char const *idx_B,
         F          beta,
         char const *idx_C);

```

The object `C` when it's called with
```cpp
C.contract(...TODO)
```
does the following operation

```
 C[idx_C] = beta*C[idx_C] + alpha * A[idx_A] * B[idx_B]
```

### The sum functions

```cpp
void
sum(F          alpha,
   Tensor     &A,
   char const *idx_A,
   F          beta,
   char const *idx_B,
   std::vector< std::vector<int> > nonZeroA,
   std::vector< std::vector<int> > nonZeroB,
   std::function<F(const F)>       &fseq,
   bool verbose=false);

```
```cpp
void
sum(F          alpha,
   Tensor     &A,
   char const *idx_A,
   F          beta,
   char const *idx_B,
   std::vector< std::vector<int> > nonZeroA,
   std::vector< std::vector<int> > nonZeroB,
   bool verbose=false);

```
```cpp
void
sum(F          alpha,
   Tensor     &A,
   char const *idx_A,
   F          beta,
   char const *idx_B,
   std::function<F(const F)>    &fseq,
   bool verbose=false);

```
```cpp
void
sum(F          alpha,
   Tensor     &A,
   char const *idx_A,
   F          beta,
   char const *idx_B,
   bool verbose=false);

```


When the object `B` is called, the operation

```
B[idx_B] = beta*B[idx_B] + fseq(alpha*A[idx_A])
```

is done. If the framework does not support block sparsity the nonZeroConditions are redundant.

### Read from file

```cpp
void 
read_dense_from_file(MPI_File &file);
```
Object reads data from MPI File.

### Write

```cpp
void
write(int64_t         npair,
           int64_t const * global_idx,
           F const       * data);
```

Object reads ```npair``` elements from the adress ```data```. The data is written according to the provided index map ```global_idx```. The global index ```g=i+j*m+k*m*N+l*m*n*p``` is associated with the tensor element (i,j,k,l) with tensor edge lenghts {m,n,p,q}.

### Read all

```cpp
void
read_all(F *data);
```
Object writes the whole tensor data to pointer ```data```. Note: this is done on every rank (i.e. not memory scalable).

### The slice function 

```cpp
void
slice(std::vector<int64_t> const offsets,
      std::vector<int64_t> const ends,
      F                          beta,
      tensor               const &A,
      std::vector<int64_t> const offsets_A,
      std::vector<int64_t> const ends_A,
      F                          alpha );
```

The object ```C``` when calling 
```cpp
C.slice( ... )
```
adds to a slice of this tensor 
```
C[offsets,ends)=beta*B[offsets,ends) + alpha*A[offsets_A,ends_A) 
```

