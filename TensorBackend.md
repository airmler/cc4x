# Description of the tensor API

In the following the type `F` is the tensor template type
and `Tensor` is a short form for `cc4x::Tensor`.

## Constructors

A tensor should be a class providing
the following constructors

```cpp
Tensor(int                                      order,
       std::vector<int64_t> const               len,
       std::vector< std::vector<int> > const    nonZero);
```

- `order` number of dimensions of tensor
- `len` edge lengths of tensor
- `nonZero` TODO


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

### The sum function

```cpp
void
sum(F          alpha,
   Tensor     &A,
   char const *idx_A,
   F          beta,
   char const *idx_B);

```

When the object `B` is called, the operation

```
B[idx_B] = beta*B[idx_B] + alpha*A[idx_A]
```

is done



### The slice function (optional)

TODO
