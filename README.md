# Scope of the project 

cc4x is a canoncial coupled cluster calculations code. It uses tensor
contraction libraries for the occurrent tensor operations. Currently, only CTF
(cyclops) is supported.  However, adding other contraction libraries is
straightforward (see TensorBackend.md). The offered functionalities are
restricted! Please check out [cc4s][https://github.com/cc4s/cc4s] for general
survey. 

# Dependencies
You will need [yaml-cpp][https://github.com/jbeder/yaml-cpp/] and
[cyclops][https://github.com/cyclops-community/ctf]
