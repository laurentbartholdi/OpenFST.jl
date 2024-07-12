# OpenFST

Implements finite state transducers via the C++ library [OpenFST](https://www.openfst.org). Its syntax follows closely that of the C++ library.

```julia
julia> using Pkg; Pkg.add("OpenFST"); Pkg.build("OpenFST")

julia> # the sheep language: ba, baa, baaa, ...

julia> # start state is 0, accept states are [3], alphabet is {a,b}

julia> sheep = StdVectorFst(0, [(0, 1, :b), (1, 2, :a), (2, 3, :a), (3, 3, :a)], [3], Dict(:a => 1, :b => 2))
StdVectorFst(4 state(s)…)

julia> OpenFST.growth(sheep) # number of words
OpenFST.growth(sheep)
10-element Vector{BigInt}:
 0
 0
 0
 1
 1
 1
 1
 1
 1
 1

julia> sheepsheep = concat(sheep,sheep) # concatenation: baba, baaba, babaa, ...
StdVectorFst(8 state(s)…)

julia> OpenFST.growth(sheepsheep)
10-element Vector{BigInt}:
 0
 0
 0
 0
 0
 0
 1
 2
 3
 4

julia> sheep∪sheepsheep |> determinize |> minimize
StdVectorFst(9 state(s)…)

julia> (initial, transitions, accepting) = unpack(ans)
(0, [(0, 1, :ϵ), (1, 2, :b), (2, 3, :a), (3, 4, :a), (4, 5, :ϵ), (4, 4, :a), (5, 6, :b), (6, 7, :a), (7, 8, :a), (8, 8, :a)], [4, 8])
```
