module OpenFST

using CxxWrap, Test, SparseArrays

@wrapmodule(joinpath(@__DIR__,"../deps/lib","libfst_jl"))

function __init__()
    @initcxx
end

const runtests = true

export SymbolTable

export StdVectorFst, addarc!, setfinal!, setstart!, isfinal, setinputsymbols!, setoutputsymbols!, inputsymbols, outputsymbols, unpack, properties, numstates, addstate!, start, final, StateIterator, ArcIterator

export isexpanded, ismutable, iserror, isacceptor, isdeterministic, hasepsilons, islabelsorted, isweighted, iscyclic, isinitialcyclic, istopsorted, isaccessible, iscoaccessible, isstring, hasweightedcycles

export arcsort!, closure, closure!, compose, concat, concat!, connect, connect!, determinize, difference, disambiguate, epsnormalize, equivalent, isomorphic, minimize, minimize_rt, minimize!, minimize_rt!, project, push, rmepsilon, rmepsilon!, synchronize, topsort!, verify, encode, encode!, decode, decode!

export canonize!, omegawords!, state, states, cantorencode, cantorencode!, cantordecode, cantordecode!

# SymbolTable
addsymbol!(s::SymbolTable, name::Symbol, id::Integer) = addsymbol!(s, String(name), id)
addsymbol!(s::SymbolTable, name::Symbol) = addsymbol!(s, String(name))
Base.setindex!(s::SymbolTable, id::Integer, name::Union{String,Symbol}) = addsymbol!(s, name, id)
find(s::SymbolTable, name::Symbol) = find(s, String(name))
find(s::SymbolTable, id::Integer) = Symbol(_find(s, id))
Base.getindex(s::SymbolTable, id::Integer) = find(s, id)
Base.getindex(s::SymbolTable, name::Union{String,Symbol}) = (id = find(s, name); id < 0 ? throw("Unknown symbol $name") : id)
Base.length(s::SymbolTable) = Int(numsymbols(s))
function Base.iterate(s::SymbolTable, state = 0)
    if state < length(s)
        key = getnthkey(s,state)
        (find(s,key) => key, state+1)
    else
        nothing
    end
end
Base.show(io::IO, s::SymbolTable) = print(io, "SymbolTable($(length(s)) symbol(s)…)")
function SymbolTable(generator; no_epsilon = false)
    s = SymbolTable()
    no_epsilon || addsymbol!(s, "ϵ", 0)
    for g=generator
        addsymbol!(s, g...)
    end
    s
end
SymbolTable(generator...; no_epsilon = false) = SymbolTable(generator, no_epsilon = no_epsilon)

runtests && @testset "Symbol tables" begin
    sym = SymbolTable(:a => 3, :b => 2)
    @test Dict(sym) == Dict(:a => 3, :b => 2, :ϵ => 0)
    @test length(sym) == 3
    @test sym[0] == :ϵ
    @test sym[3] == :a
    @test sym[:a] == 3
    @test length(SymbolTable(Dict(:a => 3, :b => 2))) == 3
    @test length(SymbolTable((:a => 3,), no_epsilon = true)) == 1
end

# StdVectorFst
Base.length(f::StdVectorFst) = Int(numstates(f))
addarc!(f::StdVectorFst, from::Integer, to::Integer, ilabel::Integer, olabel::Integer = ilabel, weight::Float32 = 0.f0) = _addarc!(f, from, to, ilabel, olabel, weight)
addarc!(f::StdVectorFst, from::Integer, to::Integer, ilabel::Symbol, olabel::Symbol = ilabel, weight::Float32 = 0.f0) = _addarc!(f, from, to, inputsymbols(f)[ilabel], outputsymbols(f)[olabel], weight)
setfinal!(f::StdVectorFst, to::Integer, weight::Bool = true) = setfinal!(f, to, weight ? 0.f0 : Inf32)
isfinal(f::StdVectorFst, from::Integer) = !isinf(final(f,from))
setinputsymbols!(f::StdVectorFst, s::SymbolTable) = _setinputsymbols!(f, ConstCxxPtr(s))
inputsymbols(f::StdVectorFst) = (s = _inputsymbols(f); isnull(s) ? nothing : s[])
setoutputsymbols!(f::StdVectorFst, s::SymbolTable) = _setoutputsymbols!(f, ConstCxxPtr(s))
outputsymbols(f::StdVectorFst) = (s = _outputsymbols(f); isnull(s) ? nothing : s[])

function isdeterministic(f::StdVectorFst, side::Symbol = :INPUT)
    if side == :INPUT
        return isideterministic(f)
    elseif side == :OUTPUT
        return isodeterministic(f)
    end
    error("side should be :INPUT or :OUTPUT")
end

function hasepsilons(f::StdVectorFst, side::Symbol)
    if side == :INPUT
        return hasiepsilons(f)
    elseif side == :OUTPUT
        return hasoepsilons(f)
    end
    error("side should be :INPUT or :OUTPUT")
end

function islabelsorted(f::StdVectorFst, side::Symbol = :INPUT)
    if side == :INPUT
        return isilabelsorted(f)
    elseif side == :OUTPUT
        return isolabelsorted(f)
    end
    error("side should be :INPUT or :OUTPUT")
end

function Base.iterate(state::StateIterator, nonce = nothing)
    if Done(state)
        return nothing
    end
    v = Value(state)
    Next(state)
    (v, state)
end
function Base.iterate(state::ArcIterator, nonce = nothing)
    if Done(state)
        return nothing
    end
    v = Value(state)
    Next(state)
    (splat(v[]), state)
end
Base.IteratorSize(::ArcIterator) = Base.SizeUnknown()

function StdVectorFst(start::Integer, arcs::Union{Vector{Tuple{Int,Int,Symbol}},Vector{Tuple{Int,Int,Symbol,Symbol}}}, final, isymbols::SymbolTable, osymbols = isymbols)
    f = StdVectorFst()
    maxstate = max([a[1] for a in arcs]...,[a[2] for a in arcs]...,0)
    for i=0:maxstate
        addstate!(f)
    end
    setstart!(f, start)
    setinputsymbols!(f, isymbols)
    setoutputsymbols!(f, osymbols)
    for a=arcs
        addarc!(f, a...)
    end
    for s=final
        setfinal!(f, s)
    end
    f
end
StdVectorFst(start, arcs, final, isymbols::Dict{Symbol,Int}) = StdVectorFst(start, arcs, final, SymbolTable(isymbols))
StdVectorFst(start, arcs, final, isymbols::Dict{Symbol,Int}, osymbols::Dict{Symbol,Int}) = StdVectorFst(start, arcs, final, SymbolTable(isymbols), SymbolTable(osymbols))

# unpack fst into arguments to be fed to StdVectorFst (except for the symbol table at the end)
function unpack(f::StdVectorFst; acceptor::Union{Bool,Nothing} = nothing)
    init = start(f)
    finals = Int[]
    arcs = Tuple[]
    isyms = inputsymbols(f)
    osyms = outputsymbols(f)
    accept = isa(acceptor,Nothing) ? isacceptor(f) : acceptor
    
    for state=StateIterator(f)
        for (ilabel,olabel,weight,next)=ArcIterator(f,state)
            newarc = (state, next, isyms == nothing ? ilabel : isyms[ilabel])
            if !accept
                newarc = (newarc..., osyms == nothing ? olabel : osyms[olabel])
            end
            push!(arcs, newarc)
        end
        isfinal(f,state) && push!(finals, state)
    end
    if isyms != nothing && accept
        (init, convert(Vector{Tuple{Int,Int,Symbol}},arcs), finals)
    elseif isyms != nothing && osyms != nothing && !accept
        (init, convert(Vector{Tuple{Int,Int,Symbol,Symbol}},arcs), finals)
    else
        (init, arcs, finals)
    end
end

"""A hack to print an FST. Uses external commands."""
function Base.show(io::IO, f::StdVectorFst)
    if !(isdefined(Main, :IJulia) && Main.IJulia.inited)
        print(io, "StdVectorFst(", numstates(f)," state(s)…)")
        return
    end
        
    mktemp() do fname, fh
        _write(fname, f)
        close(fh)
        flags = ``
        acceptor = isacceptor(f) ? `--acceptor` : ``
        img = read(pipeline(fname, `fstdraw $flags $acceptor --portrait`, `dot -Tsvg`), String)
        display("image/svg+xml", img)
    end
end

runtests && @testset "StdVectorFst" begin
    f = StdVectorFst()
    @test unpack(f) == (-1, [], [])
    init = addstate!(f)
    setstart!(f, init)
    @test unpack(f) == (init, [], [])
    setfinal!(f, init)
    addarc!(f, 0, 0, 0, 0)
    setinputsymbols!(f, SymbolTable(Dict(:a => 1)))
    setoutputsymbols!(f, SymbolTable(Dict(:z => 1)))
    addarc!(f, 0, 0, :a, :z)
    @test unpack(f, acceptor=false) == (init, [(0, 0, :ϵ, :ϵ), (0, 0, :a, :z)], [0])
end

# Methods for StdVectorFst
function arcsort!(f::StdVectorFst, side = :StdILabelCompare)
    if side == :StdILabelCompare
        _arcsort_i!(f)
    elseif side == :StdOLabelCompare
        _arcsort_o!(f)
    else
        error("side must be :StdILabelCompare or :StdOLabelCompare")
    end
    testerror(f)
end

function closure!(f::StdVectorFst, closure_type::Symbol = :*)
    if closure_type == :*
        _closure_s!(f)
    elseif closure_type == :+
        _closure_p!(f)
    else
        error("closure_type should be :* or :+")
    end
    testerror(f)
end
closure(f::StdVectorFst, closure_type...) = closure!(copy(f), closure_type...)

compose(fsts::StdVectorFst...) = reduce(_compose,fsts)
Base.:∘(fsts::StdVectorFst...) = compose(fsts...)

function concat(fsts::StdVectorFst...)
    f = copy(fsts[end])
    for i=length(fsts)-1:-1:1
        concat!(fsts[i], f)
    end
    testerror(f)
end
Base.:*(fsts::StdVectorFst...) = concat(fsts...)
Base.:^(f::StdVectorFst, n::Integer) = (n == 0 ? StdVectorFst(0,[],[0]) : n == 1 ? f : n < 0 ? error("^(::StdVectorFst,n::Int) needs n≥0") : isodd(n) ? (f*f)^(n÷2)*f : (f*f)^(n÷2))

connect(f::StdVectorFst) = connect!(copy(f))

Base.:-(fst1::StdVectorFst, fst2::StdVectorFst) = difference(fst1, fst2)

function epsnormalize(f::StdVectorFst, side::Symbol = :INPUT)
    if side == :INPUT
        return _epsnormalize_i(f)
    elseif side == :OUTPUT
        return _epsnormalize_o(f)
    end
    error("side should be :INPUT or :OUTPUT")
end

Base.:≈(f::StdVectorFst, g::StdVectorFst) = equivalent(f, g)

Base.intersect(fsts::StdVectorFst...) = reduce(_intersect, fsts)

function Base.intersect!(f::StdVectorFst, fsts::StdVectorFst...)
    for g in fsts
        f = _intersect(f, g)
        testerror(f)
    end
    f
end

Base.issubset(f::StdVectorFst, g::StdVectorFst) = equivalent(f,f∩g)
              
Base.inv(f::StdVectorFst) = inv!(copy(f))

≅(f::StdVectorFst, g::StdVectorFst) = isomorphic(f, g)

minimize!(f::StdVectorFst) = minimize!(f, false)
minimize_rt!(f::StdVectorFst) = minimize_rt!(f, false)
minimize(f::StdVectorFst, allow_nondet...) = minimize!(copy(f), allow_nondet...)
minimize_rt(f::StdVectorFst, allow_nondet...) = minimize_rt!(copy(f), allow_nondet...)

function project(f::StdVectorFst, side::Symbol = :INPUT)
    if side == :INPUT
        return _project_i(f)
    elseif side == :OUTPUT
        return _project_o(f)
    end
    error("side should be :INPUT or :OUTPUT")
end

function pushlabels(f::StdVectorFst, side::Symbol = :REWEIGHT_TO_INITIAL)
    if side == :REWEIGHT_TO_INITIAL
        return _push_i(f)
    elseif side == :REWEIGHT_TO_FINAL
        return _push_f(f)
    end
    error("side should be :REWEIGHT_TO_INITIAL or :REWEIGHT_TO_FINAL")
end

rmepsilon(f::StdVectorFst) = rmepsilon!(copy(f))

topsort(f::StdVectorFst) = topsort!(copy(f))

function Base.union!(f::StdVectorFst, fsts::StdVectorFst...)
    for g in fsts
        _union!(f, g)
    end
    f
end
Base.union(fsts::StdVectorFst...) = union!(copy(fsts[1]), fsts[2:end]...)

encode(f::StdVectorFst) = (g = copy(f); e = encode!(g); (g,e))

decode(f::StdVectorFst, encoder) = decode!(copy(f), encoder)

function cantorencode!(f::StdVectorFst; symbols::Bool = false)
    transducer2acceptor!(f)
    symbols && transducer2acceptor_st!(f)
    f
end
cantorencode(f::StdVectorFst; args...) = cantorencode!(copy(f); args...)

function cantordecode!(f::StdVectorFst; symbols::Bool = false)
    acceptor2transducer!(f)
    symbols && acceptor2transducer_st!(f)
    f
end
cantordecode(f::StdVectorFst; args...) = cantordecode!(copy(f); args...)

function Base.hash(f::StdVectorFst, inithashkey::UInt...)
    hashkey = hash(numstates(f),inithashkey...)
    for state=StateIterator(f), arc=ArcIterator(f,state)
        hashkey = hash(arc,hashkey)
    end
    hashkey
end

# Missing:
#@ Prune
#@ RandEquivalent
#@ RandGen
#@ Relabel
#@ Replace
#@ Reweight
#@ ShortestDistance
#@ ShortestPath
#@ StateMap

"""Some missing commands"""
function cleanedup(fst::StdVectorFst)
    result = copy(fst)
    rmepsilon!(result)
    result = determinize(result)
    ns = numstates(result)
    while true
        rmepsilon!(result)
        minimize!(result)
        ns2 = numstates(result)
        ns == ns2 && break
        ns = ns2
    end
    arcsort!(result)
    result
end

function prefixes(fst::StdVectorFst)
    result = copy(fst)
    new = addstate!(result)
    setstart!(result, new)
    for state=StateIterator(result)
        addarc!(result, new, state, :ϵ, :ϵ)
        setfinal!(result, state)
    end
    result
end

function states(f::StdVectorFst, level = 1)
    if level == 0
        return reshape([f],1,1)
    end
    din = length(inputsymbols(f))-1
    dout = length(outputsymbols(f))-1
    topstates = [state(f,i,o) for o=1:dout, i=1:din] # transposed!
    substates = [states(g,level-1) for g=topstates]
    hvcat(dout,substates...)
end
    
runtests && @testset "Immediate operations" begin
    triple = StdVectorFst(0, [(0, 0, :𝟎, :𝟎), (0, 1, :𝟏, :𝟏), (1, 2, :𝟏, :𝟎), (1, 0, :𝟎, :𝟏), (2, 1, :𝟎, :𝟎), (2, 2, :𝟏, :𝟏)], [0, 1, 2], Dict(:𝟎 => 1, :𝟏 => 2))

    arcsort!(triple)

    @test triple |> unpack == (0, [(0, 0, :𝟎, :𝟎), (0, 1, :𝟏, :𝟏), (1, 0, :𝟎, :𝟏), (1, 2, :𝟏, :𝟎), (2, 1, :𝟎, :𝟎), (2, 2, :𝟏, :𝟏)], [0, 1, 2])

    sheep = StdVectorFst(0, [(0, 1, :b), (1, 2, :a), (2, 3, :a), (3, 3, :a)], [3], Dict(:a => 1, :b => 2))

    @test closure!(copy(sheep)) |> unpack == (closure(sheep,:*) |> unpack) == (4, [(0, 1, :b), (1, 2, :a), (2, 3, :a), (3, 3, :a), (3, 0, :ϵ), (4, 0, :ϵ)], [3, 4])
    @test closure(sheep,:+) |> unpack == (0, [(0, 1, :b), (1, 2, :a), (2, 3, :a), (3, 3, :a), (3, 0, :ϵ)], [3])
    @test sheep |> unpack == (0, [(0, 1, :b), (1, 2, :a), (2, 3, :a), (3, 3, :a)], [3])

    @test numstates(compose(triple,triple)) == 9

    concat(sheep,sheep,sheep)
    sheep^2 |> unpack == (4, [(0, 1, :b), (1, 2, :a), (2, 3, :a), (3, 3, :a), (4, 5, :b), (5, 6, :a), (6, 7, :a), (7, 0, :ϵ), (7, 7, :a)], [3])

    sheeep = rmepsilon!(sheep∪sheep)
    @test disambiguate(sheeep) ≈ sheep
    @test numstates(determinize(sheep∪sheep)) == 5

    sheep_s = determinize(rmepsilon!(closure(sheep,:*)))
    sheep_p = determinize(rmepsilon!(closure(sheep,:+)))

    @test difference(sheep_s,sheep_p) |> unpack == (0, [], [0])

    @test numstates(epsnormalize(sheep∪sheep)) == 7

    @test sheep_s∩sheep_p == sheep_p

    @test inv(triple)∘triple |> isacceptor

    @test project(triple) |> unpack == (0, [(0, 0, :𝟎), (0, 1, :𝟏), (1, 0, :𝟎), (1, 2, :𝟏), (2, 1, :𝟎), (2, 2, :𝟏)], [0, 1, 2])
    @test project(triple,:OUTPUT) |> unpack == (0, [(0, 0, :𝟎), (0, 1, :𝟏), (1, 0, :𝟏), (1, 2, :𝟎), (2, 1, :𝟎), (2, 2, :𝟏)], [0, 1, 2])

    p = arcsort!(inv(triple)∘triple)
    @test numstates(minimize(p)) == 1

    @test equivalent(p,minimize(p))

    # push
    # synchronize
    # topsort
    # verify

    encodee, encoder = encode(triple)
    @test decode(encodee, encoder) == triple

    # canonize    
end

"""Growth series of languages."""
function growth(f::StdVectorFst, maxdegree::Int = 10)
    init, arcs, accept = unpack(rmepsilon(f))
    maxstate = max([a[1] for a in arcs]...,[a[2] for a in arcs]...,0)+1
    mat = sparse([a[1]+1 for a in arcs], [a[2]+1 for a in arcs], [BigInt(1) for a in arcs], maxstate, maxstate)
    vecfinal = sparsevec(accept.+1, [BigInt(1) for i in accept],maxstate)

    v = sparsevec([init+1],[BigInt(1)], maxstate)'
    l = BigInt[]
    for i=1:maxdegree
        push!(l, v*vecfinal)
        v *= mat
    end
    l
end

function balls(gens::Vector,n::Int)
    ball = [Set(gens)]
    for i=2:n
        push!(ball,typeof(ball[1])())
        for x=ball[i÷2], y=ball[(i+1)÷2]
            push!(ball[i],x*y)
        end
    end
    ball
end

runtests && @testset "Growth" begin
    @test growth(StdVectorFst(0,[(0,1,:c),(1,2,:c)],[2],Dict(:c=>1))) == [0,0,1,0,0,0,0,0,0,0]
end

include("buchi.jl")

end
