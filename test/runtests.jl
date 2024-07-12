using OpenFST
using Test

@testset "Symbol tables" begin
    sym = SymbolTable(:a => 3, :b => 2)
    @test Dict(sym) == Dict(:a => 3, :b => 2, :ϵ => 0)
    @test length(sym) == 3
    @test sym[0] == :ϵ
    @test sym[3] == :a
    @test sym[:a] == 3
    @test length(SymbolTable(Dict(:a => 3, :b => 2))) == 3
    @test length(SymbolTable((:a => 3,), no_epsilon = true)) == 1
end

@testset "StdVectorFst" begin
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

@testset "Immediate operations" begin
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

@testset "Growth" begin
    @test growth(StdVectorFst(0,[(0,1,:c),(1,2,:c)],[2],Dict(:c=>1))) == [0,0,1,0,0,0,0,0,0,0]
end
