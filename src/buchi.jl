# Buchi transducers, as (semi)group elements

export BuchiElement

struct BuchiElement
    t::StdVectorFst
    function BuchiElement(f::StdVectorFst)
        rmepsilon!(f)
        e = encode!(f)
        f = determinize(f)
        arcsort!(f)
        omegawords!(f)
        minimize!(f)
        decode!(f,e)
        OpenFST.canonize!(f)
        arcsort!(f)
        new(f)
    end
end

Base.:*(a::BuchiElement,b::BuchiElement) = BuchiElement(a.t∘b.t)
Base.inv(a::BuchiElement) = BuchiElement(inv(a.t))
Base.isone(a::BuchiElement) = isacceptor(a.t)
Base.one(a::BuchiElement) = BuchiElement(project(a.t))
Base.zero(a::BuchiElement) = BuchiElement(StdVectorFst(-1,Tuple{Int64, Int64, Symbol}[],0:-1,inputsymbols(a.t)))
Base.iszero(a::BuchiElement) = numstates(a.t)==0
Base.:^(a::BuchiElement, n::Int) = (n < 0 ? inv(a^(-n)) : n == 0 ? one(a) : n == 1 ? a : (b = a^(n÷2); isodd(n) ? b*b*a : b*b))
Base.:(==)(a::BuchiElement, b::BuchiElement) = a.t == b.t
Base.hash(a::BuchiElement) = hash(a.t)
Base.show(io::IO, a::BuchiElement) = show(io, a.t)

state(a::BuchiElement, i, o) = BuchiElement(state(a.t, i, o))
states(a::BuchiElement, args...) = states(a.t, args...) .|> BuchiElement
Base.copy(a::BuchiElement) = BuchiElement(copy(a.t))
cantorencode(a::BuchiElement; args...) = cantorencode(a.t; args...)
cantorencode!(a::BuchiElement; args...) = cantorencode!(a.t; args...)
cantordecode(a::BuchiElement; args...) = cantordecode(a.t; args...)
cantordecode!(a::BuchiElement; args...) = cantordecode!(a.t; args...)

Base.intersect(a::BuchiElement, args...) = BuchiElement(intersect(a.t, (a.t for a=args)...))
Base.union(a::BuchiElement, args...) = BuchiElement(union(a.t, (a.t for a=args)...))
Base.issubset(a::BuchiElement, b::BuchiElement) = issubset(a.t, b.t)
