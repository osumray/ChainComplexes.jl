abstract type AbstractChainComplex{T} end

struct ChainComplex{T} <: AbstractChainComplex{T}
    differentials :: Vector{Matrix{T}}
end

struct AlignedChainComplex{T} <: AbstractChainComplex{T}
    filteredbases::Dict{Int, Vector{Matrix{T}}}
    inversebases::Dict{Int, Matrix{T}}
    bettinumbers::Vector{Int}
end

Base.eltype(cc::ChainComplex) = (eltype ∘ first)(cc.differentials)
Base.eltype(acc::AlignedChainComplex) = (eltype ∘ first ∘ last ∘ first)(acc.filteredbases)
Base.length(cc::ChainComplex)::Integer = length(cc.differentials) + 1
Base.length(acc::AlignedChainComplex)::Integer = length(acc.filteredbases) + 1


function _bettinumbers(filteredbases, topdim)
    bettinumbers = [size(filteredbases[i][2])[2] for i in 0:topdim]
    while ~(isempty(bettinumbers) || nonzero(last(bettinumbers)))
        pop!(bettinumbers)
    end
    return bettinumbers
end

function AlignedChainComplex(cc::ChainComplex)
    T = eltype(cc)
    return AlignedChainComplex{T}(cc)
end

function AlignedChainComplex{T}(cc::ChainComplex{T}) where T
    topdim = length(cc.differentials)
    cycles = Dict{Int, Matrix{T}}()
    boundaries = Dict{Int, Matrix{T}}()
    filteredbases = Dict{Int, Vector{Matrix{T}}}()
    inversebases = Dict{Int, Matrix{T}}()
    for pairs in enumerate(cc.differentials)
        i, D = pairs
        ker, img = kernel_image(D)
        cycles[i] = ker
        boundaries[i - 1] = img
    end
    cycles[0] = (standard_basis ∘ first)(cc.differentials)
    boundaries[topdim] = zeros(T, ((last ∘ size ∘ last)(cc.differentials), 0))
    for i in 0:topdim
        filtration = extend_basis(boundaries[i], cycles[i])
        filteredbases[i] = filtration
        inversebases[i] = inverse(reduce(hcat, filtration))
    end
    bettinumbers = _bettinumbers(filteredbases, topdim)
    return AlignedChainComplex(filteredbases, inversebases, bettinumbers)
end

struct ChainComplexMorphism{T} <: Morphism{ChainComplex{T}}
    source::ChainComplex{T}
    target::ChainComplex{T}
    verticals::Vector{Matrix{T}}
end

struct AlignedChainComplexMorphism{T} <: Morphism{AlignedChainComplex{T}}
    source::AlignedChainComplex{T}
    target::AlignedChainComplex{T}
    verticals::Vector{Matrix{T}}
end

function AlignedChainComplexMorphism(m::ChainComplexMorphism)
    s = AlignedChainComplex(m.source)
    t = AlignedChainComplex(m.target)
    return AlignedChainComplexMorphism(s, t, m)
end

Base.eltype(cc::ChainComplexMorphism) = eltype(cc.source)
Base.eltype(acc::AlignedChainComplexMorphism) = eltype(acc.source)

function AlignedChainComplexMorphism{T}(m::ChainComplexMorphism{T},
                                     s::AlignedChainComplex{T},
                                     t::AlignedChainComplex{T}) where T
    alignedverticals = Vector{Matrix{T}}()
    for pair in enumerate(m.verticals)
        i, v = pair
        i -= 1
        newv = t.inversebases[i] * v * reduce(hcat, s.filteredbases[i])
        newv = eltype(m) == Bool ? Matrix{Bool}(newv .% 2) : newv
        push!(alignedverticals, newv)
    end
    # if eltype(m) == Bool
    #     for pair in enumerate(alignedverticals)
    #         i, v = pair
    #         alignedverticals[i] = Matrix{Bool}(v .% 2)
    #     end
    # end
    return AlignedChainComplexMorphism(s, t, alignedverticals)
end

function isverticaliso(m::AlignedChainComplexMorphism, degree::Int)::Bool
    if min(map(length, (m.source.bettinumbers, m.target.bettinumbers))...) < degree + 1
        return true
    else
        if m.source.bettinumbers[degree + 1] != m.target.bettinumbers[degree + 1]
            return false
        else
            s_bases = m.source.filteredbases[degree]
            t_bases = m.target.filteredbases[degree]
            _, s_dim_boundaries = size(s_bases[1])
            _, s_dim_cycles = size(s_bases[2])
            _, t_dim_boundaries = size(t_bases[1])
            _, t_dim_cycles = size(t_bases[2])
            v = m.verticals[degree + 1]
            source_homology = (s_dim_boundaries + 1):(s_dim_boundaries + s_dim_cycles)
            target_homology = (t_dim_boundaries + 1):(t_dim_boundaries + t_dim_cycles)
            homology_map = v[target_homology, source_homology]
            iso = issquare(homology_map) && isfullrank(homology_map)
            return iso
        end
    end
end

function isquasi_iso(m::AlignedChainComplexMorphism, upto)::Bool
    upto = min(upto, length(m.verticals) + 1)
    for degree in 0:(upto - 1)
        if !isverticaliso(m, degree)
            return false
        end
    end
    return true
end

function isquasi_iso(m::AlignedChainComplexMorphism)::Bool
    upto = length(m.verticals + 1)
    return isquasi_iso(m, upto)
end

const ≅ = isquasi_iso
