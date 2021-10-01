primitive type F2 <: Signed 8 end

F2(x::Int8) = reinterpret(F2, x)
F2(x::Number) = F2(Int8(x))
Core.Int8(x::F2) = reinterpret(Int8, x)
Core.Int32(x::F2) = Int32(Int8(x))
Core.Int64(x::F2) = Int64(Int8(x))
Core.Bool(x::F2) = ~iszero(x)
Base.show(io::IO, x::F2) = print(io, Int8(x))
#Base.:+(a::F2, b::F2) = F2(mod(Int8(a) + Int8(b), Int8(2)))
Base.:+(a::F2, b::F2) = F2(Int8(a) ⊻ Int8(b))
Base.:-(a::F2, b::F2) = a + b
Base.:*(a::F2, b::F2) = F2(Int8(a) & Int8(b))

function standard_basis(T, dim::Integer)
    return Matrix{T}(I, (dim, dim))
end

function standard_basis(example::Array)
    T = eltype(example)
    dim, _ = size(example)
    return standard_basis(T, dim)
end

function swaprow!(x, i, j)
    for k in axes(x, 2)
        x[i, k], x[j, k] = x[j, k], x[i, k]
    end
    return x
end

function swapcol!(x, i, j)
    for k in axes(x, 1)
        x[k, i], x[k, j] = x[k, j], x[k, i]
    end
    return x
end

nonzero = (!) ∘ iszero

function smith_normal_form(M::Matrix{F2})
    @inbounds begin
    M = copy(M)
    nrows, ncols = size(M)
    U = standard_basis(F2, nrows)
    V = standard_basis(F2, ncols)
    for c in 1:min(nrows, ncols)
        found = false
        k = c
        l = c
        while !found
            if nonzero(M[k, l])
                swaprow!(M, c, k)
                swaprow!(U, c, k)
                swapcol!(M, c, l)
                swapcol!(V, c, l)
                found = true
            elseif l == ncols && k == nrows
                # no nonzeros in row or column
                found = true
            elseif l == ncols
                l = c
                k += 1
            else
                l += 1
            end
        end
        for i in findall(nonzero, view(M, c+1:nrows, c)) .+ c
            M[i, :] += M[c, :]
            U[i, :] += U[c, :]
        end
        for j in findall(nonzero, view(M, c, c+1:ncols)) .+ c
            M[:, j] += M[:, c]
            V[:, j] += V[:, c]
        end
    end
    end
    return U, M, V
end

function smith_normal_form(M::Union{BitMatrix, Matrix{Bool}})
    @inbounds begin
    M = copy(M)
    nrows, ncols = size(M)
    U = standard_basis(Bool, nrows)
    V = standard_basis(Bool, ncols)
    for c in 1:min(nrows, ncols)
        found = false
        k = c
        l = c
        while !found
            if nonzero(M[k, l])
                swaprow!(M, c, k)
                swaprow!(U, c, k)
                swapcol!(M, c, l)
                swapcol!(V, c, l)
                found = true
            elseif l == ncols && k == nrows
                # no nonzeros in row or column
                found = true
            elseif l == ncols
                l = c
                k += 1
            else
                l += 1
            end
        end
        for i in findall(identity, view(M, c+1:nrows, c)) .+ c
            for j in c:ncols
                M[i, j] = M[i, j] ⊻ M[c, j]
            end
            U[i, :] = U[i, :] .⊻ U[c, :]
        end
        for j in findall(identity, view(M, c, c+1:ncols)) .+ c
            for i in c:nrows
                M[i, j] = M[i, j] ⊻ M[i, c]
            end
            V[:, j] = V[:, j] .⊻ V[:, c]
        end
    end
    end
    return U, M, V
end

function rowreduction(M::Matrix{F2})
    @inbounds begin
    M = copy(M)
    nrows, ncols = size(M)
    U = standard_basis(F2, nrows)
    pivotcols = []
    pivotcol = 1 :: Int64
    for r in 1:nrows
        if pivotcol > ncols # check
            return U, M, pivotcols
        end
        i = r :: Int64
        while iszero(M[i, pivotcol])
            i += 1
            if i > nrows
                i = r
                pivotcol += 1
                if pivotcol > ncols #check
                    return U, M, pivotcols
                end
            end
        end
        swaprow!(M, i, r)
        swaprow!(U, i, r)
        # findall(nonzero, view(M, :, pivotcol))
        for j in 1:nrows
            if nonzero(M[j, pivotcol]) && j != r
                for k in pivotcol:ncols
                    M[j, k] = M[j, k] + M[r, k] #check
                end
                U[j, :] += U[r, :]
            end
        end
        push!(pivotcols, pivotcol)
        pivotcol += 1
    end
    end
    return U, M, pivotcols
end

function rowreduction(M::Union{BitMatrix, Matrix{Bool}})
    @inbounds begin
    M = copy(M)
    nrows, ncols = size(M)
    U = standard_basis(Bool, nrows)
    pivotcols = []
    pivotcol = 1 :: Int64
    for r in 1:nrows
        if pivotcol > ncols # check
            return U, M, pivotcols
        end
        i = r :: Int64
        while !M[i, pivotcol]
            i += 1
            if i > nrows
                i = r
                pivotcol += 1
                if pivotcol > ncols #check
                    return U, M, pivotcols
                end
            end
        end
        swaprow!(M, i, r)
        swaprow!(U, i, r)
        for j in findall(view(M, :, pivotcol))
            if j != r
                for k in pivotcol:ncols
                    M[j, k] = M[j, k] ⊻ M[r, k] #check
                end
                U[j, :] = U[j, :] .⊻ U[r, :]
            end
        end
        push!(pivotcols, pivotcol)
        pivotcol += 1
    end
    end
    return U, M, pivotcols
end

function issquare(M)
    nrows, ncols = size(M)
    square = nrows == ncols
    return square
end

function isfullrank(M)
    _, echelonform, _ = rowreduction(M)
    fullrank = isfullrank(M, echelonform)
    return fullrank
end

function isfullrank(M, echelonform)
    fullrank = all(any(map(nonzero, echelonform), dims=2))
    return fullrank
end

function inverse(M)
    inverse, echelonform, _ = rowreduction(M)
    square = issquare(M)
    fullrank = isfullrank(M, echelonform)
    if ~square || ~fullrank
        throw(ErrorException("Matrix is not invertible"))
    else
        return inverse
    end
end

function kernel_image(M)
    U, M, V = smith_normal_form(M)
    kernel = V[:, [i[2] for i in findall(!, any(nonzero.(M), dims=1))]]
    image = inverse(U)[:, [i[1] for i in findall(any(nonzero.(M), dims=2))]]
    return kernel, image
end

function extend_basis_relative(U, V)#::FilteredBasis
    dim, nelements = size(U)
    all_vectors = hcat(U, V)
    transform, echelonform, pivotcols = rowreduction(all_vectors)
    extended_basis = all_vectors[:, pivotcols]
    return [U, extended_basis[:, (nelements + 1):end]]
end

function extend_basis(U)#::FilteredBasis
    return extend_basis_relative(U, standard_basis(U))
end

function extend_basis(Us...)#::FilteredBasis
    Us = [Us...]
    ambient_basis = standard_basis(first(Us))
    push!(Us, ambient_basis)
    descending = reverse(Us)
    subspace = pop!(descending)
    filtered_bases = [subspace]
    while ~isempty(descending)
        superspace = pop!(descending)
        _, extension = extend_basis_relative(subspace, superspace)
        push!(filtered_bases, extension)
        subspace = hcat(subspace, extension)
    end
    return filtered_bases
end
