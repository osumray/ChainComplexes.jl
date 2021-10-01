# ChainComplexes

## Overview
A package for working with [chain complexes](https://en.wikipedia.org/wiki/Chain_complex) over the finite field with two elements and their morphisms.

The homology groups of chain complexes can be calculated and whether a morphism of chain complexes is a [quasi-isomorphism](https://en.wikipedia.org/wiki/Quasi-isomorphism).

The type `ChainComplex{T}` contains a single field `differentials`, a vector of matrices of type `T`.
The ith element of the `differentials` field represents the ith differential in the chain complex.

The type `AlignedChainComplex{T}` can be constructed from a `ChainComplex{T}` type.
It finds a filtered basis for each underlying vector space of the chain complex such that
the differentials of this chain complex can be written in block form:

![](./images/blockform.svg)

and thus the ith [Betti number](https://ncatlab.org/nlab/show/Betti+number) can be found
by the number of zero columns of the ith differential minus the number of non-zero rows of the (i-1)th
differential. 


