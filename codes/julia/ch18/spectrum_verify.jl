# spectrum_verify.jl
#
# Chapter 18: Linear Spectral Eigenproblems --- shared utility.
#
# Implements the drift-with-N diagnostic of Boyd (2000, Ch. 7.5), which
# compares two spectra computed at different truncations and decides, mode
# by mode, which eigenvalues are believable.
#
# Given spectra lambda_N1 and lambda_N2 (with N2 > N1 typically), the
# routine computes three diagnostics per low-resolution mode j:
#
#   1. intermodal separation                sigma_j             (Boyd 7.19)
#   2. scaled ordinal drift  delta_ord_j = |lam1[j] - lam2[j]|    / sigma_j
#   3. scaled nearest drift  delta_nst_j = min_k |lam1[j] - lam2[k]| / sigma_j
#
# A mode is flagged `trusted` if max(delta_ord_j, delta_nst_j) < tol and
# both delta_ord_j and delta_nst_j are finite. The chapter recommends
# plotting the reciprocal 1/delta_j so trusted modes float to the top of
# a semilog plot.
#
# This utility is used from Etude 18.5 onward and is the chapter's
# central reusable artefact.
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

module SpectrumVerify

export intermodal_separation, ordinal_drift, nearest_drift, verify_spectrum,
       SpectrumReport

"""
    intermodal_separation(lam)

Compute sigma_j per Boyd Eq. 7.19 for a real-sorted spectrum.

    sigma_1     = |lam_1 - lam_2|
    sigma_j     = (|lam_j - lam_{j-1}| + |lam_{j+1} - lam_j|) / 2   for 1 < j < N
    sigma_N     = |lam_N - lam_{N-1}|

The separation guards against division by zero when an eigenvalue is
accidentally close to zero (|lam_j| would be a poor scale) and is more
appropriate than |lam_j| for problems where neighbouring eigenvalues
matter more than absolute magnitude.
"""
function intermodal_separation(lam::AbstractVector{<:Real})
    n = length(lam)
    sigma = zeros(n)
    if n == 1
        sigma[1] = abs(lam[1])  # degenerate fallback
        return sigma
    end
    sigma[1] = abs(lam[1] - lam[2])
    for j in 2:n-1
        sigma[j] = 0.5 * (abs(lam[j] - lam[j-1]) + abs(lam[j+1] - lam[j]))
    end
    sigma[n] = abs(lam[n] - lam[n-1])
    # Protect against degenerate-eigenvalue zero separations: fall back to |lam|
    for j in 1:n
        if sigma[j] < 1e-14
            sigma[j] = max(abs(lam[j]), 1e-14)
        end
    end
    return sigma
end

"""
    ordinal_drift(lam1, lam2)

Compute scaled ordinal drift delta_j,ordinal = |lam1[j] - lam2[j]| / sigma_j
per Boyd Eq. 7.20, where sigma_j is the intermodal separation of lam1.

The returned vector has length length(lam1); if lam2 is shorter than
lam1 only the first length(lam2) entries are meaningful and the rest are
Inf.
"""
function ordinal_drift(lam1::AbstractVector{<:Real}, lam2::AbstractVector{<:Real})
    n1 = length(lam1)
    n2 = length(lam2)
    sigma = intermodal_separation(lam1)
    delta = fill(Inf, n1)
    for j in 1:min(n1, n2)
        delta[j] = abs(lam1[j] - lam2[j]) / sigma[j]
    end
    return delta
end

"""
    nearest_drift(lam1, lam2)

Compute scaled nearest drift delta_j,nearest = min_k |lam1[j] - lam2[k]| / sigma_j
per Boyd Eq. 7.21.

Nearest matching is essential when mode ordering is not invariant
between resolutions (e.g. Laplace tidal equations where Rossby and
gravity waves interlace). For well-ordered problems the nearest and
ordinal drifts agree.
"""
function nearest_drift(lam1::AbstractVector{<:Real}, lam2::AbstractVector{<:Real})
    n1 = length(lam1)
    sigma = intermodal_separation(lam1)
    delta = fill(Inf, n1)
    for j in 1:n1
        # minimum |lam1[j] - lam2[k]| over all k
        mn = Inf
        for k in 1:length(lam2)
            d = abs(lam1[j] - lam2[k])
            if d < mn
                mn = d
            end
        end
        delta[j] = mn / sigma[j]
    end
    return delta
end

"""
    SpectrumReport

Container for the output of `verify_spectrum`. Fields:
    - lam1            : low-resolution spectrum (input, possibly sorted)
    - lam2            : high-resolution spectrum (input, possibly sorted)
    - sigma           : intermodal separation of lam1
    - delta_ordinal   : scaled ordinal drift (Boyd 7.20)
    - delta_nearest   : scaled nearest drift (Boyd 7.21)
    - trusted         : Boolean vector; trusted[j] = true iff mode j believable
    - tol             : tolerance used
    - n_trusted       : sum(trusted)
"""
struct SpectrumReport
    lam1::Vector{Float64}
    lam2::Vector{Float64}
    sigma::Vector{Float64}
    delta_ordinal::Vector{Float64}
    delta_nearest::Vector{Float64}
    trusted::Vector{Bool}
    tol::Float64
    n_trusted::Int
end

"""
    verify_spectrum(lam1, lam2; tol=1e-3, sort_input=true)

Drift-with-N diagnostic. Given two spectra at different N, return a
SpectrumReport.

A mode j is flagged trusted when BOTH the ordinal and nearest drifts
fall below `tol`. Using only the ordinal drift misses problems with
mode reordering; using only nearest misses problems where two bad
eigenvalues happen to lie close by. Demanding both is conservative and
matches Boyd's pedagogical recommendation.

When lam1 or lam2 is complex, pass the sorted real parts (or |lam|) --
the utility requires real, sorted input and does not silently convert.
"""
function verify_spectrum(lam1::AbstractVector{<:Real},
                         lam2::AbstractVector{<:Real};
                         tol::Real = 1e-3,
                         sort_input::Bool = true)
    l1 = sort_input ? sort(collect(float.(lam1))) : collect(float.(lam1))
    l2 = sort_input ? sort(collect(float.(lam2))) : collect(float.(lam2))
    sigma = intermodal_separation(l1)
    d_ord = ordinal_drift(l1, l2)
    d_nst = nearest_drift(l1, l2)
    trusted = (d_ord .< tol) .& (d_nst .< tol) .& isfinite.(d_ord) .& isfinite.(d_nst)
    return SpectrumReport(l1, l2, sigma, d_ord, d_nst, trusted,
                          float(tol), count(trusted))
end

end # module
