
"""
Computes the worst case distribution with a bounded total deviation `ξ`
from the underlying probability distribution `p̄` for the random variable `z`.

Efficiently computes the solution of:
min_p   p^T * z
s.t.    || p - p̄ ||_1  ≤ ξ
        1^T p = 1
        p ≥ 0

Notes
-----
This implementation works in O(n log n) time because of the sort. Using
quickselect to choose the right quantile would work in O(n) time.

This function does not check whether the provided probability distribution sums
to 1.

Returns
-------
Optimal solution `p` and the objective value
"""
function worstcase_l1(z::Vector{Float64}, p̄::Vector{Float64}, ξ::Float64)
    (maximum(p̄) <= 1 + 1e-9 && minimum(p̄) >= -1e-9)  || "values must be between 0 and 1"
    ξ >= 0 || "ξ must be nonnegative"
    (length(z) > 0 && length(z) == length(p̄)) ||
            "z's values needs to be same length as p̄'s values"
    
    ξ = clamp(ξ, 0, 2)
    size = length(z)
    sorted_ind = sortperm(z)

    out = copy(p̄)       #duplicate it
    k = sorted_ind[1]   #index begins at 1

    ϵ = min(ξ / 2, 1 - p̄[k])
    out[k] += ϵ
    i = size -1

    while ϵ > 0 && i > 0
        k = sorted_ind[i]
        i -= 1
        difference = min(ϵ, out[k])
        out[k] -= difference
        ϵ -= difference
    end

    return out, out'*z
end
using LinearAlgebra
using DataStructures
using Printf
using Infinity

struct GradientsL1_w
    grads :: Vector{Float64}
    donors :: Vector{Int}
    receivers::Vector{Int}
    donor_greater::Vector{Bool}
    sorted::Vector{Float64}
end

function GradientsL1_w(z::Vector{Float64}, w::Vector{Float64})
    epsilon = 1e-8
    element_count = Int(length(z))

    @assert length(w) == element_count

    grads = Float64[]
    donors = Int[]
    receivers = Int[]
    donor_greater = Bool[]

    # Identifing possible receivers
    z_increasing = sortperm(z)
    possible_receivers = Int[]
    smallest_w = Inf

    for iz in z_increasing
        @assert w[iz] > epsilon
        if w[iz] < smallest_w
            push!(possible_receivers, iz)
            smallest_w = w[iz]
        end
    end

    # Computing grads for donor-receiver pairs
    for i = 1:element_count
        for j in possible_receivers
            if z[i] <= z[j]
                continue
            end
            # Case a: donor ≤ pbar value
            grad = (-z[i] + z[j]) / (w[i] + w[j])
            push!(grads, grad < -epsilon ? grad : 0)
            push!(donors, i)
            push!(receivers, j)
            push!(donor_greater, false)
        end
    end

    # Case b: donor > pbar value
    for i in possible_receivers
        for j in possible_receivers
            if z[i] <= z[j]
                continue
            end
            if abs(w[i] - w[j]) > epsilon && w[i] < w[j]
                grad = (-z[i] + z[j]) / (-w[i] + w[j])
                push!(grads, grad < -epsilon ? grad : 0)
                push!(donors, i)
                push!(receivers, j)
                push!(donor_greater, true)
            end
        end
    end

    sorted = sortperm(grads)

    return GradientsL1_w(grads, donors, receivers, donor_greater, sorted)
end

function steepest_solution(gradients::GradientsL1_w, index::Int)
    @assert index >= 1 && index <= Int(length(gradients.sorted))
    e = Int(gradients.sorted[index])
    return gradients.grads[e], gradients.donors[e], gradients.receivers[e], gradients.donor_greater[e]
end

function worstcase_l1_w(z::Vector{Float64}, pbar::Vector{Float64}, w::Vector{Float64}, xi::Float64)
    @assert maximum(pbar) <= 1 + 1e-9 && minimum(pbar) >= -1e-9 "values must be between 0 and 1"
    @assert xi >= 0 "xi must be nonnegative"
    @assert length(z) > 0 && length(z) == length(pbar) "z's values needs to be same length as pbar's values"
    @assert sum(pbar) == 1 "Values of pbar must sum to one"

    epsilon = 1e-10
    p = copy(pbar)
    xi_rest = xi
    grad_epsilon = 1e-5
    grad_que = [] #tuple

    gradients = GradientsL1_w(z, w)

    for k in 1:length(z)
        push!(grad_que, steepest_solution(gradients, Int(k)))

        while length(grad_que) > 1 && grad_que[1][1] < grad_que[end][1] - grad_epsilon
            popfirst!(grad_que)
        end

        for g in grad_que

            _, donor, receiver, donor_greater = g

            if receiver == 0 continue end

            if donor_greater && pbar[donor] <= pbar[donor] + epsilon continue end

            if !donor_greater && pbar[donor] > pbar[donor] + epsilon continue end

            if pbar[donor] < epsilon continue end

            weight_change = donor_greater ? (-w[donor] + w[receiver]) : (w[donor] + w[receiver])
            @assert weight_change > 0

            donor_step = min(xi_rest / weight_change, 
                             pbar[donor] > pbar[donor] + epsilon ? (pbar[donor] - pbar[donor]) : pbar[donor])
            pbar[donor] -= donor_step
            pbar[receiver] += donor_step
            xi_rest -= donor_step * weight_change

            if xi_rest < epsilon break end
        end
        if xi_rest < epsilon break end
    end

    objective = dot(pbar, z)
    return (pbar, objective)
end
