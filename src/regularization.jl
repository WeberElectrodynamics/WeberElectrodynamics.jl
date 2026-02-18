@inline function _compute_pair_distances!(
    rb::RegularizationBuffers,
    q::Vector{Float64},
)
    dims = rb.dims
    @inbounds for k = 1:rb.n_pairs
        i = rb.pair_i[k]
        j = rb.pair_j[k]
        qi_start = (i - 1) * dims
        qj_start = (j - 1) * dims
        r2 = 0.0
        for d = 1:dims
            dq = q[qi_start+d] - q[qj_start+d]
            r2 += dq * dq
        end
        rb.pair_distance[k] = sqrt(r2)
    end
    return nothing
end

@inline function _find_min_pair(rb::RegularizationBuffers)::Tuple{Int,Float64}
    min_idx = 1
    min_r = rb.pair_distance[1]
    @inbounds for k = 2:rb.n_pairs
        r = rb.pair_distance[k]
        if r < min_r
            min_r = r
            min_idx = k
        end
    end
    return min_idx, min_r
end

@inline function _build_adjacency!(rb::RegularizationBuffers, threshold::Float64)
    n = rb.n_particles
    @inbounds begin
        fill!(rb.adjacency, false)
        for i = 1:n
            rb.adjacency[i, i] = true
        end
        for k = 1:rb.n_pairs
            if rb.pair_distance[k] <= threshold
                i = rb.pair_i[k]
                j = rb.pair_j[k]
                rb.adjacency[i, j] = true
                rb.adjacency[j, i] = true
            end
        end
    end
    return nothing
end

@inline function _build_component_from_anchor!(
    rb::RegularizationBuffers,
    anchor_i::Int,
    anchor_j::Int,
)
    n = rb.n_particles
    adjacency = rb.adjacency
    visited = rb.visited
    queue = rb.queue
    nodes = rb.component_nodes

    @inbounds begin
        fill!(visited, false)

        head = 1
        tail = 0

        tail += 1
        queue[tail] = anchor_i
        visited[anchor_i] = true

        if !visited[anchor_j]
            tail += 1
            queue[tail] = anchor_j
            visited[anchor_j] = true
        end

        count = 0
        while head <= tail
            node = queue[head]
            head += 1

            count += 1
            nodes[count] = node

            for nbr = 1:n
                if adjacency[node, nbr] && !visited[nbr]
                    visited[nbr] = true
                    tail += 1
                    queue[tail] = nbr
                end
            end
        end

        fill!(rb.component_mask, false)
        for idx = 1:count
            rb.component_mask[nodes[idx]] = true
        end

        rb.active_count = count
        for idx = 1:count
            rb.active_nodes[idx] = nodes[idx]
        end
    end

    return rb.active_count
end

@inline function _component_min_distance(rb::RegularizationBuffers)::Float64
    if rb.active_count <= 1
        return Inf
    end
    mask = rb.component_mask
    min_r = Inf
    @inbounds for k = 1:rb.n_pairs
        i = rb.pair_i[k]
        j = rb.pair_j[k]
        if mask[i] && mask[j]
            r = rb.pair_distance[k]
            if r < min_r
                min_r = r
            end
        end
    end
    return min_r
end

@inline function _component_max_distance(rb::RegularizationBuffers)::Float64
    if rb.active_count <= 1
        return 0.0
    end
    mask = rb.component_mask
    max_r = 0.0
    @inbounds for k = 1:rb.n_pairs
        i = rb.pair_i[k]
        j = rb.pair_j[k]
        if mask[i] && mask[j]
            r = rb.pair_distance[k]
            if r > max_r
                max_r = r
            end
        end
    end
    return max_r
end

@inline function _build_chain_order!(
    rb::RegularizationBuffers,
    q::Vector{Float64},
)
    n_active = rb.active_count
    if n_active == 0
        return 0
    end

    dims = rb.dims
    nodes = rb.active_nodes
    chain = rb.chain_order
    used = rb.chain_used

    @inbounds begin
        fill!(used, false)

        min_node = nodes[1]
        for idx = 2:n_active
            node = nodes[idx]
            if node < min_node
                min_node = node
            end
        end

        chain[1] = min_node
        used[min_node] = true

        for pos = 2:n_active
            prev = chain[pos-1]
            prev_start = (prev - 1) * dims
            best_node = 0
            best_r2 = Inf

            for idx = 1:n_active
                cand = nodes[idx]
                if !used[cand]
                    cand_start = (cand - 1) * dims
                    r2 = 0.0
                    for d = 1:dims
                        dq = q[prev_start+d] - q[cand_start+d]
                        r2 += dq * dq
                    end
                    if r2 < best_r2 || (r2 == best_r2 && cand < best_node)
                        best_r2 = r2
                        best_node = cand
                    end
                end
            end

            chain[pos] = best_node
            used[best_node] = true
        end
    end

    return n_active
end

@inline function _component_omega(rb::RegularizationBuffers)::Float64
    mask = rb.component_mask
    omega = 0.0
    @inbounds for k = 1:rb.n_pairs
        i = rb.pair_i[k]
        j = rb.pair_j[k]
        if mask[i] && mask[j]
            r = rb.pair_distance[k]
            omega += 1.0 / max(r, eps(Float64))
        end
    end
    return omega
end

@inline function _extract_pair_relative_state!(
    rb::RegularizationBuffers,
    q::Vector{Float64},
    p::Vector{Float64},
    masses::Vector{Float64},
    i::Int,
    j::Int,
)
    dims = rb.dims
    qi_start = (i - 1) * dims
    qj_start = (j - 1) * dims

    mi = masses[i]
    mj = masses[j]
    mu = mi * mj / (mi + mj)

    rel_q = rb.rel_q
    rel_p = rb.rel_p

    @inbounds for d = 1:dims
        rel_q[d] = q[qi_start+d] - q[qj_start+d]
        vi_d = p[qi_start+d] / mi
        vj_d = p[qj_start+d] / mj
        rel_p[d] = mu * (vi_d - vj_d)
    end

    return mu
end

@inline function _lc_lift!(
    rb::RegularizationBuffers,
)
    x = rb.rel_q[1]
    y = rb.rel_q[2]
    px = rb.rel_p[1]
    py = rb.rel_p[2]

    r = sqrt(x * x + y * y)
    if r <= eps(Float64)
        rb.lc_u[1] = 0.0
        rb.lc_u[2] = 0.0
        rb.lc_U[1] = 0.0
        rb.lc_U[2] = 0.0
        return 0.0
    end

    u1 = sqrt(max((r + x) * 0.5, 0.0))
    u2 = 0.0
    if u1 > 100 * eps(Float64)
        u2 = y / (2 * u1)
    else
        u2 = sqrt(max((r - x) * 0.5, 0.0))
        if y < 0
            u2 = -u2
        end
        if abs(u2) > 100 * eps(Float64)
            u1 = y / (2 * u2)
        else
            u1 = 0.0
        end
    end

    rb.lc_u[1] = u1
    rb.lc_u[2] = u2

    rb.lc_U[1] = 2 * (u1 * px + u2 * py)
    rb.lc_U[2] = 2 * (-u2 * px + u1 * py)

    return r
end

@inline function _lc_project!(
    q_rel::Vector{Float64},
    p_rel::Vector{Float64},
    u::Vector{Float64},
    U::Vector{Float64},
)
    u1 = u[1]
    u2 = u[2]

    q_rel[1] = u1 * u1 - u2 * u2
    q_rel[2] = 2 * u1 * u2

    r = u1 * u1 + u2 * u2
    if r <= eps(Float64)
        p_rel[1] = 0.0
        p_rel[2] = 0.0
        return 0.0
    end

    p_rel[1] = (u1 * U[1] - u2 * U[2]) / (2 * r)
    p_rel[2] = (u2 * U[1] + u1 * U[2]) / (2 * r)

    return r
end

@inline function _ks_jacobian!(J::Matrix{Float64}, u::Vector{Float64})
    u1, u2, u3, u4 = u

    @inbounds begin
        J[1, 1] = 2 * u1
        J[1, 2] = -2 * u2
        J[1, 3] = -2 * u3
        J[1, 4] = 2 * u4

        J[2, 1] = 2 * u2
        J[2, 2] = 2 * u1
        J[2, 3] = -2 * u4
        J[2, 4] = -2 * u3

        J[3, 1] = 2 * u3
        J[3, 2] = 2 * u4
        J[3, 3] = 2 * u1
        J[3, 4] = 2 * u2
    end

    return nothing
end

@inline function _ks_constraint(u::Vector{Float64}, U::Vector{Float64})::Float64
    return u[4] * U[1] - u[3] * U[2] + u[2] * U[3] - u[1] * U[4]
end

@inline function _ks_project_constraint!(
    U::Vector{Float64},
    u::Vector{Float64},
    n::Vector{Float64},
)
    @inbounds begin
        n[1] = u[4]
        n[2] = -u[3]
        n[3] = u[2]
        n[4] = -u[1]
    end

    n2 = n[1] * n[1] + n[2] * n[2] + n[3] * n[3] + n[4] * n[4]
    if n2 <= eps(Float64)
        return abs(_ks_constraint(u, U))
    end

    psi = _ks_constraint(u, U)
    alpha = psi / n2
    @inbounds for k = 1:4
        U[k] -= alpha * n[k]
    end

    return abs(_ks_constraint(u, U))
end

@inline function _ks_lift!(
    rb::RegularizationBuffers,
)
    x1 = rb.rel_q[1]
    x2 = rb.rel_q[2]
    x3 = rb.rel_q[3]
    p1 = rb.rel_p[1]
    p2 = rb.rel_p[2]
    p3 = rb.rel_p[3]

    r = sqrt(x1 * x1 + x2 * x2 + x3 * x3)
    u = rb.ks_u
    U = rb.ks_U

    if r <= eps(Float64)
        fill!(u, 0.0)
        fill!(U, 0.0)
        return 0.0
    end

    if r + x1 > 100 * eps(Float64)
        den = sqrt(2 * (r + x1))
        u[1] = 0.5 * den
        u[2] = x2 / den
        u[3] = x3 / den
        u[4] = 0.0
    else
        den = sqrt(max(2 * (r - x1), 0.0))
        u[1] = 0.0
        u[2] = 0.5 * den
        if den > 100 * eps(Float64)
            u[3] = x2 / den
            u[4] = -x3 / den
        else
            u[3] = 0.0
            u[4] = 0.0
        end
    end

    J = rb.ks_J
    _ks_jacobian!(J, u)

    @inbounds begin
        U[1] = J[1, 1] * p1 + J[2, 1] * p2 + J[3, 1] * p3
        U[2] = J[1, 2] * p1 + J[2, 2] * p2 + J[3, 2] * p3
        U[3] = J[1, 3] * p1 + J[2, 3] * p2 + J[3, 3] * p3
        U[4] = J[1, 4] * p1 + J[2, 4] * p2 + J[3, 4] * p3
    end

    return r
end

@inline function _detect_regularization_component!(
    rb::RegularizationBuffers,
    q::Vector{Float64},
    chain_enabled::Bool,
)
    _compute_pair_distances!(rb, q)

    if rb.n_pairs == 0
        rb.is_active = false
        rb.active_mode = REG_MODE_NONE
        rb.active_count = 0
        return false, REG_MODE_NONE, Inf
    end

    min_pair_idx, min_r = _find_min_pair(rb)

    if rb.is_active
        _build_adjacency!(rb, rb.r_off)
        _build_component_from_anchor!(rb, rb.active_anchor_i, rb.active_anchor_j)
        comp_max_r = _component_max_distance(rb)
        if comp_max_r > rb.r_off || rb.active_count < 2
            rb.is_active = false
            rb.active_mode = REG_MODE_NONE
            rb.active_count = 0
        end
    end

    if !rb.is_active
        if min_r <= rb.r_on
            rb.is_active = true
            rb.active_anchor_i = rb.pair_i[min_pair_idx]
            rb.active_anchor_j = rb.pair_j[min_pair_idx]

            _build_adjacency!(rb, rb.r_on)
            _build_component_from_anchor!(rb, rb.active_anchor_i, rb.active_anchor_j)
        end
    end

    if rb.is_active
        if rb.active_count == 2
            rb.active_mode = REG_MODE_PAIR
        elseif rb.active_count > 2 && chain_enabled
            rb.active_mode = REG_MODE_CHAIN
            _build_chain_order!(rb, q)
        else
            rb.active_mode = REG_MODE_NONE
        end
    end

    return rb.is_active, rb.active_mode, min_r
end
