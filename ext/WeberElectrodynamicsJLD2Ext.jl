module WeberElectrodynamicsJLD2Ext

using JLD2
using WeberElectrodynamics

function _diagnostics_to_archive(d::WeberElectrodynamics.RegularizationDiagnostics)
    return (
        enabled = d.enabled,
        requested_backend = d.requested_backend,
        used_backend = d.used_backend,
        activation_count = d.activation_count,
        deactivation_count = d.deactivation_count,
        active_steps = d.active_steps,
        pair_steps = d.pair_steps,
        adaptive_pair_steps = d.adaptive_pair_steps,
        lifted_pair_steps = d.lifted_pair_steps,
        chain_steps = d.chain_steps,
        unregularized_steps = d.unregularized_steps,
        backend_fallback_steps = d.backend_fallback_steps,
        total_substeps = d.total_substeps,
        max_substeps_used = d.max_substeps_used,
        max_constraint_violation = d.max_constraint_violation,
        ks_constraint_projection_count = d.ks_constraint_projection_count,
        ks_constraint_iteration_count = d.ks_constraint_iteration_count,
        min_encounter_distance = d.min_encounter_distance,
        mode_history = copy(d.mode_history),
    )
end

function _diagnostics_from_archive(data, n_steps::Int)
    d = WeberElectrodynamics.RegularizationDiagnostics(
        data.enabled,
        n_steps,
        data.requested_backend,
        data.used_backend,
    )
    d.activation_count = data.activation_count
    d.deactivation_count = data.deactivation_count
    d.active_steps = data.active_steps
    d.pair_steps = data.pair_steps
    d.adaptive_pair_steps = data.adaptive_pair_steps
    d.lifted_pair_steps = data.lifted_pair_steps
    d.chain_steps = data.chain_steps
    d.unregularized_steps = data.unregularized_steps
    d.backend_fallback_steps = data.backend_fallback_steps
    d.total_substeps = data.total_substeps
    d.max_substeps_used = data.max_substeps_used
    d.max_constraint_violation = data.max_constraint_violation
    d.ks_constraint_projection_count = data.ks_constraint_projection_count
    d.ks_constraint_iteration_count = data.ks_constraint_iteration_count
    d.min_encounter_distance = data.min_encounter_distance
    n = min(length(data.mode_history), length(d.mode_history))
    if n > 0
        copyto!(d.mode_history, 1, data.mode_history, 1, n)
    end
    return d
end

function WeberElectrodynamics.save_solution(
    path::AbstractString,
    sol::WeberElectrodynamics.HamiltonianSolution;
    metadata = NamedTuple(),
)
    prob = sol.prob
    terms = WeberElectrodynamics.term_names(prob.system)
    terms == [:weber, :zollner] || throw(
        ArgumentError(
            "save_solution currently archives default Weber/Zollner systems only; got terms=$terms",
        ),
    )

    archive = (
        format_version = 1,
        problem = (
            n_particles = WeberElectrodynamics.n_particles(prob),
            dims = WeberElectrodynamics.dims(prob),
            tspan = prob.tspan,
            q_initial = copy(prob.q_initial),
            p_initial = copy(prob.p_initial),
            params = copy(WeberElectrodynamics.params(prob)),
            kappas = copy(WeberElectrodynamics.kappas(prob)),
            dt = prob.dt,
            convergence_tolerance = prob.convergence_tolerance,
            maximum_iterations = prob.maximum_iterations,
            term_names = terms,
        ),
        solution = (
            t = copy(sol.t),
            q = [copy(q) for q in sol.q],
            p = [copy(p) for p in sol.p],
            retcode = sol.retcode,
            regularization = _diagnostics_to_archive(sol.regularization),
        ),
        metadata = metadata,
    )

    JLD2.jldsave(path; archive = archive)
    return path
end

function WeberElectrodynamics.load_solution(path::AbstractString)
    data = JLD2.load(path)
    haskey(data, "archive") ||
        throw(ArgumentError("archive does not contain an `archive` entry"))
    archive = data["archive"]
    archive.format_version == 1 ||
        throw(ArgumentError("unsupported solution archive format version"))

    pdat = archive.problem
    pdat.term_names == [:weber, :zollner] ||
        throw(ArgumentError("cannot reconstruct archived custom Hamiltonian system"))

    system = WeberElectrodynamics.HamiltonianSystem(pdat.n_particles, pdat.dims)
    prob = WeberElectrodynamics.HamiltonianProblem(
        system,
        pdat.tspan,
        Vector{Float64}(pdat.q_initial),
        Vector{Float64}(pdat.p_initial),
        Vector{Float64}(pdat.params),
        Vector{Float64}(pdat.kappas),
        Float64(pdat.dt),
        Float64(pdat.convergence_tolerance),
        Int(pdat.maximum_iterations),
    )

    sdat = archive.solution
    t = Vector{Float64}(sdat.t)
    q = [Vector{Float64}(state) for state in sdat.q]
    p = [Vector{Float64}(state) for state in sdat.p]
    diagnostics =
        _diagnostics_from_archive(sdat.regularization, length(sdat.regularization.mode_history))

    return WeberElectrodynamics.HamiltonianSolution(
        t,
        q,
        p,
        prob,
        sdat.retcode,
        diagnostics,
    )
end

end
