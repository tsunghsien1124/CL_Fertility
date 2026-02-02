###############  Minimal EGM vs Traditional plot (ret-2) ###############
using Optim
using Plots
using Parameters: @unpack

# -------------------------
# Helpers: age index/label
# -------------------------
ageidx_retminus2(p) = p.age_ret - p.age_min - 1
agelabel(p, age_i::Int) = p.age_min + age_i - 1

# -------------------------
# Compute EV_next / c_next_CE at ret-2 (N==4 branch in your fill_EV_Euc!)
#   V_next keeps n dimension, inf fixed at 2
# -------------------------
function compute_EV_CE_retminus2!(variables, parameters; age_i::Int=ageidx_retminus2(parameters))
    V_next = @view variables.V[:, :, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, :, :, :, 2, age_i+1]
    fill_EV_Euc!(variables.EGM_ws.EV_next, variables.EGM_ws.c_next_CE, V_next, c_next, parameters, age_i)
    return nothing
end

# -------------------------
# ret-2 feasibility threshold (a_star):
# minimal a s.t. at ap=a_min we can afford e_bar + c_floor:
#   (a*R + w_bar) - a_min >= e_bar + c_floor
# -------------------------
function retminus2_a_star(parameters; ϵ_i::Int, ν_i::Int, n_i::Int)
    age_i = ageidx_retminus2(parameters)
    @unpack R, a_min, c_floor, w_grid, q_bar_P_grid = parameters
    w_bar = w_grid[ϵ_i, ν_i, age_i]
    e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    a_star = (e_bar + c_floor + a_min - w_bar) / R
    return a_star, w_bar, e_bar, age_i
end

# -------------------------
# Traditional solver at ret-2 for one state (given M, n, P, e_bar):
# maximize u(c,e)+β*EV_next(ap) over ap in [a_min, M-e_bar-c_floor]
# where e is solved by solve_e_bisect when n>0, or e=0 if n=0
# -------------------------
function opt_state_retminus2(
    M::Float64,
    n::Float64,
    P::Float64,
    e_bar::Float64,
    EV_next_a::AbstractVector{Float64},
    parameters;
    Ngrid::Int=120,
    refine::Bool=true,
)
    @unpack a_grid, a_min, β, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ = parameters
    @unpack ψ, γ, κ, c_floor, V_penalty = parameters

    ap_lo = a_min
    ap_hi = M - e_bar - c_floor
    if !(ap_hi > ap_lo) || !isfinite(ap_hi)
        return (ok=false, V=V_penalty, c=c_floor, ap=a_min, e=e_bar)
    end

    # minimization objective: f(ap) = -V(ap) for feasible, 1e30 otherwise
    function f(ap::Float64)
        if ap < ap_lo || ap > ap_hi
            return 1e30
        end

        m = M - ap
        if !isfinite(m)
            return 1e30
        end

        # continuation value
        Vp = lininterp(a_grid, EV_next_a, ap)
        if !isfinite(Vp) || Vp <= V_penalty
            return 1e30
        end

        if n == 0.0
            c = m
            if !isfinite(c) || c <= c_floor
                return 1e30
            end
            V = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β*Vp
            return -V
        else
            if m <= e_bar + c_floor
                return 1e30
            end
            e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
            if !isfinite(e) || e < e_bar
                return 1e30
            end
            c = m - e
            if !isfinite(c) || c <= c_floor
                return 1e30
            end
            scale = (n / P)^one_m_κ
            V = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
                         c_floor=c_floor, V_penalty=V_penalty) + β*Vp
            return -V
        end
    end

    # coarse grid search (minimize f)
    bestf  = 1e30
    bestap = ap_lo
    for k in 0:Ngrid
        ap = ap_lo + (ap_hi - ap_lo) * (k / Ngrid)
        fv = f(ap)
        if fv < bestf
            bestf  = fv
            bestap = ap
        end
    end
    if !(bestf < 1e29)  # still basically infeasible everywhere
        return (ok=false, V=V_penalty, c=c_floor, ap=a_min, e=e_bar)
    end

    ap_star = bestap
    f_star  = bestf

    if refine
        step = (ap_hi - ap_lo) / max(Ngrid, 10)
        lo = max(ap_lo, bestap - 5step)
        hi = min(ap_hi, bestap + 5step)
        res = optimize(f, lo, hi, Brent())
        ap_star = Optim.minimizer(res)
        f_star  = Optim.minimum(res)
    end

    # recover policies at optimum
    m = M - ap_star
    if n == 0.0
        e_star = 0.0
        c_star = m
    else
        e_star = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
        c_star = m - e_star
    end

    V_star = -f_star
    ok = isfinite(V_star) && V_star > V_penalty/2
    return (ok=ok, V=V_star, c=c_star, ap=ap_star, e=e_star)
end


# -------------------------
# Main plot: EGM vs Traditional a'(a) at ret-2 slice
# Includes:
#   - vline at feasibility a_star
#   - vline where EGM leaves binding a'=a_min (optional)
#   - scatter infeasible points (plotted at a low y-level)
# -------------------------
using Optim
using Plots
using Parameters: @unpack

# age helpers
ageidx_retminus2(p) = p.age_ret - p.age_min - 1
agelabel(p, age_i::Int) = p.age_min + age_i - 1

"""
plot_retminus2_compare_all

Produces 3 plots for a given (ϵ_i, ν_i, n_val) at ret-2 slice:
  1) a'(a)
  2) e(a)
  3) V(a)

Also shows:
  - feasibility threshold a_star
  - a_bind where EGM leaves binding a'=a_min
  - infeasible points as scatter at low y level
"""
function plot_retminus2_compare_all(
    variables, parameters;
    ϵ_i::Int=4,
    ν_i::Int=3,
    n_val::Float64=1.0,
    Ngrid::Int=200,
    refine::Bool=true,
    tol_ap::Float64=1e-10,
    saveprefix::Union{Nothing,String}=nothing,  # e.g. joinpath(pwd(),"ret2_e4v3_n2")
    do_display::Bool=true,
)
    # ensure EV_next exists (needed by traditional solver)
    compute_EV_CE_retminus2!(variables, parameters)

    @unpack a_grid, aR_grid, a_min, c_floor, n_grid, w_grid, P_grid, q_bar_P_grid, V_penalty = parameters
    age_i   = ageidx_retminus2(parameters)
    age_lbl = agelabel(parameters, age_i)

    n_i = findfirst(==(n_val), n_grid)
    n_i === nothing && error("n_val=$n_val not in n_grid=$n_grid")

    # slice constants
    w_bar = w_grid[ϵ_i, ν_i, age_i]
    P     = P_grid[n_i, ϵ_i, ν_i, age_i]
    e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    n     = n_grid[n_i]

    # feasibility wrt minimal required spending at ap=a_min:
    # M - a_min >= e_bar + c_floor
    Mvec = aR_grid .+ w_bar
    feas = (Mvec .- a_min) .>= (e_bar .+ c_floor)

    # EGM slices
    ap_egm = collect(@view variables.policy_a_p[:, n_i, ϵ_i, ν_i, 2, age_i])
    e_egm  = collect(@view variables.policy_e[:,   n_i, ϵ_i, ν_i, 2, age_i])
    V_egm  = collect(@view variables.V[:,          n_i, ϵ_i, ν_i, 2, age_i])
    c_egm  = collect(@view variables.policy_c[:,   n_i, ϵ_i, ν_i, 2, age_i])

    # binding points (feasible & ap=a_min)
    bind = falses(length(a_grid))
    @inbounds for i in eachindex(a_grid)
        bind[i] = feas[i] && isfinite(ap_egm[i]) && (abs(ap_egm[i] - a_min) <= tol_ap)
    end

    # a_star line
    a_star, _, _, _ = retminus2_a_star(parameters; ϵ_i=ϵ_i, ν_i=ν_i, n_i=n_i)

    # first a where feasible and EGM leaves binding (ap > a_min)
    a_bind = NaN
    for i in eachindex(a_grid)
        if feas[i] && isfinite(ap_egm[i]) && (ap_egm[i] > a_min + tol_ap)
            a_bind = a_grid[i]
            break
        end
    end

    # continuation EV (as function of ap grid)
    EV_next_a = @view variables.EGM_ws.EV_next[:, n_i, ϵ_i, 2]

    # Traditional slices
    ap_trad = fill(NaN, length(a_grid))
    e_trad  = fill(NaN, length(a_grid))
    V_trad  = fill(NaN, length(a_grid))

    for idx in eachindex(a_grid)
        if !feas[idx]
            continue
        end
        M = Mvec[idx]
        sol = opt_state_retminus2(M, n, P, e_bar, EV_next_a, parameters; Ngrid=Ngrid, refine=refine)
        if sol.ok
            ap_trad[idx] = sol.ap
            e_trad[idx]  = sol.e
            V_trad[idx]  = sol.V
        end
    end

    # Plot-ready EGM series: NaN out infeasible to avoid misleading plateaus
    ap_egm_plot = copy(ap_egm); ap_egm_plot[.!feas] .= NaN
    e_egm_plot  = copy(e_egm);  e_egm_plot[.!feas]  .= NaN
    V_egm_plot  = copy(V_egm);  V_egm_plot[.!feas]  .= NaN

    title_base = "ret-2 (age=$age_lbl), n=$n_val, (ϵ,ν)=($ϵ_i,$ν_i)"

    # helper: draw infeasible points at low y
    function add_infeasible_scatter!(p, x, feasmask, series_list; label="infeasible", ms=3)
        infeas_idx = findall(!, feasmask)
        isempty(infeas_idx) && return p

        finite_vals = Float64[]
        for s in series_list
            append!(finite_vals, filter(isfinite, s))
        end
        ymax = isempty(finite_vals) ? 1.0 : maximum(finite_vals)
        y_low = -0.05 * ymax
        scatter!(p, x[infeas_idx], fill(y_low, length(infeas_idx)); label=label, ms=ms)
        return p
    end

    # -----------------
    # Plot 1: a'(a)
    # -----------------
    p_ap = plot(a_grid, ap_egm_plot; label="EGM a'(a)", lw=3,
                title="$title_base: a'(a)", xlabel="a", ylabel="a'")
    plot!(p_ap, a_grid, ap_trad; label="Traditional (Optim) a'(a)", lw=3, ls=:dash)
    vline!(p_ap, [a_star]; label="feasibility a*", ls=:dot, lw=2)
    if isfinite(a_bind)
        vline!(p_ap, [a_bind]; label="EGM leaves binding", ls=:dashdot, lw=2)
    end
    scatter!(p_ap, a_grid[bind], ap_egm_plot[bind]; label="feasible & binding (a'=a_min)", ms=3)
    add_infeasible_scatter!(p_ap, a_grid, feas, [ap_egm_plot, ap_trad])

    # -----------------
    # Plot 2: e(a)
    # -----------------
    p_e = plot(a_grid, e_egm_plot; label="EGM e(a)", lw=3,
               title="$title_base: e(a)", xlabel="a", ylabel="e")
    plot!(p_e, a_grid, e_trad; label="Traditional (Optim) e(a)", lw=3, ls=:dash)
    hline!(p_e, [e_bar]; label="e_bar", ls=:dot, lw=2)
    vline!(p_e, [a_star]; label="feasibility a*", ls=:dot, lw=2)
    if isfinite(a_bind)
        vline!(p_e, [a_bind]; label="EGM leaves binding", ls=:dashdot, lw=2)
    end
    add_infeasible_scatter!(p_e, a_grid, feas, [e_egm_plot, e_trad])

    # -----------------
    # Plot 3: V(a)
    # -----------------
    p_V = plot(a_grid, V_egm_plot; label="EGM V(a)", lw=3,
               title="$title_base: V(a)", xlabel="a", ylabel="V")
    plot!(p_V, a_grid, V_trad; label="Traditional (Optim) V(a)", lw=3, ls=:dash)
    vline!(p_V, [a_star]; label="feasibility a*", ls=:dot, lw=2)
    if isfinite(a_bind)
        vline!(p_V, [a_bind]; label="EGM leaves binding", ls=:dashdot, lw=2)
    end
    add_infeasible_scatter!(p_V, a_grid, feas, [V_egm_plot, V_trad])

    println("== ret-2 compare-all diagnostics ==")
    println("age=$age_lbl, n=$n_val, (ϵ,ν)=($ϵ_i,$ν_i)")
    println("w_bar=$w_bar, e_bar=$e_bar, a_star=$a_star, a_bind=$a_bind")
    println("feasible count = $(count(feas)) / $(length(feas))")
    println("binding count  = $(count(bind))")

    if saveprefix !== nothing
        savefig(p_ap, saveprefix * "_ap.png")
        savefig(p_e,  saveprefix * "_e.png")
        savefig(p_V,  saveprefix * "_V.png")
        println("Saved: $(saveprefix)_ap.png, $(saveprefix)_e.png, $(saveprefix)_V.png")
    end

    if do_display
        display(p_ap); display(p_e); display(p_V)
    end

    return (p_ap=p_ap, p_e=p_e, p_V=p_V,
            ap_trad=ap_trad, e_trad=e_trad, V_trad=V_trad,
            feas=feas, bind=bind, a_star=a_star, a_bind=a_bind)
end

out = plot_retminus2_compare_all(variables, parameters;
    ϵ_i=1, ν_i=3, n_val=2.0,
    Ngrid=200, refine=true,
    saveprefix=joinpath(pwd(), "ret2_e4v3_n2"),
    do_display=true
)
