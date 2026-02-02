###############  EGM diagnostics toolbox (copy-paste) ###############
using Optim
using Plots

# -------------------------
# Helpers: age indices
# -------------------------
ageidx_retminus1(p) = p.age_ret - p.age_min
ageidx_retminus2(p) = p.age_ret - p.age_min - 1
agelabel(p, age_i::Int) = p.age_min + age_i - 1

# -------------------------
# Compute EV_next / c_next_CE at ret-1 (N==3 branch of fill_EV_Euc!)
#   V_next uses n=1 slice after retirement
# -------------------------
function compute_EV_CE_retminus1!(variables, parameters)
    age_i = ageidx_retminus1(parameters)
    EV_next   = variables.EGM_ws.EV_next
    c_next_CE = variables.EGM_ws.c_next_CE

    V_next = @view variables.V[:, 1, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, 1, :, :, 2, age_i+1]

    fill_EV_Euc!(EV_next, c_next_CE, V_next, c_next, parameters, age_i)
    return nothing
end

# -------------------------
# Compute EV_next / c_next_CE at ret-2 (N==4 branch of fill_EV_Euc!)
#   V_next keeps n dimension, inf fixed at 2
# -------------------------
function compute_EV_CE_retminus2!(variables, parameters)
    age_i = ageidx_retminus2(parameters)
    EV_next   = variables.EGM_ws.EV_next
    c_next_CE = variables.EGM_ws.c_next_CE

    V_next = @view variables.V[:, :, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, :, :, :, 2, age_i+1]

    fill_EV_Euc!(EV_next, c_next_CE, V_next, c_next, parameters, age_i)
    return nothing
end

# -------------------------
# Diagnostic 1: ret-1 EV/CE should be n-invariant (because after retirement V_next has no n)
# -------------------------
function check_retminus1_EV_n_invariance(variables, parameters; verbose=true)
    compute_EV_CE_retminus1!(variables, parameters)
    EV_next   = variables.EGM_ws.EV_next
    c_next_CE = variables.EGM_ws.c_next_CE
    @unpack ϵ_size, n_size = parameters

    if verbose
        println("== Check n-invariance of EV and CE at age_ret-1 ==")
    end

    for ϵ_i in 1:ϵ_size
        baseEV = @view EV_next[:, 1, ϵ_i, 2]
        baseCE = @view c_next_CE[:, 1, ϵ_i, 2]

        maxΔEV = 0.0
        maxΔCE = 0.0
        for n_i in 2:n_size
            dEV = maximum(abs.(@view(EV_next[:, n_i, ϵ_i, 2]) .- baseEV))
            dCE = maximum(abs.(@view(c_next_CE[:, n_i, ϵ_i, 2]) .- baseCE))
            maxΔEV = max(maxΔEV, dEV)
            maxΔCE = max(maxΔCE, dCE)
        end

        if verbose
            println("ϵ_i=$ϵ_i: max |ΔEV| across n = $maxΔEV, max |ΔCE| across n = $maxΔCE")
        end
    end
    return nothing
end

# -------------------------
# Diagnostic 2: budget residuals at ret-1 (show ALL vs VALID)
# VALID means not penalty point
# -------------------------
function budget_check_retminus1(variables, parameters; ϵ_i=4, ν_i=3, tol=1e-9)
    age_i = ageidx_retminus1(parameters)
    @unpack a_size, n_size, aR_grid, w_grid, c_floor, V_penalty = parameters

    V_cur  = @view variables.V[:, :, ϵ_i, ν_i, 2, age_i]
    c_cur  = @view variables.policy_c[:, :, ϵ_i, ν_i, 2, age_i]
    ap_cur = @view variables.policy_a_p[:, :, ϵ_i, ν_i, 2, age_i]
    e_cur  = @view variables.policy_e[:, :, ϵ_i, ν_i, 2, age_i]  # for n=0 should be 0

    w_bar = w_grid[ϵ_i, ν_i, age_i]

    max_all = -Inf
    arg_all = (0,0)
    max_valid = -Inf
    arg_valid = (0,0)
    n_valid_bad = 0
    n_penalty = 0

    for a_i in 1:a_size, n_i in 1:n_size
        M = aR_grid[a_i] + w_bar
        c  = c_cur[a_i, n_i]
        ap = ap_cur[a_i, n_i]
        e  = e_cur[a_i, n_i]
        res = M - ap - c - e

        absres = abs(res)
        if absres > max_all
            max_all = absres
            arg_all = (a_i, n_i)
        end

        V = V_cur[a_i, n_i]
        is_penalty = (V <= (V_penalty/2)) || (!isfinite(V)) || (c <= c_floor) || (!isfinite(c))
        if is_penalty
            n_penalty += 1
        else
            if absres > max_valid
                max_valid = absres
                arg_valid = (a_i, n_i)
            end
            if absres > tol
                n_valid_bad += 1
            end
        end
    end

    println("== age=$(agelabel(parameters, age_i)), (ϵ,ν)=($ϵ_i,$ν_i) ret-1 budget check ==")
    println("max |budget residual| (ALL points)    = $max_all at (a_i,n_i)=$arg_all")
    println("max |budget residual| (VALID points)  = $max_valid at (a_i,n_i)=$arg_valid")
    println("count of VALID points with |res|>tol  = $n_valid_bad")
    println("count of NON-VALID (penalty) points   = $n_penalty")

    # print worst ALL point detail
    a0, n0 = arg_all
    M0  = aR_grid[a0] + w_bar
    println("\n[Worst ALL] (a_i=$a0,n_i=$n0): M=$M0, c=$(c_cur[a0,n0]), ap=$(ap_cur[a0,n0]), e=$(e_cur[a0,n0]), V=$(V_cur[a0,n0])")
    return nothing
end

# -------------------------
# ret-2 feasibility threshold:
# necessary & sufficient for existence of (ap>=a_min, e>=e_bar, c>=c_floor):
#   M(a) - a_min >= e_bar + c_floor
# => a_star = (e_bar + c_floor + a_min - w_bar)/R
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
# Traditional solver at ret-2 for one state (a,n,ϵ,ν):
# maximize u(c,e)+β*EV_next(ap) over ap in [a_min, M-e_bar-c_floor]
# where e is solved by your solve_e_bisect (n>0), or e=0 if n=0.
# -------------------------
using Optim

# -------------------------
# Traditional solver at ret-2 for one state (a,n,ϵ,ν):
# maximize u(c,e)+β*EV_next(ap) over ap in [a_min, M-e_bar-c_floor]
# with HARD feasibility for continuation value:
#   ap < ap_floor  -> V_penalty (no interpolation across penalty region)
# -------------------------
function opt_state_retminus2(
    M::Float64,
    n::Float64,
    P::Float64,
    e_bar::Float64,
    EV_next_a::AbstractVector{Float64},
    parameters;
    Ngrid::Int=80,
    refine::Bool=true,
    # optional: if you also want to require next-period c policy to be valid, pass it in
    c_next_a::Union{Nothing,AbstractVector{Float64}}=nothing,
)
    @unpack a_grid, a_min, β, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, ψ, γ, κ, c_floor, V_penalty = parameters

    # 1) Upper bound from within-period feasibility
    ap_hi = M - e_bar - c_floor
    if !(ap_hi > a_min) || !isfinite(ap_hi)
        return (ok=false, V=V_penalty, c=c_floor, ap=a_min, e=e_bar)
    end

    # 2) Build "valid continuation" mask (match your EGM skipping logic as much as possible)
    valid = isfinite.(EV_next_a) .& (EV_next_a .> V_penalty)
    if c_next_a !== nothing
        valid .&= isfinite.(c_next_a) .& (c_next_a .> c_floor)
    end

    i0 = findfirst(valid)
    if i0 === nothing
        # no feasible continuation anywhere
        return (ok=false, V=V_penalty, c=c_floor, ap=a_min, e=e_bar)
    end

    ap_floor = a_grid[i0]
    # IMPORTANT: effective lower bound is max(a_min, ap_floor)
    ap_lo = max(a_min, ap_floor)

    if !(ap_hi > ap_lo)
        return (ok=false, V=V_penalty, c=c_floor, ap=ap_lo, e=e_bar)
    end

    # 3) Interpolate continuation ONLY on valid tail
    aV = @view a_grid[i0:end]
    VV = @view EV_next_a[i0:end]

    @inline function Vcont(ap::Float64)
        # hard cut: do NOT interpolate across penalty region
        if ap < ap_floor
            return V_penalty
        end
        Vp = lininterp(aV, VV, ap)
        if !isfinite(Vp) || Vp <= V_penalty
            return V_penalty
        end
        return Vp
    end

    # 4) objective in ap (return very negative finite value, avoid -Inf in Brent)
    NEG = -1.0e300
    function obj(ap::Float64)
        if !(ap_lo <= ap <= ap_hi)
            return NEG
        end
        m = M - ap

        if n == 0.0
            c = m
            if !isfinite(c) || c <= c_floor
                return NEG
            end
            Vp = Vcont(ap)
            if Vp <= V_penalty
                return NEG
            end
            return u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β*Vp
        else
            if m <= e_bar + c_floor
                return NEG
            end
            e = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
            c = m - e
            if !isfinite(e) || e < e_bar || !isfinite(c) || c <= c_floor
                return NEG
            end
            Vp = Vcont(ap)
            if Vp <= V_penalty
                return NEG
            end
            scale = (n / P)^one_m_κ
            return u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
                            c_floor=c_floor, V_penalty=V_penalty) + β*Vp
        end
    end

    # 5) coarse grid search
    bestV  = NEG
    bestap = ap_lo
    for k in 0:Ngrid
        ap = ap_lo + (ap_hi - ap_lo) * (k / Ngrid)
        V = obj(ap)
        if V > bestV
            bestV = V
            bestap = ap
        end
    end
    if !(bestV > V_penalty) || !isfinite(bestV)
        return (ok=false, V=V_penalty, c=c_floor, ap=ap_lo, e=e_bar)
    end

    # 6) refine with Brent in a LOCAL bracket (ensure endpoints are feasible-ish)
    ap_star = bestap
    V_star  = bestV
    if refine
        step = (ap_hi - ap_lo) / max(Ngrid, 10)
        lo = max(ap_lo, bestap - 5step)
        hi = min(ap_hi, bestap + 5step)

        # safety: ensure obj(lo), obj(hi) are finite (not NEG)
        # if not, shrink around bestap until they are
        for _ in 1:8
            if obj(lo) > NEG/10 && obj(hi) > NEG/10
                break
            end
            lo = max(ap_lo, bestap - 0.5*(bestap-lo))
            hi = min(ap_hi, bestap + 0.5*(hi-bestap))
        end

        res = optimize(ap -> -obj(ap), lo, hi, Brent())
        ap_star = Optim.minimizer(res)
        V_star  = -Optim.minimum(res)
    end

    # 7) recover c,e at optimum
    m = M - ap_star
    if n == 0.0
        c_star = m
        e_star = 0.0
    else
        e_star = solve_e_bisect(m, n, P, ψ, γ, κ, one_m_κ, e_bar; c_floor=c_floor)
        c_star = m - e_star
    end

    ok = isfinite(V_star) && (V_star > V_penalty)
    return (ok=ok, V=V_star, c=c_star, ap=ap_star, e=e_star, ap_floor=ap_floor, ap_lo=ap_lo, ap_hi=ap_hi)
end

# -------------------------
# Compare EGM vs Traditional for ret-2, one (a_i,n_i,ϵ_i,ν_i)
# -------------------------
function compare_one_state_retminus2(variables, parameters; ϵ_i=4, ν_i=3, a_i=22, n_val=1.0, Ngrid=120, refine=true)
    compute_EV_CE_retminus2!(variables, parameters) # ensure EV_next filled

    age_i = ageidx_retminus2(parameters)
    @unpack a_grid, aR_grid, n_grid, w_grid, P_grid, q_bar_P_grid, V_penalty = parameters

    n_i = findfirst(==(n_val), n_grid)
    n_i === nothing && error("n_val=$n_val not in n_grid=$n_grid")

    # EGM objects at ret-2 (inf fixed=2)
    V_egm  = variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    c_egm  = variables.policy_c[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    ap_egm = variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    e_egm  = variables.policy_e[a_i, n_i, ϵ_i, ν_i, 2, age_i]

    # state info
    w_bar = w_grid[ϵ_i, ν_i, age_i]
    M     = aR_grid[a_i] + w_bar
    P     = P_grid[n_i, ϵ_i, ν_i, age_i]
    e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    n     = n_grid[n_i]

    # continuation (from workspace)
    EV_next_a = @view variables.EGM_ws.EV_next[:, n_i, ϵ_i, 2]

    sol = opt_state_retminus2(M, n, P, e_bar, EV_next_a, parameters; Ngrid=Ngrid, refine=refine)

    # diffs
    dV  = abs(V_egm - sol.V)
    dc  = abs(c_egm - sol.c)
    dap = abs(ap_egm - sol.ap)
    de  = abs(e_egm - sol.e)

    ok = sol.ok && (V_egm > V_penalty/2)

    println("== age=$(agelabel(parameters, age_i)), n=$n_val, (ϵ,ν,a)=($ϵ_i,$ν_i,$a_i) ret-2 ==")
    println("ok = $ok")
    println("|ΔV|  = $dV")
    println("|Δc|  = $dc")
    println("|Δap| = $dap")
    println("|Δe|  = $de")

    return (ok=ok, dV=dV, dc=dc, dap=dap, de=de, egm=(V_egm,c_egm,ap_egm,e_egm), trad=sol)
end

# -------------------------
# Build traditional a'(a) for ret-2 slice at given (ϵ,ν,n), and plot with:
# - NaN out infeasible region
# - vline at a_star
# -------------------------
ageidx_retminus2(par) = (par.age_ret - par.age_min - 1)
agelabel(par, age_i::Int) = par.age_min + age_i - 1

function compute_EV_CE_retminus2!(variables, parameters; age_i::Int=ageidx_retminus2(parameters))
    V_next = @view variables.V[:, :, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, :, :, :, 2, age_i+1]
    fill_EV_Euc!(variables.EGM_ws.EV_next, variables.EGM_ws.c_next_CE, V_next, c_next, parameters, age_i)
    return nothing
end

function retminus2_a_star(parameters; ϵ_i::Int, ν_i::Int, n_i::Int)
    age_i = ageidx_retminus2(parameters)
    w_bar = parameters.w_grid[ϵ_i, ν_i, age_i]
    e_bar = parameters.q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    # minimal a such that: a*R + w_bar - a_min >= e_bar + c_floor
    a_star = (e_bar + parameters.c_floor + parameters.a_min - w_bar) / parameters.R
    return a_star, w_bar, e_bar, age_i
end

function plot_retminus2_ap_compare(
    variables, parameters;
    ϵ_i::Int=4,
    ν_i::Int=3,
    n_val::Float64=1.0,
    Ngrid::Int=120,
    refine::Bool=true,
    show_infeasible::Bool=true,
    savepath::Union{Nothing,String}=nothing,   # <- 新增：要存檔就給路徑
    do_display::Bool=false                      # <- 新增：要不要 display
)
    compute_EV_CE_retminus2!(variables, parameters)

    age_i = ageidx_retminus2(parameters)
    age_lbl = agelabel(parameters, age_i)

    @unpack a_grid, aR_grid, a_min, c_floor, n_grid, w_grid, P_grid, q_bar_P_grid = parameters

    n_i = findfirst(==(n_val), n_grid)
    n_i === nothing && error("n_val=$n_val not in n_grid=$n_grid")

    # EGM policy (ret-2, inf fixed=2)
    ap_egm = collect(@view variables.policy_a_p[:, n_i, ϵ_i, ν_i, 2, age_i])

    # constants for this slice
    w_bar = w_grid[ϵ_i, ν_i, age_i]
    P     = P_grid[n_i, ϵ_i, ν_i, age_i]
    e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    n     = n_grid[n_i]

    # feasibility
    Mvec = aR_grid .+ w_bar
    feas = (Mvec .- a_min) .>= (e_bar .+ c_floor)

    # a_star vertical line
    a_star, _, _, _ = retminus2_a_star(parameters; ϵ_i=ϵ_i, ν_i=ν_i, n_i=n_i)

    # continuation values
    EV_next_a = @view variables.EGM_ws.EV_next[:, n_i, ϵ_i, 2]

    # traditional ap(a)
    ap_trad = fill(NaN, length(a_grid))
    for idx in eachindex(a_grid)
        if !feas[idx]
            continue
        end
        M = Mvec[idx]
        sol = opt_state_retminus2(M, n, P, e_bar, EV_next_a, parameters; Ngrid=Ngrid, refine=refine)
        ap_trad[idx] = sol.ap
    end

    # NaN out infeasible for EGM too (avoid misleading 0 line)
    ap_egm_plot = copy(ap_egm)
    ap_egm_plot[.!feas] .= NaN

    title_txt = "age=ret-2 (age=$age_lbl), n=$n_val, (ϵ,ν)=($ϵ_i,$ν_i): a' policy"

    p = plot(a_grid, ap_egm_plot, label="EGM  a'(a)", lw=3,
             title=title_txt, xlabel="a", ylabel="a'")
    plot!(p, a_grid, ap_trad, label="Traditional (Optim)  a'(a)", lw=3, ls=:dash)
    vline!(p, [a_star], label="min feasible a (M-a_min ≥ e_bar+c_floor)", ls=:dot, lw=2)

    if show_infeasible
        y0 = fill(a_min, length(a_grid))
        y0[feas] .= NaN
        scatter!(p, a_grid, y0, label="infeasible states", ms=3)
    end

    println("== ret-2 slice summary ==")
    println("age=$age_lbl, n=$n_val, (ϵ,ν)=($ϵ_i,$ν_i)")
    println("w_bar=$w_bar, e_bar=$e_bar, a_star=$a_star")
    println("feasible count = $(count(feas)) / $(length(feas))")

    if savepath !== nothing
        savefig(p, savepath)
        println("Saved figure to: $savepath")
    end
    if do_display
        display(p)
    end

    return (p=p, ap_trad=ap_trad, feas=feas, a_star=a_star)
end

function plot_retminus2_ap_compare_v2(
    variables, parameters;
    ϵ_i::Int=4,
    ν_i::Int=3,
    n_val::Float64=1.0,
    Ngrid::Int=120,
    refine::Bool=true,
    tol_ap::Float64=1e-10,
    savepath::Union{Nothing,String}=joinpath(pwd(), "ret2_ap_policy_compare.png"),
    do_display::Bool=false
)
    compute_EV_CE_retminus2!(variables, parameters)

    @unpack a_grid, aR_grid, a_min, c_floor, n_grid, w_grid, P_grid, q_bar_P_grid = parameters
    age_i = ageidx_retminus2(parameters)
    age_lbl = agelabel(parameters, age_i)

    n_i = findfirst(==(n_val), n_grid)
    n_i === nothing && error("n_val=$n_val not in n_grid=$n_grid")

    # slice constants
    w_bar = w_grid[ϵ_i, ν_i, age_i]
    P     = P_grid[n_i, ϵ_i, ν_i, age_i]
    e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]

    # feasibility wrt minimal required spending (c_floor + e_bar) at ap=a_min
    Mvec = aR_grid .+ w_bar
    feas = (Mvec .- a_min) .>= (e_bar .+ c_floor)

    # EGM policies
    ap_egm = collect(@view variables.policy_a_p[:, n_i, ϵ_i, ν_i, 2, age_i])

    # classify: feasible but binding ap=a_min
    bind = falses(length(a_grid))
    @inbounds for i in eachindex(a_grid)
        bind[i] = feas[i] && (abs(ap_egm[i] - a_min) <= tol_ap)
    end

    # a_star line (feasibility threshold)
    a_star, _, _, _ = retminus2_a_star(parameters; ϵ_i=ϵ_i, ν_i=ν_i, n_i=n_i)

    # "leave-the-boundary" threshold: first a where feasible and ap_egm > a_min
    a_bind = NaN
    for i in eachindex(a_grid)
        if feas[i] && (ap_egm[i] > a_min + tol_ap)
            a_bind = a_grid[i]
            break
        end
    end

    # continuation EV (as function of ap grid)
    EV_next_a = @view variables.EGM_ws.EV_next[:, n_i, ϵ_i, 2]

    # traditional ap(a) from your optimizer/grid method
    ap_trad = fill(NaN, length(a_grid))
    for idx in eachindex(a_grid)
        if !feas[idx]
            continue
        end
        M = Mvec[idx]
        sol = opt_state_retminus2(M, n_val, P, e_bar, EV_next_a, parameters; Ngrid=Ngrid, refine=refine)
        ap_trad[idx] = sol.ap
    end

    # plot-ready series (avoid drawing EGM on infeasible)
    ap_egm_plot = copy(ap_egm)
    ap_egm_plot[.!feas] .= NaN

    title_txt = "age=ret-2 (age=$age_lbl), n=$n_val, (ϵ,ν)=($ϵ_i,$ν_i): a' policy"

    p = plot(a_grid, ap_egm_plot, label="EGM  a'(a)", lw=3,
             title=title_txt, xlabel="a", ylabel="a'")
    plot!(p, a_grid, ap_trad, label="Traditional (Optim)  a'(a)", lw=3, ls=:dash)

    vline!(p, [a_star], label="feasibility a* (M-a_min ≥ e_bar+c_floor)", ls=:dot, lw=2)
    if isfinite(a_bind)
        vline!(p, [a_bind], label="EGM leaves a'=a_min", ls=:dashdot, lw=2)
    end

    # mark binding points (feasible but ap=a_min)
    scatter!(p, a_grid[bind], ap_egm_plot[bind], label="feasible & binding (a'=a_min)", ms=3)

    # put infeasible points BELOW axis to avoid confusion with ap=0 plateau
    infeas_idx = findall(!, feas)
    if !isempty(infeas_idx)
        y_low = -0.05 * maximum(skipmissing(vcat(ap_egm_plot, ap_trad)))
        scatter!(p, a_grid[infeas_idx], fill(y_low, length(infeas_idx)), label="infeasible", ms=3)
    end

    println("== ret-2 plot diagnostics ==")
    println("age=$age_lbl, n=$n_val, (ϵ,ν)=($ϵ_i,$ν_i)")
    println("w_bar=$w_bar, e_bar=$e_bar, a_star=$a_star, a_bind=$a_bind")
    println("feasible count = $(count(feas)) / $(length(feas))")
    println("binding count  = $(count(bind))")

    if savepath !== nothing
        savefig(p, savepath)
        println("Saved figure to: $savepath")
    end
    if do_display
        display(p)
    end

    return (p=p, ap_trad=ap_trad, feas=feas, bind=bind, a_star=a_star, a_bind=a_bind)
end


# -------------------------
# Convenience: run the “old checks” + make your ret-2 plot in one go
# -------------------------
function run_egm_checks!(variables, parameters; ϵ_i=4, ν_i=3, n_val=1.0)
    check_retminus1_EV_n_invariance(variables, parameters; verbose=true)
    budget_check_retminus1(variables, parameters; ϵ_i=ϵ_i, ν_i=ν_i)
    plot_retminus2_ap_compare(variables, parameters; ϵ_i=ϵ_i, ν_i=ν_i, n_val=n_val)
    return nothing
end

###############  End diagnostics toolbox ###############

out = plot_retminus2_ap_compare(variables, parameters;
    ϵ_i=5, ν_i=3, n_val=4.0,
    Ngrid=120, refine=true,
    savepath=joinpath(pwd(), "ret2_ap_policy_e4v3_n1.png"), 
    do_display=true
)



