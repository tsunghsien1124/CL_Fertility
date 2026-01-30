using JLD2
using Plots

old = load("workspace_benchmark.jld2")
new = load("workspace_benchmark_new.jld2")

@assert all(old["V"][:,1,:,:,2,end] .≈ new["V"][:,1,:,:,2,end])

n_ind = 1
ϵ_ind = 4
ν_ind = 3
inf_ind = 2
age_ret_ind = old["parameters"].age_ret - old["parameters"].age_min + 1
# age_ret_ind = old["parameters"].age_max - old["parameters"].age_min

plot(old["parameters"].a_grid, old["parameters"].a_grid[old["policy_a_p"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_ret_ind]], label="old")
plot!(new["parameters"].a_grid, new["policy_a_p"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_ret_ind], label="new")

plot(old["parameters"].a_grid, old["V"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_ret_ind], label="old")
plot!(new["parameters"].a_grid, new["V"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_ret_ind], label="new")

n_ind = 2
ϵ_ind = 4
ν_ind = 3
inf_ind = 2
# age_inf_ind = old["parameters"].age_inf - old["parameters"].age_min + 1
age_inf_ind = old["parameters"].age_ret - old["parameters"].age_min

plot(old["parameters"].a_grid, old["parameters"].a_grid[old["policy_a_p"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_inf_ind]], label="old")
plot!(new["parameters"].a_grid, new["policy_a_p"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_inf_ind], label="new")

plot(old["parameters"].a_grid, old["V"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_inf_ind], label="old")
plot!(new["parameters"].a_grid, new["V"][:,n_ind,ϵ_ind,ν_ind,inf_ind,age_inf_ind], label="new")

w_old = exp(old["parameters"].ϵ_grid[ϵ_ind] + old["parameters"].ν_grid[ν_ind] + old["parameters"].h_grid[age_inf_ind])
w_new = exp(new["parameters"].ϵ_grid[ϵ_ind] + new["parameters"].ν_grid[ν_ind] + new["parameters"].h_grid[age_inf_ind])
@assert w_old == w_new

a_ind = 22 # old["parameters"].a_size -1 
a_old = old["parameters"].a_grid[a_ind]
a_p_old = old["parameters"].a_grid[old["policy_a_p"][a_ind,n_ind,ϵ_ind,ν_ind,inf_ind,age_inf_ind]]
x_old = old["parameters"].x_grid[old["policy_x"][a_ind,n_ind,ϵ_ind,ν_ind,inf_ind,age_inf_ind]]
l_old = old["parameters"].l_grid[old["policy_l"][a_ind,n_ind,ϵ_ind,ν_ind,inf_ind,age_inf_ind]]
CoH_old = a_old * (1.0 + old["parameters"].r) + (1.0 - l_old) * w_old - a_p_old - x_old
n_old = old["parameters"].n_grid[n_ind]

function quality_function(x::Real, l::Real, n::Integer, μ::Real, θ::Real, ψ_1::Real, ψ_2::Real)
    """
    child quality function
        x - time input
        l - monetary input
        n - number of kids
    """
    if n > 0
        return (μ * (x / (n^ψ_1))^θ + (1.0 - μ) * (l / (n^ψ_2))^θ)^(1.0 / θ)
    else
        return 0.0
    end
end

@assert quality_function(x_old, l_old, n_old, old["parameters"].μ, old["parameters"].θ, old["parameters"].ψ_1, old["parameters"].ψ_2) > old["parameters"].q_bar

function budget_check_retminus1_feasible(variables::Mutable_Variables, parameters::NamedTuple;
    ϵ_i::Int=4, ν_i::Int=3, tol::Float64=1e-9)

    @unpack age_ret, age_min, a_size, n_size = parameters
    @unpack aR_grid, a_min, w_grid, P_grid, q_bar_P_grid = parameters
    @unpack c_floor, V_penalty, one_m_κ, ψ, κ, γ = parameters
    @unpack n_grid = parameters

    age_i = age_ret - age_min
    w_bar = w_grid[ϵ_i, ν_i, age_i]

    V_cur  = @view variables.V[:, :, ϵ_i, ν_i, 2, age_i]
    c_cur  = @view variables.policy_c[:, :, ϵ_i, ν_i, 2, age_i]
    ap_cur = @view variables.policy_a_p[:, :, ϵ_i, ν_i, 2, age_i]
    e_cur  = @view variables.policy_e[:, :, ϵ_i, ν_i, 2, age_i]

    max_budget_valid = 0.0
    max_budget_all   = 0.0
    bad_valid = 0
    bad_all   = 0

    argmax_all = (0,0)
    argmax_valid = (0,0)

    for n_i in 1:n_size
        P     = P_grid[n_i, ϵ_i, ν_i, age_i]
        e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
        n     = n_grid[n_i]
        scale = (n / P)^one_m_κ

        for a_i in 1:a_size
            M  = aR_grid[a_i] + w_bar
            V  = V_cur[a_i, n_i]
            c  = c_cur[a_i, n_i]
            ap = ap_cur[a_i, n_i]
            e  = e_cur[a_i, n_i]

            # budget residual (你原本的定義)
            bud = (n_i == 1) ? (M - (c + ap)) : (M - (c + ap + e))

            # 先算 all-points 最大（會被 penalty 點主宰）
            if abs(bud) > max_budget_all
                max_budget_all = abs(bud)
                argmax_all = (a_i, n_i)
            end

            # 判定 “這個點是否是你真的解出來的可行最適解”
            # (1) 不用 V == V_penalty 的點
            # (2) c 必須 > c_floor
            # (3) n>0 時，e 必須 >= e_bar（給一點容忍）
            is_valid = isfinite(V) && (V > V_penalty/2) && isfinite(c) && (c > c_floor)
            if n_i > 1
                is_valid &= isfinite(e) && (e >= e_bar - 1e-10)
                # 若你也想順便看內點 FOC
                # if e > e_bar + 1e-10
                #     intra = abs(ψ*scale*e^(-κ) - c^(-γ))
                # end
            end

            if is_valid
                if abs(bud) > max_budget_valid
                    max_budget_valid = abs(bud)
                    argmax_valid = (a_i, n_i)
                end
                if abs(bud) > tol
                    bad_valid += 1
                end
            else
                bad_all += 1
            end
        end
    end

    println("== age=$(age_ret-1), ϵ_i=$ϵ_i, ν_i=$ν_i ==")
    println("max |budget residual| (ALL points)    = $max_budget_all at (a_i,n_i)=$argmax_all")
    println("max |budget residual| (VALID points)  = $max_budget_valid at (a_i,n_i)=$argmax_valid")
    println("count of VALID points with |res|>tol  = $bad_valid")
    println("count of NON-VALID (penalty) points   = $bad_all")

    # 印出造成最大殘差的點看一下到底發生什麼
    let (ai, ni) = argmax_all
        P     = P_grid[ni, ϵ_i, ν_i, age_i]
        e_bar = q_bar_P_grid[ni, ϵ_i, ν_i, age_i]
        M  = aR_grid[ai] + w_bar
        c  = c_cur[ai, ni]
        ap = ap_cur[ai, ni]
        e  = e_cur[ai, ni]
        println("\n[Worst ALL] a_i=$ai n_i=$ni: M=$M, c=$c, ap=$ap, e=$e, e_bar=$e_bar, V=$(V_cur[ai,ni])")
        println("Feasibility check at ap=a_min: M-a_min = $(M-a_min) vs e_bar+c_floor = $(e_bar + c_floor)")
    end
end

function build_EV_retminus1!(EV_tmp, cCE_tmp, variables, parameters)
    @unpack age_ret, age_min = parameters
    age_i = age_ret - age_min  # age_ret-1 index

    V_next = @view variables.V[:, 1, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, 1, :, :, 2, age_i+1]
    fill_EV_Euc!(EV_tmp, cCE_tmp, V_next, c_next, parameters, age_i)
    return age_i
end

function brute_force_retminus1_state(
    M::Float64, n::Float64, P::Float64, e_bar::Float64,
    EV_next_a::AbstractVector{Float64},
    parameters::NamedTuple;
    apN::Int = 400, eN::Int = 400,
)
    @unpack a_grid, a_min, a_max, β = parameters
    @unpack one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ = parameters
    @unpack c_floor, V_penalty = parameters

    # feasible upper bound for a'
    ap_hi = min(a_max, M - c_floor - (n > 0 ? e_bar : 0.0))
    if ap_hi < a_min
        return (V_penalty, c_floor, a_min, (n > 0 ? e_bar : 0.0), false)
    end

    bestV  = V_penalty
    bestc  = c_floor
    bestap = a_min
    beste  = (n > 0 ? e_bar : 0.0)
    found  = false

    if n == 0.0
        for ap in LinRange(a_min, ap_hi, apN)
            Vp = lininterp(a_grid, EV_next_a, ap)
            if !isfinite(Vp) || Vp <= V_penalty
                continue
            end
            c = M - ap
            if c <= c_floor || !isfinite(c)
                continue
            end
            V = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β*Vp
            if V > bestV
                bestV = V; bestc = c; bestap = ap; found = true
            end
        end
        return (bestV, bestc, bestap, 0.0, found)
    end

    # n > 0
    scale = (n / P)^one_m_κ
    for ap in LinRange(a_min, ap_hi, apN)
        Vp = lininterp(a_grid, EV_next_a, ap)
        if !isfinite(Vp) || Vp <= V_penalty
            continue
        end

        e_hi = M - c_floor - ap
        if e_hi < e_bar
            continue
        end

        for e in LinRange(e_bar, e_hi, eN)
            c = M - ap - e
            if c <= c_floor || !isfinite(c)
                continue
            end
            V = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
                         c_floor=c_floor, V_penalty=V_penalty) + β*Vp
            if V > bestV
                bestV = V; bestc = c; bestap = ap; beste = e; found = true
            end
        end
    end

    return (bestV, bestc, bestap, beste, found)
end

function compare_EGM_vs_bruteforce_retminus1(
    variables::Mutable_Variables, parameters::NamedTuple;
    ϵ_i::Int=4, ν_i::Int=3,
    apN::Int=400, eN::Int=400,
    print_top::Int=10
)
    @unpack a_size, n_size, aR_grid, w_grid, P_grid, q_bar_P_grid, n_grid = parameters
    @unpack age_ret, age_min, V_penalty, c_floor = parameters

    # build EV at age_ret-1
    EV_tmp  = fill(V_penalty, a_size, n_size, parameters.ϵ_size, parameters.inf_size)
    cCE_tmp = fill(NaN,       a_size, n_size, parameters.ϵ_size, parameters.inf_size)
    age_i = build_EV_retminus1!(EV_tmp, cCE_tmp, variables, parameters)

    w_bar = w_grid[ϵ_i, ν_i, age_i]

    V_EGM  = @view variables.V[:, :, ϵ_i, ν_i, 2, age_i]
    c_EGM  = @view variables.policy_c[:, :, ϵ_i, ν_i, 2, age_i]
    ap_EGM = @view variables.policy_a_p[:, :, ϵ_i, ν_i, 2, age_i]
    e_EGM  = @view variables.policy_e[:, :, ϵ_i, ν_i, 2, age_i]

    # collect worst cases
    diffs = Vector{Tuple{Float64,Int,Int,Float64,Float64,Float64,Float64}}()
    # (ΔV, a_i, n_i, Δc, Δap, Δe, V_BF)

    for n_i in 1:n_size
        EV_next_a = @view EV_tmp[:, n_i, ϵ_i, 2]
        P     = P_grid[n_i, ϵ_i, ν_i, age_i]
        e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
        n     = n_grid[n_i]

        for a_i in 1:a_size
            M = aR_grid[a_i] + w_bar

            (Vbf, cbf, apbf, ebf, ok) =
                brute_force_retminus1_state(M, n, P, e_bar, EV_next_a, parameters; apN=apN, eN=eN)

            Vegm = V_EGM[a_i, n_i]
            if ok && isfinite(Vegm) && Vegm > V_penalty/2 && cbf > c_floor
                ΔV  = abs(Vegm - Vbf)
                Δc  = abs(c_EGM[a_i, n_i]  - cbf)
                Δap = abs(ap_EGM[a_i, n_i] - apbf)
                Δe  = abs(e_EGM[a_i, n_i]  - ebf)
                push!(diffs, (ΔV, a_i, n_i, Δc, Δap, Δe, Vbf))
            end
        end
    end

    sort!(diffs, by = x -> x[1], rev=true)

    println("== Compare EGM vs brute-force at age=$(age_ret-1), ϵ_i=$ϵ_i, ν_i=$ν_i ==")
    println("Checked feasible points: ", length(diffs))
    for k in 1:min(print_top, length(diffs))
        (ΔV, a_i, n_i, Δc, Δap, Δe, Vbf) = diffs[k]
        println("Top#$k: (a_i=$a_i, n_i=$n_i)  |ΔV|=$ΔV, |Δc|=$Δc, |Δap|=$Δap, |Δe|=$Δe")
    end

    return diffs
end

using Optim

function opt_retminus1_state(
    M::Float64, n::Float64, P::Float64, e_bar::Float64,
    EV_next_a::AbstractVector{Float64},
    parameters::NamedTuple;
)
    @unpack a_grid, a_min, a_max, β = parameters
    @unpack one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ = parameters
    @unpack c_floor, V_penalty = parameters

    # feasible upper bound for a'
    ap_hi = min(a_max, M - c_floor - (n > 0 ? e_bar : 0.0))
    if ap_hi < a_min
        return (V_penalty, c_floor, a_min, (n > 0 ? e_bar : 0.0), false)
    end

    # ---------- n = 0 : 1D over ap ----------
    if n == 0.0
        res_ap = optimize(ap -> begin
                Vp = lininterp(a_grid, EV_next_a, ap)
                if !isfinite(Vp) || Vp <= V_penalty
                    return 1e30
                end
                c = M - ap
                if !isfinite(c) || c <= c_floor
                    return 1e30
                end
                V = u_CRRA(c, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β*Vp
                -V
            end,
            a_min, ap_hi, Brent()
        )

        ap_star = Optim.minimizer(res_ap)
        V_star  = -Optim.minimum(res_ap)
        c_star  = M - ap_star
        return (V_star, c_star, ap_star, 0.0, true)
    end

    # ---------- n > 0 : nested 1D (ap outer, e inner) ----------
    scale = (n / P)^one_m_κ

    # inner maximize in e, given (m,Vp)
    inner_minimum(m, Vp) = begin
        e_hi = m - c_floor
        if e_hi < e_bar
            return 1e30
        end
        res_e = optimize(e -> begin
                c = m - e
                if !isfinite(c) || c <= c_floor || !isfinite(e) || e < e_bar
                    return 1e30
                end
                V = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
                             c_floor=c_floor, V_penalty=V_penalty) + β*Vp
                -V
            end,
            e_bar, e_hi, Brent()
        )
        Optim.minimum(res_e)  # negative value
    end

    res_ap = optimize(ap -> begin
            Vp = lininterp(a_grid, EV_next_a, ap)
            if !isfinite(Vp) || Vp <= V_penalty
                return 1e30
            end
            m = M - ap
            inner_minimum(m, Vp)
        end,
        a_min, ap_hi, Brent()
    )

    ap_star = Optim.minimizer(res_ap)

    # recover e*, c*, V*
    Vp = lininterp(a_grid, EV_next_a, ap_star)
    m  = M - ap_star
    e_hi = m - c_floor
    if e_hi < e_bar || !isfinite(Vp) || Vp <= V_penalty
        return (V_penalty, c_floor, ap_star, e_bar, false)
    end

    res_e = optimize(e -> begin
            c = m - e
            if !isfinite(c) || c <= c_floor || !isfinite(e) || e < e_bar
                return 1e30
            end
            V = u_CRRA_e(c, e, one_m_γ, inv_one_m_γ, one_m_κ, ψ_inv_one_m_κ, scale;
                         c_floor=c_floor, V_penalty=V_penalty) + β*Vp
            -V
        end,
        e_bar, e_hi, Brent()
    )

    e_star = Optim.minimizer(res_e)
    V_star = -Optim.minimum(res_e)
    c_star = m - e_star
    return (V_star, c_star, ap_star, e_star, true)
end

function compare_one_state_retminus1(variables, parameters; ϵ_i=4, ν_i=3, a_i=50, n_i=4)
    @unpack age_ret, age_min, aR_grid, w_grid, P_grid, q_bar_P_grid, n_grid = parameters

    age_i = age_ret - age_min

    # build EV at age_ret-1 (same as your solver)
    EV_tmp  = fill(parameters.V_penalty, parameters.a_size, parameters.n_size, parameters.ϵ_size, parameters.inf_size)
    cCE_tmp = fill(NaN, parameters.a_size, parameters.n_size, parameters.ϵ_size, parameters.inf_size)
    V_next = @view variables.V[:, 1, :, :, 2, age_i+1]
    c_next = @view variables.policy_c[:, 1, :, :, 2, age_i+1]
    fill_EV_Euc!(EV_tmp, cCE_tmp, V_next, c_next, parameters, age_i)

    w_bar = w_grid[ϵ_i, ν_i, age_i]
    M = aR_grid[a_i] + w_bar

    P     = P_grid[n_i, ϵ_i, ν_i, age_i]
    e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    n     = n_grid[n_i]
    EV_next_a = @view EV_tmp[:, n_i, ϵ_i, 2]

    (Vopt, copt, apopt, eopt, ok) = opt_retminus1_state(M, n, P, e_bar, EV_next_a, parameters)

    Vegm  = variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    cegm  = variables.policy_c[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    apegm = variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    eegm  = variables.policy_e[a_i, n_i, ϵ_i, ν_i, 2, age_i]

    println("ok = $ok")
    println("|ΔV|  = ", abs(Vegm - Vopt))
    println("|Δc|  = ", abs(cegm - copt))
    println("|Δap| = ", abs(apegm - apopt))
    println("|Δe|  = ", abs(eegm - eopt))
end

function compare_slice_retminus1_allan(variables, parameters; ϵ_i=4, ν_i=3)
    @unpack age_ret, age_min, aR_grid, w_grid, P_grid, q_bar_P_grid, n_grid = parameters
    @unpack a_size, n_size, V_penalty = parameters

    age_i = age_ret - age_min

    EV_tmp  = fill(V_penalty, a_size, n_size, parameters.ϵ_size, parameters.inf_size)
    cCE_tmp = fill(NaN, a_size, n_size, parameters.ϵ_size, parameters.inf_size)
    V_next  = @view variables.V[:, 1, :, :, 2, age_i+1]
    c_next  = @view variables.policy_c[:, 1, :, :, 2, age_i+1]
    fill_EV_Euc!(EV_tmp, cCE_tmp, V_next, c_next, parameters, age_i)

    w_bar = w_grid[ϵ_i, ν_i, age_i]
    maxdV = 0.0; maxdc = 0.0; maxdap = 0.0; maxde = 0.0
    argmax = (0,0)

    for n_i in 1:n_size, a_i in 1:a_size
        Vegm = variables.V[a_i, n_i, ϵ_i, ν_i, 2, age_i]
        if !isfinite(Vegm) || Vegm <= V_penalty
            continue
        end

        M     = aR_grid[a_i] + w_bar
        P     = P_grid[n_i, ϵ_i, ν_i, age_i]
        e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
        n     = n_grid[n_i]
        EV_next_a = @view EV_tmp[:, n_i, ϵ_i, 2]

        (Vopt, copt, apopt, eopt, ok) = opt_retminus1_state(M, n, P, e_bar, EV_next_a, parameters)
        ok || continue

        cegm  = variables.policy_c[a_i, n_i, ϵ_i, ν_i, 2, age_i]
        apegm = variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i]
        eegm  = variables.policy_e[a_i, n_i, ϵ_i, ν_i, 2, age_i]

        dV  = abs(Vegm - Vopt)
        dc  = abs(cegm - copt)
        dap = abs(apegm - apopt)
        de  = abs(eegm - eopt)

        if dV > maxdV
            maxdV = dV; maxdc = dc; maxdap = dap; maxde = de
            argmax = (a_i, n_i)
        end
    end

    println("max |ΔV|  = $maxdV at (a_i,n_i)=$argmax")
    println("corr |Δc|  = $maxdc")
    println("corr |Δap| = $maxdap")
    println("corr |Δe|  = $maxde")
end

function debug_retminus1_n0(variables, parameters; ϵ_i=4, ν_i=3, a_i=22)
    @unpack age_ret, age_min, aR_grid, w_grid, P_grid, q_bar_P_grid, n_grid = parameters
    @unpack a_size, n_size, ϵ_size, inf_size, V_penalty = parameters
    @unpack a_grid, β, one_m_γ, inv_one_m_γ, c_floor = parameters

    age_i = age_ret - age_min

    EV_tmp  = fill(V_penalty, a_size, n_size, ϵ_size, inf_size)
    cCE_tmp = fill(NaN,       a_size, n_size, ϵ_size, inf_size)
    V_next  = @view variables.V[:, 1, :, :, 2, age_i+1]
    c_next  = @view variables.policy_c[:, 1, :, :, 2, age_i+1]
    fill_EV_Euc!(EV_tmp, cCE_tmp, V_next, c_next, parameters, age_i)

    n_i = 1
    n   = n_grid[n_i]  # = 0.0
    w_bar = w_grid[ϵ_i, ν_i, age_i]
    M = aR_grid[a_i] + w_bar

    EV_next_a = @view EV_tmp[:, n_i, ϵ_i, 2]

    # EGM policy
    ap_egm = variables.policy_a_p[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    c_egm  = variables.policy_c[a_i, n_i, ϵ_i, ν_i, 2, age_i]
    Vp_egm = lininterp(a_grid, EV_next_a, ap_egm)
    V_egm_eval = u_CRRA(c_egm, one_m_γ, inv_one_m_γ; c_floor=c_floor, V_penalty=V_penalty) + β*Vp_egm

    # Optim policy
    P     = P_grid[n_i, ϵ_i, ν_i, age_i]
    e_bar = q_bar_P_grid[n_i, ϵ_i, ν_i, age_i]
    Vopt, copt, apopt, eopt, ok = opt_retminus1_state(M, n, P, e_bar, EV_next_a, parameters)

    println("EGM : ap=$ap_egm, c=$c_egm, V(eval)=$V_egm_eval")
    println("OPT : ap=$apopt, c=$copt, V      =$Vopt, ok=$ok")
    println("|ΔV(eval)|=", abs(V_egm_eval - Vopt), "  |Δap|=", abs(ap_egm-apopt), "  |Δc|=", abs(c_egm-copt))
end
