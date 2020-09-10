module FourBarLinkages

# export compute_generic_solutions, save_generic_solutions, load_generic_solutions,
#     solve_instance, four_bars, FourBar, animate, configurations
export start_demo


import HomotopyContinuation
const HC = HomotopyContinuation
using LinearAlgebra
# using HomotopyContinuation, LinearAlgebra, DynamicPolynomials
# using Interpolations
# using StaticArrays
# using Makie
using Parameters: @unpack
using AbstractPlotting, GLMakie
using AbstractPlotting.MakieLayout
using ColorSchemes, Printf

# import JLD2, FileIO
#
# include("sp_homotopy.jl")
#

# This follows the description in
#  Wampler, C. W., Morgan, A. P., and Sommese, A. J. (March 1, 1992).
#  "Complete Solution of the Nine-Point Path Synthesis Problem for Four-Bar Linkages."
#  ASME. J. Mech. Des. March 1992; 114(1): 153–159. https://doi.org/10.1115/1.2916909
struct FourBar
    x::ComplexF64
    a::ComplexF64
    y::ComplexF64
    b::ComplexF64
    P₀::ComplexF64 # coupler point
end
function Base.getproperty(F::FourBar, field::Symbol)
    if field === :A
        getfield(F, :P₀) + getfield(F, :a)
    elseif field === :B
        getfield(F, :P₀) + getfield(F, :b)
    elseif field === :C
        getfield(F, :P₀) + getfield(F, :y)
    elseif field === :D
        getfield(F, :P₀) + getfield(F, :x)
    elseif field === :u
        getfield(F, :x) - getfield(F, :a)
    elseif field === :v
        getfield(F, :y) - getfield(F, :b)
    else
        getfield(F, field)
    end
end

function FourBar(;
    A::ComplexF64,
    B::ComplexF64,
    C::ComplexF64,
    D::ComplexF64,
    P₀ = (D + C) / 2,
)
    x = D - P₀
    a = A - P₀
    y = C - P₀
    b = B - P₀
    FourBar(x, a, y, b, P₀)
end

function FourBar(F::FourBar, angles)
    FourBar(
        A = F.A,
        B = F.B,
        C = (F.B + F.v * angles.μ),
        D = (F.A + F.u * angles.λ),
        P₀ = (F.A + F.u * angles.λ - F.x * angles.θ),
    )
end

function four_bar_from_lengths(; A::ComplexF64, B::ComplexF64, BC, CD, AD)
    HC.@var a[1:2] b[1:2] c[1:2] d[1:2] l[1:5]
    F = HC.System(
        [
            sum((b - c) .^ 2) - BC^2,
            sum((c - d) .^ 2) - CD^2,
            sum((a - d) .^ 2) - AD^2,
            sum(l .* [c; d; 1]),
        ],
        [c; d],
        [l; a; b],
    )

    for i = 1:5
        res =
            HC.solve(F, target_parameters = [randn(5); real(A); imag(A); real(B); imag(B)])
        S = HC.real_solutions(res)
        if !isempty(S)
            s = S[1]
            C = complex(s[1], s[2])
            D = complex(s[3], s[4])
            return FourBar(A = A, B = B, C = C, D = D)
        end
    end
    error("Could not find a valid configuartion.")
end

# function FourBar(solution::Vector{ComplexF64}, P₀::ComplexF64)
#     x, a, y, b = solution
#     FourBar(x, a, y, b, P₀)
# end
#
Base.show(io::IO, F::FourBar) = print_fieldnames(io, F)
#
"""
     print_fieldnames(io::IO, obj)

 A better default printing for structs.
 """
function print_fieldnames(io::IO, obj)
    println(io, typeof(obj), ":")
    for name in fieldnames(typeof(obj))
        if getfield(obj, name) !== nothing
            println(io, " • ", name, " → ", getfield(obj, name))
        end
    end
end


is_conjugated_pair(u, v, tol) = abs(u - conj(v)) < tol
is_valid_loop_solution(μ, θ, τ, τ̂; tol::Float64 = 1e-12) =
    is_conjugated_pair(τ, τ̂, tol) && abs(abs(μ) - 1) < tol && abs(abs(θ) - 1) < tol

function loop_homotopy()
    params = HC.@var x a y b x̂ â ŷ b̂
    HC.@var λ μ θ τ τ̂
    HC.Homotopy(
        [
            (x - a) * λ - x * θ - τ + a,
            (x̂ - â) * θ - x̂ * λ - (τ̂ - â) * λ * θ,
            (y - b) * μ - y * θ - τ + b,
            (ŷ - b̂) * θ - ŷ * μ - (τ̂ - b̂) * μ * θ,
        ],
        [μ, θ, τ, τ̂],
        λ,
        collect(params),
    )
end

#
# function loop_equations(F::FourBar)
#     @unpack x, a, y, b = F
#     x̂, â, ŷ, b̂ = conj.((x, a, y, b))
#     vars = @var λ μ θ τ τ̂ # and offsets
#     System([
#         (x - a) * λ - x * θ - τ + a,
#         (x̂ - â) * θ - x̂ * λ - (τ̂ - â) * λ * θ,
#         (y - b) * μ - y * θ - τ + b,
#         (ŷ - b̂) * θ - ŷ * μ - (τ̂ - b̂) * μ * θ,
#     ], vars
# end

function trace_curve(F::FourBar; Δt = 1e-2, max_steps = 1_000, accuracy = 1e-8)
    H = loop_homotopy()
    params = let
        @unpack x, a, y, b = F
        x̂, â, ŷ, b̂ = conj.((x, a, y, b))
        [x, a, y, b, x̂, â, ŷ, b̂]
    end
    # compute once solutions for our loop homotopy at λ₀
    λ_gen = randn(ComplexF64)
    F = HC.System(H.expressions, H.variables, [H.t; H.parameters])
    generic_solutions = HC.solutions(HC.solve(F, target_parameters = [λ_gen; params]))


    # a valid solution is
    φ = φ₀ = 0.0
    μ, θ, τ, τ̂ = [1.0, 1, 0, 0im]

    angles = [(λ = cis(φ), μ = μ, θ = θ)]
    # loop, (λ, _) = loop_equations(F)
    # φ = φ₀ = angle(λ₀)
    # angles = [(cis(φ₀), x₀[1], x₀[2])]
    # μ₀, θ₀ = angle(x₀[1]), angle(x₀[2])
    tracker = HC.Tracker(HC.fix_parameters(H, params); options = (max_steps = 200,))

    # tracker = coretracker(SPHomotopy(loop, λ), [randn(ComplexF64, 4)]; accuracy = accuracy)
    HC.init!(tracker, [1.0, 1, 0, 0im], cis(φ), cis(φ + Δt))
    x = tracker.state.x
    y = copy(x)
    Δ = Δt
    for i = 2:max_steps
        if φ * (φ + Δ) < 0 # jump over zero
            φ′ = 0.0
        else
            φ′ = φ + Δ
        end
        retcode = HC.track!(tracker, y, cis(φ), cis(φ′))
        if retcode != HC.TrackerCode.success
            @warn "PathTracker failed with $retcode"
            break
        end
        μ, θ, τ, τ̂ = x
        if is_valid_loop_solution(μ, θ, τ, τ̂)
            if abs(angle(angles[end].μ) - angle(μ)) < 2Δt
                φ = φ′
                Δ = clamp(2Δ, -Δt, Δt)
                push!(angles, (λ = cis(φ), μ = μ, θ = θ))
                y .= x
                # reduce step size since μ jump is too large
            else
                Δ /= 2
                continue
            end
        else
            # jump to different branch
            branch_solutions = map(generic_solutions) do s
                HC.solution(HC.track(tracker, s, λ_gen, cis(φ)))
            end
            y .= branch_solutions[last(findmax(norm.(branch_solutions .- Ref(y))))]
            Δ = -Δ
            μ, θ = y[1], y[2]
            push!(angles, (λ = cis(φ), μ = μ, θ = θ))
        end
        if φ > π
            φ -= 2π
        elseif φ < -π
            φ += 2π
        end
        if abs(φ - φ₀) < 0.5Δt && abs(1 - y[1]) < 1e-3 && abs(1 - y[2]) < 1e-3
            pop!(angles)
            break
        end
    end
    angles
end



fourbar_to_points(F) = Point.(reim.((F.A, F.B, F.C, F.D)))


function energy(
    fourbar;
    resting_lengths::Vector{Float64},
    elasticities::Vector{Float64},
    E::ComplexF64,
    F::ComplexF64,
)
    c1, c2 = elasticities
    r1, r2 = resting_lengths
    c1 * (sqrt(abs2(fourbar.C - E)) - r1)^2 + c2 * (sqrt(abs2(fourbar.D - F)) - r2)^2
end

function energy_system(;
    bar_lengths::Vector,
    resting_lengths,
    elasticities,
    A::ComplexF64,
    B::ComplexF64,
)
    HC.@var p[1:2, 1:6]

    bars = [(2, 3), (3, 4), (4, 1)]
    cables = [(3, 5), (2, 6)]

    HC.@var l²[1:3] δ[1:2] λ[1:5]

    G = [
        [sum((p[:, i] - p[:, j]) .^ 2) .- l²[k] for (k, (i, j)) in enumerate(bars)]
        [sum((p[:, i] - p[:, j]) .^ 2) .- δ[k]^2 for (k, (i, j)) in enumerate(cables)]
    ]
    Q = sum((δ .- resting_lengths) .^ 2)
    L = Q + λ'G

    X = [p[:, 3]; p[:, 4]]
    # fix constants
    L′ = HC.subs(
        L,
        p[:, 1] => [real(A), imag(A)],
        p[:, 2] => [real(B), imag(B)],
        l² => bar_lengths .^ 2,
    )
    Y = [p[:, 5]; p[:, 6]]
    ∇L = HC.differentiate(L′, [X; δ; λ])
    F = HC.System(∇L, variables = [X; δ; λ], parameters = Y)
end

function local_minimum_checker(;
    bar_lengths::Vector,
    resting_lengths,
    elasticities,
    A::ComplexF64,
    B::ComplexF64,
)
    HC.@var p[1:2, 1:6]


    bars = [(2, 3), (3, 4), (4, 1)]
    cables = [(3, 5), (2, 6)]

    HC.@var l²[1:3] r[1:2] δ[1:2] λ[1:5]

    G = HC.subs(
        [
            [sum((p[:, i] - p[:, j]) .^ 2) .- l²[k] for (k, (i, j)) in enumerate(bars)]
            [sum((p[:, i] - p[:, j]) .^ 2) .- δ[k]^2 for (k, (i, j)) in enumerate(cables)]
        ],
        p[:, 1] => [real(A), imag(A)],
        p[:, 2] => [real(B), imag(B)],
        l² => bar_lengths .^ 2,
    )

    Q = sum((δ .- resting_lengths) .^ 2)
    L = Q + λ'G

    X = [p[:, 3]; p[:, 4]]
    # fix constants
    L′ = HC.subs(
        L,
        p[:, 1] => [real(A), imag(A)],
        p[:, 2] => [real(B), imag(B)],
        l² => bar_lengths .^ 2,
    )
    Y = [p[:, 5]; p[:, 6]]
    ∇L = HC.differentiate(L′, [X; δ])
    HL = HC.CompiledSystem(HC.System(∇L, [X; δ], [λ; Y]))
    dG = HC.CompiledSystem(HC.System(G, [X; δ], [λ; Y]))

    (s, params) -> begin
        v = s[1:6]
        q = [s[7:end]; params]
        W = HC.jacobian!(zeros(6, 6), HL, v, q)
        V = nullspace(HC.jacobian!(zeros(5, 6), dG, v, q))
        all(e -> e ≥ 1e-14, eigvals(V' * W * V))
    end
end




function comp_min_energy_positions(;
    A::ComplexF64,
    B::ComplexF64,
    E::ComplexF64,
    F::ComplexF64,
    bar_lengths,
    elasticities,
    resting_lengths,
)

    energy_sys = energy_system(
        A = A,
        B = B,
        bar_lengths = bar_lengths,
        elasticities = elasticities,
        resting_lengths = resting_lengths,
    )
    gen_params = randn(ComplexF64, 4)
    gen_res = HC.solve(energy_sys; target_parameters = gen_params)


    min_checker = local_minimum_checker(
        A = A,
        B = B,
        bar_lengths = bar_lengths,
        elasticities = elasticities,
        resting_lengths = resting_lengths,
    )
    res = HC.solve(
        energy_sys,
        gen_res;
        start_parameters = gen_params,
        target_parameters = [reim(F)..., reim(E)...],
    )
    min_energy_sols = filter(
        s -> s[5] > 0 && s[6] > 0 && min_checker(s, [reim(F)..., reim(E)...]),
        HC.real_solutions(res),
    )
    min_energy_positions = map(min_energy_sols) do s
        FourBar(; A = A, B = B, C = complex(s[1], s[2]), D = complex(s[3], s[4]))
    end
end

# coordinate transformation between two rects
function transferrects(pos, rectfrom, rectto)
    fracpos = (pos .- rectfrom.origin) ./ rectfrom.widths
    fracpos .* rectto.widths .+ rectto.origin
end

function add_control_node!(ax, n; kwargs...)
    selected = Ref(false)
    plt = scatter!(ax, n; kwargs...)
    lift(events(ax.scene).mouseposition) do pos
        x, y = transferrects(pos, ax.scene.px_area[], ax.limits[])
        if AbstractPlotting.is_mouseinside(ax.scene) && selected[]
            p = Point2f0(x, y)
            n[] = p
        end
        nothing
    end
    on(events(ax.scene).mousedrag) do drag
        if ispressed(ax.scene, Mouse.left)
            if drag == Mouse.down
                plot, _idx = mouse_selection(ax.scene)
                if plot == plt
                    selected[] = true
                end
            end
        elseif drag == Mouse.up || !AbstractPlotting.is_mouseinside(ax.scene)
            selected[] = false
        end
    end
    n
end



function start_demo(;
A = -1.0 + 0im,
    B = 1.0 + 0im,
    E = 5.0 + 0im,
    F = -4.0 - 4im,
    bar_lengths = [3.0, 2.0, 1.5],
    elasticities = [1.0, 4.0],
    resting_lengths = [0.1, 0.1])

    scene, layout = layoutscene(30, resolution = ((1200, 1200)))
    ax = (layout[1:3, 1] = LAxis(scene, title = "Control plane"))
    energy_ax = (layout[4, 1] = LAxis(scene, title = "Energy"))
    energy_ax.ytickformat = xs -> [@sprintf("%05.1f", x) for x in xs]

    min_energy_positions = comp_min_energy_positions(;
        A = A,
        B = B,
        E = E,
        F = F,
        bar_lengths = bar_lengths,
        elasticities = elasticities,
        resting_lengths = resting_lengths,
    )

    fourbar = min_energy_positions[1]
    angles = trace_curve(fourbar; Δt = 1e-1, max_steps = 10_000)
    loop = FourBar.(Ref(fourbar), angles)


    loop_index = Node(1)
    NABCD = lift(k -> fourbar_to_points(loop[k]), loop_index)
    NF = Node(Point(reim(F)))

    N_energy = lift(NF) do f
        energy.(
            loop;
            E = E,
            F = complex(f[1], f[2]),
            elasticities = elasticities,
            resting_lengths = resting_lengths,
        )
    end

    on(N_energy) do energy
        n = length(energy)
        k = loop_index[]
        # let's go to the right
        # println(mod1(k - 1, n), " ", k, " ", mod1(k + 1, n))
        # println(energy[mod1(k - 1, n)], " ", energy[k], " ", energy[mod1(k + 1, n)])
        while energy[mod1(k + 1, n)] < energy[k] * (1 + 1e-8)
            k = mod1(k + 1, n)
        end
        # @show energy[mod1(k + 1, n)] < energy[k]
        # if we didn't do a step, go to the left
        if k == loop_index[]
            while energy[mod1(k - 1, n)] < energy[k] * (1 + 1e-8)
                k = mod1(k - 1, n)
            end
        end
        loop_index[] = k
    end
    # scene = Scene(limits = limits, resolution = (1200, 1200), scale_plot = false)
    linesegments!(ax, lift(NABCD) do (A, B, C, D)
        [B => C, C => D, D => A]
    end, linewidth = 2, color = :black)
    linesegments!(ax, lift(NABCD, NF) do (A, B, C, D), F
        [D => F]
    end, color = :teal)
    linesegments!(ax, lift(NABCD) do (A, B, C, D)
        [C => Point(reim(E))]
    end, color = :teal)

    add_control_node!(ax, NF; markersize = 20px, marker = :cross, color = :black)
    scatter!(ax, lift(collect, NABCD), color = :black)
    scatter!(ax, Point(reim(E)), marker = :rect, color = :black)


    # scatter!(ax, Point(reim(F)), marker=:cross)
    lines!(
        ax,
        map(fb -> Point(reim(fb.P₀)), loop),
        linestyle = :dash,
    )


    NABCD[] = fourbar_to_points(first(loop))

    xlims!(ax, [-10, 10]) # as vector
    ylims!(ax, [-10, 10]) # as vector
    ax.aspect = AxisAspect(1)

    scatter!(energy_ax, 1:length(N_energy[]), N_energy)
    scatter!(energy_ax, lift((i, E) -> Point(i, E[i]), loop_index, N_energy), color = :red)
    autolimits!(energy_ax)
    on(N_energy) do _
        autolimits!(energy_ax)
    end
    display(scene)
end

# function δ_angles_pairs(F::FourBar, δ)
#     loop, (_, _, _, τ, τ̂) = loop_equations(F)
#     pairs = Tuple{ComplexF64,NTuple{3,ComplexF64}}[]
#     for δᵢ in δ
#         # do only 1:3 since overdetermined sometimes doesn't work
#         # I think the randomization is then "bad".
#         # Which seems to happen for these equations quite often.
#         sols = solutions(solve([subs(f, τ => δᵢ, τ̂ => conj(δᵢ)) for f in loop][1:3]))
#         # filter out non physical solutions
#         for (λ, μ, θ) in sols
#             if isapprox(abs(λ), 1; atol = 1e-6) &&
#                isapprox(abs(μ), 1; atol = 1e-6) &&
#                isapprox(abs(θ), 1; atol = 1e-6)
#                 push!(pairs, (δᵢ, (λ, μ, θ)))
#             end
#         end
#     end
#     pairs
# end

# function equations()
#     @polyvar x a y b x̂ â ŷ b̂
#     @polyvar γ[1:8] γ̂[1:8] δ[1:8] δ̂[1:8]
#     #system of polynomials
#     D1 = [(â*x-δ̂[i]*x)*γ[i]+(a*x̂-δ[i]*x̂)*γ̂[i]+(â-x̂)*δ[i]+(a-x)*δ̂[i]-δ[i]*δ̂[i] for i in 1:8]
#     D2 = [(b̂*y-δ̂[i]*y)*γ[i]+(b*ŷ-δ[i]*ŷ)*γ̂[i]+(b̂-ŷ)*δ[i]+(b-y)*δ̂[i]-δ[i]*δ̂[i] for i in 1:8]
#     D3 = [γ[i]*γ̂[i]+γ[i]+γ̂[i] for i in 1:8]
#     (F=[D1;D2;D3], xayb=[x, a, y, b, x̂, â, ŷ, b̂], γ=γ, γ̂=γ̂, δ=δ, δ̂=δ̂)
# end
#
# function compute_start_pair()
#     eqs = equations()
#     α = rand(8)
#     Γ, Γ̂ = cis.(α.*2π) .- 1, cis.(α.*(-2).*π) .- 1
#     xayb₀ = randn(ComplexF64,8)
#
#     start_help = map(1:8) do i
#         Dᵢ = [eqs.F[i], eqs.F[i+8]]
#         Gᵢ = [subs(f, eqs.xayb => xayb₀, eqs.γ=>Γ,eqs.γ̂=>Γ̂) for f in Dᵢ];
#         first(solutions(solve(Gᵢ)))
#     end
#
#     δ₀, δ̂₀ = first.(start_help), last.(start_help)
#     [xayb₀;Γ;Γ̂], [δ₀; δ̂₀]
# end
#
# function compute_γ_γ̂(x, a, y, b, δ, x̂, â, ŷ, b̂, δ̂)
#     γ = @MVector zeros(eltype(x), 8)
#     γ̂ = @MVector zeros(eltype(x), 8)
#     for j in 1:8
#         Aⱼ = @SMatrix [(â-δ̂[j])*x (a-δ[j])*x̂; (b̂-δ̂[j])*y (b-δ[j])*ŷ]
#         cⱼ = -@SVector [δ[j]*(â-x̂)+δ̂[j]*(a-x)-δ[j]*δ̂[j],
#                        δ[j]*(b̂-ŷ)+δ̂[j]*(b-y)-δ[j]*δ̂[j]]
#         γⱼ, γ̂ⱼ = Aⱼ \ cⱼ
#         γ[j] = γⱼ
#         γ̂[j] = γ̂ⱼ
#     end
#     SVector(γ), SVector(γ̂)
# end
#
# r(x,a,y,b) = ((x-a)*y/(x-y), (b*x - a*y)/(x-y), a-x, a)
# function roberts_cognates(v, δ, δ̂)
#     x, a, y, b, x̂, â, ŷ, b̂ = v
#     xayb₁ = r(x, a, y, b)
#     x̂âŷb̂₁ = r(x̂, â, ŷ, b̂)
#     γ₁, γ̂₁ = compute_γ_γ̂(xayb₁..., δ, x̂âŷb̂₁..., δ̂)
#     xayb₂ = r(xayb₁...)
#     x̂âŷb̂₂ = r(x̂âŷb̂₁...)
#     γ₂, γ̂₂ = compute_γ_γ̂(xayb₂..., δ, x̂âŷb̂₂..., δ̂)
#     v₁ = SVector(xayb₁..., x̂âŷb̂₁..., γ₁..., γ̂₁...)
#     v₂ = SVector(xayb₂..., x̂âŷb̂₂..., γ₂..., γ̂₂...)
#     v₁, v₂
# end
# # switch roles of (x,a) and (y,b)
# # We have the variable ordering x a y b x̂ â ŷ b̂
# # So we need to switch to y b x a ŷ b̂ x̂ â
# symmetry(s) = [[s[3],s[4],s[1],s[2],s[7],s[8],s[5],s[6]]; s[9:24]]
#
# compute_generic_solutions(;kwargs...) = compute_generic_solutions(compute_start_pair()...;kwargs...)
#
# jld2_file(filename) = endswith(filename, ".jld2") ? filename : filename * ".jld2"
# function compute_generic_solutions(x₀, p₀; filename=nothing)
#     eqs = equations()
#     δ₀, δ̂₀ = p₀[1:8], p₀[9:16]
#     group_actions = GroupActions(symmetry, s -> roberts_cognates(s, δ₀, δ̂₀))
#     @info "Computing all 1442 generic solutions..."
#     result = monodromy_solve(eqs.F, x₀, p₀;
#                 parameters=[eqs.δ;eqs.δ̂],
#                 group_actions=group_actions,
#                 target_solutions_count=1442,
#                 equivalence_classes=true)
#     δ₀, δ̂₀ = result.parameters[1:8], result.parameters[9:16]
#     data = Dict(["δ₀" => δ₀, "δ̂₀" => δ̂₀, "solutions" => reduce(hcat, result.solutions)])
#     if filename !== nothing
#         FileIO.save(jld2_file(filename), data)
#     end
#     data
# end
#
# function save_generic_solutions(result::MonodromyResult, filename::String)
#     δ₀, δ̂₀ = result.parameters[1:8], result.parameters[9:16]
#     FileIO.save(endswith(filename, ".jld2") ? filename : filename * ".jld2",
#                 "δ₀", δ₀, "δ̂₀", δ̂₀,
#                 "solutions", reduce(hcat, result.solutions))
# end
#
# const GENERIC_SOLUTIONS_PATH = joinpath(@__DIR__, "..", "data", "four_bar_start_solutions.jld2")
#
# function load_generic_solutions(filename=GENERIC_SOLUTIONS_PATH)
#     data = FileIO.load(endswith(filename, ".jld2") ? filename : filename * ".jld2")
#     (solutions=data["solutions"], δ₀=data["δ₀"], δ̂₀=data["δ̂₀"])
# end
#
# function solve_instance(δ::Vector{ComplexF64}, filename::String=GENERIC_SOLUTIONS_PATH; kwargs...)
#     eqs = equations()
#     δ̂ = conj.(δ)
#     generic_sols, δ₀, δ̂₀ = load_generic_solutions(filename)
#     start_sols = [view(generic_sols, 1:24, i) for i in 1:size(generic_sols,2)]
#     res = solve(eqs.F, start_sols;
#             parameters = [eqs.δ; eqs.δ̂],
#             start_parameters=[δ₀;δ̂₀],
#             target_parameters=[δ;δ̂], kwargs...)
# end
#
# is_conjugated_pair(u, v, tol) = abs(u - conj(v)) < tol
# is_physical_solution(s, tol) = all(j -> is_conjugated_pair(s[j], s[j+4], tol), 1:4)
# function physical_four_bars(solutions; tol=1e-10)
#     filter(s -> is_physical_solution(s, tol), solutions)
# end
#
# function four_bars(points::Vector{<:Complex}, filename::String=GENERIC_SOLUTIONS_PATH; kwargs...)
#     @assert length(points) == 9 "Expected 9 points"
#     P₀ = points[1]
#     δ = points[2:9] .- P₀
#     result = solve_instance(δ, filename)
#     four_bars(result, points; kwargs...)
# end
# function four_bars(result, points; real_tol=1e-10)
#     P₀ = points[1]
#     fourbars = FourBar[]
#     for s in solutions(result; only_nonsingular=true)
#         if is_physical_solution(s, real_tol)
#             push!(fourbars, FourBar(s, P₀))
#         end
#     end
#     fourbars
# end
#
# function is_valid_loop_solution(r)
#     is_conjugated_pair(r[3], r[4], 1e-4) || return false
#     abs(abs(r[1])-1) < 1e-4 && abs(abs(r[2]) - 1) < 1e-4
# end
#
#
# function loop_equations(F::FourBar)
#     @unpack x, a, y, b = F
#     x̂, â, ŷ, b̂ = conj.((x,a,y,b))
#     vars = @polyvar λ μ θ τ τ̂ # and offsets
#     [(x-a)*λ-x*θ-τ+a,
#      (x̂-â)*θ-x̂*λ-(τ̂-â)*λ*θ,
#      (y-b)*μ-y*θ-τ+b,
#      (ŷ-b̂)*θ-ŷ*μ-(τ̂-b̂)*μ*θ], vars
# end
#
# function trace_points(F::FourBar, λ₀, x₀; Δt = 1e-3, max_steps=20_000, accuracy=1e-8)
#     loop, (λ, _) = loop_equations(F)
#     φ = φ₀ = angle(λ₀)
#     angles = [(cis(φ₀), x₀[1], x₀[2])]
#     μ₀, θ₀ = angle(x₀[1]), angle(x₀[2])
#     tracker = coretracker(SPHomotopy(loop, λ), [randn(ComplexF64, 4)]; accuracy=accuracy)
#     HC.init!(tracker, x₀, cis(φ), cis(φ+Δt))
#     x = current_x(tracker)
#     y = copy(x)
#     for i in 2:max_steps
#         retcode = track!(tracker, y, cis(φ), cis(φ+Δt))
#
#         if retcode != HC.CoreTrackerStatus.success
#             @warn "PathTracker failed with $retcode"
#             break
#         end
#
#         if is_valid_loop_solution(x)
#             φ += Δt
#             push!(angles, (cis(φ), x[1], x[2]))
#             y .= current_x(tracker)
#         else
#             # jump to different branch
#             branch_solutions = solutions(solve([subs(f, λ => cis(φ)) for f in loop]))
#             y .= branch_solutions[last(findmax(norm.([s - y for s in branch_solutions])))]
#             Δt = -Δt
#             push!(angles, (cis(φ), y[1], y[2]))
#         end
#         if φ > π
#             φ -= 2π
#         elseif φ < -π
#             φ += 2π
#         end
#         if abs(φ - φ₀) < 0.5Δt &&
#            abs(x₀[1] - y[1]) < 1e-3 &&
#            abs(x₀[2] - y[2]) < 1e-3
#                 break
#         end
#     end
#     angles
# end
#
# function δ_angles_pairs(F::FourBar, δ)
#     loop, (_, _, _, τ, τ̂) = loop_equations(F)
#     pairs = Tuple{ComplexF64, NTuple{3, ComplexF64}}[]
#     for δᵢ in δ
#         # do only 1:3 since overdetermined sometimes doesn't work
#         # I think the randomization is then "bad".
#         # Which seems to happen for these equations quite often.
#         sols = solutions(solve([subs(f, τ => δᵢ, τ̂ => conj(δᵢ)) for f in loop][1:3]))
#         # filter out non physical solutions
#         for (λ, μ, θ) in sols
#             if isapprox(abs(λ), 1; atol=1e-6) &&
#                isapprox(abs(μ), 1; atol=1e-6) &&
#                isapprox(abs(θ), 1; atol=1e-6)
#                 push!(pairs, (δᵢ, (λ, μ, θ)))
#             end
#         end
#     end
#     pairs
# end
#
# function missed_coupler_points(δ_angles_pairs, angles)
#     filter(δ_angles_pairs) do (δ, (λ, μ, θ))
#         α = SVector(λ, μ, θ)
#         for s in angles
#             if norm(α - SVector(s)) < 1e-1
#                 return false
#             end
#         end
#         true
#     end
# end
#
# function configurations(F::FourBar, coupler_points::Vector{ComplexF64}; kwargs...)
#     δ = coupler_points[2:9] .- coupler_points[1]
#     pairs = δ_angles_pairs(F, δ)
#     # by the computation is (λ, μ, θ) = (0,0,0) a valid angle configuration
#     angles₁ = trace_points(F, 0.0, [1.0, 1.0, 0.0, 0im]; kwargs...)
#     curves = [angles₁]
#     pairs₂ = missed_coupler_points(pairs, angles₁)
#     if !isempty(pairs₂)
#         δ, (λ, μ, θ) = pairs₂[1]
#         angles₂ = trace_points(F, λ, [μ, θ, δ, conj(δ)]; kwargs...)
#         push!(curves, angles₂)
#     end
#     curves
# end
#
#
# to_point(z::ComplexF64) = Point2f0(reim(z)...)
#
# function four_bar_positions(F::FourBar, angles)
#     A = to_point(F.A)
#     B = to_point(F.B)
#     map(angles) do ((λ, μ, θ))
#         (A=A, B=B, C=to_point(F.A + F.u * λ), D=to_point(F.B + F.v * μ),
#             P=to_point(F.A + F.u * λ - F.x * θ))
#     end
# end
#
# function compute_limits(positions)
#     xmin = ymin = Float32(Inf)
#     xmax = ymax = -Float32(Inf)
#     for pos in positions
#         xmin = min(xmin, pos.A[1], pos.B[1], pos.C[1], pos.D[1], pos.P[1])
#         xmax = max(xmax, pos.A[1], pos.B[1], pos.C[1], pos.D[1], pos.P[1])
#         ymin = min(ymin, pos.A[2], pos.B[2], pos.C[2], pos.D[2], pos.P[2])
#         ymax = max(ymax, pos.A[2], pos.B[2], pos.C[2], pos.D[2], pos.P[2])
#     end
#     w = (xmax - xmin)
#     h = (ymax - ymin)
#     FRect(xmin - 0.1w, ymin-0.1h, 1.1w, 1.1h)
# end
#
# function fix_scene_limits!(scene, limits)
#     O = Point2f0(limits.origin...)
#     w, h = limits.widths
#     O1 = O + Point2f0(0,w)
#     O2 = O1 + Point2f0(h,0)
#     O3 = O + Point2f0(h,0)
#     Makie.lines!(scene, [O,O1,O2,O3, O], color = :transparent, show_axis=false);
# end
#
# function static!(scene, F::FourBar, coupler_points; markersize=1.0)
#     Makie.scatter!(scene, [to_point(F.A), to_point(F.B)], color=:BLACK,
#                 markersize=markersize, marker='▲', show_axis=false)
#     if coupler_points !== nothing
#         Makie.scatter!(scene, to_point.(coupler_points), marker=:x, markersize=markersize,
#             color=:INDIANRED, show_axis=false)
#     end
#     scene
# end
#
#
# function animate(F::FourBar, coupler_points; Δt=1e-3, kwargs...)
#     angles = configurations(F, coupler_points; Δt=Δt)
#     Makie.animate(F, angles, coupler_points; kwargs...)
# end
# function animate(F::FourBar, angles::Vector{<:Vector}, coupler_points::Vector{ComplexF64}; kwargs...)
#     Makie.animate(F, angles..., coupler_points; kwargs...)
# end
# function animate(F::FourBar,
#         angles::Vector{NTuple{3,ComplexF64}},
#         coupler_points::Union{Nothing,Vector{ComplexF64}}=nothing;
#         show_axis=true,
#         fps=24, seconds=5,
#         color=:CADETBLUE,
#         color2=color,
#         loop::Union{Int,Bool} = false,
#         filename::Union{String,Nothing}=nothing)
#
#     positions = four_bar_positions(F, angles)
#     P = map(pos -> pos.P, positions)
#     #convert vectors to points
#     #plot Fourbar
#     limits = compute_limits(positions);
#     markersize = max(limits.widths...)/50
#     scene = Scene(limits=limits,  );
#     fix_scene_limits!(scene, limits)
#
#     static!(scene, F, coupler_points; markersize=markersize)
#
#
#     source, loop_closed_ref = add_mechanism!(scene, positions;
#             color=color,
#             show_axis=show_axis, markersize=markersize)
#     itp = interpolate_curve(positions)
#     N = seconds * fps
#     if filename !== nothing
#         record(scene, filename, 0:N; framerate=fps) do k
#             push!(source, round(Int, itp(k/N)))
#         end
#     else
#         display(scene)
#
#         if loop === false
#             nloops = 1
#         elseif loop === true
#             nloops = typemax(Int64)
#         else
#             nloops = loop::Int
#         end
#
#         for _ in 1:nloops
#             for k in 0:N
#                 push!(source, round(Int, itp(k/N)))
#                 sleep(1/24)
#             end
#             loop_closed_ref[] = true
#         end
#     end
# end
#
#
#
# function animate(F::FourBar,
#         angles1::Vector{NTuple{3,ComplexF64}},
#         angles2::Vector{NTuple{3,ComplexF64}},
#         coupler_points::Union{Nothing,Vector{ComplexF64}}=nothing;
#         show_axis=true,
#         color=:CADETBLUE,
#         color2=:DODGERBLUE,
#         fps=24, seconds=5,
#         loop::Union{Int,Bool} = false,
#         filename::Union{String,Nothing}=nothing)
#
#     positions1 = four_bar_positions(F, angles1)
#     positions2 = four_bar_positions(F, angles2)
#     limits = compute_limits([positions1; positions2])
#     markersize = max(limits.widths...) / 50
#
#     scene = Scene(limits=limits, resolution = (1500,1500), scale_plot=false);
#
#     # Draw mechanism ankers and coupler points
#     static!(scene, F, coupler_points; markersize=markersize)
#
#     source1, loop_closed_ref1 = add_mechanism!(scene, positions1; show_axis=show_axis,
#                 color=color, markersize=markersize)
#     source2, loop_closed_ref2 = add_mechanism!(scene, positions2;
#                 color=color2, show_axis=show_axis, markersize=markersize)
#
#
#     itp1 = interpolate_curve(positions1)
#     itp2 = interpolate_curve(positions2)
#
#     N = seconds * fps
#     if filename !== nothing
#         record(scene, filename, 1:N; framerate=fps) do k
#             push!(source1, round(Int, itp1(k/N)))
#             push!(source2, round(Int, itp2(k/N)))
#         end
#     else
#         display(scene)
#         if loop === false
#             nloops = 1
#         elseif loop === true
#             nloops = typemax(Int64)
#         else
#             nloops = loop::Int
#         end
#
#         for _ in 1:nloops
#             for k in 0:N
#                 push!(source1, round(Int, itp1(k/N)))
#                 push!(source2, round(Int, itp2(k/N)))
#                 sleep(1/fps)
#             end
#             loop_closed_ref1[] = true
#             loop_closed_ref2[] = true
#         end
#     end
# end
#
#
# function interpolate_curve(pos)
#     n = length(pos)
#     partials = [0.0]
#     l = 0.0
#     for i in 2:n
#         l += norm(pos[i].P - pos[i-1].P)
#         push!(partials, l)
#     end
#     # normalize partials to length 1
#     partials ./= partials[end]
#     itp = interpolate((partials,), 1:n, Gridded(Linear()))
#     itp
# end
#
# function add_mechanism!(scene, positions; color=:CADETBLUE, markersize=1.0, show_axis=show_axis)
#     loop_closed = Ref(false)
#     source = Node(1)
#     P = map(pos -> pos.P, positions)
#     curve_at(t) = loop_closed[] ? (@view P[1:end]) : view(P, 1:t)
#     if show_axis == false
#         Makie.lines!(scene, P, color = :transparent, show_axis=show_axis);
#     end
#     Makie.lines!(scene, lift(curve_at, source),
#                 color=color, linewidth=5, show_axis=show_axis);
#     fourbar_at = t -> begin
#         A, B, C, D, Pᵢ = positions[t]
#         [A,C,Pᵢ,D,B,D,C]
#     end
#     lines!(scene, lift(fourbar_at, source), color=:black,
#                 linewidth = 3, show_axis=show_axis)
#     scatter!(scene, lift(i -> [P[i]], source), color=color,
#         markersize=0.5markersize,
#         show_axis=show_axis)
#     scatter!(scene, lift(i -> [positions[i].D, positions[i].C], source), color=:black,
#         markersize=0.25markersize,
#         show_axis=show_axis)
#     source, loop_closed
# end

end # module
