### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ da59901a-f1f6-11f0-aa2c-097a16ba3480
begin
    using DrWatson
    HOME = ENV["HOME"]
    DrWatson.quickactivate(
        joinpath(HOME, "projects", "bac", "CavityEquilibria"),
        "CavityEquilibria",
    )

    using IntervalArithmetic
    using StaticArrays
    using IntervalRootFinding

    using JLD2
    using Plots
    using ProgressLogging

    using CavityEquilibria
    using CavityEquilibria.Parameters
    using CavityEquilibria.Laser
    using CavityEquilibria.RootFinding
    using CavityEquilibria.Util
end

# ╔═╡ f5cb8094-8bfa-4bd9-a29b-a8dcea76a7d8
const params = DEFAULT_PARAMS

# ╔═╡ 76e62065-9561-4601-ba2b-9b3bc9b91c6d
Ω = vcat(
    range(1, 400, length = 100) .* 1e3*2π,  #.|> x-> 10^x,
) |> sort

# ╔═╡ 6eb1ad31-f974-4c5e-82f9-d13a9aa96a13
begin
    const N = 1
    const ivz = [interval(10^-2, 5e1) for _ = 1:N]
    const log10ivz = log10.(ivz)
    const d = 0 #params.R*100
    const κ = 18e4 * 2π
    const vd = [0]
end;

# ╔═╡ 1f1935f8-dabf-4d77-957f-121eecaea100
v, err, w = let
    tol = 1e-8
    w = []
    v = []
    err = []
    vphi = [0]
    @progress for Δ in Ω
        #Δ = prepare_freq(Δ_; params=params)

        tmp = find_roots_log10(log10ivz, vd, vphi; Δ = Δ, κ = κ, params = params)
        v_ = tmp .|> (x -> mid(x))
        err_ = tmp .|> x -> IntervalArithmetic.radius(x)
        append!(w, [Δ for _ in v_])
        append!(v, v_)
        append!(err, err_)
    end
    #v, w

    v, err, w
end

# ╔═╡ 3d1cb251-f943-42ce-918b-f85981aaaaa3
let p = datadir("N-Particle", "P1")
    mkpath(p)
    p = joinpath(p, "NP$(now_nodots() )-log.jld2")
    @save p v err w
    p
end

# ╔═╡ ff1ffcd4-4634-46e7-a1d1-ccebced9ad52
scatter(w/2π/1e3, v, yaxis = :log, yerror = err)

# ╔═╡ Cell order:
# ╠═da59901a-f1f6-11f0-aa2c-097a16ba3480
# ╠═f5cb8094-8bfa-4bd9-a29b-a8dcea76a7d8
# ╠═76e62065-9561-4601-ba2b-9b3bc9b91c6d
# ╠═6eb1ad31-f974-4c5e-82f9-d13a9aa96a13
# ╠═1f1935f8-dabf-4d77-957f-121eecaea100
# ╠═3d1cb251-f943-42ce-918b-f85981aaaaa3
# ╠═ff1ffcd4-4634-46e7-a1d1-ccebced9ad52
