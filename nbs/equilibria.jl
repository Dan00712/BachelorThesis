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

    using Random

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
    const single_z = let
        samples = 50

        rng = Random.seed!(0)
        logs = [-2, log10(5)+1]
        span = logs[2] - logs[1]
        foo = 10 .^ (rand(rng, samples) .* span .+ logs[1])
        #vcat(-1 .* foo, foo)
        foo
    end
    const Pz = Iterators.product([single_z for _ = 1:N]...)
    const d = 0 #params.R*100
    const κ = 18e4 * 2π
    const vd = let
        l = N/2 - 1
        [-l+k for k = 1:N]
    end
end;

# ╔═╡ 1f1935f8-dabf-4d77-957f-121eecaea100
w, z = let
    tol = 1e-8
    w = Float64[]
    z_parts = []
    vphi = [0]
    @progress for Δ in Ω
        for pz in Pz
            r = find_roots(collect(pz), vd, vphi, Δ = Δ, κ = κ, params = params)
            if length(z_parts) == 0 || !(([r ≈ z_parts[i] for i = 1:N] |> any))
                push!(z_parts, r)
                push!(w, Δ)
            end
        end
    end

    # map z :: Vector{Vector{...} l N} l M ↦ z :: Matrix M x N
    z = reduce(hcat, z_parts)'

    w, z
end

# ╔═╡ 3d1cb251-f943-42ce-918b-f85981aaaaa3
let
    datadir("N-Particle", "P1-test") |> mkpath
    @save datadir("N-Particle", "P1-test", "optim.jld2") w z
    p = datadir("N-Particle", "P1")
    mkpath(p)
    @save joinpath(p, "last.jld2") w z
    p = joinpath(p, "NP-$(now_nodots() ).jld2")
    @save p w z
    p
end

# ╔═╡ ff1ffcd4-4634-46e7-a1d1-ccebced9ad52
let
    condition(x) = x > 0 && x < 35
    w_ = [w[i] for (i, zi) in enumerate(collect(z)[:, 1]) if condition(zi)]
    scatter(
        w_ ./ 2π .* 1e-3,
        filter(collect(z)[:, 1]) do z
            condition(z)
        end;
        yaxis = :log,
        ylabel = "z/μm",
        xlabel = "Δ /2πHz",
    )
end

# ╔═╡ Cell order:
# ╠═da59901a-f1f6-11f0-aa2c-097a16ba3480
# ╠═f5cb8094-8bfa-4bd9-a29b-a8dcea76a7d8
# ╠═76e62065-9561-4601-ba2b-9b3bc9b91c6d
# ╠═6eb1ad31-f974-4c5e-82f9-d13a9aa96a13
# ╠═1f1935f8-dabf-4d77-957f-121eecaea100
# ╠═3d1cb251-f943-42ce-918b-f85981aaaaa3
# ╠═ff1ffcd4-4634-46e7-a1d1-ccebced9ad52
