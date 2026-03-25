### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 99ed7992-fb8d-11f0-8324-8f2ec916a976
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
if isinteractive()
    plotlyjs()
end

using CavityEquilibria
using CavityEquilibria.Parameters
using CavityEquilibria.Laser
using CavityEquilibria.RootFinding
using CavityEquilibria.Util
end

# ╔═╡ a82a841c-4f3f-4d42-b30f-4b32a63963eb
const Φ = range(-π/2, π, 50)

# ╔═╡ 4faf0797-817b-4987-9007-b0d26f171d4e
begin 
const params = DEFAULT_PARAMS


const N = 2
const sz = let
    VALS = 15
    Random.seed!(0)

    #range is [-2, 2]
    logsz = rand(VALS) .* 4 .- 2
    1e-6 .* vcat(exp10.(logsz), -1 .* exp10.(logsz))
end
const vz = Iterators.product([sz for _ = 1:N]...)
const d = 1.5e-6    # m === 1.5 μm
const κ = 18e4 * 2π
const Δ = let
    # kc * d/2 = 2 π
    kc = 4π/d
    # c = ωc/kc
    ωc = params.c * kc
    
    ωc - params.ω0
end
const vd = [-d/2, d/2]
end;

# ╔═╡ 5726ef17-12a4-4f42-baa6-8e7e2538e6cb
Δ

# ╔═╡ eb0b1209-5d62-4672-b37f-1b6f03e3cf89
round_better(y; f=10) = round(y*f)/f

# ╔═╡ 3ea86749-c38c-43e2-a585-480308d6d3ef
# ╠═╡ disabled = true
#=╠═╡
let
a = @animate for ϕ_ in Φ
	z = []
    for zg in vz
        try
            r = find_roots(collect(zg), vd, [0, ϕ_]; Δ = Δ, κ = κ, params = params)

            if (length(z) == 0 || !any(v -> isapprox(v, r), z)) && all(abs.(r) .< 1)
                push!(z, r)
            end
        catch err
            if !(err isa ConvergenceError)
                rethrow(err)
            end
        end
    end
	z = Iterators.reduce(hcat, z)

	p = plot(;
			 xlims=(-1, 1),
			 ylims=(-1, 1)
	)
	scatter!(p,
		z[1, :] .* 1e6,
		z[2, :] .* 1e6,
		label="ϕ=$(round_better(ϕ_))"
	)
end
	gif(a, fps=10)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═99ed7992-fb8d-11f0-8324-8f2ec916a976
# ╠═a82a841c-4f3f-4d42-b30f-4b32a63963eb
# ╠═4faf0797-817b-4987-9007-b0d26f171d4e
# ╠═5726ef17-12a4-4f42-baa6-8e7e2538e6cb
# ╠═eb0b1209-5d62-4672-b37f-1b6f03e3cf89
# ╠═3ea86749-c38c-43e2-a585-480308d6d3ef
