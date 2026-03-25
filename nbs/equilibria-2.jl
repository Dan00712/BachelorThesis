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
const Φ = range(0, 2π, 50)

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
const κ = 18e4 * 2π
const λc = 1549.999999e-9 # 1500nm
const Δ = let
	ωc = 2π * params.c/λc

	ωc - params.ω0
end
d = let
	n = 2

	n * λc
end
const vd = [-d, d]
end;

# ╔═╡ 7f5823ac-fe8c-40c5-8f27-a54a2feff0fc
params.c, params.λ0, λc

# ╔═╡ e6360f46-1a3a-4dda-8cab-871df54030ce
Δ, d

# ╔═╡ eb0b1209-5d62-4672-b37f-1b6f03e3cf89
round_better(y; f=10) = round(y*f)/f

# ╔═╡ 3ea86749-c38c-43e2-a585-480308d6d3ef
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
			 xlims=(-25, 25),
			 ylims=(-25, 25),
	)
	scatter!(p,
		z[1, :] .* 1e6,
		z[2, :] .* 1e6,
		label="ϕ=$(round_better(ϕ_))",
		xlabel="z1/μm",
		ylabel="z2/μm",
		legend=:topright
	)
end
	gif(a, fps=2)
end

# ╔═╡ Cell order:
# ╠═99ed7992-fb8d-11f0-8324-8f2ec916a976
# ╠═a82a841c-4f3f-4d42-b30f-4b32a63963eb
# ╠═4faf0797-817b-4987-9007-b0d26f171d4e
# ╠═7f5823ac-fe8c-40c5-8f27-a54a2feff0fc
# ╠═e6360f46-1a3a-4dda-8cab-871df54030ce
# ╠═eb0b1209-5d62-4672-b37f-1b6f03e3cf89
# ╠═3ea86749-c38c-43e2-a585-480308d6d3ef
