#### Lastest change: 10.5.2024 by P. Öffner 



### Used to calculate the Schrödinger example! The main is SchrPics.jl

using DifferentialEquations
using Plots
using SummationByPartsOperators




function SchrödingerNT(Op, u0func; cfl = 0.03, tmax = 10.0)
	tspan = (0.0, tmax)
	dx = (Op.xR - Op.xL)/Op.N
	dt = cfl*dx
	u0ar = zeros(2, Op.N)
	Grid = Op.cpoints[:]
	for n=1:N
		u0ar[:, n] = u0func(Grid[n])
	end
	tp = TotalProp(Op, u0ar)
	u0ar = u0ar / sqrt(tp)
	p = Op	
	println("Solving the Problem")
	prob = ODEProblem(SchrödingeR!, u0ar, tspan, p)
	sol = solve(prob, SSPRK33(), dt = dt)
	return p, sol
end

function SchrödingeR!(du, u, Op, t)
	du[1, :] = -Op.Dx*Op.Dx*u[2, :] + Op.cpoints.^2 .* u[2, :] + Op.Pi*Op.Bx*Op.Dx*u[2, :]
	du[2, :] = Op.Dx*Op.Dx*u[1, :] - Op.cpoints.^2 .* u[1, :] - Op.Pi*Op.Bx*Op.Dx*u[1, :]
end

function TotalProp(Op, u)
	return dot(u[1, :], Op.P*u[1, :]) + dot(u[2, :], Op.P*u[2, :])
end

## Produce a video!

function AnimateSol(sol, fname; xlabel, ylabel)
	kmax = length(sol.t)
	Pr = Prop(sol(0.0))
	ymin, ymax = minimum(Pr), maximum(Pr)
	anim = @animate for k=1:kmax
		Plots.plot(Op.cpoints, Prop(sol(sol.t[k])), ylims = (ymin, ymax), xlabel = xlabel, ylabel = ylabel)
	end
	gif(anim, fname, fps=30)
end

function TotalPropOT(Op, sol)
	kmax = length(sol.t)
	Pr = zeros(kmax)
	for k=1:kmax
		Pr[k] = TotalProp(Op, sol(sol.t[k]))
	end
	return Pr
end

Prop(u) = u[1, :].^2 + u[2, :].^2
