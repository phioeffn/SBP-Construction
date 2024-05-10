include("bandedpoly.jl")

#include("SchOp.jl")

using DifferentialEquations
using Plots
#using SummationByPartsOperators
#D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=4, xmin=-10.0, xmax=10.0, N=100)
#Dx = Matrix(D)
#M = mass_matrix(D)

function SchrödingerNT(Op, u0func; cfl = 0.01, tmax = 10.0)
	tspan = (0.0, tmax)
	dx = (Op.xR - Op.xL)/Op.N
	dt = cfl*dx
	u0ar = zeros(Complex, Op.N)
	Grid = Op.cpoints[:]
	for n=1:N
		u0ar[n] = u0func(Grid[n])
	end
	S = Op.P*Op.Dx
	V = diagm(Op.cpoints)
	V = V^2
	p = (-im*Op.Pi*transpose(S)*Op.Dx-im*V)
#	p = (-im*inv(M)*transpose(M*Dx)*Dx - im*V)
	println("Solving the Problem")
	prob = ODEProblem(Schrödinger!, u0ar, tspan, p)
	sol = solve(prob, SSPRK33(), dt = dt)
	return p, sol
end

function Schrödinger!(du, u, p, t)
	du[:] = p*u
end
