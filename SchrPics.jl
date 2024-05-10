#### Lastest change: 10.5.2024 by P. Öffner 

###Main of the the Schrödinger example

include("HarmOsziR.jl")
include("SchpOp.jl")
include("SchOp.jl")
using LaTeXStrings


p, sol = SchrödingerNT(Op, x->[exp(-(x-2.5)^2), 0.0], tmax = 100.0, cfl=0.03)
p2, sol2 = SchrödingerNT(pOp, x->[exp(-(x-2.5)^2), 0.0], tmax = 100.0, cfl=0.03)
p3, sol3 = SchrödingerNT(pOp2, x->[exp(-(x-2.5)^2), 0.0], tmax = 100.0, cfl=0.03)


pot = TotalPropOT(p, sol)
pot2 = TotalPropOT(p2, sol2)

for k=0:4
	Plots.plot(Op.cpoints, Prop(sol(pi/8*k)), ylims = (0.0, 1.0), ylabel = L"u_1^2 + u_2^2", xlabel = "x", label = "FSBP")
	Plots.plot!(pOp.cpoints, Prop(sol2(pi/8*k)), ylims = (0.0, 1.0), ylabel = L"u_1^2 + u_2^2", xlabel = "x", label = "SBP")
	Plots.savefig("Harmtpi" * string(k) *".png")
end

Plots.plot(Op.cpoints, Prop(sol(100.0)), ylims = (0.0, 1.0), ylabel = L"u_1^2 + u_2^2", xlabel = "x", label = "FSBP @ N=100")
Plots.plot!(pOp.cpoints, Prop(sol2(100.0)), ylims = (0.0, 1.0), ylabel = L"u_1^2 + u_2^2", xlabel = "x", label = "SBP @ N=100")
Plots.plot!(pOp2.cpoints, Prop(sol3(100.0)), ylims = (0.0, 1.0), ylabel = L"u_1^2 + u_2^2", xlabel = "x", label = "SBP @ N=150")

Plots.savefig("Harmt100.png")


Plots.plot(sol.t, pot, ylabel = L"\int u_1^2 + u_2^2 \mathrm{d} x", xlabel = L"t", label = "FSBP")
Plots.plot!(sol2.t, pot2, ylabel = L"\int u_1^2 + u_2^2 \mathrm{d} x", xlabel = L"t", label ="SBP")

Plots.savefig("Harmpot.png")
#AnimateSol(sol, "HarmonicOsziPoly.gif", ylabel = L"u_1^2 + u_2^2", xlabel = L"x")



