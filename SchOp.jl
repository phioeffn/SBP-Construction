#### Philipp Öffner 
#### Lastest change: 10.5.2024

###Function for the Schrödinger example

include("bfs.jl")
include("structs.jl")
include("LFSBPcon2.jl")

using Plots
using FastGaussQuadrature


xl = -10.0
xr = 10.0

N = 100


cpoints = range(xl, xr, length=N)
bwidth = 100

funcs = [fpow(0), fpow(1)]
dxfuncs = [dfpow(0), dfpow(1)]



for n=0:10
	global funcs = vcat(funcs, [x->Phi(n, x)])
	global dxfuncs = vcat(dxfuncs, [x->dPhi(n, x)])
end

Op = FSBP(xl, xr, N, cpoints, bwidth, funcs, dxfuncs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
S, P, B, Dx, V, dV, fsl, Tm = MakeBFSBPOpti(Op)

