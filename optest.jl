include("LFSBcon.jl")
include("LFSBPcon2.jl")
include("structs.jl")
include("bfs.jl")

using Plots
using LinearAlgebra
using SummationByPartsOperators
using Latexify

##Comparison of FSBP Operators Part I 

xL = -1.0
xR = 1.0
N = 9
cpoints = range(xL, xR, length=N)
bwidth = 10
funcs = [fpow(0), fpow(1), fpow(2), fpow(3), fpow(4)]#, fpow(5), fpow(6), fpow(7)]#, fpow(8)]
dxfuncs = [dfpow(0), dfpow(1), dfpow(2), dfpow(3), dfpow(4)]#, dfpow(5), dfpow(6), dfpow(7)]#, dfpow(8)]



xm = collect(-1.0:0.1:1.0)



Op = FSBP(xL, xR, N, cpoints, bwidth, funcs, dxfuncs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
S, P, B, Dx, V, dV, fsl, Tm = MakeBFSBPOpti(Op)


DxC =  derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=4, xmin=-1.0, xmax=1.0, N=9)
Pc = mass_matrix(DxC)
Qc = Pc*Matrix(DxC)

fc = open("fsbpP.tex", "w")
write(fc, latexify(P, fmt = FancyNumberFormatter(4)))
close(fc)

fc = open("fsbpQ.tex", "w")
write(fc, latexify(P*Dx, fmt = FancyNumberFormatter(4)))
close(fc)

cc = open("csbpP.tex", "w")
write(cc, latexify(Matrix(Pc), fmt = FancyNumberFormatter(4)))
close(cc)

cc = open("csbpQ.tex", "w")
write(cc, latexify(Pc*Matrix(DxC), fmt = FancyNumberFormatter(4)))
close(cc)




