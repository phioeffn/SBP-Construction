#### Lastest change: 10.5.2024 by P. Öffner 



include("LFSBcon.jl")
include("LFSBPcon2.jl")
include("structs.jl")
include("bfs.jl")

using Plots

xL = -10.0
xR = 10.0
N = 100
cpoints = range(xL, xR, length=N)
bwidth = 100
funcs = [fpow(0), fpow(1), fpow(2), fpow(3)]#, fpow(4), fpow(5), fpow(6), fpow(7)]#, fpow(8)]
dxfuncs = [dfpow(0), dfpow(1), dfpow(2), dfpow(3)]#, dfpow(4), dfpow(5), dfpow(6), dfpow(7)]#, dfpow(8)]
d = 3.0


pOp = FSBP(xL, xR, N, cpoints, bwidth, funcs, dxfuncs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
N2 = 150 #150
cpoints2 = range(xL, xR, length=N2)
pOp2 = FSBP(xL, xR, N2, cpoints2, bwidth, funcs, dxfuncs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


DxC =  derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=4, xmin=-10.0, xmax=10.0, N=N)
Pc = mass_matrix(DxC)
Qc = Pc*Matrix(DxC)
pOp.Dx = Matrix(DxC)
pOp.P = Matrix(Pc)
pOp.Pi = inv(Matrix(Pc))
pOp.Bx = Qc + transpose(Qc)

DxC =  derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=4, xmin=-10.0, xmax=10.0, N=N2)
Pc = mass_matrix(DxC)
Qc = Pc*Matrix(DxC)
pOp2.Dx = Matrix(DxC)
pOp2.P = Matrix(Pc)
pOp2.Pi = inv(Matrix(Pc))
pOp2.Bx = Qc + transpose(Qc)



