#### Lastest change: 10.5.2024 by P. Öffner 



include("LFSBcon.jl")
include("LFSBPcon2.jl")
include("structs.jl")
include("bfs.jl")


using GaussQuadrature
using Plots



xL = -1.0 # Boundary elements normally -1 is always used execpt for one example 
xR = 1.0
N = 20# #20 # Have to be changed for higher dimensions 

### Different point distributions 

# Uniform points: 
cpoints = range(xL, xR, length=N)

#Gauss-Lobatto nodes: 
#cpoints, w = legendre(5,both) #for the test with basis 6.... so it is really the case ...

## Chebyshev-Lobatto nodes:

#cpoints, w = chebyshev(20,1,both)  

# Random point selection 

#Random number points

 # vv= zeros(N)
 # vv[1]=xL
 # vv[N]=xR
 # ra = rand(N-2).*2 
 # ra = ra.-1.
 # ra= sort(ra)
#for i in 2:N-1
#vv[i]=ra[i-1]
#end
#cpoints = vv


### First exammple
bwidth = 100
#funcs = [fpow(0), fpow(1), fpow(2), fpow(3), fpow(4), fpow(5), fpow(6), fpow(7)]#, fpow(8)]
#dxfuncs = [dfpow(0), dfpow(1), dfpow(2), dfpow(3), dfpow(4), dfpow(5), dfpow(6), dfpow(7)]#, dfpow(8)]
### Test not included in the paper 
#d = 3.0
#funcs = [ fpow(0), fpow(1),GBF(0, 5, d), GBF(1, 5, d), GBF(2, 5, d), GBF(3, 5, d), GBF(4, 5, d), GBF(5, 5, d)]
#dxfuncs = [dfpow(0), dfpow(1), dGBF(0, 5, d), dGBF(1, 5, d), dGBF(2, 5, d), dGBF(3, 5, d), dGBF(4, 5, d), dGBF(5, 5, d)]

### Test from the paper 
d = 3.0
funcs = [  fpow(1), GBF(0, 5, d), GBF(1, 5, d)]
dxfuncs = [ dfpow(1), dGBF(0, 5, d), dGBF(1, 5, d)]

### Test not included in the paper 

#funcs = [  fpow(0), fpow(1), GBF(0, 5, d), GBF(1, 5, d)]
#dxfuncs = [ dfpow(0), dfpow(1), dGBF(0, 5, d), dGBF(1, 5, d)]

#### Basis number six: Polynomial basis to test if we have Gauss-Lobatto is exactly the Gauss-Quadrature SBP 
#### They have to be extended for higher degrees. Additional test from the paper

#funcs = [fpow(0), fpow(1), fpow(2), fpow(3)]
#dxfuncs = [dfpow(0), dfpow(1), dfpow(2), dfpow(3)]


#### Test for checking the example for SINUM Paper: 

#funcs = [fpow(0), fpow(1), fexp(1)]
#dxfuncs = [dfpow(0), dfpow(1), dfexp(1)]

#### Further Tests not included in the paper 


#for n=-5:15
#	global funcs = vcat(funcs, [GBF(n, 5, d)])
#	global dxfuncs = vcat(dxfuncs, [dGBF(n, 5, d)])
#end

#funcs = [fleg(0), fleg(1), fleg(2), fleg(3), fleg(4), fleg(5), fleg(6), fleg(7)]
#dxfuncs = [dfleg(0), dfleg(1), dfleg(2), dfleg(3), dfleg(4), dfleg(5), dfleg(6), dfleg(7)]

#funcs = [fleg(0)]
#dxfuncs = [dfleg(0)]

xm = collect(-1.0:0.1:1.0)


Op = FSBP(xL, xR, N, cpoints, bwidth, funcs, dxfuncs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
S, P, B, Dx, V, dV, fsl, Tm = MakeBFSBPOpti(Op)




### Safe stuff in external file  and compare the operators. 
# fc = open("fsbpP_compare_4.tex", "w")
# write(fc, latexify(P, fmt = FancyNumberFormatter(4)))
# close(fc)

# fc = open("fsbpQ_compare_4.tex", "w")
# write(fc, latexify(P*Dx, fmt = FancyNumberFormatter(4)))
# close(fc)

