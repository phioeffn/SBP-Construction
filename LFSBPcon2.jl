#### Lastest change: 13.5.2024 by P. Öffner 

include("LFSBcon.jl")
using Optim
#using Enzyme
#using ReverseDiff
using Zygote


function sig(x)
	return 1.0/(1.0 + exp(-x))
end

function softmax(x)
	return (1.0E-6 .+ sig.(x))./sum(1.0E-6 .+ sig.(x))
end

function softmax_2(x)
	return (1.0E-6 .+ sig.(x))
end

###This is our solver to calculate the FSBP operators with the new construciton procedure. 


function MakeBFSBPOpti(Op)
	M = DSOBM(Op.funcs, Op.dxfuncs, Op.cpoints)
	#M = diagm(ones(length(Op.funcs)))
	Op.dob = Array{Any}(undef, length(Op.funcs))
        Op.dobx = Array{Any}(undef, length(Op.funcs))
        for k = 1:length(Op.funcs)
                Op.dob[k] = x->DOBF(M, Op.funcs, k, x)
                Op.dobx[k] = x->DOBFx(M, Op.dxfuncs, k, x)
        end

        # Creation of all variables
        V, dV = MakeVander(Op), MakeDxVander(Op)
        W = vcat(V, -dV)
        S = zeros(Op.N, Op.N)
        P = diagm(ones(Op.N)) *(Op.xR - Op.xL)/N
        X = hcat(S, P)
	C = ones(1, Op.N)               # Constant one
        B = zeros(Op.N, Op.N)
	B[1, 1] = -1.0
        B[end, end] = 1.0
        R = B*V/2
	mu = (Op.xR - Op.xL)
        #mu = 1   #If we set something not the the length. Have to included if we do not requ9re the data..

        Tm = zeros(Op.N, Op.N) # A "Template" matrix storing the non-zero elements allowed in S
        for k=1:Op.N
                for l=1:Op.N
                        if abs(k-l) < Op.bwidth
                                Tm[k, l] = 1.0
                        end
                end
        end 

	s = zeros(Op.N, Op.N)
	p = zeros(Op.N)
	function f(s, p)
		# Projection onto Diagonal quadrature
	     	w = mu*softmax(p) # The version where we include that quadrature is exact for constants
                # w = mu*softmax_2(p) # Version where we are not exact for constants! Have to be included if we like to test them. 
                # w = exp.(p)
                # Skew symmetrie
                S = (Tm.*s - transpose(Tm.*s)) / 2.0
		return norm(S*V - diagm(w)*dV + R)^2
	end

        #fsl(s, p) = norm((Tm.*s- transpose(Tm.*s))/2.0*V - diagm(mu*softmax_2(p))*dV + R)^2 # andere Formel //  Version where we are not exact for constants
                # w = exp.(p)
	fsl(s, p) = norm((Tm.*s- transpose(Tm.*s))/2.0*V - diagm(mu*softmax(p))*dV + R)^2 # The version where we include that quadrature is exact for constants
        
  
	fres(S, P) = norm((Tm.*S- transpose(Tm.*S))/2.0*V - P*dV + R)^2
	fsv(x) = fsl(x[:, 1:end-1], x[:, end])
	x = hcat(s, p)
	fg = fsv'


        res = optimize(fsv, fg, x, LBFGS(), Optim.Options(g_tol=1.0E-16, iterations=50000); inplace = false)
        print(res)
	Y = res.minimizer
	S, P = Y[:, 1:end-1], diagm(mu*softmax(Y[:, end]))
        #S, P = Y[:, 1:end-1], diagm(softmax_2(Y[:, end])) # Version not exact for constants not working, or maybe
	S = (S - transpose(S)) / 2.0
	println("Resulting objective value ", fsv(Y), " ", fres(S, P))
	Q = S + B/2
        Dx = pinv(P)*Q
        Op.Dx = Dx
        Op.P = P
        Op.Pi = inv(P)
        Op.Bx = B
        return S, P, B, Dx, V, dV, fres, Tm
end

