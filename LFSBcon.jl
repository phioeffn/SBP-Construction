#### Lastest change: 10.5.2024 by P. Öffner 


using LinearAlgebra

# Returns a Vandermonde matrix for the basis and collocation points in Op
function MakeVander(Op)
        nfunc = length(Op.dob)
        V = zeros(Op.N, nfunc)
        for k=1:Op.N
                for l=1:nfunc
                        V[k, l] = Op.dob[l](Op.cpoints[k])
                end
        end
        return V
end

# Returns the Vandermonde matrix for the derivatives in x direction
# of the basis and collocation points in Op.
function MakeDxVander(Op)
        nfunc = length(Op.dobx)
        dV = zeros(Op.N, nfunc)
        for k=1:Op.N
                for l=1:nfunc
                        dV[k, l] = Op.dobx[l](Op.cpoints[k])
                end
        end
        return dV
end

# Returns locally stable vandermonde matrices. 
# Several linear dependend columns are added to allow
# the condition number of restrictions to be lower.
function MakeVanderStable(Op)
	nfunc = length(Op.funcs)
	# we copy the basis functions NCopy times. The new copys are
	# locally orthogonal
	NCopy = floor(Int, Op.N/Op.bwidth)
	# coeffs saves the linear combination coefficients used in every copy
	coeffs = zeros(NCopy, nfunc, nfunc)
	# The local vandermonde matrix on every stencil
	lv = zeros(Op.bwidth, nfunc)

	# Calculation of the needed coefficients
	for n = 1:NCopy
		# Calculation of the local vandermonde matrix
		for k = 1:Op.bwidth
			for l = 1:nfunc
				lv[k, l] = Op.funcs[l](Op.cpoints[k + (n-1)*Op.bwidth])
			end
		end
		vsvd = svd(lv)
		coeffs[n, :, :] = vsvd.Vt#vsvd.U*diagm(vsvd.S)
	end

	# Definition of the new functions
	cfuncs = Array{Any}(undef, nfunc*NCopy)
	for n=1:NCopy
		for k = 1:nfunc
			cfuncs[(n-1)*nfunc+k] = function(x)
				val = 0.0
				for l=1:nfunc
					val = val + Op.funcs[l](x)*coeffs[n,l, k]
				end
				return val
			end
		end
	end
	
	# Calculation of the new Vandermonde Matrix
	V = zeros(Op.N, nfunc*NCopy)
	for k = 1:Op.N
		for l=1:nfunc*NCopy
			V[k, l] = cfuncs[l](Op.cpoints[k])
		end
	end
	return V	
end
