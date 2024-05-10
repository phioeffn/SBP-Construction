
#### Lastest change: 10.5.2024 by P. Öffner 


# naive basis functions /polynomials 
fpow(k) = x->x^k
dfpow(k) = x->k > 0 ? k*x^(k-1) : 0.0


# naive basis functions /exponential functions
fexp(k) = x->exp(k*x)
dfexp(k) = x->k > 0 ? k*exp(x*k) : 0.0

# Legendre bais polynomials
# at point x, counting from 0
function Pleg(n, x)
        if n == 0
                return 1.0
        elseif n == 1
                return x
        else
                return ((2*n-1)*x*Pleg(n-1, x) - (n-1)*Pleg(n-2, x))/n
        end
end

# recursively gives the derivative of
# Lagrange polynomial $n$, counting from zero,
# at position x
function dPleg(n, x)
        if n == 0
                return 0.0
        elseif n == 1
                return 1.0
        else
                return ((2*n-1)*(Pleg(n-1, x) + x*dPleg(n-1, x)) - (n-1)*dPleg(n-2, x))/n
        end
end

fleg(k) = x->Pleg(k, x)
dfleg(k) = x->dPleg(k, x)

fleg(k, xm, w) = x->Pleg(k, (x-xm)/w)
dfleg(k, xm, w) = x->dPleg(k, (x-xm)/w)/w

# A radial basis function, formed by a gaussian centered at -1 + 2*k/N in the reference element
GBF(k, N, d) = x->exp(-(-1.0 + 2*k/N - x)^2/d^2)
dGBF(k, N, d) = x->exp(-(-1.0 + 2*k/N - x)^2/d^2)*2/d^2*(-1.0 + 2*k/N - x)


# Hermite polynomials
function H(n, x)
	if n == 0
		return 1
	elseif n== 1
		return 2*x
	else
		return 2*x*H(n-1, x)-2*(n-1)*H(n-2, x)
	end
end

function dH(n, x)
	if n > 0
		return 2*n*H(n-1, x)
	else 
		return 0.0
	end
end

function Phi(n, x)
	return exp(-x^2/2)/sqrt(2^n*factorial(n)*sqrt(pi))*H(n, x)
end

function dPhi(n, x)
	return (-exp(-x^2/2)*x*H(n, x) + exp(-x^2/2)*dH(n, x))/sqrt(2^n*factorial(n)*sqrt(pi))
end

# functions that return a set of DOBs from a set of basis functions.

# Calculates a discrete orthogonal basis from the monomials, saved in a matrix M
# Every row column in M contains the coefficients for the monomial basis.
function DOBM(fcts, pts)
                n = length(pts)
		m = length(fcts)
                M = zeros(m, m)
                # space for the new orthogonal basis
		O = zeros(n, m)
		a = zeros(n)

                M .= 0.0
		O .= 0.0
                # loop over all basis functions
                for i=1:m
			M[i, i] = 1.0
                        # evaluate the function that should be evaluated
                        for j=1:n
				a[j] = fcts[i](pts[j])
                        end
                        # loop over all previous new base vectors
                        for j=1:i-1
				M[:, i] = M[:, i] - dot(a, O[:, j])/dot(O[:, j], O[:, j])*M[:, j]
                        end
                        # save the new base vectors
                        for j=1:n
				O[j, i] = DOBF(M, fcts, i, pts[j])
                        end
                        # Normalize
			r = norm(O[:, i])
			O[:, i] = O[:, i] / r
			M[:, i] = M[:, i] / r

                end
		return M
end

# Calculates a discrete orthogonal basis from the monomials, saved in a matrix M
# Every row column in M contains the coefficients for the monomial basis.
function DSOBM(fcts, dxfcts, pts)
                n = length(pts)
                m = length(fcts)
                M = zeros(m, m)
                # space for the new orthogonal basis
                O = zeros(n, m)
		P = zeros(n, m)
                a = zeros(n)
		b = zeros(n)

                M .= 0.0
                O .= 0.0
		P .= 0.0
                # loop over all basis functions
                for i=1:m
                        M[i, i] = 1.0
                        # evaluate the function that should be evaluated
                        for j=1:n
                                a[j] = fcts[i](pts[j])
				b[j] = -dxfcts[i](pts[j])
                        end
                        # loop over all previous new base vectors
                        for j=1:i-1
				M[:, i] = M[:, i] - (dot(a, O[:, j]) + dot(b, P[:, j]))/(dot(O[:, j], O[:, j]) + dot(P[:, j], P[:, j]))*M[:, j]
                        end
                        # save the new base vectors
                        for j=1:n
                                O[j, i] = DOBF(M, fcts, i, pts[j])
				P[j, i] = -DOBF(M, dxfcts, i, pts[j])
                        end
                        # Normalize
			r = sqrt(dot(O[:, i], O[:, i]) + dot(P[:, i], P[:, i]))
                        O[:, i] = O[:, i] / r
                        M[:, i] = M[:, i] / r
			P[:, i] = P[:, i] / r

                end
                return M
end


# Evaluates the Discrete Orthogonal basis function k at x
# k is one-indexed
function DOBF(M, fcts, k, x)
        val = 0.0
        for i=1:k
        	val = val + fcts[i](x)*M[i, k]
        end
	return val
end

# Evaluates the k-th Discrete Orthogonal basis functions derivative  at x
# k is one-indexed
 function DOBFx(M, dxfcts, k, x)
        val = 0.0
        for i=1:k
                val = val + dxfcts[i](x)*M[i, k]
        end
        return val
end

