#### Lastest change: 10.5.2024 by P. Öffner 

using DifferentialEquations
include("bandedpoly.jl")

u0func(x) = sin(pi*x)

# Performs a numerical test for the linear advection equation. Op is the used operator.
# N and K are the number of elements in x and y direction. v is the advection direction and dt the
# timestep.
function LinAdvNT(Op; N = 10, CFL = 0.01, tmax = 10.0)
        tspan = (0.0, tmax)                                                                     # the timespan for the time integrator
	dx = 2.0/N
	dt = CFL*dx/Op.N
        u0ar = zeros(N, Op.N)                                                    # The tensor initialized with the initial condition.
        Grid = zeros(N, Op.N)                                                 # A tensor with coordinates of all nodes.

        #initialize grid
        for n=1:N
		for i=1:Op.N
			Grid[n, i] = (n-0.5)*dx + (dx/2)*Op.cpoints[i]
                        u0ar[n, i] = u0func(Grid[n,i])
                end
        end

	# Periodic
	ul(u, x, t) = u[end, end]
	ur(u, x, t) = u[1, 1]

	p = LinAdv(0.0, 2.0, Grid, dx, dt, Op, ul, ur, N) 
                                                                                                # Solve the problem
        println("Solving the Problem")
        prob = ODEProblem(LinAdvLLF!, u0ar, tspan, p)
        sol = solve(prob, SSPRK33(), dt=dt)
        return p, sol                                                                           # and return the result
end

fLin(u) = u
# The local Lax-Friedrichs flux for a linear advection with speed v = (v_1, v_2), 
# internal value ui, external value uo and outer surface normal n.
fLinLLF(ul, ur) = 0.5*(fLin(ul) + fLin(ur) + (ul - ur))


function LinAdvLLF!(du, u, p, t)
	for n=1:p.N
        	du[n, :] = -p.Op.Dx*u[n, :]*2/p.dx
	end

	ffluxdiff = zeros(Op.N)
	for n=2:p.N-1
		
		ffluxdiff[1] = fLin(u[n, 1]) - fLinLLF(u[n-1,end], u[n,1])
		ffluxdiff[end] = fLin(u[n, end]) - fLinLLF(u[n,end], u[n+1,1])	
		du[n, :] += p.Op.Pi*(p.Op.Bx*ffluxdiff)*2/p.dx
	end

	n = 1
	ffluxdiff[1] = fLin(u[n, 1]) - fLinLLF(p.ul(u, p.Grid[1, 1], t), u[n,1])
	ffluxdiff[end] = fLin(u[n, end]) - fLinLLF(u[n,end], u[n+1,1])
	du[n, :] += p.Op.Pi*(p.Op.Bx*ffluxdiff)*2/p.dx

	n = p.N
	ffluxdiff[1] = fLin(u[n, 1]) - fLinLLF(u[n-1,end], u[n,1])
	ffluxdiff[end] = fLin(u[n, end]) - fLinLLF(u[n,end], p.ur(u, p.Grid[end, end], t))
	du[n, :] += p.Op.Pi*(p.Op.Bx*ffluxdiff)*2/p.dx

end
