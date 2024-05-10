#### Lastest change: 10.5.2024 by P. Öffner 

###Advection in different settings! Pictures have to be adapted!

using Plots
include("LinAdv.jl")
#For the error convergence

Nmax = 6
Narr = collect(2:Nmax)
errarr = zeros(length(Narr))
for n=1:length(Narr)
	p, sol = LinAdvNT(Op, N=Narr[n])
	delta = sol(0.0) - sol(10.0)
	errarr[n] = sqrt(dot(delta, delta)/Narr[n])
end


#errarr1=errarr
#plot!(Narr, errarr1,  xaxis=:log, yaxis=:log, label="Error rand", lw=3, ls=:dot)
plot(Narr, errarr,  xaxis=:log, yaxis=:log, label="Error equi", lw=2)
plot!([Narr[1], Narr[end]], [errarr[end]*(Narr[end]/Narr[1])^2, errarr[end]], label="Order "*string(2))
plot!([Narr[1], Narr[end]], [errarr[end]*(Narr[end]/Narr[1])^3, errarr[end]], label="Order "*string(3))
plot!([Narr[1], Narr[end]], [errarr[end]*(Narr[end]/Narr[1])^4, errarr[end]], label="Order "*string(4))
xlabel!("N")
ylabel!("error")

#ylims!(0.008, 0.5) #for the polynomial exact....

savefig("Random_points_20_basis6.png")   #### Safefig_adapted



for p=6:8
	plot!([Narr[1], Narr[end]], [errarr[end]*(Narr[end]/Narr[1])^p, errarr[end]], label="Order "*string(p))
end

