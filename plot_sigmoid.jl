using Plots
logistic(x)=1/(1+exp(-x)) 
glogistic(x,aa)=(1+exp(-x))^(-aa) 
x = range(-10, 10, length=1000)
y1 = logistic.(x)
y2 = glogistic.(x,2)
y3 = glogistic.(x,0.5)

plot(x, [y1], title="Sigmoid Functions", label="Logistic.", ls=:dot , linewidth=3)
plot!(x, [y2 y3], title="Sigmoid Functions", label=[ "Gen. log.(2)" "Gen. log.(0.5)"] , linewidth=3)

savefig("Sigmoid.png")      # saves the CURRENT_PLOT as a .png
