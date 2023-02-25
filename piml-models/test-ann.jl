# regression ann flux

using Flux, Plots

x = hcat(collect(Float32, -3:0.1:3)...)

f(x) = @. 3x + 2;

y = f(x)

x = x .* reshape(rand(Float32, 61), (1, 61));

plot(vec(x), vec(y), lw = 3, seriestype = :scatter, label = "", title = "Generated data", xlabel = "x", ylabel= "y");





