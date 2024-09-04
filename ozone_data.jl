using DifferentialEquations
using Plots

# Defining parameters
k1 = 0.01  # ppm^-1 day^-1
k2 = 0.005  # day^-1
α = 0.1  # strength of the complex term
β = 2π/24  # rad/day, representing a daily cycle
γ = 50  # ppm, baseline ozone concentration

# Defining concentrations for NOx and VOC
NOx = 30  # ppm
VOC = 0.4  # ppm

# Initial ozone concentration 
u0 = [0.1] # ppm

# Defining the complex term f(t, [O3])
f(t, O3) = α * (β*t) * (O3 - γ)^2

# Defining the ODE function
function ozone_ode!(du, u, p, t)
    O3 = u[1]
    du[1] = k1 * NOx * VOC - k2 * O3 + f(t, O3)
end

# ODE problem
N_days = 30
tspan = (0.0, Float64(N_days)) 
datasize = N_days
t = range(tspan[1],tspan[2],length=datasize)

prob = ODEProblem(ozone_ode!, u0, tspan)

# Solve the ODE
sol = Array(solve(prob, Tsit5(), u0=u0, saveat=t))  

# Ploting the results
O3 = sol[1,:]
plot( O3, xlabel="Time (days)", ylabel="Ozone Concentration (ppm)", 
     title="Simulated Ozone Concentration", legend=false)
