using WaterLily, ParametricBodies # weymouth
using StaticArrays, ForwardDiff # required for simulation
using Statistics, CSV, DataFrames, Plots # data management and calculation

include("plot.jl")
include("util.jl")

function run_singlefoilsim(;T=Float64,mem=Array) # runs simulation for a single foil ; you can call this function with no inputs.

 ########################### INPUTS###########################
    coeff = SA[.3,.2,.3]; # CST coefficients
    L = 32; # chord length
    St = 0.3; # Strouhal number - > determines heave frequency
    U = 1; # free stream velocity
    Re = 5000; # Reynolds number - > determines viscosity
    n = 10; # number of chordwise lengths
    m = 8; # number of vertical lengths
    nNose = 4.; # x location of nose
    cycles = 3; # number of cycles to run simulation for
    metrics = false; # plot metrics

    nose = SA[nNose*L,0.5f0m*L] # calculated nose location
    ω = T(π*St*U/L); # calculated heave frequency
    h₀ = T(L); # heave amplitude
    θ₀ = T(20pi/180); # pitch amplitude

    h(t) = h₀*sin(ω*t) # heave function 
    θfunc(t) = θ₀*cos(ω*t) # pitch angle function

 #########################################################

 ################### OUTPUTS #############################

 # will output a gif vizualizing the motion, a CSV file holding the non-dimensionalized forces, moment, and power.
 # they are nondimensionalized the standard way, F/.5ρU^2L, M/.5ρU^2L^2, P/.5ρU^3L (here ρ = 1, simulation is created by giving Re and U; viscosity is calculated based on the reynolds.)
 # if you toggle metrics to true in the inputs it will output the gif along with force Plots

 #########################################################

 ############### MAKE SIMULATION ################
    function map(x,t) 
        hev = SA[0,h(t)]
        θv = θ(t) 
        R = SA[cos(θv) -sin(θv); sin(θv) cos(θv)];# rotation matrix 
        ξ = R*(x-nose-hev) #transform coodrinates; move to origin and align with x-axis
        return SA[ξ[1],abs(ξ[2])]    # reflect to positive y

    end

    function CSTgen(s,t)
        zetaflag = 1;
        psi = (1-s)^2
        ncoeff = size(coeff)[1]-1
        Classfunc = (psi^.5 )*(1. - psi)^1.
        Shapefunc = 0; 
        for i in 1:ncoeff+1 
            K = factorial(ncoeff)/(factorial(i-1)*factorial(ncoeff-(i-1)))
            Shapefunc = Shapefunc + (coeff[i]*K)*psi^(i-1)*(1 - psi)^(ncoeff-(i-1))
        end
        zeta = Classfunc*Shapefunc
        return L*SA[psi,zeta*zetaflag]
    end

    function get_force(sim,t) # takes nd time ; pretty sure this is the BIGGEST time waster
        sim_step!(sim,t,remeasure=true,verbose=true) # Non-dimensional time
        #return -1*my∮nds(sim.flow.p,sim.flow.V,sim.body,t*sim.L/sim.U), -1*sum(moment∮nds(sim.flow.p,sim.flow.V,sim.body,SA[T(nNose*L),0.5f0m*L + h(t*sim.L/sim.U)],t*sim.L/sim.U)) # Real time thrust coefficent
        return  -1.0 .* foil∮nds(sim.flow.p,sim.flow.V,sim.body,SA[T(nNose*L),0.5f0m*L + h(t*sim.L/sim.U)],t*sim.L/sim.U)
    end

    body = ParametricBody(CSTgen,(0,1);map,T,mem) # create body
    sim = Simulation((n*L,m*L),(U,0),L;Δt= .01,ν=U*L/Re,body,T,mem) # create simulation

    data = DataFrame(t=[],Thrust=[],Side_Force=[],Moment=[],power=[],ypos=[],theta=[]) # create data frame
    dat = sim.flow.σ[inside(sim.flow.σ)] |> Array; # create array for holding vorticity data

    anim = @animate for t ∈ range(0,cycles * 2/St,length = 240) # main loop

        forces,moment = get_force(sim,t) # get raw forces and moment
        
        hdot = ForwardDiff.derivative(h,t*sim.L/sim.U) # heave velocity
        θdot = ForwardDiff.derivative(θfunc,t*sim.L/sim.U) # pitch velocity

        power = -1*forces[2]*hdot + moment*θdot # calculate power

        ypos = h(t*sim.L/sim.U)
        theta = θfunc(t*sim.L/sim.U)*180/pi

        push!(data,[t,forces[1]/(.5*sim.U^2*sim.L),forces[2]/(.5*sim.U^2*sim.L),moment/(.5*sim.U^2*sim.L^2),power/(.5*sim.U^3*sim.L),ypos,theta])

        get_omega!(dat,sim); #plot_vorticity(dat;clims=(-5,5))   

        if metrics == true
            l = @layout [
                a{} [grid(6,1)]
            ];

            p1 = plot_vorticity(dat,t,;clims=(-5,5)); body_plot!(sim); # plot bodyz
            p2 = plot(data.t,data.Thrust,legend=false); xlims!(0,cycles * 2/St); ylabel!("Fd")
            p3 = plot(data.t,data.Side_Force,legend=false); xlims!(0,cycles * 2/St); ylabel!("Fs")
            p4 = plot(data.t,data.Moment,legend=false); xlims!(0,cycles * 2/St); ylabel!("M")
            p5 = plot(data.t,data.power,legend=false); ylabel!("P"); xlims!(0,cycles * 2/St)
            p6 = plot(data.t,data.ypos,legend=false,c = :orange); ylabel!("h"); xlims!(0,cycles * 2/St)
            p7 = plot(data.t,data.theta,legend=false, c = :orange); ylabel!("θ"); xlims!(0,cycles * 2/St); xlabel!("Time")

            plot(p1,p2,p3,p4,p5,p6,p7,layout = l,size = (800,800))
        else
            plot_vorticity(dat,t,;clims=(-5,5)); body_plot!(sim); # plot bodyz
        end

    end

    gif(anim, "foil.gif", fps=24)

    efficency = -1*(mean(data.Thrust))/(mean(data.power))

    CSV.write("L="*string(L)*"_"*"St="*string(St)*"_"*string(coeff)*".csv",  data, writeheader=true)

end


### example usage

run_singlefoilsim()


