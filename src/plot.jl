using StaticArrays, Plots, DataFrames, Statistics, WaterLily

function get_omega!(dat,sim)
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u) * sim.L / sim.U
    copyto!(dat,sim.flow.σ[inside(sim.flow.σ)])
end

function get_p!(dat,sim)

    test = copy(sim.flow.p);
    test[sim.flow.μ₀[:,:,1] .==0] .=0

    copyto!(dat,test[inside(sim.flow.p)])
end

function plot_vorticity(ω,t; clims=()) 
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. ω=min(clims[2],max(clims[1],ω))
    else
        clims = (minimum(ω),maximum(ω))
    end
    Plots.contourf(ω',dpi=300,
    color=palette(:RdBu_11), clims=clims, linewidth=0, levels = 10,
    aspect_ratio=:equal, legend=false, border=:none, title = "t = "*string(round(t,digits=2)))
end

function plot_vorticity(ω; clims=()) 
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. ω=min(clims[2],max(clims[1],ω))
    else
        clims = (minimum(ω),maximum(ω))
    end
    Plots.contourf(ω',dpi=300,
    color=palette(:RdBu_11), clims=clims, linewidth=0, levels = 10,
    aspect_ratio=:equal, legend=false, border=:none)
end



function plot_pressure(p; clims=()) 
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. p=min(clims[2],max(clims[1],p))
    else
        clims = (minimum(p),maximum(p))
    end

    c = cgrad(:RdBu_11, rev = true)
    p[abs.(p).< 0.1] .= 0

    Plots.contourf(p',dpi=300,
    color=c, clims=clims, linewidth=0, levels = 101,
    aspect_ratio=:equal, legend=true, border=:none)
end


function body_plot(sim;levels=[0],color=:turbo,R=inside(sim.flow.p))
    WaterLily.measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim))
    Plots.contour(sim.flow.σ[R]';levels=levels,color, cbar=false, lw=1, aspect_ratio=:equal)
end

function body_plot!(sim;levels=[0],lines=:black,R=inside(sim.flow.p))
    WaterLily.measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim))
    #Plots.contourf(sim.flow.σ[R]';levels,lines)
    Plots.contour!(sim.flow.σ[R]';levels,lines)
end

function body_plotdark!(sim;levels=[0,-100],lines=:orange,R=inside(sim.flow.p))
    WaterLily.measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim))
    Plots.contourf!(sim.flow.σ[R]';levels,lines)
    #Plots.contour!(sim.flow.σ[R]';levels,lines)
end

macro makegif(gif::Bool, expr::Expr)
    if gif
        return esc(:(Plots.@gif $expr))
    else
        return esc(:($expr))
    end
end
