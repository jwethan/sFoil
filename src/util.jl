using WaterLily,LinearAlgebra

#### MULTIFOIL FUNCTIONS ####
function multipar∮nds(p::AbstractArray{T,N},df::AbstractArray{T},body::AbstractBody,t=0,L=32,s=32,nfin = 1) where {T,N}
    multiparnds!(df,body,t)
    for i in 1:N
        WaterLily.@loop df[I,i] = df[I,i]*p[I] over I ∈ inside(p)
    end

    spacing = L + s
    tiplocations = [nfin*spacing for nfin in range(0,nfin-1,nfin)] .+ 2L

    forces = []

    for i in eachindex(tiplocations)
        if tiplocations[i] == tiplocations[end]
            push!(forces,reshape(sum(df[round(Int,tiplocations[i])  : end,:,:],dims=1:N),N))
        else
            push!(forces,reshape(sum(df[round(Int,tiplocations[i]) : round(Int, tiplocations[i]+spacing),:,:],dims=1:N),N))
        end
    end

    return forces
end
multiparnds!(a,body,t=0) = apply!(a) do i,x
    d = ParametricBodies.sdf(body,x,t)
    n = ForwardDiff.gradient(y -> ParametricBodies.sdf(body,y,t), x)
    n[i]*WaterLily.kern(clamp(d,-1,1))
end

#### OLD SINGLE FOIL FUNCTIONS ####
function par∮nds(p::AbstractArray{T,N},df::AbstractArray{T},body::AbstractBody,t=0) where {T,N}
    parnds!(df,body,t)
    for i in 1:N
        WaterLily.@loop df[I,i] = df[I,i]*p[I] over I ∈ inside(p)
    end
    reshape(sum(df,dims=1:N),N) |> Array
end
parnds!(a,body,t=0) = apply!(a) do i,x
    d = ParametricBodies.sdf(body,x,t)
    n = ForwardDiff.gradient(y -> ParametricBodies.sdf(body,y,t), x)
    n[i]*WaterLily.kern(clamp(d,-1,1))
end

function par∮ndsmoment(p::AbstractArray{T,N},df::AbstractArray{T},body::AbstractBody,m,L,St,U,t=0) where {T,N}
    parndsmoment!(df,body,m,L,St,U,t)
    for i in 1:N
        WaterLily.@loop df[I,i] = df[I,i]*p[I] over I ∈ inside(p)
    end
    sum(reshape(sum(df,dims=1:N),N))
end
parndsmoment!(a,body,m,L,St,U,t=0) = apply!(a) do i,x
    d = ParametricBodies.sdf(body,x,t)
    n = ForwardDiff.gradient(y -> ParametricBodies.sdf(body,y,t), x)

    h₀ = L
    ω = Float64(π*St*U/h₀)

    dfromnose = x - SA[4L,.5f0m*L+h₀*sin(ω*t)]

    if i == 1
       -1*n[i]*WaterLily.kern(clamp(d,-1,1))*dfromnose[2]
    else
        n[i]*WaterLily.kern(clamp(d,-1,1))*dfromnose[1]
    end
end


#### 2.0 SINGLE FOIL FUNCTIONS ####
function my∮nds(p::AbstractArray{T,N},df::AbstractArray{T},body::AbstractBody,t=0) where {T,N}
    WaterLily.@loop df[I,:] .= p[I]*nds(body,loc(0,I,T),t) over I ∈ inside(p)
    [sum(@inbounds(df[inside(p),i])) for i ∈ 1:N] |> Array
end
@inline function nds(body::AbstractBody,x,t)
    d,n,_ = measure(body,x,t)
    n*WaterLily.kern(clamp(d,-1,1))
end

function moment∮nds(p::AbstractArray{T,N},df::AbstractArray{T},body::AbstractBody,noseloc::StaticArray,t=0) where {T,N}
    WaterLily.@loop df[I,:] .= p[I]*momentnds(body,loc(0,I,T),noseloc,t) over I ∈ inside(p) 
    [sum(@inbounds(df[inside(p),i])) for i ∈ 1:N] |> Array
end
@inline function momentnds(body::AbstractBody,x,noseloc::StaticArray,t)
    d,n,_ = measure(body,x,t)
    r = x - noseloc  #[dx,dy] ;distance vector from tip to x, used to grab perpendicular distances for moment calc
    dvec = reverse(r).*SA[-1.,1.] #[-dy,dx] 
    @. n*dvec*WaterLily.kern(clamp(d,-1,1)) #[nx,ny]*[-dy,dx]   #  ---> M = p(-nx*dy+ny*dx)
end
#Example usage
#force = -1*WaterLily.∮nds(sim.flow.p,sim.flow.V,sim.body,t);
#moment = -1*sum(moment∮nds(sim.flow.p,sim.flow.V,sim.body,SA[nNose*L,0.5f0m*L + h(t)],t));  

### 3.0 SINGLE FOIL FUNCTIONS ###
function foil∮nds(p::AbstractArray{T,N},df::AbstractArray{T},body::AbstractBody,xₒ::StaticArray,t=0) where {T,N}
    """
    Calculates the force and moment on a body due to a pressure field, using sim.flow.p and sim.body. Df is any matrix of the same sized used to hold calculations and xₒ is the location you
    want to take the moment about.
    """
    
    WaterLily.@loop df[I,:] .= p[I]*nds(body,loc(0,I,T),t) over I ∈ inside(p)
    force = [sum(@inbounds(df[inside(p),i])) for i ∈ 1:N] |> Array

    M = similar(p)
    WaterLily.@loop M[I] = p[I]*Mnds(body,loc(0,I,T),xₒ,t) over I ∈ inside(p) 
    moment = sum(@inbounds(M[inside(p)]))

    return force,moment
end
@inline function nds(body::AbstractBody,x,t=0)
    d,n,_ = measure(body,x,t)
    n*WaterLily.kern(clamp(d,-1,1))
end
@inline function Mnds(body::AbstractBody,x,xₒ,t=0)    
    d,n,_ = measure(body,x,t)
    r = x - xₒ 
    cross(r,n*WaterLily.kern(clamp(d,-1,1))) 
end

function CSTgen(s,t)
    zetaflag = 1; psi = (1-s)^2
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


function polyarea(x::Vector{Float64}, y::Vector{Float64}; flag = true)
    """
    polyarea(x::Vector{Float64}, y::Vector{Float64})
    
    Compute the area of a polygon defined by the vertices (x, y).

    If ρ = 1 then area = Mass and the second moment of area is equal to the moment of inertia (or second moment of mass).

    essesntially dA = dm, if its uniform density than ρdA = dM

    #set to automatically correct for a centroid shifted from the origin, if you want the values without the correction set flag to false. (returns parallelaxis theorm equiv.)
    """

    # need to able to center the centroid at the origin for covenience

    area = 0.; Ix = 0.; Iy = 0.; xcen = 0.; ycen = 0.
    N = length(x)

    for i in 1:N-1
        j = i + 1
        da = x[i]*y[j] - x[j]*y[i]
        area += da
        xcen += (x[i]+x[j])*da
        ycen += (y[i]+y[j])*da
    end

    xcen = xcen/(3. *area)
    ycen = ycen/(3. *area)

    if flag == true # correct for shifted centroid
        x = x .- xcen
        y = y .- ycen
    end

    for i in 1 : N-1
        j = i + 1
        da = x[i]*y[j] - x[j]*y[i]
        Ix +=  1/12. * (y[i]^2 + y[i]*y[j] + y[j]^2.) * da
        Iy +=  1/12. * (x[i]^2 + x[i]*x[j] + x[j]^2.) * da
    end

    I0 = Ix + Iy

    return 1/2. * area, Ix, Iy, I0, (xcen,ycen)

end

function foilMetrics(L::Int,coeff::AbstractVector,flag)

    """ 
    Wrapper for calculating the centroid,area and second moments of area of a foil defined by CST coefficients. 
    """

    function CSTgen(s,t)
        psi = (1-s)^2
        ncoeff = size(coeff)[1]-1
        Classfunc = (psi^.5 )*(1. - psi)^1.
        Shapefunc = 0; 
        for i in 1:ncoeff+1 
            K = factorial(ncoeff)/(factorial(i-1)*factorial(ncoeff-(i-1)))
            Shapefunc = Shapefunc + (coeff[i]*K)*psi^(i-1)*(1 - psi)^(ncoeff-(i-1))
        end
        zeta = Classfunc*Shapefunc
        return L*SA[psi,zeta]
    end

    x = []; y = []
    for s in 0:.1:1
        pnt = CSTgen(s,0)
        push!(x,pnt[1])
        push!(y,pnt[2])
    end

    x  = [x[1:end];reverse(x[1:end-1])] ; y  = [y[1:end];-1*reverse(y[1:end-1])]
    x = Float64.(x);y = Float64.(y)

    return polyarea(x,y;flag = flag)[1:4]

end


function lowpass(x::Float64, xfₒ::Float64, ωcut::Float64)
    """
    lowpass(x::Float64,xₒ::Float64, ωcut::Float64)
    
    A simple lowpass filter for a single value x where x is a new measurement, xfₒ is the previous filtered value of x, and ωcut is the cutoff frequency.
    """
  
    y = 1-cos(ωcut);      
    a = -y + (y^ 2. + 2. *y)^.5; #low pass filter ωᵥ is the cut off frequency
    xf =  a*x + (1. - a)*xfₒ
    return xf

end