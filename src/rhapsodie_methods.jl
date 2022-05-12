#
# grad_tools.jl
#
# Provide the RHAPSODIE methods to reconstruct the polarimetric 
# parameters from a dataset using:
#
# - VMLM-B (when "mixed" polarimetric parameters are used)
# [Thiébaut, 2002] Thiébaut, E. (2002). Optimization issues in 
# blind deconvolution algorithms. In Astronomical Data Analysis II,
# volume 4847, pages 174–183. International Society for Optics and Photonics.
#
# - Forward-Backward with backtracking (when "stokes" parameters are used)
# [Beck and Teboulle, 2009] Beck, A. and Teboulle, M. (2009). 
# A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems. 
# SIAM J. Imaging Sci., 2(1) :183–202. (TODO: clean implementation)
#
#
# ------------------------------------------------
#
# This file is part of Rhapsodie
#
#
# Copyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

#------------------------------------------------
struct hyperparameters{T<: AbstractFloat}
    λ::T # Weight of the regularization
    ρ::T # Ratio between the weight and the threshold μ = λ/(2*ρ);
end

```
    apply_rhapsodie(x0::PolarimetricMap, A::D, d::Array{data_table,1}, par::Array{T,1}; mem=3, maxeval=50, maxiter=50, xtol=(1e-3,1e-8), gtol=(1e-3,1e-8), ftol=(1e-3,1e-8))

where : 
    - x0 is an initial reconstruction
    - A is the operator of the convolution by the PSF
    - par is a set of hyper parameters given as follow : [λ_I, λ_Q+U, ρ_I, ρ_Q+U]
      where ρ is the ratio between the contribution λ and the threshold μ = λ/2ρ of the hypersmooth regularization;

```

function apply_rhapsodie(x0::PolarimetricMap, A::D, d::Array{data_table,1}, par::Array{T,1}; mem=3, maxeval=50, maxiter=50, xtol=(1e-3,1e-8), gtol=(1e-3,1e-8), ftol=(1e-3,1e-8)) where {T <: AbstractFloat, D <:Mapping}
    #par =  

    n1,n2 = size(x0)
    X0 = convert(Array{T,3},x0);
    μ=[hyperparameters(par[1], par[3]); 
       hyperparameters(par[2], par[4])];
       #[(HyperbolicEdgePreserving(10. ^par[3], (0.5,0.5)),10. ^par[1]);
       #(HyperbolicEdgePreserving(10. ^par[4],(0.5,0.5,0.)),10. ^par[2])];
       
    lower_born=vcreate(X0);
    vfill!(view(lower_born,:,:,1),0.0)
    vfill!(view(lower_born,:,:,2:3),-Inf)
   
    g=vcreate(X0);
    fg!(x,g)=apply_gradient!(PolarimetricMap(x0.parameter_type,x), A, g, d, μ)
    x = vmlmb(fg!, X0, mem=mem, maxeval=maxeval, maxiter=maxiter, lower=lower_born, xtol=xtol,  gtol=gtol, ftol=ftol);
    return PolarimetricMap(x0.parameter_type,x)
end



function apply_gradient!(X::PolarimetricMap, A::D, g::Array{T,3}, d::Array{data_table,1}, μ::Array{hyperparameters{T},1}) where {T <: AbstractFloat, D <:Mapping}

    n1, n2, n3 = size(g)
    @assert (n1,n2) == size(X)
    @assert n3 == 3
    if X.parameter_type == "intensities"
        error("Global reconstruction not implemented on Iu, Ip, and θ. 
               Only use 'stokes' or 'mixed' parameters.")
    end
    
    Ay = cat(A*X.I[:,:], A*X.Q[:,:], A*X.U[:,:], dims=3)
    # Compute data fidelity term and gradient. (As gradient is initially set to
    # zero, we can recycle it between x and y.)
    @assert size(g) == size(Ay)
    vfill!(g, 0)
    local f::Float64 = 0.0;
    for k = 1:length(d)
        f += fg!(Ay, g, d[k])
    end

    # Convert gradient w.r.t. y into gradient w.r.t. x.  Nothing has to be done
    # for the 2nd and 3rd fields (Q and U) or if Ip = 0.
    @assert size(g) == (n1,n2,3)
    @inbounds for i3 in 1:3
        g[:,:,i3].= A'*view(g,:,:,i3)[:,:];
    end

    if X.parameter_type == "mixed"
        @inbounds for i2 in 1:n2
            for i1 in 1:n1
                if X.Ip[i1,i2] > 0
                    g[i1,i2,1] += (X.Q[i1,i2]*g[i1,i2,2] +
                                   X.U[i1,i2]*g[i1,i2,3])/X.Ip[i1,i2]
                end
            end
        end 
 	    #f+=cost!(μ[1][2] , μ[1][1], X.Iu[:,:], view(g,:,:,1), false);
 	    f+=apply_edge_preserving_smoothing!(X.Iu[:,:], view(g,:,:,1), μ[1].λ, μ[1].ρ);
     elseif X.parameter_type == "stokes"
 	    f+=apply_edge_preserving_smoothing!(X.I[:,:], view(g,:,:,1), μ[1].λ, μ[1].ρ)
 	    #f+=cost!(μ[1][2] , μ[1][1], X.I[:,:], view(g,:,:,1), false);    
     end
    f+=apply_edge_preserving_smoothing!(cat(X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:3), μ[2].λ, μ[2].ρ)
	#f+=cost!(μ[2][2] , μ[2][1], cat(X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:3), false);

	return f
end
   

function apply_edge_preserving_smoothing!(x::M,
                                   g::M,
                                   λ::Real, 
                                   ρ::Real) where {T <: AbstractFloat, M<:AbstractArray{T,2}}
    m,n = size(x)                               
    f = zero(T);
    r = zero(T);
    μ = λ/(2*ρ);
    x1 = zero(T);
    x2 = zero(T);
        
    for i=1:m-1
        for j=1:n-1
            x1= (x[i,j] - x[i+1,j])/2
            x2= (x[i,j] - x[i,j+1])/2
            
            ndx = x1^2 + x2^2;
            r =  ndx + μ^2;
            ## Cost functon ##
            f += λ*(√r -  μ);
            if r>0
                ## Gradient in x ##
                ∂r=2*√r;
                g[i,j] .+= λ*(x1 + x2)/∂r;
                g[i+1,j] .-= λ*x1/∂r; 
                g[i,j+1] .-= λ*x2/∂r; 
                
             end
        end
    end
    
    return f
end

function apply_edge_preserving_smoothing!(x::M,
                                   g::M,
                                   λ::Real, 
                                   ρ::Real) where {T <: AbstractFloat, M<:AbstractArray{T,3}}
    m,n = size(x)                               
    f = zero(T);
    r = zero(T);
    μ = λ/(2*ρ);
    xQ1 = zero(T);
    xQ2 = zero(T);
    xU1 = zero(T);
    xU2 = zero(T);
        
    for i=1:m-1
        for j=1:n-1
            xQ1= (x[i,j,1] - x[i+1,j,1])/2
            xQ2= (x[i,j,1] - x[i,j+1,1])/2
            xU1= (x[i,j,2] - x[i+1,j,2])/2
            xU2= (x[i,j,2] - x[i,j+1,2])/2
            
            ndx = xQ1^2 + xQ2^2 + xU1^2 + xU2^2;
            r =  ndx + μ^2;
            ## Cost functon ##
            f += λ*(√r -  μ);
            if r>0
                ## Gradient in x ##
                ∂r=2*√r;
                g[i,j,1] .+= λ*(xQ1 + xQ2)/∂r;
                g[i+1,j,1] .-= λ*xQ1/∂r; 
                g[i,j+1,1] .-= λ*xQ2/∂r; 
                
                g[i,j,2] .+= λ*(xU1 + xU2)/∂r;
                g[i+1,j,2] .-= λ*xU1/∂r; 
                g[i,j+1,2] .-= λ*xU2/∂r; 
                
             end
        end
    end
    
    return f
end
         
         
function apply_edge_preserving_smoothing!(x::AbstractArray{<:AbstractFloat},
                                   g::AbstractArray{<:AbstractFloat}, 
                                   hx::AbstractArray{<:AbstractFloat}, 
                                   hλ::AbstractArray{<:AbstractFloat}, 
                                   λ::Real, 
                                   ρ::Real)
    m,n = size(x)                               
    f = zero(T);
    r = zero(T);
    ndx= zero(T);
    ∂r= zero(T);
    μ = λ/(2*ρ);
    ∂μ = 2* μ^2; #is actually λ * ∂(μ²)/∂λ
    x1 = zero(T);
    x2 = zero(T);
        
    for i=1:m-1
        for j=1:n-1
            x1= (x[i,j] - x[i+1,j])/2
            x2= (x[i,j] - x[i,j+1])/2
            
            ndx = x1^2 + x2^2;
            r =  ndx + μ^2;
            ## Cost functon ##
            f += λ*(√r -  μ);
            if r>0
                ## Gradient in x ##
                ∂r=2*√r;
                g[i,j] += λ*(x1 + x2)/∂r;
                g[i+1,j] -= λ*x1/∂r; 
                g[i,j+1] -= λ*x2/∂r; 
                
                ## Hessian in x ##
                ∂r *= 2*r;
                
                
                hx[(j-1)*m+i,(j-1)*m+i] +=  λ*(2*μ^2 + (x1-x2)^2)/∂r; #∂g_{i,j}/∂x_{i,j} 
                
                hx[j*m+i, j*m+i] +=  λ*(μ^2 + x1^2)/∂r; #∂g_{i+1,j}/∂x_{i+1,j} 
                
                hx[(j-1)*m+i+1,(j-1)*m+i+1] +=  λ*(μ^2 + x2^2)/∂r; #∂g_{i,j+1}/∂x_{i,j+1} 
                
                hx[j*m+i,(j-1)*m+i] +=  λ*(-μ^2 - x1^2 +x1*x2)/∂r;
                hx[(j-1)*m+i+1,(j-1)*m+i] +=  λ*(-μ^2 - x2^2 +x1*x2)/∂r;
                
                hx[(j-1)*m+i,j*m+i] +=  λ*(-μ^2 - x1^2 +x1*x2)/∂r;                
                hx[(j-1)*m+i,(j-1)*m+i+1] +=  λ*(-μ^2 - x2^2 +x1*x2)/∂r;
                
                hx[j*m+i,(j-1)*m+i+1] +=  λ*(-x1*x2)/∂r;
                hx[(j-1)*m+i+1,j*m+i] +=  λ*(-x1*x2)/∂r;

                ## Hessian in λ ##
                
                ∂ndx = ndx + ∂μ;
                ∂r /=∂ndx*2;
                
                hλ[i,j] +=  (x1 + x2)/∂r;
                hλ[i+1,j] -= x1/∂r ;
                hλ[i,j+1] -= x2/∂r ;
             end
        end
    end
    
    return f
end

