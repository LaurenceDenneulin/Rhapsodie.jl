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

"""
    apply_rhapsodie('x0', 'A' , 'data', 'hyperparameters'; 'kwd')

where : 
    - x0 is an initial reconstruction
    - A is the operator of the convolution by the PSF
    - par is a set of hyper parameters given as follow : λ_I, λ_Q+U, ρ_I, ρ_Q+U
      where ρ is the ratio between the contribution λ and the threshold μ = λ/2ρ of the hypersmooth regularization;

"""

function apply_rhapsodie(x0::TPolarimetricMap, A::D, d::Array{Tdata_table,1}, par::Array{T,1}; mem=3, maxeval=50, maxiter=50, α::Real, xtol=(1e-3,1e-8), gtol=(1e-3,1e-8), ftol=(1e-3,1e-8)) where {T <: AbstractFloat, D <:Mapping}
    n1,n2 = size(x0)
    X0 = convert(Array{T,3},x0);
    μ=[hyperparameters(par[1], par[3]); 
       hyperparameters(par[2], par[4])];
    lower_born=vcreate(X0);
    vfill!(view(lower_born,:,:,1:2),0.0)
    vfill!(view(lower_born,:,:,3:4),-Inf)
   
    g=vcreate(X0);
    fg!(x,g) = apply_gradient!(TPolarimetricMap(x0.parameter_type, x), A, g, d, μ, α)
    x = vmlmb(fg!, X0, mem=mem, maxeval=maxeval, maxiter=maxiter, lower=lower_born, xtol=xtol,  gtol=gtol, ftol=ftol, verb=true);
    return TPolarimetricMap(x0.parameter_type, x)
end



function apply_gradient!(X::TPolarimetricMap, A::D, g::Array{T,3}, d::Array{Tdata_table,1}, μ::Array{hyperparameters{T},1}, α::Real) where {T <: AbstractFloat, D <:Mapping}

    n1, n2, n3 = size(g)
    @assert (n1,n2) == size(X)
    @assert n3 == 4
    

    if X.parameter_type == "intensities"
        error("Global reconstruction not implemented on Iu, Ip, and θ. 
               Only use 'stokes' or 'mixed' parameters.")
    end
    
    Ay = cat(X.I_star[:,:], A*X.I_disk[:,:], A*X.Q[:,:], A*X.U[:,:], dims=3)
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
    @assert size(g) == (n1,n2,4)
    @inbounds for i3 in 2:4
        g[:,:,i3].= A'*view(g,:,:,i3)[:,:];
    end

    if X.parameter_type == "mixed"
        @inbounds for i2 in 1:n2
            for i1 in 1:n1
                if X.Ip_disk[i1,i2] > 0
                    g[i1,i2,2] += (X.Q[i1,i2]*g[i1,i2,3] +
                                  X.U[i1,i2]*g[i1,i2,4])/X.Ip_disk[i1,i2]
                    # New try
                    # g[i1, i2, 3] += g[i1, i2, 2] * (X.Q[i1, i2] / X.Ip_disk[i1, i2])
                    # g[i1, i2, 4] += g[i1, i2, 2] * (X.U[i1, i2] / X.Ip_disk[i1, i2])
                end
            end
        end 
 	    #f+=cost!(μ[1][2] , μ[1][1], X.Iu[:,:], view(g,:,:,1), false);
 	    f+=apply_tikhonov!(X.Iu_star[:,:], view(g,:,:,1), μ[1].λ / (2 * μ[1].ρ));
        f+=apply_edge_preserving_smoothing!(cat(X.Iu_disk[:,:], X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:4), μ[2].λ, μ[2].ρ, α)
     elseif X.parameter_type == "stokes"
 	    f+=apply_tikhonov!(X.I_star[:,:], view(g,:,:,1), μ[1].λ / (2 * μ[1].ρ));
        f+=apply_edge_preserving_smoothing!(cat(X.I_disk[:,:], X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:4), μ[2].λ, μ[2].ρ, α)
 	    #f+=cost!(μ[1][2] , μ[1][1], X.I[:,:], view(g,:,:,1), false);    
     end
	#f+=cost!(μ[2][2] , μ[2][1], cat(X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:3), false);

	return f
end

function apply_tikhonov!(x::AbstractArray{T,2},
        g::AbstractArray{T,2},
        λ::Real) where {T <: AbstractFloat}
    m,n = size(x)                               
    f = zero(T);
    r = zero(T);
    x1 = zero(T);
    x2 = zero(T);

    for i=1:m-1
        for j=1:n-1
            x1= (x[i,j] - x[i+1,j])/2
            x2= (x[i,j] - x[i,j+1])/2
            r = x1^2 + x2^2;
            ## Cost functon ##
            f += λ * r / 2;
                ## Gradient in x ##
                g[i,j] += λ*(x1 + x2);
                g[i+1,j] -= λ*x1; 
                g[i,j+1] -= λ*x2; 
        end
    end
    return f
end

function apply_edge_preserving_smoothing!(x::AbstractArray{T,3},
                                   g::AbstractArray{T,3},
                                   λ::Real, 
                                   ρ::Real,
                                   α::Real) where {T <: AbstractFloat}
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
            xIu1 = (x[i,j,1] - x[i+1,j,1])/2;
            xIu2 = (x[i,j,1] - x[i,j+1,1])/2;
            xQ1= (x[i,j,2] - x[i+1,j,2])/2
            xQ2= (x[i,j,2] - x[i,j+1,2])/2
            xU1= (x[i,j,3] - x[i+1,j,3])/2
            xU2= (x[i,j,3] - x[i,j+1,3])/2
            
            ndx = α * (xIu1^2 + xIu2^2) + xQ1^2 + xQ2^2 + xU1^2 + xU2^2;
            r =  ndx + μ^2;
            ## Cost functon ##
            f += λ*(√r -  μ);
            if r>0
                ## Gradient in x ##
                ∂r=2*√r;
                
                g[i,j,1] += λ* α *(xIu1 + xIu2)/∂r;
                g[i+1,j,1] -= λ * α * xIu1/∂r; 
                g[i,j+1,1] -= λ * α * xIu2/∂r; 
                
                g[i,j,2] += λ*(xQ1 + xQ2)/∂r;
                g[i+1,j,2] -= λ*xQ1/∂r; 
                g[i,j+1,2] -= λ*xQ2/∂r; 
                
                g[i,j,3] += λ*(xU1 + xU2)/∂r;
                g[i+1,j,3] -= λ*xU1/∂r; 
                g[i,j+1,3] -= λ*xU2/∂r; 
                
             end
        end
    end
    
    return f
end