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


function apply_gradient!(X::PolarimetricMap, A::D, g::Array{T,3}, d::Array{data_table,1}, μ::Array{Tuple{H,T},1}) where {T <: AbstractFloat, H<:AbstractCost, D <:Mapping}

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
 	    f+=cost!(μ[1][2] , μ[1][1], X.Iu[:,:], view(g,:,:,1), false);
     elseif X.parameter_type == "stokes"
 	    f+=cost!(μ[1][2] , μ[1][1], X.I[:,:], view(g,:,:,1), false);    
     end

	f+=cost!(μ[2][2] , μ[2][1], cat(X.Q[:,:], X.U[:,:], dims=3), view(g,:,:,2:3), false);

	return f
end
   
function apply_rhapsodie(x0::PolarimetricMap, A::D, d::Array{data_table,1}, par::Array{T,1}; mem=3, maxeval=50, maxiter=50, xtol=(1e-3,1e-8), gtol=(1e-3,1e-8), ftol=(1e-3,1e-8)) where {T <: AbstractFloat, D <:Mapping}
    n1,n2 = size(x0)
    X0 = convert(Array{T,3},x0);
    μ=[(HyperbolicEdgePreserving(10. ^par[3], (0.5,0.5)),10. ^par[1]);
       (HyperbolicEdgePreserving(10. ^par[4],(0.5,0.5,0.)),10. ^par[2])];
       
    lower_born=vcreate(X0);
    vfill!(view(lower_born,:,:,1),0.0)
    vfill!(view(lower_born,:,:,2:3),-Inf)
   
    g=vcreate(X0);
    fg!(x,g)=apply_gradient!(PolarimetricMap(x0.parameter_type,x), A, g, d, μ)
    x = vmlmb(fg!, X0, mem=mem, maxeval=maxeval, maxiter=maxiter, lower=lower_born, xtol=xtol,  gtol=gtol, ftol=ftol);
    return PolarimetricMap(x0.parameter_type,x)
end

      
