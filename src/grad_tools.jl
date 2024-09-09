#
# grad_tools.jl
#
# Provide tools for the calculus of the RHAPSODIE data fidelity term. 
#
# ------------------------------------------------
#
# This file is part of Rhapsodie
#
#
# Copyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

#------------------------------------------------


struct parameters_table{T<: AbstractFloat}
     cols::NTuple{3, Int64} #Input size (reconstruction)
     rows::NTuple{2, Int64} #Output size (data)
     dataset_length::Int64
     Nframe::Int64
     Nrot::Int64
     Nangle::Int64
     v::Vector{NTuple{2, NTuple{3, T}}}
     indices::Array{Int64, 2}
     center::Array{T, 1}
     psf_center::NTuple{2,Array{T, 1}}
     epsilon::Vector{NTuple{2,Array{T, 1}}} 
     derotang::Vector{T}
end


struct FieldTransformOperator{T<:AbstractFloat, L<:Mapping, R<:Mapping} <: LinearMapping
    cols::NTuple{3, Int} 
    rows::NTuple{2, Int} 
    v_l::NTuple{3, T}
    v_r::NTuple{3, T}
    H_l::L              
    H_r::R
end

struct data_table{T<:AbstractFloat, A<:FieldTransformOperator{T}} # Other types
    data::Array{T, 2};
    weights::Array{T, 2};
    H::A
end


function reset_instrument(V::Vector{NTuple{2, NTuple{3, T}}})  where {T <: AbstractFloat}
    push!(Parameters, parameters_table(get_par().cols, 
                                       get_par().rows, 
                                       get_par().dataset_length, 
                                       get_par().Nframe, 
                                       get_par().Nrot, 
                                       get_par().Nangle, 
                                       V, 
                                       get_par().indices, 
                                       get_par().center, 
                                       get_par().psf_center,
                                       get_par().epsilon, 
                                       get_par().derotang)); 

    popfirst!(Parameters);
end



function TransRotate(A::AffineTransform2D{Float64}, EPSILON, ANGLE, CENTER, NEWCENTER)
	return translate(rotate(translate(CENTER[1]-EPSILON[1],CENTER[2]-EPSILON[2], A), ANGLE),-NEWCENTER[1], -NEWCENTER[2])
end

function bbox_size(inp_dims::NTuple{2,Integer},
                  A::AffineTransform2D{Float64})
                    
	xmin=typemax(Float64);
	ymin=typemax(Float64);
	xmax=typemin(Float64);
	ymax=typemin(Float64);
	
	xmin_int=typemin(Int64);
	ymin_int=typemin(Int64);
	xmax_int=typemin(Int64);
	ymax_int=typemin(Int64);
	
	width=typemin(Float64);
	height=typemin(Float64);
	
	Ind=[repeat(1:inp_dims[1], inner=inp_dims[2]) repeat(1:inp_dims[2], outer=inp_dims[1])]
	for it=1:(inp_dims[1]*inp_dims[2])
		(x,y)=A(Ind[it,1], Ind[it,2]);
		xmin=min(x, xmin);
		ymin=min(y, ymin);
		xmax=max(x, xmax);
		ymax=max(y, ymax);
	end
	xmin_int=floor(Int64,xmin);
	ymin_int=floor(Int64,ymin);
	xmax_int=ceil(Int64,xmax);
	ymax_int=ceil(Int64,ymax);

	width=xmax_int-xmin_int+1+50;
	height=ymax_int-ymin_int+1+50;
    out_dims=(width, height);

	return(out_dims, (xmin_int, xmax_int, ymin_int, ymax_int))
end

function set_fft_op(PSF::AbstractArray{T,2}, PSFCenter::AbstractArray{T,1}) where {T <: AbstractFloat}
 	MapSize=get_par().cols[1:2]
	MapCenter=floor.(MapSize./2).+1
	MAP=zeros(MapSize)
	ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)

	Id = AffineTransform2D{Float64}()
	Recentering=translate(-(MapCenter[1]-PSFCenter[1]), -(MapCenter[2]-PSFCenter[2]), Id)

	LazyAlgebra.apply!(MAP, ker, Recentering, PSF)
	MAP./=sum(MAP)
	F=FFTOperator(MAP)
	FFT=F\Diag(F*ifftshift(MAP)) .*F;
	
	return FFT
end

function FieldTransformOperator(cols::NTuple{3,Int64},
                                rows::NTuple{2,Int64},
                                v_l::NTuple{3,T},
                                v_r::NTuple{3,T},
                                T_l::TwoDimensionalTransformInterpolator{T},
                                T_r::TwoDimensionalTransformInterpolator{T},
                                A::M) where {T <: AbstractFloat, M <: Mapping}
    FieldTransformOperator(cols, rows, v_l, v_r, T_l*A, T_r*A)
end

function vcreate(::Type{LazyAlgebra.Direct}, A::FieldTransformOperator{T},
                 x::AbstractArray{T,3}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.cols
    Array{T,2}(undef, A.rows)
end

function vcreate(::Type{LazyAlgebra.Adjoint}, A::FieldTransformOperator{T},
                 x::AbstractArray{T,2}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.rows
    Array{T,3}(undef, A.cols)
end


function apply!(α::Real,
                ::Type{LazyAlgebra.Direct},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,3},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,2}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.cols
    @assert size(dst) == R.rows
    n = R.rows[2]
    @assert iseven(n)

    z_l = zeros(R.cols[1:2]);
    z_r = zeros(R.cols[1:2]);
    for i=1:length(R.v_l)
        # FIXME: use Generalized Matrix-Vector Multiplication of LazyAlgebra here.
        x_i = view(src,:,:,i)
	z_l .+= R.v_l[i]*x_i;
	z_r .+= R.v_r[i]*x_i;
    end
    dst[:, 1:(n÷2)] = R.H_l*z_l;
    dst[:, (n÷2)+1:n] = R.H_r*z_r;
    return dst
end

function apply!(α::Real,
                ::Type{LazyAlgebra.Adjoint},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,2},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,3}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.rows
    @assert size(dst) == R.cols
    n = R.rows[2]
    @assert iseven(n)
    y_l = R.H_l'*view(src, :, 1:(n÷2))
    y_r = R.H_r'*view(src, :, (n÷2)+1:n)
    for i=1:length(R.v_l)
        #dst[:,:,i] = R.v_l[i]*y_l + R.v_r[i]*y_r;
        # FIXME: The following should be equivalent and much faster:
         vcombine!(view(dst,:,:,i), R.v_l[i], y_l, R.v_r[i], y_r)
    end
    return dst;
end


#=
function apply!(α::Real,
                ::Type{LazyAlgebra.Direct},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,3},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,2}) where {T<:AbstractFloat}
                
    @assert size(src) == R.cols
    @assert size(dst) == R.rows
    @assert β==0 && α==1
    vfill!(dst,0);
    
    dst_left=zeros(R.cols[1:2]);
    dst_right=zeros(R.cols[1:2]);
    
    for i=1:length(R.v_l)                                                     
        dst_left .+=(R.v_l[i]*src[:,:,i]); 
        dst_right .+=(R.v_r[i]*src[:,:,i]);   
    end   
    dst .+= [R.T_l*R.A_l*dst_left R.T_r*R.A_r*dst_right];
    return dst                                                           
end

function apply!(α::Real,
                ::Type{LazyAlgebra.Adjoint},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,2},              
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,3}) where {T<:AbstractFloat}
    @assert β==0 && α==1                                                                                  
    @assert size(src) == R.rows
    @assert size(dst) == R.cols 
    for i=1:length(R.v_l)                                                         
        dst_left=R.A_l'*R.T_l'*(R.v_l[i]*src[:,1:end÷2]);
        dst_right=R.A_r'*R.T_r'*(R.v_r[i]*src[:,end÷2+1:end]);
        dst[:,:,i] = dst_left + dst_right;
    end
    return dst;
end
=#
      
function mydot(x::AbstractArray{Tx,N},
               y::AbstractArray{Ty,N}) where {Tx<:AbstractFloat,Ty<:AbstractFloat,N}
    @assert axes(x) == axes(y)
    s = zero(promote_type(Tx,Ty))
    @inbounds @simd for i in eachindex(x,y)
       s += x[i]*y[i]
    end
    return s
end
             
function fg!(x::AbstractArray{T,3},g::AbstractArray{T,3}) where {T<:AbstractFloat}
    local f::Float64 = 0.0;
    vfill!(g,0.0)
    for k=1:length(dataset)
        if sum(dataset[k].weights) !=0
            f += fg!(x, g, dataset[k])
        end
    end
    #crop!(g)
    return f
end       

        
function fg!(x::AbstractArray{T,3},g::AbstractArray{T,3}, d::data_table{T}) where {T<:AbstractFloat}
    r = d.H*x - d.data;
    wr = d.weights .* r;
    g .+= d.H'*wr;
    #f = vdot(Float64, r,wr);
    f = mydot(r,wr);
    #crop!(g)
    return Float64(f)/2
end       
        
      

  
#--------------------------------------------UTILS------------------------------------------------
      
function SetCropOperator()
    DSIZE=get_par().rows[1]
    MASK=ones(get_par().cols[1:2]);
    for k=1:length(Trans_Table)
    X1=Trans_Table[k][1](1,1);
    X2=Trans_Table[k][1](1,DSIZE);
    X3=Trans_Table[k][1](DSIZE, DSIZE);
    X4=Trans_Table[k][1](DSIZE, 1);
    Mask1=zeros(get_par().cols[1:2]);

    for i=1:get_par().cols[1]
	    for j=1:get_par().cols[2]	
		    if ((i-X1[1])*(X2[1]-X1[1]) +(j-X1[2])*(X2[2]-X1[2]) >0)&&((i-X2[1])*(X3[1]-X2[1]) +(j-X2[2])*(X3[2]-X2[2]) >0)&&((i-X3[1])*(X4[1]-X3[1]) +(j-X3[2])*(X4[2]-X3[2]) >0) &&((i-X4[1])*(X1[1]-X4[1]) +(j-X4[2])*(X1[2]-X4[2]) >0)	
		    Mask1[i,j,:] .=1;
		    end
	    end
    end
    X1=Trans_Table[k][2](1,1);
    X2=Trans_Table[k][2](1,DSIZE);
    X3=Trans_Table[k][2](DSIZE, DSIZE);
    X4=Trans_Table[k][2](DSIZE, 1);
    Mask2=zeros(get_par().cols[1:2]);

    for i=1:get_par().cols[1]
	    for j=1:get_par().cols[2]	
		    if ((i-X1[1])*(X2[1]-X1[1]) +(j-X1[2])*(X2[2]-X1[2]) >0)&&((i-X2[1])*(X3[1]-X2[1]) +(j-X2[2])*(X3[2]-X2[2]) >0)&&((i-X3[1])*(X4[1]-X3[1]) +(j-X3[2])*(X4[2]-X3[2]) >0) &&((i-X4[1])*(X1[1]-X4[1]) +(j-X4[2])*(X1[2]-X4[2]) >0)	
		    Mask2[i,j,:] .=1;
		    end
	    end
    end
    MASK .*= Mask1.*Mask2
    end
    #write(FitsFile, "MASK.fits", MASK, overwrite=true)
    push!(MASK_save, MASK)
end     

        
function crop(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    #@assert size(X) .==   get_par().cols
    xymin_1=findfirst(get_MASK() .!= 0.)
    xymin_2=findfirst(get_MASK()' .!= 0.)
    xymax_1=findlast(get_MASK() .!= 0.)
    xymax_2=findlast(get_MASK()' .!= 0.)
    
    xmin=min(xymin_1[1], xymin_2[2]);
    ymin=min(xymin_1[2], xymin_2[1]);
    xmax=max(xymax_1[1], xymax_2[2]);
    ymax=max(xymax_1[2], xymax_2[1]);
    
    Y=X[xmin:xmax, ymin:ymax];
    Y[.!isfinite.(Y)].=0    
    return Y     
end        

function crop!(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    #@assert size(X) .==   get_par().cols 
    X[.!isfinite.(X)].=0   
    X .*= get_MASK()        
end        
        
function crop(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,3}}
    #@assert size(X) .==   get_par().cols
    Y=copy(X);
    return crop!(Y)    
end        

function crop!(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,3}}
    #@assert size(X) .==   get_par().cols
    for k in size(X)[3] 
        crop!(view(X,:,:,k))
    end     
end        
   
                
function crop(X::PolarimetricMap{T}) where {T<:AbstractFloat}  
    #@assert size(X) .==   get_par().cols
    return PolarimetricMap(X.parameter_type,
                           crop(view(X.I,:,:)),
                           crop(view(X.Q,:,:)),
                           crop(view(X.U,:,:)),
                           crop(view(X.Iu,:,:)),
                           crop(view(X.Ip,:,:)),        
                           crop(view(X.θ,:,:)))   
end        

function crop!(X::PolarimetricMap{T})  where {T<:AbstractFloat}
        crop!(view(X.I,:,:));
        crop!(view(X.Q,:,:));
        crop!(view(X.U,:,:));
        crop!(view(X.Iu,:,:));
        crop!(view(X.Ip,:,:));        
        crop!(view(X.θ,:,:));
end        
        
function pad(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    X_size = size(X);
    center_diff = X_size./2 .- get_par().cols[1:2]./2;
    ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)    
    Id = AffineTransform2D{Float64}()
    center_change = translate(center_diff[1],center_diff[2], Id)   
    PAD=TwoDimensionalTransformInterpolator(get_par().cols[1:2], X_size, ker, ker, center_change)
    return PAD*X;
end         
        
function pad(X::PolarimetricMap{T}) where {T<:AbstractFloat}  
    #@assert size(X) .==   get_par().cols
    return PolarimetricMap(X.parameter_type,
                           pad(view(X.I,:,:)),
                           pad(view(X.Q,:,:)),
                           pad(view(X.U,:,:)),
                           pad(view(X.Iu,:,:)),
                           pad(view(X.Ip,:,:)),        
                           pad(view(X.θ,:,:)))   
end