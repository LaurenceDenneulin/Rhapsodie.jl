#
# grad_tools.jl
#
# Provide tools to load the instrumental parameters and the 
# calibrated data, and for the calculus of the RHAPSODIE data fidelity term. 
#
# ------------------------------------------------
#
# This file is part of Rhapsodie
#
#
# Copyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

#------------------------------------------------

const ImageInterpolator{T<:AbstractFloat, K<:Kernel{T}} = TwoDimensionalTransformInterpolator{T,K,K}
const MyKer = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)
#const MyKer = LinearInterpolators.RectangularSpline(Float64, LinearInterpolators.Flat)
#const MyKer = LinearInterpolators.LinearSpline(Float64, LinearInterpolators.Flat)

struct parameters_table{T<: AbstractFloat}
     cols::NTuple{3,Int64} #Imput size (reconstruction)
     rows::NTuple{2,Int64} #Output size (data)
     dataset_length::Int64
     Nframe::Int64
     Nrot::Int64
     Nangle::Int64
     v::Vector{NTuple{2, NTuple{3, T}}}
     indices::Array{Int64,2}
     center::Array{T,1}
     psf_center::NTuple{2,Array{T,1}}
     epsilon::Vector{NTuple{2,Array{T,1}}} 
     derotang::Vector{T}
end


struct FieldTransformOperator{T<:AbstractFloat, L<:Mapping, R<:Mapping} <: LinearMapping
    cols::NTuple{3,Int} 
    rows::NTuple{2,Int} 
    v_l::NTuple{3,T}
    v_r::NTuple{3,T}
    H_l::L              
    H_r::R      
end

struct data_table{T<:AbstractFloat,A<:FieldTransformOperator{T}} #autres types
    data::Array{T, 2};
    weights::Array{T, 2};
    H::A
end


const Trans_Table = Vector{NTuple{2,AffineTransform2D}}(); #Contains all the affine transform used 
const Parameters = parameters_table[];
get_par()::parameters_table = Parameters[1];
const dataset = data_table[];

const PSF_save = Vector{Array{Float64,2}}();
get_PSF()=PSF_save[1];
const EPSILON_save = Array{Float64,1}(undef,1);
function set_epsilon(epsilon)
    EPSILON_save[1]=epsilon
end
get_epsilon()::Float64=EPSILON_save[1];
const PRECOND_SAVE= Vector{Any}();
U()=PRECOND_SAVE[1];

const MASK_save = Vector{Array{Float64,2}}();
get_MASK()=MASK_save[1];

function load_parameters(size_data::NTuple{3,Int64}, 
                         Nframe::Int64,
                         Nrot::Int64, 
                         Nangle::Int64, 
                         center::Array{Float64,1}, 
                         psf_center::NTuple{2,Array{Float64,1}},
                         epsilon::Vector{NTuple{2,Array{Float64,1}}})
    sets_indices=Indices(size_data[3],Nrot,Nframe)
    sets_v=Set_Vi(sets_indices)
    bbox_output=(0,0)
    Id = AffineTransform2D{Float64}()
    for k=1:size_data[3]
        A_left=translate( epsilon[k][1][1], epsilon[k][1][2], Id)
        A_right=translate( epsilon[k][2][1], epsilon[k][2][2], Id) 
         
        out_left=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
        out_right=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
               
        xmax=max(bbox_output[1],out_left[1],out_right[1]);
        ymax=max(bbox_output[2],out_left[2],out_right[2]);
        
        bbox_output =(xmax,ymax);  
    end
    bbox_output = bbox_output .+ 4;
    push!(Parameters, parameters_table((bbox_output[1],bbox_output[2],3), 
                                       (size_data[1], size_data[2]), 
                                       size_data[3], 
                                       Nframe, 
                                       Nrot, 
                                       Nangle, 
                                       sets_v, 
                                       sets_indices, 
                                       center, 
                                       psf_center,
                                       epsilon, 
                                       Float64[])); 
    
        newcenter= (bbox_output .+1)./2;
        centerdiff=newcenter .- (center[1], center[2])
        for k=1:size_data[3]
        A_left=translate( epsilon[k][1][1] + centerdiff[1], epsilon[k][1][2] + centerdiff[2], Id)
        A_right=translate( epsilon[k][2][1] + centerdiff[1], epsilon[k][2][2] + centerdiff[2], Id)  
        push!(Trans_Table, (A_left, A_right))
    end
end



function load_parameters(size_data::NTuple{3,Int64}, 
                         Nframe::Int64,
                         Nrot::Int64, 
                         Nangle::Int64, 
                         center::Array{Float64,1}, 
                         psf_center::NTuple{2,Array{Float64,1}}, 
                         epsilon::Vector{NTuple{2,Array{Float64,1}}}, 
                         derotang::Vector{Float64})
    sets_indices=Indices(size_data[3],Nrot,Nframe)
    sets_v=Set_Vi(sets_indices)
    bbox_output=(0,0)
    Id = AffineTransform2D{Float64}()
    for k=1:size_data[3]
        A_left=TransRotate(Id, epsilon[k][1], derotang[k], center, center)
        A_right=TransRotate(Id, epsilon[k][2], derotang[k], center, center)   
        
        out_left=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
        out_right=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
                
        xmax=max(bbox_output[1],out_left[1],out_right[1]);
        ymax=max(bbox_output[2],out_left[2],out_right[2]);
        
        bbox_output =(xmax,ymax);                               
    end
    bbox_output = bbox_output .+ 4;
    push!(Parameters, parameters_table((bbox_output[1],bbox_output[2],3), 
                                       (size_data[1], size_data[2]), 
                                       size_data[3], 
                                       Nframe, 
                                       Nrot, 
                                       Nangle, 
                                       sets_v, 
                                       sets_indices, 
                                       center, 
                                       psf_center,
                                       epsilon, 
                                       derotang)); 
    
    newcenter= (bbox_output .+1)./2
    for k=1:size_data[3]
        A_left=inv(TransRotate(Id, epsilon[k][1], derotang[k], center, newcenter))
        A_right=inv(TransRotate(Id, epsilon[k][2], derotang[k], center, newcenter))   
        push!(Trans_Table, (A_left, A_right))
    end
end

function load_parameters(size_data::NTuple{3,Int64}, 
                         Nframe::Int64,
                         Nrot::Int64, 
                         Nangle::Int64, 
                         center::Array{Float64,1}, 
                         psf_center::NTuple{2,Array{Float64,1}},
                         epsilon::Vector{NTuple{2,Array{Float64,1}}},
                         v::Array{Float64,2})
    sets_indices=Indices(size_data[3],Nrot,Nframe)
    sets_v=Set_Vi(Nframe, size_data[3], v)
    bbox_output=(0,0)
    Id = AffineTransform2D{Float64}()
    for k=1:size_data[3]
        A_left=translate( epsilon[k][1][1], epsilon[k][1][2], Id)
        A_right=translate( epsilon[k][2][1], epsilon[k][2][2], Id) 
         
        out_left=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
        out_right=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
               
        xmax=max(bbox_output[1],out_left[1],out_right[1]);
        ymax=max(bbox_output[2],out_left[2],out_right[2]);
        
        bbox_output =(xmax,ymax);  
    end
    bbox_output = bbox_output .+ 4;
    push!(Parameters, parameters_table((bbox_output[1],bbox_output[2],3), 
                                       (size_data[1], size_data[2]), 
                                       size_data[3], 
                                       Nframe, 
                                       Nrot, 
                                       Nangle, 
                                       sets_v, 
                                       sets_indices, 
                                       center, 
                                       psf_center,
                                       epsilon, 
                                       Float64[])); 
    
        newcenter= (bbox_output .+1)./2;
        centerdiff=newcenter .- (center[1], center[2])
        for k=1:size_data[3]
        A_left=translate( epsilon[k][1][1] + centerdiff[1], epsilon[k][1][2] + centerdiff[2], Id)
        A_right=translate( epsilon[k][2][1] + centerdiff[1], epsilon[k][2][2] + centerdiff[2], Id)  
        push!(Trans_Table, (A_left, A_right))
    end
end



function load_parameters(size_data::NTuple{3,Int64}, 
                         Nframe::Int64,
                         Nrot::Int64, 
                         Nangle::Int64, 
                         center::Array{Float64,1}, 
                         psf_center::NTuple{2,Array{Float64,1}}, 
                         epsilon::Vector{NTuple{2,Array{Float64,1}}}, 
                         derotang::Vector{Float64},
                         v::Array{Float64,2})
    sets_indices=Indices(size_data[3],Nrot,Nframe)
    sets_v=Set_Vi(Nframe, size_data[3], v)
    bbox_output=(0,0)
    Id = AffineTransform2D{Float64}()
    for k=1:size_data[3]
        A_left=TransRotate(Id, epsilon[k][1], derotang[k], center, center)
        A_right=TransRotate(Id, epsilon[k][2], derotang[k], center, center)   
        
        out_left=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
        out_right=bbox_size((size_data[1], Int64(size_data[2]/2)), A_left)[1]
                
        xmax=max(bbox_output[1],out_left[1],out_right[1]);
        ymax=max(bbox_output[2],out_left[2],out_right[2]);
        
        bbox_output =(xmax,ymax);                               
    end
    bbox_output = bbox_output .+ 4;
    push!(Parameters, parameters_table((bbox_output[1],bbox_output[2],3), 
                                       (size_data[1], size_data[2]), 
                                       size_data[3], 
                                       Nframe, 
                                       Nrot, 
                                       Nangle, 
                                       sets_v, 
                                       sets_indices, 
                                       center, 
                                       psf_center,
                                       epsilon, 
                                       derotang)); 
    
    newcenter= (bbox_output .+1)./2
    for k=1:size_data[3]
        A_left=inv(TransRotate(Id, epsilon[k][1], derotang[k], center, newcenter))
        A_right=inv(TransRotate(Id, epsilon[k][2], derotang[k], center, newcenter))   
        push!(Trans_Table, (A_left, A_right))
    end
end




function Load_Data(name_data, name_weight, name_psf)
	data=read(FitsArray, name_data);
	weight=read(FitsArray, name_weight);
	psf=read(FitsArray, name_psf);
	A_l=set_fft_op((psf[1:end÷2,:]'[:,:]), get_par().psf_center[1]);
	A_r=set_fft_op((psf[end÷2+1:end,:]'[:,:]), get_par().psf_center[1]);
	ker= MyKer;
	SetCropOperator()
    for k=1:size(data)[3]   
        output_size=(get_par().rows[1], Int64(get_par().rows[2]/2));
        input_size= get_par().cols[1:2];
    	T1=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Trans_Table[k][1])
    	T2=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Trans_Table[k][2])
	    push!(dataset, data_table((data[:,:,k]')[:,:], 
	                              (weight[:,:,k]')[:,:], 
	                              FieldTransformOperator(get_par().cols, 
	                                                     get_par().rows, 
	                                                     get_par().v[k][1],
	                                                     get_par().v[k][2],
	                                                     T1,
	                                                     T2,
	                                                     A_l)));#,
	                                                     #A_r)));
	end
end

function Load_Data(name_data, name_weight)
	data=read(FitsArray, name_data);
	weight=read(FitsArray, name_weight);
	ker= MyKer;
	SetCropOperator()
    for k=1:size(data)[3]   
        output_size=(get_par().rows[1], Int64(get_par().rows[2]/2));
        input_size= get_par().cols[1:2];
    	T1=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Trans_Table[k][1])
    	T2=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Trans_Table[k][2])
    	promote
	    push!(dataset, data_table((data[:,:,k]')[:,:], 
	                              (weight[:,:,k]')[:,:], 
	                              FieldTransformOperator(get_par().cols, 
	                                                     get_par().rows, 
	                                                     get_par().v[k][1],
	                                                     get_par().v[k][2],
	                                                     T1,
	                                                     T2,
	                                                     LazyAlgebra.Id)));
	end
end



function Indices(S::Int64,Nrot::Int64,Nframe::Int64)
#S = le nombre total de frame gauche ou droite
#Nrot : Nombre de positions du derotateur
#Nframe : Nombre de frame par positions de la lame demi-onde
	
	Nangle=Int32.(S/(Nrot*Nframe*4)) #Nombre de rotations de la lame demi-onde
	INDICES=zeros(4,Nframe*Nangle*Nrot)
	for i=1:4
		ind=repeat(range(0,stop=4*Nframe*(Nrot*Nangle-1),length=Nrot*Nangle), inner=Nframe)+(Nframe*i .-mod.(range(1,stop=Nframe*Nangle*Nrot,length=Nframe*Nangle*Nrot),Nframe))
		INDICES[i,:]=ind
	end
	INDICES=round.(Int64, INDICES);
	return INDICES
end


function Set_Vi(Indices::AbstractArray{Int64,2}; alpha=[0, pi/4, pi/8, 3*pi/8], psi=[0,pi/2])
	J(a)=[cos(2*a) sin(2*a); sin(2*a) -cos(2*a)]
	P(a)=[cos(a), -sin(a)]
    v=Vector{NTuple{2, NTuple{3, Float64}}}(undef, length(Indices));
	for i=1:4
		v1=J(alpha[i])*P(psi[1])
		vn1=norm(v1)^2;
		v2=J(alpha[i])*P(psi[2])
		vn2=norm(v2)^2;
		for k=1:length(Indices[i,:])		
		v[Indices[i,k]] =((round(vn1/2, digits=4), round((abs(v1[1])^2-abs(v1[2])^2)/2, digits=4), round(real(v1[1]*v1[2]), digits=4)),
		                  (round(vn2/2, digits=4), round((abs(v2[1])^2-abs(v2[2])^2)/2, digits=4), round(real(v2[1]*v2[2]), digits=4)));
        end
	end
	return v
end

function Set_Vi(Nframe::Int, dataset_length::Int, mueller_instru::AbstractArray{Float64,2})
    @assert dataset_length == Nframe * size(mueller_instru)[1]
    
    v=Vector{NTuple{2, NTuple{3, Float64}}}(undef, dataset_length);
    for k=1:size(mueller_instru)[1]
        for l=1:Nframe
            v[Nframe*(k-1) + l] =((mueller_instru[k,1], mueller_instru[k,2], mueller_instru[k,3]),
                                  (mueller_instru[k,4], mueller_instru[k,5], mueller_instru[k,6]));
        end
    end
	return v
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
 	MapSize=get_par().cols[1:2];
	MapCenter=floor.(MapSize./2).+1
	MAP=zeros(MapSize);
	ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)

	Id = AffineTransform2D{Float64}()
	Recentering=translate(-(MapCenter[1]-PSFCenter[1]), -(MapCenter[2]-PSFCenter[2]), Id)

	LazyAlgebra.apply!(MAP, ker, Recentering, PSF);
	MAP./=sum(MAP);
	push!(PSF_save,MAP);
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
    dst[:, 1:(n÷2)  ] = R.H_l*z_l;
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
    #f = vdot(Float64, r,wr); #FIXME : mettre une issue !
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
    xymin_1=findfirst(get_Mask() != 0.)
    xymin_2=findfirst(get_Mask()' != 0.)
    xymax_1=findlast(get_Mask() != 0.)
    xymax_2=findlast(get_Mask()' != 0.)
    
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
        


