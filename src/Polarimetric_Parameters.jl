#
# Polarimetric_Parameters.jl
#
# Provide the type polarimetric parameter.
#
# ------------------------------------------------
#
# This file is part of RHAPSODIE
#
#
# Copyright (c) 2017-2021 Laurence Denneulin


#------------------------------------------------
# Structure definition
     
    #For pixles
    struct PolarimetricPixel{T<: AbstractFloat} 
        parameter_type::AbstractString   #either "stokes", "intensities" or "mixed"
        I::T
        Q::T
        U::T
        Iu::T
        Ip::T
        θ::T
    end
    
   #For map of pixels
    struct PolarimetricMap{T<: AbstractFloat} 
        parameter_type::AbstractString   #either "stokes", "intensities" or "mixed"
        I::Array{T,2};
        Q::Array{T,2};
        U::Array{T,2};
        Iu::Array{T,2};
        Ip::Array{T,2};
        θ::Array{T,2};
    end
    
    PolarimetricMap(parameter_type::AbstractString,        
                    I::AbstractArray{T,2},
                    Q::AbstractArray{T,2},
                    U::AbstractArray{T,2},
                    Iu::AbstractArray{T,2},
                    Ip::AbstractArray{T,2},
                    θ::AbstractArray{T,2}) where {T<:AbstractFloat} =  PolarimetricMap(parameter_type,        
                                                         convert(Array{T},I),
                                                         convert(Array{T},Q),
                                                         convert(Array{T},U),
                                                         convert(Array{T},I),
                                                         convert(Array{T},Ip),
                                                         convert(Array{T},θ))

                                     
#------------------------------------------------
# Constructors 
"""
    PolarimetricPixel(parameter_type, x) -> PolarimetricPixel
    PolarimetricMap(parameter_type, x) -> PolarimetricMap
    
create an object of type PolarimetricParameter from either:
    - Parameters I, Q, U (i.e. parameter_type = 'stokes')
    - Parameters Iu, Ip and θ (i.e. parameter_type = 'intensities')
    - Parameters Iu, Q, U (i.e. parameter_type = 'mixed')
    
Each parameter can be called from the structur. For exemple with a 
construction from Stokes parameters S=(I,Q,U):

    using RHAPSODIE
    X = PolarimetricParameter(S, 'stokes');
    X.I #yields the Stokes parameter I
    X.Ip #yields the polarized intensity Ip
    X[1,1] #yields a PolarimetricPix at the CartesianIndex (1,1);
    X[1,1].I #yields the Stokes parameter I at the CartesianIndex (1,1); 

"""    
    function PolarimetricPixel(parameter_type::AbstractString, 
                             x1::T, 
                             x2::T, 
                             x3::T) where {T<:AbstractFloat}
        if parameter_type == "stokes"
            I = x1            # Stokes parameter I (total light intensity)
            Q = x2            # Stokes parameter Q
            U = x3            # Stokes parameter U
            Ip = sqrt(Q^2 + U^2)    # intensity of linearly polarized light
            Iu = I-Ip    # intensity of unpolarized light
            θ = atan(U, Q)/2 # angle of linearly polarized light
        elseif parameter_type == "intensities"
            Iu = x1          # intensity of unpolarized light
            Ip = x2          # intensity of linearly polarized light
            θ = x3       # angle of linearly polarized light
            I = Iu + Ip   # Stokes parameter I (total light intensity)
            Q = Ip *cos(2*θ);    # Stokes parameter Q
            U = Ip *sin(2*θ);    # Stokes parameter U
        elseif parameter_type == "mixed" 
            Iu = x1          # intensity of unpolarized light   
            Q = x2            # Stokes parameter Q
            U = x3            # Stokes parameter U  
            Ip = sqrt(Q^2 + U^2)    # intensity of linearly polarized light
            I = Iu + Ip     # Stokes parameter I (total light intensity
            θ = atan(U, Q)/2 # angle of linearly polarized light
        else
            error("unkown type, only known types : stokes, intensities and mixed")
        end
        PolarimetricPixel(parameter_type,I,Q,U,Iu,Ip,θ)
    end

    function PolarimetricPixel(parameter_type::AbstractString, 
                             x::N) where {N<:AbstractArray{Float64,1}}
        @assert length(x) == 3 
        PolarimetricPixel(parameter_type, 
                        x[1],
                        x[2],
                        x[3]);
    end 
    
    function PolarimetricMap(x::Array{PolarimetricPixel,2}) 
        T = Float64;
        n1, n2 = size(x)
        par_type=x[1].parameter_type;          
        I = Array{T}(undef, n1, n2)    # Stokes parameter I (total light intensity)
        Q = Array{T}(undef, n1, n2)    # Stokes parameter Q
        U = Array{T}(undef, n1, n2)    # Stokes parameter U
        Iu = Array{T}(undef, n1, n2)    # intensity of unpolarized light
        Ip = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
        θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1   
                pix= x[i1,i2];
                if par_type == pix.parameter_type
                   I[i1,i2] = pix.I;
                   Q[i1,i2] = pix.Q;
                   U[i1,i2] = pix.U;
                   Iu[i1,i2] = pix.Iu;
                   Ip[i1,i2] = pix.Ip;
                   θ[i1,i2] = pix.θ;
                else
                   error("polarimetric types must match")
                end
            end 
        end
        PolarimetricMap(par_type,I,Q,U,Iu,Ip,θ)
    end
 
    function PolarimetricMap(parameter_type::AbstractString, 
                             x1::A, 
                             x2::A, 
                             x3::A) where {T<:AbstractFloat, A<:AbstractArray{T,2}}                   
        n1, n2 = size(x1)
        @assert ((n1,n2) == size(x2)) && ((n1,n2) == size(x3)) 
            if parameter_type == "stokes"
            I = x1            # Stokes parameter I (total light intensity)
            Q = x2            # Stokes parameter Q
            U = x3            # Stokes parameter U
            Iu = Array{T}(undef, n1, n2)    # intensity of unpolarized light
            Ip = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
            θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
            @inbounds for i2 in 1:n2
                @simd for i1 in 1:n1
                    Ip[i1,i2] = sqrt(Q[i1,i2]^2 + U[i1,i2]^2);
                    Iu[i1,i2] = I[i1, i2]-Ip[i1,i2];
                    θ[i1,i2] = atan(U[i1,i2], Q[i1,i2])/2;
                end
            end
        elseif parameter_type == "intensities"
            I = Array{T}(undef, n1, n2)    # Stokes parameter I (total light intensity)
            Q = Array{T}(undef, n1, n2)    # Stokes parameter Q
            U = Array{T}(undef, n1, n2)    # Stokes parameter U
            Iu = x1          # intensity of unpolarized light
            Ip = x2          # intensity of linearly polarized light
            θ = x3       # angle of linearly polarized light
             @inbounds for i2 in 1:n2
                @simd for i1 in 1:n1
                    I[i1,i2] = Iu[i1,i2] + Ip[i1,i2];
                    Q[i1,i2] = Ip[i1,i2] .*cos.(2*θ[i1,i2]);
                    U[i1,i2] = Ip[i1,i2] .*sin.(2*θ[i1,i2]);
                end
            end
        elseif parameter_type == "mixed" 
            I = Array{T}(undef, n1, n2)     # Stokes parameter I (total light intensity
            Q = x2            # Stokes parameter Q
            U = x1            # Stokes parameter U
            Iu = x3          # intensity of unpolarized light     
            Ip = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
            θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
            @inbounds for i2 in 1:n2
                @simd for i1 in 1:n1
                    Ip[i1,i2] = sqrt(Q[i1,i2]^2 + U[i1,i2]^2);
                    I[i1,i2] = Iu[i1,i2] + Ip[i1,i2];
                    θ[i1,i2] = atan(U[i1,i2], Q[i1,i2])/2;
                end
            end

        else
            error("unkown type, only known types : stokes, intensities and mixed")
        end
        PolarimetricMap(parameter_type,I,Q,U,Iu,Ip,θ)
    end

    function PolarimetricMap(parameter_type::AbstractString, 
                             x::Array{T,3}) where {T<:AbstractFloat}
        n1, n2, n3 = size(x)
        @assert n3 == 3 
        PolarimetricMap(parameter_type, 
                        view(x,:,:,1),
                        view(x,:,:,2),
                        view(x,:,:,3));
    end
    
#------------------------------------------------
# Base fonction redefinitions
    
    Base.size(A::PolarimetricMap) = size(A.I)
    Base.length(A::PolarimetricMap) =prod(size(A))*3
    Base.getindex(X::PolarimetricMap, i::CartesianIndex{2}) where {N} =
    PolarimetricPixel(X.parameter_type, X.I[i], X.Q[i], X.U[i], X.Iu[i], X.Ip[i], X.θ[i])
    Base.getindex(X::PolarimetricMap, i::Int) = 
    PolarimetricPixel(X.parameter_type, X.I[i], X.Q[i], X.U[i], X.Iu[i], X.Ip[i], X.θ[i])
    Base.getindex(X::PolarimetricMap, i::Int, j::Int) where {T<:Tuple} = 
    getindex(X, CartesianIndex(i,j))

    +(x::PolarimetricMap, y::PolarimetricMap) = PolarimetricMap(x.parameter_type,
                                                                x.I + y.I, 
                                                                x.Q + y.Q, 
                                                                x.U + y.U,
                                                                x.Iu + y.Iu,
                                                                x.Ip + y.Ip, 
                                                                x.θ + y.θ)
                                                                

     -(x::PolarimetricMap, y::PolarimetricMap) = PolarimetricMap(x.parameter_type,
                                                                x.I - y.I, 
                                                                x.Q - y.Q, 
                                                                x.U - y.U,
                                                                x.Iu - y.Iu,
                                                                x.Ip - y.Ip, 
                                                                x.θ - y.θ)
      
     vcopy(x::PolarimetricMap) = PolarimetricMap(x.parameter_type,
                                                 x.I, 
                                                 x.Q, 
                                                 x.U,
                                                 x.Iu,
                                                 x.Ip, 
                                                 x.θ)
     function vcreate(x::PolarimetricMap)
     
        @assert (x.parameter_type == "stokes") | 
                (x.parameter_type == "intensities") | 
                (x.parameter_type == "mixed")

         return PolarimetricMap(x.parameter_type,
                                Array{Float64,2}(undef,size(x)),
                                Array{Float64,2}(undef,size(x)),
                                Array{Float64,2}(undef,size(x)),
                                Array{Float64,2}(undef,size(x)),
                                Array{Float64,2}(undef,size(x)),
                                Array{Float64,2}(undef,size(x)))
     end     
                                        
     function +(x::PolarimetricMap, y::Array{T,3}) where {T<:AbstractFloat} 
        @assert size(y)[1:2] == size(x)       
        if x.parameter_type == "stokes"
           I=x.I + view(y,:,:,1);
           Q=x.Q + view(y,:,:,2);
           U=x.U + view(y,:,:,3);
           return PolarimetricMap("stokes", I, Q, U)
        elseif x.parameter_type == "intensities"
           Iu=x.Iu + view(y,:,:,1);
           Ip=x.Ip + view(y,:,:,2);
           θ=x.θ + view(y,:,:,3);
           return PolarimetricMap("intensities", Iu, Ip, θ)
        elseif x.parameter_type == "mixed"
           Iu=x.Iu + view(y,:,:,1);
           Q=x.Q + view(y,:,:,2);
           U=x.U + view(y,:,:,3);
           return PolarimetricMap("mixed", Iu, Q, U)
        else
            error("unknown parameter type")
        end
     end
     
     +(y::Array{T,3}, x::PolarimetricMap) where {T<:AbstractFloat} = x + y
     -(x::PolarimetricMap, y::Array{T,3}) where {T<:AbstractFloat} = x + (-y)
     
     
     function convert(::Type{Array{T,3}}, x::PolarimetricMap{T}) where {T <:AbstractFloat}
         if x.parameter_type == "stokes"
           return cat(x.I, x.Q, x.U, dims=3)
        elseif x.parameter_type == "intensities"
           return cat(x.Iu, x.Ip, x.θ, dims=3)
        elseif x.parameter_type == "mixed"
           return cat(x.Iu, x.Q, x.U, dims=3)
        else
            error("unknown parameter type")
        end
     end
    
#------------------------------------------------
# Writting function to save PolarimetricMap in fits file
"""
    write(X,'filename.fits') 
    
where X is a PolarimetricMap, write a fitsfile

"""

function write(X::PolarimetricMap, filename::AbstractString)
    data=cat(X.Iu', X.Ip', X.θ', X.I', X.Q', X.U',dims=3)
    header=FitsHeader()#;MAPORDER = (2, "Iu, Ip, theta, I, Q, U"))
    fitsdata=FitsImage(data)#, hdr=header)

    write!(filename, fitsdata)
end

"""
    read('parameter_type','filename.fits') -> PolarimetricMap
    
create an object of type PolarimetricMap from a fits file with:
    - Parameters I, Q, U (i.e. parameter_type = 'stokes')
    - Parameters Iu, Ip and θ (i.e. parameter_type = 'intensities')
    - Parameters Iu, Q, U (i.e. parameter_type = 'mixed')
   
"""


function read(parameter_type::AbstractString, filename::AbstractString)
    X=read(FitsArray, filename);
    return PolarimetricMap(parameter_type, 
                           view(X,:,:,4)', 
                           view(X,:,:,5)', 
                           view(X,:,:,6)', 
                           view(X,:,:,1)', 
                           view(X,:,:,2)', 
                           view(X,:,:,3)')
end

