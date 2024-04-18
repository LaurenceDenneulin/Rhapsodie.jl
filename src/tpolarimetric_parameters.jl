using EasyFITS


#------------------------------------------------
# TRAX Struct definition
struct TPolarimetricPixel{T<: AbstractFloat} 
    parameter_type::AbstractString   #either "stokes", "intensities" or "mixed"
    I::T
    Q::T
    U::T
    Iu_star::T
    Iu_disk::T
    Ip_disk::T # We supposed that the star's light is unpolarized
    θ::T
end

#For map of pixels
struct TPolarimetricMap{T<: AbstractFloat} 
    parameter_type::AbstractString   #either "stokes", "intensities" or "mixed"
    I::Array{T,2};
    Q::Array{T,2};
    U::Array{T,2};
    Iu_star::Array{T,2};
    Iu_disk::Array{T,2};
    Ip_disk::Array{T,2}; # We supposed that the star's light is unpolarized
    θ::Array{T,2};
end

TPolarimetricMap(parameter_type::AbstractString,        
                    I::AbstractArray{T,2},
                    Q::AbstractArray{T,2},
                    U::AbstractArray{T,2},
                    Iu_star::AbstractArray{T,2},
                    Iu_disk::AbstractArray{T,2},
                    Ip_disk::AbstractArray{T,2},
                    θ::AbstractArray{T,2}) where {T<:AbstractFloat} = 
                        TPolarimetricMap(parameter_type,        
                        convert(Array{T},I),
                        convert(Array{T},Q),
                        convert(Array{T},U),
                        convert(Array{T},Iu_star),
                        convert(Array{T},Iu_disk),
                        convert(Array{T},Ip_disk),
                        convert(Array{T},θ))

#------------------------------------------------
# Constructors 
"""
    TPolarimetricPixel(parameter_type, x) -> TPolarimetricPixel
    TPolarimetricMap(parameter_type, x) -> TPolarimetricMap
    
create an object of type PolarimetricParameter from either:
    - Parameters I, Q, U (i.e. parameter_type = 'stokes')
    - Parameters Iu, Ip and θ (i.e. parameter_type = 'intensities')
    - Parameters Iu, Q, U (i.e. parameter_type = 'mixed')
    
Each parameter can be called from the structure. For exemple with a 
construction from Stokes parameters S=(I,Q,U):

    using Rhapsodie
    X = PolarimetricParameter(S, 'stokes');
    X.I #yields the Stokes parameter I
    X.Ip #yields the polarized intensity Ip
    X[1,1] #yields a PolarimetricPix at the CartesianIndex (1,1);
    X[1,1].I #yields the Stokes parameter I at the CartesianIndex (1,1); 

    TPolarimetricMap(parameter_type, n1, n2) -> TPolarimetricMap
    
yields an empty    
"""    
    function TPolarimetricPixel(parameter_type::AbstractString, 
                             x1::T, 
                             x2::T, 
                             x3::T) where {T<:AbstractFloat}
        if parameter_type == "stokes"
            I = x1            # Stokes parameter I (total light intensity)
            Q = x2            # Stokes parameter Q
            U = x3            # Stokes parameter U
            Ip_disk = sqrt(Q^2 + U^2)    # intensity of linearly polarized light
            Iu_star = I-Ip_disk    # intensity of unpolarized light
            Iu_disk = I-Ip_disk # TODO: Where do we get the Iu_disk and Iu_star var ?
            θ = atan(U, Q)/2 # angle of linearly polarized light
        elseif parameter_type == "intensities"
            Iu_star = x1          # TODO: Where do we get the Iu_disk and Iu_star var ?
            Iu_disk = x1          # intensity of unpolarized light
            Ip_disk = x2          # intensity of linearly polarized light
            θ = x3       # angle of linearly polarized light
            I = Iu_star + Iu_disk + Ip_disk   # Stokes parameter I (total light intensity)
            Q = Ip_disk *cos(2*θ);    # Stokes parameter Q
            U = Ip_disk *sin(2*θ);    # Stokes parameter U
        elseif parameter_type == "mixed" 
            Iu_star = x1          # intensity of unpolarized light   
            Iu_disk = x1          # TODO: Where do we get the Iu_disk and Iu_star var ?
            Q = x2            # Stokes parameter Q
            U = x3            # Stokes parameter U  
            Ip_disk = sqrt(Q^2 + U^2)    # intensity of linearly polarized light
            I = Iu_star + Iu_disk + Ip_disk     # Stokes parameter I (total light intensity)
            θ = atan(U, Q)/2 # angle of linearly polarized light
        else
            error("unkown type, known types are : 'stokes', 'intensities' and 'mixed'")
        end
        TPolarimetricPixel(parameter_type, I, Q, U, Iu_star, Iu_disk, Ip_disk, θ)
    end

    function TPolarimetricPixel(parameter_type::AbstractString, 
                             x::N) where {N<:AbstractArray{Float64,1}}
        @assert length(x) == 3
        TPolarimetricPixel(parameter_type, 
                        x[1],
                        x[2],
                        x[3]);
    end 
    # Question: Why having not having only one func ? No need for the above constructor ? (difference is array in parameter)
    
    function TPolarimetricMap(x::Array{TPolarimetricPixel,2}) 
        T = Float64;
        n1, n2 = size(x)
        par_type=x[1].parameter_type;          
        I = Array{T}(undef, n1, n2)    # Stokes parameter I (total light intensity)
        Q = Array{T}(undef, n1, n2)    # Stokes parameter Q
        U = Array{T}(undef, n1, n2)    # Stokes parameter U
        Iu_star = Array{T}(undef, n1, n2)    # TODO
        Iu_disk = Array{T}(undef, n1, n2)    # intensity of unpolarized light
        Ip_disk = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
        θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1   
                pix= x[i1,i2];
                if par_type == pix.parameter_type
                   I[i1,i2] = pix.I;
                   Q[i1,i2] = pix.Q;
                   U[i1,i2] = pix.U;
                   Iu_star[i1,i2] = pix.Iu_star;
                   Iu_disk[i1,i2] = pix.Iu_disk;
                   Ip_disk[i1,i2] = pix.Ip_disk;
                   θ[i1,i2] = pix.θ;
                else
                   error("polarimetric types must match")
                end
            end 
        end
        TPolarimetricMap(par_type, I, Q, U, Iu_star, Iu_disk, Ip_disk, θ)
    end
 
    function TPolarimetricMap(parameter_type::AbstractString, 
                             x1::A, 
                             x2::A, 
                             x3::A) where {T<:AbstractFloat, A<:AbstractArray{T,2}}                   
        n1, n2 = size(x1)
        @assert ((n1,n2) == size(x2)) && ((n1,n2) == size(x3)) 
            if parameter_type == "stokes"
            I = x1            # Stokes parameter I (total light intensity)
            Q = x2            # Stokes parameter Q
            U = x3            # Stokes parameter U
            Iu_star = Array{T}(undef, n1, n2)    # intensity of unpolarized light
            Iu_disk = Array{T}(undef, n1, n2)    # TODO 
            Ip_disk = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
            θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
            @inbounds for i2 in 1:n2
                @simd for i1 in 1:n1
                    Ip_disk[i1,i2] = sqrt(Q[i1,i2]^2 + U[i1,i2]^2);
                    Iu_disk[i1,i2] = I[i1, i2]-Ip_disk[i1,i2];
                    Iu_star[i1,i2] = I[i1, i2]-Ip_disk[i1,i2]; # TODO
                    θ[i1,i2] = atan(U[i1,i2], Q[i1,i2])/2;
                end
            end
        elseif parameter_type == "intensities"
            I = Array{T}(undef, n1, n2)    # Stokes parameter I (total light intensity)
            Q = Array{T}(undef, n1, n2)    # Stokes parameter Q
            U = Array{T}(undef, n1, n2)    # Stokes parameter U
            Iu_star = x1          # intensity of unpolarized light
            Iu_disk = x1          # TODO
            Ip_disk = x2          # intensity of linearly polarized light
            θ = x3       # angle of linearly polarized light
             @inbounds for i2 in 1:n2
                # STOPPED HERE
                @simd for i1 in 1:n1
                    I[i1,i2] = Iu_star[i1,i2] + Iu_disk[i1,i2] + Ip_disk[i1,i2];
                    Q[i1,i2] = Ip[i1,i2] .*cos.(2*θ[i1,i2]);
                    U[i1,i2] = Ip[i1,i2] .*sin.(2*θ[i1,i2]);
                end
            end
        elseif parameter_type == "mixed" 
            I = Array{T}(undef, n1, n2)     # Stokes parameter I (total light intensity
            Q = x2            # Stokes parameter Q
            U = x3            # Stokes parameter U
            Iu = x1          # intensity of unpolarized light     
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
        TPolarimetricMap(parameter_type,I,Q,U,Iu,Ip,θ)
    end

    function TPolarimetricMap(parameter_type::AbstractString, 
                             x::Array{T,3}) where {T<:AbstractFloat}
        n1, n2, n3 = size(x)
        @assert n3 == 3 
        TPolarimetricMap(parameter_type, 
                        view(x,:,:,1),
                        view(x,:,:,2),
                        view(x,:,:,3));
    end
    
    function TPolarimetricMap(parameter_type::AbstractString, n1::Int, n2::Int)
        return TPolarimetricMap(parameter_type,
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2),
                               Array{Float64,2}(undef, n1, n2))
    end
    
#------------------------------------------------
# Base fonction redefinitions
    
    Base.size(A::TPolarimetricMap) = size(A.I)
    Base.length(A::TPolarimetricMap) =prod(size(A))*3
    Base.length(A::TPolarimetricPixel) =3
    Base.getindex(X::TPolarimetricMap, i::CartesianIndex{2}) where {N} =
    TPolarimetricPixel(X.parameter_type, X.I[i], X.Q[i], X.U[i], X.Iu[i], X.Ip[i], X.θ[i])
    Base.getindex(X::TPolarimetricMap, i::Int) = 
    TPolarimetricPixel(X.parameter_type, X.I[i], X.Q[i], X.U[i], X.Iu[i], X.Ip[i], X.θ[i])
    Base.getindex(X::TPolarimetricMap, i::Int, j::Int) where {T<:Tuple} = 
    getindex(X, CartesianIndex(i,j))

    function Base.setindex!(X::TPolarimetricMap{Float64}, x::TPolarimetricPixel{Float64}, i::Int64, j::Int64)
        X.I[i,j]=x.I;
        X.Q[i,j]=x.Q;
        X.U[i,j]=x.U;
        X.Iu[i,j]=x.Iu;
        X.Ip[i,j]=x.Ip;
        X.θ[i,j]=x.θ;
    end
        

    +(x::TPolarimetricMap, y::TPolarimetricMap) = TPolarimetricMap(x.parameter_type,
                                                                x.I + y.I, 
                                                                x.Q + y.Q, 
                                                                x.U + y.U,
                                                                x.Iu + y.Iu,
                                                                x.Ip + y.Ip, 
                                                                x.θ + y.θ)
                                                                

     -(x::TPolarimetricMap, y::TPolarimetricMap) = TPolarimetricMap(x.parameter_type,
                                                                x.I - y.I, 
                                                                x.Q - y.Q, 
                                                                x.U - y.U,
                                                                x.Iu - y.Iu,
                                                                x.Ip - y.Ip, 
                                                                x.θ - y.θ)
      
     vcopy(x::TPolarimetricMap) = TPolarimetricMap(x.parameter_type,
                                                 x.I, 
                                                 x.Q, 
                                                 x.U,
                                                 x.Iu,
                                                 x.Ip, 
                                                 x.θ)
     function vcreate(x::TPolarimetricMap)
     
        @assert (x.parameter_type == "stokes") | 
                (x.parameter_type == "intensities") | 
                (x.parameter_type == "mixed")
         n1,n2=size(x);
         return TPolarimetricMap(x.parameter_type, n1, n2)
     end     
                                        
     function +(x::TPolarimetricMap, y::Array{T,3}) where {T<:AbstractFloat} 
        @assert size(y)[1:2] == size(x)       
        if x.parameter_type == "stokes"
           I=x.I + view(y,:,:,1);
           Q=x.Q + view(y,:,:,2);
           U=x.U + view(y,:,:,3);
           return TPolarimetricMap("stokes", I, Q, U)
        elseif x.parameter_type == "intensities"
           Iu=x.Iu + view(y,:,:,1);
           Ip=x.Ip + view(y,:,:,2);
           θ=x.θ + view(y,:,:,3);
           return TPolarimetricMap("intensities", Iu, Ip, θ)
        elseif x.parameter_type == "mixed"
           Iu=x.Iu + view(y,:,:,1);
           Q=x.Q + view(y,:,:,2);
           U=x.U + view(y,:,:,3);
           return TPolarimetricMap("mixed", Iu, Q, U)
        else
            error("unknown parameter type")
        end
     end
     
     +(y::Array{T,3}, x::TPolarimetricMap) where {T<:AbstractFloat} = x + y
     -(x::TPolarimetricMap, y::Array{T,3}) where {T<:AbstractFloat} = x + (-y)
     
     
     function convert(::Type{Array{T,3}}, x::TPolarimetricMap{T}) where {T <:AbstractFloat}
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
    
const IU_HEADER_POS = 1
const IP_HEADER_POS = 2
const THETA_HEADER_POS = 3
const I_HEADER_POS = 4
const Q_HEADER_POS = 5
const U_HEADER_POS = 6
#------------------------------------------------
# Writting function to save TPolarimetricMap in fits file
"""
    write(X,'filename.fits') 
    
where X is a TPolarimetricMap, write a fitsfile

"""
function write_polar_map(X::TPolarimetricMap, filename::AbstractString; overwrite::Bool = false)
    data = cat(X.Iu', X.Ip', X.θ', X.I', X.Q', X.U',dims=3)
    writefits(filename,
    ["D" => ("Ok", "")],
    data, overwrite=overwrite)
end

"""
    read('parameter_type','filename.fits') -> TPolarimetricMap
    
create an object of type TPolarimetricMap from a fits file with:
    - Parameters I, Q, U (i.e. parameter_type = 'stokes')
    - Parameters Iu, Ip and θ (i.e. parameter_type = 'intensities')
    - Parameters Iu, Q, U (i.e. parameter_type = 'mixed')
   
"""

function read_and_fill_polar_map(parameter_type::AbstractString, filename::AbstractString)
    X=readfits(filename);
    return TPolarimetricMap(parameter_type,
                           view(X,:,:,I_HEADER_POS)',
                           view(X,:,:,Q_HEADER_POS)',
                           view(X,:,:,U_HEADER_POS)',
                           view(X,:,:,IU_HEADER_POS)',
                           view(X,:,:,IP_HEADER_POS)',
                           view(X,:,:,THETA_HEADER_POS)')
end
