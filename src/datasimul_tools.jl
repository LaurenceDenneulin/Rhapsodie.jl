#
# datasimul_tools.jl
#
# Provide tools to simulate synthetic parameters and dataset.
#
# ------------------------------------------------
#
# This file is part of Rhapsodie
#
#
# Copyright (c) 2017-2021 Laurence Denneulin


#------------------------------------------------

function data_generator(model::AbstractArray{T,N}, weights::AbstractArray{T,N};bad=zero(T), seed = 1992) where {T<:AbstractFloat,N}   
    seed === nothing ||  Random.seed!(seed);
    
    data = Array{T}(undef, size(model));
    @inbounds for i in eachindex(data, weights)
        w=weights[i] 
        (isfinite(w) && w >= 0 ) || error("invalid weights")
        if w >0            
            data[i] = model[i]  +randn()/sqrt(w)    
        elseif w ==0 
            data[i]=bad;
        end
    end
    return data
end

function generate_model(S::PolarimetricMap, A::Mapping)
    @assert size(S) == get_par().cols[1:2];
    
    M=Array{Float64,3}(undef, get_par().rows[1], 
                              get_par().rows[2], 
                              get_par().dataset_length)
	
	ker= LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)    
    for k=1:get_par().dataset_length
        output_size=(get_par().rows[1], Int64(get_par().rows[2]/2));
        input_size= get_par().cols[1:2];
    	T1=TwoDimensionalTransformInterpolator(output_size, 
    	                                       input_size, 
    	                                       ker, 
    	                                       ker, 
    	                                       Trans_Table[k][1])
    	T2=TwoDimensionalTransformInterpolator(output_size, 
    	                                       input_size, 
    	                                       ker, 
    	                                       ker, 
    	                                       Trans_Table[k][2])
	    
	    F=FieldTransformOperator(get_par().cols, 
	                                        get_par().rows, 
	                                        get_par().v[k][1],
	                                        get_par().v[k][2],
	                                        T1,
	                                        T2,
	                                        A)
	    
	    M[:,:,k].=F * cat(S.I, S.Q, S.U, dims=3);	    
	end
    return M
end
    
function data_simulator(Good_Pix, tau, A::Mapping; ro_noise=8.5)
    map_size=get_par().cols[1:2];
    
    S = generate_parameters(map_size, tau);

    M=generate_model(S,A);
    
    VAR=max.(M,zero(eltype(M))) .+ro_noise^2
	W=Good_Pix ./ VAR
	D=data_generator(M, W)
	
	check_MSE(M,D,W);
	
    CS=PolarimetricMap("intensities", A*S.Iu, A*S.Ip, A*S.θ)
	return D,W,S,CS
end

function data_simulator(Good_Pix, tau, A::Mapping, S::PolarimetricMap; ro_noise=8.5)
    @assert size(S) == get_par().cols[1:2];
   
    M=generate_model(S,A);
    
    VAR=max.(M,zero(eltype(M))) .+ro_noise^2
	W=Good_Pix ./ VAR
	D=data_generator(M, W)
	
	check_MSE(M,D,W);
    CS=PolarimetricMap("intensities", A*S.Iu, A*S.Ip, A*S.θ)
	return D,W,CS
end

function check_MSE(model, data, weights)
	MSE = vdot(data-model, weights.*(data-model)) ;
	N=count(weights .> 0);
	println("MSE=$MSE, N=$N, MSE/N=$(MSE/N)");
end

function generate_parameters(map_size, tau)
	Ip=zeros(map_size);
	Iu=zeros(map_size);
	θ=zeros(map_size);
	STAR1=zeros(map_size);
	STAR2=zeros(map_size);
	
	for i=1:map_size[1]
    	for j=1:map_size[2]
    		r1=(map_size[1]+1)/2-i;
    		r2=(map_size[2]+1)/2-j;
    		if (r1^2+r2^2<=20^2)
        		Iu[i,j]=1000;
        		Ip[i,j]=tau*Iu[i,j]/(1-tau);
    		end 
    		if ((r1^2+r2^2>=25^2)&&(r1^2+r2^2<=27^2))
        		Iu[i,j]=1000;
        		Ip[i,j]=tau*Iu[i,j]/(1-tau);
    		end
    		if ((r1^2+r2^2>=32^2)&&(r1^2+r2^2<=40^2))
        		Iu[i,j]=1000;
        		Ip[i,j]=tau*Iu[i,j]/(1-tau);
    		end
    		θ[i,map_size[2]+1-j]=atan(j-map_size[2]/2,i-map_size[1]/2);
			STAR1[i,j]=200*exp(-((i-map_size[1]/2)^2+(j-map_size[2]/2)^2)/(2*75^2))
			STAR2[i,j]=100000*exp(-((i-map_size[1]/2)^2+(j-map_size[2]/2)^2)/(2*7^2))
			if ((((map_size[1]+1)/2-i)^2+((map_size[2]+1)/2-j)^2)<=10^2)
        		STAR2[i,j]=800;
        		Iu[i,j]=0;
        		Ip[i,j]=0;	
    		end
			if ((((map_size[1]+1)/2-i)^2+((map_size[2]+1)/2-j)^2)<=70^2)
        		STAR1[i,j]=50;		
    		end
		end
	end    
	θ=θ.*(Ip.!=0);
	STAR=STAR1+STAR2
	STAR[round(Int64,10*map_size[1]/16)-3,round(Int64,10*map_size[2]/16)]=20000.0;
	STAR[round(Int64,10*map_size[1]/16),round(Int64,10*map_size[2]/16)-3]=100000.0;

    return PolarimetricMap("intensities", Iu+STAR, Ip, θ);
end
