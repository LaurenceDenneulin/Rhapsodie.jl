include("../src/Rhapsodie.jl")
include("test_data_simulation.jl")
#-----------------------------------------------
# Loading the dataset parameters
par=Rhapsodie.readdlm("test_results/Parameters.txt")

DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);		
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
DerotAng=Vector{Float64}();
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
    push!(DerotAng, deg2rad(-par[7]));
end
		
psf_center=Rhapsodie.readdlm("../data/PSF_centers_Airy.txt");

Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, DerotAng)


for tau in [0.01, 0.03, 0.07, 0.1, 0.15, 0.25, 0.5]
#----------------------------------------------
# Loading the pre-processed data (i.e. bad pixels have been interpolated	
    Rhapsodie.Load_Data("test_results/DATA_$tau-$DSIZE.fits", 
                         "test_results/WEIGHT_$tau-$DSIZE.fits")
	Sdim=length(Rhapsodie.dataset)
    
    DATA=zeros(Rhapsodie.get_par().cols[1], Rhapsodie.get_par().cols[2], Sdim,2);
    WEIGHT=zeros(Rhapsodie.get_par().cols[1], Rhapsodie.get_par().cols[2], Sdim,2);

    ker= Rhapsodie.CatmullRomSpline(Float64, Rhapsodie.Flat)       
    input_size=(Rhapsodie.get_par().rows[1], Rhapsodie.get_par().rows[2]÷2);
    output_size= Rhapsodie.get_par().cols[1:2];

# Pre processing
    for i=1:Sdim
        T1=Rhapsodie.TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, inv(Rhapsodie.Trans_Table[i][1]))
        T2=Rhapsodie.TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, inv(Rhapsodie.Trans_Table[i][2]))
    
        I1=T1*Rhapsodie.dataset[i].data[:,1:end÷2]
        I2=T2*Rhapsodie.dataset[i].data[:,end÷2+1:end]
        W1=T1*Rhapsodie.dataset[i].weights[:,1:end÷2]
        W2=T2*Rhapsodie.dataset[i].weights[:,end÷2+1:end]

        DATA[:,:,i,1]=I1;
        DATA[:,:,i,2]=I2;
        WEIGHT[:,:,i,1]=W1;
        WEIGHT[:,:,i,2]=W2;
    end

#-----------------------------------------------------
# Reconstruction
# Non-Linear Separable inverse method
    #MI=Rhapsodie.NonLinear_Method(DATA, WEIGHT);
    #Rhapsodie.write(MI, "test_results/Results_Separable_NonLinear_$tau-$DSIZE.fits")

# Linear Separable inverse method
    ML=Rhapsodie.Linear_Method(DATA, WEIGHT);
    Rhapsodie.write(ML, "test_results/Results_Separable_Linear_$tau-$DSIZE.fits")

# Double Ratio    
    DR=Rhapsodie.Double_Ratio(DATA);
    Rhapsodie.write(DR, "test_results/Results_Separable_DoubleRatio_$tau-$DSIZE.fits")
    
# Double Difference   
    DD=Rhapsodie.Double_Difference(DATA);
    Rhapsodie.write(DD, "test_results/Results_Separable_DoubleDifference_$tau-$DSIZE.fits")
		
    empty!(Rhapsodie.dataset);
end

