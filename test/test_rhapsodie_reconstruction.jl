include("../src/RHAPSODIE.jl")
#include("test_data_simulation.jl")

tau=0.1;

par=RHAPSODIE.readdlm("test_results/Parameters.txt")
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
		
psf_center=RHAPSODIE.readdlm("../data/PSF_centers_Airy.txt");

RHAPSODIE.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, DerotAng)

RHAPSODIE.Load_Data("test_results/DATA_$tau-$DSIZE.fits", 
                    "test_results/WEIGHT_$tau-$DSIZE.fits")
                         
PSF=RHAPSODIE.read(RHAPSODIE.FitsArray, "../data/PSF_parametered_Airy.fits");
const A=RHAPSODIE.set_fft_op(PSF[1:end÷2,:]'[:,:],psf_center[1:2]);


X0 = RHAPSODIE.read("stokes", "test_results/Results_Separable_Linear_$tau-$DSIZE.fits");

regularisation_parameters = [0.5 , 0. , -1., -3.]; #(in log10)

@time x = RHAPSODIE.apply_rhapsodie(X0, A, RHAPSODIE.dataset, regularisation_parameters, 
                          maxeval=3000, maxiter=3000);

RHAPSODIE.write(x, "test_results/RHAPSODIE_nonlinearresults_$tau-$DSIZE.fits")  

