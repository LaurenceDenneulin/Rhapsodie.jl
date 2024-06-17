using Rhapsodie
using DelimitedFiles
using EasyFITS

#include("test_separable_reconstruction.jl")

par=readdlm("data_for_demo/Parameters.txt")
DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
DerotAng = deg2rad.(readdlm("data_for_demo/pds70_angles.txt", Float64)[1:64])
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)

Rhapsodie.load_data("test_results/DATA.fits", 
          "test_results/WEIGHT.fits")

PSF=readfits("data_for_demo/PSF_parametered_Airy.fits");
const A=set_fft_op(PSF[1:end÷2,:]'[:,:],psf_center[1:2]);

X0 = TPolarimetricMap("mixed", zeros(Rhapsodie.get_par().cols));
regularisation_parameters = 10 .^[0.5 , -1. , 3, -3.]; #(in log10) star, disk


@time x = apply_rhapsodie(X0, A, Rhapsodie.dataset, regularisation_parameters,
                          maxeval=1000, maxiter=100);

crop!(x)
write_polar_map(x, "test_results/RHAPSODIE_non_linear_results.fits", overwrite=true)