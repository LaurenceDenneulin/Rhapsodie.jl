using Rhapsodie
using DelimitedFiles
using EasyFITS

#include("test_separable_reconstruction.jl")

tau=0.25;

max_iter = 100
α=10^-5
par=readdlm("data/Parameters.txt")
DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
# DerotAng = deg2rad.(readdlm("data_for_demo/pds70_angles.txt", Float64)[1:64])
max_angle = 64
DerotAng = [deg2rad(i) for i in range(1, max_angle, length=64)]
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
    push!(DerotAng, deg2rad(-par[7]));
end
		
psf_center=readdlm("data/PSF_centers_Airy.txt");

Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, DerotAng)


const I_HEADER_POS = 1
const I_STAR_HEADER_POS = 2
const I_DISK_HEADER_POS = 3
const Q_HEADER_POS = 4
const U_HEADER_POS = 5
const IU_HEADER_POS = 6
const IP_DISK_HEADER_POS = 9
const THETA_HEADER_POS = 10

true_file = readfits("test_results/TRUE.fits")
true_polar_map = PolarimetricMap("mixed",
                                view(true_file, :, :, I_STAR_HEADER_POS)' .+ view(true_file, :, :, I_DISK_HEADER_POS)',
                                view(true_file, :, :, Q_HEADER_POS)',
                                view(true_file, :, :, U_HEADER_POS)', 
                                view(true_file, :, :, IU_HEADER_POS)',
                                view(true_file, :, :, IP_DISK_HEADER_POS)',
                                view(true_file, :, :, THETA_HEADER_POS)')

Load_Data("test_results/DATA.fits", 
          "test_results/WEIGHT.fits")

PSF=readfits("data/PSF_parametered_Airy.fits");
const A=set_fft_op(PSF[1:end÷2,:]'[:,:],psf_center[1:2]);

X0 = PolarimetricMap("mixed", zeros(Rhapsodie.get_par().cols));
regularisation_parameters = 10 .^[0.5 , -1. , -1., -3.]; #(in log10)
regularisation_parameters[1] = 0

regularisation_parameter_list = [10^i for i in range(-5, 1, length=16)]

mse_list = Vector{Vector{Float64}}(undef, length(regularisation_parameter_list))
for k=1:length(regularisation_parameter_list)
    println("------Iteration: ", k, "------")
    regularisation_parameters[4] = regularisation_parameter_list[k]
    @time x = Rhapsodie.apply_rhapsodie(X0, A, Rhapsodie.dataset, regularisation_parameters,
                            maxeval=1000, maxiter=1000);

    crop!(x)
    mse_list[k] = Rhapsodie.MSE_object(x, true_polar_map)
    println("MSE: ", mse_list[k])
end
writedlm("test_results/rhapsodie_method_results/mse_list.txt", mse_list)