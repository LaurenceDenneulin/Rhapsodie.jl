using Rhapsodie
using DelimitedFiles
using EasyFITS

max_iter = 700
α=10^-5
par=readdlm("data_for_demo/Parameters.txt")
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
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)
mse_list = Vector{Vector{Float64}}()

for k in range(-3, 0, step=0.5)
    parameter_type = "stokes"
    base_path = joinpath("test_results", "contrast_10e$(k)")

    x_true = read_and_fill_polar_map(parameter_type, joinpath(pwd(), base_path, "TRUE.fits"))
    x_est = read_and_fill_polar_map(parameter_type, joinpath(base_path, "rhapsodie_method_results", "max_iter_700", "RHAPSODIE_non_linear_results_regul_param_0.21877616239495523.fits"))
    
    res = MSE_object(x_est, x_true)

    println("I: ", res[1], " | I_star: ", res[2], " | I_disk: ", res[3], " | Q: ", res[4], " | U: ", res[5], " | Iu: ", res[6], " | Iu_star: ", res[7], " | Iu_disk: ", res[8], " | Ip: ", res[9], " | Theta: ", res[10])
    push!(mse_list, res)
end

writedlm("test_results/mse_list.txt", mse_list)