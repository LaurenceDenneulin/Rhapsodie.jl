using Rhapsodie
using DelimitedFiles
using EasyFITS
using Dates

if prod(readdir() .!= "test_results")     
    mkdir("test_results")
end

DSIZE=256;
NTOT=64;
Nframe=2;
Nrot=NTOT
Nangle=NTOT ÷ (Nframe*4)
Center=[DSIZE/2, DSIZE/2];
DerotAng = deg2rad.(readdlm("data_for_demo/pds70_angles.txt", Float64)[1:64])
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],[10.7365 , -1.39344]));
end

psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");
Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)

writedlm("data_for_demo/Parameters.txt", [DSIZE; NTOT; Nframe; Nrot; Center; 300; [10.7365 , -1.39344]]);

psf = readfits("data_for_demo/PSF_parametered_Airy.fits");
const A=set_fft_op((psf[1:end÷2,:]'), get_par().psf_center[1]);


ddit_fits = readfits("data_for_demo/ddit_simulated_data.fits");
I_tot = view(ddit_fits, :, :, 1)
max_I_tot = maximum(I_tot)
I_tot *= 100 / max_I_tot

Ip_disk = ddit_fits[:,:,2]
Ip_disk *= 100 / max_I_tot

scattering = ddit_fits[:,:,3]
star_map = Rhapsodie.generate_parameters(size(Ip_disk), 0)
Iu_star = star_map.I_star
Iu_disk = I_tot - Ip_disk


X0 = Rhapsodie.TPolarimetricMap("intensities", Iu_star, Iu_disk, Ip_disk, scattering)

GoodPixMap = rand(0.0:1e-16:1.0,(DSIZE, 2*DSIZE)).< 0.9;

data, weight, S, S_convolved = Rhapsodie.ddit_data_simulator(GoodPixMap, A, X0, ro_noise=8.5);
writefits("test_results/DATA.fits",
["DATE" => (now(), "date of creation")],
mapslices(transpose,data,dims=[1,2]), overwrite=true)

writefits("test_results/WEIGHT.fits",
["DATE" => (now(), "date of creation")],
mapslices(transpose,weight,dims=[1,2]), overwrite=true)

Rhapsodie.write_polar_map(S_convolved, "test_results/TRUE_convolved.fits", overwrite=true)
Rhapsodie.write_polar_map(S, "test_results/TRUE.fits", overwrite=true)