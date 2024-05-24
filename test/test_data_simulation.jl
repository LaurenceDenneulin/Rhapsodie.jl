using Rhapsodie
using DelimitedFiles
using EasyFITS
using Dates

if prod(readdir() .!= "test_results")     
    mkdir("test_results")
end

DSIZE=128;
NTOT=64;
Nframe=2;
Nrot=1
Nangle=Int64.(NTOT/(Nrot*Nframe*4))
Center=[DSIZE/2;  DSIZE/2];
DerotAng=Vector{Float64}();
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],[10.7365 , -1.39344]));
    push!(DerotAng, deg2rad(-300.0));
end

psf_center=readdlm("data/PSF_centers_Airy.txt");
load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, DerotAng)

writedlm("test_results/Parameters.txt", [DSIZE; NTOT; Nframe; Nrot; Center; 300; [10.7365 , -1.39344]]);

psf = readfits("data/PSF_parametered_Airy.fits");
const A=set_fft_op((psf[1:end÷2,:]'), get_par().psf_center[1]);

BadPixMap = rand(0.0:1e-16:1.0,(DSIZE, 2*DSIZE)).< 0.9;

for tau in [0.03, 0.07, 0.1, 0.15, 0.25, 0.5]

    data, weight, S, S_convolved = data_simulator(BadPixMap, tau, A);
    writefits("test_results/DATA_$tau-$DSIZE.fits",
    ["DATE" => (now(), "date of creation")],
    mapslices(transpose,data,dims=[1,2]), overwrite=true)

    writefits("test_results/WEIGHT_$tau-$DSIZE.fits",
    ["DATE" => (now(), "date of creation")],
    mapslices(transpose,weight,dims=[1,2]), overwrite=true)

    write_polar_map(S_convolved, "test_results/TRUE_$tau-$DSIZE.fits", overwrite=true)
end