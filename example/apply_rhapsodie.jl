using Revise
using RhapsodieDirect
using Rhapsodie
using DelimitedFiles
using AstroFITS
using InterpolationKernels

T=Float32
ker=CatmullRomSpline(T, Flat)

mkpath("test_results")


data = mapslices(transpose,readfits(Array{T,3},"pds70_data/DATA_processed_coro.fits"),dims=(1,2))
weight = mapslices(transpose,readfits(Array{T,3},"pds70_data/WEIGHT_processed_coro.fits"),dims=(1,2))

par=readdlm("pds70_data/Parameters.txt",T)
dit=readdlm("pds70_data/Ditering.txt",T)
polar_val=readdlm("pds70_data/instruments_values_with_crosstalk.txt",T)


object_params=ObjectParameters{T,Int64}((150,150),(75.,75.))
object_type = "mixed"

data_params=DatasetParameters{T,Int64}((par[1],2*par[1]), par[2], par[3],par[2]÷(par[3]*4), (par[6],par[5]))

indices=get_indices_table(data_params)
polar_params=set_default_polarisation_coefficients(T,indices)

field_params=FieldTransformParameters[]
ndit = Int((data_params.frames_total/size(dit)[1]) )
for i=1:data_params.frames_total
    it = (i-1)÷ndit +  1
    push!(field_params, FieldTransformParameters(ker,
                                                deg2rad(-par[6 + i]-1.75),
                                                (0.,0.) .+ (dit[it,2], dit[it,1]),
                                                (-par[end] , -par[end-1]).+ (dit[it,2], dit[it,1]),
                                                Tuple(polar_val[i,1:3]),
                                                Tuple(polar_val[i,4:6])))
end

field_transforms=load_field_transforms(object_params,
                                       data_params,
                                       field_params)

	
psf_center=readdlm("data/PSF_centers_Airy.txt",T)
psf=readfits(Array{T,2},"data/PSF_parametered_Airy.fits")
blur=set_fft_operator(object_params,(psf[1:end÷2,:]'), psf_center[1:2]; pad_size=50)[1]


data_cube,weights_cube = pre_processing(data,weight,object_params, data_params, field_params)
Sl=Linear_Method(data_cube,weights_cube,field_params)
S0=PolarimetricMap(object_type, Sl.Iu, Sl.Q, Sl.U)
write(S0, "test_results/init.fits")

#S1=PolarimetricMap("intensities", blur'*(blur*Sl.Iu), blur'*(blur*Sl.Ip), Sl.θ )

#H = LinearDirectModel(size(S0), (data_params.size...,data_params.frames_total),S0.parameter_type,field_transforms,blur)



λ_Iu = T(4.)
for λ_QU in collect(T,2.:4.)
    for mu in collect(T,-3.:0)

        reg_par = 10 .^[λ_Iu, λ_QU, λ_Iu + 5, λ_QU + mu]

        H = LinearDirectModel(object_params.size, size(data),object_type,field_transforms,blur)
           
        D=Dataset(data,weight,H)

        S=apply_rhapsodie(S0, D, reg_par; mem=3, maxeval=1000, maxiter=1000)

        write(S, "test_results/$(object_type)_invertion_mu=$(reg_par[2])_rho=$(reg_par[4]).fits")

    end
end

