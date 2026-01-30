#
# Rhapsodie.jl
#
# Package for the Reconstruction of High-contrAst Polarized
# SOurces and Deconvolution for cIrcumstellar Environments (Rhapsodie)
#
#----------------------------------------------------------
#
# 
# Copiyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

module Rhapsodie

    export
        apply_rhapsodie,
        apply_edge_preserving_smoothing!,
        Double_Difference,
        Double_Ratio,
        Linear_Method,
        mse_intensities,
        NonLinear_Method

    using OptimPackNextGen
    using AstroFITS
    using DelimitedFiles
    using RhapsodieDirect
    using Statistics
    using StaticArrays

    include("separable_methods.jl") #FIXME : PolarimetricPixels does not exist anymore
    include("rhapsodie_methods.jl")
    include("sure_tools.jl")
end

