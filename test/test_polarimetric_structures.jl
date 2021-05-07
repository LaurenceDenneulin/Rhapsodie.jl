include("../src/RHAPSODIE.jl")
S=abs.(randn(5,5,3));
X=RHAPSODIE.PolarimetricMap("stokes", S);
x=RHAPSODIE.PolarimetricPixel("stokes", S[1,1,:]);
X[1,1]
X[1,1].I
RHAPSODIE.write(X,"test_results/rand_test.fits")
RHAPSODIE.read("mixed","test_results/rand_test.fits")
