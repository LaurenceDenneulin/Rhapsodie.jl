using Rhapsodie

x = PolarimetricPixel("intensities", 100., 1., pi/8);
n1,n2 = (500, 499);
X = PolarimetricMap("intensities", n1,n2)

S=abs.(randn(5,5,3));
X=Rhapsodie.PolarimetricMap("stokes", S);
x=Rhapsodie.PolarimetricPixel("stokes", S[1,1,:]);
X[1,1]
X[1,1].I
Rhapsodie.write(X,"test_results/rand_test.fits")
Rhapsodie.read("mixed","test_results/rand_test.fits")
