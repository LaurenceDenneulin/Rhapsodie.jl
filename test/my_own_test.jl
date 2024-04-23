using Rhapsodie

model = Array{Float64, 2}(undef, 2, 3)
weights = ones(Float64, 2, 3)

data = data_generator(model, weights)

println(data)


atan(cos(6),sin(6))/2