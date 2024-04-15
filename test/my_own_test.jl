f = open("test/filename.txt", "r")  # Open file in read mode
content = read(f, String)      # Read content of the file into a string
close(f)
println(content)