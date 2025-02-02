function read_input_file(filename)
    open(filename, "r") do file
        lines = readlines(file)
        
        # Parse number of vertices
        n = parse(Int, split(lines[1], "=")[2])
        
        # Parse travel time matrix t
        t = zeros(n, n)
        t_lines = split(strip(split(lines[2], "=")[2]), ";")
        for i in 1:n
            t[i, :] = parse.(Int, split(strip(t_lines[i])))
        end
        
        # Parse nominal travel times th
        th = parse.(Int, split(strip(split(lines[3], "=")[2]), ","))
        
        # Parse threshold T
        T = parse(Int, split(lines[4], "=")[2])
        
        # Parse demands d
        d = parse.(Int, split(strip(split(lines[5], "=")[2]), ","))
        
        # Parse capacity C
        C = parse(Int, split(lines[6], "=")[2])
        
        # Define sets
        V = 1:n
        A = [(i, j) for i in V, j in V if i != j]  # Set of arcs
        
        return V, A, th, t, d, C, T
    end
end
