include("SpinGlass.jl")
using .SpinGlass
using Random
using FileIO
using MPI

# Function to parse the input data
function parse_input(input_data)
    # Split the input data into lines
    lines = split(input_data, '\n')

    # Initialize variables
    judge = ""
    L = 0
    Dimension = 0
    T = 0.0
    P = 0.0
    chi = 0
    Ndisorder = 0
    Nsample = 0
    Ntherm = 0
    Nreplic = 0
    seed = 0

    # Process each line of the input data
    for line in lines
        # Skip empty lines
        if isempty(line)
            continue
        end

        # Split the line into key and value
        key, value = split(line)

        # Assign values to variables based on the key
        if key == "Type"
            judge = value
        elseif key == "P"
            P = parse(Float64, value)
        elseif key == "L"
            L = parse(Int, value)
        elseif key == "Dimension"
            Dimension = parse(Int, value)
        elseif key == "T"
            T = parse(Float64, value)
        elseif key == "chi"
            chi = parse(Int, value)
        elseif key == "Ndisorder"
            Ndisorder = parse(Int, value)
        elseif key == "Nsample"
            Nsample = parse(Int, value)
        elseif key == "Ntherm"
            Ntherm = parse(Int, value) 
        elseif key == "Nreplic"
            Nreplic = parse(Int, value) 
        elseif key == "seed"
            seed = parse(Int, value)
        end
    end

    # Return the parsed values as a tuple
    return [T, P], [L, Dimension, chi, Ndisorder, Nsample, Ntherm, Nreplic, seed]
end

MPI.Init()
comm = MPI.COMM_WORLD
my_rank = MPI.Comm_rank(comm)
num_procs = MPI.Comm_size(comm)
root = 0
T = Array{Float64}(undef, 2)
parsed_params = Array{Int64}(undef, 8)
if my_rank == root
    input_data = read(stdin, String)
    T, parsed_params = parse_input(input_data) 
end
T = MPI.Bcast!(T, root, comm)
parsed_params = MPI.Bcast!(parsed_params, root, comm)
parsed_params[end] += my_rank
pa=Parament("Tensor",T...,parsed_params...,num_procs)
sp=Spin(pa)

st=Statistics(pa)
da=DataFile(pa)
te=Tensor(pa,Array{Any,2}(undef,pa.L,pa.L))
Random.seed!(parsed_params[end]+my_rank)
SimulationTNMH(pa,sp,te,st,da,comm) 
#SimulationMetro(pa,sp,st,da,comm)
MPI.Finalize()

