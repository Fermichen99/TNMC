module SpinGlass
using MPI
include("Para.jl")
include("Spin.jl")
include("Metro.jl")
include("wrapping.jl")
include("Obs/Struct.jl")
include("Stat.jl")
include("Obs/Pdf.jl")
include("Obs/C4q.jl")
include("Obs/Lnk.jl")
include("Tensor.jl")
include("Simu.jl")


export SimulationMetro, SimulationTNMH

end