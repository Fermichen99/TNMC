function generate_sequence(n::Int)
    sequence = Int[]
    push!(sequence, 1)
    current_value = 2
    while current_value <= n
        push!(sequence, current_value)
        current_value *= 2
    end
    return sequence
end

struct C4q
    r::Vector{Int}; q::Vector{Float16}
    lenr::Int;  lenq::Int
    Quan::Array{Float64, 3}
    Ave::Array{Float64, 2}
    Dev::Array{Float64, 2}

    function C4q(pa::Parament)
        r = generate_sequence(pa.L÷2)
        q = range(0, stop=1, length=21)
        lenr = length(r);  lenq = length(q)
        Quan = zeros(Float64, lenr, lenq,  pa.Ndisorder)
        Ave = zeros(Float64, lenr, lenq)
        Dev = zeros(Float64, lenr, lenq)
        
        new(r, q, lenr, lenq, Quan, Ave, Dev)
    end
end
export C4q

struct Lnk
    q::Vector{Float16}
    lenq::Int
    Quan::Array{Float64, 2}
    Ave::Array{Float64, 1}
    Dev::Array{Float64, 1}

    function Lnk(pa::Parament)
        q = range(0, stop=1, length=21)
        lenq = length(q)
        Quan = zeros(Float64, lenq,  pa.Ndisorder)
        Ave = zeros(Float64, lenq)
        Dev = zeros(Float64, lenq)
        
        new(q, lenq, Quan, Ave, Dev)
    end
end
export Lnk

struct Pdf
    Npdf::Int
    PdfCoe::Float64
    Quan::Array{Float64, 2}
    Ave::Vector{Float64}
    Dev::Vector{Float64}
    Qva::Vector{Float64}

    function Pdf(pa::Parament)
        Npdf = 41
        PdfCoe = sqrt(pa.Vol/ (2 *pi))
        Quan = zeros(Float64, Npdf, pa.Ndisorder)
        Ave = zeros(Float64, Npdf)
        Dev = zeros(Float64, Npdf)
        Qva = range(-1, stop=1, length=Npdf)
        
        new(Npdf, PdfCoe, Quan, Ave, Dev, Qva)
    end
end
export Pdf


