
struct Parament
    judge::String
    L::Int
    Dimension::Int
    chi::Int
    Nblck::Int
    Ndisorder::Int
    Nsample::Int
    Ntherm::Int
    Nreplic::Int
    Nmpi::Int
    seed::Int
    # judge::Any
    P::Float64
    T::Float64
    SWP::Float64
    Beta::Float64
    dn :: Float64

    Vol::Int
    nnb::Int
    Barray::Array{Float64, 2}
    BAntiarray::Array{Float64, 2}
    field::Array{Float64, 1}
    Antifield::Array{Float64, 1}
    I1::Array{Float64, 1}
    I2::Array{Float64, 2}
    I3::Array{Float64, 3}
    I4::Array{Float64, 4}
    BondTensor::Array{Any, 1}
    
    # Constructor for Parament
    function Parament(judge, T, P, L, Dimension, chi, Ndisorder, Nsample, Ntherm, Nreplic, seed, Nmpi)
        Vol = L^Dimension
        SWP = 1 - exp(-2/T)
        # Ntherm = div(round(Int, Nsample), 5)
        Nblck = Ndisorder
        nnb = Dimension * 2
        Beta = 1.0 / T
        
        dn = Dimension*(L-1)/L
        Barray = exp.(Beta * [1.0 -1.0; -1.0 1.0])
        BAntiarray = exp.(-Beta * [1.0 -1.0; -1.0 1.0])
        field = vec(exp.(Beta * [1.0 -1.0]))
        Antifield = vec(exp.(-Beta * [1.0 -1.0]))
        I1 = vec([1.0 1.0])
        I2 = zeros(2, 2)
        I3 = zeros(2, 2, 2)
        I4 = zeros(2, 2, 2, 2)
        for i in 1:2
            I2[i ,i] = 1.0
            I3[i, i, i] = 1.0
            I4[i, i, i, i] = 1.0
        end
        BondTensor = [reshape([1], 1, 1), deepcopy(Barray), deepcopy(BAntiarray)]
        new(judge, L, Dimension, chi, Nblck,Ndisorder, Nsample, Ntherm, Nreplic, Nmpi, seed, P, T, SWP, Beta, dn, Vol, nnb, 
            Barray, BAntiarray, field, Antifield, I1, I2, I3, I4, BondTensor)
    end
end
export Parament