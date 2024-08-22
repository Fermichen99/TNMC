
mutable struct Spin
    Dx::Vector{Int}
    Dy::Vector{Int}
    # spinnew::Array{Int, 3}
    # spinold::Array{Int, 3}
    # lnprnew::Array{Float64,1}
    # lnprold::Array{Float64,1}
    bond::Array{Int, 3}

    function Spin(pa::Parament)
        Dx=[-1,0,1,0]
        Dy=[0,1,0,-1]
        # spinnew = zeros(Int, pa.L, pa.L, pa.Nreplic)
        # spinold = zeros(Int, pa.L, pa.L, pa.Nreplic)
        # lnprnew = zeros(Float64, pa.Nreplic)
        # lnprold = zeros(Float64, pa.Nreplic)
        bond = zeros(Int, pa.L, pa.L, 4)

        new(Dx, Dy, bond)
    end
end
export Spin

function InitBondFBC(pa::Parament,sp::Spin)
    for y in 1:pa.L, x in 1:pa.L
        for k in 1:div(round(Int, pa.nnb), 2)
            jx,jy=x+sp.Dx[k],y+sp.Dy[k]
            if 0<jx<=pa.L && 0<jy<=pa.L 
                # jy = mod1(y + sp.Dy[k], pa.L)
                # jx = mod1(x + sp.Dx[k], pa.L)
                if rand() < pa.P
                    bond = -1
                else 
                    bond = 1
                end
                sp.bond[y,x,k] = bond
                sp.bond[jy,jx,k+2] = bond
            end
        end
    end
end

function InitSpin(pa::Parament,sp::Spin)
    for re in 1:2, y in 1:pa.L, x in 1:pa.L
        sp.spin[re,y,x] = rand([1,-1])
    end
end

        
export InitBondFBC