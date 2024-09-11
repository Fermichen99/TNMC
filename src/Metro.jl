function Metropolis(pa::Parament, sp::Spin, spin::Array{Int, 3}, Beta::Float64)
    Acc = 0 
    for i in 1:2
        for y in 1:pa.L, x in 1:pa.L
            ns = spin[y, x, i]
            Nei = 0
            for k in 1:pa.nnb
                ky = (y + sp.Dy[k] + pa.L - 1) % pa.L + 1
                kx = (x + sp.Dx[k] + pa.L - 1) % pa.L + 1
                # ky = mod1(y + sp.Dy[k], pa.L)
                # kx = mod1(x + sp.Dx[k], pa.L)
                Nei += spin[ky, kx, i] * sp.bond[y, x, k]
            end
            DeltaE = 2 * ns * Nei * Beta 
            if DeltaE < 0 || rand() < exp(-DeltaE)
                spin[y, x, i] *= -1
                Acc += 1
            end
        end
    end
    Acc /= 2*pa.Vol
    return Acc
end

function Metropolis_tensor(pa::Parament, sp::Spin, spin::Array{Int, 3}, Beta::Float64)
    for i in 1:2
        for y in 1:pa.L, x in 1:pa.L
            ns = spin[y, x, i]
            Nei = 0
            for k in 1:pa.nnb
                ky = (y + sp.Dy[k] + pa.L - 1) % pa.L + 1
                kx = (x + sp.Dx[k] + pa.L - 1) % pa.L + 1
                # ky = mod1(y + sp.Dy[k], pa.L)
                # kx = mod1(x + sp.Dx[k], pa.L)
                Nei += spin[ky, kx, i] * sp.bond[y, x, k]
            end
            DeltaE = 2 * ns * Nei * Beta 
            if DeltaE < 0 || rand() < exp(-DeltaE)
                spin[y, x, i] *= -1
            end
        end
    end
    return spin
end


# function Metropolis(pa::Parament, sp::Spin, spin::Array{Int, 2}, Beta::Float64, lnpr::Float64)
#     for y in 1:pa.L, x in 1:pa.L
#         ns = spin[y, x]
#         Nei = 0
#         for k in 1:pa.nnb
#             ky = (y + sp.Dy[k] + pa.L - 1) % pa.L + 1
#             kx = (x + sp.Dx[k] + pa.L - 1) % pa.L + 1
#             Nei += spin[ky, kx] * sp.bond[y, x, k]
#         end
#         DeltaE = 2 * ns * Nei * Beta 
#         if DeltaE < 0 || rand() < exp(-DeltaE)
#             spin[y, x] *= -1
#         end
#         pr = 1.0/(1.0 + exp(-2 * Nei * Beta))
#         if spin[y,x] == 1
#             lnpr += log(pr)
#         else
#             lnpr += log(1-pr)
#         end
#     end
#     return spin, lnpr
# end

function GetEnergy(pa::Parament,sp::Spin,spin::Array{Int, 2})
    Energy = 0.0
    for y in 1:pa.L, x in 1:pa.L
        ns = spin[y, x]
        Nei = 0
        for k in 1:pa.nnb
            # jx,jy=x+sp.Dx[k],y+sp.Dy[k]
            # if 0<jx<=pa.L && 0<jy<=pa.L 
                jy = mod1(y + sp.Dy[k], pa.L)
                jx = mod1(x + sp.Dx[k], pa.L)
                Nei += spin[jy, jx] * sp.bond[y, x, k]
            # end
        end
        Energy += -ns * Nei/2
    end
    return Energy
end

