using Random, Graphs

# nnb = 4; nnbh = nnb÷2
# Dx = [0, -1, 1, 0]
# Dy = [-1, 0, 0, 1]

nnb = 6; nnbh = nnb÷2
Dx = [-1, 0, -1, 1, 0, 1]
Dy = [-1, -1, 0, 0, 1, 1]
# Initialize directions based on lattice connectivity
# function initialize_directions(nnb::Int)
#     if nnb == 4
#         Dx = [0, -1, 1, 0]
#         Dy = [-1, 0, 0, 1]
#     elseif nnb == 6
#         Dx = [-1, 0, -1, 1, 0, 1]
#         Dy = [-1, -1, 0, 0, 1, 1]
#     end
#     Dx, Dy
# end

function make_bonds(sp::Spin, spin_lattice::Array{Int,2}, p::Float64, M::Float64)
    N = size(spin_lattice, 1)
    bondmaj = SimpleGraph(N^2)
    bondmin = SimpleGraph(N^2)
    spin = M > 0 ? 1 : -1

    for y in 1:N
        for x in 1:N
            for b in nnbh+1:nnb
                ry = y + Dy[b]
                jy = ry > N ? 1 : (ry < 1 ? N : ry)
                rx = x + Dx[b]
                jx = rx > N ? 1 : (rx < 1 ? N : rx)
                if spin_lattice[jy, jx] == spin_lattice[y, x] && rand() < p
                    if spin_lattice[y, x] == spin
                        add_edge!(bondmaj, (y - 1) * N + x, (jy - 1) * N + jx)
                    else
                        add_edge!(bondmin, (y - 1) * N + x, (jy - 1) * N + jx)
                    end
                end
            end
        end
    end
    return bondmaj, bondmin
end 


function check_wrapping(N::Int, component, bonds)
    ID = zeros(Int, N^2); visited = zeros(Int, N^2)
    for grid in component
        ID[grid] = 1
    end
    id = component[1]; queue = []; push!(queue, id)

    rl0 = 0; rx0 = []; ry0 = []; newlp = 0
    Tx  = zeros(Int, N^2); Ty = zeros(Int, N^2)
    while !isempty(queue)
        current = popfirst!(queue)
        visited[current] = 1
        y = (current-1)÷N + 1;  x = (current-1)%N+1
        for b in 1:nnb
            ry = y + Dy[b];  jy = ry > N ? 1 : (ry < 1 ? N : ry)
            rx = x + Dx[b];  jx = rx > N ? 1 : (rx < 1 ? N : rx)
            index = (jy-1)*N + jx
            cx = Tx[current]+Dx[b]; cy = Ty[current]+Dy[b]
            if ID[index] == 0 continue end
            if !has_edge(bonds, current, index) continue end
            if visited[index] == 0 
                push!(queue, index); visited[index] = 1
                Tx[index] = cx;  Ty[index] = cy
            elseif cx != Tx[index] || cy != Ty[index]
                jkx = abs(cx-Tx[index])÷N; jky = abs(cy-Ty[index])÷N
                if rl0 == 0
                    rl0 = 1; push!(rx0,jkx); push!(ry0,jky)
                else
                    newlp = 1
                    for i in 1:rl0
                        if rx0[i] == jkx && ry0[i] == jky
                            newlp = 0; break
                        end
                    end
                    if newlp == 1
                        rl0 += 1; push!(rx0,jkx); push!(ry0,jky)
                    end
                end
            end
        end
    end
    if rl0 == 0
        return 0
    elseif rl0 == 1
        return 1
    elseif rl0 > 1
        return 2
    end
end


# 检测Percolation函数
function Wrapping(bonds::SimpleGraph, N::Int)
    labels = connected_components(bonds)
    component = last(sort(labels, by=length))
    R = check_wrapping(N, component, bonds) 

    R0 = (R==0) ? 1 : 0
    R1 = (R==1) ? 1 : 0
    R2 = (R==2) ? 1 : 0

    return R0, R1, R2
end

# 主函数
function Wrapping_R(pa::Parament,sp::Spin,spin_lattice::Array{Int,2},P::Float64, M::Float64)
    bondmaj,bondmin = make_bonds(sp,spin_lattice, P, M)
    R0maj,R1maj,R2maj = Wrapping(bondmaj, pa.L)
    R0min,R1min,R2min = Wrapping(bondmin, pa.L)
    R2tot = (R2maj + R2min > 0) ? 1 : 0
    R0tot = (R0maj + R0min == 2) ? 1 : 0
    R1tot = (R2tot==0 && R0tot==0) ? 1 : 0
    return R0tot,R1tot,R2tot,R0maj,R1maj,R2maj,R0min,R1min,R2min
end

pl = [i for i in 0.1:0.1:1]
function GetR(pa::Parament,sp::Spin,spin_lattice::Array{Int,2},M::Float64)
    #pl = [0.5]
    RP = []
    for p in pl
        push!(RP,Wrapping_R(pa,sp,spin_lattice,p,M))
    end
    return RP
end








