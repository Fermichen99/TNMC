using LinearAlgebra
using BenchmarkTools
using OMEinsum:@ein_str
import Base.GC: gc
export InitTensor, TNMH
BLAS.set_num_threads(1) 


struct MpsState
    Tensor::Array{Any}  # Allow for different array types
    lnZ::Array{Float64, 1}
end

mutable struct Tensor
    len::Int
    StoreTensor::Array{Any,2}
    HorizontTensorsDw::MpsState
    VerticalTensorsLe::MpsState
    VerticalTensorsRi::MpsState
    
    function Tensor(pa::Parament, StoreTensor::Array{Any,2}) 
        len = size(StoreTensor)[1]
        HorizontTensorsDw = MpsState(Array{Any,2}(undef, pa.L, pa.L), zeros(Float64, pa.L))
        VerticalTensorsLe = MpsState(Array{Any,1}(undef, pa.L), zeros(Float64, pa.L))
        VerticalTensorsRi = MpsState(Array{Any,1}(undef, pa.L), zeros(Float64, pa.L))
        new(len, StoreTensor, HorizontTensorsDw, VerticalTensorsLe, VerticalTensorsRi)
    end
end

export Tensor

function safe_svd(A::Matrix{Float64}; tol::Float64=1e-12)
    try
        F = svd(A)
        return F
    catch e
        println("SVD failed with error: $e. Checking matrix properties.")
        if any(isnan, A) || any(isinf, A)
            println("Matrix contains NaN or Inf. Replacing with small value.")
            A = replace(A, x -> isnan(x) || isinf(x) ? eps() : x)
        else
            # 检查矩阵的条件数是否过大
            cond_number = cond(A)
            if cond_number > 1e12
                println("Matrix is near singular. Adding small value to diagonal.")
            end
             # 检查矩阵的秩是否过低
            rank_A = rank(A)
            println("Matrix rank: $rank_A / $(min(size(A)...))")
            if rank_A < min(size(A)...)
                println("Matrix is low-rank. Adding small value to diagonal.")
            end
            # 计算矩阵A的最小维度
            min_dim = min(size(A)...)
            for i in 1:min_dim
                A[i, i] += tol
            end
        end
        F = svd(A)
        return F  
    end
end

function GetAntiPara(sp::Spin,nnb::Int, y::Int, x::Int, ju::Int)
    if sp.bond[y, x, nnb] == 0 
        if  nnb == 4
            return 4
        else
            return 1
        end
    else
        if ju == 1
            if sp.bond[y, x, nnb] == 1
                return 2
            else
                return 3
            end 
        elseif ju == 0
            return 4
        end
    end
end

const STORE_CACHE = Dict{Tuple{Vararg{Int}},Any}()

function OutTensorGrid(pa::Parament, sp::Spin, te::Tensor)
    bond_indices = [(1,0), (3,1), (2,1), (4,0)]
    bond_configs = [Tuple([GetAntiPara(sp, dir, y, x, param) for (dir, param) in bond_indices])
                    for y in 1:pa.L, x in 1:pa.L]
    code = ein"ijkl,ai,bj,ck,dl->abcd"
    @inbounds for idx in CartesianIndices(te.StoreTensor)
        current_config = bond_configs[idx]
        current_bonds = [pa.BondTensor[i] for i in current_config]
        cached_value = get!(STORE_CACHE, current_config) do
            code(pa.I4, current_bonds...)
        end
        te.StoreTensor[idx] = cached_value
    end

end           

function reshape_tensor(pa::Parament,Ftensor)
    for x in 1:pa.L
        dims = size(Ftensor[pa.L,x])
        Ftensor[pa.L,x] = reshape(Ftensor[pa.L,x], (dims[1:2]...,dims[4]...))
    end
end

function InitTensor(pa::Parament, sp::Spin, te::Tensor)
    OutTensorGrid(pa, sp, te)
    reshape_tensor(pa,te.StoreTensor)
    StorelnZ(pa, te)
end


function StorelnZ(pa::Parament,te::Tensor)
    lnZ=0
    te.HorizontTensorsDw.lnZ[pa.L]=0 
    te.HorizontTensorsDw.Tensor[pa.L,:] .= te.StoreTensor[pa.L,:]
    for y in pa.L:-1:2
        res,te.HorizontTensorsDw.Tensor[y-1,:] = Compress(Eat(pa,te.HorizontTensorsDw.Tensor[y,:] ,te.StoreTensor[y-1,:]),pa.chi)
        lnZ += res
        te.HorizontTensorsDw.lnZ[y-1] = lnZ
    end
end

function Eat(pa::Parament, mps::Array{T,S}, mpo::Array{T,O}) where {T,S,O}
    result = similar(mps)  # 预分配同类型内存
    for i in 1:pa.L
        tmp = ein"ikj,abjc->aibkc"(mps[i], mpo[i])
        result[i] = reshape(tmp, size(mps[i],1)*size(mpo[i],1), :, 2)
    end
    return result
end

function Compress(mps::Array{T,S},chi::Int) where {T,S}
    residual=0
    len=length(mps)
    for i in 1:len 
        mps[i]=permutedims(mps[i],(1,3,2))
    end
    for i in 1:len-1
        F=qr(reshape(mps[i],size(mps[i],1)*2,:))
        Q=Matrix(F.Q)
        mps[i] = reshape(Q,size(mps[i],1),2,:)
        mps[i+1] = reshape(ein"ij,jab->iab"(F.R,mps[i+1]), size(F.R,1), 2, size(mps[i+1],3))  
        if mod(i,20) == 0
            tnorm = norm(mps[i+1])
            mps[i+1] /= tnorm
            residual += log(tnorm)
        end
    end   
    for i in len:-1:2
        F=safe_svd(reshape(ein"ijk,kab->ijab"(mps[i-1],mps[i]), size(mps[i-1],1)*2,size(mps[i],3)*2))     
        if size(F.S,1)>chi
            Vt=F.Vt[1:chi,:]; U=F.U[:,1:chi]
            mps[i]=reshape(Vt, : , 2, size(mps[i],3))
            mps[i-1]=reshape(U*Diagonal(F.S[1:chi]), size(mps[i-1],1),2,:)
        else 
            mps[i]=reshape(F.Vt, :, 2, size(mps[i],3)) 
            mps[i-1]=reshape(F.U*Diagonal(F.S), size(mps[i-1],1),2,:)            
        end   
    end       
    tnorm=norm(mps[1])
    mps[1] /= tnorm
    residual += log(tnorm)
    for i in 1:len
        mps[i]=permutedims(mps[i],(1,3,2))
    end
    return residual,mps
end

function Contruction(pa::Parament, te::Tensor, J::Array{Int}, mps::Array{T,S}) where {T,S}
    te.VerticalTensorsRi.Tensor[pa.L] = mps[pa.L]
    for i in pa.L:-1:2
        mps[i-1] = ein"ijk,k,lin->ljn"(mps[i],OutField(J[i],pa),mps[i-1])
        te.VerticalTensorsRi.Tensor[i-1] = mps[i-1]
    end
end


function OutField(ju::Int,pa::Parament)
    return ju == 1 ? pa.field :  ju == -1 ? pa.Antifield : ones(2)
end 

function OutConditionalZ(x::Int,y::Int,ju::Int,pa::Parament,te::Tensor)
    if x == 1
        return log.(abs.(te.VerticalTensorsRi.Tensor[1][:]))
    else
        tensor =  ein"ij,jkl->ikl"(te.VerticalTensorsLe.Tensor[x-1], te.VerticalTensorsRi.Tensor[x])
        return log.(abs.(tensor[:]))
    end
end


function GetlnZ(pa::Parament, sp::Spin, x::Int, y::Int, spin::Array{Int,2})
    lnZ = 0; nnb = 4
    jx=x+sp.Dx[nnb]; jy=y+sp.Dy[nnb]
    if 0<jx<=pa.L && 0<jy<=pa.L
        lnZ += spin[jy,jx] *  sp.bond[y, x, nnb] *pa.Beta
    end
    return [lnZ,-lnZ]
end


function ChangeVerticalTensorsLe(x::Int,y::Int,s::Int,pa::Parament,te::Tensor) 
    vecs = s == 1 ? [1,0] : [0,1] 
    C = ein"k,ijk->ij"(vecs,te.HorizontTensorsDw.Tensor[y,x])
    if x == 1
        te.VerticalTensorsLe.Tensor[x] = C
    else
        te.VerticalTensorsLe.Tensor[x] = te.VerticalTensorsLe.Tensor[x-1]*C
    end
end

"""
TNMH(pa::Parament, sp::Spin, te::Tensor)

This function performs the TNMH algorithm for a 2D RBIM system.

# Arguments
- `pa::Parament`: The parameters of the RBIM system.
- `sp::Spin`: The spin configuration of the RBIM system.
- `te::Tensor`: The tensor representation of the RBIM system.

# Returns
- `spin`: The updated spin configuration after the TNMH algorithm.
- `lnpr`: The logarithm of the acceptance probability.
- `lnZU[1]`: The logarithm of the partition function.

"""
function TNMH(pa::Parament, sp::Spin, te::Tensor)
    lnZt =  0
    #####################
    lnZpast = 0
    spin = zeros(Int, pa.L, pa.L)
    lnpr = 0; lnZq = zeros(2); lnZU = zeros(2)
    for y in 1:pa.L
        J = spin[max(y-1,1),:].*sp.bond[y,:,4]
        Contruction(pa,te,J,te.HorizontTensorsDw.Tensor[y,:])
        for x in 1:pa.L
            ZU = GetlnZ(pa,sp,x,y,spin)
            lnZU .+= ZU
            lnZq = OutConditionalZ(x,y,spin[y,max(x-1,1)]*sp.bond[y,x,1],pa,te)
            lnZq .+= lnZU
            if y!=pa.L
                lnZq .+=  te.HorizontTensorsDw.lnZ[y]
            end
            pi= 1.0/(1.0+exp(lnZq[2]-lnZq[1]))
    
            # if lnZq[2]+log(1+exp(lnZq[1]-lnZq[2]))-lnZpast > 1e-4
                # println(lnZq[2],' ',lnZq[1],' ',1.0/(1.0+exp(lnZq[2]-lnZq[1])),' ',lnZq[2]+log(1+exp(lnZq[1]-lnZq[2]))-lnZpast,' ',lnZq[2]+log(1+exp(lnZq[1]-lnZq[2])),' ',y,' ',x)
            # end

            ####--
            if rand()<pi 
                spin[y,x]=1
                lnZU[2] -= ZU[2]; lnZU[2] += ZU[1]
                ChangeVerticalTensorsLe(x,y,1,pa,te)                  
                lnZpast = lnZq[1]
                lnpr += log(pi)
            else
                spin[y,x]=-1
                lnZU[1] -= ZU[1]; lnZU[1] += ZU[2]
                ChangeVerticalTensorsLe(x,y,-1,pa,te)
                lnZpast = lnZq[2]
                lnpr += log(1-pi)
            end
        end
        for x in 1:pa.L-1
            jx=x+sp.Dx[3]
            if 0<jx<=pa.L 
                lnZU .+= spin[y,x] * spin[y,jx] *  sp.bond[y, x, 3] *pa.Beta 
            end
        end
    end
    return spin,lnpr,lnZU[1]
end

function Accept(pa::Parament,sp::Spin,st::Statistics,da::DataFile, spinold::Array{Int,2},spinnew::Array{Int,2},lnprold::Float64,lnprnew::Float64,Energyold::Float64,Energynew::Float64)
    deltaE = Energynew-Energyold
    spinnext = zeros(Int, pa.L, pa.L)
    pi = exp(deltaE+lnprold-lnprnew)
    # @show deltaE, lnprold-lnprnew
    st.Acc[2] += 1.0 
    if pi>=1
        spinnext .= spinnew
        prnext= lnprnew
        Ennext= Energynew
        st.Acc[1] += 1.0 
        print(da.temdoc["Acc"],  @sprintf("%8i",1), " ")
    else
        if rand()<pi
            spinnext .= spinnew
            prnext = lnprnew
            Ennext = Energynew
            st.Acc[1] += 1.0 
            print(da.temdoc["Acc"],  @sprintf("%8i",1), " ")
        else
            spinnext  .= spinold
            prnext= lnprold
            Ennext= Energyold
            print(da.temdoc["Acc"],  @sprintf("%8i",0), " ")
        end
    end
    return spinnext,prnext,Ennext
end 
