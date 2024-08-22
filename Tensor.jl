
using LinearAlgebra:svd, Diagonal, qr, norm
using TensorOperations
using OMEinsum
using BenchmarkTools
import Base.GC: gc
export InitTensor, TNMH


struct MpsState
    Tensor::Array{Any}  # Allow for different array types
    lnZ::Array{Float64, 1}
end

mutable struct Tensor
    StoreTensor::Array{Any,2}
    ChangTensor::Array{Any,2}
    HorizontTensorsDw::MpsState
    VerticalTensorsLe::MpsState
    VerticalTensorsRi::MpsState
    
    function Tensor(pa::Parament, StoreTensor::Array{Any,2}) 
        ChangTensor = similar(StoreTensor)
        HorizontTensorsDw = MpsState(Array{Any,2}(undef, pa.L, pa.L), zeros(Float64, pa.L))
        VerticalTensorsLe = MpsState(Array{Any,1}(undef, pa.L), zeros(Float64, pa.L))
        VerticalTensorsRi = MpsState(Array{Any,1}(undef, pa.L), zeros(Float64, pa.L))
        new(StoreTensor, ChangTensor, HorizontTensorsDw, VerticalTensorsLe, VerticalTensorsRi)
    end
end

export Tensor

function GetAntiPara(sp::Spin,nnb::Int, y::Int, x::Int)
    return sp.bond[y, x, nnb] == 0 ? 1 : (sp.bond[y, x, nnb] == 1 ? 2 : 3)
end

function OutTensorGrid(pa::Parament, sp::Spin, te::Tensor)
    # tensor=Array{Any, 2}(undef,pa.L,pa.L)
    for y in 1:pa.L, x in 1:pa.L
        if y == 1  || y == pa.L
            if x == 1
                A = reshape(pa.I2, 1, 2, 2)
            elseif x == pa.L
                A = reshape(pa.I2, 2, 1, 2)
            else  
                A = pa.I3
            end
        else         
            if x == 1
                A = reshape(pa.I3, 1, 2, 2, 2)
            elseif x == pa.L
                A = reshape(pa.I3, 2, 1, 2, 2)
            else  
                A = pa.I4
            end
        end
        if y == pa.L
            A = ein"ijk,lj->ilk"(A, pa.BondTensor[GetAntiPara(sp, 3, y, x)])
        elseif y == 1
            A = ein"ijk,lj,ka->ila"(A, pa.BondTensor[GetAntiPara(sp, 3, y, x)],
            pa.BondTensor[GetAntiPara(sp, 2, y, x)])         
        else
            A = ein"ijkl,nj,ka->inal"(A, pa.BondTensor[GetAntiPara(sp, 3, y, x)],
            pa.BondTensor[GetAntiPara(sp, 2, y, x)])  
        end   
        te.StoreTensor[y,x] = A
    end
end           


function InitTensor(pa::Parament, sp::Spin, te::Tensor)
    OutTensorGrid(pa, sp, te)
    StorelnZ(pa, te)
end


function StorelnZ(pa::Parament,te::Tensor)
    lnZ=0
    te.HorizontTensorsDw.lnZ[pa.L]=0 
    te.HorizontTensorsDw.Tensor[pa.L,:] .= te.StoreTensor[pa.L,:]
    for y in pa.L:-1:3
        # res,te.HorizontTensorsDw.Tensor[y-1,:] = Compress(Eat(pa,te.HorizontTensorsDw.Tensor[y,:] ,te.StoreTensor[y-1,:]),pa.chi)

        res,te.HorizontTensorsDw.Tensor[y-1,:] = Compress(Eat(pa,te.HorizontTensorsDw.Tensor[y,:] ,te.StoreTensor[y-1,:]),pa.chi)
        lnZ += res
        te.HorizontTensorsDw.lnZ[y-1] = lnZ
    end
end

function Contruction(pa::Parament, te::Tensor, MpsUp::Array{T,S}, MpsDw::Array{T,S}) where {T,S}
    lnZ = 0
    len = length(MpsUp)
    mps = Array{Any,1}(undef, len)
    mps[len] = ein"ikj,abj->iakb"(MpsUp[len],MpsDw[len])
    tnorm = norm(mps[pa.L])
    mps[pa.L] /= tnorm
    te.VerticalTensorsRi.Tensor[len] = mps[len]
    lnZ += log(tnorm)
    te.VerticalTensorsRi.lnZ[len] = lnZ
    for i in pa.L:-1:2
        mps[i-1] =  ein"abij,nak,lbk->nlij"(mps[i],MpsUp[i-1],MpsDw[i-1])
        tnorm = norm(mps[i-1])
        mps[i-1] /= tnorm
        te.VerticalTensorsRi.Tensor[i-1] = mps[i-1]
        lnZ += log(tnorm)
        te.VerticalTensorsRi.lnZ[i-1] = lnZ
        mps[i]=nothing
    end
end

function Contruction(pa::Parament, te::Tensor, mps::Array{T,S}) where {T,S}
    lnZ = 0
    tnorm = norm(mps[pa.L])
    mps[pa.L] /= tnorm
    te.VerticalTensorsRi.Tensor[pa.L] = mps[pa.L]
    lnZ += log(tnorm)
    te.VerticalTensorsRi.lnZ[pa.L] = lnZ
    for i in pa.L:-1:3
        mps[i-1] =  ein"ji,kj->ki"(mps[i], mps[i-1])
        tnorm = norm(mps[i-1])
        mps[i-1] /= tnorm
        te.VerticalTensorsRi.Tensor[i-1] = mps[i-1]
        lnZ += log(tnorm)
        te.VerticalTensorsRi.lnZ[i-1] = lnZ
        mps[i]=nothing
    end
end

function Eat(pa::Parament,mps::Array{T,S},mpo::Array{T,O}) where {T,S,O}
    return [reshape(ein"ikj,abjc->aibkc"(mps[i],mpo[i]),size(mps[i],1)*size(mpo[i],1),:,2) for i in 1:pa.L]
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
        F=svd(reshape(ein"ijk,kab->ijab"(mps[i-1],mps[i]), size(mps[i-1],1)*2,size(mps[i],3)*2))     
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

function OutField(J::Int,pa::Parament,s::Int)
    ju = J*s
    return ju == 1 ? pa.field : pa.Antifield
end 

function InputTensor(x::Int,y::Int,pa::Parament,sp::Spin,te::Tensor,spin::Int)
    B = nothing; C = nothing
    for nnb in [2,3]
        jx=x+sp.Dx[nnb]; jy=y+sp.Dy[nnb]
        if 0<jx<=pa.L && 0<jy<=pa.L
            Field = OutField(sp.bond[y,x,nnb], pa, spin)
            if nnb == 3
                if y == pa.L
                    B = reshape(ein"ij,i->j"(te.VerticalTensorsRi.Tensor[x+1] , Field), 1, :)
                else
                    B = ein"ijkl,i->jkl"(te.VerticalTensorsRi.Tensor[x+1], Field)
                    B = reshape(B, size(B,1), :)
                end
            else
                C = ein"ijl,l->ij"(te.HorizontTensorsDw.Tensor[jy,jx], Field)
                if x != 1
                    C = ein"ij,jk->ik"(te.VerticalTensorsLe.Tensor[x-1], C)
                end
            end
        end
    end
    return B, C
end

function GetlnZ(x::Int,y::Int,pa::Parament,B::T,C::O) where {T,O}
    lnZ = 0
    if y == pa.L
        if x == pa.L
            lnZ += 0
        else
            tnorm = norm(B)
            lnZ += log(tnorm) 
        end
    else
        if x == pa.L
            tnorm = norm(C)
            lnZ += log(tnorm)
        else
            End = ein"ij,jk->ik"(C, B)
            tnorm = norm(End) 
            lnZ += log(tnorm) 
        end
    end
    return lnZ
end

function GetlnZ(pa::Parament, sp::Spin, x::Int, y::Int, spin::Array{Int,2})
    lnZ = 0
    for nnb in [1,2,3,4]
        jx=x+sp.Dx[nnb]; jy=y+sp.Dy[nnb]
        if 0<jx<=pa.L && 0<jy<=pa.L
            lnZ += spin[jy,jx] *  sp.bond[y, x, nnb] *pa.Beta
        end
    end
    return lnZ, -lnZ
end

function ChangeGridTensor(y::Int,pa::Parament,sp::Spin,te::Tensor,spin::Array{Int,2})
    for x in 1:pa.L
        te.ChangTensor[y,x] = nothing 
        jx=x+sp.Dx[2]; jy=y+sp.Dy[2]
        if 0<jx<=pa.L && 0<jy<=pa.L
            Field = OutField(sp.bond[y,x,2], pa, spin[y,x])
            if jy != pa.L
                te.ChangTensor[jy,jx] = ein"ijkl,l->ijk"(te.StoreTensor[jy,jx], Field)
            else
                te.ChangTensor[jy,jx] = ein"ijl,l->ij"(te.StoreTensor[jy,jx], Field)     
            end
        end
    end
end

function ChangeVerticalTensorsLe(pa::Parament,x::Int,te::Tensor,C::T) where {T}
    tnorm = norm(C)
    C ./= tnorm
    te.VerticalTensorsLe.Tensor[x] = C
    if x == 1
        te.VerticalTensorsLe.lnZ[x] = log(tnorm)
    else
        te.VerticalTensorsLe.lnZ[x] = log(tnorm) + te.VerticalTensorsLe.lnZ[x-1]
    end
    if x != 1 te.VerticalTensorsLe.Tensor[x-1] = nothing end
    if x != pa.L te.VerticalTensorsRi.Tensor[x+1] = nothing end
end

function TNMH(pa::Parament,sp::Spin,te::Tensor)
    #####################
    # teCheck = CheckTensor(deepcopy(te.StoreTensor))
    lnZt =  0
    #####################
    lnZpast = 0
    spin = zeros(Int, pa.L, pa.L)
    lnpr = 0; lnZ1 = 0; lnZ2= 0; lnZup=0; lnZdw=0
    for y in 1:pa.L
        if y == 1 te.ChangTensor[1,:] .= te.StoreTensor[1,:]  end
        if y == pa.L
            Contruction(pa,te,te.ChangTensor[y,:])
        else
            Contruction(pa,te,te.ChangTensor[y,:],te.HorizontTensorsDw.Tensor[y+1,:])
        end
        for x in 1:pa.L
            ZU,ZD = GetlnZ(pa,sp,x,y,spin)
            lnZup+=ZU;   lnZdw+=ZD
            spin[y,x] = 1
            B,C1 = InputTensor(x,y,pa,sp,te,spin[y,x])
            lnZ1 = GetlnZ(x,y,pa,B,C1) + lnZup
            if x == 1 
                lnZ1 += te.VerticalTensorsRi.lnZ[x+1]
            elseif x == pa.L
                lnZ1 += te.VerticalTensorsLe.lnZ[x-1]
            else
                lnZ1 += te.VerticalTensorsLe.lnZ[x-1] + te.VerticalTensorsRi.lnZ[x+1]
            end
            if y!=pa.L
                lnZ1 +=  te.HorizontTensorsDw.lnZ[y+1]
            end
            ################
            # tensors1 = InputCheckTensor(pa,sp,teCheck,y,x,1)
            # lnZ1Check = GetChecklnZ(pa,y,x,tensors1)+ lnZup
            ###############            
            spin[y,x] = -1
            B,C2 = InputTensor(x,y,pa,sp,te,spin[y,x])
            lnZ2 = GetlnZ(x,y,pa,B,C2) + lnZdw 
            if x == 1 
                lnZ2 += te.VerticalTensorsRi.lnZ[x+1]
            elseif x == pa.L
                lnZ2 += te.VerticalTensorsLe.lnZ[x-1]
            else
                lnZ2 += te.VerticalTensorsLe.lnZ[x-1] + te.VerticalTensorsRi.lnZ[x+1]
            end
            if y!=pa.L
                lnZ2 +=  te.HorizontTensorsDw.lnZ[y+1]
            end
            ##############
            # tensors2 = InputCheckTensor(pa,sp,teCheck,y,x,-1)
            # lnZ2Check = GetChecklnZ(pa,y,x,tensors2)+ lnZdw
            # @show lnZ2, lnZ2Check
            #############
            pi=1.0/(1.0+exp(lnZ2-lnZ1))
            # pi=1.0/(1.0+exp(lnZ2Check-lnZ1Check))
            ####--
            # @show pi, 1.0/(1.0+exp(lnZ2Check-lnZ1Check)),  y, x
            if y == 1 && x == 1
                lnZt = log(exp(lnZ1)+exp(lnZ2))
            end
            println(lnZ1,' ',lnZ2,' ',pi,' ',lnZ2+log(1+exp(lnZ1-lnZ2))-lnZpast,' ',y,' ',x,' ',lnZup,' ',lnZdw)
            # open("1.txt", "a") do file 
            #     println(file,lnZ1,' ',lnZ2,' ',pi,' ',lnZ2+log(1+exp(lnZ1-lnZ2))-lnZpast,' ',y,' ',x)
            # end
            # if lnZ2+log(1+exp(lnZ1-lnZ2))-lnZpast>1e-4
            #     println(lnZ2+log(1+exp(lnZ1-lnZ2))-lnZpast,' ',y,' ',x)
            #     println(pi,' ',exp(lnZ2-lnZ1),' ',log(pi))
            # end
        
            ####--
            if rand()<pi 
                lnZdw -= ZD; lnZdw += ZU
                spin[y,x]=1
                if y != pa.L
                    ChangeVerticalTensorsLe(pa,x,te,C1)
                end
                ####Check
                # teCheck.StoreTensor[:,:] .= tensors1[:,:]
                lnZpast = lnZ1 
                #### 
                lnpr += log(pi)
            else
                lnZup -= ZU; lnZup += ZD
                if y != pa.L
                    ChangeVerticalTensorsLe(pa,x,te,C2)
                end
                ####Check
                # teCheck.StoreTensor[:,:] .= tensors2[:,:]
                lnZpast = lnZ2 
                #### 
                lnpr += log(1-pi)
            end
        end
        ####
        ####
        if y != pa.L
            ChangeGridTensor(y,pa,sp,te,spin)
        end
        if y == pa.L-1 
            te.VerticalTensorsLe.lnZ[:] .= 0
        end
    end
    return spin,lnpr
end

function Accept(pa::Parament,sp::Spin,st::Statistics,da::DataFile, spinold::Array{Int,2},spinnew::Array{Int,2},lnprold::Float64,lnprnew::Float64)
    Energyold = GetEnergy(pa,sp,spinold); Energynew = GetEnergy(pa,sp,spinnew)
    deltaE = Energynew-Energyold
    # @show -pa.Beta*deltaE, lnprold-lnprnew
    pi = exp(-pa.Beta*deltaE+lnprold-lnprnew)
    @show pi, -pa.Beta*deltaE, lnprold-lnprnew
    # @show Energyold, Energynew, lnprold, lnprnew
    error("ting")
    st.Acc[2] += 1.0 
    if pi>=1
        spinnext  = spinnew
        prnext= lnprnew
        st.Acc[1] += 1.0 
        print(da.temdoc["Acc"],  @sprintf("%8i",1), " ")
    else
        if rand()<pi
            spinnext = spinnew
            prnext = lnprnew
            st.Acc[1] += 1.0 
            print(da.temdoc["Acc"],  @sprintf("%8i",1), " ")
        else
            spinnext  = spinold
            prnext= lnprold
            print(da.temdoc["Acc"],  @sprintf("%8i",0), " ")
        end
    end
    return spinnext,prnext
end 


#Check
struct CheckTensor
    StoreTensor::Array{Any,2}
end

function InputCheckTensor(pa::Parament, sp::Spin, te::CheckTensor,  y::Int, x::Int, spin::Int)
    tensors = deepcopy(te.StoreTensor)
    for nnb in [2,3]
        jx=x+sp.Dx[nnb]; jy=y+sp.Dy[nnb]
        if 0<jx<=pa.L && 0<jy<=pa.L 
            Field = OutField(sp.bond[y,x,nnb], pa, spin)
            if nnb == 3
                if y == pa.L
                    tensors[jy,jx] = ein"ij,i->j"(tensors[jy,jx] , Field)
                    tensors[jy,jx] = reshape(tensors[jy,jx], 1, size(tensors[jy,jx],1))
                else
                    tensors[jy,jx] = ein"ijk,i->jk"(tensors[jy,jx], Field)
                    tensors[jy,jx] = reshape(tensors[jy,jx], 1, size(tensors[jy,jx],1), :)
                end
            else
                if y == pa.L-1
                    tensors[jy,jx] = ein"ijl,l->ij"(tensors[jy,jx], Field)
                else
                    tensors[jy,jx] = ein"ijkl,l->ijk"(tensors[jy,jx], Field)
                end                    
            end
        end
    end
    return tensors
end

function GetChecklnZ(pa::Parament, y::Int, x::Int, tensorsA::Array{Any,2})
    lnZ=0
    tensors = deepcopy(tensorsA)
    if y == pa.L-1
        for i in x+1:pa.L
            tensors[y+1,i] = reshape(ein"ijk,abk->iajb"(tensors[y,i], tensors[y+1,i]), size(tensors[y,i],1)*size(tensors[y+1,i],1),:)
        end        
    elseif y < pa.L-1
        for i in x+1:pa.L
            tensors[y+1,i] = reshape(ein"ijk,abck->iajbc"(tensors[y,i], tensors[y+1,i]), size(tensors[y,i],1)*size(tensors[y+1,i],1),:,size(tensors[y+1,i],4))
        end
    end
    for head in pa.L:-1:y+3
        tensors[head-1,:] = Eat(pa,tensors[head,:],tensors[head-1,:])
        lnZ += 0
    end
    if y < pa.L-1
        res = full_contraction(pa,tensors[y+1,:],tensors[y+2,:])
    elseif y == pa.L-1
        for i in pa.L-1:-1:1
            tensors[y+1,i] = ein"ij,jk->ik"(tensors[y+1,i],tensors[y+1,i+1])
        end
        res = log(norm(tensors[y+1,1]))
    else
        for i in x+1:pa.L-1
            tensors[y,i+1] = ein"ij,jk->ik"(tensors[y,i],tensors[y,i+1])
        end
        res = log(norm(tensors[y,pa.L]))
    end
    lnZ+=res
    return lnZ
end

function full_contraction(pa::Parament,mps1::Array{T,S},mps2::Array{T,S}) where {T,S}
    mps = [reshape(ein"ikj,abj->iakb"(mps1[i],mps2[i]),size(mps1[i],1)*size(mps2[i],1),:) for i in 1:pa.L]
    for i in 1:pa.L-1
        mps[i+1] = reshape(ein"ij,jk->ik"(mps[i],mps[i+1]),size(mps[i],1),size(mps[i+1],2))
    end
    scale=norm(mps[pa.L])
    return log(scale)
end
