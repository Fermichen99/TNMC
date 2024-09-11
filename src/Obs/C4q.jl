

function NormC4q(pa::Parament, c4q::C4q, idis::Int)
    nor = 1/pa.Nsample
    for r in 1:c4q.lenr
        for q in 1:c4q.lenq
            c4q.Quan[r,q,idis] *= nor
        end
    end
end

function StatC4q(pa::Parament, c4q::C4q, pdf::Pdf, Nblck::Int, comm::MPI.Comm)
    eps = 1e-15; tol = 0.5
    MPI.Allreduce!(c4q.Quan, MPI.SUM, comm)
    for r in 1:c4q.lenr
        for q in 1:c4q.lenq
            c4q.Ave[r,q] = sum(c4q.Quan[r,q,:])/pa.Ndisorder
        end
    end
    iN = pdf.Npdf÷2
    for r in 1:c4q.lenr
        for q in 1:c4q.lenq
            c4q.Ave[r,q] /= pdf.Ave[q+iN]
        end
    end

    BlckSize = pa.Ndisorder÷Nblck
    if Nblck < pa.Ndisorder
        for r in 1:c4q.lenr
            for q in 1:c4q.lenq
                k0 = 1 
                for k in 1:Nblck
                    c4q.Quan[r,q,k] = sum(c4q.Quan[r,q,k0:k0+BlckSize-1])/BlckSize
                    c4q.Quan[r,q,k] /= pdf.Quan[q+iN,k]
                    k0 = k0 + BlckSize
                end
            end
        end
    end

    nor = 1/Nblck
    for r in 1:c4q.lenr
        for q in 1:c4q.lenq
            devp=0; c4q.Dev[r,q]=0;
            for k in 1:Nblck
                devn = c4q.Quan[r,q,k]-c4q.Ave[r,q]
                c4q.Dev[r,q] += devn^2
            end
            c4q.Dev[r,q] *= nor
            c4q.Dev[r,q] = sqrt(c4q.Dev[r,q]/(Nblck-1))
        end
    end
end  

function C4qV(pa::Parament, c4q::C4q, q::Array{Int,2}, idis::Int, PV::Vector{Float64})
    indexr = 0; iN = c4q.lenq-1; link = 0
    for r in c4q.r
        sum = 0; indexr += 1 
        for y in 1:pa.L, x in 1:pa.L
            before = (0 < (y - r) <= pa.L) ? q[y-r, x] : 0
            behind = (0 < (y + r) <= pa.L) ? q[r+r, x] : 0
            left = (0 < (x - r) <= pa.L) ? q[y, x-r] : 0
            right = (0 < (x + r) <= pa.L) ? q[y, x+r] : 0
            sum += q[y,x]*(left+right+before+behind)
        end
        if r==1 link = sum/(pa.L*(pa.L-1)) end
        for q in 1:c4q.lenq
            c4q.Quan[indexr, q, idis] += sum/pa.Vol*PV[q+iN]
        end
    end
    return link
end
