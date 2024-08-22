
function NormLnk(pa::Parament, lnk::Lnk, idis::Int)
    nor = 1/pa.Nsample
    for q in 1:lnk.lenq
        lnk.Quan[q,idis] *= nor
    end
end

function StatLnk(pa::Parament, lnk::Lnk, pdf::Pdf, Nblck::Int, comm::MPI.Comm)
    eps = 1e-15; tol = 0.5
    MPI.Allreduce!(lnk.Quan, MPI.SUM, comm)
    for q in 1:lnk.lenq
        lnk.Ave[q] = sum(lnk.Quan[q,:])/pa.Ndisorder
    end
    iN = pdf.Npdf÷2
    for q in 1:lnk.lenq
        lnk.Ave[q] /= pdf.Ave[q+iN]
    end

    BlckSize = pa.Ndisorder÷Nblck
    if Nblck < pa.Ndisorder
        for q in 1:lnk.lenq
            k0 = 1 
            for k in 1:Nblck
                lnk.Quan[q,k] = sum(lnk.Quan[q,k0:k0+BlckSize-1])/BlckSize
                lnk.Quan[q,k] /= pdf.Quan[q+iN,k]
                k0 = k0 + BlckSize
            end
        end
    end

    nor = 1/Nblck
    for q in 1:lnk.lenq
        devp=0; lnk.Dev[q]=0;
        for k in 1:Nblck
            devn = lnk.Quan[q,k]-lnk.Ave[q]
            lnk.Dev[q] += devn^2
        end
        lnk.Dev[q] *= nor
        lnk.Dev[q] = sqrt(lnk.Dev[q]/(Nblck-1))
    end
end  

function Link(Qlink::Float64, lnk::Lnk, idis::Int, PV::Vector{Float64})
    iN = lnk.lenq-1
    for q in 1:lnk.lenq
        lnk.Quan[q, idis] += Qlink*PV[q+iN]
    end
end
