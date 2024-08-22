
function NormPdf(pa::Parament, pdf::Pdf, idis::Int)
    nor = 1/pa.Nsample
    for i in 1:pdf.Npdf
        pdf.Quan[i,idis] *= nor
    end
end

function StatPdf(pa::Parament, pdf::Pdf, Nblck::Int, comm::MPI.Comm)
    eps = 1e-15; tol = 0.5
    MPI.Allreduce!(pdf.Quan, MPI.SUM, comm)
    for i in 1:pdf.Npdf
        pdf.Ave[i] = sum(pdf.Quan[i,:])/pa.Ndisorder
    end

    BlckSize = pa.Ndisorder÷Nblck
    if Nblck < pa.Ndisorder
        for i in 1:pdf.Npdf
            k0 = 1 
            for k in 1:Nblck
                pdf.Quan[i,k] = sum(pdf.Quan[i,k0:k0+BlckSize-1])/BlckSize
                k0 = k0 + BlckSize
            end
        end
    end

    nor = 1/Nblck
    for i in 1:pdf.Npdf
        devp=0; pdf.Dev[i]=0;
        for k in 1:Nblck
            devn = pdf.Quan[i,k]-pdf.Ave[i]
            pdf.Dev[i] += devn^2
        end
        pdf.Dev[i] *= nor
        pdf.Dev[i] = sqrt(pdf.Dev[i]/(Nblck-1))
    end
end  

function qEAF(pa::Parament, pdf::Pdf, data::Vector{Float64})
    # Find the left and right edges of the peak 
    idx=argmax(data)
    left_idx = idx
    right_idx = idx            
    while left_idx > 1 && data[left_idx] > data[idx]/2.0
        left_idx -= 1
    end       
    while right_idx < pdf.Npdf && data[right_idx] > data[idx]/2.0
        right_idx += 1
    end
    qEA = pdf.Qva[idx]
    half_widths = pdf.Qva[right_idx]-pdf.Qva[left_idx]
    return abs(qEA), abs(half_widths)
end

function cal_QEA(pa::Parament,st::Statistics,pdf::Pdf)
    qEA, Hwi = qEAF(pa,pdf,pdf.Ave[:])
    st.Ave["QeaD"] = qEA;  st.Ave["HwiD"] = Hwi
    for k in 1:st.Nblck
        qEA, Hwi = qEAF(pa,pdf,pdf.Quan[:,k])
        st.Quan["QeaD"][k] = qEA;  st.Quan["HwiD"][k] = Hwi
    end
end