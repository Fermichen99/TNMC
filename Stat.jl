include("DicDat.jl")
using Printf

listpa=["Ene1","Ene2","ChiE","Mag1","Mag2","Mag4","Magcor","Magsid","Magbul","ChiM","BinM","Qvr1","Qvr2","Qvr4","ChiQ","BinQ",
        "Qnk1","Qnk2","ChiK","BinK","CtdE","CtdM","CtdQ","CtdK","BtdM","BtdQ","ChdM","ChdQ","ChQQ","CtQQ","BiQQ"]
listPdf = ["Pdf","C4q","Lnk"]
listtem = ["Ene1","Qnk1","Mag1","Qvr1","Magcor","Magsid","Magbul","Acc","time","memory"]
single = String[]
complx = String[]
if "Pdf" in listPdf push!(listpa, "QeaD"); push!(listpa, "HwiD") end
# Iterate through listpa and classify keys into "single" and "complx"
for key in listpa
    if get(Dat, key, "") == "single"
        push!(single, key)
    elseif get(Dat, key, "") == "complx"
        push!(complx, key)
    end
end

mutable struct Statistics
    NObs_b::Int
    NObs_c::Int
    NObs::Int
    Quan::Dict{String, Vector{Float64}}
    Ave::Dict{String, Float64}
    Dev::Dict{String, Float64}
    Cor::Dict{String, Float64}
    Nblck::Int
    Acc::Array{Float64, 1}
    pdf::Union{Pdf, Nothing}  
    c4q::Union{C4q, Nothing}  
    lnk::Union{Lnk, Nothing}  

    RW::Array{Float64, 2}

    function Statistics(pa::Parament)
        NObs_b = length(single)
        NObs_c = length(complx)
        NObs = NObs_b + NObs_c 

        # Initialize dictionaries
        Quan = Dict{String, Vector{Float64}}()
        Ave = Dict{String, Float64}()
        Dev = Dict{String, Float64}()
        Cor = Dict{String, Float64}()

        for param in listpa
            Quan[param] = zeros(Float64, pa.Ndisorder)
            Ave[param] = 0.0
            Dev[param] = 0.0
            Cor[param] = 0.0
        end
        Nblck = pa.Ndisorder
        Acc = zeros(Float64, 2)

        pdf = nothing
        c4q = nothing
        lnk = nothing
        if "Pdf" in listPdf  pdf=Pdf(pa)  end
        if "C4q" in listPdf  c4q=C4q(pa)  end
        if "Lnk" in listPdf  lnk=Lnk(pa)  end

        RW = zeros(Float64, length(pl), 9)
        new(NObs_b, NObs_c, NObs, Quan, Ave, Dev, Cor, Nblck, Acc, pdf, c4q, lnk, RW)
    end
end
export Statistics


struct DataFile
    SpinGlass_dir :: String
    folder_name::String
    temdoc :: Dict{String, IO}
    Rdoc   :: Dict{Float64, IO}

    function DataFile(pa::Parament)
        script_dir = dirname(@__FILE__)
        SpinGlass_dir = joinpath(script_dir) # "..", ".."
        Pstr=@sprintf("%.5f",pa.P)
        Lstr=@sprintf("%.0f",pa.L)
        folder_name = joinpath(SpinGlass_dir, "dat", string("L=",Lstr) ,string("P=", Pstr), string("seed=",pa.seed))
        if !ispath(folder_name)
            mkpath(folder_name)
        end

        Rfloder = joinpath(folder_name, "RW")
        if !isdir(Rfloder)
            mkdir(Rfloder)
        end
        Rdoc = Dict{Float64, IO}()
        for p in pl
            Rdoc[p] = open(joinpath(Rfloder, "p=$(p)"), "a")
            print(Rdoc[p], @sprintf("%6s", "RW"), @sprintf("%8i",pa.L), " ", @sprintf("%10.6f",pa.T), " ", @sprintf("%10.6f",pa.P), " ", @sprintf("%6i", pa.Nsample), " ", @sprintf("%6i", pa.Ndisorder), " ", @sprintf("%6i", pa.Ntherm), " ", @sprintf("%6i", pa.chi), " ", @sprintf("%6i", pa.seed), " ", @sprintf("%10.6f", p), "\n")
        end

        # obsdoc = open(joinpath(folder_name, "$(pa.Dimension)Dobs"),"a") 
        # pdfdoc = open(joinpath(folder_name, "$(pa.Dimension)Dpdf"),"a")  
        # c4qdoc = open(joinpath(folder_name, "$(pa.Dimension)Dc4q"),"a") 
        # lnkdoc = open(joinpath(folder_name, "$(pa.Dimension)Dlnk"),"a")     
        
        temdoc = Dict{String, IO}()
        for para in listtem
            temdoc[para] = open(joinpath(folder_name, "$(pa.Dimension)D$(para)"), "a")
            print(temdoc[para], @sprintf("%6s", para), @sprintf("%8i",pa.L), " ", @sprintf("%10.6f",pa.T), " ", @sprintf("%10.6f",pa.P), " ", @sprintf("%6i", pa.Nsample), " ", @sprintf("%6i", pa.Ndisorder), " ", @sprintf("%6i", pa.Ntherm), " ", @sprintf("%6i", pa.chi), " ", @sprintf("%6i", pa.seed), "\n")
        end
        new(SpinGlass_dir, folder_name, temdoc, Rdoc)
    end
end
export DataFile


function measure(pa::Parament,sp::Spin,st::Statistics,da::DataFile,idis::Int,Spin::Array{Int,3}) 
    Ene = 0;  Mag = 0; RPh = zeros(Float64, length(pl), 9)
    Magcor = 0;  Magsid = 0; Magbul = 0; qtL = pa.L÷4+1; haL = pa.L-qtL+1
    nor = 1/pa.Nreplic
    for i in 1:pa.Nreplic
        E = abs(GetEnergy(pa,sp,Spin[:,:,i])/pa.Vol)
        M = sum(Spin[:,:,i])/pa.Vol
        Magcor = (Spin[1,1,i]+Spin[1,pa.L,i]+Spin[pa.L,1,i]+Spin[pa.L,pa.L,i])/4.0
        Magsid = (sum(Spin[qtL:haL,1,i])+sum(Spin[qtL:haL,pa.L,i])+sum(Spin[pa.L,qtL:haL,i])+sum(Spin[1,qtL:haL,i]))/(2.0*pa.L)
        Magbul =  sum(Spin[qtL:haL,qtL:haL,i])*4/pa.Vol
        RP = GetR(pa,sp,Spin[:,:,i],M)
        print(da.temdoc["Ene1"],  @sprintf("%10.6f",E), " ")
        print(da.temdoc["Mag1"],  @sprintf("%10.6f",M), " ")
        print(da.temdoc["Magcor"],  @sprintf("%10.6f",Magcor), " ")
        print(da.temdoc["Magsid"],  @sprintf("%10.6f",Magsid), " ")
        print(da.temdoc["Magbul"],  @sprintf("%10.6f",Magbul), " ")
        M = abs(M)
        Ene += E; Mag += M;
        for i in 1:length(pl)
            RPh[i,:] .+= RP[i][:]
        end
    end
    Ene *= nor;  Mag *= nor; RPh .*= nor
    st.Quan["Ene1"][idis] += Ene
    st.Quan["Mag1"][idis] += Mag
    st.RW .+= RPh

    if "Ene2" in listpa    
        st.Quan["Ene2"][idis] += sum(Ene.^2)*nor
    end        
    if "Ene4" in listpa    
        st.Quan["Ene4"][idis] += sum(Ene.^4)*nor
    end  
    if "Mag2" in listpa    
        st.Quan["Mag2"][idis] += sum(Mag.^2)*nor
    end  
    if "Mag4" in listpa    
        st.Quan["Mag4"][idis] += sum(Mag.^4)*nor
    end  
    n = 0
    for i in 1:pa.Nreplic-1
        for j in i+1:pa.Nreplic
            n += 1
            MeasureComposite(pa,st,da,idis,i,j,Spin)  
        end
    end
    for para in listtem
        flush(da.temdoc[para])
    end
end

function MeasureComposite(pa::Parament, st::Statistics, da::DataFile, idis::Int, i1::Int, i2::Int, spin::Array{Int,3})
    Overlap = spin[:,:,i1] .* spin[:,:,i2]
    Qlink = 0
    if "Qvr1" in listpa
        q = sum(Overlap)/pa.Vol
        st.Quan["Qvr1"][idis] += q  
        print(da.temdoc["Qvr1"],  @sprintf("%10.6f",q), " ")
    end
    if "Qvr2" in listpa    
        st.Quan["Qvr2"][idis] += q^2
    end  
    if "Qvr4" in listpa    
        st.Quan["Qvr4"][idis] += q^4
    end  
    if "Pdf" in listPdf
        PV = zeros(Float64, st.pdf.Npdf)
        for i in 1:st.pdf.Npdf
            PV[i] =  st.pdf.PdfCoe*exp(-pa.Vol*(st.pdf.Qva[i]-q)^2/2.0)
            st.pdf.Quan[i,idis] += PV[i]
        end
        if "C4q" in listPdf
            Qlink = C4qV(pa,st.c4q,Overlap,idis,PV)
        end
        if "Lnk" in listPdf
            Link(Qlink,st.lnk,idis,PV)
        end
    end 
    if "Qnk1" in listpa
        print(da.temdoc["Qnk1"],  @sprintf("%10.6f",Qlink), " ")
        st.Quan["Qnk1"][idis] += Qlink
    end
    if "Qnk2" in listpa    
        st.Quan["Qnk2"][idis] += Qlink^2
    end  
    if "Qnk4" in listpa    
        st.Quan["Qnk4"][idis] += Qlink^4
    end  
end

function NormSample(pa::Parament,st::Statistics,da::DataFile,idis::Int)
    nor = 1/pa.Nsample
    for key in single
        st.Quan[key][idis] *= nor
    end
    st.RW .*= nor
    for (i,p) in enumerate(pl)
        for j in 1:9
            print(da.Rdoc[p],  @sprintf("%10.6f",st.RW[i,j]), " ")
        end
        print(da.Rdoc[p], "\n")
        flush(da.Rdoc[p])
    end
    st.RW .= 0

    for para in listtem
        print(da.temdoc[para], "\n")
        flush(da.temdoc[para])
    end
    # if "Pdf" in listPdf NormPdf(pa,st.pdf,idis) end
    # if "C4q" in listPdf NormC4q(pa,st.c4q,idis) end
    # if "Lnk" in listPdf NormLnk(pa,st.lnk,idis) end
end

function substract_one_if_odd(number::Int)
    if number % 2 == 1  # Check if the number is odd
        number -= 1     # Add one to the number
    end
    return number
end

function stat_analy(pa::Parament,st::Statistics,comm::MPI.Comm)
    eps = 1e-15; tol = 0.5
    for key in single
        MPI.Allreduce!(st.Quan[key], MPI.SUM, comm)
        st.Ave[key] = sum(st.Quan[key])/pa.Ndisorder
    end

    Nblck = pa.Ndisorder; nor = 1/Nblck;
    # while true
    #     prt = true   
        for key in single
            devp=0; st.Cor[key]=0; st.Dev[key]=0;
            for k in 1:Nblck
                devn = st.Quan[key][k]-st.Ave[key]
                st.Dev[key] += devn^2
                st.Cor[key] += devn*devp
                devp = devn 
            end
            st.Dev[key] *= nor; st.Cor[key] *= nor 
            if st.Dev[key] > eps  st.Cor[key] /= st.Dev[key] end
            st.Dev[key] = sqrt(st.Dev[key]/(Nblck-1))
            if abs(st.Cor[key])>tol   prt = false end
        end
        st.Nblck = Nblck

        # if prt      break  end
        # if Nblck<=4  
        #   prt = false;     
        #   println("your data(Qvr1) has too big correlation, please check it")
        #   break 
        # end
        # Nblck = substract_one_if_odd(Nblck)
        # Nblck /= 2; nor = 1/Nblck; Nblck = Int(Nblck)
        # #-- coarsen blocking ------------------------------------------
        # for key in single
        #     k0 = 1 
        #     for k in 1:Nblck
        #         st.Quan[key][k] = (st.Quan[key][k0] + st.Quan[key][k0+1])*0.5
        #         k0 = k0 + 2
        #     end
        # end
    # end

    if "Pdf" in listPdf 
        StatPdf(pa, st.pdf, Nblck, comm) 
        cal_Obs_pdf(pa,st)
    end

    if "C4q" in listPdf StatC4q(pa,st.c4q,st.pdf,Nblck,comm) end
    if "Lnk" in listPdf StatLnk(pa,st.lnk,st.pdf,Nblck,comm) end

    cal_Obs_comp(pa,st)

    for key in complx
        devp=0; st.Cor[key]=0; st.Dev[key]=0;
        for k in 1:Nblck
            devn = st.Quan[key][k]-st.Ave[key]
            st.Dev[key] += devn^2
            st.Cor[key] += devn*devp
            devp = devn 
        end
        st.Dev[key] *= nor; st.Cor[key] *= nor 
        if st.Dev[key] > eps  st.Cor[key] /= st.Dev[key] end
        st.Dev[key] = sqrt(st.Dev[key]/(Nblck-1))
    end
end    

function substract_one_if_odd(number::Int)
    if number % 2 == 1  # Check if the number is odd
        number -= 1     # Add one to the number
    end
    return number
end

function cal_Obs_pdf(pa::Parament,st::Statistics) 
    if "QeaD" in complx cal_QEA(pa,st,st.pdf) end
end

function cal_Obs_comp(pa::Parament,st::Statistics) 
    if "ChiE" in complx cal_Chi(pa,st,"ChiE", "Ene1", "Ene2") end
    if "ChiM" in complx cal_Chi(pa,st,"ChiM", "Mag1", "Mag2") end
    if "ChiQ" in complx cal_Chi(pa,st,"ChiQ", "Qvr1", "Qvr2") end
    if "ChiK" in complx cal_Chi(pa,st,"ChiK", "Qnk1", "Qnk2") end 
    if "BinE" in complx cal_Q(pa,st,"BinE", "Ene2", "Ene4") end
    if "BinM" in complx cal_Q(pa,st,"BinM", "Mag2", "Mag4") end
    if "BinQ" in complx cal_Q(pa,st,"BinQ", "Qvr2", "Qvr4") end    
    if "BtdM" in complx cal_Qtd(pa,st,"BtdM", "Mag2", "Mag4") end
    if "BtdQ" in complx cal_Qtd(pa,st,"BtdQ", "Qvr2", "Qvr4") end 
    if "ChdM" in complx cal_Chd(pa,st,"ChdM", "Mag1", "Mag2") end
    if "ChdQ" in complx cal_Chd(pa,st,"ChdQ", "Qvr1", "Qvr2") end 
    if "CtdQ" in complx cal_Ctd(pa,st,"CtdQ", "Qvr1", "Qvr2") end 
    if "CtdE" in complx cal_Ctd(pa,st,"CtdE", "Ene1", "Ene2") end
    if "CtdM" in complx cal_Ctd(pa,st,"CtdM", "Mag1", "Mag2") end
    if "CtdK" in complx cal_Ctd(pa,st,"CtdK", "Qnk1", "Qnk2") end
    if "BinK" in complx cal_Qlink(pa,st,"BinK", "ChiK", "CtdK") end   
    if "ChQQ" in complx cal_Chi(pa,st,"ChQQ", "Qvr2", "Qvr4") end
    if "CtQQ" in complx cal_Ctd(pa,st,"CtQQ", "Qvr2", "Qvr4") end 
    if "BiQQ" in complx cal_Qlink(pa,st,"BiQQ", "CtQQ", "ChQQ") end   
end

function cal_Chi(pa::Parament,st::Statistics,jb::String,b1::String,b2::String)
    st.Ave[jb] = st.Ave[b2]-st.Ave[b1]^2
    for k in 1:st.Nblck
        st.Quan[jb][k] = st.Quan[b2][k]-st.Quan[b1][k]^2
    end
end

function cal_Q(pa::Parament,st::Statistics,jb::String,b1::String,b2::String)
    tmp = st.Ave[b1]^2 
    if abs(tmp)>1e-15  tmp = st.Ave[b2]/tmp end
    st.Ave[jb] = tmp
    for k in 1:st.Nblck
        tmp = st.Quan[b1][k]^2
        if abs(tmp)>1e-15  tmp = st.Quan[b2][k]/tmp end
        st.Quan[jb][k] = tmp
    end
end

function cal_Chd(pa::Parament,st::Statistics,jb::String,b1::String,b2::String)
    st.Ave[jb] = st.Ave[b2]-st.Ave[b1]
    for k in 1:st.Nblck
        st.Quan[jb][k] = st.Quan[b2][k]-st.Quan[b1][k]
    end
end

function cal_Ctd(pa::Parament,st::Statistics,jb::String,b1::String,b2::String)
    nor = 1/st.Nblck
    for i in 1:st.Nblck
        st.Quan[jb][i] = st.Quan[b2][i] - st.Quan[b1][i]^2
    end
    st.Ave[jb] = sum(st.Quan[jb][1:st.Nblck])*nor
end

function cal_Qtd(pa::Parament,st::Statistics,jb::String,b1::String,b2::String)
    nor = 1/st.Nblck
    for i in 1:st.Nblck
        st.Quan[jb][i] = st.Quan[b2][i] / st.Quan[b1][i]^2
    end
    st.Ave[jb] = sum(st.Quan[jb][1:st.Nblck])*nor
end

function cal_Qlink(pa::Parament,st::Statistics,jb::String,b1::String,b2::String)
    tmp = st.Ave[b2] 
    if abs(tmp)>1e-15  tmp = tmp/st.Ave[b1] end
    st.Ave[jb] = tmp
    for k in 1:st.Nblck
        tmp = st.Quan[b2][k]
        if abs(tmp)>1e-15  tmp = tmp/st.Quan[b1][k] end
        st.Quan[jb][k] = tmp
    end
end

function write2file(pa::Parament,st::Statistics,da::DataFile)
    print(da.obsdoc, pa.L, " ", pa.T, " ", pa.P, " ", pa.Nsample, " ", pa.Ndisorder, " ", pa.Ntherm, " ", pa.chi, " ", pa.seed, "\n")
    for key in listpa
        print(da.obsdoc, @sprintf("%-10s", key), @sprintf("%-20.15f", st.Ave[key]), @sprintf("%-20.15f", st.Dev[key]), @sprintf("%-20.15f", st.Cor[key]), "\n")
    end
    close(da.obsdoc)

    if "Pdf" in listPdf 
        print(da.pdfdoc, pa.L, " ", pa.T, " ", pa.P, " ", pa.Nsample, " ", pa.Ndisorder, " ", pa.Ntherm, " ", pa.chi, " ", pa.seed, "\n")
        for j in 1:st.pdf.Npdf
            print(da.pdfdoc, @sprintf("%-20.15f", st.pdf.Qva[j]), @sprintf("%-20.15f", st.pdf.Ave[j]), @sprintf("%-20.15f", st.pdf.Dev[j]), "\n")
        end
        close(da.pdfdoc)
    end

    if "C4q" in listPdf 
        print(da.c4qdoc, pa.L, " ", pa.T, " ", pa.P, " ", pa.Nsample, " ", pa.Ndisorder, " ", pa.Ntherm, " ", pa.chi, " ", pa.seed, "\n")
        for r in 1:st.c4q.lenr
            for q in 1:st.c4q.lenq
                print(da.c4qdoc, @sprintf("%d", st.c4q.r[r]), "  ", @sprintf("%-10.6f", st.c4q.q[q]), @sprintf("%-20.15f", st.c4q.Ave[r,q]), @sprintf("%-20.15f", st.c4q.Dev[r,q]), "\n")
            end
        end
        close(da.c4qdoc)
    end

    if "Lnk" in listPdf 
        print(da.lnkdoc, pa.L, " ", pa.T, " ", pa.P, " ", pa.Nsample, " ", pa.Ndisorder, " ", pa.Ntherm, " ", pa.chi, " ", pa.seed, "\n")
        for q in 1:st.lnk.lenq
            print(da.lnkdoc, @sprintf("%-10.6f", st.lnk.q[q]), @sprintf("%-20.15f", st.lnk.Ave[q]), @sprintf("%-20.15f", st.lnk.Dev[q]), "\n")
        end
        close(da.lnkdoc)
    end
    for para in listtem
        close(da.temdoc[para])
    end
    for p in pl
        close(da.Rdoc[p])
    end
end

function closedoc(da::DataFile)
    for para in listtem
        close(da.temdoc[para])
    end
    for p in pl
        close(da.Rdoc[p])
    end
end


export write2file




