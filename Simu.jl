using JLD2
using ProgressMeter
using Dates
using Profile

function SimulationMetro(pa::Parament,sp::Spin,st::Statistics, da::DataFile,comm::MPI.Comm)
    # @showprogress 1 "Processing" for i in 1:pa.Nblck
    start_time = now()
    st.Acc[:] .= 0.0
    rank = MPI.Comm_rank(comm)
    for idis in rank*(pa.Ndisorder÷pa.Nmpi)+1:(rank+1)*(pa.Ndisorder÷pa.Nmpi)
        InitBondFBC(pa,sp)

        spin = [rand([-1, 1]) for _ in 1:pa.L, _ in 1:pa.L, _ in 1:2]

        for i in 1:pa.Ntherm
            Metropolis(pa,sp,spin,pa.Beta)
        end

        for isamp in 1:pa.Nsample
            Acc = Metropolis(pa,sp,spin,pa.Beta)
            print(da.temdoc["Acc"],  @sprintf("%10.6f",Acc), " ")
            measure(pa,sp,st,da,idis,spin)
        end
        NormSample(pa,st,da,idis)
    end
    try
        MPI.Barrier(comm)
    catch err
        if rank % pa.Nmpi == 0
            println("Error occurred on process $rank: $err. Excluding this process.")
        end
        return
    end
    # stat_analy(pa,st,comm)
    # write2file(pa,st,da)

    @show st.Acc[1]/st.Acc[2], st.Nblck, pa.L, pa.chi
    end_time = now()
    elapsed_time = end_time - start_time
    println("Elapsed time (seconds): ", Dates.value(elapsed_time) / 1e3)
end

function SimulationTNMH(pa::Parament,sp::Spin,te::Tensor,st::Statistics,da::DataFile,comm::MPI.Comm)
    # @showprogress 1 "Processing" for i in 1:pa.Nblck
    start_time = now()
    st.Acc[:] .= 0.0
    rank = MPI.Comm_rank(comm)

    spinold = zeros(Int, pa.L, pa.L,pa.Nreplic); spinnew = zeros(Int, pa.L, pa.L,pa.Nreplic)
    lnprold = zeros(Float64, pa.Nreplic); lnprnew = zeros(Float64, pa.Nreplic)
    for idis in rank*(pa.Ndisorder÷pa.Nmpi)+1:(rank+1)*(pa.Ndisorder÷pa.Nmpi)        
        InitBondFBC(pa,sp)
        # time0 = now()
        InitTensor(pa,sp,te)
        # time00 = now()
        # @show pa.L, pa.chi
        # println("Pertime (seconds): ", Dates.value(time00-time0) / 1e3)

        for i in 1:pa.Nreplic
            spinold[:,:,i],lnprold[i]=TNMH(pa,sp,te)
        end

        for itherm in 1:pa.Ntherm
            for i in 1:pa.Nreplic
                # time1 = now()
                spinnew[:,:,i],lnprnew[i]=TNMH(pa,sp,te)                      
                # time2 = now()
                # timen = time2 - time1
                # println("Pertime (seconds): ", Dates.value(timen) / 1e3)
                spinold[:,:,i],lnprold[i]=Accept(pa,sp,st,da,spinold[:,:,i],spinnew[:,:,i],lnprold[i],lnprnew[i])  
            end
            if pa.L >= 96 gc() end
        end
        for isamp in 1:pa.Nsample
            for i in 1:pa.Nreplic
                time1 = now()
                spinnew[:,:,i],lnprnew[i]=TNMH(pa,sp,te)
                spinold[:,:,i],lnprold[i]=Accept(pa,sp,st,da,spinold[:,:,i],spinnew[:,:,i],lnprold[i],lnprnew[i])
                time2 = now()
                elapsed_seconds = Dates.value(time2-time1) / 1e3
                print(da.temdoc["time"],  @sprintf("%10.6f",elapsed_seconds), " ")
                print(da.temdoc["memory"],  @sprintf("%12.4f",round(get_memory_usage_mb(), digits=2)), " ") 
            end
            measure(pa,sp,st,da,idis,Metropolis_tensor(pa,sp,spinold[:,:,:],pa.Beta))
            if pa.L >= 96 gc() end
        end
        NormSample(pa,st,da,idis)
    end
    try
        MPI.Barrier(comm)
    catch err
        if rank % pa.Nmpi == 0
            println("Error occurred on process $rank: $err. Excluding this process.")
        end
        return
    end
    # stat_analy(pa,st,comm)
    # write2file(pa,st,da)
    closedoc(da)

    @show st.Acc[1]/st.Acc[2], st.Nblck, pa.L, pa.chi
    end_time = now()
    elapsed_time = end_time - start_time
    println("Elapsed time (seconds): ", Dates.value(elapsed_time) / 1e3)
end

# Convert memory usage to MB or GB
function bytes_to_human_readable(bytes)
    if bytes < 1024
        return "$bytes bytes"
    elseif bytes < 1024^2
        return "$(bytes / 1024) KB"
    elseif bytes < 1024^3
        return "$(bytes / 1024^2) MB"
    else
        return "$(bytes / 1024^3) GB"
    end
end

function get_memory_usage_mb()
    pid = getpid()  # 获取当前进程ID
    command = `ps -o rss= -p $pid`  # 构造获取常驻集大小的命令
    rss = read(command, String)  # 执行命令并读取输出
    rss = strip(rss)  # 去除前后空格
    memory_usage_kb = parse(Int, rss)  # 常驻集大小（KB）
    memory_usage_mb = memory_usage_kb / 1024  # 转换为 MB
    return memory_usage_mb
end

function total_memory_usage()
    total_memory = 0
    
    # 获取所有全局变量的内存占用
    for var in names(Main)
        total_memory += Base.summarysize(Main.eval(Symbol(var)))
    end
    
    # 获取当前函数的内存占用
    total_memory += Base.summarysize(stacktrace()[1].func)
    
    return total_memory
end
