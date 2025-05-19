module DCON

using Fortran90Namelists
using OrderedCollections
using NCDatasets
using UUIDs
using FilePathsBase
using FilePaths
using Printf
using AbstractTrees
export Input, run_dcon, run_rdcon, run_stride, run_stability

include("plottings.jl")
include("runner.jl")
include("namelist_tree.jl")


execdir = joinpath(@__DIR__,"..", "executables")
inputdir = joinpath(@__DIR__, "..","inputs")


function Input()
    check_cpu()
    println(inputdir)
    dcon_in = Fortran90Namelists.readnml(joinpath(inputdir,"dcon.in"));
    rdcon_in = Fortran90Namelists.readnml(joinpath(inputdir,"rdcon.in"));
    stride_in = Fortran90Namelists.readnml(joinpath(inputdir,"stride.in"));
    vac_in = Fortran90Namelists.readnml(joinpath(inputdir,"vac.in"));
    equil_in = Fortran90Namelists.readnml(joinpath(inputdir,"equil.in"));
    

    return InputFiles(
        dcon = parse_dcon_namelist(dcon_in),
        rdcon = parse_rdcon_namelist(rdcon_in),
        stride = parse_stride_namelist(stride_in),
        vac = parse_vac_namelist(vac_in),
        equil = parse_equil_namelist(equil_in),
        ideal_n = [1, 2, 3],
        resistive_n = Int64[1],
        ballooning = true,
        directory = nothing,
        remove = true
    )
end

function check_cpu()
    if Sys.ARCH == :x86_64 && Sys.islinux()
        @info "DCON.jl has binary file for your CPU architecture" 
    else
        error("DCON.jl does not have a binary for your CPU architecture. Sorry.")
    end
end


function check_directory(ini::InputFiles)
    if ini.directory === nothing
        return joinpath("/tmp", "dcon_$(uuid4())")
    else
        if !isdir(ini.directory)
            error("Directory :  $(ini.directory) Does not exist")
        end
        return ini.directory
    end
end

function remove_directory(tmp_path::String)
    if isdir(tmp_path)
        rm(tmp_path; force=true, recursive=true)
        @info "Removed Directory !" path=tmp_path
    else
        @warn "No such directory" path=tmp_path
    end
end

function save_input(ini::InputFiles; tmp_path=nothing , save_list::Vector{String} = ["dcon", "equil", "vac", "stride", "rdcon"])
    
    if tmp_path === nothing 
        tmp_path = check_directory(ini)
        @info "Directory is not specified. Create new working directory" tmp_path
    end

    if !isdir(tmp_path)
        mkpath(tmp_path)
    end

    # 전체 가능한 구조체 목록
    all_namelists = Dict(
        "dcon"   => (filename = "dcon.in",   structure = ini.dcon),
        "equil"  => (filename = "equil.in",  structure = ini.equil),
        "vac"    => (filename = "vac.in",    structure = ini.vac),
        "stride" => (filename = "stride.in", structure = ini.stride),
        "rdcon"  => (filename = "rdcon.in",  structure = ini.rdcon)
    )

    # 선택된 목록만 저장
    for name in save_list
        if haskey(all_namelists, name)
            fname, structure = all_namelists[name]
            dict_for_saving = struct_to_dict(structure)
            Fortran90Namelists.writenml(joinpath(tmp_path, fname), dict_for_saving)
            run(`chmod +x $(joinpath(tmp_path, fname))`)
        else
            @warn "Requested namelist \"$name\" is not recognized"
        end
    end
end


function save_executables(ini::InputFiles; tmp_path=nothing , save_list::Vector{String} = ["dcon", "equil", "vac", "stride", "rdcon"])
    for f in save_list
        cp(joinpath(execdir, f), joinpath(tmp_path, f); force=true)
    end
end


function run_dcon(ini::InputFiles, gfile::String)
    tmp_dir = check_directory(ini)
    println("tmp_dir은 :", tmp_dir)
    save_input(ini,tmp_path = tmp_dir,save_list = ["dcon", "equil", "vac"])
    save_executables(ini,tmp_path = tmp_dir, save_list = ["dcon"])
    cp(gfile, joinpath(tmp_dir, "gfile"); force=true)

    δW = Dict{Int, Union{Float64, Nothing}}()
    dcon_results = Dict{Int, Union{DCONOutputNetCDF, Nothing}}()  # 결과는 DCONOutputNetCDF 또는 nothing
    
    if ini.ideal_n === nothing
        ini.ideal_n = [ini.dcon.DCON_CONTROL.nn]
    end

    for n in ini.ideal_n
        try
            ini.dcon.DCON_CONTROL.nn = n
            save_input(ini; tmp_path = tmp_dir, save_list = ["dcon"])
            @info "Running DCON for n = $n"

            cd(tmp_dir) do
                run_code("./dcon")
            end


            result_file = joinpath(tmp_dir, "dcon_output_n$(n).nc")
            result = load_dcon_output_nc(result_file)
            
            dcon_results[n] = result
            δW[n] = result.W_t_eigenvalue[1,1]

        catch e
            @warn "DCON failed for n=$n" exception = e
            dcon_results[n] = nothing
            δW[n] = nothing
        end
    end
    
    if ini.remove 
        remove_directory(tmp_dir)
    end

    return dcon_results, δW
end



function run_rdcon(ini::InputFiles, gfile::String)
    tmp_dir = check_directory(ini)
    println("tmp_dir은 :", tmp_dir)
    save_input(ini,tmp_path = tmp_dir,save_list = ["rdcon","equil", "vac"])
    save_executables(ini,tmp_path = tmp_dir, save_list = ["rdcon"])
    cp(gfile, joinpath(tmp_dir, "gfile"); force=true)
    
    rdcon_results = Dict{Int, Union{RdconOutputNetCDF, Nothing}}()
    delta_prime = Dict{Int, Union{Float64, Nothing}}()
    
    if ini.resistive_n === nothing
        ini.resistive_n = [ini.rdcon.RDCON_CONTROL.nn]
    end

    for n in ini.resistive_n
        try
            ini.rdcon.RDCON_CONTROL.nn = n
            save_input(ini; tmp_path = tmp_dir, save_list = ["rdcon"])
            @info "Running RDCON for n = $n"

            cd(tmp_dir) do
                run_code("./rdcon")
            end


            result_file = joinpath(tmp_dir, "rdcon_output_n$(n).nc")
            result = load_rdcon_output_nc(result_file)
            
            rdcon_results[n] = result
            delta_prime[n] = result.Delta_prime[1,1,1]

        catch e
            @warn "RDCON failed for n=$n" exception = e
            rdcon_results[n] = nothing
            delta_prime[n] = nothing
        end
    end
    
    if ini.remove 
        remove_directory(tmp_dir)
    end

    return rdcon_results, delta_prime
end

function run_stride(ini::InputFiles, gfile::String)

    tmp_dir = check_directory(ini)
    println("tmp_dir은 :", tmp_dir)
    save_input(ini,tmp_path = tmp_dir,save_list = ["stride","equil", "vac"])
    save_executables(ini,tmp_path = tmp_dir, save_list = ["stride"])
    cp(gfile, joinpath(tmp_dir, "gfile"); force=true)
    
    stride_results = Dict{Int, Union{StrideOutputNetCDF, Nothing}}()
    delta_prime = Dict{Int, Union{Float64, Nothing}}()
    
    if ini.resistive_n === nothing
        ini.resistive_n = [ini.stride.stride_control.nn]
    end

    for n in ini.resistive_n
        try
            ini.stride.stride_control.nn = n
            save_input(ini; tmp_path = tmp_dir, save_list = ["stride"])
            @info "Running STRIDE for n = $n"

            cd(tmp_dir) do
                run_code("./stride")
            end


            result_file = joinpath(tmp_dir, "stride_output_n$(n).nc")
            result = load_stride_output_nc(result_file)
            
            stride_results[n] = result
            delta_prime[n] = result.Delta_prime[1,1,1]

        catch e
            @warn "STRIDE failed for n=$n" exception = e
            stride_results[n] = nothing
            delta_prime[n] = nothing
        end
    end
    
    if ini.remove 
        remove_directory(tmp_dir)
    end

    return stride_results, delta_prime
end


function run_ballooning(ini::InputFiles, gfile::String)

    if ini.ballooning
        @info "Ballooning Flag is on!"
    else
        @warn("Ballooning Flag is off! turning on ballooning flag")
    end

    ini.dcon.DCON_CONTROL.bal_flag = true
    ini.dcon.DCON_CONTROL.nn = 1
    tmp_dir = check_directory(ini)
    println("tmp_dir은 :", tmp_dir)
    save_input(ini,tmp_path = tmp_dir,save_list = ["dcon", "equil", "vac"])
    save_executables(ini,tmp_path = tmp_dir, save_list = ["dcon"])
    cp(gfile, joinpath(tmp_dir, "gfile"); force=true)
    
    ca1 = nothing

    try
        save_input(ini; tmp_path = tmp_dir, save_list = ["dcon"])
        @info "Running DCON for n = ∞ ( ballooning mode)"

        cd(tmp_dir) do
            run_code("./dcon")
        end

        result_file = joinpath(tmp_dir, "dcon_output_n1.nc")
        result = load_dcon_output_nc(result_file)
        ca1 = result.ca1
    catch e
        @warn "DCON failed for  n = ∞ ( ballooning mode)" exception = e
        ca1 = nothing
    end

    if ini.remove 
        remove_directory(tmp_dir)
    end

    ini.dcon.DCON_CONTROL.bal_flag = false
    is_all_positive = all(≥(0.0), ca1)

    return ca1, is_all_positive

end





function run_stability(ini::InputFiles, gfile::String ; plot_flag = true)

    dcon_result, dcon_δW = run_dcon(ini, gfile)
    rdcon_result, rdcon_delta_prime = run_rdcon(ini, gfile)
    stride_result, stride_delta_prime = run_stride(ini, gfile)
    ballooning_result = nothing
    ballooning_stable = false

    if ini.ballooning
        ballooning_result, ballooning_stable = run_ballooning(ini, gfile)
    end

    if plot_flag
        for (n, result) in rdcon_result
            if result !== nothing
                plot(result; log=true)
            else
                @warn "No result for n = $n, skipping plot."
            end
        end

        for (n, result) in stride_result
            if result !== nothing
                plot(result; log=true)
            else
                @warn "No result for n = $n, skipping plot."
            end
        end
    end

    println("\n========== Stability Analysis Summary ==========\n")

    # DCON results
    println("🔹 DCON ΔW")
    for (mode, val) in dcon_δW
        println("  - Mode (n = $mode): δW = $(round(val, sigdigits=5))")
    end

    # RDCON results
    println("\n🔹 RDCON Δ′")
    for (mode, val) in rdcon_delta_prime
        println("  - Mode (n = $mode): Δ′ = $(round(val, sigdigits=5))")
    end

    # STRIDE results
    println("\n🔹 STRIDE Δ′")
    for (mode, val) in stride_delta_prime
        println("  - Mode (n = $mode): Δ′ = $(round(val, sigdigits=5))")
    end

    # Ballooning 결과
    if ini.ballooning
        println("\n🔹 Ballooning Stability")
        status = ballooning_stable ? "Stable ✅" : "Unstable ❌"
        println("  - Ballooning Mode (n = ∞): $status")
    end

    return dcon_result,rdcon_result, stride_result, ballooning_result, dcon_δW, rdcon_delta_prime, stride_delta_prime, ballooning_stable
end



end