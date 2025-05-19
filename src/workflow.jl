


import Pkg
Pkg.activate("..")  # src에서 상위 폴더의 Project.toml 사용
Pkg.instantiate()
include("DCON.jl")
using .DCON

ini = DCON.Input();
ini
gfile_path = "./../g000000.00000"
ini.ideal_n = [1]

dcon_result, δW = DCON.run_dcon(ini, gfile_path)
rdcon_result, Δ = DCON.run_rdcon(ini, gfile_path)
stride_result, Δ = DCON.run_stride(ini, gfile_path)
dcon, rdcon, stride, ballooning, _, _, _, _ = DCON.run_stability(ini, gfile_path)

ini.ballooning = true
ca1 = DCON.run_ballooning(ini, gfile_path)


dcon_result[1].ca1

dcon_result[1]
rdcon_result[1].q_rational

rdcon[1]

rdcon[1].Delta_prime[:,:,1]
typeof(rdcon[1])

#rdcon
#stride
#run
#plotting


# dcon_in = Fortran90Namelists.readnml(joinpath("./inputs","dcon.in"))
# dcon = parse_dcon_namelist(dcon_in)





gfile_path = "g000000.00000"
tmp_path = joinpath("/tmp", "dcon_$(uuid4())")

# tmp_path = "/home/aspire1019/Stability/"
mkpath(tmp_path)


for f in ["dcon.in", "equil.in", "rdcon.in", "stride.in", "vac.in"]
    cp(joinpath(inputdir, f), joinpath(tmp_path, f); force=true)
end


dcon_in = Fortran90Namelists.readnml(joinpath(inputdir,"dcon.in"))
dcon_struct = parse_dcon_namelist(dcon_in)
typeof(dcon_struct)
dcon_struct.DCON_CONTROL.nn = 2
dict_for_saving = struct_to_dict(dcon_struct)
Fortran90Namelists.writenml(joinpath(tmp_path,"dcon.in"),dict_for_saving)



equil_in = Fortran90Namelists.readnml(joinpath(inputdir,"equil.in"))
dcon_struct = parse_equil_namelist(equil_in)
dict_for_saving = struct_to_dict(dcon_struct)
Fortran90Namelists.writenml(joinpath(tmp_path,"equil.in"),dict_for_saving)


vac_in = Fortran90Namelists.readnml("vac.in")
dcon_struct = parse_vac_namelist(vac_in)
dict_for_saving = struct_to_dict(dcon_struct)
Fortran90Namelists.writenml("hello",dict_for_saving)


stride_in = Fortran90Namelists.readnml("./inputs/stride.in")
dcon_struct = parse_stride_namelist(stride_in)
dict_for_saving = struct_to_dict(dcon_struct)
Fortran90Namelists.writenml("hello",dict_for_saving)






for f in ["dcon", "rdcon", "stride"]
    cp(joinpath(execdir, f), joinpath(tmp_path, f); force=true)
end
cp(gfile_path, joinpath(tmp_path, "gfile"); force=true)

for f in ["dcon", "rdcon", "stride"]
    run(`chmod +x $(joinpath(tmp_path, f))`)
end


tmp_path = "/home/jmlmir/Stability/wth"

cd(tmp_path) do
    run_code("./dcon");
end
cd(tmp_path) do
    run_code("./rdcon");
end
cd(tmp_path) do
    run_code("./stride");
end

dcon_result = load_dcon_output_nc(joinpath(tmp_path,"dcon_output_n1.nc"))
rdcon_result = load_rdcon_output_nc(joinpath(tmp_path,"rdcon_output_n1.nc"))
stride_result = load_stride_output_nc(joinpath(tmp_path,"stride_output_n1.nc"))
typeof(dcon_result)



#deltaW
dcon_result.W_t_eigenvalue[1,1]
rdcon_result.Delta_prime[:,:,1]











plot(rdcon_result, log=true)

rm(tmp_path; force=true, recursive=true)

