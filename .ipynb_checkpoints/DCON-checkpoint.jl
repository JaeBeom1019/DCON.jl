using Revise
include("DCON_file_IO.jl")
input = Dict()


for arr in dcon.input_file_list
    dcon, comments  = read_dcon_input("$arr")
end

input["dcon.in"]["DCON_OUTPUT"]["out_bal1"]
input
write_dcon_input("dcon.in",input["dcon.in"])




# 확인
println(input["dcon"]["DCON_OUTPUT"]["bin_bal2"])  # f (false)
println(input["dcon"]["DCON_CONTROL"]["nn"])       # 1 (정수)
println(input["dcon"]["DCON_CONTROL"]["psiedge"])  # 0.988 (실수)

# 사용법
write_dcon_input("dcon_new.in", input["dcon"])  # 새로운 파일로 저장

include("galsol.jl")

directory = "."  # 파일이 위치한 디렉토리 경로
galsol = GalsolData(directory)
for (psi, value) in galsol.left[1][-12]
    println("ψ = $psi, 값 = $value")
end

dcon = DCON()
load_dcon_input(dcon)
dcon.gfile="/home/aspire1019/code/etc/OpenFUSIONToolkit/src/examples/TokaMaker/DIIID/g192185.02440"
DCON.load_gfile(dcon)
run(`bash -c "dcon"`)





#=
자료구조
DCON - gfile, dcon_input_dir, input_file_list, δW, Δ'
=#


# 사용 예시
input = Dict{String, DconInput}()
file_path = "rdcon.in"
input["rdcon.in"], comments = read_dcon_input(file_path)

# 값 접근 테스트
println("sing1_flag 값: ", input["rdcon.in"]["RDCON_CONTROL"].sing1_flag)
println("gal_tol 값: ", input["rdcon.in"]["GAL_INPUT"].gal_tol)

# 파일 쓰기
write_dcon_input("rdcon_new.in", input, comments)