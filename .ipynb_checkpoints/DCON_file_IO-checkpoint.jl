using Printf

export DCON, load_dcon_input, load_gfile

mutable struct DCON
    dcon_input_dir::String
    input_file_list::Vector{String}
    gfile::String
    Δ::Matrix{Float64}
    _Δ::Float64

    #dcon output
    mpsi::Int
    mtheta::Int
    psilow::Float64
    psihigh::Float64
    rs_right::Float64
    rs_left::Float64
    zs_top::Float64
    zs_bot::Float64
    amean::Float64
    rmean::Float64
    aratio::Float64
    kappa::Float64
    delta1::Float64
    delta2::Float64
    volume::Float64
    li1::Float64
    li2::Float64
    li3::Float64
    ro::Float64
    zo::Float64
    psio::Float64
    betap1::Float64
    betap2::Float64
    betat::Float64
    betan::Float64
    betaj::Float64
    bwall::Float64
    bt0::Float64
    q0::Float64
    qmin::Float64
    q95::Float64
    qmax::Float64
    qa::Float64
    crnt::Float64
    I_over_aB::Float64
    δW::Vector{Dict{String, Float64}}
    _δW::Float64


    function DCON(dcon_input_dir="/home/aspire1019/code/DCON", 
                input_file_list=["dcon.in", "vac.in", "equil.in", "rdcon.in"])
        
        obj = new(dcon_input_dir, input_file_list, "")
        try
            println("💪 DCON/RDCON 모듈이 load 되었습니다.")
            #환경변수 julia와 연결
            env_output = read(`bash -c "module load gpec/dev && env"`, String)
            for line in split(env_output, "\n")
                parts = split(line, "=", limit=2)
                if length(parts) == 2
                    ENV[parts[1]] = parts[2]
                end
            end
        catch e
            println("⚠️ 모듈을 로드하는데 실패했습니다: ", e)
        end
        println("📂 Working directory:", pwd())
        g_files = filter(f -> startswith(f, "g"), readdir(pwd()))
        obj.gfile = !isempty(g_files) ? g_files[1] : ""
        return obj
    end
end


"""
📂 DCON 입력 파일을 GPEC input directory에서 가져옵니다.
"""
function load_dcon_input(path::Union{DCON,String})
    if typeof(path) == DCON
        source_dir = path.dcon_input_dir
    else
        #그대로 PATH 사용하기
        source_dir = path
    end

    destination_dir = pwd()  # 현재 디렉토리
    input_list = path.input_file_list

    println("📂 Source directory: ", source_dir)
    println("📂 Destination directory: ", destination_dir)
    # source_dir 안의 모든 파일 목록 가져오기
    for file in readdir(source_dir)
        if file in input_list
            try
                full_source_path = joinpath(source_dir, file)
                full_dest_path = joinpath(destination_dir, file)
    
                if isfile(full_source_path)
                    cp(full_source_path, full_dest_path; force=true)
                    println("✅ $(file) 을 성공적으로 이동했습니다.")
                else
                    println("🚨 $(file) 은 존재하지 않습니다.")
                end
            catch e
                println("⚠️ $(file)을 가져올 수 없습니다. 오류: ", e)
            end
        end
    end

end


"""
Geqdsk를 load합니다.
DCON을 넣었을 경우에는 DCON.gfile이 load 되며,
아무것도 안 적었을 시, 현재 디렉토리에서 g로 시작하는 파일이 gfile로 인식됩니다.

# 예제
``` julia
julia > load_gfile("절대경로")
julia > load_gfile("dcon::DCON")
julia > load_gfile()

"""
function load_gfile(path=nothing)

    g_file = nothing  # 기본값 설정

    if isa(path, DCON)  # path가 DCON 객체일 경우
        g_file = path.gfile
        println("🗂️ `dcon.gfile` 사용: ", g_file)
    elseif isa(path, String)  # path가 문자열(경로)일 경우
        g_file = path
    else  # path가 아무것도 없을 경우, 현재 디렉토리에서 g로 시작하는 파일 탐색
        g_files = filter(f -> startswith(f, "g"), readdir(pwd()))
        g_file = !isempty(g_files) ? g_files[1] : error("🚨 `g`로 시작하는 파일이 없습니다.")
        println("✅ `g` 파일 자동 탐색: ", g_file)
    end

    println("✅ gfile의 경로는 : ", g_file)

end


function write_equil_in(dcon::DCON)

    #equil.in 수정하기
    lines = readlines("equil.in")

    # eq_filename 수정
    updated_lines = map(lines) do line
        if occursin("eq_filename=", line)
            return """    eq_filename="$(dcon.gfile)" ! Path to input file (geqdsk, etc.) Must specified (no default)."""
        else
            return line
        end
    end

    # 변경된 내용 저장
    open(joinpath(pwd(),"equil.in"), "w") do f
        for line in updated_lines
            println(f, line)
        end
    end

    println("✅ `equil.in`의 gfile 경로가 $(dcon.gfile) 로 업데이트되었습니다!")
end


function parse_dcon_file(dcon::DCON)

    filename = "dcon.out"
    lines = readlines(filename)
    
    # 1D 데이터 파싱
    mpsi_data = split(lines[8])
    shape_data = split(lines[13])
    li_data = split(lines[18])
    beta_data = split(lines[23])
    q_data = split(lines[28])
    
    # Total Energy Eigenvalues 파싱
    eig_values = Vector{Dict{String, Float64}}()
    start_reading = false
    for line in lines
        if startswith(line, " Total Energy Eigenvalues:")
            start_reading = true
            continue
        elseif startswith(line, "   isol") || isempty(strip(line))
            continue
        elseif start_reading
            if occursin("E", line)
                vals = split(strip(line))
                if length(vals) >= 4
                    push!(eig_values, Dict(
                        "isol" => parse(Float64, vals[1]),
                        "plasma" => parse(Float64, vals[2]),
                        "vacuum" => parse(Float64, vals[3]),
                        "total" => parse(Float64, vals[4])
                    ))
                else
                    break
                end
            else
                break
            end
        end
    end

    dcon.mpsi = parse(Int, mpsi_data[1])
    dcon.mtheta = parse(Int, mpsi_data[2])
    dcon.psilow = parse(Float64, mpsi_data[3])
    dcon.psihigh = parse(Float64, mpsi_data[4])
    dcon.rs_right = parse(Float64, mpsi_data[5])
    dcon.rs_left = parse(Float64, mpsi_data[6])
    dcon.zs_top = parse(Float64, mpsi_data[7])
    dcon.zs_bot = parse(Float64, mpsi_data[8])
    dcon.amean = parse(Float64, shape_data[1])
    dcon.rmean = parse(Float64, shape_data[2])
    dcon.aratio = parse(Float64, shape_data[3])
    dcon.kappa = parse(Float64, shape_data[4])
    dcon.delta1 = parse(Float64, shape_data[5])
    dcon.delta2 = parse(Float64, shape_data[6])
    dcon.volume = parse(Float64, shape_data[7])
    dcon.li1 = parse(Float64, li_data[1])
    dcon.li2 = parse(Float64, li_data[2])
    dcon.li3 = parse(Float64, li_data[3])
    dcon.ro = parse(Float64, li_data[4])
    dcon.zo = parse(Float64, li_data[5])
    dcon.psio = parse(Float64, li_data[6])
    dcon.betap1 = parse(Float64, beta_data[1])
    dcon.betap2 = parse(Float64, beta_data[2])
    dcon.betat = parse(Float64, beta_data[3])
    dcon.betan = parse(Float64, beta_data[4])
    dcon.betaj = parse(Float64, beta_data[5])
    dcon.bwall = parse(Float64, beta_data[6])
    dcon.bt0 = parse(Float64, beta_data[7])
    dcon.q0 = parse(Float64, q_data[1])
    dcon.qmin = parse(Float64, q_data[2])
    dcon.q95 = parse(Float64, q_data[3])
    dcon.qmax = parse(Float64, q_data[4])
    dcon.qa = parse(Float64, q_data[5])
    dcon.crnt = parse(Float64, q_data[6])
    dcon.I_over_aB = parse(Float64, q_data[7])
    dcon.δW = eig_values
    dcon._δW = minimum(dcon.δW[i]["total"] for i in eachindex(dcon.δW))
    
    return dcon
end


function parse_delta_gw_file(filename::String)
    lines = readlines(filename)
    
    # Delta matrix 데이터를 저장할 7x7 배열 초기화
    arr = Array{Complex{Float64}, 2}(undef, 7, 7)
    
    # 데이터 시작 위치 찾기
    start_idx = findfirst(line -> startswith(line, "  1l"), lines)
    if start_idx === nothing
        error("Delta matrix data not found in file")
    end
    
    # 7줄의 데이터 파싱 (1l부터 7r까지)
    for i in 1:7
        line = lines[start_idx + i - 1]
        values = split(line)[2:end]  # 첫 번째 열 (인덱스 라벨) 제외
        
        for j in 1:7
            # 각 요소는 4개의 복소수 (re/im 쌍)으로 구성됨
            re1 = parse(Float64, values[(j-1)*8 + 1])  # re 1l
            im1 = parse(Float64, values[(j-1)*8 + 2])  # im 1l
            re2 = parse(Float64, values[(j-1)*8 + 3])  # re 1r
            im2 = parse(Float64, values[(j-1)*8 + 4])  # im 1r
            re3 = parse(Float64, values[(j-1)*8 + 5])  # re 2l
            im3 = parse(Float64, values[(j-1)*8 + 6])  # im 2l
            re4 = parse(Float64, values[(j-1)*8 + 7])  # re 2r
            im4 = parse(Float64, values[(j-1)*8 + 8])  # im 2r
            
            # 복소수 계산: re1 + im1*im - (re2 + im2*im) - (re3 + im3*im) + (re4 + im4*im)
            arr[i, j] = (re1 + im1 * im) - (re2 + im2 * im) - (re3 + im3 * im) + (re4 + im4 * im)
        end
    end
    
    return arr
end

# # 사용 예시
# filename = "delta_gw.out"
# arr = parse_delta_gw_file(filename)

# # 결과 확인
# println("arr[1][1]: ", arr[1, 1])
# println("arr[7][7]: ", arr[7, 7])




# function read_dcon_input(file_path)
#     dcon_dict = Dict()
#     current_section = ""

#     for line in eachline(file_path)
#         line = strip(line)
#         if isempty(line) || startswith(line, "!")  # 빈 줄 또는 주석 제외
#             continue
#         end

#         if startswith(line, "&")  # 섹션 시작
#             current_section = strip(line[2:end])  # "&DCON_CONTROL" → "DCON_CONTROL"
#             dcon_dict[current_section] = Dict()
#         elseif line == "/"  # 섹션 끝
#             current_section = ""
#         elseif current_section != ""  # 키-값 저장
#             parts = split(line, "="; limit=2)
#             if length(parts) == 2
#                 key = strip(parts[1])
#                 value = strip(parts[2])
#                 value = replace(value, "!" => "")  # 주석 제거
#                 value = strip(split(value, " ")[1])  # 숫자나 논리값만 남김

#                 # 숫자로 변환 가능한 경우 변환
#                 if occursin(".", value) || occursin("e", value)  # 실수 변환
#                     value = tryparse(Float64, value)
#                 elseif value in ["t", "f"]  # 논리값 변환
#                     value = value == "t"
#                 else  # 정수 변환
#                     value = tryparse(Int, value) !== nothing ? parse(Int, value) : value
#                 end

#                 dcon_dict[current_section][key] = value
#             end
#         end
#     end

#     return dcon_dict
# end


# function write_dcon_input(file_path, dcon_dict)
#     open(file_path, "w") do f
#         for (section, values) in dcon_dict
#             println(f, "&$section")  # 섹션 헤더 출력

#             for (key, value) in values
#                 formatted_value = value
#                 if value isa Bool  # Boolean 값 변환
#                     formatted_value = value ? "t" : "f"
#                 elseif value isa Float64  # 실수 변환
#                     formatted_value = @sprintf("%.6E", value)  # 지수 표기법
#                 end

#                 println(f, "    $key = $formatted_value")
#             end
#             println(f, "/\n")  # 섹션 종료
#         end
#     end
#     println("✅ `$file_path` 파일이 성공적으로 업데이트되었습니다!")
# end
