using AbstractTrees
using AbstractTrees: print_tree
using Crayons  # 컬러 출력용

# 🎨 Crayon 스타일 정의
const color_key    = Crayon(foreground=:light_red, bold=false)         # key = 필드 이름
const color_value  = Crayon(foreground=:light_gray, bold=false)                     # 값
const color_type   = Crayon(foreground=:light_blue, bold=false) # 일반 타입명 출력용

# 🌲 구조체 트리 탐색을 위한 children 정의
function AbstractTrees.children(x)
    T = typeof(x)
    if isstructtype(T)
        return [field => getfield(x, field) for field in fieldnames(T)]
    else
        return []
    end
end



# Pair 타입 (field => value) 노드 출력
function AbstractTrees.printnode(io::IO, x::Pair)
    # print directly to the IO stream without pre-stringifying ANSI codes
    print(io, color_key, string(x.first)) # Print field name in green bold
    print(io, " = ")                     # Print separator
    print(io, color_value, repr(x.second)) # Print value in blue
end

# 기본 타입의 노드 출력 (예: Vector, Nothing 등)
function AbstractTrees.printnode(io::IO, x)
    print(io, color_type, string(typeof(x))) # Print type name in light gray
end


# ===== show() 자동 오버로드 =====
# 출력 깊이를 제한한 show 함수
function show_shallow(T::Type)
    @eval Base.show(io::IO, mime::MIME"text/plain", x::$T) = print_tree(io, x; maxdepth=1, indicate_truncation=false)
end

# 트리로 보고 싶은 구조체 목록
tree_types = [
    InputFiles,
    DCONNamelist, DCONControl, DCONOutput,
    EquilNamelist, EquilControl, EquilOutput,
    VacNamelist, Modes, Debugs, Vacdat, Shape, Diagns, Sprk,
    RdconNamelist, RdconControl, RdconOutput, GalInput, GalOutput, UaDiagnoseList,
    StrideNamelist, StrideControl, StrideOutput, StrideParams
]

# 모두 적용
foreach(show_shallow, tree_types)