using AbstractTrees
using AbstractTrees: print_tree
using Crayons  # 컬러 출력용

# 🎨 Crayon 스타일 정의
const color_key    = Crayon(foreground=:light_red, bold=false)       
const color_value  = Crayon(foreground=:light_gray, bold=false) 
const color_type   = Crayon(foreground=:light_blue, bold=false) 

# Define children
function AbstractTrees.children(x)
    T = typeof(x)
    if isstructtype(T)
        return [field => getfield(x, field) for field in fieldnames(T)]
    else
        return []
    end
end



# Pair type (field => value)
function AbstractTrees.printnode(io::IO, x::Pair)
    # print directly to the IO stream without pre-stringifying ANSI codes
    print(io, color_key, string(x.first)) # Print field name in green bold
    print(io, " = ")                     # Print separator
    print(io, color_value, repr(x.second)) # Print value in blue
end

# basic type node print (ex: Vector, Nothing)
function AbstractTrees.printnode(io::IO, x)
    print(io, color_type, string(typeof(x))) # Print type name in light gray
end


# ===== show() overload =====
function show_shallow(T::Type)
    @eval Base.show(io::IO, mime::MIME"text/plain", x::$T) = print_tree(io, x; maxdepth=1, indicate_truncation=false)
end

# Structure list
tree_types = [
    InputFiles,
    DCONNamelist, DCONControl, DCONOutput,
    EquilNamelist, EquilControl, EquilOutput,
    VacNamelist, Modes, Debugs, Vacdat, Shape, Diagns, Sprk,
    RdconNamelist, RdconControl, RdconOutput, GalInput, GalOutput, UaDiagnoseList,
    StrideNamelist, StrideControl, StrideOutput, StrideParams
]

foreach(show_shallow, tree_types)