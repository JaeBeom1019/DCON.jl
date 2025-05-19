using CairoMakie, LaTeXStrings
using AbstractTrees, Terming

include("namelist_structure.jl")

function plot(result::Union{DCON.RdconOutputNetCDF,DCON.StrideOutputNetCDF}; log = true)

    fig = Figure(size=(600, 500))
    mat = result.Delta_prime[:, :, 1]
    codename = result isa DCON.RdconOutputNetCDF ? "Rdcon" : "Stride"
    qvals = result.q_rational

    if log == true
        mat_log = sign.(mat) .* log10.(abs.(mat) .+ 1e-20)  # log(0) 방지

        ax = Axis(fig[1, 1];
            title = latexstring("\\text{Log scaled } \\Delta' \\text{ matrix (n=1, ", codename, ")}"),
            xlabel = L"q",
            ylabel = L"q",
            yreversed = true,
            xaxisposition = :top,
            aspect = DataAspect(),
            titlesize = 25,
            xlabelsize = 20,
            ylabelsize = 20,
        )

        hm = CairoMakie.heatmap!(
            ax,
            qvals, qvals, mat_log;
            colormap = Reverse(:RdBu),
            colorrange = (-maximum(abs, mat_log), maximum(abs, mat_log))
        )


        Colorbar(fig[1, 2], hm;
            label = "sign(Δ′) × log₁₀(|Δ′|)",
            labelsize = 18,
            width = 15,
            ticklabelsize = 17,
        )
    else
        ax = Axis(fig[1, 1];
            title = latexstring("\\Delta' \\text{ matrix (n = 1, ", codename, ")}"),
            xlabel = L"m",
            ylabel = L"m",
            yreversed = true,
            xaxisposition = :top,
            aspect = DataAspect(),
            titlesize = 25,
            xlabelsize = 20,
            ylabelsize = 20,
        )

        hm = CairoMakie.heatmap!(
            ax,
            qvals, qvals, mat;
            colormap = Reverse(:RdBu),
            colorrange = (-maximum(abs, mat), maximum(abs, mat))
        )

        Colorbar(fig[1, 2], hm;
            width = 15,
            ticklabelsize = 15,
        )
    end

    display(fig)

end



# Extend `AbstractTrees.printnode` for custom node display
function AbstractTrees.printnode(io::IO, x::T) where {T}
    print(io, Terming.bold_on(), Terming.rgb(100, 150, 255), nameof(T), Terming.color_off(), Terming.bold_off())
end

# Method for fields of a struct (non-struct values)
function AbstractTrees.printnode(io::IO, x::Pair{Symbol, Any})
    field_name = x.first
    field_value = x.second
    field_type = typeof(field_value)

    print(io, Terming.bold_on(), Terming.rgb(0, 200, 200), field_name, Terming.color_off(), Terming.bold_off())
    print(io, Terming.rgb(200, 200, 200), "::", field_type, Terming.color_off())
    print(io, " = ", Terming.rgb(255, 150, 0), repr(field_value), Terming.color_off())
end

# Extend `AbstractTrees.children` for structs
AbstractTrees.children(x) = [] # Default for non-structs or leaf nodes

function AbstractTrees.children(x::T) where {T}
    if isstructtype(T)
        child_nodes = []
        for field in fieldnames(T)
            field_val = getfield(x, field)
            # If the field value is itself a struct, add it directly
            if isstructtype(typeof(field_val))
                push!(child_nodes, field_val)
            else
                # Otherwise, represent it as a Pair{Symbol, Any}
                push!(child_nodes, field => field_val)
            end
        end
        return child_nodes
    else
        return []
    end
end

# Custom `show` method for InputFiles to enable `ini` in REPL
Base.show(io::IO, mime::MIME"text/plain", x::InputFiles) = print_tree(io, x)

# Helper function to print tree with a max depth
function print_tree(io::IO, x; max_depth::Int = -1)
    if max_depth == -1
        AbstractTrees.print_tree(io, x)
    else
        AbstractTrees.print_tree(io, x, max_depth = max_depth)
    end
end