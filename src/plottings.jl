using CairoMakie, LaTeXStrings


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
