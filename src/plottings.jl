using Plots, LaTeXStrings

include("namelist_structure.jl")

function plot(result::Union{DCON.RdconOutputNetCDF,DCON.StrideOutputNetCDF}; log = true)
    mat = result.Delta_prime[:, :, 1]
    codename = result isa DCON.RdconOutputNetCDF ? "Rdcon" : "Stride"
    qvals = result.q_rational

    plt = if log
        mat_log = sign.(mat) .* log10.(abs.(mat) .+ 1e-20)
        heatmap(
            qvals, qvals, mat_log;
            title = latexstring("\\text{Log scaled } \\Delta' \\text{ matrix (n=1, ", codename, ")}"),
            xlabel = L"q",
            ylabel = L"q",
            colorbar_title = "sign(Δ′) × log₁₀(|Δ′|)",
            colormap = :RdBu,
            clim = (-maximum(abs, mat_log), maximum(abs, mat_log)),
            aspect_ratio = :equal,
            titlefontsize = 14,
            xlabelfontsize = 12,
            ylabelfontsize = 12,
            colorbar_titlefontsize = 10,
            tickfontsize = 10,
            yflip = true,
            xaxis = :top,
        )
    else
        heatmap(
            qvals, qvals, mat;
            title = latexstring("\\Delta' \\text{ matrix (n = 1, ", codename, ")}"),
            xlabel = L"q",
            ylabel = L"q",
            colormap = :RdBu,
            clim = (-maximum(abs, mat), maximum(abs, mat)),
            aspect_ratio = :equal,
            titlefontsize = 14,
            xlabelfontsize = 12,
            ylabelfontsize = 12,
            tickfontsize = 10,
            yflip = true,
            xaxis = :top,
        )
    end

    display(plt)   # 중요!
    return plt     # return 도 해줘야 REPL이나 스크립트에서 이어쓸 수 있음
end