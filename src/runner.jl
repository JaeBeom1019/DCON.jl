function run_code(exe_path::String)
    
    env_limited = merge(ENV, Dict(
        "OMP_NUM_THREADS" => "1",
        "OPENBLAS_NUM_THREADS" => "1",
        "MKL_NUM_THREADS" => "1",
        "NUMEXPR_NUM_THREADS" => "1"
    ))

    run(Cmd(`$exe_path`; env = env_limited))
end
