using NCDatasets
using OrderedCollections
# ================================
# DCON_CONTROL 구조체 정의
# ================================


Base.@kwdef mutable struct DCONControl
    bal_flag::Bool
    mat_flag::Bool
    ode_flag::Bool
    vac_flag::Bool
    mer_flag::Bool
    sas_flag::Bool
    dmlim::Float64
    psiedge::Float64
    qlow::Float64
    qhigh::Float64
    sing_start::Int
    reform_eq_with_psilim::Bool
    nn::Int
    delta_mlow::Int
    delta_mhigh::Int
    delta_mband::Int
    mthvac::Int
    thmax0::Float64
    kin_flag::Bool
    con_flag::Bool
    kinfac1::Float64
    kinfac2::Float64
    kingridtype::Int
    passing_flag::Bool
    trapped_flag::Bool
    ktanh_flag::Bool
    ktc::Float64
    ktw::Float64
    ion_flag::Bool
    electron_flag::Bool
    dcon_kin_threads::Int
    tol_nr::Float64
    tol_r::Float64
    crossover::Float64
    singfac_min::Float64
    ucrit::Float64
    termbycross_flag::Bool
    use_classic_splines::Bool
    wv_farwall_flag::Bool
end

# ================================
# DCON_OUTPUT 구조체 정의
# ================================
Base.@kwdef mutable struct DCONOutput
    out_fund::Bool
    crit_break::Bool
    ahb_flag::Bool
    msol_ahb::Int
    mthsurf0::Int
    bin_euler::Bool
    euler_stride::Int
    out_bal1::Bool
    bin_bal1::Bool
    out_bal2::Bool
    bin_bal2::Bool
    netcdf_out::Bool
end

# ================================
# 전체 구조체
# ================================
Base.@kwdef struct DCONNamelist
    DCON_CONTROL::DCONControl
    DCON_OUTPUT::DCONOutput
end

function parse_dcon_namelist(d::OrderedDict)::DCONNamelist
    return DCONNamelist(
        DCON_CONTROL = DCONControl(; d[:DCON_CONTROL]...),
        DCON_OUTPUT = DCONOutput(; d[:DCON_OUTPUT]...)
    )
end

function struct_to_dict(control::DCONControl)
    return OrderedDict{Symbol, Any}(field => getfield(control, field) for field in fieldnames(DCONControl))
end

function struct_to_dict(output::DCONOutput)
    return OrderedDict{Symbol, Any}(field => getfield(output, field) for field in fieldnames(DCONOutput))
end

function struct_to_dict(nml::DCONNamelist)
    return OrderedDict(
        :DCON_CONTROL => struct_to_dict(nml.DCON_CONTROL),
        :DCON_OUTPUT  => struct_to_dict(nml.DCON_OUTPUT)
    )
end

#EquilControl

Base.@kwdef mutable struct EquilControl
    eq_type::String
    eq_filename::String
    jac_type::String
    power_bp::Int
    power_b::Int
    power_r::Int
    grid_type::String
    psilow::Float64
    psihigh::Float64
    mpsi::Int
    mtheta::Int
    newq0::Int
    etol::Float64
    use_classic_splines::Bool
    input_only::Bool
    use_galgrid::Bool
end

Base.@kwdef mutable struct EquilOutput
    gse_flag::Bool
    out_eq_1d::Bool
    bin_eq_1d::Bool
    out_eq_2d::Bool
    bin_eq_2d::Bool
    out_2d::Bool
    bin_2d::Bool
    dump_flag::Bool
end

Base.@kwdef struct EquilNamelist
    EQUIL_CONTROL::EquilControl
    EQUIL_OUTPUT::EquilOutput
end

function parse_equil_namelist(d::OrderedDict)::EquilNamelist
    return EquilNamelist(
        EQUIL_CONTROL = EquilControl(; d[:EQUIL_CONTROL]...),
        EQUIL_OUTPUT  = EquilOutput(; d[:EQUIL_OUTPUT]...)
    )
end

function struct_to_dict(c::EquilControl)
    return OrderedDict{Symbol, Any}(f => getfield(c, f) for f in fieldnames(EquilControl))
end

function struct_to_dict(c::EquilOutput)
    return OrderedDict{Symbol, Any}(f => getfield(c, f) for f in fieldnames(EquilOutput))
end

function struct_to_dict(nml::EquilNamelist)
    return OrderedDict(
        :EQUIL_CONTROL => struct_to_dict(nml.EQUIL_CONTROL),
        :EQUIL_OUTPUT  => struct_to_dict(nml.EQUIL_OUTPUT)
    )
end

# ========== VAC.IN STRUCTURES ========== #
Base.@kwdef mutable struct Modes
    mth::Int
    lsymz::Bool
    leqarcw::Int
    lzio::Int
    lgato::Int
    lrgato::Int
end

Base.@kwdef mutable struct Debugs
    checkd::Bool
    check1::Bool
    check2::Bool
    checke::Bool
    checks::Bool
    wall::Bool
    lkplt::Int
    verbose_timer_output::Bool
end

Base.@kwdef mutable struct Vacdat
    ishape::Int
    aw::Float64
    bw::Float64
    cw::Float64
    dw::Float64
    tw::Float64
    nsing::Int
    epsq::Float64
    noutv::Int
    idgt::Int
    idot::Int
    idsk::Int
    delg::Float64
    delfac::Float64
    cn0::Int
end

Base.@kwdef mutable struct Shape
    ipshp::Int
    xpl::Float64
    apl::Float64
    a::Float64
    b::Float64
    bpl::Float64
    dpl::Float64
    r::Float64
    abulg::Float64
    bbulg::Float64
    tbulg::Float64
    qain::Float64
end

Base.@kwdef mutable struct Diagns
    lkdis::Bool
    ieig::Int
    iloop::Int
    lpsub::Int
    nloop::Int
    nloopr::Int
    nphil::Int
    nphse::Int
    xofsl::Float64
    ntloop::Int
    aloop::Float64
    bloop::Float64
    dloop::Float64
    rloop::Float64
    deloop::Float64
    mx::Int
    mz::Int
    nph::Int
    nxlpin::Int
    nzlpin::Int
    epslp::Float64
    xlpmin::Float64
    xlpmax::Float64
    zlpmin::Float64
    zlpmax::Float64
    linterior::Int
end

Base.@kwdef mutable struct Sprk
    nminus::Int
    nplus::Int
    mphi::Int
    lwrt11::Int
    civ::Float64
    sp2sgn1::Int
    sp2sgn2::Int
    sp2sgn3::Int
    sp2sgn4::Int
    sp2sgn5::Int
    sp3sgn1::Int
    sp3sgn2::Int
    sp3sgn3::Int
    sp3sgn4::Int
    sp3sgn5::Int
    lff::Int
    ff::Float64
    fv::Vector{Float64}
end

Base.@kwdef struct VacNamelist
    modes::Modes
    debugs::Debugs
    vacdat::Vacdat
    shape::Shape
    diagns::Diagns
    sprk::Sprk
end

function parse_vac_namelist(d::OrderedDict)::VacNamelist
    return VacNamelist(
        modes = Modes(; d[:modes]...),
        debugs = Debugs(; d[:debugs]...),
        vacdat = Vacdat(; d[:vacdat]...),
        shape = Shape(; d[:shape]...),
        diagns = Diagns(; d[:diagns]...),
        sprk = Sprk(; d[:sprk]...)
    )
end

function struct_to_dict(v::VacNamelist)
    return OrderedDict(
        :modes => OrderedDict{Symbol, Any}(field => getfield(v.modes, field) for field in fieldnames(Modes)),
        :debugs => OrderedDict{Symbol, Any}(field => getfield(v.debugs, field) for field in fieldnames(Debugs)),
        :vacdat => OrderedDict{Symbol, Any}(field => getfield(v.vacdat, field) for field in fieldnames(Vacdat)),
        :shape => OrderedDict{Symbol, Any}(field => getfield(v.shape, field) for field in fieldnames(Shape)),
        :diagns => OrderedDict{Symbol, Any}(field => getfield(v.diagns, field) for field in fieldnames(Diagns)),
        :sprk => OrderedDict{Symbol, Any}(field => getfield(v.sprk, field) for field in fieldnames(Sprk))
    )
end

# ==================================
# RDCON_CONTROL 섹션 구조체 정의
# ==================================
Base.@kwdef mutable struct RdconControl
    bal_flag::Bool
    mat_flag::Bool
    ode_flag::Bool
    vac_flag::Bool
    gal_flag::Bool
    sas_flag::Bool
    dmlim::Float64
    sing_start::Int
    nn::Int
    delta_mlow::Int
    delta_mhigh::Int
    delta_mband::Int
    mthvac::Int
    thmax0::Float64
    tol_nr::Float64
    tol_r::Float64
    crossover::Float64
    singfac_min::Float64
    ucrit::Float64
    cyl_flag::Bool
    sing1_flag::Bool
    sing_order::Int
    sing_order_ceiling::Bool
    regrid_flag::Bool
end

# ==================================
# RDCON_OUTPUT 섹션 구조체 정의
# ==================================
Base.@kwdef mutable struct RdconOutput
    crit_break::Bool
    ahb_flag::Bool
    msol_ahb::Int
    mthsurf0::Int
    bin_euler::Bool
    euler_stride::Int
    out_bal1::Bool
    bin_bal1::Bool
    out_bal2::Bool
    bin_bal2::Bool
end

# ==================================
# GAL_INPUT 구조체
# ==================================
Base.@kwdef mutable struct GalInput
    nx::Int
    pfac::Float64
    gal_tol::Float64
    dx1dx2_flag::Bool
    dx0::Float64
    dx1::Float64
    dx2::Float64
    cutoff::Int
    solver::String
    nq::Int
    diagnose_integrand::Bool
    diagnose_map::Bool
    diagnose_grid::Bool
    diagnose_lsode::Bool
    ndiagnose::Int
end

# ==================================
# GAL_OUTPUT 구조체
# ==================================
Base.@kwdef mutable struct GalOutput
    interp_np::Int
    restore_uh::Bool
    restore_us::Bool
    restore_ul::Bool
    bin_delmatch::Bool
    out_galsol::Bool
    bin_galsol::Bool
    b_flag::Bool
end

# ==================================
# UA_DIAGNOSE_LIST 구조체 (비어있을 수 있음)
# ==================================
Base.@kwdef mutable struct UaDiagnoseList
end

# ==================================
# 전체 RDCON NAMELIST 구조체
# ==================================
Base.@kwdef struct RdconNamelist
    GAL_INPUT::GalInput
    GAL_OUTPUT::GalOutput
    RDCON_CONTROL::RdconControl
    RDCON_OUTPUT::RdconOutput
    UA_DIAGNOSE_LIST::UaDiagnoseList
end

function parse_rdcon_namelist(d::OrderedDict)::RdconNamelist
    return RdconNamelist(
        GAL_INPUT = GalInput(; d[:GAL_INPUT]...),
        GAL_OUTPUT = GalOutput(; d[:GAL_OUTPUT]...),
        RDCON_CONTROL = RdconControl(; d[:RDCON_CONTROL]...),
        RDCON_OUTPUT = RdconOutput(; d[:RDCON_OUTPUT]...),
        UA_DIAGNOSE_LIST = UaDiagnoseList()  # 비어있음
    )
end

function struct_to_dict(x::GalInput)
    return OrderedDict(Symbol(f) => getfield(x, f) for f in fieldnames(GalInput))
end

function struct_to_dict(x::GalOutput)
    return OrderedDict(Symbol(f) => getfield(x, f) for f in fieldnames(GalOutput))
end

function struct_to_dict(x::RdconControl)
    return OrderedDict(Symbol(f) => getfield(x, f) for f in fieldnames(RdconControl))
end

function struct_to_dict(x::RdconOutput)
    return OrderedDict(Symbol(f) => getfield(x, f) for f in fieldnames(RdconOutput))
end

function struct_to_dict(x::UaDiagnoseList)
    return OrderedDict{Symbol, Any}()
end

function struct_to_dict(nml::RdconNamelist)
    return OrderedDict(
        :GAL_INPUT => struct_to_dict(nml.GAL_INPUT),
        :GAL_OUTPUT => struct_to_dict(nml.GAL_OUTPUT),
        :RDCON_CONTROL => struct_to_dict(nml.RDCON_CONTROL),
        :RDCON_OUTPUT => struct_to_dict(nml.RDCON_OUTPUT),
        :UA_DIAGNOSE_LIST => struct_to_dict(nml.UA_DIAGNOSE_LIST),
    )
end


Base.@kwdef mutable struct StrideControl
    bal_flag::Bool
    mat_flag::Bool
    ode_flag::Bool
    vac_flag::Bool
    mer_flag::Bool
    sas_flag::Bool
    dmlim::Float64
    qlow::Float64
    qhigh::Float64
    sing_start::Int
    nn::Int
    delta_mlow::Int
    delta_mhigh::Int
    delta_mband::Int
    mthvac::Int
    thmax0::Float64
    tol_nr::Float64
    tol_r::Float64
    crossover::Float64
    singfac_min::Float64
    ucrit::Float64
    sing_order::Int
    use_classic_splines::Bool
    use_notaknot_splines::Bool
end

Base.@kwdef mutable struct StrideOutput
    crit_break::Bool
    ahb_flag::Bool
    msol_ahb::Int
    mthsurf0::Int
    bin_euler::Bool
    euler_stride::Int
    out_bal1::Bool
    bin_bal1::Bool
    out_bal2::Bool
    bin_bal2::Bool
    netcdf_out::Bool
end

Base.@kwdef mutable struct StrideParams
    nThreads::Int
    fourfit_metric_parallel::Bool
    vac_parallel::Bool
    nIntervalsTot::Int
    grid_packing::String
    axis_mid_pt_skew::Float64
    asymp_at_sing::Bool
    kill_big_soln_for_ideal_dW::Bool
    calc_delta_prime::Bool
    calc_dp_with_vac::Bool
    big_soln_err_tol::Float64
    integrate_riccati::Bool
    riccati_bounce::Bool
    riccati_match_hamiltonian_evals::Bool
    verbose_riccati_output::Bool
    ric_dt::Float64
    ric_tol::Float64
    verbose_performance_output::Bool
end

Base.@kwdef struct StrideNamelist
    stride_control::StrideControl
    stride_output::StrideOutput
    stride_params::StrideParams
end

function parse_stride_namelist(d::OrderedDict)::StrideNamelist
    return StrideNamelist(
        stride_control = StrideControl(; d[:stride_control]...),
        stride_output  = StrideOutput(; d[:stride_output]...),
        stride_params  = StrideParams(; d[:stride_params]...)
    )
end

function struct_to_dict(s::StrideControl)
    return OrderedDict{Symbol, Any}(field => getfield(s, field) for field in fieldnames(StrideControl))
end
function struct_to_dict(s::StrideOutput)
    return OrderedDict{Symbol, Any}(field => getfield(s, field) for field in fieldnames(StrideOutput))
end
function struct_to_dict(s::StrideParams)
    return OrderedDict{Symbol, Any}(field => getfield(s, field) for field in fieldnames(StrideParams))
end

function struct_to_dict(nml::StrideNamelist)
    return OrderedDict(
        :stride_control => struct_to_dict(nml.stride_control),
        :stride_output  => struct_to_dict(nml.stride_output),
        :stride_params  => struct_to_dict(nml.stride_params)
    )
end


#InputFiles
Base.@kwdef mutable struct InputFiles
    dcon::DCONNamelist
    rdcon::RdconNamelist
    stride::StrideNamelist
    vac::VacNamelist
    equil::EquilNamelist
    ideal_n::Union{Nothing, Vector{Int}} = nothing
    resistive_n::Union{Nothing, Vector{Int}} = nothing
    ballooning::Bool
    directory::Union{Nothing,String} = nothing
    remove::Bool
end


Base.@kwdef struct DCONOutputNetCDF
    # Dimensions
    i::Vector{Int}
    m::Vector{Int}
    mode::Vector{Int}
    psi_n::Vector{Float64}

    # 1D Profiles
    f::Vector{Float64}
    mu0p::Vector{Float64}
    dvdpsi::Vector{Float64}
    q::Vector{Float64}
    di::Vector{Float64}
    dr::Vector{Float64}
    ca1::Vector{Float64}

    # Eigenvectors and Eigenvalues
    W_p_eigenvector::Array{Float64, 3}
    W_p_eigenvalue::Array{Float64, 2}

    W_v_eigenvector::Array{Float64, 3}
    W_v_eigenvalue::Array{Float64, 2}

    W_t_eigenvector::Array{Float64, 3}
    W_t_eigenvalue::Array{Float64, 2}

    W_t::Array{Float64, 3}

    # Global attributes
    title::String
    jacobian::String
    power_bp::Int
    power_b::Int
    power_r::Int
    mpsi::Int
    mtheta::Int
    mlow::Int
    mhigh::Int
    mpert::Int
    mband::Int
    psilow::Float64
    amean::Float64
    rmean::Float64
    aratio::Float64
    kappa::Float64
    delta1::Float64
    delta2::Float64
    li1::Float64
    li2::Float64
    li3::Float64
    ro::Float64
    zo::Float64
    psio::Float64
    betap1::Float64
    betap2::Float64
    betap3::Float64
    betat::Float64
    betan::Float64
    bt0::Float64
    q0::Float64
    qmin::Float64
    qmax::Float64
    qa::Float64
    crnt::Float64
    q95::Float64
    qlim::Float64
    psilim::Float64
    shot::Int
    time::Int
    n::Int
    version::String
    cpu_time::Float64
    wall_time::Float64
end


Base.@kwdef struct StrideOutputNetCDF
    i::Vector{Int32}
    m::Vector{Int32}
    mode::Vector{Int32}
    psi_n::Vector{Float64}
    lr_index::Vector{Int32}
    lr_prime::Vector{Int32}
    r::Vector{Int32}
    r_prime::Vector{Int32}
    psi_n_rational::Vector{Float64}
    q_rational::Vector{Float64}
    f::Vector{Float64}
    mu0p::Vector{Float64}
    dvdpsi::Vector{Float64}
    q::Vector{Float64}
    di::Vector{Float64}
    dr::Vector{Float64}
    ca1::Vector{Float64}
    W_p_eigenvector::Array{Float64,3}
    W_p_eigenvalue::Matrix{Float64}
    W_v_eigenvector::Array{Float64,3}
    W_v_eigenvalue::Matrix{Float64}
    W_t_eigenvector::Array{Float64,3}
    W_t_eigenvalue::Matrix{Float64}
    Delta::Array{Float64,3}
    A_prime::Array{Float64,3}
    B_prime::Array{Float64,3}
    Gamma_prime::Array{Float64,3}
    Delta_prime::Array{Float64,3}

    # global attributes
    title::String
    jacobian::String
    power_bp::Int
    power_b::Int
    power_r::Int
    mpsi::Int
    mtheta::Int
    mlow::Int
    mhigh::Int
    mpert::Int
    mband::Int
    psilow::Float64
    amean::Float64
    rmean::Float64
    aratio::Float64
    kappa::Float64
    delta1::Float64
    delta2::Float64
    li1::Float64
    li2::Float64
    li3::Float64
    ro::Float64
    zo::Float64
    psio::Float64
    betap1::Float64
    betap2::Float64
    betap3::Float64
    betat::Float64
    betan::Float64
    bt0::Float64
    q0::Float64
    qmin::Float64
    qmax::Float64
    qa::Float64
    crnt::Float64
    q95::Float64
    qlim::Float64
    psilim::Float64
    shot::Int
    time::Int
    n::Int
    version::String
end



function load_dcon_output_nc(path::String)::DCONOutputNetCDF
    ds = Dataset(path, "r")
    output = DCONOutputNetCDF(
        i = ds["i"][:],
        m = ds["m"][:],
        mode = ds["mode"][:],
        psi_n = ds["psi_n"][:],
        f = ds["f"][:],
        mu0p = ds["mu0p"][:],
        dvdpsi = ds["dvdpsi"][:],
        q = ds["q"][:],
        di = ds["di"][:],
        dr = ds["dr"][:],
        ca1 = ds["ca1"][:],
        W_p_eigenvector = ds["W_p_eigenvector"][:, :, :],
        W_p_eigenvalue = ds["W_p_eigenvalue"][:, :],
        W_v_eigenvector = ds["W_v_eigenvector"][:, :, :],
        W_v_eigenvalue = ds["W_v_eigenvalue"][:, :],
        W_t_eigenvector = ds["W_t_eigenvector"][:, :, :],
        W_t_eigenvalue = ds["W_t_eigenvalue"][:, :],
        W_t = ds["W_t"][:, :, :],

        # global attrs
        title = ds.attrib["title"],
        jacobian = ds.attrib["jacobian"],
        power_bp = ds.attrib["power_bp"],
        power_b = ds.attrib["power_b"],
        power_r = ds.attrib["power_r"],
        mpsi = ds.attrib["mpsi"],
        mtheta = ds.attrib["mtheta"],
        mlow = ds.attrib["mlow"],
        mhigh = ds.attrib["mhigh"],
        mpert = ds.attrib["mpert"],
        mband = ds.attrib["mband"],
        psilow = ds.attrib["psilow"],
        amean = ds.attrib["amean"],
        rmean = ds.attrib["rmean"],
        aratio = ds.attrib["aratio"],
        kappa = ds.attrib["kappa"],
        delta1 = ds.attrib["delta1"],
        delta2 = ds.attrib["delta2"],
        li1 = ds.attrib["li1"],
        li2 = ds.attrib["li2"],
        li3 = ds.attrib["li3"],
        ro = ds.attrib["ro"],
        zo = ds.attrib["zo"],
        psio = ds.attrib["psio"],
        betap1 = ds.attrib["betap1"],
        betap2 = ds.attrib["betap2"],
        betap3 = ds.attrib["betap3"],
        betat = ds.attrib["betat"],
        betan = ds.attrib["betan"],
        bt0 = ds.attrib["bt0"],
        q0 = ds.attrib["q0"],
        qmin = ds.attrib["qmin"],
        qmax = ds.attrib["qmax"],
        qa = ds.attrib["qa"],
        crnt = ds.attrib["crnt"],
        q95 = ds.attrib["q95"],
        qlim = ds.attrib["qlim"],
        psilim = ds.attrib["psilim"],
        shot = ds.attrib["shot"],
        time = ds.attrib["time"],
        n = ds.attrib["n"],
        version = ds.attrib["version"],
        cpu_time = ds.attrib["cpu_time"],
        wall_time = ds.attrib["wall_time"],
    )
    close(ds)
    return output
end

function load_stride_output_nc(path::String)::StrideOutputNetCDF
    ds = Dataset(path, "r")

    result = StrideOutputNetCDF(
        i = ds["i"][:],
        m = ds["m"][:],
        mode = ds["mode"][:],
        psi_n = ds["psi_n"][:],
        lr_index = ds["lr_index"][:],
        lr_prime = ds["lr_prime"][:],
        r = ds["r"][:],
        r_prime = ds["r_prime"][:],
        psi_n_rational = ds["psi_n_rational"][:],
        q_rational = ds["q_rational"][:],
        f = ds["f"][:],
        mu0p = ds["mu0p"][:],
        dvdpsi = ds["dvdpsi"][:],
        q = ds["q"][:],
        di = ds["di"][:],
        dr = ds["dr"][:],
        ca1 = ds["ca1"][:],
        W_p_eigenvector = ds["W_p_eigenvector"][:, :, :],
        W_p_eigenvalue = ds["W_p_eigenvalue"][:, :],
        W_v_eigenvector = ds["W_v_eigenvector"][:, :, :],
        W_v_eigenvalue = ds["W_v_eigenvalue"][:, :],
        W_t_eigenvector = ds["W_t_eigenvector"][:, :, :],
        W_t_eigenvalue = ds["W_t_eigenvalue"][:, :],
        Delta = ds["Delta"][:, :, :],
        A_prime = ds["A_prime"][:, :, :],
        B_prime = ds["B_prime"][:, :, :],
        Gamma_prime = ds["Gamma_prime"][:, :, :],
        Delta_prime = ds["Delta_prime"][:, :, :],

        title = ds.attrib["title"],
        jacobian = ds.attrib["jacobian"],
        power_bp = ds.attrib["power_bp"],
        power_b = ds.attrib["power_b"],
        power_r = ds.attrib["power_r"],
        mpsi = ds.attrib["mpsi"],
        mtheta = ds.attrib["mtheta"],
        mlow = ds.attrib["mlow"],
        mhigh = ds.attrib["mhigh"],
        mpert = ds.attrib["mpert"],
        mband = ds.attrib["mband"],
        psilow = ds.attrib["psilow"],
        amean = ds.attrib["amean"],
        rmean = ds.attrib["rmean"],
        aratio = ds.attrib["aratio"],
        kappa = ds.attrib["kappa"],
        delta1 = ds.attrib["delta1"],
        delta2 = ds.attrib["delta2"],
        li1 = ds.attrib["li1"],
        li2 = ds.attrib["li2"],
        li3 = ds.attrib["li3"],
        ro = ds.attrib["ro"],
        zo = ds.attrib["zo"],
        psio = ds.attrib["psio"],
        betap1 = ds.attrib["betap1"],
        betap2 = ds.attrib["betap2"],
        betap3 = ds.attrib["betap3"],
        betat = ds.attrib["betat"],
        betan = ds.attrib["betan"],
        bt0 = ds.attrib["bt0"],
        q0 = ds.attrib["q0"],
        qmin = ds.attrib["qmin"],
        qmax = ds.attrib["qmax"],
        qa = ds.attrib["qa"],
        crnt = ds.attrib["crnt"],
        q95 = ds.attrib["q95"],
        qlim = ds.attrib["qlim"],
        psilim = ds.attrib["psilim"],
        shot = ds.attrib["shot"],
        time = ds.attrib["time"],
        n = ds.attrib["n"],
        version = ds.attrib["version"]
    )

    close(ds)
    return result
end

Base.@kwdef struct RdconOutputNetCDF
    # dimensions
    i::Vector{Int32}
    m::Vector{Int32}
    mode::Vector{Int32}
    psi_n::Vector{Float64}
    lr_index::Vector{Int32}
    lr_prime::Vector{Int32}
    r::Vector{Int32}
    r_prime::Vector{Int32}

    # 1D profiles
    psi_n_rational::Vector{Float64}
    q_rational::Vector{Float64}
    f::Vector{Float64}
    mu0p::Vector{Float64}
    dvdpsi::Vector{Float64}
    q::Vector{Float64}
    di::Vector{Float64}
    dr::Vector{Float64}
    ca1::Vector{Float64}

    # 2D / 3D matrices
    W_p_eigenvector::Array{Float64, 3}
    W_p_eigenvalue::Array{Float64, 2}
    W_v_eigenvector::Array{Float64, 3}
    W_v_eigenvalue::Array{Float64, 2}
    W_t_eigenvector::Array{Float64, 3}
    W_t_eigenvalue::Array{Float64, 2}
    Delta::Array{Float64, 3}
    A_prime::Array{Float64, 3}
    B_prime::Array{Float64, 3}
    Gamma_prime::Array{Float64, 3}
    Delta_prime::Array{Float64, 3}
    A_prime_sym::Array{Float64, 3}
    B_prime_sym::Array{Float64, 3}
    Gamma_prime_sym::Array{Float64, 3}
    Delta_prime_sym::Array{Float64, 3}

    # global attributes
    title::String
    jacobian::String
    power_bp::Int
    power_b::Int
    power_r::Int
    mpsi::Int
    mtheta::Int
    mlow::Int
    mhigh::Int
    mpert::Int
    mband::Int
    psilow::Float64
    amean::Float64
    rmean::Float64
    aratio::Float64
    kappa::Float64
    delta1::Float64
    delta2::Float64
    li1::Float64
    li2::Float64
    li3::Float64
    ro::Float64
    zo::Float64
    psio::Float64
    betap1::Float64
    betap2::Float64
    betap3::Float64
    betat::Float64
    betan::Float64
    bt0::Float64
    q0::Float64
    qmin::Float64
    qmax::Float64
    qa::Float64
    crnt::Float64
    q95::Float64
    qlim::Float64
    psilim::Float64
    shot::Int
    time::Int
    n::Int
    version::String
end

function load_rdcon_output_nc(path::String)::RdconOutputNetCDF
    ds = Dataset(path, "r")
    output = RdconOutputNetCDF(
        i = ds["i"][:],
        m = ds["m"][:],
        mode = ds["mode"][:],
        psi_n = ds["psi_n"][:],
        lr_index = ds["lr_index"][:],
        lr_prime = ds["lr_prime"][:],
        r = ds["r"][:],
        r_prime = ds["r_prime"][:],
        psi_n_rational = ds["psi_n_rational"][:],
        q_rational = ds["q_rational"][:],
        f = ds["f"][:],
        mu0p = ds["mu0p"][:],
        dvdpsi = ds["dvdpsi"][:],
        q = ds["q"][:],
        di = ds["di"][:],
        dr = ds["dr"][:],
        ca1 = ds["ca1"][:],
        W_p_eigenvector = ds["W_p_eigenvector"][:, :, :],
        W_p_eigenvalue = ds["W_p_eigenvalue"][:, :],
        W_v_eigenvector = ds["W_v_eigenvector"][:, :, :],
        W_v_eigenvalue = ds["W_v_eigenvalue"][:, :],
        W_t_eigenvector = ds["W_t_eigenvector"][:, :, :],
        W_t_eigenvalue = ds["W_t_eigenvalue"][:, :],
        Delta = ds["Delta"][:, :, :],
        A_prime = ds["A_prime"][:, :, :],
        B_prime = ds["B_prime"][:, :, :],
        Gamma_prime = ds["Gamma_prime"][:, :, :],
        Delta_prime = ds["Delta_prime"][:, :, :],
        A_prime_sym = ds["A_prime_sym"][:, :, :],
        B_prime_sym = ds["B_prime_sym"][:, :, :],
        Gamma_prime_sym = ds["Gamma_prime_sym"][:, :, :],
        Delta_prime_sym = ds["Delta_prime_sym"][:, :, :],
        title = ds.attrib["title"],
        jacobian = ds.attrib["jacobian"],
        power_bp = ds.attrib["power_bp"],
        power_b = ds.attrib["power_b"],
        power_r = ds.attrib["power_r"],
        mpsi = ds.attrib["mpsi"],
        mtheta = ds.attrib["mtheta"],
        mlow = ds.attrib["mlow"],
        mhigh = ds.attrib["mhigh"],
        mpert = ds.attrib["mpert"],
        mband = ds.attrib["mband"],
        psilow = ds.attrib["psilow"],
        amean = ds.attrib["amean"],
        rmean = ds.attrib["rmean"],
        aratio = ds.attrib["aratio"],
        kappa = ds.attrib["kappa"],
        delta1 = ds.attrib["delta1"],
        delta2 = ds.attrib["delta2"],
        li1 = ds.attrib["li1"],
        li2 = ds.attrib["li2"],
        li3 = ds.attrib["li3"],
        ro = ds.attrib["ro"],
        zo = ds.attrib["zo"],
        psio = ds.attrib["psio"],
        betap1 = ds.attrib["betap1"],
        betap2 = ds.attrib["betap2"],
        betap3 = ds.attrib["betap3"],
        betat = ds.attrib["betat"],
        betan = ds.attrib["betan"],
        bt0 = ds.attrib["bt0"],
        q0 = ds.attrib["q0"],
        qmin = ds.attrib["qmin"],
        qmax = ds.attrib["qmax"],
        qa = ds.attrib["qa"],
        crnt = ds.attrib["crnt"],
        q95 = ds.attrib["q95"],
        qlim = ds.attrib["qlim"],
        psilim = ds.attrib["psilim"],
        shot = ds.attrib["shot"],
        time = ds.attrib["time"],
        n = ds.attrib["n"],
        version = get(ds.attrib, "version", "")
    )
    close(ds)
    return output
end