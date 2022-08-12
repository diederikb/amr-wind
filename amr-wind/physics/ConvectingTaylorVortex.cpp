#include "amr-wind/physics/ConvectingTaylorVortex.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"

namespace amr_wind {
namespace ctv {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real UExactFV::operator()(
    const amrex::Real u0,
    const amrex::Real v0,
    const amrex::Real alpha,
    const amrex::Real beta,
    const amrex::Real A,
    const amrex::Real nu,
    const amrex::Real dx,
    const amrex::Real dy,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real t) const
{   
    const amrex::Real x_L = x - 0.5 * dx;
    const amrex::Real x_R = x + 0.5 * dx;
    const amrex::Real y_B = y - 0.5 * dy;
    const amrex::Real y_T = y + 0.5 * dy;
    return u0 + A / alpha * ((std::sin(alpha * (x_R - u0 * t)) - std::sin(alpha * (x_L - u0 * t))) * (std::cos(beta * (y_T - v0 * t)) - std::cos(beta * (y_B - v0 * t)))) * std::exp(-(alpha * alpha + beta * beta) * nu * t) / (dx * dy) ;

    //return u0 - A * beta * std::cos(alpha * (x - u0 * t)) *
    //                std::sin(beta * (y - v0 * t)) *
    //                std::exp(-(alpha * alpha + beta * beta) * nu * t);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real VExactFV::operator()(
    const amrex::Real u0,
    const amrex::Real v0,
    const amrex::Real alpha,
    const amrex::Real beta,
    const amrex::Real A,
    const amrex::Real nu,
    const amrex::Real dx,
    const amrex::Real dy,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real t) const
{
    const amrex::Real x_L = x - 0.5 * dx;
    const amrex::Real x_R = x + 0.5 * dx;
    const amrex::Real y_B = y - 0.5 * dy;
    const amrex::Real y_T = y + 0.5 * dy;
    return v0 - A / beta * ((std::cos(alpha * (x_R - u0 * t)) - std::cos(alpha * (x_L - u0 * t))) * (std::sin(beta * (y_T - v0 * t)) - std::sin(beta * (y_B - v0 * t)))) * std::exp(-(alpha * alpha + beta * beta) * nu * t) / (dx * dy) ;
    //return v0 + A * alpha * std::sin(alpha * (x - u0 * t)) *
    //                std::cos(beta * (y - v0 * t)) *
    //                std::exp(-(alpha * alpha + beta * beta) * nu * t);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real WExactFV::operator()(
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/) const
{
    return 0.0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real GpxExactFV::operator()(
    const amrex::Real u0,
    const amrex::Real /*unused*/,
    const amrex::Real alpha,
    const amrex::Real beta,
    const amrex::Real A,
    const amrex::Real nu,
    const amrex::Real dx,
    const amrex::Real dy,
    const amrex::Real x,
    const amrex::Real /*unused*/,
    const amrex::Real t) const
{
    return 0.5 * A * A * alpha * beta * beta *
           std::sin(2.0 * alpha * (x - u0 * t)) *
           std::exp(-2.0 * (alpha * alpha + beta * beta) * nu * t);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real GpyExactFV::operator()(
    const amrex::Real /*unused*/,
    const amrex::Real v0,
    const amrex::Real alpha,
    const amrex::Real beta,
    const amrex::Real A,
    const amrex::Real nu,
    const amrex::Real dx,
    const amrex::Real dy,
    const amrex::Real /*unused*/,
    const amrex::Real y,
    const amrex::Real t) const
{
    return 0.5 * A * A * alpha * alpha * beta *
           std::sin(2.0 * beta * (y - v0 * t)) *
           std::exp(-2.0 * (alpha * alpha + beta * beta) * nu * t);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real GpzExactFV::operator()(
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/) const
{
    return 0.0;
}

} // namespace

ConvectingTaylorVortex::ConvectingTaylorVortex(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_gradp(sim.repo().get_field("gp"))
    , m_density(sim.repo().get_field("density"))
    , m_mesh_mapping(sim.has_mesh_mapping())
{
    {
        amrex::ParmParse pp("CTV");
        pp.query("density", m_rho);
        pp.query("u0", m_u0);
        pp.query("v0", m_v0);
        pp.query("alpha", m_alpha);
        pp.query("beta", m_beta);
        pp.query("A", m_A);
        pp.query("activate_pressure", m_activate_pressure);
        pp.query("error_log_file", m_output_fname);
    }
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity", m_nu);
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        const int maxlevel = sim.mesh().maxLevel();
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time";
        for (int lev = 0; lev < maxlevel + 1; ++lev) {
            f << std::setw(m_w - 1) << "L2_u_lvl_" << lev
              << std::setw(m_w - 1) << "L2_v_lvl_" << lev
              << std::setw(m_w - 1) << "L2_w_lvl_" << lev
              << std::setw(m_w - 1) << "L2_gpx_lvl_" << lev
              << std::setw(m_w - 1) << "L2_gpy_lvl_" << lev
              << std::setw(m_w - 1) << "L2_gpz_lvl_" << lev;
        }
        f << std::setw(m_w) << "L2_u_all_lvl" << std::setw(m_w) << "L2_v_all_lvl" 
          << std::setw(m_w) << "L2_w_all_lvl" << std::setw(m_w) << "L2_gpx_all_lvl" 
          << std::setw(m_w) << "L2_gpy_all_lvl" << std::setw(m_w) << "L2_gpz_all_lvl";
        f << std::endl;
        f.close();
    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void ConvectingTaylorVortex::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    using namespace utils;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();

    const auto u0 = m_u0;
    const auto v0 = m_v0;
    const auto alpha = m_alpha;
    const auto beta = m_beta;
    const auto A = m_A;
    const auto nu = m_nu;
    const bool activate_pressure = m_activate_pressure;
    const auto mesh_mapping = m_mesh_mapping;

    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& pressure = m_repo.get_field("p")(level);
    auto& gradp = m_repo.get_field("gp")(level);
    Field const* nu_coord_cc =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    Field const* nu_coord_nd =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_nd")) : nullptr;

    density.setVal(m_rho);

    UExactFV u_exact;
    VExactFV v_exact;
    WExactFV w_exact;
    GpxExactFV gpx_exact;
    GpyExactFV gpy_exact;
    GpzExactFV gpz_exact;

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        auto vel = velocity.array(mfi);
        auto gp = gradp.array(mfi);
        amrex::Array4<amrex::Real const> nu_cc =
            mesh_mapping ? ((*nu_coord_cc)(level).const_array(mfi))
                         : amrex::Array4<amrex::Real const>();

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = mesh_mapping ? (nu_cc(i, j, k, 0))
                                             : (prob_lo[0] + (i + 0.5) * dx[0]);
                amrex::Real y = mesh_mapping ? (nu_cc(i, j, k, 1))
                                             : (prob_lo[1] + (j + 0.5) * dx[1]);

                vel(i, j, k, 0) = u_exact(u0, v0, alpha, beta, A, nu, dx[0], dx[1], x, y, 0.0);
                vel(i, j, k, 1) = v_exact(u0, v0, alpha, beta, A, nu, dx[0], dx[1], x, y, 0.0);
                vel(i, j, k, 2) = w_exact(u0, v0, alpha, beta, A, nu, dx[0], dx[1], x, y, 0.0);

                if (activate_pressure) {
                    gp(i, j, k, 0) = gpx_exact(u0, v0, alpha, beta, A, nu, dx[0], dx[1], x, y, 0.0);
                    gp(i, j, k, 1) = gpy_exact(u0, v0, alpha, beta, A, nu, dx[0], dx[1], x, y, 0.0);
                    gp(i, j, k, 2) = gpz_exact(u0, v0, alpha, beta, A, nu, dx[0], dx[1], x, y, 0.0);
                }
            });

        if (activate_pressure) {
            const auto& nbx = mfi.nodaltilebox();
            auto pres = pressure.array(mfi);
            amrex::Array4<amrex::Real const> nu_nd =
                mesh_mapping ? ((*nu_coord_nd)(level).const_array(mfi))
                             : amrex::Array4<amrex::Real const>();

            amrex::ParallelFor(
                nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real x = mesh_mapping ? (nu_nd(i, j, k, 0))
                                                 : (prob_lo[0] + i * dx[0]);
                    amrex::Real y = mesh_mapping ? (nu_nd(i, j, k, 1))
                                                 : (prob_lo[1] + j * dx[1]);

                    pres(i, j, k, 0) =
                        -0.25 * A * A * (beta * beta * std::cos(2.0 * alpha * x) +
                                 alpha * alpha * std::cos(2.0 * beta * y));
                });
        }
    }
}

template <typename T>
amrex::Real ConvectingTaylorVortex::compute_error_squared(int lev, const Field& field)
{
    amrex::Real error_squared = 0.0;
    const amrex::Real time = m_time.new_time();
    const auto u0 = m_u0;
    const auto v0 = m_v0;
    const auto alpha = m_alpha;
    const auto beta = m_beta;
    const auto A = m_A;
    const auto nu = m_nu;
    T f_exact;
    const auto comp = f_exact.m_comp;
    const auto mesh_mapping = m_mesh_mapping;

    Field const* nu_coord_cc =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    Field const* mesh_fac_cc =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;

    const int nlevels = m_repo.num_active_levels();
    if (lev <= nlevels - 1) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        if (m_sim.has_overset()) {
            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();

                const auto& iblank_arr =
                    m_repo.get_int_field("iblank_cell")(lev).array(mfi);
                const auto& imask_arr = level_mask.array(mfi);
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (iblank_arr(i, j, k) < 1) {
                            imask_arr(i, j, k) = 0;
                        }
                    });
            }
        }

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& prob_lo = m_mesh.Geom(lev).ProbLoArray();

        const auto& fld = field(lev);
        auto const& fld_arr = fld.const_arrays();
        auto const& mask_arr = level_mask.const_arrays();
        amrex::MultiArray4<amrex::Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();
        amrex::MultiArray4<amrex::Real const> nu_cc =
            mesh_mapping ? ((*nu_coord_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();

        error_squared += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real>{}, fld, amrex::IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& fld_bx = fld_arr[box_no];
                auto const& mask_bx = mask_arr[box_no];

                amrex::Real x = mesh_mapping ? (nu_cc[box_no](i, j, k, 0))
                                             : (prob_lo[0] + (i + 0.5) * dx[0]);
                amrex::Real y = mesh_mapping ? (nu_cc[box_no](i, j, k, 1))
                                             : (prob_lo[1] + (j + 0.5) * dx[1]);
                amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                const amrex::Real u = fld_bx(i, j, k, comp);
                const amrex::Real u_exact = f_exact(u0, v0, alpha, beta, A, nu, dx[0], dx[1], x, y, 0.0);
                const amrex::Real cell_vol =
                    dx[0] * fac_x * dx[1] * fac_y * dx[2] * fac_z;



                return cell_vol * mask_bx(i, j, k) * (u - u_exact) *
                       (u - u_exact);
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error_squared);

    return error_squared;
}

void ConvectingTaylorVortex::output_error()
{
    const int maxlevel = m_sim.mesh().maxLevel();

    // TODO: gradp analytical solution has not been adjusted for mesh mapping
    amrex::Real u_err_sq_lvl;
    amrex::Real v_err_sq_lvl;
    amrex::Real w_err_sq_lvl;
    amrex::Real gpx_err_sq_lvl;
    amrex::Real gpy_err_sq_lvl;
    amrex::Real gpz_err_sq_lvl;
    amrex::Real u_err_sq_all_lvl(0.0);
    amrex::Real v_err_sq_all_lvl(0.0);
    amrex::Real w_err_sq_all_lvl(0.0);
    amrex::Real gpx_err_sq_all_lvl(0.0);
    amrex::Real gpy_err_sq_all_lvl(0.0);
    amrex::Real gpz_err_sq_all_lvl(0.0);

    const amrex::Real total_vol = m_mesh.Geom(0).ProbDomain().volume();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.new_time();
        f.close();
    }
    for (int lev = 0; lev < maxlevel + 1; ++lev) {
        u_err_sq_lvl = compute_error_squared<UExactFV>(lev,m_velocity);
        v_err_sq_lvl = compute_error_squared<VExactFV>(lev,m_velocity);
        w_err_sq_lvl = compute_error_squared<WExactFV>(lev,m_velocity);
        gpx_err_sq_lvl = compute_error_squared<GpxExactFV>(lev,m_gradp);
        gpy_err_sq_lvl = compute_error_squared<GpyExactFV>(lev,m_gradp);
        gpz_err_sq_lvl = compute_error_squared<GpzExactFV>(lev,m_gradp);
        
        u_err_sq_all_lvl += u_err_sq_lvl;
        v_err_sq_all_lvl += v_err_sq_lvl;
        w_err_sq_all_lvl += w_err_sq_lvl;
        gpx_err_sq_all_lvl += gpx_err_sq_lvl;
        gpy_err_sq_all_lvl += gpy_err_sq_lvl;
        gpz_err_sq_all_lvl += gpz_err_sq_lvl;

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream f;
            f.open(m_output_fname.c_str(), std::ios_base::app);
            f << std::setw(m_w) << std::sqrt(u_err_sq_lvl / total_vol) << std::setw(m_w) << std::sqrt(v_err_sq_lvl / total_vol)
              << std::setw(m_w) << std::sqrt(w_err_sq_lvl / total_vol) << std::setw(m_w) << std::sqrt(gpx_err_sq_lvl / total_vol)
              << std::setw(m_w) << std::sqrt(gpy_err_sq_lvl / total_vol) << std::setw(m_w) << std::sqrt(gpz_err_sq_lvl / total_vol);
            f.close();
        }
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setw(m_w) << std::sqrt(u_err_sq_all_lvl / total_vol) << std::setw(m_w) << std::sqrt(v_err_sq_all_lvl / total_vol)
          << std::setw(m_w) << std::sqrt(w_err_sq_all_lvl / total_vol) << std::setw(m_w) << std::sqrt(gpx_err_sq_all_lvl / total_vol)
          << std::setw(m_w) << std::sqrt(gpy_err_sq_all_lvl / total_vol) << std::setw(m_w) << std::sqrt(gpz_err_sq_all_lvl / total_vol)
          << std::endl;
        f.close();
    }
}

void ConvectingTaylorVortex::post_init_actions() { output_error(); }

void ConvectingTaylorVortex::post_advance_work() { output_error(); }

} // namespace ctv
} // namespace amr_wind
