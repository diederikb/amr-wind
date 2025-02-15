#include "amr-wind/utilities/sampling/FreeSurface.H"
#include "amr-wind/utilities/io_utils.H"
#include <AMReX_MultiFabUtil.H>
#include <utility>
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace free_surface {

FreeSurface::FreeSurface(CFDSim& sim, std::string label)
    : m_sim(sim), m_label(std::move(label)), m_vof(sim.repo().get_field("vof"))
{
#ifdef AMREX_USE_GPU
    amrex::Print() << "WARNING: FreeSurface: Running on GPUs..." << std::endl;
#endif
}

FreeSurface::~FreeSurface() = default;

void FreeSurface::initialize()
{
    BL_PROFILE("amr-wind::FreeSurface::initialize");

    {
        amrex::ParmParse pp(m_label);
        pp.query("output_frequency", m_out_freq);
        pp.query("output_format", m_out_fmt);
        // Load parameters of freesurface sampling
        pp.query("num_instances", m_ninst);
        pp.query("search_direction", m_coorddir);
        pp.getarr("num_points", m_npts_dir);
        pp.getarr("start", m_start);
        pp.getarr("end", m_end);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_start.size()) == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_end.size()) == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == 2);

        switch (m_coorddir) {
        case 0: {
            m_gc1 = 1;
            m_gc2 = 2;
            break;
        }
        case 1: {
            m_gc1 = 0;
            m_gc2 = 2;
            break;
        }
        case 2: {
            m_gc1 = 0;
            m_gc2 = 1;
            break;
        }
        default: {
            amrex::Abort(
                "FreeSurface: Invalid coordinate search direction "
                "encountered");
            break;
        }
        }
    }

    // Calculate total number of points
    m_npts = m_npts_dir[0] * m_npts_dir[1];

    // Turn parameters into 2D grid
    m_locs.resize(m_npts);
    m_out.resize(m_npts * m_ninst);

    // Get size of sample grid spacing
    amrex::Vector<amrex::Real> dx = {0.0, 0.0};
    dx[0] = (m_end[m_gc1] - m_start[m_gc1]) / amrex::max(m_npts_dir[0] - 1, 1);
    dx[1] = (m_end[m_gc2] - m_start[m_gc2]) / amrex::max(m_npts_dir[1] - 1, 1);

    // Store locations
    int idx = 0;
    for (int j = 0; j < m_npts_dir[1]; ++j) {
        for (int i = 0; i < m_npts_dir[0]; ++i) {
            // Initialize output values to 0.0
            for (int ni = 0; ni < m_ninst; ++ni) {
                m_out[idx * m_ninst + ni] = m_start[m_coorddir];
            }
            // Grid direction 1
            m_locs[idx][0] = m_start[m_gc1] + dx[0] * i;
            // Grid direction 2
            m_locs[idx][1] = m_start[m_gc2] + dx[1] * j;

            ++idx;
        }
    }

    if (m_out_fmt == "netcdf") {
        prepare_netcdf_file();
    }
}

void FreeSurface::post_advance_work()
{

    BL_PROFILE("amr-wind::FreeSurface::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) {
        return;
    }

    // Zero data in output array
    for (int n = 0; n < m_npts * m_ninst; n++) {
        m_out[n] = 0.0;
    }
    // Set up device vector of outputs, initialize to above phi0
    const auto& plo0 = m_sim.mesh().Geom(0).ProbLoArray();
    const auto& phi0 = m_sim.mesh().Geom(0).ProbHiArray();
    amrex::Gpu::DeviceVector<amrex::Real> dout(m_npts, phi0[2] + 1.0);
    auto* dout_ptr = dout.data();

    // Sum of interface locations at each point (assumes one interface only)
    const int finest_level = m_vof.repo().num_active_levels() - 1;

    // Loop instances
    for (int ni = 0; ni < m_ninst; ++ni) {
        for (int lev = 0; lev <= finest_level; lev++) {

            // Use level_mask to identify smallest volume
            amrex::iMultiFab level_mask;
            if (lev < finest_level) {
                level_mask = makeFineMask(
                    m_sim.mesh().boxArray(lev),
                    m_sim.mesh().DistributionMap(lev),
                    m_sim.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
            } else {
                level_mask.define(
                    m_sim.mesh().boxArray(lev),
                    m_sim.mesh().DistributionMap(lev), 1, 0, amrex::MFInfo());
                level_mask.setVal(1);
            }

            const auto& vof = m_vof(lev);
            const auto& geom = m_sim.mesh().Geom(lev);
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
                geom.CellSizeArray();
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxi =
                geom.InvCellSizeArray();
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo =
                geom.ProbLoArray();

            // c = captured
            const int dir = m_coorddir;
            const int gc1 = m_gc1;
            const int gc2 = m_gc2;

            // Loop points in 2D grid
            for (int n = 0; n < m_npts; ++n) {
                amrex::GpuArray<amrex::Real, 2> loc;
                for (int d = 0; d < 2; ++d) {
                    loc[d] = m_locs[n][d];
                }

                m_out[ni * m_npts + n] = amrex::max(
                    m_out[ni * m_npts + n],
                    amrex::ReduceMax(
                        vof, level_mask, 0,
                        [=] AMREX_GPU_HOST_DEVICE(
                            amrex::Box const& bx,
                            amrex::Array4<amrex::Real const> const& vof_arr,
                            amrex::Array4<int const> const& mask_arr)
                            -> amrex::Real {
                            amrex::Real height_fab = 0.0;

                            amrex::Loop(
                                bx,
                                [=, &height_fab](int i, int j, int k) noexcept {
                                    // Initialize height measurement
                                    amrex::Real ht = plo[dir];
                                    // Cell location
                                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
                                        xm;
                                    xm[0] = plo[0] + (i + 0.5) * dx[0];
                                    xm[1] = plo[1] + (j + 0.5) * dx[1];
                                    xm[2] = plo[2] + (k + 0.5) * dx[2];
                                    // (1) Check that cell height is below
                                    // previous instance. (2) Check if cell
                                    // contains 2D grid point: complicated
                                    // conditional is to avoid double-counting
                                    // and includes exception for lo boundary.
                                    // (3) Check if cell is obviously
                                    // multiphase, then check if cell might have
                                    // interface at top or bottom
                                    if ((dout_ptr[n] >
                                         xm[dir] + 0.5 * dx[dir]) &&
                                        (((plo[gc1] == loc[0] &&
                                           xm[gc1] - loc[0] == 0.5 * dx[gc1]) ||
                                          (xm[gc1] - loc[0] < 0.5 * dx[gc1] &&
                                           loc[0] - xm[gc1] <=
                                               0.5 * dx[gc1])) &&
                                         ((plo[gc2] == loc[1] &&
                                           xm[gc2] - loc[1] == 0.5 * dx[gc2]) ||
                                          (xm[gc2] - loc[1] < 0.5 * dx[gc2] &&
                                           loc[1] - xm[gc2] <=
                                               0.5 * dx[gc2])))) {

                                        int ip = i;
                                        int jp = j;
                                        int kp = k;
                                        int im = i;
                                        int jm = j;
                                        int km = k;
                                        amrex::Real mx = 0.0;
                                        amrex::Real my = 0.0;
                                        amrex::Real mz = 0.0;
                                        amrex::Real alpha = 1.0;
                                        // Get modified indices for checking up
                                        // and down and orient normal in search
                                        // direction
                                        switch (dir) {
                                        case 0:
                                            ip += 1;
                                            im -= 1;
                                            mx = 1.0;
                                            break;
                                        case 1:
                                            jp += 1;
                                            jm -= 1;
                                            my = 1.0;
                                            break;
                                        case 2:
                                            kp += 1;
                                            km -= 1;
                                            mz = 1.0;
                                            break;
                                        }
                                        // If cell is full of single phase
                                        // (accounts for when interface is at
                                        // intersection of cells but lower one
                                        // is not single-phase)
                                        bool calc_flag = false;
                                        if ((ni % 2 == 0 &&
                                             vof_arr(i, j, k) >= 1.0 - 1e-12) ||
                                            (ni % 2 == 1 &&
                                             vof_arr(i, j, k) <= 1e-12)) {
                                            // put bdy at top
                                            alpha = 1.0;
                                            if (ni % 2 == 1) {
                                                mx *= -1.0;
                                                my *= -1.0;
                                                mz *= -1.0;
                                                alpha *= -1.0;
                                            }
                                            calc_flag = true;
                                        }
                                        // Multiphase cell case
                                        if ((vof_arr(i, j, k) < (1.0 - 1e-12) &&
                                             vof_arr(i, j, k) > 1e-12)) {
                                            // Get interface reconstruction
                                            multiphase::fit_plane(
                                                i, j, k, vof_arr, mx, my, mz,
                                                alpha);
                                            calc_flag = true;
                                        }

                                        if (calc_flag) {
                                            // Reassign slope coefficients
                                            amrex::Real mdr = 0.0;
                                            amrex::Real mg1 = 0.0;
                                            amrex::Real mg2 = 0.0;
                                            switch (dir) {
                                            case 0:
                                                mdr = mx;
                                                mg1 = my;
                                                mg2 = mz;
                                                break;
                                            case 1:
                                                mdr = my;
                                                mg1 = mx;
                                                mg2 = mz;
                                                break;
                                            case 2:
                                                mdr = mz;
                                                mg1 = mx;
                                                mg2 = my;
                                                break;
                                            }
                                            // Get height of interface
                                            if (mdr == 0) {
                                                // If slope is undefined in z,
                                                // use middle of cell
                                                ht = xm[dir];
                                            } else {
                                                // Intersect 2D point with plane
                                                ht = (xm[dir] - 0.5 * dx[dir]) +
                                                     (alpha -
                                                      mg1 * dxi[gc1] *
                                                          (loc[0] -
                                                           (xm[gc1] -
                                                            0.5 * dx[gc1])) -
                                                      mg2 * dxi[gc2] *
                                                          (loc[1] -
                                                           (xm[gc2] -
                                                            0.5 * dx[gc2]))) /
                                                         (mdr * dxi[dir]);
                                            }
                                            // If interface is below lower
                                            // bound, continue to look
                                            if (ht < xm[dir] - 0.5 * dx[dir]) {
                                                ht = plo[dir];
                                            }
                                            // If interface is above upper
                                            // bound, limit it
                                            if (ht >
                                                xm[dir] + 0.5 * dx[dir] *
                                                              (1.0 + 1e-8)) {
                                                ht = xm[dir] + 0.5 * dx[dir];
                                            }
                                        }
                                    }
                                    // Offset by removing lo and contribute to
                                    // whole
                                    height_fab = amrex::max(
                                        height_fab,
                                        mask_arr(i, j, k) * (ht - plo[dir]));
                                });
                            return height_fab;
                        }));
            }
        }

        // Loop points in 2D grid
        for (int n = 0; n < m_npts; n++) {
            amrex::ParallelDescriptor::ReduceRealMax(m_out[ni * m_npts + n]);
        }
        // Add problo back to heights, making them absolute, not relative
        for (int n = 0; n < m_npts; n++) {
            m_out[ni * m_npts + n] += plo0[m_coorddir];
        }
        // Copy last m_out to device vector
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, &m_out[ni * m_npts],
            &m_out[(ni + 1) * m_npts - 1] + 1, dout.begin());
    }

    process_output();
}

void FreeSurface::process_output()
{
    if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("FreeSurface: Invalid output format encountered");
    }
}

void FreeSurface::write_ascii()
{
    BL_PROFILE("amr-wind::FreeSurface::write_ascii");
    amrex::Print()
        << "WARNING: FreeSurface: ASCII output will impact performance"
        << std::endl;

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    const std::string fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        //
        // Have I/O processor open file and write everything.
        //
        std::ofstream File;

        File.open(fname.c_str(), std::ios::out | std::ios::trunc);

        if (!File.good()) {
            amrex::FileOpenFailed(fname);
        }

        std::string str1 = "x";
        std::string str2 = "y";
        switch (m_coorddir) {
        case 0: {
            str1 = "y";
            str2 = "z";
            break;
        }
        case 1: {
            str1 = "x";
            str2 = "z";
            break;
        }
        default: {
            // Do nothing, initial settings are good
            break;
        }
        }

        // Metadata
        File << m_npts << '\n';
        File << str1 << ' ' << m_npts_dir[0] << ' ' << str2 << ' '
             << m_npts_dir[1] << '\n';

        // Points in grid (x, y, z0, z1, ...)
        for (int n = 0; n < m_npts; ++n) {
            File << m_locs[n][0] << ' ' << m_locs[n][1];
            for (int ni = 0; ni < m_ninst; ++ni) {
                File << ' ' << m_out[ni * m_npts + n];
            }
            File << '\n';
        }

        File.flush();

        File.close();

        if (!File.good()) {
            amrex::Abort("FreeSurface::write_ascii(): problem writing file");
        }
    }
}

void FreeSurface::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    m_ncfile_name = post_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string ngp_name = "num_grid_points";
    const std::string ninst_name = "num_instances";
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind data sampling output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(ngp_name, m_npts);
    ncf.def_dim(ninst_name, m_ninst);
    ncf.def_dim("ndim2", 2);
    ncf.def_var("time", NC_DOUBLE, {nt_name});

    // Metadata related to the 2D grid used to sample
    const std::vector<int> ijk{m_npts_dir[0], m_npts_dir[1]};
    ncf.put_attr("ijk_dims", ijk);
    ncf.put_attr("start", m_start);
    ncf.put_attr("end", m_end);

    // Set up array of data for locations in 2D grid
    ncf.def_var("coordinates2D", NC_DOUBLE, {ngp_name, "ndim2"});

    // Set up array for height outputs
    ncf.def_var("heights", NC_DOUBLE, {nt_name, ninst_name, ngp_name});

    ncf.exit_def_mode();

    // Populate data that doesn't change
    const std::vector<size_t> start{0, 0};
    std::vector<size_t> count{0, 2};
    count[0] = m_npts;
    auto xy = ncf.var("coordinates2D");
    xy.put(&m_locs[0][0], start, count);

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please "
        "recompile or "
        "use native format");
#endif
}

void FreeSurface::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF

    if (!amrex::ParallelDescriptor::IOProcessor()) return;
    auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of the next timestep
    const size_t nt = ncf.dim(nt_name).len();
    {
        auto time = m_sim.time().new_time();
        ncf.var("time").put(&time, {nt}, {1});
    }

    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, 0};

    count[1] = 1;
    count[2] = m_npts;
    auto var = ncf.var("heights");
    for (int ni = 0; ni < m_ninst; ++ni) {
        var.put(&m_out[ni * m_npts], start, count);
        ++start[1];
    }

    ncf.close();
#endif
}

} // namespace free_surface
} // namespace amr_wind
