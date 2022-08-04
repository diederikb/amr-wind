#ifndef ConvectingTaylorVortex_H
#define ConvectingTaylorVortex_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {
namespace ctv {

namespace {
struct UExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real /*u0*/,
        const amrex::Real /*v0*/,
        const amrex::Real /*alpha*/,
        const amrex::Real /*beta*/,
        const amrex::Real /*A*/,
        const amrex::Real /*nu*/,
        const amrex::Real /*x*/,
        const amrex::Real /*y*/,
        const amrex::Real /*t*/) const;
    const int m_comp{0};
};

struct VExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real /*u0*/,
        const amrex::Real /*v0*/,
        const amrex::Real /*alpha*/,
        const amrex::Real /*beta*/,
        const amrex::Real /*A*/,
        const amrex::Real /*nu*/,
        const amrex::Real /*x*/,
        const amrex::Real /*y*/,
        const amrex::Real /*t*/) const;
    const int m_comp{1};
};

struct WExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/) const;
    const int m_comp{2};
};

struct GpxExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real /*u0*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*alpha*/,
        const amrex::Real /*beta*/,
        const amrex::Real /*A*/,
        const amrex::Real /*nu*/,
        const amrex::Real /*x*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*t*/) const;
    const int m_comp{0};
};

struct GpyExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real /*unused*/,
        const amrex::Real /*v0*/,
        const amrex::Real /*alpha*/,
        const amrex::Real /*beta*/,
        const amrex::Real /*A*/,
        const amrex::Real /*nu*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*y*/,
        const amrex::Real /*t*/) const;
    const int m_comp{1};
};

struct GpzExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/) const;
    const int m_comp{2};
};

} // namespace

/** Convecting Taylor Vortex physics
 *  \ingroup physics
 */
class ConvectingTaylorVortex : public Physics::Register<ConvectingTaylorVortex>
{
public:
    static std::string identifier() { return "ConvectingTaylorVortex"; }

    explicit ConvectingTaylorVortex(const CFDSim& sim);

    ~ConvectingTaylorVortex() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    template <typename T>
    amrex::Real compute_error(const Field& /*field*/);

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

    amrex::Real get_u0() const { return m_u0; }

    amrex::Real get_v0() const { return m_v0; }

    amrex::Real get_nu() const { return m_nu; }

private:
    const amr_wind::SimTime& m_time;
    const amr_wind::CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    Field& m_velocity;
    Field& m_gradp;
    Field& m_density;

    void output_error();

    //! initial density value
    amrex::Real m_rho{1.0};

    //! mean x-velocity value
    amrex::Real m_u0{1.0};

    //! mean y-velocity value
    amrex::Real m_v0{1.0};

    //! x-direction wave number
    amrex::Real m_alpha{utils::pi()};

    //! y-direction wave number
    amrex::Real m_beta{utils::pi()};

    //! amplitude
    amrex::Real m_A{1.0/utils::pi()};

    //! viscosity
    amrex::Real m_nu;

    //! output precision
    const int m_w{18};

    //! error log file
    std::string m_output_fname{"ctv.log"};

    bool m_activate_pressure{false};
    bool m_mesh_mapping{false};
};
} // namespace ctv
} // namespace amr_wind

#endif /* ConvectingTaylorVortex_H */
