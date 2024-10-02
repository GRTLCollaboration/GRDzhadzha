/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "FRWScalarLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "FRW.hpp"
#include "InitialScalarData.hpp"
#include "MatterEvolution.hpp"
#include "ScalarField.hpp"
#include "ScalarPotential.hpp"

// Initial data for field and metric variables
void FRWScalarLevel::initialData()
{
    CH_TIME("FRWScalarLevel::initialData");
    if (m_verbosity)
        pout() << "FRWScalarLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for fields
    SetValue set_zero(0.0);
    FRW frw_bg(m_p.bg_params, m_dx, m_time);
    InitialScalarData initial_sf(m_p.initial_params, m_dx);
    auto compute_pack = make_compute_pack(set_zero, frw_bg);
    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);
    BoxLoops::loop(initial_sf, m_state_new, m_state_new, FILL_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void FRWScalarLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                     const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    ScalarPotential potential(m_p.initial_params);
    ScalarFieldWithPotential scalar_field(potential);
    FRW frw_bg(m_p.bg_params, m_dx, m_time);
    MatterEvolution<ScalarFieldWithPotential, FRW> my_matter(
        scalar_field, frw_bg, m_p.sigma, m_dx, m_p.center);
    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);
}

void FRWScalarLevel::specificPostTimeStep()
{
    // Check for nans on every level
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());

    // At any level, but after the min_level timestep
    int min_level = 0;
    bool calculate_diagnostics = at_level_timestep_multiple(min_level);
    if (calculate_diagnostics)
    {
        fillAllGhosts();
        ScalarPotential potential(m_p.initial_params);
        ScalarFieldWithPotential scalar_field(potential);
        FRW frw_bg(m_p.bg_params, m_dx, m_time);
        BoxLoops::loop(frw_bg, m_state_new, m_state_diagnostics,
                       SKIP_GHOST_CELLS);
    }

    // write out the integral after each timestep on minimum level
}

void FRWScalarLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                             const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion);
}
