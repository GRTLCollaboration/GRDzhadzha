/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FRWSCALARLEVEL_HPP_
#define FRWSCALARLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "ScalarField.hpp"
#include "ScalarPotential.hpp"

//!  A class for the evolution of a scalar field, minimally coupled to gravity
/*!
     The class takes some initial data for a scalar field (variables phi and Pi)
     and evolves it using the CCZ4 equations. It is possible to specify an
   initial period of relaxation for the conformal factor chi, for non analytic
   initial conditions (for example, a general field configuration at a moment of
   time symmetry assuming conformal flatness). \sa MatterCCZ4(),
   ConstraintsMatter(), ScalarField(), RelaxationChi()
*/
class FRWScalarLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<FRWScalarLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // Typedef for scalar field
    typedef ScalarField<ScalarPotential> ScalarFieldWithPotential;

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    //! To do after each timestep
    virtual void specificPostTimeStep();

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#endif /* FRWSCALARLEVEL_HPP_ */
