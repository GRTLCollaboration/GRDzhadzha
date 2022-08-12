/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOOSTEDBHSCALARLEVEL_HPP_
#define BOOSTEDBHSCALARLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "ComplexScalarField.hpp"
#include "ComplexScalarPotential.hpp"

//!  A class for the evolution of a scalar field, minimally coupled to gravity
/*!
     The class takes some initial data for a scalar field (variables phi and Pi)
     and evolves it on a fixed metric background.
*/
class BoostedBHScalarLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<BoostedBHScalarLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // Typedef for scalar field
    typedef ComplexScalarField<ComplexScalarPotential> ScalarFieldWithPotential;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance();

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    // to do post each time step on every level
    virtual void specificPostTimeStep();

    // Any actions that should happen just before plot files output
    virtual void prePlotLevel();

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#endif /* BOOSTEDBHSCALARLEVEL_HPP_ */
