/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOOSTEDBHPROCALEVEL_HPP_
#define BOOSTEDBHPROCALEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "ProcaField.hpp"

//!  A class for the evolution of a proca field, minimally coupled to gravity
/*!
     The class takes some initial data for a proca field (variables phi and Pi)
     and evolves it on a fixed metric background.
*/
class BoostedBHProcaLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<BoostedBHProcaLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // Typedef for proca field
    typedef ProcaField ProcaField;

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    // to do post each time step on every level
    virtual void specificPostTimeStep();

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#endif /* BOOSTEDBHPROCALEVEL_HPP_ */
