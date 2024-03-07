/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALPROCADATA_HPP_
#define INITIALPROCADATA_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoostedBH.hpp"
#include "TensorAlgebra.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "ProcaField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class InitialProcaData
{
  protected:
    const double m_amplitude;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const BoostedBH::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t> using Vars = ProcaField::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialProcaData(const double a_amplitude, 
                     const std::array<double, CH_SPACEDIM> a_center,
                     const BoostedBH::params_t a_bg_params, const double a_dx)
        : m_dx(a_dx), m_amplitude(a_amplitude), m_center(a_center),
          m_bg_params(a_bg_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        BoostedBH boosted_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        boosted_bh.compute_metric_background(metric_vars, coords);
        
        //Compute contravariant spatial metric
        Tensor<2,data_t> gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);

        // set the field variables to constant values
        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);

        // set the vector values
        vars.Avec[0] = m_amplitude;
        vars.Avec[1] = 0.0;
        vars.Avec[2] = 0.0;
        
        // set auxiliary field values
        vars.Zvec = 0.0;

        // set electric field values consistent with Avec and Avec0
        FOR1(i)
        {
          vars.Evec[i] = 0.0;
          FOR3(j,k,l)
          {
            vars.Evec[i] += gamma_UU[i][k] * gamma_UU[j][l] * vars.Avec[j] * metric_vars.K_tensor[k][l];
          }
        }

        // set scalar part of Proca field
        vars.Avec0 = 0.0;

        // Store the initial values of the variables
        current_cell.store_vars(vars);
    }
};

#endif /* INITIALPROCADATA_HPP_ */
