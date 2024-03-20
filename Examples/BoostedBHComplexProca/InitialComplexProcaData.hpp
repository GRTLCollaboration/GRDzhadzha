/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALPROCADATA_HPP_
#define INITIALPROCADATA_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoostedBH.hpp"
#include "Cell.hpp"
#include "ComplexProcaField.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class InitialProcaData
{
  protected:
    const double m_amplitude;
    const double m_mass;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const BoostedBH::params_t m_bg_params;
    const std::string m_proca_initial_data_profile;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t> using Vars = ComplexProcaField::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialProcaData(
        const double a_amplitude, const double a_mass,
        const std::array<double, CH_SPACEDIM> a_center,
        const BoostedBH::params_t a_bg_params, const double a_dx,
        const std::string a_proca_initial_data_profile = "uniform-x")
        : m_dx(a_dx), m_amplitude(a_amplitude), m_mass{a_mass},
          m_center(a_center), m_bg_params(a_bg_params),
          m_proca_initial_data_profile(a_proca_initial_data_profile)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        /*
        We choose Avec0 = 0 and set Evec=0 to satisfy gauss constraint
        */
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        BoostedBH boosted_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        boosted_bh.compute_metric_background(metric_vars, coords);

        // Compute contravariant spatial metric
        Tensor<2, data_t> gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);

        // set the field variables default values to zero
        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);

        // set the vector values
        if (m_proca_initial_data_profile == "uniform-x")
        {
            vars.Avec_Re[0] = m_amplitude;
        }
        else if (m_proca_initial_data_profile == "uniform-y")
        {
            vars.Avec_Re[1] = m_amplitude;
        }
        else if (m_proca_initial_data_profile == "uniform-z")
        {
            vars.Avec_Re[2] = m_amplitude;
        }

        // Since Evec is minus the conjugate momentum to Avec, we set E = -
        // dAvec/dt
        //  Assuming Avec is a phasor, i.e. Avec[i] ~ e^(-i*w*t)
        FOR1(i)
        {
            vars.Evec_Im[i] = 0;

            FOR1(j)
            {
                vars.Evec_Im[i] += -gamma_UU[i][j] * m_mass * vars.Avec_Re[j];
            }
        }

        // Store the initial values of the variables
        current_cell.store_vars(vars);
    }
};

#endif /* INITIALPROCADATA_HPP_ */
