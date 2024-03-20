/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROCACONSTRAINT_HPP_
#define PROCACONSTRAINT_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoostedBH.hpp"
#include "Cell.hpp"
#include "ComplexProcaField.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class ProcaConstraint
{
  protected:
    const double m_dx;
    const double m_mu;
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const std::array<double, CH_SPACEDIM> m_center;
    const BoostedBH::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t> using Vars = ComplexProcaField::Vars<data_t>;

  public:
    //! The constructor for the class
    ProcaConstraint(const std::array<double, CH_SPACEDIM> a_center,
                    const BoostedBH::params_t a_bg_params, const double a_mu,
                    const double a_dx)
        : m_center(a_center), m_bg_params(a_bg_params), m_mu(a_mu), m_dx(a_dx),
          m_deriv(a_dx)
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

        // calculate full spatial christoffel symbols
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        // The matter vars
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        // The constraint sets Avec0 from the derivs of Evec
        data_t Avec0_Re = 0.0;
        data_t Avec0_Im = 0.0;
        FOR1(i)
        {
            Avec0_Re += -d1.Evec_Re[i][i];
            Avec0_Im += -d1.Evec_Im[i][i];

            FOR1(j)
            {
                Avec0_Re += -chris_phys.ULL[i][i][j] * vars.Evec_Re[j];
                Avec0_Im += -chris_phys.ULL[i][i][j] * vars.Evec_Im[j];
            }
        }
        Avec0_Re = Avec0_Re / m_mu / m_mu;
        Avec0_Im = Avec0_Im / m_mu / m_mu;

        // Store the initial values of the variables
        current_cell.store_vars(Avec0_Re, c_Avec0_Re);
        current_cell.store_vars(Avec0_Im, c_Avec0_Im);
    }
};

#endif /* PROCACONSTRAINTS_HPP_ */
