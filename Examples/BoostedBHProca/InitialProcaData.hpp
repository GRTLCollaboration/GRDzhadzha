/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALPROCADATA_HPP_
#define INITIALPROCADATA_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoostedBH.hpp"
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
    const double m_mu;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const BoostedBH::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t> using Vars = ProcaField::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialProcaData(const double a_amplitude, const double a_mu,
                     const std::array<double, CH_SPACEDIM> a_center,
                     const BoostedBH::params_t a_bg_params, const double a_dx)
        : m_dx(a_dx), m_amplitude(a_amplitude), m_mu(a_mu), m_center(a_center),
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
        // const auto metric_vars = current_cell.template
        // load_vars<MetricVars>();
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        // useful coordinate quantities
        const double M = m_bg_params.mass;
        const double z = coords.z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;

        // the peak of the bound state for the n=l=0, j=m=1 state
        // See arXiv : 1704.05081 eqns 9 and 10
        double alpha = M * m_mu;
        double r_0 = 1.0 / (m_mu * alpha);

        // set the field variable to approx profile
        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        data_t Avec = m_amplitude * (exp(-r / r_0)) / det_gamma;
        double coswt = 1.0;
        double sinwt = 0.0;

        // set the vector values
        vars.Avec[0] = -Avec * coswt;
        vars.Avec[1] = -Avec * sinwt;
        vars.Avec[2] = 0.0;

        // Store the initial values of the variables
        current_cell.store_vars(vars);
    }
};

#endif /* INITIALPROCADATA_HPP_ */
