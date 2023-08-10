/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ENERGYCONSERVATION_HPP_
#define ENERGYCONSERVATION_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the energy integrands and fluxes with type matter_t
//! and writes them to the grid, as in https://arxiv.org/pdf/2104.13420.pdf
//! NB This assumes a background with a time like Killing Vector
template <class matter_t, class background_t> class EnergyConservation
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;                        //!< The matter object
    const double m_dx;                              //!< The grid spacing
    const background_t m_background;                //!< The metric background
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center

  public:
    EnergyConservation(matter_t a_matter, background_t a_background,
                       double a_dx, std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

        // get the metric vars from the background
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        // some useful quantities
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
        const auto lapse = metric_vars.lapse;
        const auto shift = metric_vars.shift;
        const data_t det_gamma = compute_determinant_sym(metric_vars.gamma);

        const data_t R = coords.get_radius();
        data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t r2sintheta = sqrt(rho2) * R;

        // the unit coordinate vector in the radial direction
        Tensor<1, data_t> si_L;
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;

        // See eqn (10) - note the sign is reversed to give a positive rho
        data_t rhoEnergy = emtensor.rho * metric_vars.lapse;
        FOR1(i) { rhoEnergy += -emtensor.Si[i] * metric_vars.shift[i]; }
        rhoEnergy *= sqrt(det_gamma);

        // See eqn (11) - note sign reversal
        data_t fluxEnergy = 0.0;
        FOR1(i)
        {
            fluxEnergy += -lapse * si_L[i] * emtensor.rho * shift[i];
            FOR1(j)
            {
                fluxEnergy +=
                    si_L[i] * emtensor.Si[j] *
                    (shift[i] * shift[j] + lapse * lapse * gamma_UU[i][j]);
                FOR1(k)
                {
                    fluxEnergy += -si_L[i] * lapse * gamma_UU[i][j] * shift[k] *
                                  emtensor.Sij[j][k];
                }
            }
        }

        // Add the volume factor to account for the spherical surface and
        // normal vector proper lengths
        fluxEnergy *= sqrt(det_gamma);

        // Store the values on the grid
        current_cell.store_vars(rhoEnergy, c_rhoEnergy);
        current_cell.store_vars(fluxEnergy, c_fluxEnergy);
    }
};

#endif /* ENERGYCONSERVATION_HPP_ */
