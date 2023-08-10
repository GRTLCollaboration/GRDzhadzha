/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ANGULARMOMCONSERVATION_HPP_
#define ANGULARMOMCONSERVATION_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the angular momentum densities and fluxes and writes them to
//! the grid, assumes the axis of rotation is the z axis.
//! See https://arxiv.org/pdf/2203.13845.pdf
//! This also assumes a Killing vector in the background spacetime in phi
//! such that there is no source term
template <class matter_t, class background_t> class AngularMomConservation
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
    AngularMomConservation(matter_t a_matter, background_t a_background,
                           double a_dx,
                           std::array<double, CH_SPACEDIM> a_center)
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

        // Some useful quantities
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
        const data_t det_gamma = compute_determinant_sym(metric_vars.gamma);
        const data_t R = coords.get_radius();
        data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t r2sintheta = sqrt(rho2) * R;

        // the unit vector in the radial direction
        Tensor<1, data_t> si_L;
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;

        Tensor<1, data_t> dxdphi;
        dxdphi[0] = -coords.y;
        dxdphi[1] = coords.x;
        dxdphi[2] = 0;

        // angular momentum density
        data_t rhoAngMom = 0.0;
        FOR1(i) { rhoAngMom += emtensor.Si[i] * dxdphi[i]; }
        rhoAngMom *= sqrt(det_gamma);

        // flux density on surface
        data_t fluxAngMom = 0.0;
        FOR2(i, j)
        {
            fluxAngMom +=
                -si_L[i] * metric_vars.shift[i] * emtensor.Si[j] * dxdphi[j];
            FOR1(k)
            {
                fluxAngMom += metric_vars.lapse * emtensor.Sij[i][j] *
                              dxdphi[j] * gamma_UU[i][k] * si_L[k];
            }
        }

        // Add the volume factor to account for the spherical surface and
        // normal vector proper lengths
        fluxAngMom *= sqrt(det_gamma);

        // store values on the grid
        current_cell.store_vars(rhoAngMom, c_rhoAngMom);
        current_cell.store_vars(fluxAngMom, c_fluxAngMom);
    }
};

#endif /* ANGULARMOMCONSERVATION_HPP_ */
