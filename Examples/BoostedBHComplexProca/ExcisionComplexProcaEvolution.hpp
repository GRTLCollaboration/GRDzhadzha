/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONPROCAEVOLUTION_HPP_
#define EXCISIONPROCAEVOLUTION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"

//! Does excision for fixed BG BH solutions
template <class matter_t, class background_t> class ExcisionProcaEvolution
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;

  protected:
    const double m_dx;             //!< The grid spacing
    const double m_z_exc_rad_mult; //!< Multiplier for radius of z excision.
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;

  public:
    ExcisionProcaEvolution(const double a_dx, const double a_z_exc_rad_mult,
                           const std::array<double, CH_SPACEDIM> a_center,
                           background_t a_background)
        : m_dx(a_dx), m_z_exc_rad_mult(a_z_exc_rad_mult), m_deriv(m_dx),
          m_center(a_center), m_background(a_background)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        bool is_excised = m_background.check_if_excised(coords);
        if (is_excised)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars vars;
            VarsTools::assign(vars, 0.0);
            // assign values of rhs or vars in output box
            current_cell.store_vars(vars);
        } // else do nothing

        // calculate outer horizon size
        const double bh_mass = m_background.m_params.mass;
        const double bh_vel = m_background.m_params.velocity;
        const double boost = pow(1.0 - bh_vel * bh_vel, -0.5);

        double r_plus = 0.5 * bh_mass;
        double x = coords.x;
        double y = coords.y;
        double z = coords.z;
        double radius = sqrt(x * x * boost * boost + y * y +
                             z * z); // boosted coordinate radius

        // Excise auxiliary Z field directly at horizon, since it can drive
        // errors
        if (radius < m_z_exc_rad_mult * r_plus)
        {
            current_cell.store_vars(0.0, c_Zvec_Re);
            current_cell.store_vars(0.0, c_Zvec_Im);
        }
    }
};

#endif /* EXCISIONPROCAEVOLUTION_HPP_ */
