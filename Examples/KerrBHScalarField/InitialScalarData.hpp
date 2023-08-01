/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "ScalarField.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a constant scalar field given params for initial
//! matter config
class InitialScalarData
{
  public:
    struct params_t
    {
        double mass;
        double amplitude;
        std::array<double, CH_SPACEDIM> center;
    };

    //! The constructor for the class
    InitialScalarData(const params_t a_params, const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        ScalarField<>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.);

        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t radius = coords.get_radius();

        // set the field vars
        vars.phi = m_params.amplitude * exp((radius - 30.0) * (radius - 30.0) /
                                            0.1 / m_params.mass / m_params.mass);
        vars.Pi = 0;

        current_cell.store_vars(vars);
    }

  protected:
    const double m_dx;
    const params_t m_params;
};

#endif /* INITIALSCALARDATA_HPP_ */
