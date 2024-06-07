#ifndef ASQUARED_HPP_
#define ASQUARED_HPP_


#include "ADMFixedBGVars.hpp"
#include "BoostedBH.hpp"
#include "Cell.hpp"
#include "ComplexProcaField.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"



template <class background_t> class ProcaSquared
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const std::array<double, CH_SPACEDIM> m_center;
    const background_t m_background;

    // Use the variable definition in ADMVars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    template <class data_t> using MatterVars = ComplexProcaField::template Vars<data_t>;

  public:
    ProcaSquared(double a_dx, const std::array<double, CH_SPACEDIM> a_center,
                 const background_t a_background)
        : m_dx(a_dx), m_deriv{m_dx}, m_center{a_center}, m_background{
                                                             a_background} {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // compute background variables
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        // load variables from Chombo grid
        const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto matter_vars_d1 =
            m_deriv.template diff1<MatterVars>(current_cell);

        // compute contravariant conformal spatial metric
        const auto gamma_UU{
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma)};

        data_t Asquared_re;
        Asquared_re = -matter_vars.Avec0_Re * matter_vars.Avec0_Re;
        FOR2(i, j)
        {
            Asquared_re +=
                gamma_UU[i][j] * matter_vars.Avec_Re[i] * matter_vars.Avec_Re[j];
        };

        data_t Asquared_im;
        Asquared_im = -matter_vars.Avec0_Im * matter_vars.Avec0_Im;
        FOR2(i, j)
        {
            Asquared_im +=
                gamma_UU[i][j] * matter_vars.Avec_Im[i] * matter_vars.Avec_Im[j];
        };
        
        data_t Asquared = Asquared_re + Asquared_im;

        current_cell.store_vars(Asquared, c_Asquared);
    };
};


#endif /* ASQUARED_HPP_ */