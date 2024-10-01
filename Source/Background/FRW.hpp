/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FRW_HPP_
#define FRW_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes an FRW background


class FRW
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double rho0; // initial energy density
        double omega; // equation of state parameter = 0 (MD), 1/3 (RD)
        std::array<double, CH_SPACEDIM> center; 

    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;


    const params_t m_params;

    const double m_time;
    const double m_dx;
 

   FRW(params_t a_params, double a_dx, double a_time) : m_params(a_params), m_dx(a_dx), m_time(a_time){}


    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

       

        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;

     
        compute_metric_background(metric_vars, coords);
        



       //Reconstruct FRW quantities to check as diagnostics
       data_t scalefac = sqrt(metric_vars.gamma[0][0]);
       data_t hubbleparam = -metric_vars.K_tensor[0][0];

       current_cell.store_vars(scalefac, c_scalefac);
       current_cell.store_vars(hubbleparam, c_hubbleparam);
        

       
    }

    
    //Analytical background
     template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Coordinates<data_t> &coords) const
    {

        data_t scalefac; 
        data_t hubparam; 

        // Zero components - to do this in level?
        vars.K = 0.0;

        FOR(i) vars.shift[i] = 0.0;
        FOR(i,j) vars.gamma[i][j] = 0.0;
        FOR(i,j) vars.K_tensor[i][j] = 0.0;
        
        
        FOR(i,j, k) vars.d1_gamma[i][j][k] = 0.0;
        FOR(i) vars.d1_lapse[i] = 0.0;
        FOR(i,j) vars.d1_shift[i][j] = 0.0;

        // ds^2 = -dt^2 + a^2 dx_i dx^i assumin flat - TODO: curvature, conformal time

        vars.lapse = 1.0; // Conformal time - a^2

        if (m_params.rho0 == 0.0){

            scalefac = 1.0;

            FOR(i) vars.gamma[i][i] = scalefac; // Minkowski

        }
        
        
        else{

        data_t K0 = -sqrt(24.0*M_PI*m_params.rho0); //Ham constraint
            
        data_t t0 = -1.0/((1.0+m_params.omega)*K0); //Defining t0 from time evolution of scale factor
        
        
        
        scalefac = pow((m_time/t0)+1, 2/(3*(1+m_params.omega)));  //Scale factor evolution

        hubparam = 2.0/(3.0*(1+m_params.omega)*(m_time+t0)); //Hubble parameter evolution
        
        
        
        
        vars.K = -3.0*hubparam; //K = -3H (not used but included for completeness)
        
        FOR(i) vars.gamma[i][i] = pow(scalefac, 2.0); //gamma_ij = a^2

        FOR(i) vars.K_tensor[i][i] = -hubparam; //Kij = -H gamma_ij
        
        }
       

    }


   
    



};

#endif /* FRW_HPP_ */
