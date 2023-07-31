/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMFIXEDBGVARS_HPP_
#define ADMFIXEDBGVARS_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

//! A structure for the decomposed elements of the Energy Momentum Tensor in
//! 3+1D - in GRChombo this is in CCZ4Geometry so it needs to be reproduced here
//! (Ideally it should just be moved to TensorAlgebra)
template <class data_t> struct emtensor_t
{
    Tensor<2, data_t> Sij; //!< S_ij = T_ij
    Tensor<1, data_t> Si;  //!< S_i = T_ia_n^a
    data_t S;              //!< S = S^i_i
    data_t rho;            //!< rho = T_ab n^a n^b
};

/// Namespace for ADM vars for fixed BG evolution
namespace ADMFixedBGVars
{
/// Vars object for ADM vars used in FixedBG evolution
template <class data_t> struct Vars
{
    // ADM vars needed in matter only rhs (ok for Proca and SF)
    Tensor<2, data_t> gamma;
    Tensor<2, data_t> K_tensor;
    data_t K;
    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<2, Tensor<1, data_t>> d1_gamma;
    Tensor<1, data_t> d1_lapse;
    Tensor<2, data_t> d1_shift;
};

} // namespace ADMFixedBGVars

#endif /* ADMFIXEDBGVARS_HPP_ */
