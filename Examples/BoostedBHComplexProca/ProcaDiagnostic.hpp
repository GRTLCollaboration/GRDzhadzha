#ifndef PROCAGAUSSDIAGNOSTIC_HPP_INCLUDED
#define PROCAGAUSSDIAGNOSTIC_HPP_INCLUDED

#include "ADMFixedBGVars.hpp"
#include "BoostedBH.hpp"
#include "Cell.hpp"
#include "ComplexProcaField.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

class ProcaGaussDiagnostic
{
  protected:
    template <class data_t>
    using MatterVars =
        ComplexProcaField::Vars<data_t>; // type alias the matter variables

  public:
    ProcaGaussDiagnostic(){}; // explicit default constructor

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // Load the auxiliary Z field and save it to diagnostic variable
        // c_Zvec_out. The auxiliary Z field is equal to minus the Gauss
        // constraint. so instead of calculating the Gauss constraint directly,
        // we just have to load the Z field

        MatterVars<data_t> matter_vars{current_cell.template load_vars<
            MatterVars>()}; // load the matter variables from the Chombo grid
        current_cell.store_vars(matter_vars.Zvec_Re, c_Zvec_Re_out);
        current_cell.store_vars(matter_vars.Zvec_Im, c_Zvec_Im_out);
    };
};

#endif /* PROCAGAUSSDIAGNOSTIC_HPP_INCLUDED */