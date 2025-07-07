/**
 * @file 
*/

#pragma once

#include "General.h"

class Turbulent_LES
{
private:
    static constexpr double smag = 0.04; // Smagorinsky constant
    static constexpr double sq_smag = SQ(smag);

public:
    static D_Phy_real calculate_stress_tensor(const D_Phy_DDF *ddf, const D_Phy_DDF *feq);
    static D_Phy_real smag_subgrid_tau(D_Phy_DDF tau, D_Phy_DDF stressNorm);
};