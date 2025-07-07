/**
 * @file les.cpp
 * @brief the implementation of the les turbulence model
 * @version 0.1
 * @date 2024-02-23
 * 
 */

#include "lbm_kernel/les.hpp"
#include "Lat_Manager.hpp"

inline D_Phy_real Turbulent_LES::calculate_stress_tensor(const D_Phy_DDF *ddf, const D_Phy_DDF *feq) 
{
    D_Phy_real mag = 0.0;
    D_Phy_DDF elem[9] = {0};

    for (D_int q = 0; q < C_Q; ++q) {
        D_Phy_DDF delta = ddf[q] - feq[q];
        elem[0] += ex[q] * ex[q] * delta;
        elem[1] += ex[q] * ey[q] * delta;
        elem[2] += ex[q] * ez[q] * delta;
        elem[3] += ey[q] * ex[q] * delta;
        elem[4] += ey[q] * ey[q] * delta;
        elem[5] += ey[q] * ez[q] * delta;
        elem[6] += ez[q] * ex[q] * delta;
        elem[7] += ez[q] * ey[q] * delta;
        elem[8] += ez[q] * ez[q] * delta;
    }

    for (int i = 0; i < 9; ++i) {
        mag += SQ(elem[i]);
    }

    return sqrt(mag);


}

inline D_Phy_real Turbulent_LES::smag_subgrid_tau(D_Phy_real tau, D_Phy_real stressNorm) 
{
    D_Phy_real vis = (tau - 0.5) / 3.0;
    D_Phy_real stress = (sqrt(SQ(vis) + 18 * sq_smag * stressNorm) - vis) / (6.0 * sq_smag);
    D_Phy_real tau0 = 3 * (vis + sq_smag * stress) + 0.5;
    return tau0;
}