/**
 * @file Collision_Model.cpp
 * @brief  And implement some simple collision model, such as LBGK, MRT.
 *         Cumulant collision model is implemented in the Collision_Cumulant.cpp
 * @version 0.1
 * @date 2024-02-19
 * 
 */

#include "lbm_kernel/Collision_Manager.hpp"
#include "Lat_Manager.hpp"

void Collision_LBGK::collide(D_int i_level)
{
    // D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
    // D_Phy_real tau_alias = tau_[i_level];

    // auto collision_function = [i_level, &tau_alias] (D_mapLat* lattice_ptr,User* user, _s_DDF* f_ptr)
    auto collision_function = [i_level, this] (D_mapLat* lattice_ptr)
    {
        _s_DDF* f_ptr = &(user_[i_level].df);
        for (auto lat_iter : *lattice_ptr)
        {
            D_morton lat_code = lat_iter.first;

            D_Phy_DDF feq[C_Q];
            calculateFeq(user_[i_level].density.at(lat_code), user_[i_level].velocity.at(lat_code), feq);
            D_Phy_real tau0 = tau_[i_level];

            #ifdef LES
                D_Phy_real stressNorm = Turbulent_LES::calculate_stress_tensor(ddf, feq);
                tau0 = Turbulent_LES::smag_subgrid_tau(omega, stressNorm);
            #endif

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                f_ptr->fcol[i_q].at(lat_code) = f_ptr->f[i_q].at(lat_code) - (f_ptr->f[i_q].at(lat_code) - feq[i_q]) / tau0;
            }
        }
    };

    collision_function(&(Lat_Manager::pointer_me->lat_f.at(i_level)));

    collision_function(&(Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level)));
    
    for (auto lat_iter : Lat_Manager::pointer_me->lat_overlap_C2F.at(i_level))
    {
        D_morton lat_code = lat_iter.first;
        for (D_int i_q = 0; i_q < C_Q; ++i_q)
        {
            user_[i_level].df.fcol[i_q].at(lat_code) = user_[i_level].df.f[i_q].at(lat_code);
        }
    }
}

/**
 * @todo KBC collision
*/