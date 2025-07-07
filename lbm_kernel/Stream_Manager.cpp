/**
 * @file Stream_Model.h
 * @brief implement the stream model
 * @version 0.1
 * @date 2024-02-25
 * 
 */

#include "lbm_kernel/Stream_Manager.hpp"
#include "Lat_Manager.hpp"

/**
 * @note Read fCol write f, for respective stream and collision
*/
void ABPatternPush::stream(D_int level)
{
    D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
    _s_DDF* f_ptr = &(user_[level].df);

    for (auto lat_iter : *lat_ptr)
    {
        D_morton lat_code = lat_iter.first;
        f_ptr->f[0].at(lat_code) = f_ptr->fcol[0].at(lat_code);
        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, level, ex[i_q], ey[i_q], ez[i_q]);
            f_ptr->f[i_q].at(neighbor_code) = f_ptr->fcol[i_q].at(lat_code);
        }
    }

    for (auto lat_iter : Lat_Manager::pointer_me->lat_overlap_C2F.at(level)) {
        D_morton lat_code = lat_iter.first;
        f_ptr->f[0].at(lat_code) = f_ptr->fcol[0].at(lat_code);
        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, level, ex[i_q], ey[i_q], ez[i_q]);
            f_ptr->f[i_q].at(neighbor_code) = f_ptr->fcol[i_q].at(lat_code);
        }
    }
}

/**
 * @note Read fCol write f, for respective stream and collision
*/
void ABPatternPull::stream(D_int level)
{
    D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
    _s_DDF* f_ptr = &(user_[level].df);

    for (auto lat_iter : *lat_ptr)
    {
        D_morton lat_code = lat_iter.first;
        f_ptr->f[0].at(lat_code) = f_ptr->fcol[0].at(lat_code);
        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, level, ex[i_q], ey[i_q], ez[i_q]);
            f_ptr->f[i_q].at(lat_code) = f_ptr->fcol[i_q].at(neighbor_code);
        }
    }

    for (auto lat_iter : Lat_Manager::pointer_me->lat_overlap_C2F.at(level)) {
        D_morton lat_code = lat_iter.first;
        f_ptr->f[0].at(lat_code) = f_ptr->fcol[0].at(lat_code);
        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, level, ex[i_q], ey[i_q], ez[i_q]);
            f_ptr->f[i_q].at(lat_code) = f_ptr->fcol[i_q].at(neighbor_code);
        }
    }
}

/**
 * @note for Stream collision fusion, need to swap pointers between f and fCol
*/
// void ABCollisionStreamFusionPush::stream(D_int level)
// {
//     D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
//     _s_DDF* f_ptr = &(user_[level].df);

//     for (auto lat_iter : *lat_ptr)
//     {
//         D_morton lat_code = lat_iter.first;

//         // Collision part
//         D_Phy_Rho temp_rho = user_[level].density.at(lat_code);
        
//     }
// }

/**
 * @note for Stream collision fusion, need to swap pointers between f and fCol
*/
// class ABCollisionStreamFusionPull::stream(D_int level)
// {
//     D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
//     _s_DDF* f_ptr = &(user_[level].df);
// }

// void 