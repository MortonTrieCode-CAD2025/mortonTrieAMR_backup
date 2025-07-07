// /**
//  * @file LBM_Partition.cpp
//  * @brief 
//  */

// #include "LBM_Partition.h"
// #include "General.h"
// #include "Grid_Manager.h"
// #include "LBM_Manager.h"
// #include "Lat_Manager.h"
// #include "user.h"

// DEPRECATED void LBM_Partition_1D::lbgkCollisionStreamFusion(D_int refine_level)
// {
    
// }

// DEPRECATED void LBM_Partition_1D::ABPatternStream(D_int level)
// {

// }

// // void LBM_Partition_1D::ibMethod(const std::string applyPosition)
// DEPRECATED void LBM_Partition_1D::ibMethod()
// {

// }

// DEPRECATED void LBM_Partition_1D::AMR_transDDF_coalescence(D_int level)
// {

// }

// DEPRECATED void LBM_Partition_1D::AMR_transDDF_explosion(D_int level)
// {

// }

// DEPRECATED void LBM_Partition_2D::lbgkCollisionStreamFusion(D_int level)
// {
//     D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
//     _s_DDF* f_ptr = &(user_[level].df);

//     for (auto lat_iter : *lat_ptr)
//     {
//         D_morton lat_code = lat_iter.first;

//         // Collision part
//         D_Phy_Rho temp_rho = user_[level].density.at(lat_code);
//         D_uvw temp_velocity =  user_[level].velocity.at(lat_code);
//         D_Phy_Velocity temp_u = temp_velocity.x, temp_v = temp_velocity.y, temp_w = temp_velocity.z;
//         D_Phy_real euv, feq;
//         D_Phy_real tau0 = tau_[level];
//         #ifdef LES
//             D_Phy_Rho CS = 0.16; // Smagorinsky constant
//             D_Phy_real Qxx = f_ptr->fcol[1].at(lat_code) + f_ptr->fcol[2].at(lat_code) + f_ptr->fcol[7].at(lat_code) + f_ptr->fcol[8].at(lat_code) + f_ptr->fcol[9].at(lat_code) + f_ptr->fcol[10].at(lat_code) + f_ptr->fcol[11].at(lat_code) + f_ptr->fcol[12].at(lat_code) + f_ptr->fcol[13].at(lat_code) + f_ptr->fcol[14].at(lat_code) - temp_u * temp_u - temp_rho / 3.;
//             D_Phy_real Qyy = f_ptr->fcol[3].at(lat_code) + f_ptr->fcol[4].at(lat_code) + f_ptr->fcol[7].at(lat_code) + f_ptr->fcol[8].at(lat_code) + f_ptr->fcol[9].at(lat_code) + f_ptr->fcol[10].at(lat_code) + f_ptr->fcol[15].at(lat_code) + f_ptr->fcol[16].at(lat_code) + f_ptr->fcol[17].at(lat_code) + f_ptr->fcol[18].at(lat_code) - temp_v * temp_v - temp_rho / 3.;
//             D_Phy_real Qzz = f_ptr->fcol[5].at(lat_code) + f_ptr->fcol[6].at(lat_code) + f_ptr->fcol[11].at(lat_code) + f_ptr->fcol[12].at(lat_code) + f_ptr->fcol[13].at(lat_code) + f_ptr->fcol[14].at(lat_code) + f_ptr->fcol[15].at(lat_code) + f_ptr->fcol[16].at(lat_code) + f_ptr->fcol[17].at(lat_code) + f_ptr->fcol[18].at(lat_code) - temp_w * temp_w - temp_rho / 3.;
//             D_Phy_real Qxy = f_ptr->fcol[7].at(lat_code) + f_ptr->fcol[8].at(lat_code) - f_ptr->fcol[9].at(lat_code) - f_ptr->fcol[10].at(lat_code) - temp_u * temp_v;
//             D_Phy_real Qyz = f_ptr->fcol[15].at(lat_code) + f_ptr->fcol[16].at(lat_code) - f_ptr->fcol[17].at(lat_code) - f_ptr->fcol[18].at(lat_code) - temp_v * temp_w;
//             D_Phy_real Qzx = f_ptr->fcol[11].at(lat_code) + f_ptr->fcol[12].at(lat_code) - f_ptr->fcol[13].at(lat_code) - f_ptr->fcol[14].at(lat_code) - temp_w * temp_u;

//             D_Phy_real Qave = sqrt(Qxx * Qxx + Qyy * Qyy + Qzz * Qzz + 2.0 * Qxy * Qxy + 2.0 * Qyz * Qyz + 2.0 * Qzx * Qzx);
//             tau0 = sqrt(tau0 * tau0 + 18.0 * CS * CS / two_power_n(2*level) * Qave / 1.0) / 2.0 + tau0 / 2.0;
//         #endif

//         D_Phy_real uu = temp_u * temp_u + temp_v * temp_v + temp_w * temp_w;

//         for (D_int i_q = 0; i_q < C_Q; ++i_q) {
//             euv = ex[i_q] * temp_u + ey[i_q] * temp_v + ez[i_q] * temp_w;
//             feq = temp_rho * w[i_q] * (1. + 3.*euv + 4.5*euv*euv - 1.5 * uu);
//             f_ptr->fcol[i_q].at(lat_code) = f_ptr->fcol[i_q].at(lat_code) - (f_ptr->fcol[i_q].at(lat_code) - feq) / tau0;
//         }
        
//         // Stream part
//         f_ptr->f[0].at(lat_code) = f_ptr->fcol[0].at(lat_code);

//         // No need to judge whether the neighbor node is on the same refinement level with i_grid,
//         // Even if the neighbor node is on the other level, the correspoding ddf is stored in the ghost layer between differenet resolution interface.

// #ifdef PUSH
//         f_ptr->f[1].at(Morton_Assist::find_x1(lat_code, level)) = f_ptr->fcol[1].at(lat_code);
//         f_ptr->f[2].at(Morton_Assist::find_x0(lat_code, level)) = f_ptr->fcol[2].at(lat_code);
//         f_ptr->f[3].at(Morton_Assist::find_y1(lat_code, level)) = f_ptr->fcol[3].at(lat_code);
//         f_ptr->f[4].at(Morton_Assist::find_y0(lat_code, level)) = f_ptr->fcol[4].at(lat_code);
//         f_ptr->f[5].at(Morton_Assist::find_z1(lat_code, level)) = f_ptr->fcol[5].at(lat_code);
//         f_ptr->f[6].at(Morton_Assist::find_z0(lat_code, level)) = f_ptr->fcol[6].at(lat_code);
//         f_ptr->f[7].at(Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[7].at(lat_code);
//         f_ptr->f[8].at(Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[8].at(lat_code);
//         f_ptr->f[9].at(Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[9].at(lat_code);
//         f_ptr->f[10].at(Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[10].at(lat_code);
//         f_ptr->f[11].at(Morton_Assist::find_z1(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[11].at(lat_code);
//         f_ptr->f[12].at(Morton_Assist::find_z0(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[12].at(lat_code);
//         f_ptr->f[13].at(Morton_Assist::find_z0(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[13].at(lat_code);
//         f_ptr->f[14].at(Morton_Assist::find_z1(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[14].at(lat_code);
//         f_ptr->f[15].at(Morton_Assist::find_z1(Morton_Assist::find_y1(lat_code, level), level)) = f_ptr->fcol[15].at(lat_code);
//         f_ptr->f[16].at(Morton_Assist::find_z0(Morton_Assist::find_y0(lat_code, level), level)) = f_ptr->fcol[16].at(lat_code);
//         f_ptr->f[17].at(Morton_Assist::find_z0(Morton_Assist::find_y1(lat_code, level), level)) = f_ptr->fcol[17].at(lat_code);
//         f_ptr->f[18].at(Morton_Assist::find_z1(Morton_Assist::find_y0(lat_code, level), level)) = f_ptr->fcol[18].at(lat_code);
// #endif

//     } // for i_lat

//     for (D_int i_q = 0; i_q < C_Q; ++i_q)
//         std::swap(user_[level].df.fcol[i_q], user_[level].df.f[i_q]);

//     std::cout << "level " << level << " lbgkCollisionStreamFusion Finish." << std::endl;    
// }

// DEPRECATED void LBM_Partition_2D::ABPatternStream(D_int level)
// {
//     D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
//     _s_DDF* f_ptr = &(user_[level].df);
//     std::cout << "level " << level << std::endl;
//     for (auto lat_iter : *lat_ptr)
//     {
//         D_morton lat_code = lat_iter.first;

//         f_ptr->f[0].at(lat_code) = f_ptr->fcol[0].at(lat_code);

//         // No need to judge whether the neighbor node is on the same refinement level with i_grid,
//         // Even if the neighbor node is on the other level, the correspoding ddf is stored in the ghost layer between differenet resolution interface.
// #if C_Q == 19
//     #ifdef PUSH
//             f_ptr->f[1].at(Morton_Assist::find_x1(lat_code, level)) = f_ptr->fcol[1].at(lat_code);
//             f_ptr->f[2].at(Morton_Assist::find_x0(lat_code, level)) = f_ptr->fcol[2].at(lat_code);
//             f_ptr->f[3].at(Morton_Assist::find_y1(lat_code, level)) = f_ptr->fcol[3].at(lat_code);
//             f_ptr->f[4].at(Morton_Assist::find_y0(lat_code, level)) = f_ptr->fcol[4].at(lat_code);
//             f_ptr->f[5].at(Morton_Assist::find_z1(lat_code, level)) = f_ptr->fcol[5].at(lat_code);
//             f_ptr->f[6].at(Morton_Assist::find_z0(lat_code, level)) = f_ptr->fcol[6].at(lat_code);
//             f_ptr->f[7].at(Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[7].at(lat_code);
//             f_ptr->f[8].at(Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[8].at(lat_code);
//             f_ptr->f[9].at(Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[9].at(lat_code);
//             f_ptr->f[10].at(Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[10].at(lat_code);
//             f_ptr->f[11].at(Morton_Assist::find_z1(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[11].at(lat_code);
//             f_ptr->f[12].at(Morton_Assist::find_z0(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[12].at(lat_code);
//             f_ptr->f[13].at(Morton_Assist::find_z0(Morton_Assist::find_x1(lat_code, level), level)) = f_ptr->fcol[13].at(lat_code);
//             f_ptr->f[14].at(Morton_Assist::find_z1(Morton_Assist::find_x0(lat_code, level), level)) = f_ptr->fcol[14].at(lat_code);
//             f_ptr->f[15].at(Morton_Assist::find_z1(Morton_Assist::find_y1(lat_code, level), level)) = f_ptr->fcol[15].at(lat_code);
//             f_ptr->f[16].at(Morton_Assist::find_z0(Morton_Assist::find_y0(lat_code, level), level)) = f_ptr->fcol[16].at(lat_code);
//             f_ptr->f[17].at(Morton_Assist::find_z0(Morton_Assist::find_y1(lat_code, level), level)) = f_ptr->fcol[17].at(lat_code);
//             f_ptr->f[18].at(Morton_Assist::find_z1(Morton_Assist::find_y0(lat_code, level), level)) = f_ptr->fcol[18].at(lat_code);
//     #endif

//     #ifdef PULL
//             f_ptr->f[1].at(lat_code) = f_ptr->fcol[1].at(Morton_Assist::find_x0(lat_code, level));
//             f_ptr->f[2].at(lat_code) = f_ptr->fcol[2].at(Morton_Assist::find_x1(lat_code, level));
//             f_ptr->f[3].at(lat_code) = f_ptr->fcol[3].at(Morton_Assist::find_y0(lat_code, level));
//             f_ptr->f[4].at(lat_code) = f_ptr->fcol[4].at(Morton_Assist::find_y1(lat_code, level));
//             f_ptr->f[5].at(lat_code) = f_ptr->fcol[5].at(Morton_Assist::find_z0(lat_code, level));
//             f_ptr->f[6].at(lat_code) = f_ptr->fcol[6].at(Morton_Assist::find_z1(lat_code, level));
//             f_ptr->f[7].at(lat_code) = f_ptr->fcol[7].at(Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, level), level));
//             f_ptr->f[8].at(lat_code) = f_ptr->fcol[8].at(Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, level), level));
//             f_ptr->f[9].at(lat_code) = f_ptr->fcol[9].at(Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, level), level));
//             f_ptr->f[10].at(lat_code) = f_ptr->fcol[10].at(Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, level), level));
//             f_ptr->f[11].at(lat_code) = f_ptr->fcol[11].at(Morton_Assist::find_z0(Morton_Assist::find_x0(lat_code, level), level));
//             f_ptr->f[12].at(lat_code) = f_ptr->fcol[12].at(Morton_Assist::find_z1(Morton_Assist::find_x1(lat_code, level), level));
//             f_ptr->f[13].at(lat_code) = f_ptr->fcol[13].at(Morton_Assist::find_z1(Morton_Assist::find_x0(lat_code, level), level));
//             f_ptr->f[14].at(lat_code) = f_ptr->fcol[14].at(Morton_Assist::find_z0(Morton_Assist::find_x1(lat_code, level), level));
//             f_ptr->f[15].at(lat_code) = f_ptr->fcol[15].at(Morton_Assist::find_z0(Morton_Assist::find_y0(lat_code, level), level));
//             f_ptr->f[16].at(lat_code) = f_ptr->fcol[16].at(Morton_Assist::find_z1(Morton_Assist::find_y1(lat_code, level), level));
//             f_ptr->f[17].at(lat_code) = f_ptr->fcol[17].at(Morton_Assist::find_z1(Morton_Assist::find_y0(lat_code, level), level));
//             f_ptr->f[18].at(lat_code) = f_ptr->fcol[18].at(Morton_Assist::find_z0(Morton_Assist::find_y1(lat_code, level), level));
//     #endif

// #endif // end if D3Q19
//     }



//     for (auto lat_iter : Lat_Manager::pointer_me->lat_overlap_C2F.at(level))
//     {
//         D_morton lat_code = lat_iter.first;

//         f_ptr->f[0].at(lat_code) = f_ptr->fcol[0].at(lat_code);
//         #if C_Q == 19
//             #ifdef PULL
//                     f_ptr->f[1].at(lat_code) = f_ptr->fcol[1].at(Morton_Assist::find_x0(lat_code, level));
//                     f_ptr->f[2].at(lat_code) = f_ptr->fcol[2].at(Morton_Assist::find_x1(lat_code, level));
//                     f_ptr->f[3].at(lat_code) = f_ptr->fcol[3].at(Morton_Assist::find_y0(lat_code, level));
//                     f_ptr->f[4].at(lat_code) = f_ptr->fcol[4].at(Morton_Assist::find_y1(lat_code, level));
//                     f_ptr->f[5].at(lat_code) = f_ptr->fcol[5].at(Morton_Assist::find_z0(lat_code, level));
//                     f_ptr->f[6].at(lat_code) = f_ptr->fcol[6].at(Morton_Assist::find_z1(lat_code, level));
//                     f_ptr->f[7].at(lat_code) = f_ptr->fcol[7].at(Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, level), level));
//                     f_ptr->f[8].at(lat_code) = f_ptr->fcol[8].at(Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, level), level));
//                     f_ptr->f[9].at(lat_code) = f_ptr->fcol[9].at(Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, level), level));
//                     f_ptr->f[10].at(lat_code)= f_ptr->fcol[10].at(Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, level), level));
//                     f_ptr->f[11].at(lat_code)= f_ptr->fcol[11].at(Morton_Assist::find_z0(Morton_Assist::find_x0(lat_code, level), level));
//                     f_ptr->f[12].at(lat_code)= f_ptr->fcol[12].at(Morton_Assist::find_z1(Morton_Assist::find_x1(lat_code, level), level));
//                     f_ptr->f[13].at(lat_code)= f_ptr->fcol[13].at(Morton_Assist::find_z1(Morton_Assist::find_x0(lat_code, level), level));
//                     f_ptr->f[14].at(lat_code)= f_ptr->fcol[14].at(Morton_Assist::find_z0(Morton_Assist::find_x1(lat_code, level), level));
//                     f_ptr->f[15].at(lat_code)= f_ptr->fcol[15].at(Morton_Assist::find_z0(Morton_Assist::find_y0(lat_code, level), level));
//                     f_ptr->f[16].at(lat_code)= f_ptr->fcol[16].at(Morton_Assist::find_z1(Morton_Assist::find_y1(lat_code, level), level));
//                     f_ptr->f[17].at(lat_code)= f_ptr->fcol[17].at(Morton_Assist::find_z1(Morton_Assist::find_y0(lat_code, level), level));
//                     f_ptr->f[18].at(lat_code)= f_ptr->fcol[18].at(Morton_Assist::find_z0(Morton_Assist::find_y1(lat_code, level), level));
//             #endif
//         #elif C_Q == 27
//             #ifdef PULL
//                 f_ptr->f[1].at(lat_code) = f_ptr->fcol[1].at(Morton_Assist::find_x0(lat_code, level));
//                 f_ptr->f[2].at(lat_code) = f_ptr->fcol[2].at(Morton_Assist::find_x1(lat_code, level));
//                 f_ptr->f[3].at(lat_code) = f_ptr->fcol[3].at(Morton_Assist::find_y0(lat_code, level));
//                 f_ptr->f[4].at(lat_code) = f_ptr->fcol[4].at(Morton_Assist::find_y1(lat_code, level));
//                 f_ptr->f[5].at(lat_code) = f_ptr->fcol[5].at(Morton_Assist::find_z0(lat_code, level));
//                 f_ptr->f[6].at(lat_code) = f_ptr->fcol[6].at(Morton_Assist::find_z1(lat_code, level));

//                 f_ptr->f[7].at(lat_code) = f_ptr->fcol[7].at(Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, level), level));
//                 f_ptr->f[8].at(lat_code) = f_ptr->fcol[8].at(Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, level), level));
//                 f_ptr->f[9].at(lat_code) = f_ptr->fcol[9].at(Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, level), level));
//                 f_ptr->f[10].at(lat_code)= f_ptr->fcol[10].at(Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, level), level));

//                 f_ptr->f[11].at(lat_code)= f_ptr->fcol[11].at(Morton_Assist::find_z0(Morton_Assist::find_x0(lat_code, level), level));
//                 f_ptr->f[12].at(lat_code)= f_ptr->fcol[12].at(Morton_Assist::find_z1(Morton_Assist::find_x1(lat_code, level), level));
//                 f_ptr->f[13].at(lat_code)= f_ptr->fcol[13].at(Morton_Assist::find_z1(Morton_Assist::find_x0(lat_code, level), level));
//                 f_ptr->f[14].at(lat_code)= f_ptr->fcol[14].at(Morton_Assist::find_z0(Morton_Assist::find_x1(lat_code, level), level));

//                 f_ptr->f[15].at(lat_code)= f_ptr->fcol[15].at(Morton_Assist::find_z0(Morton_Assist::find_y0(lat_code, level), level));
//                 f_ptr->f[16].at(lat_code)= f_ptr->fcol[16].at(Morton_Assist::find_z1(Morton_Assist::find_y1(lat_code, level), level));
//                 f_ptr->f[17].at(lat_code)= f_ptr->fcol[17].at(Morton_Assist::find_z1(Morton_Assist::find_y0(lat_code, level), level));
//                 f_ptr->f[18].at(lat_code)= f_ptr->fcol[18].at(Morton_Assist::find_z0(Morton_Assist::find_y1(lat_code, level), level));
                
//                 f_ptr->f[19].at(lat_code)= f_ptr->fcol[19].at(Morton_Assist::find_z0(Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, level), level), level));
//                 f_ptr->f[20].at(lat_code)= f_ptr->fcol[20].at(Morton_Assist::find_z0(Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, level), level), level));
//                 f_ptr->f[21].at(lat_code)= f_ptr->fcol[21].at(Morton_Assist::find_z0(Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, level), level), level));
//                 f_ptr->f[22].at(lat_code)= f_ptr->fcol[22].at(Morton_Assist::find_z0(Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, level), level), level));
//                 f_ptr->f[23].at(lat_code)= f_ptr->fcol[23].at(Morton_Assist::find_z1(Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, level), level), level));
//                 f_ptr->f[24].at(lat_code)= f_ptr->fcol[24].at(Morton_Assist::find_z1(Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, level), level), level));
//                 f_ptr->f[25].at(lat_code)= f_ptr->fcol[25].at(Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, level), level), level));
//                 f_ptr->f[26].at(lat_code)= f_ptr->fcol[26].at(Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, level), level), level));
//             #endif
//         #endif // end if D3Q19/D3Q27/..
//     }
// }

// // void LBM_Partition_2D::ibMethod(const std::string applyPosition)
// DEPRECATED void LBM_Partition_2D::ibMethod()
// {
//     auto dis_ptr = &(Lat_Manager::pointer_me->dis);
//     D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(C_max_level));

//     for (auto dis_ib : *dis_ptr)
//     {
//         D_morton lat_code = dis_ib.first;
        
//         // neighbor 1
//         D_real distance = dis_ib.second[0];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_x0(lat_code, C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[1].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[2].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[1].at(inv_code);
//             user_[C_max_level].df.f[2].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 2
//         distance = dis_ib.second[1];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_x1(lat_code, C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[2].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[1].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[2].at(inv_code);
//             user_[C_max_level].df.f[1].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 3
//         distance = dis_ib.second[2];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_y0(lat_code, C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[3].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[4].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[3].at(inv_code);
//             user_[C_max_level].df.f[4].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 4
//         distance = dis_ib.second[3];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_y1(lat_code, C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[4].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[3].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[4].at(inv_code);
//             user_[C_max_level].df.f[3].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 5
//         distance = dis_ib.second[4];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z0(lat_code, C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[5].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[6].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[5].at(inv_code);
//             user_[C_max_level].df.f[6].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 6
//         distance = dis_ib.second[5];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z1(lat_code, C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[6].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[5].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[6].at(inv_code);
//             user_[C_max_level].df.f[5].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 7
//         distance = dis_ib.second[6];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_y0(Morton_Assist::find_x0(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[7].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[8].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[7].at(inv_code);
//             user_[C_max_level].df.f[8].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 8
//         distance = dis_ib.second[7];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_y1(Morton_Assist::find_x1(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[8].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[7].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[8].at(inv_code);
//             user_[C_max_level].df.f[7].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 9
//         distance = dis_ib.second[8];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_y1(Morton_Assist::find_x0(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[9].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[10].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[9].at(inv_code);
//             user_[C_max_level].df.f[10].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 10
//         distance = dis_ib.second[9];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_y0(Morton_Assist::find_x1(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[10].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[9].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[10].at(inv_code);
//             user_[C_max_level].df.f[9].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 11
//         distance = dis_ib.second[10];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z0(Morton_Assist::find_x0(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[11].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[12].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[11].at(inv_code);
//             user_[C_max_level].df.f[12].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 12
//         distance = dis_ib.second[11];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z1(Morton_Assist::find_x1(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[12].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[11].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[12].at(inv_code);
//             user_[C_max_level].df.f[11].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 13
//         distance = dis_ib.second[12];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z1(Morton_Assist::find_x0(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[13].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[14].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[13].at(inv_code);
//             user_[C_max_level].df.f[14].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 14
//         distance = dis_ib.second[13];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z0(Morton_Assist::find_x1(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[14].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[13].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[14].at(inv_code);
//             user_[C_max_level].df.f[13].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 15
//         distance = dis_ib.second[14];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z0(Morton_Assist::find_y0(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[15].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[16].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[15].at(inv_code);
//             user_[C_max_level].df.f[16].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 16
//         distance = dis_ib.second[15];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z1(Morton_Assist::find_y1(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[16].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[15].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[16].at(inv_code);
//             user_[C_max_level].df.f[15].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 17
//         distance = dis_ib.second[16];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z1(Morton_Assist::find_y0(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[17].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[18].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[17].at(inv_code);
//             user_[C_max_level].df.f[18].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // neighbor 18
//         distance = dis_ib.second[17];
//         if (!isnanf(distance)) {
//             D_morton inv_code = Morton_Assist::find_z0(Morton_Assist::find_y1(lat_code, C_max_level),C_max_level);
//             D_Phy_DDF fcol = user_[C_max_level].df.fcol[18].at(lat_code);
//             D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[17].at(lat_code);
//             D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[18].at(inv_code);
//             user_[C_max_level].df.f[17].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
//         }

//         // Bounce-Back
//         // for (D_int i_q = 1; i_q < C_Q; ++i_q)
//         // {
//         //     if (!std::isnan(dis_ib.second[i_q-1]))
//         //         user_[C_max_level].df.f[i_q].at(lat_code) = user_[C_max_level].df.fcol[inv[i_q]].at(lat_code);
//         // }

//     }
// }

// DEPRECATED void LBM_Partition_2D::AMR_transDDF_coalescence(D_int level)
// {
//     D_int fine_level = level + 1;
//     D_int coarse_level = level;
//     for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(coarse_level))
//     {
//         D_morton ngbr_code[8]{i_lat.first, 
//                                Morton_Assist::find_x1(i_lat.first, fine_level), 
//                                 Morton_Assist::find_y1(i_lat.first, fine_level),
//                                  Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
//                                   Morton_Assist::find_z1(i_lat.first, fine_level),
//                                    Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
//                                     Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, fine_level), fine_level),
//                                      Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level), fine_level)};

//         for (D_int i_q = 0; i_q < C_Q; ++i_q) {
//             user_[coarse_level].df.f[i_q].at(i_lat.first) = 0;
//             for (D_int i_ngbr = 0; i_ngbr < 8; ++i_ngbr) {
//                 user_[coarse_level].df.f[i_q].at(i_lat.first) += user_[fine_level].df.f[i_q].at(ngbr_code[i_ngbr]);
//             }
//             user_[coarse_level].df.f[i_q].at(i_lat.first) = user_[coarse_level].df.f[i_q].at(i_lat.first) / 8;
//         }

//     }
// }

// /**
//  * @ref Massively parallel algorithms for the LBM on nonuniform grids SIAM J SCI COMPUT
//  * @param[in] level 
//  */
// DEPRECATED void LBM_Partition_2D::AMR_transDDF_explosion(D_int level)
// {
//     D_int fine_level = level;
//     D_int coarse_level = level - 1;

//     // for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(level)) 
//     for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(coarse_level)) 
//     {
//         D_morton coarse_code = i_lat.first;
//         D_morton fine_code_x1 = Morton_Assist::find_x1(coarse_code, fine_level);
//         D_morton fine_code_y1 = Morton_Assist::find_y1(coarse_code, fine_level);
//         D_morton fine_code_xy1 = Morton_Assist::find_y1(fine_code_x1, fine_level);
//         D_morton fine_code_z1 = Morton_Assist::find_z1(coarse_code, fine_level);
//         D_morton fine_code_zx1 = Morton_Assist::find_z1(fine_code_x1, fine_level);
//         D_morton fine_code_zy1 = Morton_Assist::find_z1(fine_code_y1, fine_level);
//         D_morton fine_code_xyz1 = Morton_Assist::find_z1(fine_code_xy1, fine_level);

//         for (D_int i_q = 0; i_q < C_Q; ++i_q)
//         {
//             user_[fine_level].df.fcol[i_q].at(coarse_code)   = user_[coarse_level].df.fcol[i_q].at(coarse_code);
//             user_[fine_level].df.fcol[i_q].at(fine_code_x1)  = user_[coarse_level].df.fcol[i_q].at(coarse_code);
//             user_[fine_level].df.fcol[i_q].at(fine_code_y1)  = user_[coarse_level].df.fcol[i_q].at(coarse_code);
//             user_[fine_level].df.fcol[i_q].at(fine_code_xy1) = user_[coarse_level].df.fcol[i_q].at(coarse_code);
//             user_[fine_level].df.fcol[i_q].at(fine_code_z1)  = user_[coarse_level].df.fcol[i_q].at(coarse_code);
//             user_[fine_level].df.fcol[i_q].at(fine_code_zx1) = user_[coarse_level].df.fcol[i_q].at(coarse_code);
//             user_[fine_level].df.fcol[i_q].at(fine_code_zy1) = user_[coarse_level].df.fcol[i_q].at(coarse_code);
//             user_[fine_level].df.fcol[i_q].at(fine_code_xyz1)= user_[coarse_level].df.fcol[i_q].at(coarse_code);
//         }
//     }

//     for (auto i_lat : Lat_Manager::pointer_me->ghost_overlap_F2C.at(coarse_level)) 
//     {
//         D_morton ngbr_code[8]{i_lat.first, 
//                                Morton_Assist::find_x1(i_lat.first, fine_level), 
//                                 Morton_Assist::find_y1(i_lat.first, fine_level),
//                                  Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
//                                   Morton_Assist::find_z1(i_lat.first, fine_level),
//                                    Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
//                                     Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, fine_level), fine_level),
//                                      Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level), fine_level)};

//         for (D_int i_ngbr = 0; i_ngbr < 8; ++i_ngbr) {
//             if (Lat_Manager::pointer_me->ghost_overlap_C2F.at(fine_level).find(ngbr_code[i_ngbr]) != Lat_Manager::pointer_me->ghost_overlap_C2F.at(fine_level).end()) {
//                 for (D_int i_q = 0; i_q < C_Q; ++i_q)
//                     user_[fine_level].df.fcol[i_q].at(ngbr_code[i_ngbr]) = user_[coarse_level].df.fcol[i_q].at(ngbr_code[0]);
//             }
//         }

//     }
// }

// #if C_DIMS == 3
// DEPRECATED void LBM_Partition_3D::lbgkCollisionStreamFusion(D_int refine_level)
// {
    
// }

// DEPRECATED void LBM_Partition_3D::ABPatternStream(D_int level)
// {

// }

// // void LBM_Partition_3D::ibMethod(const std::string applyPosition)
// DEPRECATED void LBM_Partition_3D::ibMethod()
// {
    
// }

// DEPRECATED void LBM_Partition_3D::AMR_transDDF_coalescence(D_int level)
// {

// }

// DEPRECATED void LBM_Partition_3D::AMR_transDDF_explosion(D_int level)
// {

// }

// #endif

// /**
//  * @skip Collision model has been encapsulated into the <Collision_Model> Class
//  *       This function "LBM_Partition::lbgkCollision" has been deprecated.
// */
// DEPRECATED void LBM_Partition::lbgkCollision(D_int level)
// {
//     D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
//     D_Phy_real tau_alias = tau_[level];

//     auto collision_kernel = [level, &tau_alias] (D_mapLat* lattice_ptr,User* user, _s_DDF* f_ptr)
//     {
//         for (auto lat_iter : *lattice_ptr)
//         {
//             D_morton lat_code = lat_iter.first;

//             user[level].density.at(lat_code) = f_ptr->f[0].at(lat_code);
//             user[level].v_old.at(lat_code) = user[level].velocity.at(lat_code);
//             user[level].velocity.at(lat_code) = {0.f, 0.f, 0.f};

//             for (D_int i_q = 1; i_q < C_Q; ++i_q)
//             {
//                 user[level].density.at(lat_code) += f_ptr->f[i_q].at(lat_code);
//                 user[level].velocity.at(lat_code).x += f_ptr->f[i_q].at(lat_code) * ex[i_q];
//                 user[level].velocity.at(lat_code).y += f_ptr->f[i_q].at(lat_code) * ey[i_q];
//                 user[level].velocity.at(lat_code).z += f_ptr->f[i_q].at(lat_code) * ez[i_q];
//             }
//             user[level].velocity.at(lat_code) = user[level].velocity.at(lat_code) / user[level].density.at(lat_code);

//             D_Phy_Rho temp_rho = user[level].density.at(lat_code);
//             D_uvw temp_velocity = user[level].velocity.at(lat_code);
//             D_Phy_Velocity temp_u = temp_velocity.x, temp_v = temp_velocity.y, temp_w = temp_velocity.z;
//             D_Phy_real euv, feq;
//             D_Phy_real tau0 = tau_alias;

//             #ifdef LES
//                 D_Phy_Rho CS = 0.16; // Smagorinsky constant
//                 D_Phy_real Qxx = f_ptr->f[1].at(lat_code) + f_ptr->f[2].at(lat_code) + f_ptr->f[7].at(lat_code)  + f_ptr->f[8].at(lat_code)  + f_ptr->f[9].at(lat_code)  + f_ptr->f[10].at(lat_code) + f_ptr->f[11].at(lat_code) + f_ptr->f[12].at(lat_code) + f_ptr->f[13].at(lat_code) + f_ptr->f[14].at(lat_code) - temp_u * temp_u - temp_rho / 3.;
//                 D_Phy_real Qyy = f_ptr->f[3].at(lat_code) + f_ptr->f[4].at(lat_code) + f_ptr->f[7].at(lat_code)  + f_ptr->f[8].at(lat_code)  + f_ptr->f[9].at(lat_code)  + f_ptr->f[10].at(lat_code) + f_ptr->f[15].at(lat_code) + f_ptr->f[16].at(lat_code) + f_ptr->f[17].at(lat_code) + f_ptr->f[18].at(lat_code) - temp_v * temp_v - temp_rho / 3.;
//                 D_Phy_real Qzz = f_ptr->f[5].at(lat_code) + f_ptr->f[6].at(lat_code) + f_ptr->f[11].at(lat_code) + f_ptr->f[12].at(lat_code) + f_ptr->f[13].at(lat_code) + f_ptr->f[14].at(lat_code) + f_ptr->f[15].at(lat_code) + f_ptr->f[16].at(lat_code) + f_ptr->f[17].at(lat_code) + f_ptr->f[18].at(lat_code) - temp_w * temp_w - temp_rho / 3.;
//                 D_Phy_real Qxy = f_ptr->f[7].at(lat_code) + f_ptr->f[8].at(lat_code) - f_ptr->f[9].at(lat_code)  - f_ptr->f[10].at(lat_code) - temp_u * temp_v;
//                 D_Phy_real Qyz = f_ptr->f[15].at(lat_code)+ f_ptr->f[16].at(lat_code)- f_ptr->f[17].at(lat_code) - f_ptr->f[18].at(lat_code) - temp_v * temp_w;
//                 D_Phy_real Qzx = f_ptr->f[11].at(lat_code)+ f_ptr->f[12].at(lat_code)- f_ptr->f[13].at(lat_code) - f_ptr->f[14].at(lat_code) - temp_w * temp_u;

//                 D_Phy_real Qave = sqrt(Qxx * Qxx + Qyy * Qyy + Qzz * Qzz + 2.0 * Qxy * Qxy + 2.0 * Qyz * Qyz + 2.0 * Qzx * Qzx);
//                 tau0 = sqrt(tau0 * tau0 + 18.0 * CS * CS / two_power_n(2*level) * Qave / 1.0) / 2.0 + tau0 / 2.0;
//             #endif

//             D_Phy_real uu = temp_u * temp_u + temp_v * temp_v + temp_w * temp_w;

//             for (D_int i_q = 0; i_q < C_Q; ++i_q) {
//                 euv = ex[i_q] * temp_u + ey[i_q] * temp_v + ez[i_q] * temp_w;
//                 feq = temp_rho * w[i_q] * (1. + 3.*euv + 4.5*euv*euv - 1.5 * uu);
//                 f_ptr->fcol[i_q].at(lat_code) = f_ptr->f[i_q].at(lat_code) - (f_ptr->f[i_q].at(lat_code) - feq) / tau0;
//             }
//         }
//     };

//     collision_kernel(lat_ptr, user_, &(user_[level].df));
//     collision_kernel(&(Lat_Manager::pointer_me->lat_overlap_F2C.at(level)), user_, &(user_[level].df));
//     for (auto lat_iter : Lat_Manager::pointer_me->lat_overlap_C2F.at(level))
//     {
//         D_morton lat_code = lat_iter.first;
//         for (D_int i_q = 0; i_q < C_Q; ++i_q)
//         {
//             user_[level].df.fcol[i_q].at(lat_code) = user_[level].df.f[i_q].at(lat_code);
//         }
//     }

// }

// DEPRECATED void LBM_Partition::calculateMacros(const D_Phy_DDF* ddf, D_Phy_Rho& rho, D_uvw& velocity)
// {
//     LBM_Manager::calculateMacros(ddf, rho, velocity);
// }