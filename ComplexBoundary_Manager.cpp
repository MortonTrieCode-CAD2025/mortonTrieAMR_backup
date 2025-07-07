/**
 * @file SolidStream_Manager.h
 * @brief fluid stream aroung obstacle
 * @version 0.1
 * @date 2024-02-26
 * 
 */

#include "ComplexBoundary_Manager.hpp"
#include "Lat_Manager.hpp"
#include "Morton_Assist.h"
#include "LBM_Manager.hpp"

void UniformBoundary::treat()
{
    auto dis_ptr = &(Lat_Manager::pointer_me->dis);
    D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(C_max_level));

    for (auto dis_ib : *dis_ptr)
    {
        D_morton lat_code = dis_ib.first;

        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            // D_real distance = dis_ib.second[i_q-1];
            D_real distance = dis_ib.second[e_inv[i_q]-1];
            if (distance != -1) {
                D_morton inv_code = Morton_Assist::find_neighbor(lat_code, C_max_level, ex[e_inv[i_q]], ey[e_inv[i_q]], ez[e_inv[i_q]]);
                D_Phy_DDF fcol = user_[C_max_level].df.fcol[i_q].at(lat_code);
                D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[e_inv[i_q]].at(lat_code);
                D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[i_q].at(inv_code);
                user_[C_max_level].df.f[e_inv[i_q]].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
            }
            else {
                D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, C_max_level, ex[i_q], ey[i_q], ez[i_q]);
                user_[C_max_level].df.f[i_q].at(lat_code) = user_[C_max_level].df.fcol[i_q].at(neighbor_code);
            }
        }
    }

}


void SinglePointInterpolation::treat()
{
    auto dis_ptr = &(Lat_Manager::pointer_me->dis);
    D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(C_max_level));
    D_real dx_inner = Lat_Manager::pointer_me->dx.at(C_max_level);

    for (auto dis_ib : *dis_ptr)
    {
        D_morton lat_code = dis_ib.first;
        D_vec lat_center;
        Morton_Assist::compute_coordinate(lat_code, C_max_level, lat_center.x, lat_center.y, lat_center.z);
        lat_center += (D_vec){dx_inner/2., dx_inner/2., dx_inner/2.};
        D_real CF_U_ = LBM_Manager::pointer_me->CF_U;
        D_uvw u_b = user_[C_max_level].velocity.at(lat_code);
        D_Phy_Rho rho_b = user_[C_max_level].density.at(lat_code);

        user_[C_max_level].df.f[0].at(lat_code) = user_[C_max_level].df.fcol[0].at(lat_code);

        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            D_real distance = dis_ib.second[e_inv[i_q]-1];

            if (distance != -1.) {

                D_vec vecR = (D_vec){0.,0.,0} - lat_center;
                D_vec angVel(0.,0.,0);
                D_vec linVel = cross_product(angVel, vecR);
                D_vec cartVel(0.,0.,0.);
                D_uvw uw((linVel.x + cartVel.x)/CF_U_, (linVel.y + cartVel.y)/CF_U_, (linVel.z + cartVel.z)/CF_U_);

                // complex boundary treatment
                D_Phy_DDF feqi_ = LBM_Manager::calculateFeq(rho_b, uw, i_q);
                D_Phy_DDF fi = user_[C_max_level].df.f[e_inv[i_q]].at(lat_code);
                D_Phy_DDF fneqi = fi - LBM_Manager::calculateFeq(rho_b, u_b, e_inv[i_q]);
                D_Phy_DDF fcoli_ = user_[C_max_level].df.fcol[i_q].at(lat_code);

                user_[C_max_level].df.f[i_q].at(lat_code) = (feqi_ + fneqi + distance * fcoli_) / (1. + distance);
            }
            else {
                // normal stream in pull mode
                D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, C_max_level, ex[i_q], ey[i_q], ez[i_q]);
                user_[C_max_level].df.f[i_q].at(lat_code) = user_[C_max_level].df.fcol[i_q].at(neighbor_code);
            }

        } // for i_q

    } // for distance
}