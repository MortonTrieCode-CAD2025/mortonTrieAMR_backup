/**
 * @file Collision_Cumulant.h
 * @brief the implementation of the cumulant collsion
 * @version 0.1
 * @date 2024-02-13
 * 
 */

#pragma once

#include "General.h"
#include "user.h"
#include "Collision_Manager.hpp"


class Collision_Cumulant : public Collision_Model
{
private:
    // Private constructor to prevent external instantiation
    friend class Collision_Model_Factory;
    // Collision_Cumulant(User* user, const std::array<D_Phy_real,C_max_level+1> tau) : Collision_Model(user, tau) {}
    // Collision_Cumulant() : Collision_Model() {}

    void fwd_central_moment_trans(const D_Phy_DDF* ddf, const D_uvw velocity, D_Phy_DDF* k);
    void fwd_cum_trans(const D_Phy_DDF *k, const D_Phy_Rho rho, D_Phy_DDF *C);
    void cum_collide_kernel(D_Phy_DDF *C, const D_Phy_DDF *k, D_Phy_real tau, const D_uvw velocity, const D_Phy_Rho rho);
    void bwd_cum_trans(const D_Phy_DDF *C, D_Phy_DDF *k, D_Phy_Rho rho);
    void bwd_central_moment_trans(const D_Phy_DDF *k, D_Phy_DDF *ddf, const D_uvw velocity);

public:
    // Static method to get the singleton instance
    // static Collision_Cumulant* getInstance(User* user, const std::array<D_Phy_real,C_max_level+1> tau) 
    static Collision_Cumulant* getInstance() 
    {
        // static Collision_Cumulant instance(user, tau);
        static Collision_Cumulant instance;
        return &instance;
    }

    std::string print_Model_Name() override { return "Cumulant"; }
    
    void collide(D_int level) override;

};