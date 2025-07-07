/**
 * @file Collision_Model.h
 * @brief Define collision model for the LBM kernel in the singleton module
 * @version 0.1
 * @date 2024-02-13
 * 
 */

#pragma once

#include <string>
#include "General.h"
#include "user.h"
#include "Lat_Manager.hpp"
#include "LBM_Manager.hpp"

class Collision_Model 
{
protected:
    User* user_;
    const std::array<D_Phy_real,C_max_level+1> tau_;

    void calculateFeq(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq) { 
        LBM_Manager::calculateFeq(rho, velocity, feq); 
    }

    void calculateFeq_HeLuo(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq) { 
        LBM_Manager::calculateFeq_HeLuo(rho, velocity, feq); 
    }

public:
    // Collision_Model(User* user, const std::array<D_Phy_real,C_max_level+1> tau) 
    Collision_Model() 
     : user_(LBM_Manager::pointer_me->user), tau_(LBM_Manager::pointer_me->tau) {}

    virtual void collide(D_int level) = 0;

    virtual std::string print_Model_Name() = 0;

    virtual ~Collision_Model() {}
};

/**
 * @brief the basic collision model: LBGK
*/
class Collision_LBGK : public Collision_Model
{
private:
    // Private constructor to prevent external instantiation
    friend class Collision_Model_Factory;
    // Collision_LBGK(User* user, const std::array<D_Phy_real,C_max_level+1> tau) : Collision_Model(user, tau) {}
    // Collision_LBGK() : Collision_Model() {}

public:
    // Static method to get the singleton instance
    // static Collision_LBGK* getInstance(User* user, const std::array<D_Phy_real,C_max_level+1> tau) 
    static Collision_LBGK* getInstance() 
    {
        // static Collision_LBGK instance(user, tau);
        static Collision_LBGK instance;
        return &instance;
    }
    
    void collide(D_int level) override;

    std::string print_Model_Name() override { return "LBGK"; }

    ~Collision_LBGK() override {}
};

