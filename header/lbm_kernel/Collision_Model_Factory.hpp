/**
 * @file Collision_Model_Factory.h
 * @brief Define collision model for the LBM kernel in the singleton module
 * @version 0.1
 * @date 2024-03-31
 * 
 */


#pragma once

#include "Collision_Manager.hpp"
#include "Collision_Cumulant.hpp"

class Collision_Model_Factory
{
public:
    // Static method to create a Collision_Model_Factory instance based on the string modelType
    // static Collision_Model* createCollisionModel(const std::string& modelType, User *user, const std::array<D_Phy_real,C_max_level+1> tau)
    static Collision_Model* createCollisionModel(const std::string& modelType)
    {
        if (modelType == "LBGK") {
            // return Collision_LBGK::getInstance(user, tau);
            return Collision_LBGK::getInstance();
        }
        else if (modelType == "Cumulant") {
            return Collision_Cumulant::getInstance();
        }
        else {
            std::string error_msg = "Invalid collision model type: " + modelType;
            throw std::invalid_argument(error_msg);
        }
        
    }

};