/**
 * @file ZouHe_Velocity.h
 * @brief Non-equilibrum Bounce-Back (ZouHe) Velocity BC = Dirichlet BC = Flux BC (suitable for Low Re Number)
 * @date 2023-07-18
 */
#define BC_NO_GHOST
#ifndef BC_NO_GHOST
#pragma once
#include "user.h"
#include "D3Q19_BC_Manager.hpp"
#include "LBM_Manager.hpp"

class nonEquilibrumBounceBack_Velocity_BC : public D3Q19_BC_Strategy
{
// protected:
//     D_Phy_real CF_U;

public:
    // nonEquilibrumBounceBack_Velocity_BC() : CF_U(LBM_Manager::pointer_me->CF_U) { };
    virtual void applyBCStrategy() = 0;
    virtual void initialBCStrategy() = 0;
};

class nonEquilibrumBounceBack_Velocity_West : public nonEquilibrumBounceBack_Velocity_BC
{
private:
    D_uvw bc_velocity;

public:
    nonEquilibrumBounceBack_Velocity_West(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumBounceBack_Velocity_East : public nonEquilibrumBounceBack_Velocity_BC
{
private:
    D_uvw bc_velocity;

public:
    nonEquilibrumBounceBack_Velocity_East(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumBounceBack_Velocity_North : public nonEquilibrumBounceBack_Velocity_BC
{
private:
    D_uvw bc_velocity;

public:
    nonEquilibrumBounceBack_Velocity_North(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumBounceBack_Velocity_South : public nonEquilibrumBounceBack_Velocity_BC
{
private:
    D_uvw bc_velocity;

public:
    nonEquilibrumBounceBack_Velocity_South(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumBounceBack_Velocity_Top : public nonEquilibrumBounceBack_Velocity_BC
{
private:
    D_uvw bc_velocity;

public:
    nonEquilibrumBounceBack_Velocity_Top(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumBounceBack_Velocity_Bot : public nonEquilibrumBounceBack_Velocity_BC
{
private:
    D_uvw bc_velocity;

public:
    nonEquilibrumBounceBack_Velocity_Bot(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};
#endif