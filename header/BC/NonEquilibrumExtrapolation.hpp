/**
 * @file nonEquilibrumExtrapolation.h
 * @brief Non-equilibrum Extrapolation Method
 * @date 2023-07-27
 */

#pragma once
#include "user.h"
#include "D3Q19_BC_Manager.hpp"
#define eflow_stable_bc // A kind of stable nonequilibrum extrapolation bc implemented in eflow. 
                        //  ifndef "eflow_stable_bc",means adopting the original nonequilibrum extrapolation theory proposed by GuoZhaoli

class nonEquilibrumExtrapolation_BC : public D3Q19_BC_Strategy
{
public:
    virtual void applyBCStrategy() = 0;
    virtual void initialBCStrategy() = 0;
};

class nonEquilibrumExtrapolation_West : public nonEquilibrumExtrapolation_BC
{
private:
    D_uvw bc_velocity;
    static constexpr std::array<uint,9> lat_bc_idx 
    // std::array<uint,9> lat_bc_idx 
        = {Bdry_Type::BDRY_FACE_WEST,
            Bdry_Type::BDRY_EDGE_WN,
             Bdry_Type::BDRY_EDGE_WS,
              Bdry_Type::BDRY_EDGE_WT,
               Bdry_Type::BDRY_EDGE_WB,
                Bdry_Type::BDRY_CORNER_WNT,
                 Bdry_Type::BDRY_CORNER_WST,
                  Bdry_Type::BDRY_CORNER_WNB,
                   Bdry_Type::BDRY_CORNER_WSB};

public:
    nonEquilibrumExtrapolation_West(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumExtrapolation_East : public nonEquilibrumExtrapolation_BC
{
private:
    D_uvw bc_velocity;
    static constexpr std::array<uint,9> lat_bc_idx
        = {Bdry_Type::BDRY_FACE_EAST,
            Bdry_Type::BDRY_EDGE_EN,
             Bdry_Type::BDRY_EDGE_ES,
              Bdry_Type::BDRY_EDGE_ET,
               Bdry_Type::BDRY_EDGE_EB,
                Bdry_Type::BDRY_CORNER_ENT,
                 Bdry_Type::BDRY_CORNER_EST,
                  Bdry_Type::BDRY_CORNER_ENB,
                   Bdry_Type::BDRY_CORNER_ESB};

public:
    nonEquilibrumExtrapolation_East(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumExtrapolation_North : public nonEquilibrumExtrapolation_BC
{
private:
    D_uvw bc_velocity;
    static constexpr std::array<uint,3> lat_bc_idx
        = {Bdry_Type::BDRY_FACE_NORTH,
            Bdry_Type::BDRY_EDGE_NT,
             Bdry_Type::BDRY_EDGE_NB};

public:
    nonEquilibrumExtrapolation_North(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumExtrapolation_South : public nonEquilibrumExtrapolation_BC
{
private:
    D_uvw bc_velocity;
    static constexpr std::array<uint,3> lat_bc_idx
        = {Bdry_Type::BDRY_FACE_SOUTH,
            Bdry_Type::BDRY_EDGE_ST,
             Bdry_Type::BDRY_EDGE_SB};

public:
    nonEquilibrumExtrapolation_South(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumExtrapolation_Top : public nonEquilibrumExtrapolation_BC
{
private:
    D_uvw bc_velocity;
    static constexpr std::array<uint,1> lat_bc_idx = {Bdry_Type::BDRY_FACE_TOP};

public:
    nonEquilibrumExtrapolation_Top(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class nonEquilibrumExtrapolation_Bot : public nonEquilibrumExtrapolation_BC
{
private:
    D_uvw bc_velocity;
    static constexpr std::array<uint,1> lat_bc_idx = {Bdry_Type::BDRY_FACE_BOT};

public:
    nonEquilibrumExtrapolation_Bot(D_uvw input_v) : bc_velocity(input_v/CF_U_) {}
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};