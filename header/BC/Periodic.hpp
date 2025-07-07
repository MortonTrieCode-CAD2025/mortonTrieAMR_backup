/**
 * @file Periodic.h
 * @date 2023-09-22
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#define BC_NO_GHOST
#ifndef BC_NO_GHOST
#pragma once
#include "user.h"
#include "D3Q19_BC_Manager.hpp"
#include <unordered_set>

class periodic_BC : public D3Q19_BC_Strategy
{
public:
    virtual void applyBCStrategy() = 0;
    virtual void initialBCStrategy() = 0;
};

class periodic_West : public periodic_BC
{
private:
    D_map_define<D_morton> ghost_bc_x0;

public:
    periodic_West();
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class periodic_East : public periodic_BC
{
private:
    D_map_define<D_morton> ghost_bc_x1;
        

public:
    periodic_East();
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class periodic_South : public periodic_BC
{
private:
    D_map_define<D_morton> ghost_bc_y0;

public:
    periodic_South();
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class periodic_North : public periodic_BC
{
private:
    D_map_define<D_morton> ghost_bc_y1;

public:
    periodic_North();
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class periodic_Top : public periodic_BC
{
private:
    D_map_define<D_morton> ghost_bc_z1;

public:
    periodic_Top();
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};

class periodic_Bot : public periodic_BC
{
private:
    D_map_define<D_morton> ghost_bc_z0;
    
public:
    periodic_Bot();
    void applyBCStrategy() override;
    void initialBCStrategy() override;
};
#endif