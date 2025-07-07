/**
 * @file tecplot.h
 * @brief 
 * @date 2023-09-21
 */

#pragma once

#include "user.h"

class IO_TECPLOTH
{
    friend class IO_Manager;
private:
    const char* pltFileName_prefix = "SphereFlow";
public:
    void write(const D_int iter);
};