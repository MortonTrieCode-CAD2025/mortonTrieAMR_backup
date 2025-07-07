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
    const char* meshFileName_prefix = "Mesh";
    std::vector<unsigned int> grid_vlevel;
    std::vector<unsigned int> grid_npoints, grid_nelements;

public:
    void write(const D_int iter);

    // Mesh output functions
    int write_mesh_manager(std::string outfile);
    void write_mesh(std::string outfile);
    template <class T_grptr> void write_mesh_grid(T_grptr &grid_ptr, int level_index, std::ofstream &plt_file, int &zone_count);
    void write_mesh_boundaries(std::ofstream &plt_file, int &zone_count);
    void write_mesh_solid_boundaries(std::ofstream &plt_file, int &zone_count);
};