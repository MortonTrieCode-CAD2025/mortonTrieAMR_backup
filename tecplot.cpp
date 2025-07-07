/**
 * @file tecplot.cpp
 * @date 2023-09-21
 */

#include "tecplot.h"
#include "Lat_Manager.hpp"
#include "LBM_Manager.hpp"
#include "Morton_Assist.h"
#include <iostream>
#include <sstream>
#include <map>

void IO_TECPLOTH::write(const D_int iter)
{
    char fileName[1024];
    sprintf(fileName, "sphereFlow_%f.dat", iter * LBM_Manager::pointer_me->CF_T[0]);
    
    D_Phy_real Convert_U = LBM_Manager::pointer_me->CF_U;

    // calculate the node numbers
    int elemNums = 0;
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        elemNums += Lat_Manager::pointer_me->lat_f.at(i_level).size() + Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level).size();
    }
    
    D_map_define<std::array<int,2>> node_idx;
    std::vector<D_morton> mortonCode_list;
    int node_idx_count = 0;
    for (auto i_level = 0; i_level <= C_max_level; ++i_level)
    {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
        {
            std::array<D_morton, 8> vertex ={
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),
                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level),i_level)
            };
            for (auto i_vet : vertex) {
                if (node_idx.find(i_vet) == node_idx.end()) {
                    ++node_idx_count;
                    node_idx.insert(make_pair(i_vet, std::array<int,2>{node_idx_count, i_level}));
                    mortonCode_list.push_back(i_vet);
                }
            }
        }
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
        {
            std::array<D_morton, 8> vertex = {
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),
                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level), i_level)
            };
            for (auto i_vet : vertex) {
                if (node_idx.find(i_vet) == node_idx.end()) {
                    ++node_idx_count;
                    node_idx.insert(make_pair(i_vet, std::array<int,2>{node_idx_count, i_level}));
                    mortonCode_list.push_back(i_vet);
                }
            }
        }
    }
    int nodeNums = node_idx.size();
    std::ofstream plt_file(fileName);
    if (!plt_file.is_open()) {
        std::cerr << "Failed to open .plt file for writing." << std::endl;
        exit(0);
    }

    plt_file << "TITLE = \"Flow Field Example\"" << std::endl;
    plt_file << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\"" << std::endl;
    plt_file << "ZONE, NODES=" << nodeNums << ", ELEMENTS=" << elemNums << ", DATAPACKING=BLOCK, VARLOCATION=([4,5,6]=CELLCENTERED), ZONETYPE=FEBRICK" << std::endl;

    for (auto i_code : mortonCode_list)
    {
        double x,y,z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        plt_file << x << std::endl;
    }
    for (auto i_code : mortonCode_list)
    {
        double x,y,z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        plt_file << y << std::endl;
    }
    for (auto i_code : mortonCode_list)
    {
        double x,y,z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        plt_file << z << std::endl;
    }

    // u
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).x * Convert_U << std::endl;
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).x * Convert_U << std::endl;
    }

    // v
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).y * Convert_U << std::endl;
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).y * Convert_U << std::endl;
    }

    // w
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).z * Convert_U << std::endl;
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).z * Convert_U << std::endl;
    }

    // cell information
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level)) 
        {
            std::array<D_morton, 8> vertex = {
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),

                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level)
                
            };
            for (auto i_vex : vertex) {
                plt_file << node_idx.at(i_vex)[0] << " ";
            }
            plt_file << std::endl;
        }
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level)) 
        {
            std::array<D_morton, 8> vertex ={
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),

                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level)
                
            };
            for (auto i_vex : vertex) {
                plt_file << node_idx.at(i_vex)[0] << " ";
            }
            plt_file << std::endl;
        }
    }

    plt_file.close();
}