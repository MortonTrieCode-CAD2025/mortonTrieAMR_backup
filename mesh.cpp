#include "header/mortonTrie.hpp"
#include "header/utils.hpp"
#include "header/mesh.hpp"
#include "header/Morton_Assist.hpp"
#include <cmath>
#include <iostream>

/**
 * @brief Generate a tree to simulate the AMR inside a cube
*/
void Mesh::initial()
{
    generate_inner_gird();

    // for (uint i_level = max_level-2; i_level > 0; --i_level) {
    //     generate_intermediate_grid(i_level);
    // }

    // generate_outer_grid();

}

void Mesh::generate_inner_gird()
{
    Timer tmr;
    
    // computer the domain of inner grid
    domain[refine_level][0] = Nx / 4;
    domain[refine_level][1] = Ny / 4 ;
    domain[refine_level][2] = Nz / 4;
    domain[refine_level][3] = Nx - domain[refine_level][0];
    domain[refine_level][4] = Ny - domain[refine_level][1];
    domain[refine_level][5] = Nz - domain[refine_level][2];

    const uint x_start = static_cast<uint>(std::floor(domain[refine_level][0] / ddx));
    const uint x_end   = static_cast<uint>(std::ceil( domain[refine_level][3] / ddx));
    const uint y_start = static_cast<uint>(std::floor(domain[refine_level][1] / ddx));
    const uint y_end   = static_cast<uint>(std::ceil( domain[refine_level][4] / ddx));
    const uint z_start = static_cast<uint>(std::floor(domain[refine_level][2] / ddx));
    const uint z_end   = static_cast<uint>(std::ceil( domain[refine_level][5] / ddx));

    tmr.reset();
    //insert morton code to the trie
    for (uint ix = x_start; ix < x_end; ++ix) {
        for (uint iy = y_start; iy < y_end; ++iy) {
            for (uint iz = z_start; iz < z_end; ++iz) {
                amr_trie.insert(Morton_Assist::morton_encode(P_Coord(ix*ddx, iy*ddx, iz*ddx)), refine_level, ix+iy+iz);
            }
        }
    } 
    std::cout << "Trie Insert: " << tmr.elapsed() << " s." << std::endl;

    
    #if defined(USE_HASHMAP) || defined(USE_RBTREE)
    tmr.reset();
    for (uint ix = x_start; ix < x_end; ++ix) {
        for (uint iy = y_start; iy < y_end; ++iy) {
            for (uint iz = z_start; iz < z_end; ++iz) {
                amr_hashmap.insert(make_pair(Morton_Assist::morton_encode(P_Coord(ix*ddx, iy*ddx, iz*ddx)), ix+iy+iz));
            }
        }
    } 
    #endif

    #ifdef USE_HASHMAP
    std::cout << "Hash Insert: " << tmr.elapsed() << " s." << std::endl;
    #endif
    #ifdef USE_RBTREE
    std::cout << "RB-Tree Insert: " << tmr.elapsed() << " s." << std::endl;
    #endif

}

void Mesh::generate_intermediate_grid(uint level)
{
    dx[level] = 1. / two_pow(level);
    // computer the domain of interlayer grid
    // domain[level+1][0] -= dx[level] * extend_width;
    // domain[level+1][1] -= dx[level] * extend_width;
    // domain[level+1][2] -= dx[level] * extend_width;
    // domain[level+1][3] += dx[level] * extend_width;
    // domain[level+1][4] += dx[level] * extend_width;
    // domain[level+1][5] += dx[level] * extend_width;

    
}

void Mesh::generate_outer_grid()
{
    dx[0] = 1.;
    domain[0][0] = 0;
    domain[0][1] = 0;
    domain[0][2] = 0;
    domain[0][3] = Nx;
    domain[0][4] = Ny;
    domain[0][5] = Nz;

}