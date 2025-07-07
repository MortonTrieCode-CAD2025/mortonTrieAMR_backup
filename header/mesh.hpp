#pragma once

#include "mortonTrie.hpp"
#include "mortonHash.hpp"

class Mesh
{
private:
    MortonTrie& amr_trie;
    D_mapInt& amr_hashmap;

    std::array<double, level_num> dx;
    std::array<std::array<double, 2*DIM>, level_num> domain; // {x0, y0, z0, x1, y1, z1}
    // std::array<std::array<uint, 2*DIM>, level_num> int_domain; // {0, 0, 0, Nx, Ny, Nz}

    void generate_inner_gird();
    void generate_intermediate_grid(uint level);
    void generate_outer_grid();

public:
    Mesh(MortonTrie& trie, D_mapInt& hashmap) : amr_trie(trie), amr_hashmap(hashmap) {}

    void initial();
    std::array<double, 2*DIM> get_domain(uint level) const { return domain[level]; }
};