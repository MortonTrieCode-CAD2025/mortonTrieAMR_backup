#pragma once
#include <bitset>
#include <cstdint>

#define DIM 3
#define STENCIL_SIZE 27
#if (DIM == 3)
    constexpr uint8_t BRANCH = 8;
#elif (DIM == 2)
    constexpr uint8_t BRANCH = 4;
#endif

constexpr uint Nx = 400/20, Ny = 300/20, Nz = 100/20;
constexpr uint8_t refine_level = 4;
constexpr uint8_t level_num = refine_level + 1;

constexpr uint8_t extend_width = 6;

constexpr uint8_t BIT = 64;  ///< Bits to store morton code 
using D_morton = std::bitset<BIT>;  ///> Bitset to store morton code, due to the method used to find neighbours in Morton_Assist.h, the bit is limited to unsigned long long
using D_3BIT = std::bitset<DIM>;

// #define USE_HASHMAP
// #define USE_RBTREE
// #define USE_MORTONTRIE
// #define RANDOM_SEARCH_TEST

// #define TRIE_DEBUG
// #define AAPattern
#define ABPattern