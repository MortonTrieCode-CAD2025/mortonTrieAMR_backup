# MortonTrie AMR

## Overview

An AMR framework for Lattice Boltzmann Method (LBM) simulations driven by MortonTrie.


## Switch Data Structure


1. Edit `header/Constants.h`
2. Select `#define C_MAP_TYPE 1/2/3`: 1 = Hash table, 2 = Red-black tree 3 = MortonTrie
3. Recompile our project


```
// This code works with all three configurations
D_map_define<double> velocity_map;
D_map_define<int> index_map;

// Common patterns that are supported
velocity_map.insert(std::make_pair(morton_code, velocity_value));
auto it = velocity_map.find(morton_code);
velocity_map[morton_code] = new_velocity;

// Range insertion (used in LBM_Manager.cpp)
velocity_map.insert(other_map.begin(), other_map.end());
```
