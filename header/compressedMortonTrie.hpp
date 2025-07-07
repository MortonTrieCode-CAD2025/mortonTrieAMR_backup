#pragma once

#include "settings.hpp"
#include "mortonTrie.hpp"

class CompressedMortonTrie {
public:

    void init_compressed_Mtrie();

    double search(const D_morton morton, const uint8_t level) const; //Top-down search
    double search_halfway(const D_morton morton, const uint8_t level); //Bottom-up search
    double search_halfway_noUpdateStack(const D_morton morton, const uint8_t level) const;
    void init_searchStack(const D_morton morton, const uint8_t level);
    void update_searchStack(const D_morton morton, const uint8_t level);
    // void outside_searchStack(const D_morton input_morton, const uint8_t level, D_morton& outside_morton, std::vector<ulong>& search_path_stack);
    void update(const D_morton morton, const uint8_t level, const double value);
    void update_halfway(const D_morton morton, const uint8_t level, const double value);

    static CompressedMortonTrie& getInstance(const MortonTrie& org_trie) {
        static CompressedMortonTrie instance(org_trie);
        return instance;
    };

    

private:
    CompressedMortonTrie(const MortonTrie& org_trie) : original_Mtrie(org_trie), depth(org_trie.get_depth()) {
        node_index = new std::vector<std::vector<long>>(depth);
        search_path_stack = new std::vector<ulong>(depth);
    };

    const MortonTrie& original_Mtrie;
    const uint8_t depth;
    std::vector<double> data;
    #ifdef ABPattern
    std::vector<double> data2;
    #endif
    std::vector<std::vector<long>> *node_index;

    std::vector<ulong> *search_path_stack;
    D_morton stack_code = 0;

    void bfs_traversal();
    
};