#include "header/mortonTrie.hpp"
#include "header/utils.hpp"
#include <iostream>

#define LOOP_METHOD_FINDFIRST1

void MortonTrie::insert(const D_morton morton, const uint8_t level, const double value) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            node->children = new TrieNode[BRANCH];
            node->children[index].level = i_depth;
        }
        
        node = &node->children[index];
    }
    // last level
        uint8_t start_pos = DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->data) {
            node->data = new double[BRANCH];
            #ifdef ABPattern
            node->data2 = new double[BRANCH];
            #endif
            // node->children[index].level = search_depth - 1;
        }
        node->data[index] = value;

    // Update existence bitmap
    set_existence_bit(morton, level);
}

double MortonTrie::search(const D_morton morton, const uint8_t level) const {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[search] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
    }

    return node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
}

void MortonTrie::init_searchStack(const D_morton morton,const uint8_t level) {
    TrieNode* node = root; search_path_stack->clear(); search_path_stack->resize(depth);

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        node = &node->children[EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos)];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;
}

double MortonTrie::search_halfway(const D_morton morton, const uint8_t level) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return search_path_stack->at(search_depth - 2)->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[Bottom-Up search] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;

    return node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
    
}


double MortonTrie::search_halfway_noUpdateStack(const D_morton morton, const uint8_t level) const {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return search_path_stack->at(search_depth - 2)->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[search_bottom_up] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
    }

    return node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)];
}

void MortonTrie::update_searchStack(const D_morton morton, const uint8_t level) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return;
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        node =&node->children[EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos)];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;
}

void MortonTrie::update(const D_morton morton, const uint8_t level, const double value) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[Update] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
    }

    #ifdef AAPattern
    node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif
    #ifdef ABPattern
    node->data2[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif

    // Update existence bitmap (in case this is a new cell)
    set_existence_bit(morton, level);
}

void MortonTrie::update_halfway(const D_morton morton, const uint8_t level, const double value) {
    TrieNode* node = root;

    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        #ifdef AAPattern
        search_path_stack->at(search_depth - 2)->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
        #endif
        #ifdef ABPattern
        search_path_stack->at(search_depth - 2)->data2[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
        #endif
    }

    node = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        if (!node->children) {
            std::cerr << "[Bottom-Up search] " << morton << " not found" << std::endl;
            exit(-1);
        }
        
        node =&node->children[index];
        search_path_stack->at(i_depth) = node;
    }

    stack_code = morton;

    #ifdef AAPattern
    node->data[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif
    #ifdef ABPattern
    node->data2[EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)] = value;
    #endif

}

// Bitmap-based O(1) existence query implementation
bool MortonTrie::exists(const D_morton& morton, const uint8_t level) const {
    uint32_t bg_index = extract_background_index(morton);
    uint32_t ref_index = extract_refinement_index(morton, level);

    if (bg_index >= BACKGROUND_GRID_SIZE || ref_index >= BITMAP_BITS_PER_CELL) {
        return false; // Out of bounds
    }

    return existence_bitmaps[bg_index][ref_index];
}

uint32_t MortonTrie::extract_background_index(const D_morton& morton) const {
    // For simplicity, use a basic approach: extract the higher bits and convert to index
    // Background bits start at bit position: refine_level * DIM
    uint32_t bg_bits_start = refine_level * DIM;

    // Extract the background portion as a simple integer
    uint64_t bg_portion = 0;
    for (uint32_t i = bg_bits_start; i < BIT; ++i) {
        if (morton[i]) {
            bg_portion |= (1ULL << (i - bg_bits_start));
        }
    }

    // For now, use modulo to ensure we stay within bounds
    // This is a simplified approach - in production, proper Morton decoding should be used
    return static_cast<uint32_t>(bg_portion % BACKGROUND_GRID_SIZE);
}

uint32_t MortonTrie::extract_refinement_index(const D_morton& morton, const uint8_t level) const {
    // Extract refinement field from lower bits
    uint32_t ref_index = 0;
    uint32_t bits_to_extract = level * DIM;

    // Extract the lower bits that represent the refinement path
    for (uint32_t i = 0; i < bits_to_extract; ++i) {
        if (morton[i]) {
            ref_index |= (1U << i);
        }
    }

    return ref_index;
}

void MortonTrie::set_existence_bit(const D_morton& morton, const uint8_t level) {
    uint32_t bg_index = extract_background_index(morton);
    uint32_t ref_index = extract_refinement_index(morton, level);

    if (bg_index < BACKGROUND_GRID_SIZE && ref_index < BITMAP_BITS_PER_CELL) {
        existence_bitmaps[bg_index][ref_index] = true;
    }
}

