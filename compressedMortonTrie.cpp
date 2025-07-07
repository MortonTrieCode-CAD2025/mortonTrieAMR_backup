#include "header/compressedMortonTrie.hpp"
#include <queue>

void CompressedMortonTrie::init_compressed_Mtrie()
{
    bfs_traversal();

    #ifdef ABPattern
    data2.insert(data2.begin(), data.begin(), data.end());
    #endif
}
void CompressedMortonTrie::bfs_traversal()
{
    std::queue<TrieNode*> ptr_queue;

    long index_addition = 0;
    node_index->at(0).insert(node_index->at(0).end(), BRANCH, -1);

    for (int j = 0; j < BRANCH; ++j) {
        if (original_Mtrie.get_root()->children[j].children) {
            node_index->at(0).at(j) = (index_addition++) * BRANCH;
            ptr_queue.push(&original_Mtrie.get_root()->children[j]);
            node_index->at(1).insert(node_index->at(1).end(), BRANCH, -1);
        }
    }
    
    while (!ptr_queue.empty())
    {
        long level_node_size = ptr_queue.size();
        uint8_t level = ptr_queue.front()->level;
        index_addition = 0;

        for (long i = 0; i < level_node_size; ++i) {
            TrieNode* node = ptr_queue.front();
            ptr_queue.pop();

            if (!node->data) {
                for (uint8_t j = 0; j < BRANCH; ++j) {
                    if (node->children[j].children || node->children[j].data) {
                        ptr_queue.push(&node->children[j]);
                        node_index->at(level+1).at(i * BRANCH + j) = (index_addition++) * BRANCH;
                        node_index->at(level+2).insert(node_index->at(level+2).end(), BRANCH, -1);
                    }
                } 
            } else {
                for (uint8_t j = 0; j < BRANCH; ++j) {
                    data.emplace_back(node->data[j]);
                }
            }
            
        }

    }
}

double CompressedMortonTrie::search(const D_morton morton, const uint8_t level) const {
    uint8_t search_depth = depth - refine_level + level;
    long ptr_index = 0;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        ptr_index = node_index->at(i_depth).at(ptr_index + index);
        if (ptr_index == -1) {
            std::cerr << "[CMTrie search] " << morton << " not found" << std::endl;
            exit(-1);
        }
    }

    return data.at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM));
}

double CompressedMortonTrie::search_halfway(const D_morton morton, const uint8_t level) {
    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return data.at(search_path_stack->at(search_depth - 2) + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM));
    }

    long ptr_index = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        ptr_index = node_index->at(i_depth).at(ptr_index + index);
        if (ptr_index == -1) {
            std::cerr << "[CMTrie Bottom_Up search] " << morton << " not found" << std::endl;
            exit(-1);
        }

        search_path_stack->at(i_depth) = ptr_index;
    }

    stack_code = morton;

    return data.at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM));
}

double CompressedMortonTrie::search_halfway_noUpdateStack(const D_morton morton, const uint8_t level) const {
    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return data.at(search_path_stack->at(search_depth - 2) + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM));
    }

    long ptr_index = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        ptr_index = node_index->at(i_depth).at(ptr_index + index);
        if (ptr_index == -1) {
            std::cerr << "[CMTrie Bottom_Up search] " << morton << " not found" << std::endl;
            exit(-1);
        }

    }


    return data.at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM));
}

void CompressedMortonTrie::update_searchStack(const D_morton morton, const uint8_t level) {

    uint8_t search_depth = depth - refine_level + level;
    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        return;
    }

    ulong ptr_index = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        ptr_index = node_index->at(i_depth).at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos));
        search_path_stack->at(i_depth) = ptr_index;
    }

    stack_code = morton;
}

// void CompressedMortonTrie::outside_searchStack(const D_morton morton, const uint8_t level) {

//     uint8_t search_depth = depth - refine_level + level;
//     D_morton xor_diff = stack_code ^ morton;
//     uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

//     if (start_search_depth == -1) {
//         return;
//     }

//     ulong ptr_index = search_path_stack->at(start_search_depth);

//     for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
//         uint8_t start_pos = (depth - i_depth) * DIM;
//         ptr_index = node_index->at(i_depth).at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos));
//         search_path_stack->at(i_depth) = ptr_index;
//     }

//     stack_code = morton;
// }

void CompressedMortonTrie::init_searchStack(const D_morton morton, const uint8_t level)
{
    search_path_stack->clear(); search_path_stack->resize(depth);

    uint8_t search_depth = depth - refine_level + level;
    long ptr_index = 0;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        ptr_index = node_index->at(i_depth).at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos));
        search_path_stack->at(i_depth) = ptr_index;
    }

    stack_code = morton;

}

void CompressedMortonTrie::update(const D_morton morton, const uint8_t level, const double value)
{
    uint8_t search_depth = depth - refine_level + level;
    long ptr_index = 0;
    for (uint8_t i_depth = 0; i_depth < search_depth - 1; ++i_depth)
    {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        ptr_index = node_index->at(i_depth).at(ptr_index + index);
        if (ptr_index == -1) {
            std::cerr << "[Update] " << morton << " not found" << std::endl;
            exit(-1);
        }
    }

    #ifdef AAPattern
    data.at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)) = value;
    #endif
    #ifdef ABPattern
    data2.at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)) = value;
    #endif
}

void CompressedMortonTrie::update_halfway(const D_morton morton, const uint8_t level, const double value) {
    uint8_t search_depth = depth - refine_level + level;

    D_morton xor_diff = stack_code ^ morton;
    uint start_search_depth = (__builtin_clzl(xor_diff.to_ulong()) - (BIT - DIM * depth)) / DIM - 1;

    if (start_search_depth == -1) {
        #ifdef AAPattern
        data.at(search_path_stack->at(search_depth - 2) + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)) = value;
        #endif
        #ifdef ABPattern
        data2.at(search_path_stack->at(search_depth - 2) + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)) = value;
        #endif
    }

    long ptr_index = search_path_stack->at(start_search_depth);

    for (uint8_t i_depth = start_search_depth + 1; i_depth < search_depth - 1; ++i_depth) {
        uint8_t start_pos = (depth - i_depth) * DIM;
        uint8_t index = EXTRACT_MORTON_CONVERT_INDEX(morton, start_pos);
        ptr_index = node_index->at(i_depth).at(ptr_index + index);
        if (ptr_index == -1) {
            std::cerr << "[CMTrie Bottom_Up search] " << morton << " not found" << std::endl;
            exit(-1);
        }

    }

    #ifdef AAPattern
    data.at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)) = value;
    #endif
    #ifdef ABPattern
    data2.at(ptr_index + EXTRACT_MORTON_CONVERT_INDEX(morton, DIM)) = value;
    #endif
}