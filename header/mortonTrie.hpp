#pragma once

#include "settings.hpp"
#include "General.h"
#include <array>
#include <vector>
#include <utility>
#include <iterator>
#include <cmath>
#include <bitset>

class TrieNode {
public:
    TrieNode* children; //in Z-order: 000, 001, 010, 011, 100, 101, 110, 111
    uint8_t level;
    double* data;
    #ifdef ABPattern
    double* data2;
    #endif
    
    TrieNode() : children(nullptr), level(0) {}
};

using NgbrList = std::array<TrieNode*, STENCIL_SIZE-1>;

class MortonTrie {
public:

    void insert(const D_morton morton, const uint8_t level, const double value);
    double search(const D_morton morton, const uint8_t level) const; //Top-down search
    double search_halfway(const D_morton morton, const uint8_t level); //Bottom-up search
    double search_halfway_noUpdateStack(const D_morton morton, const uint8_t level) const; //Bottom-up search without updating the search path stack
    void init_searchStack(D_morton morton, uint8_t level);
    void update(const D_morton morton, const uint8_t level, const double value);
    void update_halfway(const D_morton morton, const uint8_t level, const double value);
    void update_searchStack(const D_morton morton, const uint8_t level);

    // Bitmap-based O(1) existence query
    bool exists(const D_morton& morton, const uint8_t level) const;
    
    uint8_t get_depth() const { return depth; }
    TrieNode* get_root() const { return root; }
    
    static MortonTrie& getInstance() {
        static MortonTrie instance;
        return instance;
    };


private:

    TrieNode* root;
    uint8_t depth;

    // Bitmap-based existence query data structures
    static constexpr uint32_t BITMAP_BITS_PER_CELL = 1 << (refine_level * DIM); // 2^(L*DIM) bits per background cell
    static constexpr uint32_t BACKGROUND_GRID_SIZE = Nx * Ny * Nz; // Total background grid cells
    using CellBitmap = std::bitset<BITMAP_BITS_PER_CELL>;
    std::vector<CellBitmap> existence_bitmaps; // One bitmap per background grid cell

    MortonTrie() {
        depth = static_cast<uint8_t>(log2(max_of_three(Nx, Ny, Nz)) + 1 + refine_level);
        // not containing the single root level , +1 for std::ceil(log2(max_of_three(Nx, Ny, Nz)))
        root = new TrieNode();
        search_path_stack = new std::vector<TrieNode*>(depth);

        // Initialize existence bitmaps
        existence_bitmaps.resize(BACKGROUND_GRID_SIZE);
    }
    ~MortonTrie() = default;

    // Helper methods for bitmap operations
    uint32_t extract_background_index(const D_morton& morton) const;
    uint32_t extract_refinement_index(const D_morton& morton, const uint8_t level) const;
    void set_existence_bit(const D_morton& morton, const uint8_t level);

    std::vector<TrieNode*> *search_path_stack;
    D_morton stack_code = 0;
};

// Forward declaration for MortonTrieMap
template<typename Value>
class MortonTrieMap;

// Iterator class for MortonTrieMap
template<typename Value>
class MortonTrieIterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::pair<D_morton, Value>;  // Remove const from key for simplicity
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

private:
    MortonTrieMap<Value>* map_;
    typename std::vector<std::pair<D_morton, Value>>::iterator iter_;

public:
    MortonTrieIterator(MortonTrieMap<Value>* map,
                      typename std::vector<std::pair<D_morton, Value>>::iterator iter)
        : map_(map), iter_(iter) {}

    // Const constructor for const iterators
    MortonTrieIterator(const MortonTrieMap<Value>* map,
                      typename std::vector<std::pair<D_morton, Value>>::const_iterator iter)
        : map_(const_cast<MortonTrieMap<Value>*>(map)),
          iter_(const_cast<typename std::vector<std::pair<D_morton, Value>>::iterator>(iter)) {}

    reference operator*() const { return *iter_; }
    pointer operator->() const { return &(*iter_); }

    MortonTrieIterator& operator++() { ++iter_; return *this; }
    MortonTrieIterator operator++(int) { auto tmp = *this; ++iter_; return tmp; }

    bool operator==(const MortonTrieIterator& other) const { return iter_ == other.iter_; }
    bool operator!=(const MortonTrieIterator& other) const { return iter_ != other.iter_; }
};

/**
 * @brief STL-compatible wrapper for MortonTrie that provides map-like interface
 * This class allows MortonTrie to be used as a drop-in replacement for std::unordered_map and std::map
 */
template<typename Value>
class MortonTrieMap {
public:
    using key_type = D_morton;
    using mapped_type = Value;
    using value_type = std::pair<const D_morton, Value>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using iterator = MortonTrieIterator<Value>;
    using const_iterator = MortonTrieIterator<Value>;

private:
    MortonTrie& trie_;
    uint8_t default_level_;
    std::vector<std::pair<D_morton, Value>> cached_data_;
    mutable bool cache_valid_;

    void invalidate_cache() { cache_valid_ = false; }
    void update_cache() const;

public:
    explicit MortonTrieMap(uint8_t level = refine_level)
        : trie_(MortonTrie::getInstance()), default_level_(level), cache_valid_(false) {}

    // STL container interface
    iterator begin() const;
    iterator end() const;
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

    bool empty() const;
    size_type size() const;
    void clear();

    // Map interface
    std::pair<iterator, bool> insert(const value_type& value);
    std::pair<iterator, bool> insert(const std::pair<D_morton, Value>& value);
    template<typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args);

    iterator find(const D_morton& key) const;
    size_type count(const D_morton& key) const;

    Value& operator[](const D_morton& key);
    const Value& at(const D_morton& key) const;

    // Bitmap-based O(1) existence query
    bool exists(const D_morton& key) const {
        return trie_.exists(key, default_level_);
    }

    bool exists(const D_morton& key, uint8_t level) const {
        return trie_.exists(key, level);
    }

    // Range insertion for compatibility with existing code
    template<typename InputIt>
    void insert(InputIt first, InputIt last);
};