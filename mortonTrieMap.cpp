#include "header/mortonTrie.hpp"
#include <stdexcept>
#include <algorithm>
#include <type_traits>

// Note: This is a simplified implementation that stores data in a cache
// A full implementation would need to traverse the trie structure to build iterators
// For now, we'll use a cache-based approach for compatibility

template<typename Value>
void MortonTrieMap<Value>::update_cache() const {
    if (cache_valid_) return;
    
    // For now, we'll maintain a simple cache
    // In a full implementation, this would traverse the trie structure
    cache_valid_ = true;
}

template<typename Value>
typename MortonTrieMap<Value>::iterator MortonTrieMap<Value>::begin() const {
    update_cache();
    return iterator(const_cast<MortonTrieMap<Value>*>(this),
                   const_cast<std::vector<std::pair<D_morton, Value>>&>(cached_data_).begin());
}

template<typename Value>
typename MortonTrieMap<Value>::iterator MortonTrieMap<Value>::end() const {
    update_cache();
    return iterator(const_cast<MortonTrieMap<Value>*>(this),
                   const_cast<std::vector<std::pair<D_morton, Value>>&>(cached_data_).end());
}

template<typename Value>
bool MortonTrieMap<Value>::empty() const {
    update_cache();
    return cached_data_.empty();
}

template<typename Value>
typename MortonTrieMap<Value>::size_type MortonTrieMap<Value>::size() const {
    update_cache();
    return cached_data_.size();
}

template<typename Value>
void MortonTrieMap<Value>::clear() {
    // Note: MortonTrie doesn't have a clear method, so we can't actually clear it
    // This is a limitation of the current MortonTrie implementation
    cached_data_.clear();
    invalidate_cache();
}

template<typename Value>
std::pair<typename MortonTrieMap<Value>::iterator, bool> 
MortonTrieMap<Value>::insert(const value_type& value) {
    return insert(std::make_pair(value.first, value.second));
}

template<typename Value>
std::pair<typename MortonTrieMap<Value>::iterator, bool> 
MortonTrieMap<Value>::insert(const std::pair<D_morton, Value>& value) {
    // Check if key already exists
    auto it = find(value.first);
    if (it != end()) {
        return std::make_pair(it, false);
    }
    
    // Insert into trie - we need to handle the conversion from Value to double
    // This is a limitation: MortonTrie only stores double values
    if constexpr (std::is_convertible_v<Value, double>) {
        trie_.insert(value.first, default_level_, static_cast<double>(value.second));
    } else {
        // For non-convertible types, we'll store in cache only
        // This is not ideal but maintains compatibility
    }
    
    // Add to cache
    cached_data_.emplace_back(value.first, value.second);
    invalidate_cache();
    
    return std::make_pair(find(value.first), true);
}

template<typename Value>
template<typename... Args>
std::pair<typename MortonTrieMap<Value>::iterator, bool> 
MortonTrieMap<Value>::emplace(Args&&... args) {
    value_type value(std::forward<Args>(args)...);
    return insert(value);
}

template<typename Value>
typename MortonTrieMap<Value>::iterator MortonTrieMap<Value>::find(const D_morton& key) const {
    update_cache();
    auto& mutable_cache = const_cast<std::vector<std::pair<D_morton, Value>>&>(cached_data_);
    auto it = std::find_if(mutable_cache.begin(), mutable_cache.end(),
                          [&key](const auto& pair) { return pair.first == key; });
    return iterator(const_cast<MortonTrieMap<Value>*>(this), it);
}

template<typename Value>
typename MortonTrieMap<Value>::size_type MortonTrieMap<Value>::count(const D_morton& key) const {
    return find(key) != end() ? 1 : 0;
}

template<typename Value>
Value& MortonTrieMap<Value>::operator[](const D_morton& key) {
    // Find existing entry
    auto cache_it = std::find_if(cached_data_.begin(), cached_data_.end(),
                                [&key](const auto& pair) { return pair.first == key; });

    if (cache_it != cached_data_.end()) {
        return cache_it->second;
    }

    // Insert default value into both cache and trie
    Value default_value{};
    cached_data_.emplace_back(key, default_value);

    // Also insert into the trie to maintain bitmap consistency
    if constexpr (std::is_convertible_v<Value, double>) {
        trie_.insert(key, default_level_, static_cast<double>(default_value));
    } else {
        // For non-convertible types, insert a placeholder value
        trie_.insert(key, default_level_, 0.0);
    }

    invalidate_cache();
    return cached_data_.back().second;
}

template<typename Value>
const Value& MortonTrieMap<Value>::at(const D_morton& key) const {
    auto it = find(key);
    if (it == end()) {
        throw std::out_of_range("MortonTrieMap::at: key not found");
    }
    return it->second;
}

template<typename Value>
template<typename InputIt>
void MortonTrieMap<Value>::insert(InputIt first, InputIt last) {
    for (auto it = first; it != last; ++it) {
        insert(*it);
    }
}

// Explicit template instantiations for commonly used types
template class MortonTrieMap<int>;
template class MortonTrieMap<double>;
template class MortonTrieMap<float>;
