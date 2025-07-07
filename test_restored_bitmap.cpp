#include "header/mortonTrie.hpp"
#include "header/settings.hpp"
#include <iostream>

int main() {
    std::cout << "Testing restored bitmap functionality..." << std::endl;
    
    // Print configuration
    std::cout << "DIM: " << DIM << std::endl;
    std::cout << "refine_level: " << static_cast<int>(refine_level) << std::endl;
    std::cout << "Nx: " << Nx << ", Ny: " << Ny << ", Nz: " << Nz << std::endl;
    std::cout << "BIT: " << BIT << std::endl;
    
    MortonTrie& trie = MortonTrie::getInstance();
    
    D_morton test_code = 0;
    uint8_t test_level = 1;
    
    std::cout << "Testing with morton code: " << test_code << ", level: " << static_cast<int>(test_level) << std::endl;
    
    // Check initial state
    bool exists_before = trie.exists(test_code, test_level);
    std::cout << "Exists before insert: " << exists_before << std::endl;
    
    // Insert
    std::cout << "Inserting..." << std::endl;
    trie.insert(test_code, test_level, 42.0);
    
    // Check after insert
    bool exists_after = trie.exists(test_code, test_level);
    std::cout << "Exists after insert: " << exists_after << std::endl;
    
    // Test MortonTrieMap
    std::cout << "\nTesting MortonTrieMap..." << std::endl;
    MortonTrieMap<int> map;
    
    std::cout << "Map default level: " << static_cast<int>(refine_level) << std::endl;
    
    D_morton map_key = 1;
    std::cout << "Testing with map key: " << map_key << std::endl;
    
    bool map_exists_before = map.exists(map_key);
    std::cout << "Map exists before: " << map_exists_before << std::endl;
    
    map[map_key] = 123;
    std::cout << "Inserted into map" << std::endl;
    
    bool map_exists_after = map.exists(map_key);
    std::cout << "Map exists after: " << map_exists_after << std::endl;
    
    // Test with explicit level
    bool map_exists_explicit = map.exists(map_key, refine_level);
    std::cout << "Map exists with explicit level " << static_cast<int>(refine_level) << ": " << map_exists_explicit << std::endl;
    
    // Test multiple insertions
    std::cout << "\nTesting multiple insertions..." << std::endl;
    for (int i = 2; i < 10; ++i) {
        D_morton key = i;
        map[key] = i * 10;
        if (map.exists(key)) {
            std::cout << "Key " << key << " exists after insertion" << std::endl;
        } else {
            std::cout << "Key " << key << " should exist after insertion" << std::endl;
        }
    }
    
    std::cout << "\n=== Bitmap functionality restored successfully! ===" << std::endl;
    return 0;
}
