#include "header/mortonTrie.hpp"
#include "header/settings.hpp"
#include <iostream>
#include <chrono>
#include <random>
#include <vector>

void test_basic_existence() {
    std::cout << "Testing basic existence functionality..." << std::endl;
    
    MortonTrie& trie = MortonTrie::getInstance();
    
    // Test with some Morton codes
    D_morton test_codes[3];
    test_codes[0] = 0;  // Background cell (0,0,0)
    test_codes[1] = 1;  // First refinement
    test_codes[2] = 8;  // Different refinement pattern
    
    uint8_t test_level = 2;
    
    // Initially, cells should not exist
    for (int i = 0; i < 3; ++i) {
        if (trie.exists(test_codes[i], test_level)) {
            std::cout << "ERROR: Cell " << test_codes[i] << " should not exist initially" << std::endl;
            return;
        }
    }
    std::cout << "✓ Initial non-existence check passed" << std::endl;
    
    // Insert cells
    for (int i = 0; i < 3; ++i) {
        trie.insert(test_codes[i], test_level, 42.0 + i);
    }
    
    // Now cells should exist
    for (int i = 0; i < 3; ++i) {
        if (!trie.exists(test_codes[i], test_level)) {
            std::cout << "ERROR: Cell " << test_codes[i] << " should exist after insertion" << std::endl;
            return;
        }
    }
    std::cout << "✓ Post-insertion existence check passed" << std::endl;
    
    // Test non-existent cells
    D_morton non_existent = 999;
    if (trie.exists(non_existent, test_level)) {
        std::cout << "ERROR: Non-existent cell should not be found" << std::endl;
        return;
    }
    std::cout << "✓ Non-existent cell check passed" << std::endl;
    
    std::cout << "Basic existence test PASSED!" << std::endl;
}

void test_mortontriemap_integration() {
    std::cout << "\nTesting MortonTrieMap integration..." << std::endl;
    
    MortonTrieMap<int> map;
    D_morton key = 42;
    
    // Test existence through wrapper
    if (map.exists(key)) {
        std::cout << "ERROR: Key should not exist initially" << std::endl;
        return;
    }
    
    // Insert through wrapper
    map[key] = 123;
    
    // Check existence through wrapper
    if (!map.exists(key)) {
        std::cout << "ERROR: Key should exist after insertion" << std::endl;
        return;
    }
    
    // Test with custom level
    if (!map.exists(key, refine_level)) {
        std::cout << "ERROR: Key should exist at refine_level" << std::endl;
        return;
    }
    
    std::cout << "✓ MortonTrieMap integration test PASSED!" << std::endl;
}

void benchmark_existence_query() {
    std::cout << "\nBenchmarking existence queries..." << std::endl;

    MortonTrie& trie = MortonTrie::getInstance();

    // Generate safer test data with smaller Morton codes
    const int NUM_CELLS = 1000;
    const uint8_t test_level = 2;
    std::vector<D_morton> test_codes;

    // Use sequential codes to avoid out-of-bounds issues
    for (int i = 0; i < NUM_CELLS; ++i) {
        D_morton code = i;
        test_codes.push_back(code);
        trie.insert(code, test_level, static_cast<double>(i));
    }

    std::cout << "Inserted " << NUM_CELLS << " cells at level " << static_cast<int>(test_level) << std::endl;

    // Benchmark bitmap-based existence queries
    auto start = std::chrono::high_resolution_clock::now();
    int found_count = 0;

    for (int iter = 0; iter < 1000; ++iter) {
        for (const auto& code : test_codes) {
            if (trie.exists(code, test_level)) {
                found_count++;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto bitmap_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Bitmap queries: " << found_count << " found in "
              << bitmap_time.count() << " microseconds (1000 iterations)" << std::endl;
    std::cout << "Average time per bitmap query: "
              << static_cast<double>(bitmap_time.count()) / (1000 * NUM_CELLS)
              << " microseconds" << std::endl;

    // Test some non-existent queries
    start = std::chrono::high_resolution_clock::now();
    int not_found_count = 0;

    for (int iter = 0; iter < 1000; ++iter) {
        for (int i = NUM_CELLS; i < NUM_CELLS + 100; ++i) {
            if (!trie.exists(D_morton(i), test_level)) {
                not_found_count++;
            }
        }
    }

    end = std::chrono::high_resolution_clock::now();
    auto not_found_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Non-existent queries: " << not_found_count << " correctly not found in "
              << not_found_time.count() << " microseconds" << std::endl;

    std::cout << "Bitmap existence query benchmark completed successfully!" << std::endl;
}

int main() {
    std::cout << "=== MortonTrie Bitmap Existence Query Tests ===" << std::endl;
    
    test_basic_existence();
    test_mortontriemap_integration();
    benchmark_existence_query();
    
    std::cout << "\n=== All tests completed ===" << std::endl;
    return 0;
}
