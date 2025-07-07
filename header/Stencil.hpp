#pragma once

#include "mortonTrie.hpp"
#include "mortonHash.hpp"
#include "compressedMortonTrie.hpp"
#include "mesh.hpp"

class Stencil {
public:
    void compute(uint level);
    void stream(uint level);
    Stencil(Mesh& mesh_, MortonTrie& trie, D_mapInt& hashmap, CompressedMortonTrie& cmtrie) : mesh(mesh_), amr_trie(trie), amr_hashmap(hashmap), amr_CMtrie(cmtrie) {}

private:
    MortonTrie& amr_trie;
    D_mapInt& amr_hashmap;
    CompressedMortonTrie& amr_CMtrie;
    Mesh& mesh;
};