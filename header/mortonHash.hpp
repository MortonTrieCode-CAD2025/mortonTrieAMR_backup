#pragma once

#include "Morton_Assist.hpp"
#include "settings.hpp"
#ifdef USE_HASHMAP
 #include <unordered_map>
#endif
#ifdef USE_RBTREE
 #include <map>
#endif
#ifdef USE_MORTONTRIE
 #include "mortonTrie.hpp"
#endif

struct Node
{
    uint8_t flag = 0;
};


#ifdef USE_HASHMAP
template <typename value>
using D_map_define = std::unordered_map<D_morton, value>;
#endif

#ifdef USE_RBTREE
class compare
{
public:
	bool operator()(const D_morton &a, const D_morton &b) const
	{
		return a.to_ullong() < b.to_ullong();
	}
};
template <typename value>
using D_map_define = std::map<D_morton, value, compare>;
#endif

#ifdef USE_MORTONTRIE
template <typename value>
using D_map_define = MortonTrieMap<value>;
#endif

#if defined(USE_HASHMAP) || defined(USE_RBTREE) || defined(USE_MORTONTRIE)
using D_map = D_map_define<Node>;
using D_mapInt = D_map_define<int>;
#else
using D_map = std::vector<Node>;
using D_mapInt = std::vector<int>;
#endif