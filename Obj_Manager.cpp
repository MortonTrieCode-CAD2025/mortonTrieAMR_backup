/**
* @file
* @brief Main function.
* @note .
*/
#include "Obj_Manager.h"
#include <iomanip>
#include <memory>
#include <typeinfo>
#include <fstream>
#if(C_MAP_TYPE == 2)
#include <algorithm>
#include <numeric>
#include <list>
#include <forward_list>
#endif

Grid_Manager* Grid_Manager::pointer_me;
Lat_Manager* Lat_Manager::pointer_me;
IO_Manager* IO_Manager::pointer_me;
Solid_Manager* Solid_Manager::pointer_me;

// declarations
Grid_Manager gr_manager;
Solid_Manager solid_manager;
IO_Manager io_manager;
Lat_Manager lat_manager;

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#elif defined(__unix__) || defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#include <unistd.h>
#endif

#if (C_MAP_TYPE == 1)
// 内存测量工具类
class MapMemoryMeasurer {
public:
	// 方法1: 使用bucket信息估算内存占用
	template <typename K, typename V, typename Hash = std::hash<K>, 
				typename Pred = std::equal_to<K>, 
				typename Alloc = std::allocator<std::pair<const K, V>>>
	static size_t estimateMemoryUsage(const std::unordered_map<K, V, Hash, Pred, Alloc>& map) {
		// 基本容器大小
		size_t memoryUsage = sizeof(std::unordered_map<K, V, Hash, Pred, Alloc>);
		
		// 添加bucket数组的大小
		size_t bucketCount = map.bucket_count();
		memoryUsage += bucketCount * sizeof(void*);
		
		// 添加节点的大小
		size_t nodeSize = sizeof(std::pair<const K, V>) + sizeof(void*); // 节点包含pair和指针
		memoryUsage += map.size() * nodeSize;
		
		// 对于string键，还需要考虑字符串本身的内存
		if constexpr(std::is_same_v<K, std::string>) {
			for (const auto& pair : map) {
				// 字符串的容量(通常比size大)加上一个指针和size_t (SSO优化除外)
				memoryUsage += pair.first.capacity() + 1 + sizeof(size_t) + sizeof(void*);
			}
		}
		
		// 对于string值，也需要考虑字符串内存
		if constexpr(std::is_same_v<V, std::string>) {
			for (const auto& pair : map) {
				memoryUsage += pair.second.capacity() + 1 + sizeof(size_t) + sizeof(void*);
			}
		}
		
		return memoryUsage;
	}
	
	// 方法2: 更详细的内存分析
	template <typename K, typename V, typename Hash = std::hash<K>, 
				typename Pred = std::equal_to<K>, 
				typename Alloc = std::allocator<std::pair<const K, V>>>
	static void analyzeMapMemory(const std::unordered_map<K, V, Hash, Pred, Alloc>& map) {
		// 获取基本信息
		size_t size = map.size();
		size_t bucketCount = map.bucket_count();
		float loadFactor = size / static_cast<float>(bucketCount > 0 ? bucketCount : 1);
		
		// 计算估计内存
		size_t estimatedMemory = estimateMemoryUsage(map);
		
		// 输出分析结果
		std::cout << "===== unordered_map 内存分析 =====\n";
		std::cout << "类型: unordered_map<" 
					<< typeid(K).name() << ", " 
					<< typeid(V).name() << ">\n";
		std::cout << "元素数量: " << size << "\n";
		std::cout << "桶数量: " << bucketCount << "\n";
		std::cout << "负载因子: " << loadFactor << "\n";
		std::cout << "最大负载因子: " << map.max_load_factor() << "\n";
		std::cout << "估计内存占用: " << estimatedMemory << " 字节 (" 
					<< (estimatedMemory / 1024.0) << " KB)\n";
		std::cout << "每个元素平均内存: " 
					<< (size > 0 ? estimatedMemory / static_cast<double>(size) : 0) 
					<< " 字节\n";
		
		// 分析桶分布
		size_t emptyBuckets = 0;
		size_t maxBucketSize = 0;
		size_t totalChainLength = 0;
		
		for (size_t i = 0; i < bucketCount; ++i) {
			size_t bucketSize = map.bucket_size(i);
			if (bucketSize == 0) {
				emptyBuckets++;
			}
			if (bucketSize > maxBucketSize) {
				maxBucketSize = bucketSize;
			}
			totalChainLength += bucketSize;
		}
		
		std::cout << "空桶数量: " << emptyBuckets << " (" 
					<< (bucketCount > 0 ? emptyBuckets * 100.0 / bucketCount : 0) << "%)\n";
		std::cout << "最大链长: " << maxBucketSize << "\n";
		std::cout << "平均链长: " << (bucketCount > 0 ? totalChainLength / static_cast<double>(bucketCount) : 0) << "\n";
		std::cout << "===================================\n";
	}
	
	// 方法3: 通过拷贝前后的内存差测量 (需要操作系统支持)
	template <typename MapType>
	static size_t measureByDifference(const MapType& map) {
		size_t beforeMem = getCurrentRSS();
		
		// 创建map的拷贝，这会分配新的内存
		MapType mapCopy = map;
		
		size_t afterMem = getCurrentRSS();
		
		// 返回差值，这应该接近map的实际内存占用
		return afterMem - beforeMem;
	}
	
private:
	// 获取当前进程的内存使用量 (Resident Set Size)
	static size_t getCurrentRSS() {
#if defined(_WIN32)
		// Windows实现
		PROCESS_MEMORY_COUNTERS pmc;
		if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
			return pmc.WorkingSetSize;
		}
		return 0;
#elif defined(__unix__) || defined(__linux__) || defined(__APPLE__)
		// Unix/Linux/Mac实现
		FILE* file = fopen("/proc/self/statm", "r");
		if (file) {
			long rss = 0;
			if (fscanf(file, "%*s %ld", &rss) == 1) {
				fclose(file);
				return rss * sysconf(_SC_PAGESIZE);
			}
			fclose(file);
		}
		
		// 如果上面的方法失败，尝试使用getrusage
		struct rusage usage;
		if (getrusage(RUSAGE_SELF, &usage) == 0) {
			return usage.ru_maxrss * 1024;
		}
		return 0;
#else
		return 0; // 不支持的平台
#endif
	}
};
#elif(C_MAP_TYPE == 2)
class MapMemoryMeasurer {
	public:
		// 测量std::map的内存占用
		template <typename K, typename V, typename Compare = std::less<K>,
				  typename Alloc = std::allocator<std::pair<const K, V>>>
		static size_t estimateMemoryUsage(const std::map<K, V, Compare, Alloc>& map) {
			// 基本容器大小
			size_t memoryUsage = sizeof(std::map<K, V, Compare, Alloc>);
			
			// std::map通常使用红黑树实现，每个节点包含pair、颜色标记和三个指针(父、左、右)
			size_t nodeOverhead = 3 * sizeof(void*) + sizeof(char); // 三个指针和颜色标记
			size_t nodeSize = sizeof(std::pair<const K, V>) + nodeOverhead;
			
			// 添加所有节点的内存
			memoryUsage += map.size() * nodeSize;
			
			// 处理string类型的键
			if constexpr(std::is_same_v<K, std::string>) {
				for (const auto& pair : map) {
					memoryUsage += estimateStringMemory(pair.first);
				}
			}
			
			// 处理string类型的值
			if constexpr(std::is_same_v<V, std::string>) {
				for (const auto& pair : map) {
					memoryUsage += estimateStringMemory(pair.second);
				}
			}
			
			// 处理嵌套容器
			if constexpr(is_container_v<V>) {
				for (const auto& pair : map) {
					memoryUsage += estimateContainerMemory(pair.second);
				}
			}
			
			return memoryUsage;
		}
		template <typename K, typename V, typename Compare = std::less<K>,
              typename Alloc = std::allocator<std::pair<const K, V>>>
    static void analyzeMapMemory(const std::map<K, V, Compare, Alloc>& map) {
        // 获取基本信息
        size_t size = map.size();
        
        // 计算估计内存
        size_t estimatedMemory = estimateMemoryUsage(map);
        
        // 输出分析结果
        std::cout << "===== std::map 内存分析 =====\n";
        std::cout << "类型: map<" 
                  << typeid(K).name() << ", " 
                  << typeid(V).name() << ">\n";
        std::cout << "元素数量: " << size << "\n";
        std::cout << "估计内存占用: " << estimatedMemory << " 字节 (" 
                  << (estimatedMemory / 1024.0) << " KB)\n";
        std::cout << "每个元素平均内存: " 
                  << (size > 0 ? estimatedMemory / static_cast<double>(size) : 0) 
                  << " 字节\n";
        
        // 红黑树高度估计 (log2(n))
        if (size > 0) {
            double estimatedHeight = std::log2(size);
            std::cout << "估计树高度: ~" << estimatedHeight << "\n";
        }
        
        std::cout << "===================================\n";
    }
	// 通过拷贝前后的内存差测量
    template <typename MapType>
    static size_t measureByDifference(const MapType& map) {
        // 强制GC
        std::malloc(1);
        std::free(std::malloc(1));
        
        size_t beforeMem = getCurrentRSS();
        
        // 创建map的拷贝，这会分配新的内存
        MapType mapCopy = map;
        
        // 确保所有内存被真正分配
        volatile char dummy = 0;
        for (const auto& pair : mapCopy) {
            dummy += reinterpret_cast<const char*>(&pair)[0];
        }
        
        // 强制GC
        std::malloc(1);
        std::free(std::malloc(1));
        
        size_t afterMem = getCurrentRSS();
        
        // 返回差值，这应该接近map的实际内存占用
        return afterMem - beforeMem;
    }
	private:
    // 检测类型是否为容器的trait
    template <typename T, typename = void>
    struct is_container : std::false_type {};
    
    template <typename T>
    struct is_container<T, 
        std::void_t<
            typename T::value_type,
            typename T::size_type,
            typename T::iterator,
            decltype(std::declval<T>().size()),
            decltype(std::declval<T>().begin()),
            decltype(std::declval<T>().end())
        >> : std::true_type {};
    
    template <typename T>
    static constexpr bool is_container_v = is_container<T>::value;
    
    // 估计std::string的内存占用
    static size_t estimateStringMemory(const std::string& str) {
        // 检查是否使用SSO (Small String Optimization)
        // 这个阈值因实现而异，通常在15-23字节之间
        constexpr size_t ssoThreshold = 15; // 假设SSO阈值为15
        
        if (str.size() <= ssoThreshold) {
            return 0; // 字符串已包含在string对象本身中
        } else {
            // 非SSO字符串需要额外的堆分配
            // 计算实际分配的容量(通常是2的幂次或接近的值)
            return str.capacity() + 1; // +1 for null terminator
        }
    }
    
    // 估计容器的内存占用
    template <typename Container>
    static size_t estimateContainerMemory(const Container& container) {
        using ValueType = typename Container::value_type;
        
        size_t memoryUsage = 0;
        
        // 对于vector和类似容器，考虑容量
        if constexpr(has_capacity<Container>::value) {
            memoryUsage += container.capacity() * sizeof(ValueType);
        } 

		// 对于list和类似容器，考虑节点开销
		else if constexpr(std::is_same_v<Container, std::list<ValueType>> || 
			std::is_same_v<Container, std::forward_list<ValueType>>) {
		  // 每个节点包含值和至少一个指针
		  memoryUsage += container.size() * (sizeof(ValueType) + sizeof(void*));
	  }
	  // 对于set/map等树形容器
	  else if constexpr(std::is_same_v<Container, std::set<ValueType>> || 
						std::is_same_v<Container, std::multiset<ValueType>> ||
						std::is_same_v<Container, std::map<typename Container::key_type, typename Container::mapped_type>>) {
		  // 红黑树节点通常包含值、颜色标记和三个指针
		  memoryUsage += container.size() * (sizeof(ValueType) + sizeof(char) + 3 * sizeof(void*));
	  }
	  // 默认情况
	  else {
		  memoryUsage += container.size() * sizeof(ValueType);
	  }
	  
	  // 递归处理容器内的string或嵌套容器
	  if constexpr(is_string_or_container_v<ValueType>) {
		  for (const auto& item : container) {
			  memoryUsage += estimateNestedMemory(item);
		  }
	  }
	  
	  return memoryUsage;
  }
  
  // 处理嵌套的string或容器
  template <typename T>
  static size_t estimateNestedMemory(const T& value) {
	  if constexpr(std::is_same_v<T, std::string>) {
		  return estimateStringMemory(value);
	  }
	  else if constexpr(is_container_v<T>) {
		  return estimateContainerMemory(value);
	  }
	  else {
		  return 0;
	  }
  }
  
  // 检测类型是否为string或容器
  template <typename T>
  static constexpr bool is_string_or_container_v = 
	  std::is_same_v<T, std::string> || is_container_v<T>;
  
  // 检测容器是否有capacity方法
  template <typename T, typename = void>
  struct has_capacity : std::false_type {};
  
  template <typename T>
  struct has_capacity<T, std::void_t<decltype(std::declval<T>().capacity())>> 
	  : std::true_type {};
  
  // 获取当前进程的内存使用量 (Resident Set Size)
  static size_t getCurrentRSS() {
#ifdef _WIN32
	  // Windows实现
	  PROCESS_MEMORY_COUNTERS pmc;
	  if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
		  return pmc.WorkingSetSize;
	  }
	  return 0;
#elif defined(__unix__) || defined(__linux__) || defined(__APPLE__)
	  // Unix/Linux/Mac实现
	  FILE* file = fopen("/proc/self/statm", "r");
	  if (file) {
		  long rss = 0;
		  if (fscanf(file, "%*s %ld", &rss) == 1) {
			  fclose(file);
			  return rss * sysconf(_SC_PAGESIZE);
		  }
		  fclose(file);
	  }
	  
	  // 如果上面的方法失败，尝试使用getrusage
	  struct rusage usage;
	  if (getrusage(RUSAGE_SELF, &usage) == 0) {
		  return usage.ru_maxrss * 1024;
	  }
	  return 0;
#else
	  return 0; // 不支持的平台
#endif
  }

};
#endif

void Obj_Manager::initial()
{

	// Initialize solid geometry
	std::vector<Ini_Shape> ini_shape(1); // number of solids

	// ini_shape.at(0).shape_type = geofile;
	ini_shape.at(0).shape_type = geofile_stl;
	// ini_shape.at(0).bool_moving = true;
	ini_shape.at(0).bool_moving = false;
#ifdef SOLIDCENTER
	ini_shape.at(0).length.push_back(C_xb / 2);
    ini_shape.at(0).length.push_back(C_yb / 2);
	ini_shape.at(0).length.push_back(C_zb / 2);
#endif

	// ini_shape.at(0).shape_type = circle;
	// ini_shape.at(0).bool_moving = true;
	// ini_shape.at(0).numb_nodes = 2000;
	// ini_shape.at(0).length.push_back(C_xb / 2);
	// ini_shape.at(0).length.push_back(C_yb / 2);
	// ini_shape.at(0).length.push_back(1);

	// ini_shape.at(1).shape_type = geofile_stl;
	// // ini_shape.at(0).bool_moving = true;
	// ini_shape.at(1).bool_moving = false;
	// ini_shape.at(1).length.push_back(C_xb*2 / 3);
    // ini_shape.at(1).length.push_back(C_yb*2 / 3);
	// ini_shape.at(1).length.push_back(C_zb*2 / 3);

	//ini_shape.at(1).shape_type = line_fillx;
	//ini_shape.at(1).x0 = 20.;
	//ini_shape.at(1).y0 = 20.;
	//ini_shape.at(1).numb_nodes = 200;
	//ini_shape.at(1).length.push_back(-1.);
	//ini_shape.at(1).length.push_back(20.);

	Solid_Manager::pointer_me = &solid_manager;
	solid_manager.initial(ini_shape);
	// solid_manager.output_pointCloud(ini_shape);

	std::vector <Ini_Shape>().swap(ini_shape);

	// Generate and initialize mesh
	Grid_Manager::pointer_me = &gr_manager;
	Timer tmmr;
	double st = tmmr.elapsed();
	gr_manager.initial();
	double ed = tmmr.elapsed();

	std::cout << "Grid Size: " << gr_manager.gr_inner.grid.size() << std::endl;
	std::cout << "Grid space: " << gr_manager.gr_inner.get_dx() << std::endl;
	std::cout << "Generation time: " << ed - st << " s." << std::endl;

	MapMemoryMeasurer::analyzeMapMemory(gr_manager.gr_inner.grid);

	size_t diffMemory = MapMemoryMeasurer::measureByDifference(gr_manager.gr_inner.grid);
	std::cout << "差值法测量内存: " << diffMemory << " 字节 (" 
	<< (diffMemory / 1024.0) << " KB)\n\n";


	auto now = std::chrono::system_clock::now();
	auto now_time_t = std::chrono::system_clock::to_time_t(now);
	std::tm local_tm{};
    #ifdef _WIN32
    localtime_s(&local_tm, &now_time_t);
    #else
    localtime_r(&now_time_t, &local_tm);  // POSIX标准
    #endif
	std::stringstream ss;
    ss << OUTPUT_NAME << "_dx" << gr_manager.gr_inner.get_dx() << "_"
       << std::put_time(&local_tm, "%Y%m%d_%H%M%S") << ".txt";



	std::ofstream amr_nodes_file;
	amr_nodes_file.open(ss.str(), std::ios::out);
	for (auto iter =  gr_manager.gr_inner.grid.begin(); iter != gr_manager.gr_inner.grid.end(); ++iter) {
		amr_nodes_file << iter->first << "\n";
	}
	amr_nodes_file.close();
	

	// Generate lattices from mesh (node = vertex)
	// Lat_Manager::pointer_me = &lat_manager;
	// lat_manager.voxelize();
	// lat_manager.initial(); // old version
}

void Obj_Manager::time_marching_management()
{
	D_real sum_t = 0.;
	std::array<D_mapint, C_max_level + 1> map_add_nodes, map_remove_nodes;
	while (sum_t < 3 * C_dx)
	{
		time_marching(sum_t, map_add_nodes, map_remove_nodes);
		sum_t += C_dx;
	}

	
}

void Obj_Manager::time_marching(D_real sum_t, std::array<D_mapint, C_max_level + 1>  &map_add_nodes, std::array<D_mapint, C_max_level + 1>  &map_remove_nodes)
{
	D_real dt = C_dx / static_cast<D_real> (two_power_n(C_max_level));
	std::array<unsigned int, C_max_level + 1> accumulate_t{};
	std::array<unsigned int, C_max_level + 1> flag_time_step{}; // record number of time step at ilevel

	for (std::vector<unsigned int>::iterator iter = run_order.begin(); iter != run_order.end(); ++iter)
	{
		Timer tmr;
		double t0 = tmr.elapsed();

		unsigned int ilevel = *iter;
		++flag_time_step[ilevel];

		D_real current_t;

		accumulate_t[ilevel] += two_power_n(C_max_level - ilevel);
		current_t = sum_t + dt * static_cast<D_real> (accumulate_t[ilevel]);

		if (ilevel == C_max_level)
		{
#if (C_SOLID_BOUNDARY == 2)
			// update shape information
			unsigned int numb_solids = Solid_Manager::pointer_me->numb_solids;
			for (unsigned int ishape = 0; ishape < numb_solids; ++ishape)
			{
				if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_moving)
				{
					Solid_Manager::pointer_me->renew(ishape, current_t);				
					Grid_Manager::pointer_me->update_nodes_near_solid(ishape, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));
					Grid_Manager::pointer_me->update_map_node_IB(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));

					if (Solid_Manager::pointer_me->shape_solids.at(ishape).bool_enclosed)
					{
						Grid_Manager::pointer_me->update_ghost_node(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel));
					}
				}

			}
			
			if ((flag_time_step[ilevel] % 2) == 0)
			{
				Grid_Manager::pointer_me->call_update_nodes(ilevel, map_add_nodes.at(ilevel), map_remove_nodes.at(ilevel), map_add_nodes.at(ilevel - 1), map_remove_nodes.at(ilevel - 1));
				map_add_nodes[ilevel].clear();
				map_remove_nodes[ilevel].clear();
			}

#endif			

		}
		else if (ilevel == 0)
		{

		}
		else
		{
#if (C_SOLID_BOUNDARY == 2)
			if ((flag_time_step[ilevel] % 2) == 0)
			{
				Grid_Manager::pointer_me->call_update_nodes(ilevel, map_add_nodes[ilevel], map_remove_nodes[ilevel], map_add_nodes[ilevel - 1], map_remove_nodes[ilevel - 1]);
				map_add_nodes[ilevel].clear();
				map_remove_nodes[ilevel].clear();
			}
#endif	
		}

		//double t1 = tmr.elapsed();
		//double t2 = tmr.elapsed();
		//std::cout << ilevel<< ": " << t1 - t0 << ", " << t2 - t1 << std::endl;
	}
	
}

void Obj_Manager::output()
{
	// time
	auto now = std::chrono::system_clock::now();
	auto now_time_t = std::chrono::system_clock::to_time_t(now);
	std::tm local_tm{};
    #ifdef _WIN32
    localtime_s(&local_tm, &now_time_t);
    #else
    localtime_r(&now_time_t, &local_tm);  // POSIX标准
    #endif
	std::stringstream ss;
    ss << OUTPUT_NAME << "_"
       << std::put_time(&local_tm, "%Y%m%d_%H%M%S");

	// Write flowfield
	IO_Manager::pointer_me = &io_manager;
	io_manager.method = 4;
	io_manager.outfile = ss.str();

	std::vector<unsigned int> out_vlevel = {C_max_level};
	// std::vector<unsigned int> out_vlevel = {3};
	// std::vector<unsigned int> out_vlevel;
	// for (unsigned ilevel = 0; ilevel < C_max_level + 1; ++ilevel)
	// {
	// 	out_vlevel.push_back(ilevel);
	// }
	io_manager.vlevel = out_vlevel;

	io_manager.control();
}