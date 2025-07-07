/**
* @file
* @brief  Define Morton assist Class used to store functions related to morton code calculation.
*.
*/
#ifndef MORTON_ASSIST_H
#define MORTON_ASSIST_H
#include "General.h"
class Morton_Assist
{
	friend class Grid_Background;
	friend class Grid_Intermediate;
	friend class Grid_Inner;
	friend class Grid_Manager;
	friend class Math_Manager;
private:
	static D_morton mortonx_min, mortonx_max;
	static D_morton mortony_min, mortony_max;
#if (C_DIMS == 3)
	static D_morton mortonz_min, mortonz_max;
#endif
	static std::array<D_morton, C_max_level + 1> ref_one;
	static std::array<std::array<D_morton, 2>, C_max_level + 1> xmorton_ref; ///< morton code used to find neighbours in x directions by bool operator
	static std::array<std::array<D_morton, 2>, C_max_level + 1> ymorton_ref; ///< morton code used to find neighbours in y directions by bool operator
#if (C_DIMS == 3)
	static std::array<std::array<D_morton, 2>, C_max_level + 1> zmorton_ref; ///< morton code used to find neighbours in z directions by bool operator
#endif
	static std::array<D_morton, C_max_level + 1> morton_xyz; ///< morton code where the x, y, z position at ilevel is set as 1, while others are 0
	static std::array<D_morton, C_max_level + 1> morton_xyz_flip; ///< morton code where the x, y, z position at ilevel is set as 0, while others are 1
public:
	// parameters used to find neighers of the background mesh
	static	D_morton xbk_ref1, xbk_ref0;
	static	D_morton ybk_ref1, ybk_ref0;
#if (C_DIMS==3)
	static	D_morton zbk_ref1, zbk_ref0;
#endif
	static Morton_Assist* pointer_me;    ///< pointer points to the class itself
	static unsigned int bit_background;  ///< Number of bits for numbering the backround nodes
	static unsigned int bit_otherlevel;  ///< Number of bits for other nodes, C_BIT = bit_otherlevel + bit_background
//	static D_morton bk_refine; ///< morton code with bit_background be 0 and bit_otherlevel be 1
public:
	static inline D_morton compute_ref_one(unsigned int ilevel);
//  static inline D_morton compute_current_cooridnate(unsigned int ilevel);
	// Notice: not sure which one is better, x_ref, yref as static parameters or formal parameters? Formal parameters are used here.
	static inline D_morton find_x0(const D_morton &key, const unsigned int ilevel); // find neighbour in the negative x direction
	static inline D_morton find_x1(const D_morton &key, const unsigned int ilevel); // find neighbour in the positive x direction
	static inline D_morton find_y0(const D_morton &key, const unsigned int ilevel); // find neighbour in the negative y direction
	static inline D_morton find_y1(const D_morton &key, const unsigned int ilevel); // find neighbour in the positive y direction
	static inline D_morton find_x0_check_boundary(const D_morton& key, const unsigned int ilevel);
	static inline D_morton find_x1_check_boundary(const D_morton& key, const unsigned int ilevel);
	static inline D_morton find_y0_check_boundary(const D_morton& key, const unsigned int ilevel);
	static inline D_morton find_y1_check_boundary(const D_morton& key, const unsigned int ilevel);
	static inline D_morton find_neighbor(const D_morton& key, const uint level, const int x_shift, const int y_shift, const int z_shift);
	// static inline void generate_neighbor_sequence_boxStencil(const D_morton& key, const uint level, std::array<D_morton,26>& neighbor_sequence);
	// static inline void generate_neighbor_sequence_starStencil(const D_morton& key, const uint level, std::array<D_morton,6>& neighbor_sequence);

	static void compute_ref(unsigned int ilevel);
#if (C_DIMS==2)
	static void compute_coordinate(D_morton key, unsigned int ilevel, D_real &x, D_real &y); // compute cooridate according to the morton code
	static inline D_morton find_neighbor(const D_morton key, const int level, const int x_shift, const int y_shift);
#endif
#if (C_DIMS==3)
	static inline D_morton find_z0(const D_morton &key, const unsigned int ilevel); // find neighbour in the negative z direction
	static inline D_morton find_z1(const D_morton &key, const unsigned int ilevel); // find neighbour in the positive z direction
	static inline D_morton find_z0_check_boundary(const D_morton& key, const unsigned int ilevel);
	static inline D_morton find_z1_check_boundary(const D_morton& key, const unsigned int ilevel);
	static void compute_coordinate(D_morton key, unsigned int ilevel, D_real &x, D_real &y, D_real &z); // compute cooridate according to the morton code
#endif
protected:
	void morton_initial();
public:
#if (C_DIMS==2)
	D_morton morton_encode(D_uint ix, D_uint iy);
#endif
#if (C_DIMS==3)
	D_morton morton_encode(D_uint ix, D_uint iy, D_uint iz);
#endif

};
/**
* @brief function to search neighbour using bitwise operator.
* @param[in]  ilevel   refinement level
* @return              reference Morton code only one sepecified position is true for a given refinement level
*/
inline D_morton Morton_Assist::compute_ref_one(unsigned int ilevel)
{
	D_morton ref_one = 0;
	ref_one.set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS, true); // reference value which is 1 at the sepecified position for a given level
	return ref_one;
}

///**
//* @brief function to search neighbour using bitwise operator.
//* @param[in]  ilevel   refinement level
//* @return              reference Morton code C_DIMS sepecified positions are true, used identify if the node is the first leaf
//*/
//inline D_morton Morton_Assist::compute_current_cooridnate(unsigned int ilevel)
//{
//	D_morton morton_temp = 0;
//	morton_temp.set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS, true);
//	morton_temp.set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS + 1, true);
//#if (C_DIMS==3)
//	morton_temp.set(Morton_Assist::bit_otherlevel - ilevel * C_DIMS + 2, true);
//#endif
//	return morton_temp;
//}

#if (C_DIMS==2)
/**
* @brief function to search neighbour using bitwise operator.
* @param[in]  key   Morton code of current node.
* @return           Morton code of the neighbour
*/
inline D_morton Morton_Assist::find_x0_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the negative x direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key & xmorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_x1_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the positive x direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key | xmorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y0_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the negative y direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key & ymorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y1_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the positive y direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key | ymorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)));
}

inline D_morton Morton_Assist::find_x0(const D_morton& key, const unsigned int ilevel) // find neighbour in the negative x direction
{
	D_morton temp = static_cast <D_morton> ((key & xmorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_x1(const D_morton& key, const unsigned int ilevel) // find neighbour in the positive x direction
{
	D_morton temp = static_cast <D_morton> ((key | xmorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y0(const D_morton& key, const unsigned int ilevel) // find neighbour in the negative y direction
{
	D_morton temp = static_cast <D_morton> ((key & ymorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y1(const D_morton& key, const unsigned int ilevel) // find neighbour in the positive y direction
{
	D_morton temp = static_cast <D_morton> ((key | ymorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)));
}
#endif

#if (C_DIMS==3)
/**
* @brief function to search neighbour using bitwise operator.
* @param[in]  key   Morton code of current node.
* @return           Morton code of the neighbour
*/
inline D_morton Morton_Assist::find_x0_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the negative x direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_min) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key & xmorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_x1_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the positive x direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::xmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonx_max) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key | xmorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y0_check_boundary(const D_morton &key, const unsigned int ilevel)  // find neighbour in the negative y direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_min) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key & ymorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y1_check_boundary(const D_morton &key, const unsigned int ilevel)  // find neighbour in the positive y direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::ymorton_ref.at(ilevel).at(1)) == Morton_Assist::mortony_max) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key | ymorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_z0_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the negative z direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_min) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key & zmorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & zmorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_z1_check_boundary(const D_morton &key, const unsigned int ilevel) // find neighbour in the positive z direction
{
#if (C_CHECK_MORTON_BOUNDARY == 1)
	if ((key & Morton_Assist::zmorton_ref.at(ilevel).at(1)) == Morton_Assist::mortonz_max) return key;
#endif
	D_morton temp = static_cast <D_morton> ((key | zmorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & zmorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}

inline D_morton Morton_Assist::find_x0(const D_morton& key, const unsigned int ilevel) // find neighbour in the negative x direction
{
	D_morton temp = static_cast <D_morton> ((key & xmorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_x1(const D_morton& key, const unsigned int ilevel) // find neighbour in the positive x direction
{
	D_morton temp = static_cast <D_morton> ((key | xmorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y0(const D_morton& key, const unsigned int ilevel)  // find neighbour in the negative y direction
{
	D_morton temp = static_cast <D_morton> ((key & ymorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_y1(const D_morton& key, const unsigned int ilevel)  // find neighbour in the positive y direction
{
	D_morton temp = static_cast <D_morton> ((key | ymorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & ymorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & zmorton_ref.at(ilevel).at(1)));
}

inline D_morton Morton_Assist::find_z0(const D_morton& key, const unsigned int ilevel) // find neighbour in the negative z direction
{
	D_morton temp = static_cast <D_morton> ((key & zmorton_ref.at(ilevel).at(1)).to_ullong() - ref_one.at(ilevel).to_ullong());
	return ((temp & zmorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}
inline D_morton Morton_Assist::find_z1(const D_morton& key, const unsigned int ilevel) // find neighbour in the positive z direction
{
	D_morton temp = static_cast <D_morton> ((key | zmorton_ref.at(ilevel).at(0)).to_ullong() + ref_one.at(ilevel).to_ullong());
	return ((temp & zmorton_ref.at(ilevel).at(1)) | (key & xmorton_ref.at(ilevel).at(1)) | (key & ymorton_ref.at(ilevel).at(1)));
}

inline D_morton Morton_Assist::find_neighbor(const D_morton& key, const uint level, const int x_shift, const int y_shift, const int z_shift)
{
	D_morton rst_code = key;
	if (x_shift == 1) {
		rst_code = Morton_Assist::find_x1(key, level);
	}
	else if (x_shift == -1)
	{
		rst_code = Morton_Assist::find_x0(key, level);
	}

	if (y_shift == 1) {
		rst_code = Morton_Assist::find_y1(rst_code, level);
	}
	else if (y_shift == -1)
	{
		rst_code = Morton_Assist::find_y0(rst_code, level);
	}

	if (z_shift == 1) {
		rst_code = Morton_Assist::find_z1(rst_code, level);
	}
	else if (z_shift == -1)
	{
		rst_code = Morton_Assist::find_z0(rst_code, level);
	}
	
	return rst_code;
}

// /**
//  * @brief generate the neighbor sequence of a cell
//  * @param key the morton code of the cell
//  * @param level the level of the cell
//  * @param neighbor_sequence the neighbor sequence of the cell
//  */
// inline void Morton_Assist::generate_neighbor_sequence_boxStencil(const D_morton& key, const uint level, std::array<D_morton,26>& neighbor_sequence)
// {
// 	for (int i = -1; i <= 1; i++) {
// 		for (int j = -1; j <= 1; j++) {
// 			for (int k = -1; k <= 1; k++) {
// 				if (i == 0 && j == 0 && k == 0) continue;
// 				neighbor_sequence[ (i + 1) + (j + 1) * 3 + (k + 1) * 9] = Morton_Assist::find_neighbor(key, level, i, j, k);
// 			}
// 		}
// 	}
// }

// inline void Morton_Assist::generate_neighbor_sequence_starStencil(const D_morton& key, const uint level, std::array<D_morton,6>& neighbor_sequence)
// {
// 	neighbor_sequence[0] = Morton_Assist::find_neighbor(key, level, -1, 0, 0);
// 	neighbor_sequence[1] = Morton_Assist::find_neighbor(key, level, 1, 0, 0);
// 	neighbor_sequence[2] = Morton_Assist::find_neighbor(key, level, 0, -1, 0);
// 	neighbor_sequence[3] = Morton_Assist::find_neighbor(key, level, 0, 1, 0);
// 	neighbor_sequence[4] = Morton_Assist::find_neighbor(key, level, 0, 0, -1);
// 	neighbor_sequence[5] = Morton_Assist::find_neighbor(key, level, 0, 0, 1);
// }

#endif

#endif