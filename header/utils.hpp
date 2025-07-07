#pragma once

#include "settings.hpp"
#include <iostream>

#define two_pow(n) (1 << n)
#define eight_pow(n) (1 << (3*n))
constexpr double ddx = 1./two_pow(refine_level);


inline uint log2 (uint x) { 
  static const unsigned char log_2[256] = { 
    0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, 
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, 
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, 
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, 
    8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 
    8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 
    8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 
    8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8 
  }; 
  int l = -1; 
  while (x >= 256) { l += 8; x >>= 8; } 
  return l + log_2[x]; 
} 


#define EXTRACT_MORTON_CONVERT_INDEX(bitset, start_pos) \
    (static_cast<uint8_t>((bitset >> (start_pos - DIM)).to_ulong() & 0x7))

#define MAX_2(a,b) ((a) > (b) ? (a) : (b))
#define MAX_3(a,b,c) (((a) > (b)) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))

#include <chrono>
class Timer
{
public:
	Timer() : beg_(clock_::now()) {}
	void reset() { beg_ = clock_::now(); }
	double elapsed() const {
		return std::chrono::duration_cast<second_>
			(clock_::now() - beg_).count();
	}

private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> beg_;
};


#if (DIM == 3)
 #if (STENCIL_SIZE == 19)
/* D3Q19             
                         ___________________              ________15_________  
                        /:                 /|            /:                 /|
                       / :      5         / |          14 :               11 |
  z                   /__:______|_3_____ /  |         /________18________ /  |
 /|\ y               |   :      |/      |   |         |  10              |   7 
  | /                | 2-:------0-------|-1 |         |   :      .       |   |
  |/                 |   :     /|       |   |         |   :              |   |
  o------->x         |   :----/-|-------|---|         8   :-------17-----9---|
                     |  /    4  |       |  /          |  /               |   /
                     | /        6       | /           | 12               | 13  
                     |/_________________|/            |/_______16________|/    
*/
constexpr int ex[STENCIL_SIZE] = {
    0,  
    1, -1,  0,  0,  0,  0,  
    1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
constexpr int ey[STENCIL_SIZE] = {
    0,  
    0,  0,  1, -1,  0,  0,  
    1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1};
constexpr int ez[STENCIL_SIZE] = {
    0,  
    0,  0,  0,  0,  1, -1,  
    0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1};

constexpr unsigned int e_inv[STENCIL_SIZE] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};

 #elif (STENCIL_SIZE == 27)
/* D3Q27             
                             ___________________              ________15_________              25________________19 
                            /:                 /|            /:                 /|            /:                 /| 
                           / :      5         / |          14 :               11 |           / :                / | 
      z                   /__:______|_3_____ /  |         /________18________ /  |         23__:______________21  | 
     /|\ y               |   :      |/      |   |         |  10              |   7         |   :              |   | 
      | /                | 2-:------0-------|-1 |         |   :      .       |   |         |   :              |   | 
      |/                 |   :     /|       |   |         |   :              |   |         |   :              |   | 
      o------->x         |   :----/-|-------|---|         8   :-------17-----9---|         |   22-------------|--24 
                         |  /    4  |       |  /          |  /               |   /         |  /               |  /  
                         | /        6       | /           | 12               | 13          | /                | /   
                         |/_________________|/            |/_______16________|/            20_________________26    
*/
           
constexpr int ex[STENCIL_SIZE] = {
    0, 
    1, -1,  0,  0,  0,  0, 
    1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0, 
    1, -1,  1, -1, -1,  1, -1,  1};
constexpr int ey[STENCIL_SIZE] = {
    0, 
    0,  0,  1, -1, 0, 0, 
    1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,
    1, -1, -1,  1, -1,  1,  1, -1};
constexpr int ez[STENCIL_SIZE] = {
    0,  
    0,  0,  0,  0,  1, -1,  
    0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,
    1, -1,  1, -1,  1, -1,  1, -1};

constexpr unsigned int e_inv[STENCIL_SIZE] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};

 #endif
#endif

// constexpr D_morton xmorton_ref_0N = 17129119497016008704ULL;
// constexpr D_morton xmorton_ref_0P = 1317624576693538816ULL;
// constexpr D_morton xmorton_ref_1N = 17129119497016011776ULL;
// constexpr D_morton xmorton_ref_1P = 1317624576693539328ULL;
// constexpr D_morton xmorton_ref_2N = 17129119497016012160ULL;
// constexpr D_morton xmorton_ref_2P = 1317624576693539392ULL;
// constexpr D_morton xmorton_ref_3N = 17129119497016012208ULL;
// constexpr D_morton xmorton_ref_3P = 1317624576693539400ULL;
// constexpr D_morton xmorton_ref_4N = 17129119497016012214ULL;
// constexpr D_morton xmorton_ref_4P = 1317624576693539401ULL;

// constexpr D_morton ymorton_ref_0N = 15811494920322469888ULL;
// constexpr D_morton ymorton_ref_0P = 2635249153387077632ULL;
// constexpr D_morton ymorton_ref_1N = 15811494920322472448ULL;
// constexpr D_morton ymorton_ref_1P = 2635249153387078656ULL;
// constexpr D_morton ymorton_ref_2N = 15811494920322472768ULL;
// constexpr D_morton ymorton_ref_2P = 2635249153387078784ULL;
// constexpr D_morton ymorton_ref_3N = 15811494920322472808ULL;
// constexpr D_morton ymorton_ref_3P = 2635249153387078800ULL;
// constexpr D_morton ymorton_ref_4N = 15811494920322472813ULL;
// constexpr D_morton ymorton_ref_4P = 2635249153387078802ULL;

// constexpr D_morton zmorton_ref_0N = 13176245766935392256ULL;
// constexpr D_morton zmorton_ref_0P = 5270498306774155264ULL;
// constexpr D_morton zmorton_ref_1N = 13176245766935393792ULL;
// constexpr D_morton zmorton_ref_1P = 5270498306774157312ULL;
// constexpr D_morton zmorton_ref_2N = 13176245766935393984ULL;
// constexpr D_morton zmorton_ref_2P = 5270498306774157568ULL;
// constexpr D_morton zmorton_ref_3N = 13176245766935394008ULL;
// constexpr D_morton zmorton_ref_3P = 5270498306774157600ULL;
// constexpr D_morton zmorton_ref_4N = 13176245766935394011ULL;
// constexpr D_morton zmorton_ref_4P = 5270498306774157604ULL;

// constexpr D_morton ref_one_0 = 4096ULL;
// constexpr D_morton ref_one_1 = 512ULL;
// constexpr D_morton ref_one_2 = 64ULL;
// constexpr D_morton ref_one_3 = 8ULL;
// constexpr D_morton ref_one_4 = 1ULL;