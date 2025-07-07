/**
 * @file LBM_Manager.cpp
 * @brief
 * @version 0.1
 * @date 2023-03-27
 *
 */

#include "General.h"
#include "Constants.h"
#include "user.h"
#include "LBM_Manager.hpp"
// #include "Grid_Manager.h"
#include "Lat_Manager.hpp"
#include "Morton_Assist.h"
#include "Grid_Manager.h"
// #include <cmath>
// #include "Obj_Manager.h"
#include "lbm_kernel/Collision_Model_Factory.hpp"
#include "lbm_kernel/Stream_Manager.hpp"
#include "ComplexBoundary_Manager.hpp"
#include "BC/D3Q19_BC_Manager.hpp"
// #include "BC/ZouHe_Velocity.hpp"
// #include "BC/Periodic.hpp"
#include "BC/NonEquilibrumExtrapolation.hpp"
#include <chrono>

/**
 * @brief Construct a new LBM_Manager::LBM_Manager object
 */
LBM_Manager::LBM_Manager()
{
    LBM_Manager::pointer_me = this;

    // dermine run order of each refinement level
    defineRunOrder();
    setUnit("Cumulant");

    user = new User[C_max_level+1];
    allocMacrosAndDDF();

    // collideModel = Collision_Model_Factory::createCollisionModel("LBGK");
    collideModel = Collision_Model_Factory::createCollisionModel("Cumulant");
    streamModel = new ABPatternPull(user);
    // solidBCModel = new UniformBoundary(user);
    solidBCModel = new SinglePointInterpolation(user);

    bc_manager_ = new D3Q19_BoundaryCondtion_Manager;
    // std::unique_ptr<D3Q19_BC_Strategy> westBC_param_ptr = 
    //     std::unique_ptr<nonEquilibrumBounceBack_Velocity_West> (new nonEquilibrumBounceBack_Velocity_West(D_uvw(U_u0, 0., 0.)));
    std::unique_ptr<D3Q19_BC_Strategy> westBC_param_ptr = 
        std::unique_ptr<nonEquilibrumExtrapolation_West> (new nonEquilibrumExtrapolation_West(D_uvw(U_u0, U_v0, U_w0)));
    std::unique_ptr<D3Q19_BC_Strategy> eastBC_param_ptr = 
        std::unique_ptr<nonEquilibrumExtrapolation_East> (new nonEquilibrumExtrapolation_East(D_uvw(U_u0, U_v0, U_w0)));
    std::unique_ptr<D3Q19_BC_Strategy> southBC_param_ptr = 
        std::unique_ptr<nonEquilibrumExtrapolation_South> (new nonEquilibrumExtrapolation_South(D_uvw(U_u0, U_v0, U_w0)));
    std::unique_ptr<D3Q19_BC_Strategy> northBC_param_ptr = 
        std::unique_ptr<nonEquilibrumExtrapolation_North> (new nonEquilibrumExtrapolation_North(D_uvw(U_u0, U_v0, U_w0)));
    std::unique_ptr<D3Q19_BC_Strategy> topBC_param_ptr = 
        std::unique_ptr<nonEquilibrumExtrapolation_Top> (new nonEquilibrumExtrapolation_Top(D_uvw(U_u0, U_v0, U_w0)));
    std::unique_ptr<D3Q19_BC_Strategy> botBC_param_ptr = 
        std::unique_ptr<nonEquilibrumExtrapolation_Bot> (new nonEquilibrumExtrapolation_Bot(D_uvw(U_u0, U_v0, U_w0)));

    // Another style with same effect
    // std::unique_ptr<nonEquilibrumBounceBack_Velocity_West> westBC_param_ptr =
    //     std::make_unique<nonEquilibrumBounceBack_Velocity_West>(D_uvw(6.751,U_v0,U_w0));

    // bc_manager_->setup("West", std::move(westBC_param_ptr));
    //  bc_manager_->setup("East", std::make_unique<periodic_East>());
    //   bc_manager_->setup("South", std::make_unique<periodic_South>());
    //    bc_manager_->setup("North", std::make_unique<periodic_North>());
    //     bc_manager_->setup("Bot", std::make_unique<periodic_Bot>());
    //      bc_manager_->setup("Top", std::make_unique<periodic_Top>());
    bc_manager_->setup("West", std::move(westBC_param_ptr));
     bc_manager_->setup("East", std::move(eastBC_param_ptr));
      bc_manager_->setup("South", std::move(southBC_param_ptr));
       bc_manager_->setup("North", std::move(northBC_param_ptr));
        bc_manager_->setup("Bot", std::move(botBC_param_ptr));
         bc_manager_->setup("Top", std::move(topBC_param_ptr));

        //   std::cout << "bc_manager_->boundaryConditionStrategies_.size() " << bc_manager_->boundaryConditionStrategies_.size() << std::endl;
        //   std::cout << bc_manager_->density_ptr->size() << std::endl;

          /** Other LBM programs (eg. XLB or LUMA) donot treat the edge and corner
           *   Too much fuss about boundary conditions wasted my much time
           */
        //   bc_manager_->setup({"NorthEast", "NorthWest", "SouthEast", "SouthWest", 
        //                       "TopNorth", "TopSouth", "BotNorth", "BotSouth",
        //                       "TopEast", "TopWest", "BotEast", "BotWest"}, std::make_unique<nonEquilibrumBounceBack_NoSlip_Edge>());
        //    bc_manager_->setup({"All"}, std::make_unique<nonEquilibrumBounceBack_NoSlip_Corner>());

    /** @todo Shall we have a nonEquilibrumBounceBack_NoSlip_Face like Edges? */

}

void LBM_Manager::defineRunOrder()
{
    std::array<unsigned int, C_max_level + 1> accumulate_t = {};
	if (C_max_level == 1)
	{
		run_order.push_back(0);
		run_order.push_back(C_max_level);
		run_order.push_back(C_max_level);
	} else {
		unsigned int ilevel = 1;
		accumulate_t[0] = two_power_n(C_max_level);
		run_order.push_back(0);
		while (ilevel > 0)
		{
			if (accumulate_t[ilevel- 1] == accumulate_t[ilevel])
			{
				--ilevel;
			} else {
				run_order.push_back(ilevel);
				accumulate_t[ilevel] += two_power_n(C_max_level - ilevel);
				
				if (ilevel + 1 < C_max_level)
				{
					++ilevel;
				} else {
					run_order.push_back(ilevel + 1);
					run_order.push_back(ilevel + 1);
					accumulate_t[ilevel + 1] += two_power_n(C_max_level - ilevel);
					++ilevel;
				}
			}
		}
	}
}

/**
 * @param[in] rho 
 * @param[in] velocity 
 * @param[out] feq 
 */
void LBM_Manager::calculateFeq(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq)
{
#if (C_DIMS == 3)
 #if (C_Q == 19)
    #define wc (1./3.)
    #define we (1./36.)
    #define ws (1./18.)
    constexpr D_Phy_real w[C_Q] = {
        wc, 
        ws, ws, ws, ws, ws, ws, 
        we, we, we, we, we, we, we, we, we, we, we, we};

 #elif (C_Q == 27)
    #define wb (8./27.)
    #define ws (2./27.)
    #define we (1./54.)
    #define wc (1./216.)
    constexpr D_Phy_real w[C_Q] = {
        wb,
        ws, ws, ws, ws, ws, ws,
        we, we, we, we, we, we, we, we, we, we, we, we,
        wc, wc, wc, wc, wc, wc, wc, wc};
 #endif
#endif

    D_Phy_Velocity uu = dot_product(velocity, velocity);

    for (D_int i_dir = 0; i_dir < C_Q; ++i_dir)
    {
        D_Phy_Velocity euv = ex[i_dir] * velocity.x + ey[i_dir] * velocity.y + ez[i_dir] * velocity.z;
        feq[i_dir] = rho * w[i_dir] * (1. + 3.*euv + 4.5*euv*euv - 1.5*uu);
    }
}

void LBM_Manager::calculateFeq_HeLuo(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq)
{
#if (C_DIMS == 3)
 #if (C_Q == 19)
    #define wc (1./3.)
    #define we (1./36.)
    #define ws (1./18.)
    constexpr D_Phy_real w[C_Q] = {
        wc, 
        ws, ws, ws, ws, ws, ws, 
        we, we, we, we, we, we, we, we, we, we, we, we};

 #elif (C_Q == 27)
    #define wb (8./27.)
    #define ws (2./27.)
    #define we (1./54.)
    #define wc (1./216.)
    constexpr D_Phy_real w[C_Q] = {
        wb,
        ws, ws, ws, ws, ws, ws,
        we, we, we, we, we, we, we, we, we, we, we, we,
        wc, wc, wc, wc, wc, wc, wc, wc};
 #endif
#endif

    D_Phy_Velocity uu = dot_product(velocity, velocity);

    for (D_int i_dir = 0; i_dir < C_Q; ++i_dir)
    {
        D_Phy_Velocity euv = ex[i_dir] * velocity.x + ey[i_dir] * velocity.y + ez[i_dir] * velocity.z;
        feq[i_dir] = w[i_dir] * (rho + U_rho0 * (3.*euv + 4.5*euv*euv - 1.5*uu));
    }
}

D_Phy_DDF LBM_Manager::calculateFeq(const D_Phy_Rho rho, const D_uvw velocity, const int dir)
{
#if (C_DIMS == 3)
 #if (C_Q == 19)
    #define wc (1./3.)
    #define we (1./36.)
    #define ws (1./18.)
    constexpr D_Phy_real w[C_Q] = {
        wc, 
        ws, ws, ws, ws, ws, ws, 
        we, we, we, we, we, we, we, we, we, we, we, we};

 #elif (C_Q == 27)
    #define wb (8./27.)
    #define ws (2./27.)
    #define we (1./54.)
    #define wc (1./216.)
    constexpr D_Phy_real w[C_Q] = {
        wb,
        ws, ws, ws, ws, ws, ws,
        we, we, we, we, we, we, we, we, we, we, we, we,
        wc, wc, wc, wc, wc, wc, wc, wc};
 #endif
#endif

    D_Phy_Velocity uu = dot_product(velocity, velocity);
    D_Phy_Velocity euv = ex[dir] * velocity.x + ey[dir] * velocity.y + ez[dir] * velocity.z;
    return (rho * w[dir] * (1. + 3.*euv + 4.5*euv*euv - 1.5*uu));
}

D_Phy_DDF LBM_Manager::calculateFeq_HeLuo(const D_Phy_Rho rho, const D_uvw velocity, const int dir)
{
#if (C_DIMS == 3)
 #if (C_Q == 19)
    #define wc (1./3.)
    #define we (1./36.)
    #define ws (1./18.)
    constexpr D_Phy_real w[C_Q] = {
        wc, 
        ws, ws, ws, ws, ws, ws, 
        we, we, we, we, we, we, we, we, we, we, we, we};

 #elif (C_Q == 27)
    #define wb (8./27.)
    #define ws (2./27.)
    #define we (1./54.)
    #define wc (1./216.)
    constexpr D_Phy_real w[C_Q] = {
        wb,
        ws, ws, ws, ws, ws, ws,
        we, we, we, we, we, we, we, we, we, we, we, we,
        wc, wc, wc, wc, wc, wc, wc, wc};
 #endif
#endif

    D_Phy_Velocity uu = dot_product(velocity, velocity);
    D_Phy_Velocity euv = ex[dir] * velocity.x + ey[dir] * velocity.y + ez[dir] * velocity.z;
    return (w[dir] * (rho + U_rho0 * (3.*euv + 4.5*euv*euv - 1.5*uu)));
}

void LBM_Manager::initialFeq()
{
    // D_Phy_DDF feq_temp[C_Q];

    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        // D_int num_lattice = Lat_Manager::pointer_me->lat_f[i_level].size();
        // for (D_int i_latt = 0; i_latt < num_lattice; ++i_latt)
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        for (auto lat_iter : *lat_ptr)
        {
            // computeLatticeBGKFeq(user[i_level].density.at(lat_iter.first), user[i_level].velocity.at(lat_iter.first), feq_temp);

            // D_uvw _velocity = user[i_level].velocity.at(lat_iter.first);
            // D_Phy_Velocity uu = dot_product(_velocity, _velocity);

            for (D_int i_dir = 0; i_dir < C_Q; ++i_dir) {
                // D_Phy_Velocity euv = ex[i_dir] * _velocity.x + ey[i_dir] * _velocity.y + ez[i_dir] * _velocity.z;
                D_Phy_DDF feq_temp = calculateFeq(user[i_level].density.at(lat_iter.first), user[i_level].velocity.at(lat_iter.first), i_dir);
                user[i_level].df.fcol[i_dir].at(lat_iter.first) = feq_temp;
                user[i_level].df.f[i_dir].at(lat_iter.first) = feq_temp;
            }

            // for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            //     user[i_level].df.fcol[i_q].at(lat_iter.first)  = feq_temp[i_q];
            //     user[i_level].df.f[i_q].at(lat_iter.first)  = feq_temp[i_q];
            // } // for i_q
        }
    }
}

void LBM_Manager::initialBC()
{
    bc_manager_->initialBoundaryCondition("West");
     bc_manager_->initialBoundaryCondition("East");
      bc_manager_->initialBoundaryCondition("South");
       bc_manager_->initialBoundaryCondition("North");
        bc_manager_->initialBoundaryCondition("Top");
         bc_manager_->initialBoundaryCondition("Bot");
}

#ifndef BC_NO_GHOST
DEPRECATED void LBM_Manager::allocMacros()
{
    for (auto i_bc_x : Lat_Manager::pointer_me->lat_bc_x) {
        for (auto i_latbc : i_bc_x) {
            user[0].density.insert(make_pair(i_latbc.first, U_rho0/CF_Rho));
            user[0].rho_average.insert(make_pair(i_latbc.first, 0.));
            user[0].velocity.insert(make_pair(i_latbc.first, D_vec(0., 0., 0.)));
        }
    }
    for (auto i_bc_y : Lat_Manager::pointer_me->lat_bc_y) {
        for (auto i_latbc : i_bc_y) {
            user[0].density.insert(make_pair(i_latbc.first, U_rho0/CF_Rho));
            user[0].rho_average.insert(make_pair(i_latbc.first, 0.));
            user[0].velocity.insert(make_pair(i_latbc.first, D_vec(0., 0., 0.)));
        }
    }
    for (auto i_bc_z : Lat_Manager::pointer_me->lat_bc_z) {
        for (auto i_latbc : i_bc_z) {
            user[0].density.insert(make_pair(i_latbc.first, U_rho0/CF_Rho));
            user[0].rho_average.insert(make_pair(i_latbc.first, 0.));
            user[0].velocity.insert(make_pair(i_latbc.first, D_vec(0., 0., 0.)));
        }
    }

    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        for (auto i_lat : *lat_ptr)
        {
            user[i_level].density.insert(make_pair(i_lat.first, U_rho0/CF_Rho));
            user[i_level].rho_average.insert(make_pair(i_lat.first, 0.));
            user[i_level].velocity.insert(make_pair(i_lat.first, D_vec(0., 0., 0.)));
        }
    }
}
#else
DEPRECATED void LBM_Manager::allocMacros()
{
    for (auto i_bc_list : Lat_Manager::pointer_me->lat_bc) {
        D_Phy_Rho _rho0 = U_rho0/CF_Rho;
        D_map_define<D_Phy_Rho> rho0_map;
        D_map_define<D_Phy_Rho> rhoAverage_map;
        D_map_define<D_uvw> zeroUVW_map;
        for (auto bc_lat_code : i_bc_list) {
            rho0_map.insert(make_pair(bc_lat_code, _rho0));
            rhoAverage_map.insert(make_pair(bc_lat_code, 0.));
            zeroUVW_map.insert(make_pair(bc_lat_code, D_vec(0., 0., 0.)));
        }
        user[0].density.insert(rho0_map.begin(), rho0_map.end());
        user[0].rho_average.insert(rhoAverage_map.begin(), rhoAverage_map.end());
        user[0].velocity.insert(zeroUVW_map.begin(), zeroUVW_map.end());
    }

    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        D_Phy_Rho _rho0 = U_rho0/CF_Rho;
        D_map_define<D_Phy_Rho> rho0_map;
        D_map_define<D_Phy_Rho> rhoAverage_map;
        D_map_define<D_uvw> zeroUVW_map;
        for (auto i_lat : *lat_ptr) {
            rho0_map.insert(make_pair(i_lat.first, _rho0));
            rhoAverage_map.insert(make_pair(i_lat.first, 0.));
            zeroUVW_map.insert(make_pair(i_lat.first, D_vec(0., 0., 0.)));
        }
        user[i_level].density.insert(rho0_map.begin(), rho0_map.end());
        user[i_level].rho_average.insert(rhoAverage_map.begin(), rhoAverage_map.end());
        user[i_level].velocity.insert(zeroUVW_map.begin(), zeroUVW_map.end());
    }
}
#endif

#ifndef BC_NO_GHOST
void LBM_Manager::allocMacrosAndDDF()
{
    D_Phy_DDF feq_temp[C_max_level+1][C_Q];
    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level) {
        D_vec v_temp = {U_u0/CF_U, U_v0/CF_U, U_w0/CF_U};
        calculateFeq(U_rho0/CF_Rho, v_temp, feq_temp[i_level]);
    }

    for (auto i_bc_ghost : Lat_Manager::pointer_me->ghost_bc) {
        D_morton lat_code = i_bc_ghost.first;
        user[0].density.insert(make_pair(lat_code, U_rho0/CF_Rho));
        user[0].velocity.insert(make_pair(lat_code, D_vec(U_u0/CF_U, U_v0/CF_U, U_w0/CF_U)));
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[0].df.fcol[i_q].insert(make_pair(lat_code ,feq_temp[0][i_q]));
            user[0].df.f[i_q].insert(make_pair(lat_code ,feq_temp[0][i_q]));
        }
    }

    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        for (auto i_lat : *lat_ptr)
        {
            D_morton lat_code = i_lat.first;

            user[i_level].density.insert(make_pair(lat_code, U_rho0/CF_Rho));
            user[i_level].rho_average.insert(make_pair(lat_code, 0.));
            user[i_level].velocity.insert(make_pair(lat_code, D_vec(U_u0/CF_U, U_v0/CF_U, U_w0/CF_U)));
            user[i_level].v_old.insert(make_pair(lat_code, D_vec(0., 0., 0.)));

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(lat_code, feq_temp[i_level][i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(lat_code, feq_temp[i_level][i_q]));
            }
        }
    }
    
    // SOLID boundary
    for (auto solidLat_iter : Lat_Manager::pointer_me->lat_sf)
    {
        D_morton solid_code = solidLat_iter.first;

        user[C_max_level].density.insert(make_pair(solid_code, 1.));
        user[C_max_level].rho_average.insert(make_pair(solid_code, 0.));
        user[C_max_level].velocity.insert(make_pair(solid_code, D_vec(0., 0., 0.)));
        user[C_max_level].v_old.insert(make_pair(solid_code, D_vec(0., 0., 0.)));

        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[C_max_level].df.f[i_q].insert(make_pair(solid_code, 0.));
            user[C_max_level].df.fcol[i_q].insert(make_pair(solid_code, 0.));
        }
    }

    // overlap
    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        for (auto overlap_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
        {
            D_morton op_code = overlap_lat.first;

            user[i_level].density.insert(make_pair(op_code, U_rho0/CF_Rho));
            user[i_level].rho_average.insert(make_pair(op_code, 0.));
            user[i_level].velocity.insert(make_pair(op_code, D_vec(U_u0/CF_U, U_v0/CF_U, U_w0/CF_U)));
            user[i_level].v_old.insert(make_pair(op_code, D_vec(0., 0., 0.)));

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(op_code, feq_temp[i_level][i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(op_code, feq_temp[i_level][i_q]));
            }
        }

        for (auto overlap_lat : Lat_Manager::pointer_me->lat_overlap_C2F.at(i_level))
        {
            D_morton op_code = overlap_lat.first;

            user[i_level].density.insert(make_pair(op_code, U_rho0/CF_Rho));
            user[i_level].rho_average.insert(make_pair(op_code, 0.));
            user[i_level].velocity.insert(make_pair(op_code, D_vec(U_u0/CF_U, U_v0/CF_U, U_w0/CF_U)));
            user[i_level].v_old.insert(make_pair(op_code, D_vec(0., 0., 0.)));

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(op_code, feq_temp[i_level][i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(op_code, feq_temp[i_level][i_q]));
            }
        }
    }

    // ghost overlap
    for (D_int i_level = 0; i_level <= C_max_level; ++i_level) 
    {
        for (auto i_overlap_ghost : Lat_Manager::pointer_me->ghost_overlap_C2F.at(i_level)) {
            D_morton ghost_code = i_overlap_ghost.first;

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(ghost_code ,feq_temp[i_level][i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(ghost_code ,feq_temp[i_level][i_q]));
            }
        }
    }

}
#else
void LBM_Manager::allocMacrosAndDDF()
{
    D_uvw     init_v    = {U_u0/CF_U, U_v0/CF_U, U_w0/CF_U};
    D_Phy_Rho init_rho0 = U_rho0/CF_Rho;
    D_Phy_DDF init_feq[C_Q];

    Timer tmr;
    calculateFeq(init_rho0, init_v, init_feq);

    for (auto i_bc_list : Lat_Manager::pointer_me->lat_bc) 
    {
        for (auto lat_code : i_bc_list) 
        {
            user[0].density.insert(make_pair(lat_code, init_rho0));
            user[0].rho_average.insert(make_pair(lat_code, 0.));
            user[0].velocity.insert(make_pair(lat_code, init_v));
            user[0].v_old.insert(make_pair(lat_code, D_uvw(0., 0., 0.)));

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[0].df.fcol[i_q].insert(make_pair(lat_code ,init_feq[i_q]));
                user[0].df.f[i_q].insert(make_pair(lat_code ,init_feq[i_q]));
            }
        }
    }

    
    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        double t0 = tmr.elapsed();
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        for (auto i_lat : *lat_ptr)
        {
            D_morton lat_code = i_lat.first;

            user[i_level].density.insert(make_pair(lat_code, init_rho0));
            user[i_level].rho_average.insert(make_pair(lat_code, 0.));
            user[i_level].velocity.insert(make_pair(lat_code, init_v));
            user[i_level].v_old.insert(make_pair(lat_code, D_uvw(0., 0., 0.)));

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(lat_code, init_feq[i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(lat_code, init_feq[i_q]));
            }
        }
        double t1 = tmr.elapsed();
    }
    
    // SOLID boundary
    for (auto solidLat_iter : Lat_Manager::pointer_me->lat_sf)
    {
        D_morton solid_code = solidLat_iter.first;

        user[C_max_level].density.insert(make_pair(solid_code, 1.));
        user[C_max_level].rho_average.insert(make_pair(solid_code, 0.));
        user[C_max_level].velocity.insert(make_pair(solid_code, D_uvw(0., 0., 0.)));
        user[C_max_level].v_old.insert(make_pair(solid_code, D_uvw(0., 0., 0.)));

        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[C_max_level].df.f[i_q].insert(make_pair(solid_code, 0.));
            user[C_max_level].df.fcol[i_q].insert(make_pair(solid_code, 0.));
        }
    }

    // overlap
    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        for (auto overlap_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
        {
            D_morton op_code = overlap_lat.first;

            user[i_level].density.insert(make_pair(op_code, init_rho0));
            user[i_level].rho_average.insert(make_pair(op_code, 0.));
            user[i_level].velocity.insert(make_pair(op_code, init_v));
            user[i_level].v_old.insert(make_pair(op_code, D_uvw(0., 0., 0.)));

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(op_code, init_feq[i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(op_code, init_feq[i_q]));
            }
        }

        for (auto overlap_lat : Lat_Manager::pointer_me->lat_overlap_C2F.at(i_level))
        {
            D_morton op_code = overlap_lat.first;

            user[i_level].density.insert(make_pair(op_code, init_rho0));
            user[i_level].rho_average.insert(make_pair(op_code, 0.));
            user[i_level].velocity.insert(make_pair(op_code, init_v));
            user[i_level].v_old.insert(make_pair(op_code, D_uvw(0., 0., 0.)));

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(op_code, init_feq[i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(op_code, init_feq[i_q]));
            }
        }
    }

    // ghost overlap
    for (D_int i_level = 0; i_level <= C_max_level; ++i_level) 
    {
        for (auto i_overlap_ghost : Lat_Manager::pointer_me->ghost_overlap_C2F.at(i_level)) {
            D_morton ghost_code = i_overlap_ghost.first;

            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(ghost_code ,init_feq[i_q]));
                user[i_level].df.fcol[i_q].insert(make_pair(ghost_code ,init_feq[i_q]));
            }
        }
    }

}
#endif

#ifndef BC_NO_GHOST
DEPRECATED void LBM_Manager::allocDDF()
{
    for (auto i_bc_x : Lat_Manager::pointer_me->lat_bc_x) {
        for (auto i_latbc : i_bc_x) {
            D_morton lat_code = i_latbc.first;
            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[0].df.fcol[i_q].insert(make_pair(lat_code ,0.));
                user[0].df.f[i_q].insert(make_pair(lat_code ,0.));
            }
        }
    }
    for (auto i_bc_y : Lat_Manager::pointer_me->lat_bc_y) {
        for (auto i_latbc : i_bc_y) {
            D_morton lat_code = i_latbc.first;
            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[0].df.fcol[i_q].insert(make_pair(lat_code ,0.));
                user[0].df.f[i_q].insert(make_pair(lat_code ,0.));
            }
        }
    }
    for (auto i_bc_z : Lat_Manager::pointer_me->lat_bc_z) {
        for (auto i_latbc : i_bc_z) {
            D_morton lat_code = i_latbc.first;
            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[0].df.fcol[i_q].insert(make_pair(lat_code ,0.));
                user[0].df.f[i_q].insert(make_pair(lat_code ,0.));
            }
        }
    }

    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        // D_int num_lattice = Lat_Manager::pointer_me->lat_f[i_level].size();

        // Structure of Array ↑↑↑→→→↓↓↓←←←
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        for (auto i_lat : *lat_ptr)
        {
            D_morton lat_code = i_lat.first;
            for (D_int i_q = 0; i_q < C_Q; ++i_q) {
                user[i_level].df.f[i_q].insert(make_pair(lat_code ,0.));
                user[i_level].df.fcol[i_q].insert(make_pair(lat_code ,0.));
            }
        }
    }

    // not real lattice, only containing ddf for stream, not executing collison & fusion
    for (auto i_lat : Lat_Manager::pointer_me->lat_sf)
    {
        D_morton lat_code = i_lat.first;
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[C_max_level].df.f[i_q].insert(make_pair(lat_code ,0.));
            user[C_max_level].df.fcol[i_q].insert(make_pair(lat_code ,0.));
        }
    }
    
    int lat_f_fluid = 0, lat_f_ib = 0, lat_f_F2C = 0;

    for (auto i_lat : Lat_Manager::pointer_me->lat_f[C_max_level]) {
        if (i_lat.second.flag == FLUID) lat_f_fluid++;
        else if (i_lat.second.flag == IB) lat_f_ib++;
        else if (i_lat.second.flag == OVERLAP_F2C) lat_f_F2C++;
        else {
            std::cout << "Flag " << i_lat.second.flag << std::endl;
        }
    }

    std::cout << "lat_f_fluid " << lat_f_fluid;
    std::cout << " lat_f_ib " << lat_f_ib;
    std::cout << " lat_f_F2C " << lat_f_F2C << std::endl;
    std::cout << " lat_overlap_F2C " << Lat_Manager::pointer_me->lat_overlap_F2C.at(C_max_level-1).size() << std::endl;

    int lat_sf_ib = 0, lat_sf_surface = 0;
    for (auto i_latSf : Lat_Manager::pointer_me->lat_sf) {
        if (i_latSf.second.flag == IB) lat_sf_ib++;
        else if (i_latSf.second.flag == SURFACE || i_latSf.second.flag == SOLID) lat_sf_surface++;
    }
    std::cout << "lat_sf_ib " << lat_sf_ib;
    std::cout << " lat_sf_surface " << lat_sf_surface << std::endl;
}
#else
DEPRECATED void LBM_Manager::allocDDF()
{

    for (auto i_bc_list : Lat_Manager::pointer_me->lat_bc) {
        D_map_define<D_Phy_DDF> zero_ddf;
        for (auto bc_lat_code : i_bc_list) {
            zero_ddf.insert(make_pair(bc_lat_code, 0.));
        }
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[0].df.fcol[i_q].insert(zero_ddf.begin(),zero_ddf.end());
            user[0].df.f[i_q].insert(zero_ddf.begin(),zero_ddf.end());
        }
    }

    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        // D_int num_lattice = Lat_Manager::pointer_me->lat_f[i_level].size();

        // Structure of Array ↑↑↑→→→↓↓↓←←←
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        D_map_define<D_Phy_DDF> zero_ddf;
        for (auto i_lat : *lat_ptr) {
            zero_ddf.insert(make_pair(i_lat.first, 0.));
        }
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[i_level].df.f[i_q].insert(zero_ddf.begin(),zero_ddf.end());
            user[i_level].df.fcol[i_q].insert(zero_ddf.begin(),zero_ddf.end());
        }
    }

    // not real lattice, only containing ddf for stream, not executing collison & fusion
    D_map_define<D_Phy_DDF> zero_ddf;
    for (auto i_lat : Lat_Manager::pointer_me->lat_sf) {
        zero_ddf.insert(make_pair(i_lat.first, 0.));
    }
    for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[C_max_level].df.f[i_q].insert(zero_ddf.begin(),zero_ddf.end());
            user[C_max_level].df.fcol[i_q].insert(zero_ddf.begin(),zero_ddf.end());
    }
    
    int lat_f_fluid = 0, lat_f_ib = 0, lat_f_F2C = 0;

    for (auto i_lat : Lat_Manager::pointer_me->lat_f[C_max_level]) {
        if (i_lat.second.flag == FLUID) lat_f_fluid++;
        else if (i_lat.second.flag == IB) lat_f_ib++;
        else if (i_lat.second.flag == OVERLAP_F2C) lat_f_F2C++;
        else {
            std::cout << "Flag " << i_lat.second.flag << std::endl;
        }
    }

    std::cout << "lat_f_fluid " << lat_f_fluid;
    std::cout << " lat_f_ib " << lat_f_ib;
    std::cout << " lat_f_F2C " << lat_f_F2C << std::endl;
    std::cout << " lat_overlap_F2C " << Lat_Manager::pointer_me->lat_overlap_F2C.at(C_max_level-1).size() << std::endl;

    int lat_sf_ib = 0, lat_sf_surface = 0;
    for (auto i_latSf : Lat_Manager::pointer_me->lat_sf) {
        if (i_latSf.second.flag == IB) lat_sf_ib++;
        else if (i_latSf.second.flag == SURFACE || i_latSf.second.flag == SOLID) lat_sf_surface++;
    }
    std::cout << "lat_sf_ib " << lat_sf_ib;
    std::cout << " lat_sf_surface " << lat_sf_surface << std::endl;
}
#endif

void LBM_Manager::checkMach(D_Phy_Velocity u0)
{
    D_real Ma = u0 / L_Cs;
    if (Ma >= MAX_MACH) {
        std::cout << " - the MACH="<<Ma<<" is larger than "<<MAX_MACH<<", the simulation may be divergence" << std::endl;
        std::cout << " - the dt will be reset by MAX_MACH="<<MAX_MACH << std::endl;
        D_real refVel = u0 * CF_U;
        CF_U = refVel / (0.9 * MAX_MACH * L_Cs);
        for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
        {
            CF_T[i_level] = CF_L[i_level] / CF_U;
            CF_Nu[i_level] = CF_L[i_level] * CF_U;
        }
        std::cout << " - the dt is changed: " << U_dt <<" --> " << CF_T[0] << std::endl;
        U_dt = CF_T[0];
    } else {
        std::cout << " - the MACH="<<Ma<<" is satisfied the LBM convergence condition" << std::endl;
    }
}

void LBM_Manager::setUnit(std::string collision_Model)
{
    printf("------------ LBM unit conversion ------------\n");
    std::cout << " - grid spacing = " << C_dx << std::endl;
    std::cout << " - Physical Kinematic Viscosity = "<< U_kineViscosity <<" [m2/s]" << std::endl;
    std::cout << " - Physical Reference Velocity  = "<< R_velocity <<" [m/s]" << std::endl;

    D_real L_SIZE[C_max_level+1], L_TIME[C_max_level+1];
    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        L_SIZE[i_level] = 1. / two_power_n(i_level);
        L_TIME[i_level] = 1. / two_power_n(i_level);
    }

    /**
     * @brief Compute the dt range by the given space and physical parameter of fluid
     * @author ZouHang
     * @date 2023-4-24
     */
    if (collision_Model == "LBGK")
    {
        std::cout << " - MIN_TAU = " << MIN_TAU << ", MAX_MACH = " << MAX_MACH << std::endl;
        D_real max_dh = U_kineViscosity * MAX_MACH / (R_velocity * L_Cs *(MIN_TAU - 0.5));
        std::cout << " - the max space is:              " << max_dh << std::endl;
        for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
        {
            D_real dx = Lat_Manager::pointer_me->dx.at(i_level);
            D_real dt = U_dt * L_TIME[i_level];
            D_real max_dt = dx * MAX_MACH * L_Cs / R_velocity;
            D_real min_dt = SQ(dx) * (MIN_TAU - 0.5) * L_CsSq / U_kineViscosity;
            std::cout << " - [Level"<<i_level<<"] time step dt should be in:    ("<<min_dt<<", "<<max_dt<<")" << std::endl;
            if (dt < min_dt || dt > max_dt) {
                std::cout << " - your dt is out of the recommend range" << std::endl;
                std::cout << " - the dt is changed: " << U_dt <<" --> " << max_dt * 0.99 / L_TIME[i_level] << std::endl;
                U_dt = max_dt * 0.99 / L_TIME[i_level];
            }
        }
    }

    // compute the conversion factor by dh and dt
    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level) {
        CF_L[i_level] = C_dx * L_SIZE[i_level];
        CF_T[i_level] = U_dt * L_TIME[i_level];
    }
    CF_U = CF_L[0] / CF_T[0];
    CF_Rho = U_rho0 / L_Rho;
    checkMach(R_velocity/CF_U);
    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level) {
        tau[i_level] = U_kineViscosity / CF_Nu[i_level] / (L_CsSq * L_TIME[i_level]) + 0.5;
        std::cout << " - [Level"<<i_level<<"] relaxation time = " << tau[i_level] << std::endl;
    }
    std::cout << " - Re = " << R_Re << std::endl;

#ifdef DDFTRANSFACTOR
        if (i_level != C_max_level) { 
            transFactor_C2F[i_level][0] = 1./(2.*tau[i_level]) * (2.*tau[i_level] + 1./two_power_n(i_level) - 1.);
            transFactor_C2F[i_level][1] = 1./(2.*tau[i_level]) * (1./two_power_n(i_level) - 1.);
        }
        if (i_level != 0) {
            transFactor_F2C[i_level][0] = 1./(2.*tau[i_level]) * (2.*tau[i_level] + two_power_n(i_level) - 1.);
            transFactor_F2C[i_level][1] = -1./(2.*tau[i_level]) * (two_power_n(i_level) - 1.);
        }
#endif

    printf("------------ -------------------- ------------\n");
}

void LBM_Manager::infoPrint(D_int iter_no)
{
    std::cout << "[Iter] " << iter_no;
    D_real temp1, temp2;
    D_int i, j, k, id;
    D_real guerror = 0., grerr = 0., uerror = 0.;

    for (D_int i_level = 0; i_level < C_max_level+1; ++i_level)
    {
        D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
        for (auto i_lat : *lat_ptr) {
            D_morton lat_code = i_lat.first;

            if (isnanf(user[i_level].velocity.at(lat_code).x) || isnanf(user[i_level].velocity.at(lat_code).y) || isnanf(user[i_level].velocity.at(lat_code).z))
            {
                std::cout << "lat_code " << lat_code << " V=(" << user[i_level].velocity.at(lat_code).x << "," <<
                    user[i_level].velocity.at(lat_code).y << "," <<
                    user[i_level].velocity.at(lat_code).z << ")" << std::endl;
            }

            guerror += (user[i_level].velocity.at(lat_code).x - user[i_level].v_old.at(lat_code).x) * (user[i_level].velocity.at(lat_code).x - user[i_level].v_old.at(lat_code).x) + (user[i_level].velocity.at(lat_code).y - user[i_level].v_old.at(lat_code).y) * (user[i_level].velocity.at(lat_code).y - user[i_level].v_old.at(lat_code).y) + (user[i_level].velocity.at(lat_code).z - user[i_level].v_old.at(lat_code).z) * (user[i_level].velocity.at(lat_code).z - user[i_level].v_old.at(lat_code).z);
            grerr += (user[i_level].velocity.at(lat_code).x * user[i_level].velocity.at(lat_code).x + user[i_level].velocity.at(lat_code).y * user[i_level].velocity.at(lat_code).y + user[i_level].velocity.at(lat_code).z * user[i_level].velocity.at(lat_code).z);
        }
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level)) {
            D_morton lat_code = i_lat.first;

            if (isnanf(user[i_level].velocity.at(lat_code).x) || isnanf(user[i_level].velocity.at(lat_code).y) || isnanf(user[i_level].velocity.at(lat_code).z))
            {
                std::cout << "lat_code " << lat_code << " V=(" << user[i_level].velocity.at(lat_code).x << "," <<
                    user[i_level].velocity.at(lat_code).y << "," <<
                    user[i_level].velocity.at(lat_code).z << ")" << std::endl;
            }
        }
    }

    uerror = sqrt(guerror) / sqrt(grerr + 1e-16);

    std::cout << " time = " << iter_no * CF_T[0] << " Velocity error = " << uerror << std::endl;
}

void LBM_Manager::calculateMacros(const D_Phy_DDF* ddf, D_Phy_Rho& rho, D_uvw& v)
{
    rho = ddf[0];
    v = {0.,0.,0.};
    for (D_int i_q = 1; i_q < C_Q; ++i_q)
    {
        rho += ddf[i_q];
        v.x += ddf[i_q] * ex[i_q];
        v.y += ddf[i_q] * ey[i_q];
        v.z += ddf[i_q] * ez[i_q];
    }
    v = v / rho;
}

void LBM_Manager::calculateMacros(D_int level)
{
    D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(level));
    D_int nan_count = 0;

    D_Phy_DDF ddf[C_Q];

    for (auto i_lat : *lat_ptr) {
        D_morton lat_code = i_lat.first;
        load_ddf(ddf, level, lat_code);
        calculateMacros(ddf, user[level].density.at(lat_code), user[level].velocity.at(lat_code));
        if (isnanf(user[level].velocity.at(lat_code).x) || isnanf(user[level].velocity.at(lat_code).y) || isnanf(user[level].velocity.at(lat_code).z))
        {
            // 0000000000000000000000000000000000001010101001001110000110100001
            std::cout << "level " << level << std::endl;
            std::cout << "lat_code " << lat_code << " V=(" << user[level].velocity.at(lat_code).x << "," <<
                user[level].velocity.at(lat_code).y << "," <<
                user[level].velocity.at(lat_code).z << ")" << std::endl;
            for (auto iq=0; iq<C_Q; ++iq) {
                std::cout << "  iq " << iq << " ddf " << ddf[iq] << std::endl;
            }
        }
    } // for i_lat

    for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(level)) {
        D_morton lat_code = i_lat.first;
        load_ddf(ddf, level, lat_code);
        calculateMacros(ddf, user[level].density.at(lat_code), user[level].velocity.at(lat_code));
    } // for i_lat

}

void LBM_Manager::calculateMacros()
{
    for (D_int level = 0; level <= C_max_level; ++level) {
        calculateMacros(level);
    }
}

void LBM_Manager::AMR_transDDF_explosion(D_int level)
{
    D_int fine_level = level;
    D_int coarse_level = level - 1;

    // D_morton xxx_code("0000000000000000000000000000000000001010101001001110000110100001");
    // std::cout << "Before explosion [level " << level << "] " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;

    // for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(level)) 
    for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(coarse_level)) 
    {
        D_morton fine_code[8]{i_lat.first,
                               Morton_Assist::find_x1(i_lat.first, fine_level),
                                Morton_Assist::find_y1(i_lat.first, fine_level),
                                 Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
                                  Morton_Assist::find_z1(i_lat.first, fine_level),
                                   Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
                                    Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, fine_level), fine_level),
                                     Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level), fine_level)
        };

        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            for (D_int i_fineLat = 0; i_fineLat < 8; ++i_fineLat) {
                user[fine_level].df.fcol[i_q].at(fine_code[i_fineLat]) = user[coarse_level].df.fcol[i_q].at(fine_code[0]);
            }
        }
    }

    // std::cout << "After explosion lat_F2C [level " << level << "] " << user[4].df.fcol[0].at(xxx_code) << " " << user[4].df.fcol[1].at(xxx_code) << std::endl;

    for (auto i_lat : Lat_Manager::pointer_me->ghost_overlap_F2C.at(coarse_level)) 
    {
        D_morton ngbr_code[8]{i_lat.first, 
                               Morton_Assist::find_x1(i_lat.first, fine_level), 
                                Morton_Assist::find_y1(i_lat.first, fine_level),
                                 Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
                                  Morton_Assist::find_z1(i_lat.first, fine_level),
                                   Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
                                    Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, fine_level), fine_level),
                                     Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level), fine_level)
                                     };

        for (D_int i_ngbr = 0; i_ngbr < 8; ++i_ngbr) {
            if (Lat_Manager::pointer_me->ghost_overlap_C2F.at(fine_level).find(ngbr_code[i_ngbr]) != Lat_Manager::pointer_me->ghost_overlap_C2F.at(fine_level).end()) {
                for (D_int i_q = 0; i_q < C_Q; ++i_q)
                    user[fine_level].df.fcol[i_q].at(ngbr_code[i_ngbr]) = user[coarse_level].df.fcol[i_q].at(ngbr_code[0]);
            }
        }

    }

    // std::cout << "After explosion ghost_F2C [level " << level << "] " << user[4].df.fcol[0].at(xxx_code) << " " << user[4].df.fcol[1].at(xxx_code) << std::endl;

}

void LBM_Manager::AMR_transDDF_coalescence(D_int level)
{
    D_int fine_level = level + 1;
    D_int coarse_level = level;
    for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(coarse_level))
    {
        D_morton ngbr_code[8]{i_lat.first, 
                               Morton_Assist::find_x1(i_lat.first, fine_level), 
                                Morton_Assist::find_y1(i_lat.first, fine_level),
                                 Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
                                  Morton_Assist::find_z1(i_lat.first, fine_level),
                                   Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level),
                                    Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, fine_level), fine_level),
                                     Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, fine_level), fine_level), fine_level)};

        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            user[coarse_level].df.f[i_q].at(i_lat.first) = 0;
            for (D_int i_ngbr = 0; i_ngbr < 8; ++i_ngbr) {
                user[coarse_level].df.f[i_q].at(i_lat.first) += user[fine_level].df.f[i_q].at(ngbr_code[i_ngbr]);
            }
            user[coarse_level].df.f[i_q].at(i_lat.first) = user[coarse_level].df.f[i_q].at(i_lat.first) / 8;
        }

    }
}

void LBM_Manager::load_ddf(D_Phy_DDF* ddf, D_int level, D_morton code)
{
    for (D_int i_q = 0; i_q < C_Q; ++i_q) {
        ddf[i_q] = user[level].df.f[i_q].at(code);
    }
}

void LBM_Manager::store_ddf(const D_Phy_DDF* ddf, D_int level, D_morton code)
{
    for (D_int i_q = 0; i_q < C_Q; ++i_q) {
        user[level].df.f[i_q].at(code) = ddf[i_q];
    }
}