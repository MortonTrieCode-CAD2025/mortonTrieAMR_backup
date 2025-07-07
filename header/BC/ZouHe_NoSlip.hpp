// /**
//  * @file ZouHe_NoSlip.h
//  * @brief Non-equilibrum Bounce-Back boundary condition applying in non-slip wall
//  *        V = 0 for ZouHe_Velocity (suitable for Low Re Number)
//  * @date 2023-07-18
//  * @ref M. Hecht and J. Harting, “Implementation of on-site velocity boundary conditions for d3q19 lattice boltzmann simulations,” Journal of Statistical Mechanics: Theory and Experiment, vol. 2010, no. 01, 2010. doi:10.1088/1742-5468/2010/01/p01018 
//  *      M. Hecht and J. Harting, “Erratum: Implementation of on-site velocity boundary conditions for d3q19 lattice boltzmann simulations,” Journal of Statistical Mechanics: Theory and Experiment, vol. 2013, no. 02, 2013. doi:10.1088/1742-5468/2013/02/e02001 
//  */

// #pragma once
// #include "user.h"
// #include "D3Q19_BC_Manager.h"

// class nonEquilibrumBounceBack_NoSlip_BC : public D3Q19_BC_Strategy
// {
// public:
//     virtual void applyBCStrategy() = 0;
//     virtual void initialBCStrategy() = 0;
// };

// class nonEquilibrumBounceBack_NoSlip_West : public nonEquilibrumBounceBack_NoSlip_BC
// {
// private:

// public:
//     nonEquilibrumBounceBack_NoSlip_West() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;
// };

// class nonEquilibrumBounceBack_NoSlip_East : public nonEquilibrumBounceBack_NoSlip_BC
// {
// private:

// public:
//     nonEquilibrumBounceBack_NoSlip_East() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;
// };

// class nonEquilibrumBounceBack_NoSlip_North : public nonEquilibrumBounceBack_NoSlip_BC
// {
// private:

// public:
//     nonEquilibrumBounceBack_NoSlip_North() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;
// };

// class nonEquilibrumBounceBack_NoSlip_South : public nonEquilibrumBounceBack_NoSlip_BC
// {
// private:

// public:
//     nonEquilibrumBounceBack_NoSlip_South() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;
// };

// class nonEquilibrumBounceBack_NoSlip_Top : public nonEquilibrumBounceBack_NoSlip_BC
// {
// private:

// public:
//     nonEquilibrumBounceBack_NoSlip_Top() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;
// };

// class nonEquilibrumBounceBack_NoSlip_Bot : public nonEquilibrumBounceBack_NoSlip_BC
// {
// private:

// public:
//     nonEquilibrumBounceBack_NoSlip_Bot() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;
// };

// // class nonEquilibrumBounceBack_NoSlip_Edge : public nonEquilibrumBounceBack_NoSlip_BC
// // {
// // private:
// //     const std::initializer_list<std::string> edge_list;
// //     bool all_edge_same = false;
// //     void applyOnEdgeBC_NE();
// //     void applyOnEdgeBC_NW();
// //     void applyOnEdgeBC_SE();
// //     void applyOnEdgeBC_SW();

// //     void applyOnEdgeBC_TN();
// //     void applyOnEdgeBC_TS();
// //     void applyOnEdgeBC_BN();
// //     void applyOnEdgeBC_BS();

// //     void applyOnEdgeBC_TE();
// //     void applyOnEdgeBC_TW();
// //     void applyOnEdgeBC_BE();
// //     void applyOnEdgeBC_BW();
    
// // public:
// //     nonEquilibrumBounceBack_NoSlip_Edge(const std::initializer_list<std::string>& edge_name) : edge_list(edge_name) 
// //     {
// //         const std::set<std::string> predefined_edgeName = {"NorthEast", "Northeast", "northeast",
// //                                                            "NorthWest", "Northwest", "northwest",
// //                                                            "SouthWest", "Southwest", "southwest",
// //                                                            "SouthEast", "Southeast", "southeast",
// //                                                            "TopNorth", "Topnorth", "topnorth",
// //                                                            "TopSouth", "Topsouth", "topsouth",
// //                                                            "BotNorth", "Botnorth", "botnorth",
// //                                                            "BotSouth", "Botsouth", "botsouth",
// //                                                            "TopEast", "Topeast", "topeast",
// //                                                            "TopWest", "Topwest", "topwest",
// //                                                            "BotEast", "Boteast", "boteast",
// //                                                            "BotWest", "Botwest", "botwest"};

// //         if (std::set<std::string>(edge_list).count("All") > 0U || std::set<std::string>(edge_list).count("all") > 0U) {
// //             all_edge_same = true;
// //         }

// //         for (const auto& name : edge_list) {
// //             if (predefined_edgeName.find(name) == predefined_edgeName.end()) {
// //                 std::stringstream error_log;
// //                 error_log << name << " not contains in built-in BC for No-slip edge treatment" << std::endl;
// //                 log_error(error_log.str(), Log_function::logfile);
// //             }
// //         }

// //     };

// //     void applyBCStrategy() override;
// // };



// // class nonEquilibrumBounceBack_NoSlip_Corner : public nonEquilibrumBounceBack_NoSlip_BC
// // {
// // private:
// //     const std::initializer_list<std::string> corner_list;
// //     bool all_corner_same = false;
    
// // public:
// //     nonEquilibrumBounceBack_NoSlip_Corner(const std::initializer_list<std::string>& corner_name) : corner_list(corner_name) 
// //     {
// //         const std::set<std::string> predefined_cornerName = {"NorthEastTop", "Northeasttop", "northeasttop",
// //                                                              "NorthWestTop", "Northwesttop", "northwesttop",
// //                                                              "SouthWestTop", "Southwesttop", "southwesttop",
// //                                                              "SouthEastTop", "Southeasttop", "southeasttop",
// //                                                              "NorthEastBot", "Northeastbot", "northeastbot",
// //                                                              "NorthWestBot", "Northwestbot", "northwestbot",
// //                                                              "SouthWestBot", "Southwestbot", "southwestbot",
// //                                                              "SouthEastBot", "Southeastbot", "southeastbot"};

// //         if (std::set<std::string>(corner_list).count("All") > 0U || std::set<std::string>(corner_list).count("all") > 0U) {
// //             all_corner_same = true;
// //         }

// //         for (const auto& name : corner_list) {
// //             if (predefined_cornerName.find(name) == predefined_cornerName.end()) {
// //                 std::stringstream error_log;
// //                 error_log << name << " not contains in built-in BC for No-slip corner treatment" << std::endl;
// //                 log_error(error_log.str(), Log_function::logfile);
// //             }
// //         }

// //     };

    

// // };