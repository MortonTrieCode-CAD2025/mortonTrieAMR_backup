// /**
//  * @file ZouHe_NoSlip.cpp
//  * @brief Non-equilibrum Bounce-Back boundary condition applying in non-slip wall
//  *        V = 0 for ZouHe_Velocity (suitable for Low Re Number)
//  * @ref M. Hecht and J. Harting, “Implementation of on-site velocity boundary conditions for d3q19 lattice boltzmann simulations,” Journal of Statistical Mechanics: Theory and Experiment, vol. 2010, no. 01, 2010. doi:10.1088/1742-5468/2010/01/p01018 
//  *      M. Hecht and J. Harting, “Erratum: Implementation of on-site velocity boundary conditions for d3q19 lattice boltzmann simulations,” Journal of Statistical Mechanics: Theory and Experiment, vol. 2013, no. 02, 2013. doi:10.1088/1742-5468/2013/02/e02001 
//  * @date 2023-07-20
//  */

// #include "BC/ZouHe_NoSlip.h"
// #include "Lat_Manager.h"

// void nonEquilibrumBounceBack_NoSlip_West::initialBCStrategy()
// {
//     D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
//     initial_Feq(westBc_ptr);
// }

// void nonEquilibrumBounceBack_NoSlip_West::applyBCStrategy()
// {
//     D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
//     for (auto i_bc_lat : *westBc_ptr)
//     {
//         D_morton lat_code = i_bc_lat.first;
//         D_Phy_DDF f0 = f_ptr[0].at(lat_code);
//         D_Phy_DDF f2 = f_ptr[2].at(lat_code);
//         D_Phy_DDF f3 = f_ptr[3].at(lat_code);
//         D_Phy_DDF f4 = f_ptr[4].at(lat_code);
//         D_Phy_DDF f5 = f_ptr[5].at(lat_code);
//         D_Phy_DDF f6 = f_ptr[6].at(lat_code);
//         D_Phy_DDF f8 = f_ptr[8].at(lat_code);
//         D_Phy_DDF f10= f_ptr[10].at(lat_code);
//         D_Phy_DDF f12= f_ptr[12].at(lat_code);
//         D_Phy_DDF f14= f_ptr[14].at(lat_code);
//         D_Phy_DDF f15= f_ptr[15].at(lat_code);
//         D_Phy_DDF f16= f_ptr[16].at(lat_code);
//         D_Phy_DDF f17= f_ptr[17].at(lat_code);
//         D_Phy_DDF f18= f_ptr[18].at(lat_code);
//         D_Phy_DDF ny_x = 1./2. * (f3+f15+f17-f4-f18-f16);
//         D_Phy_DDF nz_x = 1./2. * (f5+f15+f18-f6-f17-f16);
//         density_ptr->at(lat_code) = f0+f3+f4+f5+f6+f15+f17+f18+f16 + 2.*(f2+f10+f8+f14+f12);
//         f_ptr[1].at(lat_code) = f2;
//         f_ptr[7].at(lat_code) = f8 - ny_x;
//         f_ptr[9].at(lat_code) = f10 + ny_x;
//         f_ptr[11].at(lat_code) = f12 - nz_x;
//         f_ptr[13].at(lat_code) = f14 + nz_x;
//     }
// }

// void nonEquilibrumBounceBack_NoSlip_East::initialBCStrategy()
// {
//     D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
//     initial_Feq(eastBc_ptr);
// }

// void nonEquilibrumBounceBack_NoSlip_East::applyBCStrategy()
// {
//     D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
//     for (auto i_bc_lat : *eastBc_ptr)
//     {
//         D_morton lat_code = i_bc_lat.first;
//         D_Phy_DDF f0 = f_ptr[0].at(lat_code);
//         D_Phy_DDF f1 = f_ptr[1].at(lat_code);
//         D_Phy_DDF f3 = f_ptr[3].at(lat_code);
//         D_Phy_DDF f4 = f_ptr[4].at(lat_code);
//         D_Phy_DDF f5 = f_ptr[5].at(lat_code);
//         D_Phy_DDF f6 = f_ptr[6].at(lat_code);
//         D_Phy_DDF f7 = f_ptr[7].at(lat_code);
//         D_Phy_DDF f9 = f_ptr[9].at(lat_code);
//         D_Phy_DDF f11= f_ptr[11].at(lat_code);
//         D_Phy_DDF f13= f_ptr[13].at(lat_code);
//         D_Phy_DDF f15= f_ptr[15].at(lat_code);
//         D_Phy_DDF f16= f_ptr[16].at(lat_code);
//         D_Phy_DDF f17= f_ptr[17].at(lat_code);
//         D_Phy_DDF f18= f_ptr[18].at(lat_code);
//         D_Phy_DDF ny_x = 1./2. * (f3+f15+f17-f4-f18-f16);
//         D_Phy_DDF nz_x = 1./2. * (f5+f15+f18-f6-f17-f16);
//         density_ptr->at(lat_code) = f0+f3+f4+f5+f6+f15+f17+f18+f16 + 2.*(f1+f7+f9+f11+f13);
//         f_ptr[2].at(lat_code) = f1;
//         f_ptr[8].at(lat_code) = f7 + ny_x;
//         f_ptr[10].at(lat_code) = f9 - ny_x;
//         f_ptr[12].at(lat_code) = f11 + nz_x;
//         f_ptr[14].at(lat_code) = f13 - nz_x;
//     }
// }

// void nonEquilibrumBounceBack_NoSlip_South::initialBCStrategy()
// {
//     D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
//     initial_Feq(southBc_ptr);
// }

// void nonEquilibrumBounceBack_NoSlip_South::applyBCStrategy()
// {
//     D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
//     for (auto i_bc_lat : *southBc_ptr)
//     {
//         D_morton lat_code = i_bc_lat.first;
//         D_Phy_DDF f0 = f_ptr[0].at(lat_code);
//         D_Phy_DDF f1 = f_ptr[1].at(lat_code);
//         D_Phy_DDF f2 = f_ptr[2].at(lat_code);
//         D_Phy_DDF f4 = f_ptr[4].at(lat_code);
//         D_Phy_DDF f5 = f_ptr[5].at(lat_code);
//         D_Phy_DDF f6 = f_ptr[6].at(lat_code);
//         D_Phy_DDF f8 = f_ptr[8].at(lat_code);
//         D_Phy_DDF f9 = f_ptr[9].at(lat_code);
//         D_Phy_DDF f11= f_ptr[11].at(lat_code);
//         D_Phy_DDF f12= f_ptr[12].at(lat_code);
//         D_Phy_DDF f13= f_ptr[13].at(lat_code);
//         D_Phy_DDF f14= f_ptr[14].at(lat_code);
//         D_Phy_DDF f16= f_ptr[16].at(lat_code);
//         D_Phy_DDF f18= f_ptr[18].at(lat_code);
//         D_Phy_DDF nx_y = 1./2. * (f1+f11+f13-f2-f14-f12);
//         D_Phy_DDF nz_y = 1./2. * (f5+f11+f14-f6-f13-f12);
//         density_ptr->at(lat_code) = f0+f1+f2+f5+f6+f11+f13+f14+f12 + 2.*(f4+f9+f8+f18+f16);
//         f_ptr[3].at(lat_code) = f4;
//         f_ptr[7].at(lat_code) = f8 - nx_y;
//         f_ptr[10].at(lat_code) = f9 + nx_y;
//         f_ptr[15].at(lat_code) = f16 - nz_y;
//         f_ptr[17].at(lat_code) = f18 + nz_y;
//     }


// }

// void nonEquilibrumBounceBack_NoSlip_North::initialBCStrategy()
// {
//     D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
//     initial_Feq(northBc_ptr);
// }

// void nonEquilibrumBounceBack_NoSlip_North::applyBCStrategy()
// {
//     D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
//     for (auto i_bc_lat : *northBc_ptr)
//     {
//         D_morton lat_code = i_bc_lat.first;
//         D_Phy_DDF f0 = f_ptr[0].at(lat_code);
//         D_Phy_DDF f1 = f_ptr[1].at(lat_code);
//         D_Phy_DDF f2 = f_ptr[2].at(lat_code);
//         D_Phy_DDF f3 = f_ptr[3].at(lat_code);
//         D_Phy_DDF f5 = f_ptr[5].at(lat_code);
//         D_Phy_DDF f6 = f_ptr[6].at(lat_code);
//         D_Phy_DDF f7 = f_ptr[7].at(lat_code);
//         D_Phy_DDF f10= f_ptr[10].at(lat_code);
//         D_Phy_DDF f11= f_ptr[11].at(lat_code);
//         D_Phy_DDF f12= f_ptr[12].at(lat_code);
//         D_Phy_DDF f13= f_ptr[13].at(lat_code);
//         D_Phy_DDF f14= f_ptr[14].at(lat_code);
//         D_Phy_DDF f15= f_ptr[15].at(lat_code);
//         D_Phy_DDF f17= f_ptr[17].at(lat_code);
//         D_Phy_DDF nx_y = 1./2. * (f1+f11+f13-f2-f14-f12);
//         D_Phy_DDF nz_y = 1./2. * (f5+f11+f14-f6-f13-f12);
//         density_ptr->at(lat_code) = f0+f1+f2+f5+f6+f11+f13+f14+f12 + 2.*(f3+f7+f10+f15+f17);
//         f_ptr[4].at(lat_code) = f3;
//         f_ptr[8].at(lat_code) = f7 + nx_y;
//         f_ptr[9].at(lat_code) = f10 - nx_y;
//         f_ptr[16].at(lat_code) = f15 + nz_y;
//         f_ptr[18].at(lat_code) = f17 - nz_y;
//     }
// }

// void nonEquilibrumBounceBack_NoSlip_Top::initialBCStrategy()
// {
//     D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
//     initial_Feq(topBc_ptr);
// }

// void nonEquilibrumBounceBack_NoSlip_Top::applyBCStrategy()
// {
//     D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
//     for (auto i_bc_lat : *topBc_ptr)
//     {
//         D_morton lat_code = i_bc_lat.first;
//         D_Phy_DDF f0 = f_ptr[0].at(lat_code);
//         D_Phy_DDF f1 = f_ptr[1].at(lat_code);
//         D_Phy_DDF f2 = f_ptr[2].at(lat_code);
//         D_Phy_DDF f3 = f_ptr[3].at(lat_code);
//         D_Phy_DDF f4 = f_ptr[4].at(lat_code);
//         D_Phy_DDF f5 = f_ptr[5].at(lat_code);
//         D_Phy_DDF f7 = f_ptr[7].at(lat_code);
//         D_Phy_DDF f8 = f_ptr[8].at(lat_code);
//         D_Phy_DDF f9 = f_ptr[9].at(lat_code);
//         D_Phy_DDF f10= f_ptr[10].at(lat_code);
//         D_Phy_DDF f11= f_ptr[11].at(lat_code);
//         D_Phy_DDF f14= f_ptr[14].at(lat_code);
//         D_Phy_DDF f15= f_ptr[15].at(lat_code);
//         D_Phy_DDF f18= f_ptr[18].at(lat_code);
//         D_Phy_DDF nx_z = 1./2. * (f1+f7+f9-f2-f10-f8);
//         D_Phy_DDF ny_z = 1./2. * (f3+f7+f10-f4-f9-f8);
//         density_ptr->at(lat_code) = f0+f1+f2+f3+f4+f7+f9+f10+f8 + 2.*(f5+f11+f14+f15+f18);
//         f_ptr[6].at(lat_code) = f5;
//         f_ptr[12].at(lat_code) = f11 + nx_z;
//         f_ptr[13].at(lat_code) = f14 - nx_z;
//         f_ptr[16].at(lat_code) = f15 + ny_z;
//         f_ptr[17].at(lat_code) = f18 - ny_z;
//     }
// }

// void nonEquilibrumBounceBack_NoSlip_Bot::initialBCStrategy()
// {
//     D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
//     initial_Feq(botBc_ptr);
// }

// void nonEquilibrumBounceBack_NoSlip_Bot::applyBCStrategy()
// {
//     D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
//     for (auto i_bc_lat : *botBc_ptr)
//     {
//         D_morton lat_code = i_bc_lat.first;
//         D_Phy_DDF f0 = f_ptr[0].at(lat_code);
//         D_Phy_DDF f1 = f_ptr[1].at(lat_code);
//         D_Phy_DDF f2 = f_ptr[2].at(lat_code);
//         D_Phy_DDF f3 = f_ptr[3].at(lat_code);
//         D_Phy_DDF f4 = f_ptr[4].at(lat_code);
//         D_Phy_DDF f6 = f_ptr[6].at(lat_code);
//         D_Phy_DDF f7 = f_ptr[7].at(lat_code);
//         D_Phy_DDF f8 = f_ptr[8].at(lat_code);
//         D_Phy_DDF f9 = f_ptr[9].at(lat_code);
//         D_Phy_DDF f10= f_ptr[10].at(lat_code);
//         D_Phy_DDF f12= f_ptr[12].at(lat_code);
//         D_Phy_DDF f13= f_ptr[13].at(lat_code);
//         D_Phy_DDF f16= f_ptr[16].at(lat_code);
//         D_Phy_DDF f17= f_ptr[17].at(lat_code);
//         D_Phy_DDF nx_z = 1./2. * (f1+f7+f9-f2-f10-f8);
//         D_Phy_DDF ny_z = 1./2. * (f3+f7+f10-f4-f9-f8);
//         density_ptr->at(lat_code) = f0+f1+f2+f3+f4+f7+f9+f10+f8 + 2.*(f6+f13+f12+f17+f16);
//         f_ptr[5].at(lat_code) = f6;
//         f_ptr[11].at(lat_code) = f12 - nx_z;
//         f_ptr[14].at(lat_code) = f13 + nx_z;
//         f_ptr[15].at(lat_code) = f16 - ny_z;
//         f_ptr[18].at(lat_code) = f17 + ny_z;
//     }
// }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyBCStrategy()
// // {
// // if (all_edge_same) {
// //     applyOnEdgeBC_NE();
// //     applyOnEdgeBC_NW();
// //     applyOnEdgeBC_SE();
// //     applyOnEdgeBC_SW();
// //     applyOnEdgeBC_TE();
// //     applyOnEdgeBC_TW();
// //     applyOnEdgeBC_BE();
// //     applyOnEdgeBC_BW();
// //     applyOnEdgeBC_TN();
// //     applyOnEdgeBC_TS();
// //     applyOnEdgeBC_BN();
// //     applyOnEdgeBC_BS();
// // }
// // else {
    
// // }

// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_NE()
// // {
// //     D_mapLat* NEedge_ptr = &(Lat_Manager::pointer_me->lat_bc_NEedge);
// //     for (auto i_bc_lat : *NEedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxy = 1./4. * (df_ptr->f5.at(lat_code) - df_ptr->f6.at(lat_code));
// //         df_ptr->f2.at(lat_code) = df_ptr->f1.at(lat_code);
// //         df_ptr->f4.at(lat_code) = df_ptr->f3.at(lat_code);
// //         df_ptr->f12.at(lat_code) = df_ptr->f7.at(lat_code);
// //         df_ptr->f13.at(lat_code) = df_ptr->f10.at(lat_code) + nxy;
// //         df_ptr->f14.at(lat_code) = df_ptr->f9.at(lat_code) - nxy;
// //         df_ptr->f17.at(lat_code) = df_ptr->f16.at(lat_code) + nxy;
// //         df_ptr->f18.at(lat_code) = df_ptr->f15.at(lat_code) - nxy;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f9.at(lat_code) +
// //                                                      df_ptr->f10.at(lat_code) +
// //                                                       df_ptr->f12.at(lat_code) +
// //                                                        df_ptr->f13.at(lat_code) +
// //                                                         df_ptr->f14.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f8.at(lat_code) = buried_link_df;
// //         df_ptr->f11.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_NW()
// // {
// //     D_mapLat* NWedge_ptr = &(Lat_Manager::pointer_me->lat_bc_NWedge);
// //     for (auto i_bc_lat : *NWedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxy = 1./4. * (df_ptr->f5.at(lat_code) - df_ptr->f6.at(lat_code));
// //         df_ptr->f1.at(lat_code) = df_ptr->f2.at(lat_code);
// //         df_ptr->f4.at(lat_code) = df_ptr->f3.at(lat_code);
// //         df_ptr->f8.at(lat_code) = df_ptr->f11.at(lat_code);
// //         df_ptr->f9.at(lat_code) = df_ptr->f14.at(lat_code) + nxy;
// //         df_ptr->f10.at(lat_code) = df_ptr->f13.at(lat_code) - nxy;
// //         df_ptr->f17.at(lat_code) = df_ptr->f16.at(lat_code) + nxy;
// //         df_ptr->f18.at(lat_code) = df_ptr->f15.at(lat_code) - nxy;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f8.at(lat_code) +
// //                                                     df_ptr->f9.at(lat_code) +
// //                                                      df_ptr->f10.at(lat_code) +
// //                                                       df_ptr->f11.at(lat_code) +
// //                                                        df_ptr->f13.at(lat_code) +
// //                                                         df_ptr->f14.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f7.at(lat_code) = buried_link_df;
// //         df_ptr->f12.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_SE()
// // {
// //     D_mapLat* SEedge_ptr = &(Lat_Manager::pointer_me->lat_bc_SEedge);
// //     for (auto i_bc_lat : *SEedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxy = 1./4. * (df_ptr->f5.at(lat_code) - df_ptr->f6.at(lat_code));
// //         df_ptr->f2.at(lat_code) = df_ptr->f1.at(lat_code);
// //         df_ptr->f3.at(lat_code) = df_ptr->f4.at(lat_code);
// //         df_ptr->f11.at(lat_code) = df_ptr->f8.at(lat_code);
// //         df_ptr->f13.at(lat_code) = df_ptr->f10.at(lat_code) + nxy;
// //         df_ptr->f14.at(lat_code) = df_ptr->f9.at(lat_code) - nxy;
// //         df_ptr->f15.at(lat_code) = df_ptr->f18.at(lat_code) + nxy;
// //         df_ptr->f16.at(lat_code) = df_ptr->f17.at(lat_code) - nxy;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f8.at(lat_code) +
// //                                                     df_ptr->f9.at(lat_code) +
// //                                                      df_ptr->f10.at(lat_code) +
// //                                                       df_ptr->f11.at(lat_code) +
// //                                                        df_ptr->f13.at(lat_code) +
// //                                                         df_ptr->f14.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f7.at(lat_code) = buried_link_df;
// //         df_ptr->f12.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_SW()
// // {
// //     D_mapLat* SWedge_ptr = &(Lat_Manager::pointer_me->lat_bc_SWedge);
// //     for (auto i_bc_lat : *SWedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxy = 1./4. * (df_ptr->f5.at(lat_code) - df_ptr->f6.at(lat_code));
// //         df_ptr->f1.at(lat_code) = df_ptr->f2.at(lat_code);
// //         df_ptr->f3.at(lat_code) = df_ptr->f4.at(lat_code);
// //         df_ptr->f7.at(lat_code) = df_ptr->f12.at(lat_code);
// //         df_ptr->f9.at(lat_code) = df_ptr->f14.at(lat_code) + nxy;
// //         df_ptr->f10.at(lat_code) = df_ptr->f13.at(lat_code) - nxy;
// //         df_ptr->f15.at(lat_code) = df_ptr->f18.at(lat_code) + nxy;
// //         df_ptr->f16.at(lat_code) = df_ptr->f17.at(lat_code) - nxy;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f9.at(lat_code) +
// //                                                      df_ptr->f10.at(lat_code) +
// //                                                       df_ptr->f12.at(lat_code) +
// //                                                        df_ptr->f13.at(lat_code) +
// //                                                         df_ptr->f14.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f8.at(lat_code) = buried_link_df;
// //         df_ptr->f11.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_TE()
// // {
// //     D_mapLat* TEedge_ptr = &(Lat_Manager::pointer_me->lat_bc_TEedge);
// //     for (auto i_bc_lat : *TEedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxz = 1./4. * (df_ptr->f3.at(lat_code) - df_ptr->f4.at(lat_code));
// //         df_ptr->f2.at(lat_code) = df_ptr->f1.at(lat_code);
// //         df_ptr->f6.at(lat_code) = df_ptr->f5.at(lat_code);
// //         df_ptr->f14.at(lat_code) = df_ptr->f9.at(lat_code);
// //         df_ptr->f11.at(lat_code) = df_ptr->f8.at(lat_code) + nxz;
// //         df_ptr->f12.at(lat_code) = df_ptr->f7.at(lat_code) - nxz;
// //         df_ptr->f16.at(lat_code) = df_ptr->f17.at(lat_code) + nxz;
// //         df_ptr->f18.at(lat_code) = df_ptr->f15.at(lat_code) - nxz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f9.at(lat_code) +
// //                                                       df_ptr->f11.at(lat_code) +
// //                                                        df_ptr->f12.at(lat_code) +
// //                                                         df_ptr->f14.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f10.at(lat_code) = buried_link_df;
// //         df_ptr->f13.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_TW()
// // {
// //     D_mapLat* TWedge_ptr = &(Lat_Manager::pointer_me->lat_bc_TWedge);
// //     for (auto i_bc_lat : *TWedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxz = 1./4. * (df_ptr->f3.at(lat_code) - df_ptr->f4.at(lat_code));
// //         df_ptr->f1.at(lat_code) = df_ptr->f2.at(lat_code);
// //         df_ptr->f6.at(lat_code) = df_ptr->f5.at(lat_code);
// //         df_ptr->f10.at(lat_code) = df_ptr->f13.at(lat_code);
// //         df_ptr->f7.at(lat_code) = df_ptr->f12.at(lat_code) + nxz;
// //         df_ptr->f8.at(lat_code) = df_ptr->f11.at(lat_code) - nxz;
// //         df_ptr->f16.at(lat_code) = df_ptr->f17.at(lat_code) + nxz;
// //         df_ptr->f18.at(lat_code) = df_ptr->f15.at(lat_code) - nxz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f10.at(lat_code) +
// //                                                       df_ptr->f11.at(lat_code) +
// //                                                        df_ptr->f12.at(lat_code) +
// //                                                         df_ptr->f13.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f9.at(lat_code) = buried_link_df;
// //         df_ptr->f14.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_BE()
// // {
// //     D_mapLat* BEedge_ptr = &(Lat_Manager::pointer_me->lat_bc_BEedge);
// //     for (auto i_bc_lat : *BEedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxz = 1./4. * (df_ptr->f3.at(lat_code) - df_ptr->f4.at(lat_code));
// //         df_ptr->f2.at(lat_code) = df_ptr->f1.at(lat_code);
// //         df_ptr->f5.at(lat_code) = df_ptr->f6.at(lat_code);
// //         df_ptr->f13.at(lat_code) = df_ptr->f10.at(lat_code);
// //         df_ptr->f11.at(lat_code) = df_ptr->f8.at(lat_code) + nxz;
// //         df_ptr->f12.at(lat_code) = df_ptr->f7.at(lat_code) - nxz;
// //         df_ptr->f15.at(lat_code) = df_ptr->f18.at(lat_code) + nxz;
// //         df_ptr->f17.at(lat_code) = df_ptr->f16.at(lat_code) - nxz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f10.at(lat_code) +
// //                                                       df_ptr->f11.at(lat_code) +
// //                                                        df_ptr->f12.at(lat_code) +
// //                                                         df_ptr->f13.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f9.at(lat_code) = buried_link_df;
// //         df_ptr->f14.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_BW()
// // {
// //     D_mapLat* BWedge_ptr = &(Lat_Manager::pointer_me->lat_bc_BWedge);
// //     for (auto i_bc_lat : *BWedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nxz = 1./4. * (df_ptr->f3.at(lat_code) - df_ptr->f4.at(lat_code));
// //         df_ptr->f2.at(lat_code) = df_ptr->f1.at(lat_code);
// //         df_ptr->f5.at(lat_code) = df_ptr->f6.at(lat_code);
// //         df_ptr->f9.at(lat_code) = df_ptr->f14.at(lat_code);
// //         df_ptr->f7.at(lat_code) = df_ptr->f12.at(lat_code) + nxz;
// //         df_ptr->f8.at(lat_code) = df_ptr->f11.at(lat_code) - nxz;
// //         df_ptr->f15.at(lat_code) = df_ptr->f18.at(lat_code) + nxz;
// //         df_ptr->f17.at(lat_code) = df_ptr->f16.at(lat_code) - nxz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f9.at(lat_code) +
// //                                                       df_ptr->f11.at(lat_code) +
// //                                                        df_ptr->f12.at(lat_code) +
// //                                                         df_ptr->f14.at(lat_code) +
// //                                                          df_ptr->f15.at(lat_code) +
// //                                                           df_ptr->f16.at(lat_code) +
// //                                                            df_ptr->f17.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f10.at(lat_code) = buried_link_df;
// //         df_ptr->f13.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_TN()
// // {
// //     D_mapLat* TNedge_ptr = &(Lat_Manager::pointer_me->lat_bc_TNedge);
// //     for (auto i_bc_lat : *TNedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nyz = 1./4. * (df_ptr->f1.at(lat_code) - df_ptr->f2.at(lat_code));
// //         df_ptr->f4.at(lat_code) = df_ptr->f3.at(lat_code);
// //         df_ptr->f6.at(lat_code) = df_ptr->f5.at(lat_code);
// //         df_ptr->f18.at(lat_code) = df_ptr->f15.at(lat_code);
// //         df_ptr->f8.at(lat_code) = df_ptr->f11.at(lat_code) + nyz;
// //         df_ptr->f12.at(lat_code) = df_ptr->f7.at(lat_code) - nyz;
// //         df_ptr->f10.at(lat_code) = df_ptr->f13.at(lat_code) + nyz;
// //         df_ptr->f14.at(lat_code) = df_ptr->f9.at(lat_code) - nyz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f9.at(lat_code) +
// //                                                       df_ptr->f10.at(lat_code) +
// //                                                        df_ptr->f11.at(lat_code) +
// //                                                         df_ptr->f12.at(lat_code) +
// //                                                          df_ptr->f13.at(lat_code) +
// //                                                           df_ptr->f14.at(lat_code) +
// //                                                            df_ptr->f15.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f16.at(lat_code) = buried_link_df;
// //         df_ptr->f17.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_TS()
// // {
// //     D_mapLat* TSedge_ptr = &(Lat_Manager::pointer_me->lat_bc_TSedge);
// //     for (auto i_bc_lat : *TSedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nyz = 1./4. * (df_ptr->f1.at(lat_code) - df_ptr->f2.at(lat_code));
// //         df_ptr->f3.at(lat_code) = df_ptr->f4.at(lat_code);
// //         df_ptr->f6.at(lat_code) = df_ptr->f5.at(lat_code);
// //         df_ptr->f16.at(lat_code) = df_ptr->f17.at(lat_code);
// //         df_ptr->f7.at(lat_code) = df_ptr->f12.at(lat_code) + nyz;
// //         df_ptr->f11.at(lat_code) = df_ptr->f8.at(lat_code) - nyz;
// //         df_ptr->f10.at(lat_code) = df_ptr->f13.at(lat_code) + nyz;
// //         df_ptr->f14.at(lat_code) = df_ptr->f9.at(lat_code) - nyz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f9.at(lat_code) +
// //                                                       df_ptr->f10.at(lat_code) +
// //                                                        df_ptr->f11.at(lat_code) +
// //                                                         df_ptr->f12.at(lat_code) +
// //                                                          df_ptr->f13.at(lat_code) +
// //                                                           df_ptr->f14.at(lat_code) +
// //                                                            df_ptr->f16.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f15.at(lat_code) = buried_link_df;
// //         df_ptr->f18.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_BN()
// // {
// //     D_mapLat* BNedge_ptr = &(Lat_Manager::pointer_me->lat_bc_BNedge);
// //     for (auto i_bc_lat : *BNedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nyz = 1./4. * (df_ptr->f1.at(lat_code) - df_ptr->f2.at(lat_code));
// //         df_ptr->f4.at(lat_code) = df_ptr->f3.at(lat_code);
// //         df_ptr->f5.at(lat_code) = df_ptr->f6.at(lat_code);
// //         df_ptr->f17.at(lat_code) = df_ptr->f16.at(lat_code);
// //         df_ptr->f8.at(lat_code) = df_ptr->f11.at(lat_code) + nyz;
// //         df_ptr->f12.at(lat_code) = df_ptr->f7.at(lat_code) - nyz;
// //         df_ptr->f9.at(lat_code) = df_ptr->f14.at(lat_code) + nyz;
// //         df_ptr->f13.at(lat_code) = df_ptr->f10.at(lat_code) - nyz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f9.at(lat_code) +
// //                                                       df_ptr->f10.at(lat_code) +
// //                                                        df_ptr->f11.at(lat_code) +
// //                                                         df_ptr->f12.at(lat_code) +
// //                                                          df_ptr->f13.at(lat_code) +
// //                                                           df_ptr->f14.at(lat_code) +
// //                                                            df_ptr->f16.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f15.at(lat_code) = buried_link_df;
// //         df_ptr->f18.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Edge::applyOnEdgeBC_BS()
// // {
// //     D_mapLat* BSedge_ptr = &(Lat_Manager::pointer_me->lat_bc_BSedge);
// //     for (auto i_bc_lat : *BSedge_ptr)
// //     {
// //         D_morton lat_code = i_bc_lat.first;
// //         D_Phy_DDF nyz = 1./4. * (df_ptr->f1.at(lat_code) - df_ptr->f2.at(lat_code));
// //         df_ptr->f3.at(lat_code) = df_ptr->f4.at(lat_code);
// //         df_ptr->f5.at(lat_code) = df_ptr->f6.at(lat_code);
// //         df_ptr->f15.at(lat_code) = df_ptr->f18.at(lat_code);
// //         df_ptr->f7.at(lat_code) = df_ptr->f12.at(lat_code) + nyz;
// //         df_ptr->f11.at(lat_code) = df_ptr->f8.at(lat_code) - nyz;
// //         df_ptr->f9.at(lat_code) = df_ptr->f14.at(lat_code) + nyz;
// //         df_ptr->f13.at(lat_code) = df_ptr->f10.at(lat_code) - nyz;
// //         // buried-link population
// //         D_Phy_DDF buried_link_df = 1./22. * (df_ptr->f1.at(lat_code) + 
// //                                               df_ptr->f2.at(lat_code) + 
// //                                                df_ptr->f3.at(lat_code) + 
// //                                                 df_ptr->f4.at(lat_code) + 
// //                                                  df_ptr->f5.at(lat_code) + 
// //                                                   df_ptr->f6.at(lat_code) +
// //                                                    df_ptr->f7.at(lat_code) +
// //                                                     df_ptr->f8.at(lat_code) +
// //                                                      df_ptr->f9.at(lat_code) +
// //                                                       df_ptr->f10.at(lat_code) +
// //                                                        df_ptr->f11.at(lat_code) +
// //                                                         df_ptr->f12.at(lat_code) +
// //                                                          df_ptr->f13.at(lat_code) +
// //                                                           df_ptr->f14.at(lat_code) +
// //                                                            df_ptr->f15.at(lat_code) +
// //                                                             df_ptr->f18.at(lat_code));
// //         df_ptr->f16.at(lat_code) = buried_link_df;
// //         df_ptr->f17.at(lat_code) = buried_link_df;
// //         df_ptr->f0.at(lat_code) = 12. * buried_link_df;
// //     }
// // }

// // void nonEquilibrumBounceBack_NoSlip_Corner::applyBCStrategy()
// // {
// //
// // }

