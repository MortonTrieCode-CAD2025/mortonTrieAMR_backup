#include "Lat_Manager.hpp"
#include "Solid_Manager.h"
#include "Morton_Assist.h"

void Lat_Manager::voxelize()
{
    update_boundary_lattices();

    update_outerLevel_lattices();

    update_innerLevel_lattices();

    update_overlap_lattices();

    // Print the statistic of Lattice Num
    std::cout << "----------- Lattice Number Count -----------" << std::endl;
    for (D_int i_level = 0; i_level <= C_max_level; ++i_level) {
        std::cout << "  Level " << i_level << " " << lat_f.at(i_level).size() << " overlap_C2F " << lat_overlap_C2F.at(i_level).size() << std::endl;
    }
    std::cout << "----------- -------------------- -----------" << std::endl;
}


void Lat_Manager::update_boundary_lattices()
{
    // Unified boundary handling for arbitrary solid positioning
    D_mapNodePtr* bdEast = &(Grid_Manager::pointer_me->bk_boundary_x.at(0).at(1));
    D_mapNodePtr* bdWest = &(Grid_Manager::pointer_me->bk_boundary_x.at(0).at(0));
    D_mapNodePtr* bdNorth = &(Grid_Manager::pointer_me->bk_boundary_y.at(0).at(1));
    D_mapNodePtr* bdSouth = &(Grid_Manager::pointer_me->bk_boundary_y.at(0).at(0));
    D_mapNodePtr* bdTop = &(Grid_Manager::pointer_me->bk_boundary_z.at(0).at(1));
    D_mapNodePtr* bdBot = &(Grid_Manager::pointer_me->bk_boundary_z.at(0).at(0));


#ifndef BC_NO_GHOST
    nNodeX = Grid_Manager::pointer_me->nx - 3;
    nNodeY = Grid_Manager::pointer_me->ny - 3;
    nNodeZ = Grid_Manager::pointer_me->nz - 3;
#else
    nNodeX = Grid_Manager::pointer_me->nx - 1;
    nNodeY = Grid_Manager::pointer_me->ny - 1;
    nNodeZ = Grid_Manager::pointer_me->nz - 1;
#endif


#ifndef BC_NO_GHOST
    for (auto bdry_node : *bdWest)
    {
        D_morton key_current = bdry_node.first;
        D_morton key_x1 = Morton_Assist::find_x1(key_current, 0);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, 0);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, 0);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, 0);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, 0);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, 0);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, 0);

        if (bdWest->find(key_yz1) == bdWest->end())
        // beyond the envelope
            continue;

        ghost_bc.insert(make_pair(key_current, Cell(Cell_Flag::BDRY_GHOST)));

        if (bdNorth->find(key_y1) != bdNorth->end() ||
             bdSouth->find(key_current) != bdSouth->end() ||
              bdTop->find(key_z1) != bdTop->end() ||
               bdBot->find(key_current) != bdBot->end())
            continue;
        else {
            lat_bc_x.at(0).insert(make_pair(key_x1, Cell(Cell_Flag::BDRY_FACE_WEST)));
            lat_f.at(0).insert(make_pair(key_x1, Cell(Cell_Flag::BDRY_FACE_WEST)));
        }
    }

    for (auto bdry_node : *bdEast) 
    {
        D_morton key_x1 = bdry_node.first;
        D_morton key_current = Morton_Assist::find_x0(key_x1, 0);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, 0);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, 0);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, 0);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, 0);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, 0);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, 0);

        // if (bdEast->find(key_xz1)==bdEast->end() || bdEast->find(key_xy1)==bdEast->end())
        if (bdEast->find(key_xyz1)==bdEast->end())
            // beyond the envelope
            continue;
        
        ghost_bc.insert(make_pair(key_current, Cell(Cell_Flag::BDRY_GHOST)));

        if (bdNorth->find(key_y1) != bdNorth->end() || 
             bdSouth->find(key_current) != bdSouth->end() ||
              bdTop->find(key_z1) != bdTop->end() ||
               bdBot->find(key_current) != bdBot->end())
            continue;
        else {
            lat_bc_x.at(1).insert(make_pair(Morton_Assist::find_x0(key_current, 0), Cell(Cell_Flag::BDRY_FACE_EAST)));
            lat_f.at(0).insert(make_pair(Morton_Assist::find_x0(key_current, 0), Cell(Cell_Flag::BDRY_FACE_EAST)));
        }
    }

    for (auto bdry_node : *bdSouth)
    {
        D_morton key_current = bdry_node.first;
        D_morton key_x1  = Morton_Assist::find_x1(key_current, 0);
        D_morton key_y1  = Morton_Assist::find_y1(key_current, 0);
        D_morton key_z1  = Morton_Assist::find_z1(key_current, 0);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, 0);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, 0);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, 0);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, 0);

        // if (bdSouth->find(key_xz1)==bdSouth->end() || bdWest->find(key_yz1)==bdWest->end() || bdEast->find(key_xyz1)==bdEast->end())
        if (bdSouth->find(key_xz1)==bdSouth->end())
            // beyond the envelope
            continue;

        ghost_bc.insert(make_pair(key_current, Cell(Cell_Flag::BDRY_GHOST)));

        if (bdEast->find(key_x1) != bdEast->end() || 
             bdWest->find(key_current) != bdWest->end() ||
              bdTop->find(key_z1) != bdTop->end() ||
               bdBot->find(key_current) != bdBot->end())
            continue;
        else {
            lat_bc_y.at(0).insert(make_pair(key_y1, Cell(Cell_Flag::BDRY_FACE_SOUTH)));
            lat_f.at(0).insert(make_pair(key_y1, Cell(Cell_Flag::BDRY_FACE_SOUTH)));
        }
    }

    for (auto bdry_node : *bdNorth)
    {
        D_morton key_y1 = bdry_node.first;
        D_morton key_current = Morton_Assist::find_y0(key_y1, 0);
        D_morton key_x1 = Morton_Assist::find_x1(key_current, 0);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, 0);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, 0);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, 0);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, 0);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, 0);

        // if (bdNorth->find(key_xyz1)==bdNorth->end() || bdWest->find(key_yz1)==bdWest->end() || bdEast->find(key_xyz1)==bdEast->end())
        if (bdNorth->find(key_xyz1)==bdNorth->end())
            // beyond the envelope
            continue;

        ghost_bc.insert(make_pair(key_current, Cell(Cell_Flag::BDRY_GHOST)));

        if (bdEast->find(key_x1) != bdEast->end() || 
             bdWest->find(key_current) != bdWest->end() ||
              bdTop->find(key_z1) != bdTop->end() ||
               bdBot->find(key_current) != bdBot->end())
            continue;
        else {
            lat_bc_y.at(1).insert(make_pair(Morton_Assist::find_y0(key_current, 0), Cell(Cell_Flag::BDRY_FACE_NORTH)));
            lat_f.at(0).insert(make_pair(Morton_Assist::find_y0(key_current, 0), Cell(Cell_Flag::BDRY_FACE_NORTH)));
        }
    }

    for (auto bdry_node : *bdBot)
    {
        D_morton key_current = bdry_node.first;
        D_morton key_x1 = Morton_Assist::find_x1(key_current, 0);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, 0);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, 0);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, 0);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, 0);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, 0);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, 0);

        if (bdBot->find(key_xy1)==bdBot->end())
            // beyond the envelope
            continue;

        // lat_bc_z.at(0).insert(make_pair(key_current, Cell(Cell_Flag::BDRY_FACE_BOT)));
        ghost_bc.insert(make_pair(key_current, Cell(Cell_Flag::BDRY_GHOST)));

        if (bdEast->find(key_x1) != bdEast->end() || 
             bdWest->find(key_current) != bdWest->end() ||
              bdNorth->find(key_y1) != bdNorth->end() ||
               bdSouth->find(key_current) != bdSouth->end())
            continue;
        else {
            lat_bc_z.at(0).insert(make_pair(key_z1, Cell(Cell_Flag::BDRY_FACE_BOT)));
            lat_f.at(0).insert(make_pair(key_z1, Cell(Cell_Flag::BDRY_FACE_BOT)));
        }
    }

    for (auto bdry_node : *bdTop)
    {
        D_morton key_z1 = bdry_node.first;
        D_morton key_current = Morton_Assist::find_z0(key_z1, 0);
        D_morton key_x1 = Morton_Assist::find_x1(key_current, 0);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, 0);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, 0);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, 0);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, 0);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, 0);

        if (bdTop->find(key_xyz1)==bdTop->end())
            // beyond the envelope
            continue;
        
        // lat_bc_z.at(1).insert(make_pair(key_current, Cell(Cell_Flag::BDRY_FACE_TOP)));
        ghost_bc.insert(make_pair(key_current, Cell(Cell_Flag::BDRY_GHOST)));

        if (bdEast->find(key_x1) != bdEast->end() || 
             bdWest->find(key_current) != bdWest->end() ||
              bdNorth->find(key_y1) != bdNorth->end() ||
               bdSouth->find(key_current) != bdSouth->end())
            continue;
        else {
            lat_bc_z.at(1).insert(make_pair(Morton_Assist::find_z0(key_current, 0), Cell(Cell_Flag::BDRY_FACE_TOP)));
            lat_f.at(0).insert(make_pair(Morton_Assist::find_z0(key_current, 0), Cell(Cell_Flag::BDRY_FACE_TOP)));
        }
    }

    std::cout << "----------- Geometry Boundary -----------\n";
    std::cout << " TOTAL BC Lats\t"   << lat_f.at(0).size()    << "\n";
    std::cout << "  - west size\t"  << lat_bc_x.at(0).size() << "\n";
    std::cout << "  - east size\t"  << lat_bc_x.at(1).size() << "\n";
    std::cout << "  - south size\t" << lat_bc_y.at(0).size() << "\n";
    std::cout << "  - north size\t" << lat_bc_y.at(1).size() << "\n";
    std::cout << "  - top size\t"   << lat_bc_z.at(0).size() << "\n";
    std::cout << "  - bot size\t"   << lat_bc_z.at(1).size() << "\n";
    std::cout << "----------- ----------------- -----------\n" << std::endl;

#else
    for (D_uint iz = 1; iz <= nNodeZ-2; ++iz) {
        for (D_uint iy = 1; iy <= nNodeY-2; ++iy) {
            lat_bc.at(Bdry_Type::BDRY_FACE_WEST).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,iy,iz)) << Morton_Assist::bit_otherlevel );
            lat_bc.at(Bdry_Type::BDRY_FACE_EAST).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,iy,iz)) << Morton_Assist::bit_otherlevel );
        }
    }
    for (D_uint iz = 1; iz <= nNodeZ-2; ++iz) {
        for (D_uint ix = 1; ix <= nNodeX-2; ++ix) {
            lat_bc.at(Bdry_Type::BDRY_FACE_SOUTH).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,0,iz)) << Morton_Assist::bit_otherlevel );
            lat_bc.at(Bdry_Type::BDRY_FACE_NORTH).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,nNodeY-1,iz)) << Morton_Assist::bit_otherlevel );
        }
    }
    for (D_uint iy = 1; iy <= nNodeY-2; ++iy) {
        for (D_uint ix = 1; ix <= nNodeX-2; ++ix) {
            lat_bc.at(Bdry_Type::BDRY_FACE_BOT).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,iy,0)) << Morton_Assist::bit_otherlevel );
            lat_bc.at(Bdry_Type::BDRY_FACE_TOP).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,iy,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
        }
    }

    for (D_uint iy = 1; iy <= nNodeY-2; ++iy) {
        lat_bc.at(Bdry_Type::BDRY_EDGE_WB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,iy,0)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_EB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,iy,0)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_WT).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,iy,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_ET).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,iy,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
    }
    for (D_uint ix = 1; ix <= nNodeX-2; ++ix) {
        lat_bc.at(Bdry_Type::BDRY_EDGE_SB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,0,0)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_NB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,nNodeY-1,0)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_ST).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,0,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_NT).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,nNodeY-1,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
    }
    for (D_uint iz = 1; iz <= nNodeZ-2; ++iz) {
        lat_bc.at(Bdry_Type::BDRY_EDGE_WS).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,iz)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_ES).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,0,iz)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_WN).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,nNodeY-1,iz)) << Morton_Assist::bit_otherlevel );
        lat_bc.at(Bdry_Type::BDRY_EDGE_EN).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,nNodeY-1,iz)) << Morton_Assist::bit_otherlevel );
    }

    lat_bc.at(Bdry_Type::BDRY_CORNER_WSB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,0)) << Morton_Assist::bit_otherlevel );
    lat_bc.at(Bdry_Type::BDRY_CORNER_WST).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
    lat_bc.at(Bdry_Type::BDRY_CORNER_WNB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,nNodeY-1,0)) << Morton_Assist::bit_otherlevel );
    lat_bc.at(Bdry_Type::BDRY_CORNER_WNT).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,nNodeY-1,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
    lat_bc.at(Bdry_Type::BDRY_CORNER_ESB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,0,0)) << Morton_Assist::bit_otherlevel );
    lat_bc.at(Bdry_Type::BDRY_CORNER_EST).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,0,nNodeZ-1)) << Morton_Assist::bit_otherlevel );
    lat_bc.at(Bdry_Type::BDRY_CORNER_ENB).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,nNodeY-1,0)) << Morton_Assist::bit_otherlevel );
    lat_bc.at(Bdry_Type::BDRY_CORNER_ENT).insert( static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(nNodeX-1,nNodeY-1,nNodeZ-1)) << Morton_Assist::bit_otherlevel );

    // std::cout << "----------- Geometry Boundary -----------\n";
    // // std::cout << " TOTAL BC Lats\t"   << lat_f.at(0).size()    << "\n";
    // std::cout << "  * west  size\t" << lat_bc.at(Bdry_Type::BDRY_FACE_WEST).size() << "\n";
    // std::cout << "  * east  size\t" << lat_bc.at(Bdry_Type::BDRY_FACE_EAST).size() << "\n";
    // std::cout << "  * south size\t" << lat_bc.at(Bdry_Type::BDRY_FACE_SOUTH).size() << "\n";
    // std::cout << "  * north size\t" << lat_bc.at(Bdry_Type::BDRY_FACE_NORTH).size() << "\n";
    // std::cout << "  * top   size\t" << lat_bc.at(Bdry_Type::BDRY_FACE_TOP).size() << "\n";
    // std::cout << "  * bot   size\t" << lat_bc.at(Bdry_Type::BDRY_FACE_BOT).size() << "\n\n";

    // std::cout << "  * bot-west   size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_BW).size() << "\n";
    // std::cout << "  * bot-east   size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_BE).size() << "\n";
    // std::cout << "  * bot-south  size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_BS).size() << "\n";
    // std::cout << "  * bot-north  size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_BN).size() << "\n";
    // std::cout << "  * top-west   size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_TW).size() << "\n";
    // std::cout << "  * top-east   size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_TE).size() << "\n";
    // std::cout << "  * top-south  size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_TS).size() << "\n";
    // std::cout << "  * top-north  size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_TN).size() << "\n";
    // std::cout << "  * west-south size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_SW).size() << "\n";
    // std::cout << "  * west-north size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_NW).size() << "\n";
    // std::cout << "  * east-south size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_SE).size() << "\n";
    // std::cout << "  * east-north size\t" << lat_bc.at(Bdry_Type::BDRY_EDGE_NE).size() << "\n\n";

    // std::cout << "  * south-west-bot size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_SWB).size() << "\n";
    // std::cout << "  * south-west-top size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_SWT).size() << "\n";
    // std::cout << "  * south-east-bot size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_SEB).size() << "\n";
    // std::cout << "  * south-east-top size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_SET).size() << "\n";
    // std::cout << "  * north-west-bot size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_NWB).size() << "\n";
    // std::cout << "  * north-west-top size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_NWT).size() << "\n";
    // std::cout << "  * north-east-bot size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_NEB).size() << "\n";
    // std::cout << "  * north-east-top size\t" << lat_bc.at(Bdry_Type::BDRY_CORNER_NET).size() << "\n";

    // std::cout << "----------- ----------------- -----------\n" << std::endl;

#endif
}


/**
 * @brief update all lattices but the inner level
 * 
 */
void Lat_Manager::update_outerLevel_lattices()
{
    D_setLat lat_bc_total;
    for (auto i_bc : lat_bc) {
        lat_bc_total.insert(i_bc.begin(),i_bc.end());
    }

    for (D_int i_level = 0; i_level < C_max_level; ++i_level)
    {
        D_mapLat* lat = &(lat_f.at(i_level));
        Grid_NIB* grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(i_level));
        for (auto iter : grid_ptr->grid)
        {
            D_morton key_current = iter.first;

#ifndef BC_NO_GHOST
            if (ghost_bc.find(key_current)!=ghost_bc.end() && i_level == 0)
                continue;
#else
            if (lat_bc_total.find(key_current)!=lat_bc_total.end() && i_level == 0)
                continue;
#endif

            D_morton key_x1 = Morton_Assist::find_x1(key_current, i_level);
            D_morton key_y1 = Morton_Assist::find_y1(key_current, i_level);
            D_morton key_z1 = Morton_Assist::find_z1(key_current, i_level);
            D_morton key_xy1 = Morton_Assist::find_y1(key_x1, i_level);
            D_morton key_xz1 = Morton_Assist::find_z1(key_x1, i_level);
            D_morton key_yz1 = Morton_Assist::find_z1(key_y1, i_level);
            D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, i_level);

            // insure all vertex at this level
            bool bool_xyz = (grid_ptr->grid.find(key_x1) != grid_ptr->grid.end()) && 
                             (grid_ptr->grid.find(key_y1) != grid_ptr->grid.end()) && 
                              (grid_ptr->grid.find(key_z1) != grid_ptr->grid.end()) && 
                               (grid_ptr->grid.find(key_xy1) != grid_ptr->grid.end()) && 
                                (grid_ptr->grid.find(key_xz1) != grid_ptr->grid.end()) && 
                                 (grid_ptr->grid.find(key_yz1) != grid_ptr->grid.end()) && 
                                  (grid_ptr->grid.find(key_xyz1) != grid_ptr->grid.end());

            if (!bool_xyz) 
                continue;
            
            int cell_flag = 999;

            if ((grid_ptr->grid.at(key_current).flag < 0 ||
                  grid_ptr->grid.at(key_x1).flag < 0 ||
                   grid_ptr->grid.at(key_y1).flag < 0 ||
                    grid_ptr->grid.at(key_z1).flag < 0 ||
                     grid_ptr->grid.at(key_xy1).flag < 0 ||
                      grid_ptr->grid.at(key_xz1).flag < 0 ||
                       grid_ptr->grid.at(key_yz1).flag < 0 ||
                        grid_ptr->grid.at(key_xyz1).flag < 0) && (i_level > 0)) {
                    lat_overlap_C2F.at(i_level).insert(make_pair(key_current, Cell(Cell_Flag::OVERLAP_C2F)));
            } else {
                lat->insert(make_pair(key_current, Cell(Cell_Flag::FLUID)));
            }
            
        } // for grid
    } // for level

}


void Lat_Manager::update_innerLevel_lattices()
{
    Grid_IB* solidgrid_ptr = &(Grid_Manager::pointer_me->gr_inner);

    //////////////////////////////////////////////////
    //  The step of constructing the inner lattice  //
    //////////////////////////////////////////////////
    //
    // lat_f[C_max_level] = FLUID + IB (need to execute the LBM simulation)
    // lat_sf = SURFACE + SOLID (no need to LBM)
    // 1. Classify the interface layer of overlapping 
    // 2. Convinced SOLID lattice  --+--> laf_sf
    // 3. Others classified into TBD --+--> lat_f[C_max_level]
    // 4. [func](update_cutting_lattices) Classify TBD of lat_f into IB or SURFACE
    //              IB --+--> lat_f       SURFACE --+--> lat_sf
    // 5. [func](propagate_fluid_lattices) Classify TBD of lat_f into FLUID or SOLID
    //            FLUID --+--> lat_f       SOLID --+--> lat_sf


    for (auto iter : solidgrid_ptr->grid)
    {
        dx.at(C_max_level) = solidgrid_ptr->get_dx();

        D_morton key_current = iter.first;
        D_morton key_x1 = Morton_Assist::find_x1(key_current, C_max_level);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, C_max_level);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, C_max_level);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, C_max_level);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, C_max_level);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, C_max_level);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, C_max_level);

        // insure all vertex at this level, not cross the inner level
        bool bool_xyz = (solidgrid_ptr->grid.find(key_x1) != solidgrid_ptr->grid.end()) && 
                         (solidgrid_ptr->grid.find(key_y1) != solidgrid_ptr->grid.end()) && 
                          (solidgrid_ptr->grid.find(key_z1) != solidgrid_ptr->grid.end()) && 
                           (solidgrid_ptr->grid.find(key_xy1) != solidgrid_ptr->grid.end()) && 
                            (solidgrid_ptr->grid.find(key_xz1) != solidgrid_ptr->grid.end()) && 
                             (solidgrid_ptr->grid.find(key_yz1) != solidgrid_ptr->grid.end()) && 
                              (solidgrid_ptr->grid.find(key_xyz1) != solidgrid_ptr->grid.end());

        if (!bool_xyz)
            continue;

        if (solidgrid_ptr->grid.at(key_current).flag < 0 &&
             solidgrid_ptr->grid.at(key_x1).flag < 0 &&
              solidgrid_ptr->grid.at(key_y1).flag < 0 &&
               solidgrid_ptr->grid.at(key_z1).flag < 0 &&
                solidgrid_ptr->grid.at(key_xy1).flag < 0 &&
                 solidgrid_ptr->grid.at(key_xz1).flag < 0 &&
                  solidgrid_ptr->grid.at(key_yz1).flag < 0 &&
                   solidgrid_ptr->grid.at(key_xyz1).flag < 0) {
            // the overlap lattices for interpolation between coarse and fine
            lat_overlap_C2F.at(C_max_level).insert(make_pair(key_current,Cell(Cell_Flag::OVERLAP_C2F)));
        }
        else if (solidgrid_ptr->grid.at(key_current).flag == Grid_Manager::pointer_me->flag_ghost ||
                  solidgrid_ptr->grid.at(key_x1).flag == Grid_Manager::pointer_me->flag_ghost ||
                   solidgrid_ptr->grid.at(key_y1).flag == Grid_Manager::pointer_me->flag_ghost ||
                    solidgrid_ptr->grid.at(key_z1).flag == Grid_Manager::pointer_me->flag_ghost ||
                     solidgrid_ptr->grid.at(key_xy1).flag == Grid_Manager::pointer_me->flag_ghost ||
                      solidgrid_ptr->grid.at(key_xz1).flag == Grid_Manager::pointer_me->flag_ghost ||
                       solidgrid_ptr->grid.at(key_yz1).flag == Grid_Manager::pointer_me->flag_ghost ||
                        solidgrid_ptr->grid.at(key_xyz1).flag == Grid_Manager::pointer_me->flag_ghost) {
            // the lattices inside the SOLID
            lat_sf.insert(make_pair(key_current,Cell(Cell_Flag::SOLID)));
        } else {
            lat_f.at(C_max_level).insert(make_pair(key_current,Cell(Cell_Flag::TBD)));
        }

    }

    update_cutting_lattices();

    // propagate_fluid_lattices();

    update_solid_lattices();

    // std::unordered_map<int, int> innerLat_flags;
    // for (auto lat_iter : lat_f.at(C_max_level)) {
    //     if (innerLat_flags.find(lat_iter.second.flag) == innerLat_flags.end()) {
    //         innerLat_flags.emplace(lat_iter.second.flag, 1);
    //     }
    //     else {
    //         innerLat_flags.at(lat_iter.second.flag)++;
    //     }
    // }
    // for (auto flag_iter : innerLat_flags) {
    //     std::cout << " flag" << flag_iter.first << " " << flag_iter.second;
    // }
    // std::cout << "\n----------- -------------------- -----------\n" << std::endl;

}


/**
 *    Real Lattice domain
 *       _________ _________ _________ ____ ____ ____ ____ __ __ __ __ 
 *      |         |         |         | F  | F  | F  | F  |__|__|__|__|
 *      |    C    |    C    |    C    |____|____|____|____|__|__|__|__|
 *      |         |         |         | F  | F  | F  | F  |__|__|__|__|
 *      |_________|_________|_________|____|____|____|____|__|__|__|__|
 *
 */
void Lat_Manager::update_overlap_lattices()
{
    auto& lat_f_alias = lat_f;
    auto& overlap_c2f_alias = lat_overlap_C2F;

    std::array<D_setLat, C_max_level+1> grid_coarse2fine;
    std::array<D_setLat, C_max_level+1> grid_fine2coarse;

    for (D_int coarse_level = 0; coarse_level < C_max_level; ++coarse_level) {
        for (auto i_overlap_layer : Grid_Manager::pointer_me->gr_NoIB.at(coarse_level).coarse2fine) {
            for (auto grid_overlapC2F : i_overlap_layer) {
                grid_coarse2fine.at(coarse_level).insert(grid_overlapC2F.first);
            }
        } 
        
        for (auto i_overlap_layer : Grid_Manager::pointer_me->gr_NoIB.at(coarse_level).fine2coarse) {
            for (auto grid_overlapF2C : i_overlap_layer) {
                grid_fine2coarse.at(coarse_level).insert(grid_overlapF2C.first);
            }
        }
    }
    for (auto i_overlap_layer : Grid_Manager::pointer_me->gr_inner.coarse2fine) {
        for (auto grid_overlapC2F : i_overlap_layer) {
            grid_coarse2fine.at(C_max_level).insert(grid_overlapC2F.first);
        }
    }
    for (auto i_overlap_layer : Grid_Manager::pointer_me->gr_inner.fine2coarse) {
        for (auto grid_overlapF2C : i_overlap_layer) {
            grid_fine2coarse.at(C_max_level).insert(grid_overlapF2C.first);
        }
    }

    // complete the Grid_Manager::pointer_me->gr_inner.fine2coarse[0]
    for (D_int coarse_level = 0; coarse_level < C_max_level; ++coarse_level) {
        D_int fine_level = coarse_level + 1;
        D_mapint *grid_f2c0;
        if (fine_level == C_max_level) {
            grid_f2c0 = &(Grid_Manager::pointer_me->gr_inner.fine2coarse.at(0));
        } else {
            grid_f2c0 = &(Grid_Manager::pointer_me->gr_NoIB.at(fine_level).fine2coarse.at(0));
        }
        for (auto i_grid_coarse : grid_coarse2fine.at(coarse_level)) {
            D_morton x1_code   = Morton_Assist::find_x1(i_grid_coarse, coarse_level); bool is_x1_inC2F   = grid_coarse2fine.at(coarse_level).find(x1_code)  != grid_coarse2fine.at(coarse_level).end();
            D_morton y1_code   = Morton_Assist::find_y1(i_grid_coarse, coarse_level); bool is_y1_inC2F   = grid_coarse2fine.at(coarse_level).find(y1_code)  != grid_coarse2fine.at(coarse_level).end();
            D_morton xy1_code  = Morton_Assist::find_y1(x1_code, coarse_level);       bool is_xy1_inC2F  = grid_coarse2fine.at(coarse_level).find(xy1_code) != grid_coarse2fine.at(coarse_level).end();
            D_morton z1_code   = Morton_Assist::find_z1(i_grid_coarse, coarse_level); bool is_z1_inC2F   = grid_coarse2fine.at(coarse_level).find(z1_code)  != grid_coarse2fine.at(coarse_level).end();
            D_morton xz1_code  = Morton_Assist::find_z1(x1_code, coarse_level);       bool is_xz1_inC2F  = grid_coarse2fine.at(coarse_level).find(xz1_code) != grid_coarse2fine.at(coarse_level).end();
            D_morton yz1_code  = Morton_Assist::find_z1(y1_code, coarse_level);       bool is_yz1_inC2F  = grid_coarse2fine.at(coarse_level).find(yz1_code) != grid_coarse2fine.at(coarse_level).end();
            D_morton xyz1_code = Morton_Assist::find_z1(xy1_code, coarse_level);      bool is_xyz1_inC2F = grid_coarse2fine.at(coarse_level).find(xyz1_code)!= grid_coarse2fine.at(coarse_level).end();

            if (is_x1_inC2F && is_y1_inC2F && is_xy1_inC2F && is_z1_inC2F && is_xz1_inC2F && is_yz1_inC2F && is_xyz1_inC2F) {
                // overlap_F2C.at(coarse_level).insert(i_grid_coarse);
                lat_overlap_F2C.at(coarse_level).insert(make_pair(i_grid_coarse, Cell(Cell_Flag::OVERLAP_F2C)));
                
                std::array<D_morton, 27> overlap_shell_code;
                overlap_shell_code[0] = i_grid_coarse;        //x0y0z0
                 overlap_shell_code[1] = Morton_Assist::find_x1(overlap_shell_code[0], fine_level); //x1y0z0
                  overlap_shell_code[2] = Morton_Assist::find_x1(overlap_shell_code[1], fine_level); //x2y0z0

                overlap_shell_code[3] = Morton_Assist::find_y1(overlap_shell_code[0], fine_level); //x0y1z0
                 overlap_shell_code[4] = Morton_Assist::find_x1(overlap_shell_code[3], fine_level); //x1y1z0
                  overlap_shell_code[5] = Morton_Assist::find_x1(overlap_shell_code[4], fine_level); //x2y1z0

                overlap_shell_code[6] = Morton_Assist::find_y1(overlap_shell_code[3], fine_level); //x0y2z0
                 overlap_shell_code[7] = Morton_Assist::find_x1(overlap_shell_code[6], fine_level); //x1y2z0
                  overlap_shell_code[8] = Morton_Assist::find_x1(overlap_shell_code[7], fine_level); //x2y2z0

                overlap_shell_code[9] = Morton_Assist::find_z1(overlap_shell_code[0], fine_level); //x0y0z1
                 overlap_shell_code[10] = Morton_Assist::find_x1(overlap_shell_code[9], fine_level); //x1y0z1
                  overlap_shell_code[11] = Morton_Assist::find_x1(overlap_shell_code[10], fine_level); //x2y0z1

                overlap_shell_code[12] = Morton_Assist::find_y1(overlap_shell_code[9], fine_level); //x0y1z1
                 overlap_shell_code[13] = Morton_Assist::find_x1(overlap_shell_code[12], fine_level); //x1y1z1
                  overlap_shell_code[14] = Morton_Assist::find_x1(overlap_shell_code[13], fine_level); //x2y1z1

                overlap_shell_code[15] = Morton_Assist::find_y1(overlap_shell_code[12], fine_level); //x0y2z1
                 overlap_shell_code[16] = Morton_Assist::find_x1(overlap_shell_code[15], fine_level); //x1y2z1
                  overlap_shell_code[17] = Morton_Assist::find_x1(overlap_shell_code[16], fine_level); //x2y2z1
                
                overlap_shell_code[18] = Morton_Assist::find_z1(overlap_shell_code[9], fine_level); //x0y0z2
                 overlap_shell_code[19] = Morton_Assist::find_x1(overlap_shell_code[18], fine_level); //x1y0z2
                  overlap_shell_code[20] = Morton_Assist::find_x1(overlap_shell_code[19], fine_level); //x2y0z2

                overlap_shell_code[21] = Morton_Assist::find_y1(overlap_shell_code[18], fine_level); //x0y1z2
                 overlap_shell_code[22] = Morton_Assist::find_x1(overlap_shell_code[21], fine_level); //x1y1z2
                  overlap_shell_code[23] = Morton_Assist::find_x1(overlap_shell_code[22], fine_level); //x2y1z2

                overlap_shell_code[24] = Morton_Assist::find_y1(overlap_shell_code[21], fine_level); //x0y2z2
                 overlap_shell_code[25] = Morton_Assist::find_x1(overlap_shell_code[24], fine_level); //x1y2z2
                  overlap_shell_code[26] = Morton_Assist::find_x1(overlap_shell_code[25], fine_level); //x2y2z2

                for (auto iter_overlap_shell : overlap_shell_code) {
                    if (grid_fine2coarse.at(fine_level).find(iter_overlap_shell) == grid_fine2coarse.at(fine_level).end()) {
                        grid_f2c0->insert(make_pair(iter_overlap_shell, Grid_Manager::pointer_me->flag_interface[0]));
                        grid_fine2coarse.at(fine_level).insert(iter_overlap_shell);
                    }
                } // for overlap_shell_code
            }
        } // for C2F
    } // for level

    // Remove overlap_f2c from lat_f[C_max_level] 
    //  (overlap_f2c collision independently, and stream implicitly by coalescence)
    for (D_int coarse_level = 0; coarse_level < C_max_level; ++coarse_level) {
        for (auto overlap_f2c : lat_overlap_F2C.at(coarse_level)) {
            lat_f.at(coarse_level).erase(overlap_f2c.first);
        }
    }

    for (D_int fine_level = 1; fine_level <= C_max_level; ++fine_level) {
        D_mapLat temp;
        temp.swap(lat_overlap_C2F.at(fine_level));
        for (auto current_code : grid_fine2coarse.at(fine_level)) {
            // D_morton current_code = grid_overlapC2F;
            D_morton x1_code   = Morton_Assist::find_x1(current_code, fine_level); bool is_x1_inC2F   = grid_fine2coarse.at(fine_level).find(x1_code)  != grid_fine2coarse.at(fine_level).end();
            D_morton y1_code   = Morton_Assist::find_y1(current_code, fine_level); bool is_y1_inC2F   = grid_fine2coarse.at(fine_level).find(y1_code)  != grid_fine2coarse.at(fine_level).end();
            D_morton xy1_code  = Morton_Assist::find_y1(x1_code, fine_level);      bool is_xy1_inC2F  = grid_fine2coarse.at(fine_level).find(xy1_code) != grid_fine2coarse.at(fine_level).end();
            D_morton z1_code   = Morton_Assist::find_z1(current_code, fine_level); bool is_z1_inC2F   = grid_fine2coarse.at(fine_level).find(z1_code)  != grid_fine2coarse.at(fine_level).end();
            D_morton xz1_code  = Morton_Assist::find_z1(x1_code, fine_level);      bool is_xz1_inC2F  = grid_fine2coarse.at(fine_level).find(xz1_code) != grid_fine2coarse.at(fine_level).end();
            D_morton yz1_code  = Morton_Assist::find_z1(y1_code, fine_level);      bool is_yz1_inC2F  = grid_fine2coarse.at(fine_level).find(yz1_code) != grid_fine2coarse.at(fine_level).end();
            D_morton xyz1_code = Morton_Assist::find_z1(xy1_code, fine_level);     bool is_xyz1_inC2F = grid_fine2coarse.at(fine_level).find(xyz1_code)!= grid_fine2coarse.at(fine_level).end();

            if (is_x1_inC2F && is_y1_inC2F && is_xy1_inC2F && is_z1_inC2F && is_xz1_inC2F && is_yz1_inC2F && is_xyz1_inC2F) {
                lat_overlap_C2F.at(fine_level).insert(make_pair(current_code, Cell(Cell_Flag::OVERLAP_C2F)));
            }
        } 
    }

    for (D_int fine_level = 1; fine_level <= C_max_level; ++fine_level) {
        for (auto coarse_overlap_lat : lat_overlap_C2F.at(fine_level)) {
            D_morton c2f_code = coarse_overlap_lat.first;

            /**
             * @todo redundancy search for inner layer of overlap
             */
            auto add_ghost_c2f = [&lat_f_alias, &overlap_c2f_alias, fine_level] (D_morton ngbr_code, std::array<D_mapLat, C_max_level+1>& ghost_overlap_C2F)
            {
                if (lat_f_alias.at(fine_level).find(ngbr_code) == lat_f_alias.at(fine_level).end() && overlap_c2f_alias.at(fine_level).find(ngbr_code) == overlap_c2f_alias.at(fine_level).end())
                    ghost_overlap_C2F.at(fine_level).insert(make_pair(ngbr_code, Cell(Cell_Flag::OVERLAP_GHOST)));
            };

            // add_ghost_c2f(Morton_Assist::find_x0(c2f_code, fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_x1(c2f_code, fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_y0(c2f_code, fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_y1(c2f_code, fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z0(c2f_code, fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z1(c2f_code, fine_level), ghost_overlap_C2F);

            // add_ghost_c2f(Morton_Assist::find_y0(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_y1(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_y0(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_y1(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);

            // add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);

            // add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_y0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_y0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_y1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            // add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_y1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);

            for (D_int i_q = 1; i_q < C_Q; ++i_q) {
                add_ghost_c2f(Morton_Assist::find_neighbor(c2f_code, fine_level, ex[i_q], ey[i_q], ez[i_q]), ghost_overlap_C2F);
            }
        }
    }


    // Generate ghost_overlap_F2C
    for (D_int coarse_level = 0; coarse_level < C_max_level; ++coarse_level) {
        for (auto coarse_overlap_lat : lat_overlap_F2C.at(coarse_level)) {
            D_morton f2c_code = coarse_overlap_lat.first;

            auto add_ghost_f2c = [&lat_f_alias, coarse_level](D_morton ngbr_code, std::array<D_mapLat, C_max_level+1>& ghost_overlap_F2C)
            {
                if (lat_f_alias.at(coarse_level).find(ngbr_code) != lat_f_alias.at(coarse_level).end())
                    ghost_overlap_F2C.at(coarse_level).insert(make_pair(ngbr_code, Cell(Cell_Flag::OVERLAP_GHOST)));
            };

            // add_ghost_f2c(Morton_Assist::find_x0(f2c_code, coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_x1(f2c_code, coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_y0(f2c_code, coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_y1(f2c_code, coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z0(f2c_code, coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z1(f2c_code, coarse_level), ghost_overlap_F2C);

            // add_ghost_f2c(Morton_Assist::find_y0(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_y1(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_y0(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_y1(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);

            // add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);

            // add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_y0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_y0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_y1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            // add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_y1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            for (D_int i_q = 1; i_q < C_Q; ++i_q) {
                add_ghost_f2c(Morton_Assist::find_neighbor(f2c_code, coarse_level, ex[i_q], ey[i_q], ez[i_q]), ghost_overlap_F2C);
            }
        }
    }

}


void Lat_Manager::update_cutting_lattices()
{
    Grid_IB* solidgrid_ptr = &(Grid_Manager::pointer_me->gr_inner);
    double dx_innerLat = dx.at(C_max_level);
    double box_halfextent[3] = {dx_innerLat*0.5f, dx_innerLat*0.5f, dx_innerLat*0.5f};

    D_setLat add_surf;

    for (D_int ishape = 0; ishape < Solid_Manager::pointer_me->numb_solids; ++ishape)
    {
        std::vector<Solid_Face>* triFace_ptr = &(Solid_Manager::pointer_me->shape_solids.at(ishape).triFace);
        D_int triFace_idx = -1;

        for (auto iTri : *triFace_ptr)
        {
            triFace_idx++;
            Solid_Face iTri1 = iTri.triFace_offset(D_vec{iTri.faceNorm.at(0),iTri.faceNorm.at(1),iTri.faceNorm.at(2)}, 1e-6);
            Solid_Face iTri2 = iTri.triFace_offset(D_vec{iTri.faceNorm.at(0),iTri.faceNorm.at(1),iTri.faceNorm.at(2)},-1e-6);

            Solid_Face triFaceSet[3] = {iTri, iTri1, iTri2};
            
            for (int count = 0; count < 3; ++count) {
                Solid_Face tri = triFaceSet[count];

                // get tri facet AABB and normal vector
                VertexBox aabb = tri.create_triFace_AABB();

                // loop the voxel on vBox, judge the flag
                for (D_int i = aabb.get_p0_x(); i < aabb.get_p1_x(); ++i) {
                    for (D_int j = aabb.get_p0_y(); j < aabb.get_p1_y(); ++j) { 
                        for (D_int k = aabb.get_p0_z(); k < aabb.get_p1_z(); ++k) {

                            D_morton cell_in_box_code = Morton_Assist::pointer_me->morton_encode(i, j, k);

                            if (solidgrid_ptr->grid.find(cell_in_box_code) != solidgrid_ptr->grid.end())
                            {
                                aabb.lat_in_aabb.insert(cell_in_box_code);
                                Cell lat_in_box = lat_f.at(C_max_level).at(cell_in_box_code);

                                // check the intersection of current cell and facet
                                D_vec lat_center;
                                Morton_Assist::pointer_me->compute_coordinate(cell_in_box_code, C_max_level, lat_center.x, lat_center.y, lat_center.z);
                                double boxcenter[3] = {lat_center.x, lat_center.y, lat_center.z};
                                double vertex_coord[3][3] = {{tri.vertex1.x, tri.vertex1.y, tri.vertex1.z},
                                                              {tri.vertex2.x, tri.vertex2.y, tri.vertex2.z},
                                                               {tri.vertex3.x, tri.vertex3.y, tri.vertex3.z}};
                                bool same_dir = dot_product((lat_center-D_vec(tri.vertex1.x,tri.vertex1.y,tri.vertex1.z)), D_vec(tri.faceNorm[0],tri.faceNorm[1],tri.faceNorm[2])) > 1e-7;
                                if (Shape::triBoxOverlap(boxcenter, box_halfextent, vertex_coord))
                                {
                                    if (same_dir) {
                                        // If the triangleFace immerses the center of the lattice 
                                        //  (i.e., the center of the lattice is in the normal direction of the triangleFace),
                                        //  the lattice maybe SURFACE or IB.
                                        // The SOLID tag will be assigned by the function "propagate_fluid_lattices".

                                        add_surf.insert(cell_in_box_code); 
                                        triFaceID_in_srf[cell_in_box_code].push_back(std::array<D_int,2>{ishape, triFace_idx});
                                    } else {
                                        lat_f.at(C_max_level).at(cell_in_box_code).flag = IB;
                                        compute_distance(tri, cell_in_box_code, lat_center);
                                    }
                                    
                                }
                                else {
                                    if (!same_dir) {
                                        lat_f.at(C_max_level).at(cell_in_box_code).flag = (lat_f.at(C_max_level).at(cell_in_box_code).flag == TBD) ? FLUID : (lat_f.at(C_max_level).at(cell_in_box_code).flag);
                                    }
                                }
                            }
                            
                        } // for k
                    } // for j
                } // for i

            } // for count[3]

        } // for iTri
    } // for ishape

    for (auto i_srf : add_surf) {
        lat_f.at(C_max_level).erase(i_srf);
        lat_sf.insert(make_pair(i_srf, Cell(Cell_Flag::SURFACE)));
    }
}


/**
 * @brief Update the fluid lattice by the means of BFS (Breadth First Search)
 * 
 */
void Lat_Manager::propagate_fluid_lattices()
{
    // Manual BFS implementation
    std::vector<D_morton> queue(lat_f.at(C_max_level).size());
    D_map_define<bool> visited;
    D_int front = -1, rear = -1;

    // Initialize visited map and find starting points
    for (auto& lat_pair : lat_f.at(C_max_level)) {
        visited.insert(std::make_pair(lat_pair.first, false));
        // Add FLUID cells as starting points for BFS
        if (lat_pair.second.flag == FLUID) {
            enqueue(lat_pair.first, front, rear, queue, visited);
        }
    }

    while (front <= rear) {
        D_morton current_code = dequeue(front, rear, queue);
        lat_f.at(C_max_level).at(current_code).flag = (lat_f.at(C_max_level).at(current_code).flag == TBD) ? FLUID : (lat_f.at(C_max_level).at(current_code).flag);

        D_morton ngbr_code[C_Q-1];
        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            ngbr_code[i_q-1] = Morton_Assist::find_neighbor(current_code, C_max_level, ex[i_q], ey[i_q], ez[i_q]);
        }

        for (D_int i_q = 0; i_q < C_Q-1; ++i_q)
        {
            if (triFaceID_in_srf.find(ngbr_code[i_q]) == triFaceID_in_srf.end()) // no TriFace inside the lattice of ngbr_code[i_q]
            {
                auto ngbr_ptr = lat_f.at(C_max_level).find(ngbr_code[i_q]);
                if (ngbr_ptr != lat_f.at(C_max_level).end()) { // ngbr_code[i_q] is fluid
                    if (ngbr_ptr->second.flag == TBD) {
                        enqueue(ngbr_code[i_q], front, rear, queue, visited);
                    }
                } else {
                    bool in_lat_sf  = lat_sf.find(ngbr_code[i_q]) == lat_sf.end();
                    bool in_lat_C2F = lat_overlap_C2F.at(C_max_level).find(ngbr_code[i_q]) == lat_overlap_C2F.at(C_max_level).end();
                    bool in_lat_F2C = lat_overlap_F2C.at(C_max_level).find(ngbr_code[i_q]) == lat_overlap_F2C.at(C_max_level).end();
                    
                    if (in_lat_sf && in_lat_C2F && in_lat_F2C) {
                        std::stringstream error;
                        error << "[propagate_fluid_lattices] lat_sf NO NEIGHBOR!" << std::endl;
                        log_error(error.str(), Log_function::logfile);
                    }
                }

           } else {
                // compute distance between FLUID lattice and SURFACE lattice
                D_vec local_coord, ngbr_coord;
                Morton_Assist::compute_coordinate(current_code, C_max_level, local_coord.x, local_coord.y, local_coord.z);
                Morton_Assist::compute_coordinate(ngbr_code[i_q], C_max_level, ngbr_coord.x, ngbr_coord.y, ngbr_coord.z);
                local_coord += D_vec(dx.at(C_max_level)/2., dx.at(C_max_level)/2., dx.at(C_max_level)/2.);
                ngbr_coord += D_vec(dx.at(C_max_level)/2., dx.at(C_max_level)/2., dx.at(C_max_level)/2.);

                D_real dist_to_latCenter;

                for (auto i_triFace : triFaceID_in_srf.at(ngbr_code[i_q]))
                {
                    D_int shape_id = i_triFace[0]; 
                    D_int triFace_id = i_triFace[1];
                    if (Shape::intersect_line_with_triangle(local_coord, ngbr_coord, 
                         Solid_Manager::pointer_me->shape_solids.at(shape_id).triFace.at(triFace_id))) 
                         {
                            if (dis.find(current_code) == dis.end()) {
                                dis.insert(make_pair(current_code, std::array<D_real, C_Q-1>()));
                                dis.at(current_code).fill(-1.);
                                dis.at(current_code)[i_q] = dist_to_latCenter / Shape::two_points_length(local_coord,ngbr_coord);
                            } else {
                                dis.at(current_code)[i_q] = min_of_two(dis.at(current_code)[i_q], dist_to_latCenter/Shape::two_points_length(local_coord,ngbr_coord));
                            }
                         }
                         
                } // for triFace intersect with the lattice

            }

        } // for i_q
    }
}


/**
 * @brief push back a cell to the BFS queue
 */
void Lat_Manager::enqueue(D_morton cur_code, D_int &front, D_int &rear, std::vector<D_morton> &queue, D_map_define<bool> &visited)
{
    if (rear == lat_f.at(C_max_level).size()-1)
        return;  // queue is full
    
    if (front == -1) {
        front = 0;
    }
    // ++rear;

    if (!visited.at(cur_code)) {
        ++rear;
        queue[rear] = cur_code;
        visited.at(cur_code) = true;
    }
}


/**
 * @brief remove and return the head cell from the queue
 */
D_morton Lat_Manager::dequeue(D_int &front, const D_int &rear, const std::vector<D_morton> &queue)
{
    if (front == -1 || front > rear) {
        return lat_f.at(C_max_level).begin()->first; // queue is empty
    }
    ++front;
    return queue[front-1];
}


void Lat_Manager::update_solid_lattices()
{
    D_setLat add_solid;

    for (auto &i_lat : lat_f.at(C_max_level))
    {
        D_morton cur_code = i_lat.first;

        if (i_lat.second.flag == TBD) {
            add_solid.insert(i_lat.first);
            continue;
        }
        
        if (i_lat.second.flag == IB || i_lat.second.flag == SURFACE) 
        {
            bool fluidNgbr = false;

            D_morton ngbr_code[26];
            ngbr_code[0] = Morton_Assist::find_x0(cur_code, C_max_level); // x0
            ngbr_code[1] = Morton_Assist::find_x1(cur_code, C_max_level); // x1
            ngbr_code[2] = Morton_Assist::find_y0(cur_code, C_max_level); // y0
            ngbr_code[3] = Morton_Assist::find_y1(cur_code, C_max_level); // y1
            ngbr_code[4] = Morton_Assist::find_z0(cur_code, C_max_level); // z0
            ngbr_code[5] = Morton_Assist::find_z1(cur_code, C_max_level); // z1
            
            ngbr_code[6] = Morton_Assist::find_y0(ngbr_code[0], C_max_level); // x0y0
            ngbr_code[7] = Morton_Assist::find_y1(ngbr_code[0], C_max_level); // x0y1
            ngbr_code[8] = Morton_Assist::find_y0(ngbr_code[1], C_max_level); // x1y0
            ngbr_code[9] = Morton_Assist::find_y1(ngbr_code[1], C_max_level); // x1y1

            ngbr_code[10]= Morton_Assist::find_z0(ngbr_code[0], C_max_level); // x0z0
            ngbr_code[11]= Morton_Assist::find_z1(ngbr_code[0], C_max_level); // x0z1
            ngbr_code[12]= Morton_Assist::find_z0(ngbr_code[1], C_max_level); // x1z0
            ngbr_code[13]= Morton_Assist::find_z1(ngbr_code[1], C_max_level); // x1z1

            ngbr_code[14]= Morton_Assist::find_z0(ngbr_code[2], C_max_level); // y0z0
            ngbr_code[15]= Morton_Assist::find_z1(ngbr_code[2], C_max_level); // y0z1
            ngbr_code[16]= Morton_Assist::find_z0(ngbr_code[3], C_max_level); // y1z0
            ngbr_code[17]= Morton_Assist::find_z1(ngbr_code[3], C_max_level); // y1z1

            ngbr_code[18]= Morton_Assist::find_z0(ngbr_code[6], C_max_level); // x0y0z0
            ngbr_code[19]= Morton_Assist::find_z1(ngbr_code[6], C_max_level); // x0y0z1
            ngbr_code[20]= Morton_Assist::find_z0(ngbr_code[7], C_max_level); // x0y1z0
            ngbr_code[21]= Morton_Assist::find_z1(ngbr_code[7], C_max_level); // x0y1z1
            ngbr_code[22]= Morton_Assist::find_z0(ngbr_code[8], C_max_level); // x1y0z0
            ngbr_code[23]= Morton_Assist::find_z1(ngbr_code[8], C_max_level); // x1y0z1
            ngbr_code[24]= Morton_Assist::find_z0(ngbr_code[9], C_max_level); // x1y1z0
            ngbr_code[25]= Morton_Assist::find_z1(ngbr_code[9], C_max_level); // x1y1z1

            for (D_int i_dir = 0; i_dir < 26; ++i_dir) {
                if (lat_f.at(C_max_level).find(ngbr_code[i_dir]) != lat_f.at(C_max_level).end()) 
                {
                    if (lat_f.at(C_max_level).at(ngbr_code[i_dir]).flag == FLUID) {
                        fluidNgbr = true;
                        break;
                    }
                } else {
                    bool in_lat_sf = lat_sf.find(ngbr_code[i_dir]) == lat_sf.end();
                    bool in_lat_C2F = lat_overlap_C2F.at(C_max_level).find(ngbr_code[i_dir]) == lat_overlap_C2F.at(C_max_level).end();
                    bool in_lat_F2C = lat_overlap_F2C.at(C_max_level).find(ngbr_code[i_dir]) == lat_overlap_F2C.at(C_max_level).end();
                    
                    if (in_lat_sf && in_lat_C2F && in_lat_F2C) {
                        std::stringstream error;
                        error << "[update_solid_lattices] lat_sf NO NEIGHBOR!" << std::endl;
                        log_error(error.str(), Log_function::logfile);
                    }
                }
                
            } // for all neighbor

            if (!fluidNgbr) {
                add_solid.insert(i_lat.first);
            }

        }
    } // for i_lat

    for (auto i_solid : add_solid) {
        lat_f.at(C_max_level).erase(i_solid);
        lat_sf.insert(make_pair(i_solid, Cell(Cell_Flag::SOLID)));
    }
}


void Lat_Manager::compute_distance(Solid_Face triFace, D_morton local_code, D_vec lat_center)
{
    double halfbox = dx.at(C_max_level) / 2;
	for (D_int i_q = 1; i_q < C_Q; ++i_q) 
    {
		D_vec ngbr_coord = lat_center + D_vec(ex[i_q]*halfbox, ey[i_q]*halfbox, ez[i_q]*halfbox);
        D_real dist_to_latCenter;
        if (Shape::intersect_line_with_triangle(lat_center, ngbr_coord, triFace, dist_to_latCenter)) {
            // dis.at(local_code)[i_q] = min_of_two(dis.at(local_code)[i_q], dist_to_latCenter / Shape::two_points_length(lat_center, ngbr_coord));
            if (dis.find(local_code) == dis.end()) {
                dis.insert(make_pair(local_code, std::array<D_real, C_Q-1>()));
                dis.at(local_code).fill(-1.);
                dis.at(local_code)[i_q] = dist_to_latCenter / Shape::two_points_length(lat_center, ngbr_coord);
            } else {
                dis.at(local_code)[i_q] = 
                    min_of_two(dis.at(local_code)[i_q], 
                     dist_to_latCenter / Shape::two_points_length(lat_center, ngbr_coord));
            }
        }
	}

}