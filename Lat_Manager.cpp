/**
 * @file Lat_Manager.cpp
 * @brief 
 * @date 2023-06-03
 * 
 */

// #include "Lat_Manager.hpp"
// #include "Grid_Manager.h"
// #include "Solid_Manager.h"
// #include "Morton_Assist.h"
// #include "General.h"
// #include "Constants.h"
// #include <set>


/**
 * @brief old version to generate lattices, has been deprecated
*/
#ifdef OLD_VERSION_LAT
DEPRECATED void Lat_Manager::initial()
{
    Grid_IB* solidgrid_ptr = &(Grid_Manager::pointer_me->gr_inner);
    Grid_NIB* grid_ptr;

    D_mapNodePtr* bdEast = &(Grid_Manager::pointer_me->bk_boundary_x.at(0).at(1));
    D_mapNodePtr* bdWest = &(Grid_Manager::pointer_me->bk_boundary_x.at(0).at(0));
    D_mapNodePtr* bdNorth = &(Grid_Manager::pointer_me->bk_boundary_y.at(0).at(1));
    D_mapNodePtr* bdSouth = &(Grid_Manager::pointer_me->bk_boundary_y.at(0).at(0));
    D_mapNodePtr* bdTop = &(Grid_Manager::pointer_me->bk_boundary_z.at(0).at(1));
    D_mapNodePtr* bdBot = &(Grid_Manager::pointer_me->bk_boundary_z.at(0).at(0));

    nNodeX = Grid_Manager::pointer_me->nx - 3;
    nNodeY = Grid_Manager::pointer_me->ny - 3;
    nNodeZ = Grid_Manager::pointer_me->nz - 3;

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

        // if (bdWest->find(key_z1)==bdWest->end() || bdWest->find(key_y1)==bdWest->end())
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
        D_morton key_x1 = Morton_Assist::find_x1(key_current, 0);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, 0);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, 0);
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

////////
////////   SIGN
////////

    for (D_int i_level = 0; i_level < C_max_level; ++i_level)
    {
        D_mapLat* lat = &(lat_f.at(i_level));
        grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(i_level));
        dx.at(i_level) = grid_ptr->get_dx();
        for (auto iter : grid_ptr->grid)
        {
            D_morton key_current = iter.first;
            if (ghost_bc.find(key_current)!=ghost_bc.end() && i_level == 0)
                continue;
            
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

            if (bool_xyz) 
            {
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
                }

                else {
                    lat->insert(make_pair(key_current, Cell(Cell_Flag::FLUID)));
                }
            }
        }
    }

    // For inner mesh

#ifndef SPHERE_TEST
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

        if (bool_xyz) {

            if (solidgrid_ptr->grid.at(key_current).flag == Grid_Manager::pointer_me->flag_near_solid && 
                 solidgrid_ptr->grid.at(key_x1).flag == Grid_Manager::pointer_me->flag_near_solid && 
                  solidgrid_ptr->grid.at(key_y1).flag == Grid_Manager::pointer_me->flag_near_solid &&
                   solidgrid_ptr->grid.at(key_z1).flag == Grid_Manager::pointer_me->flag_near_solid &&
                    solidgrid_ptr->grid.at(key_xy1).flag == Grid_Manager::pointer_me->flag_near_solid &&
                     solidgrid_ptr->grid.at(key_xz1).flag == Grid_Manager::pointer_me->flag_near_solid &&
                      solidgrid_ptr->grid.at(key_yz1).flag == Grid_Manager::pointer_me->flag_near_solid &&
                       solidgrid_ptr->grid.at(key_xyz1).flag == Grid_Manager::pointer_me->flag_near_solid) {
                // the lattices within one grid space to the solid node, now to be determined (maybe IB/SURFACE)
                lat_sf.insert(make_pair(key_current,Cell(Cell_Flag::TBD)));
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
            }
            else if (solidgrid_ptr->grid.at(key_current).flag < 0 ||
                      solidgrid_ptr->grid.at(key_x1).flag < 0 ||
                       solidgrid_ptr->grid.at(key_y1).flag < 0 ||
                        solidgrid_ptr->grid.at(key_z1).flag < 0 ||
                         solidgrid_ptr->grid.at(key_xy1).flag < 0 ||
                          solidgrid_ptr->grid.at(key_xz1).flag < 0 ||
                           solidgrid_ptr->grid.at(key_yz1).flag < 0 ||
                            solidgrid_ptr->grid.at(key_xyz1).flag < 0) {
                // the overlap lattices for interpolation between coarse and fine
                lat_overlap_C2F.at(C_max_level).insert(make_pair(key_current,Cell(Cell_Flag::OVERLAP_C2F)));
            }
            else if (solidgrid_ptr->grid.at(key_current).flag == Grid_Manager::pointer_me->flag_near_solid ||
                      solidgrid_ptr->grid.at(key_x1).flag == Grid_Manager::pointer_me->flag_near_solid ||
                       solidgrid_ptr->grid.at(key_y1).flag == Grid_Manager::pointer_me->flag_near_solid ||
                        solidgrid_ptr->grid.at(key_z1).flag == Grid_Manager::pointer_me->flag_near_solid ||
                         solidgrid_ptr->grid.at(key_xy1).flag == Grid_Manager::pointer_me->flag_near_solid ||
                          solidgrid_ptr->grid.at(key_xz1).flag == Grid_Manager::pointer_me->flag_near_solid ||
                           solidgrid_ptr->grid.at(key_yz1).flag == Grid_Manager::pointer_me->flag_near_solid ||
                            solidgrid_ptr->grid.at(key_xyz1).flag == Grid_Manager::pointer_me->flag_near_solid) {
                // part of vertexes of the lattice within one grid space to the solid node, now to be deteremined (maybe IB/SURFACE/FLUID)
                lat_f.at(C_max_level).insert(make_pair(key_current,Cell(Cell_Flag::TBD)));
            }
            else {
                // lattices of pure fluid, no surface inside the lattice
                lat_f.at(C_max_level).insert(make_pair(key_current,Cell(Cell_Flag::FLUID)));
            }
        }
    }
#endif

#ifdef SPHERE_TEST
dx.at(C_max_level) = solidgrid_ptr->get_dx();
double half_diag = sqrt(dx[C_max_level]*dx[C_max_level]/4 * 3);
int lat_f_fluid = 0, lat_f_ib = 0, lat_sf_srf = 0, lat_sf_solid = 0;

for (auto iter : solidgrid_ptr->grid)
{
    D_morton key_current = iter.first;
    D_morton key_x1 = Morton_Assist::find_x1(key_current, C_max_level);
    D_morton key_y1 = Morton_Assist::find_y1(key_current, C_max_level);
    D_morton key_z1 = Morton_Assist::find_z1(key_current, C_max_level);
    D_morton key_xy1 = Morton_Assist::find_y1(key_x1, C_max_level);
    D_morton key_xz1 = Morton_Assist::find_z1(key_x1, C_max_level);
    D_morton key_yz1 = Morton_Assist::find_z1(key_y1, C_max_level);
    D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, C_max_level);

    bool bool_xyz = (solidgrid_ptr->grid.find(key_x1) != solidgrid_ptr->grid.end()) && 
                     (solidgrid_ptr->grid.find(key_y1) != solidgrid_ptr->grid.end()) && 
                      (solidgrid_ptr->grid.find(key_z1) != solidgrid_ptr->grid.end()) && 
                       (solidgrid_ptr->grid.find(key_xy1) != solidgrid_ptr->grid.end()) && 
                        (solidgrid_ptr->grid.find(key_xz1) != solidgrid_ptr->grid.end()) && 
                         (solidgrid_ptr->grid.find(key_yz1) != solidgrid_ptr->grid.end()) && 
                          (solidgrid_ptr->grid.find(key_xyz1) != solidgrid_ptr->grid.end());

    if (!bool_xyz)
        continue;

    if (solidgrid_ptr->grid.at(key_current).flag < 0 ||
         solidgrid_ptr->grid.at(key_x1).flag < 0 ||
          solidgrid_ptr->grid.at(key_y1).flag < 0 ||
           solidgrid_ptr->grid.at(key_z1).flag < 0 ||
            solidgrid_ptr->grid.at(key_xy1).flag < 0 ||
             solidgrid_ptr->grid.at(key_xz1).flag < 0 ||
              solidgrid_ptr->grid.at(key_yz1).flag < 0 ||
               solidgrid_ptr->grid.at(key_xyz1).flag < 0) {
            lat_overlap_C2F.at(C_max_level).insert(make_pair(key_current,Cell(Cell_Flag::OVERLAP_C2F)));
        }
    else 
    {
        D_vec xyz;
        Morton_Assist::compute_coordinate(key_current, C_max_level, xyz.x, xyz.y, xyz.z);
        xyz = xyz + D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);
        D_vec sphere_center(4., 4., 4.);

        // get the (lattice_center-sphere_center) intersect point with the sphere surface
        std::vector<D_vec> intersectPoint_on_sphere;
        if ( get_intersectPoint_withSphere(xyz, sphere_center, intersectPoint_on_sphere) ) 
        {
            std::array<D_vec, 2> lattice_limit;
            lattice_limit[0] = xyz - D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);
            Morton_Assist::compute_coordinate(key_xyz1, C_max_level, lattice_limit[1].x, lattice_limit[1].y, lattice_limit[1].z);

            bool xyz_inside_sphere = Shape::two_points_length(xyz, sphere_center) < 0.5; 
            D_vec real_intersectPoint;

            if (intersectPoint_on_sphere.size() == 1) {
                real_intersectPoint = intersectPoint_on_sphere[0];
            }
            else if (intersectPoint_on_sphere.size() > 1) {
                real_intersectPoint = (Shape::two_points_length(intersectPoint_on_sphere[0], xyz) <= Shape::two_points_length(intersectPoint_on_sphere[1], xyz)) ? intersectPoint_on_sphere[0] : intersectPoint_on_sphere[1];
            }

            
            if (real_intersectPoint.x > lattice_limit[0].x && real_intersectPoint.x < lattice_limit[1].x 
             && real_intersectPoint.y > lattice_limit[0].y && real_intersectPoint.y < lattice_limit[1].y 
             && real_intersectPoint.z > lattice_limit[0].z && real_intersectPoint.z < lattice_limit[1].z) {
                if (xyz_inside_sphere) {
                    lat_sf.insert(make_pair(key_current, Cell(Cell_Flag::SURFACE)));
                }
                else {
                    lat_f.at(C_max_level).insert(make_pair(key_current, Cell(Cell_Flag::IB)));
                }
            }
            else {
                if (xyz_inside_sphere) {
                    lat_sf.insert(make_pair(key_current, Cell(Cell_Flag::SOLID)));
                }
                else {
                    lat_f.at(C_max_level).insert(make_pair(key_current, Cell(Cell_Flag::FLUID)));
                }
            }
        }
    }
}

for (auto i_lat : lat_f.at(C_max_level))
{
    D_vec xyz;
    Morton_Assist::compute_coordinate(i_lat.first, C_max_level, xyz.x, xyz.y, xyz.z);
    xyz = xyz + D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);

    std::array<D_morton, 18> nghbr_point= {Morton_Assist::find_x1(i_lat.first, C_max_level), 
                                           Morton_Assist::find_x0(i_lat.first, C_max_level),
                                           Morton_Assist::find_y1(i_lat.first, C_max_level),
                                           Morton_Assist::find_y0(i_lat.first, C_max_level),
                                           Morton_Assist::find_z1(i_lat.first, C_max_level),
                                           Morton_Assist::find_z0(i_lat.first, C_max_level),

                                           Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_y0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_y0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_y1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),

                                           Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_z0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_z0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_z1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),

                                           Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_z0(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_z0(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level),
                                           Morton_Assist::find_z1(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level)};

    std::array<D_vec, 18> nghbr_xyz;
    for (D_int i_point = 0; i_point < 18; ++i_point) {
        Morton_Assist::compute_coordinate(nghbr_point[i_point], C_max_level, nghbr_xyz[i_point].x, nghbr_xyz[i_point].y, nghbr_xyz[i_point].z);
        nghbr_xyz[i_point] = nghbr_xyz[i_point] + D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);
    }

    for (int i_dir = 0; i_dir < C_Q-1; ++i_dir) 
    {
        auto inbr = nghbr_point[i_dir];
        if (lat_f.at(C_max_level).find(inbr)==lat_f.at(C_max_level).end() && lat_overlap_C2F.at(C_max_level).find(inbr)==lat_overlap_C2F.at(C_max_level).end())
        {
            if (dis.find(i_lat.first) == dis.end())
            {
                std::array<D_real,C_Q-1> init_dis;
                init_dis.fill(NAN);
                dis.insert(make_pair(i_lat.first, init_dis));
            }
            std::vector<D_vec> intersectPoint_on_sphere;
            D_vec real_intersectPoint;
            if ( get_intersectPoint_withSphere(xyz, nghbr_xyz[i_dir], intersectPoint_on_sphere) ) 
            {
                if (intersectPoint_on_sphere.size() == 1)
                    real_intersectPoint = intersectPoint_on_sphere[0];
                else if (intersectPoint_on_sphere.size() > 1) 
                    real_intersectPoint = (Shape::two_points_length(intersectPoint_on_sphere[0], xyz) <= Shape::two_points_length(intersectPoint_on_sphere[1], xyz)) ? intersectPoint_on_sphere[0] : intersectPoint_on_sphere[1];
                
                dis.at(i_lat.first)[i_dir] = Shape::two_points_length(real_intersectPoint,xyz) / Shape::two_points_length(xyz,nghbr_xyz[i_dir]);
                if (lat_sf.find(inbr) == lat_sf.end()) {
                    lat_sf.insert(make_pair(inbr, Cell(Cell_Flag::SOLID)));
                }
            }
        }
    }
}

/*
// for (auto i_lat : lat_f.at(C_max_level))
// {
//     if (i_lat.second.flag == FLUID)
//         continue;

//     D_vec xyz;
//     Morton_Assist::compute_coordinate(i_lat.first, C_max_level, xyz.x, xyz.y, xyz.z);
//     xyz = xyz + D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);
//     std::array<D_morton, 18> nghbr_point= {Morton_Assist::find_x1(i_lat.first, C_max_level), 
//                                            Morton_Assist::find_x0(i_lat.first, C_max_level),
//                                            Morton_Assist::find_y1(i_lat.first, C_max_level),
//                                            Morton_Assist::find_y0(i_lat.first, C_max_level),
//                                            Morton_Assist::find_z1(i_lat.first, C_max_level),
//                                            Morton_Assist::find_z0(i_lat.first, C_max_level),

//                                            Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_y0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_y0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_y1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),

//                                            Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_z0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_z0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_z1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),

//                                            Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_z0(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_z0(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level),
//                                            Morton_Assist::find_z1(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level)};

//     std::array<D_vec, 18> nghbr_xyz;
//     for (D_int i_point = 0; i_point < 18; ++i_point) {
//         Morton_Assist::compute_coordinate(nghbr_point[i_point], C_max_level, nghbr_xyz[i_point].x, nghbr_xyz[i_point].y, nghbr_xyz[i_point].z);
//         nghbr_xyz[i_point] = nghbr_xyz[i_point] + D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);
//     }

//     std::array<D_real,C_Q-1> init_dis;
//     init_dis.fill(NAN);
//     dis.insert(make_pair(i_lat.first, init_dis));

//     for (D_int i_dir = 0; i_dir < 6; ++i_dir) {
//         double dist_ = distance_to_sphere(xyz, nghbr_xyz[i_dir]);
//         if (!isnan(dist_)) {
//             if (dist_ <= dx[C_max_level]/2) {
//                 dis.at(i_lat.first)[i_dir] = dist_;
//             }
//         }
//     }

//     for (D_int i_dir = 6; i_dir < C_Q-1; ++i_dir) {
//         double dist_ = distance_to_sphere(xyz, nghbr_xyz[i_dir]);
//         if (!isnan(dist_)) {
//             if (dist_ <= sqrt(dx[C_max_level]*dx[C_max_level]/2)) {
//                 dis.at(i_lat.first)[i_dir] = dist_;
//             }
//         }
//     }
    
// }

// for (auto i_lat : lat_sf)
// {
//     if (i_lat.second.flag == SURFACE) {

//         D_vec xyz;
//         Morton_Assist::compute_coordinate(i_lat.first, C_max_level, xyz.x, xyz.y, xyz.z);
//         xyz = xyz + D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);

//         std::array<D_morton, C_Q-1> nghbr_point = {
//             Morton_Assist::find_x0(i_lat.first, C_max_level), 
//             Morton_Assist::find_x1(i_lat.first, C_max_level),
//             Morton_Assist::find_y0(i_lat.first, C_max_level), 
//             Morton_Assist::find_y1(i_lat.first, C_max_level), 
//             Morton_Assist::find_z0(i_lat.first, C_max_level), 
//             Morton_Assist::find_z1(i_lat.first, C_max_level), 

//             Morton_Assist::find_y0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_y1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_y0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),

//             Morton_Assist::find_z0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level), 
//             Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_z1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_z0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),

//             Morton_Assist::find_z0(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_z1(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level),
//             Morton_Assist::find_z0(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level)
//         };

//         std::array<D_vec, 18> nghbr_xyz;
//         for (D_int i_point = 0; i_point < C_Q-1; ++i_point) {
//             Morton_Assist::compute_coordinate(nghbr_point[i_point], C_max_level, nghbr_xyz[i_point].x, nghbr_xyz[i_point].y, nghbr_xyz[i_point].z);
//             nghbr_xyz[i_point] = nghbr_xyz[i_point] + D_vec(dx[C_max_level]/2, dx[C_max_level]/2, dx[C_max_level]/2);
//         }

//         for (int i_dir = 0; i_dir < 6; ++i_dir) {
//             if (lat_f.at(C_max_level).find(nghbr_point[i_dir]) == lat_f.at(C_max_level).end())
//                 continue;
            
//             if (lat_f.at(C_max_level).at(nghbr_point[i_dir]).flag != FLUID)
//                 continue;

//             double dist_ = distance_to_sphere(nghbr_xyz[i_dir], xyz);
//             if (dist_ <= dx[C_max_level]) 
//             {
//                 if (dis.find(nghbr_point[i_dir]) == dis.end()) {
//                     std::array<D_real,C_Q-1> init_dis;
//                     init_dis.fill(NAN);
//                     dis.insert(make_pair(nghbr_point[i_dir], init_dis));
//                     dis.at(nghbr_point[i_dir])[inv[i_dir+1]-1] = dist_;
//                 }
//                 else {
//                     dis.at(nghbr_point[i_dir])[inv[i_dir+1]-1] = dist_;
//                 }
//             }
//         }
//         for (int i_dir = 6; i_dir < C_Q-1; ++i_dir) {
//             if (lat_f.at(C_max_level).find(nghbr_point[i_dir]) == lat_f.at(C_max_level).end())
//                 continue;
//             else if (lat_f.at(C_max_level).at(nghbr_point[i_dir]).flag != FLUID)
//                 continue;

//             // std::cout << "dir " << i_dir << " type " << lat_f.at(C_max_level).at(nghbr_point[i_dir]).flag <<
//             //     " (" << nghbr_xyz[i_dir].x << " , " << nghbr_xyz[i_dir].y << " , " << nghbr_xyz[i_dir].z << ")" <<
//             //     sqrt((nghbr_xyz[i_dir].x-4.) * (nghbr_xyz[i_dir].x-4.) + (nghbr_xyz[i_dir].y-4.) * (nghbr_xyz[i_dir].y-4.) + (nghbr_xyz[i_dir].z-4.) * (nghbr_xyz[i_dir].z-4.)) << std::endl;

//             double dist_ = distance_to_sphere(nghbr_xyz[i_dir], xyz);
//             // std::cout << "dist_ " << dist_ << std::endl;
//             if (dist_ <= sqrt(2*dx[C_max_level]*dx[C_max_level])) 
//             {
//                 if (dis.find(nghbr_point[i_dir]) == dis.end()) {
//                     std::array<D_real,C_Q-1> init_dis;
//                     init_dis.fill(NAN);
//                     dis.insert(make_pair(nghbr_point[i_dir], init_dis));
//                     dis.at(nghbr_point[i_dir])[inv[i_dir+1]-1] = dist_;
//                 }
//                 else {
//                     dis.at(nghbr_point[i_dir])[inv[i_dir+1]-1] = dist_;
//                 }
//             }
//         }
//     }
// }
// std::cout << "dis.size() " << dis.size() << std::endl;

*/

#endif

    // generate lat_overlap_F2C
    for (D_int fine_level = 1; fine_level < C_max_level+1; ++fine_level) {
        D_int coarse_level = fine_level - 1;
        for (auto fine_overlap_lat : lat_overlap_C2F.at(fine_level)) {
            if (lat_f.at(coarse_level).find(fine_overlap_lat.first) != lat_f.at(coarse_level).end()) {
                lat_overlap_F2C.at(coarse_level).insert(make_pair(fine_overlap_lat.first, Cell(Cell_Flag::OVERLAP_F2C)));
            }
        }
    }

    // Remove overlap_f2c from lat_f[C_max_level] 
    //  (overlap_f2c collision independently, and stream implicitly by coalescence)
    for (D_int coarse_level = 0; coarse_level < C_max_level; ++coarse_level) {
        for (auto overlap_f2c : lat_overlap_F2C.at(coarse_level)) {
            lat_f.at(coarse_level).erase(overlap_f2c.first);
        }
    }

    // int intp_count = 0;
    // for (auto i_lat : lat_f.at(C_max_level)) {

    //     std::array<D_morton, C_Q-1> nghbr_point = {
    //         Morton_Assist::find_x0(i_lat.first, C_max_level), 
    //         Morton_Assist::find_x1(i_lat.first, C_max_level),
    //         Morton_Assist::find_y0(i_lat.first, C_max_level), 
    //         Morton_Assist::find_y1(i_lat.first, C_max_level), 
    //         Morton_Assist::find_z0(i_lat.first, C_max_level), 
    //         Morton_Assist::find_z1(i_lat.first, C_max_level), 

    //         Morton_Assist::find_y0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_y1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_y0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),

    //         Morton_Assist::find_z0(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level), 
    //         Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_z1(Morton_Assist::find_x0(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_z0(Morton_Assist::find_x1(i_lat.first, C_max_level),C_max_level),

    //         Morton_Assist::find_z0(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_z1(Morton_Assist::find_y0(i_lat.first, C_max_level),C_max_level),
    //         Morton_Assist::find_z0(Morton_Assist::find_y1(i_lat.first, C_max_level),C_max_level)
    //     };
        
    //     for (auto i_nbr : nghbr_point) {
    //         if ((lat_f.at(C_max_level).find(i_nbr) == lat_f.at(C_max_level).end()) &&
    //              (lat_overlap_C2F.at(C_max_level).find(i_nbr) == lat_overlap_C2F.at(C_max_level).end()) ) {
    //                 ++intp_count;
    //                 break;
    //              }
    //     }
    // }

    // std::cout << "intp_count " << intp_count << std::endl;

    auto& lat_f_alias = lat_f;
    auto& overlap_c2f_alias = lat_overlap_C2F;

    // Generate ghost_overlap_C2F
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

            add_ghost_c2f(Morton_Assist::find_x0(c2f_code, fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_x1(c2f_code, fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_y0(c2f_code, fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_y1(c2f_code, fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z0(c2f_code, fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z1(c2f_code, fine_level), ghost_overlap_C2F);

            add_ghost_c2f(Morton_Assist::find_y0(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_y1(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_y0(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_y1(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);

            add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_x0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_x1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);

            add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_y0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_y0(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z0(Morton_Assist::find_y1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
            add_ghost_c2f(Morton_Assist::find_z1(Morton_Assist::find_y1(c2f_code,fine_level), fine_level), ghost_overlap_C2F);
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

            add_ghost_f2c(Morton_Assist::find_x0(f2c_code, coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_x1(f2c_code, coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_y0(f2c_code, coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_y1(f2c_code, coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z0(f2c_code, coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z1(f2c_code, coarse_level), ghost_overlap_F2C);

            add_ghost_f2c(Morton_Assist::find_y0(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_y1(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_y0(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_y1(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);

            add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_x0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_x1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);

            add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_y0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_y0(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z0(Morton_Assist::find_y1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
            add_ghost_f2c(Morton_Assist::find_z1(Morton_Assist::find_y1(f2c_code, coarse_level), coarse_level), ghost_overlap_F2C);
        }
    }

    // for (auto it : lat_sf) {
    //     if (it.second.flag == SURFACE) {
    //         std::cout << "it.second.flag" << std::endl;
    //     }
    // }

    for (int i_level = 0; i_level < lat_overlap_C2F.size(); ++i_level) {
        std::cout << "level " << i_level << " Overlap_C2F " << lat_overlap_C2F.at(i_level).size() << std::endl;
        std::cout << "level " << i_level << " Overlap_F2C " << lat_overlap_F2C.at(i_level).size() << std::endl;
        std::cout << "level " << i_level << " ghost_overlap_C2F " << ghost_overlap_C2F.at(i_level).size() << std::endl;
        std::cout << "level " << i_level << " ghost_overlap_F2C " << ghost_overlap_F2C.at(i_level).size() << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    D_map_define<std::vector<std::set<uint>>> triFaceIdx_in_lat; ///> Store which triFaces the Lattices contain [Lattice MortonCode within shapePoints][Shape][TriFaceID]
    // initial_triface_in_Lat(triFaceIdx_in_lat_sf);
    std::vector<Solid_Face>* triFace_ptr = &(Solid_Manager::pointer_me->shape_solids.at(0).triFace);
    double halfbox[3] = {solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2};

/*
    // for (auto iLatSf : lat_sf) {
    //     if (iLatSf.second.flag != SOLID)
    //         continue;
    //     D_vec lat_center;
    //     Morton_Assist::compute_coordinate(iLatSf.first, C_max_level, lat_center.x, lat_center.y, lat_center.z);
    //     lat_center += D_vec(solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2);
    //     double boxcenter[3] = {lat_center.x, lat_center.y, lat_center.z};
    //
    //     for (auto iFace = triFace_ptr->begin(); iFace != triFace_ptr->end(); ++iFace)
    //     {
    //         double vertex_Coord[3][3] = {{iFace->vertex1.x, iFace->vertex1.y, iFace->vertex1.z},
    //                                          {iFace->vertex2.x, iFace->vertex2.y, iFace->vertex2.z},
    //                                          {iFace->vertex3.x, iFace->vertex3.y, iFace->vertex3.z}};
    //         D_vec face2boxcenter = D_vec(lat_center.x, lat_center.y, lat_center.z) 
    //                             - D_vec((iFace->vertex1.x + iFace->vertex2.x + iFace->vertex3.x)/3,
    //                                     (iFace->vertex1.y + iFace->vertex2.y + iFace->vertex3.y)/3,
    //                                     (iFace->vertex1.z + iFace->vertex2.z + iFace->vertex3.z)/3);
    //         if (Shape::triBoxOverlap(boxcenter, halfbox, vertex_Coord))
    //         {
    //             std::cout << "solid lattice within TriFace" << std::endl;
    //         }
    //     }
    // }
*/

    int ib_in_latf = 0, surface_in_latf = 0, fluid_in_latsf = 0;

#ifndef SPHERE_TEST
    for (auto fluid_inner = lat_f.at(C_max_level).begin(); fluid_inner != lat_f.at(C_max_level).end(); ++fluid_inner) 
    {
        if (fluid_inner->second.flag == TBD)
        {
            D_vec lat_center;
            Morton_Assist::compute_coordinate(fluid_inner->first, C_max_level, lat_center.x, lat_center.y, lat_center.z);
            lat_center += D_vec(halfbox[0], halfbox[1], halfbox[2]);

            for (auto iFace = triFace_ptr->begin(); iFace != triFace_ptr->end(); ++iFace) 
            {
                double boxcenter[3] = {lat_center.x, lat_center.y, lat_center.z};
                // double vertex_Coord[3][3] = {{iFace->vertex1.x, iFace->vertex1.y, iFace->vertex1.z},
                //                              {iFace->vertex2.x, iFace->vertex2.y, iFace->vertex2.z},
                //                              {iFace->vertex3.x, iFace->vertex3.y, iFace->vertex3.z}};
                double vertex_Coord[3][3] = {{iFace->vertex1.x, iFace->vertex1.y, iFace->vertex1.z},
                                             {iFace->vertex2.x, iFace->vertex2.y, iFace->vertex2.z},
                                             {iFace->vertex3.x, iFace->vertex3.y, iFace->vertex3.z}};
                D_vec face2boxcenter = D_vec(lat_center.x, lat_center.y, lat_center.z) 
                                    - D_vec((iFace->vertex1.x + iFace->vertex2.x + iFace->vertex3.x)/3,
                                            (iFace->vertex1.y + iFace->vertex2.y + iFace->vertex3.y)/3,
                                            (iFace->vertex1.z + iFace->vertex2.z + iFace->vertex3.z)/3);

                if (Shape::triBoxOverlap(boxcenter, halfbox, vertex_Coord)) 
                {
                    if (triFaceIdx_in_lat.find(fluid_inner->first) != triFaceIdx_in_lat.end()) {
                        triFaceIdx_in_lat.at(fluid_inner->first).at(0).insert(std::distance(triFace_ptr->begin(),iFace));
                    }
                    else {
                        std::vector<std::set<uint>> a(Solid_Manager::pointer_me->numb_solids);
                        a.at(0).insert(std::distance(triFace_ptr->begin(),iFace));
                        triFaceIdx_in_lat.insert(make_pair(fluid_inner->first, a));
                    }
                    
                    bool is_exceed_centerOfLat = Shape::judge_point_within_Shape(lat_center, *iFace);
                    if (is_exceed_centerOfLat) {
                        // Solid body immersing over lattice center : SURFACE
                        fluid_inner->second.flag = SURFACE;
                        ++surface_in_latf;
                    }

                    // multiple TriFaces cutting the lattice, once the lattice is signed by surface, it should be SURFACE
                    if (!is_exceed_centerOfLat && fluid_inner->second.flag!=SURFACE) {
                        // Solid body occuping a corner of lattice
                        fluid_inner->second.flag = IB;
                        ++ib_in_latf;
                    }
                }

            } // for iFace
        }
        if (fluid_inner->second.flag==TBD)
            fluid_inner->second.flag = FLUID;
    }

    for (auto latsf_inner = lat_sf.begin(); latsf_inner != lat_sf.end(); ++latsf_inner)
    {
        D_vec lat_center;
        Morton_Assist::compute_coordinate(latsf_inner->first, C_max_level, lat_center.x, lat_center.y, lat_center.z);
        lat_center += D_vec(halfbox[0], halfbox[1], halfbox[2]);
        D_uint count_insectPoint = 0;

        for (auto iTriFace = triFace_ptr->begin(); iTriFace != triFace_ptr->end(); ++iTriFace)
        {
            double boxcenter[3] = {lat_center.x, lat_center.y, lat_center.z};
            double vertex_Coord[3][3] = {{iTriFace->vertex1.x, iTriFace->vertex1.y, iTriFace->vertex1.z},
                                         {iTriFace->vertex2.x, iTriFace->vertex2.y, iTriFace->vertex2.z},
                                         {iTriFace->vertex3.x, iTriFace->vertex3.y, iTriFace->vertex3.z}};
            D_vec face2boxcenter = D_vec(lat_center.x, lat_center.y, lat_center.z)
                                    - D_vec((iTriFace->vertex1.x + iTriFace->vertex2.x + iTriFace->vertex3.x)/3,
                                            (iTriFace->vertex1.y + iTriFace->vertex2.y + iTriFace->vertex3.y)/3,
                                            (iTriFace->vertex1.z + iTriFace->vertex2.z + iTriFace->vertex3.z)/3);

            if (Shape::triBoxOverlap(boxcenter, halfbox, vertex_Coord)) {
                if (triFaceIdx_in_lat.find(latsf_inner->first) != triFaceIdx_in_lat.end())
                    triFaceIdx_in_lat.at(latsf_inner->first).at(0).insert(std::distance(triFace_ptr->begin(),iTriFace));
                else {
                    std::vector<std::set<uint>> a(Solid_Manager::pointer_me->numb_solids);
                    a.at(0).insert(std::distance(triFace_ptr->begin(),iTriFace));
                    triFaceIdx_in_lat.insert(make_pair(latsf_inner->first, a));
                }

                bool is_exceed_centerOfLat = Shape::judge_point_within_Shape(lat_center, *iTriFace);
                if (is_exceed_centerOfLat) {
                    // Solid body immersing over lattice center : SURFACE
                    latsf_inner->second.flag = SURFACE;
                    std::cout << "latsf_inner->second.flag = SURFACE" << std::endl;
                }
                
                // multiple TriFaces cutting the lattice, once the lattice is signed by surface, it should be SURFACE
                if (!is_exceed_centerOfLat && latsf_inner->second.flag != SURFACE)
                    latsf_inner->second.flag = IB;
            }

            if (Shape::intersect_line_with_triangle(lat_center, D_vec(Solid_Manager::pointer_me->shape_solids.at(0).x0,Solid_Manager::pointer_me->shape_solids.at(0).x0,Solid_Manager::pointer_me->shape_solids.at(0).x0), *iTriFace))
                ++count_insectPoint;

        } // for iTriFace
        if (latsf_inner->second.flag == TBD) {
            if (count_insectPoint!=0 && count_insectPoint%2 == 0)
                latsf_inner->second.flag = FLUID;
        }
    }

#endif

    // std::cout << "[dis] lat_sf.size() " << lat_sf.size() << std::endl;
    // std::cout << "[dis] ib_in_latf " << ib_in_latf << std::endl;
    // std::cout << "[dis] surface_in_latf " << surface_in_latf << std::endl;
    // for (auto i_lat_sf : lat_sf) {
    //     if (i_lat_sf.second.flag != IB)
    //         std::cout << i_lat_sf.second.flag << std::endl;
    // }
    // std::cout << "[dis] fluid_in_latsf " << fluid_in_latsf << std::endl;

    // D_uint cell_index_ib = 0;
    // D_uint cell_index_fluid = 0;
    // D_uint cell_index_f2c = 0;
    // for (auto fluid_inner = lat_f.at(C_max_level).begin(); fluid_inner != lat_f.at(C_max_level).end(); ++fluid_inner) {
    //     if (fluid_inner->second.flag==IB)
    //         ++cell_index_ib;
    //     else if (fluid_inner->second.flag==FLUID)
    //         ++cell_index_fluid;
    //     else if (fluid_inner->second.flag==OVERLAP_C2F)
    //         ++cell_index_f2c;
    // }

    // std::cout << "lat_sf.size() " << lat_sf.size() << std::endl;
    
    // int latsf_tbd = 0;
    // for (auto latsf_inner = lat_sf.begin(); latsf_inner != lat_sf.end(); ++latsf_inner) {
    //     if (latsf_inner->second.flag == TBD)
    //         latsf_tbd++;
    // }

#ifndef SPHERE_TEST
    D_mapLat::iterator lat_sf_iter = lat_sf.begin();
    while(lat_sf_iter != lat_sf.end()) {
        if (lat_sf_iter->second.flag == FLUID || lat_sf_iter->second.flag == IB) {
            lat_f.at(C_max_level).insert(make_pair(lat_sf_iter->first,Cell(Cell_Flag::lat_sf_iter->second.flag)));
            lat_sf.erase(lat_sf_iter++);
        }
        else {
            ++lat_sf_iter;
        }
    }
#endif


    // int latsf_ib = 0, latsf_surface = 0, latsf_solid = 0;
    // for (auto latsf_inner = lat_sf.begin(); latsf_inner != lat_sf.end(); ++latsf_inner) {
    //     if (latsf_inner->second.flag == TBD)
    //         latsf_tbd++;
    //     else if (latsf_inner->second.flag == IB)
    //         latsf_ib++;
    //     else if (latsf_inner->second.flag == SURFACE)
    //         latsf_surface++;
    //     else if (latsf_inner->second.flag == SOLID)
    //         latsf_solid++;
    //     else {
    //         std::cout << latsf_inner->second.flag << std::endl;
    //     }
    // }
    // std::cout << "latsf_tbd " << latsf_tbd << std::endl;
    // std::cout << "latsf_ib " << latsf_ib << std::endl;
    // std::cout << "latsf_surface " << latsf_surface << std::endl;
    // std::cout << "latsf_solid " << latsf_solid << std::endl;

    // std::cout << lat_sf.size() << std::endl;
    // for (auto ilat_sf : lat_sf) {
    //     D_real xyz_temp[3];
    //     Morton_Assist::compute_coordinate(ilat_sf.first, C_max_level, xyz_temp[0], xyz_temp[1], xyz_temp[2]);
    //     std::cout << "sf xyz " << xyz_temp[0] << " , " << xyz_temp[1] << " , " << " , " << xyz_temp[2] << std::endl;
    // }

    // auto add_latsf = [&lat_f_alias] (D_morton nbgr_code, D_mapLat& latsf_)
    // {
    //     if (lat_f_alias.at(C_max_level).find(nbgr_code) == lat_f_alias.at(C_max_level).end())
    //         latsf_.insert(make_pair(nbgr_code, Cell(Cell_Flag::SOLID)));
            
    // };

#ifndef SPHERE_TEST
    auto lat_f_inner = lat_f.at(C_max_level);
    auto& lat_overlap_C2F_alisa = lat_overlap_C2F.at(C_max_level);
    for (auto iLat : lat_f_inner) {

        auto add_lat_sf = [&lat_f_inner, &iLat, &lat_overlap_C2F_alisa] (D_morton ngbr_code, D_mapLat& lat_sf)
        {
            if (lat_f_inner.find(ngbr_code) == lat_f_inner.end() && lat_overlap_C2F_alisa.find(ngbr_code)==lat_overlap_C2F_alisa.end()) {
                Cell_type ngbr_lat_type = (lat_f_inner.at(iLat.first).flag == IB) ? SOLID : SURFACE;
                lat_sf.insert(std::make_pair(ngbr_code, Cell(Cell_Flag::ngbr_lat_type)));
            }
        };

        // if (iLat.second.flag == IB) {

        add_lat_sf(Morton_Assist::find_x1(iLat.first, C_max_level), lat_sf); // Q1 right
        add_lat_sf(Morton_Assist::find_x0(iLat.first, C_max_level), lat_sf); // Q2 left
        add_lat_sf(Morton_Assist::find_y1(iLat.first, C_max_level), lat_sf); // Q3 front
        add_lat_sf(Morton_Assist::find_y0(iLat.first, C_max_level), lat_sf); // Q4 back
        add_lat_sf(Morton_Assist::find_z1(iLat.first, C_max_level), lat_sf); // Q5 top
        add_lat_sf(Morton_Assist::find_z0(iLat.first, C_max_level), lat_sf); // Q6 bottom

        add_lat_sf(Morton_Assist::find_y1(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q7
        add_lat_sf(Morton_Assist::find_y0(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q8
        add_lat_sf(Morton_Assist::find_y0(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q9
        add_lat_sf(Morton_Assist::find_y1(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q10

        add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q11
        add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q12
        add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q13
        add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q14

        add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_y1(iLat.first, C_max_level),C_max_level), lat_sf); //Q15
        add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_y0(iLat.first, C_max_level),C_max_level), lat_sf); //Q16
        add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_y1(iLat.first, C_max_level),C_max_level), lat_sf); //Q17
        add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_y0(iLat.first, C_max_level),C_max_level), lat_sf); //Q18

        // }
    }
#endif

    // std::cout << "Lattice count" << std::endl;
    // std::cout << " - " << "Fluid num:" << cell_index_fluid << std::endl;
    // std::cout << " - " << "F2C num:" << cell_index_f2c << std::endl;
    // std::cout << " - " << "Surface num: " << lat_sf.size() << std::endl;
    // std::cout << " - " << "cell_index_ib: " << cell_index_ib << std::endl;

    // Remove the case that "Solid body does not pass through the lattice center" from lat_sf,
    //   and insert it as lat_f.
    // For example, the Solid body only occupies a corner of the lattice.
/*
    // for (auto latsf_inner : lat_sf)
    // {
    //     D_vec lat_center;
    //     Morton_Assist::compute_coordinate(latsf_inner.first, C_max_level, lat_center.x, lat_center.y, lat_center.z);
    //     lat_center += D_vec(solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2);
    //     bool break_flag = false;

    //     for (uint i_shape = 0; i_shape < Solid_Manager::pointer_me->numb_solids && !break_flag; ++i_shape)
    //     {
    //         std::vector<Solid_Face>* triFace_ptr = &(Solid_Manager::pointer_me->shape_solids.at(i_shape).triFace);
    //         // if (triFaceIdx_in_lat_sf.find(latsf_inner.first) == triFaceIdx_in_lat_sf.end()) {
    //         //     latsf_inner.second.flag = SOLID;
    //         //     break_flag = true;
    //         //     break;
    //         // }
    //         // else {
    //         //     for (auto triFace_idx : triFaceIdx_in_lat_sf.at(latsf_inner.first).at(i_shape))
    //         //     {
    //         //         if (Shape::judge_point_within_Shape(lat_center, triFace_ptr->at(triFace_idx))) {
    //         //             latsf_inner.second.flag = SURFACE;
    //         //             break_flag = true;
    //         //             break;
    //         //         }
    //         //         else
    //         //             latsf_inner.second.flag = IB;
    //         //     }
    //         // }
    //         for (auto iTriFace : *triFace_ptr)
    //         {
    //             double boxcenter[3] = {lat_center.x, lat_center.y, lat_center.z};
    //             double halfbox[3] = {solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2};
    //             double vertex_Coord[3][3] = {{iTriFace.vertex1.x, iTriFace.vertex1.y, iTriFace.vertex1.z},
    //                                          {iTriFace.vertex2.x, iTriFace.vertex2.y, iTriFace.vertex2.z},
    //                                          {iTriFace.vertex3.x, iTriFace.vertex3.y, iTriFace.vertex3.z}};
    //             if (Shape::triBoxOverlap(boxcenter, halfbox, vertex_Coord)) {

    //             }
    //         } // for triFace

    //     } // for i_shape

    // } // for lat_sf
*/

// test for capture distance lattice
/*
std::vector<Solid_Face>* triFace_ptr = &(Solid_Manager::pointer_me->shape_solids.at(0).triFace);
double halfbox[3] = {solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2};
for (auto fluid_inner = lat_f.at(C_max_level).begin(); fluid_inner != lat_f.at(C_max_level).end(); ++fluid_inner) 
{
    D_vec lat_center;
    Morton_Assist::compute_coordinate(fluid_inner->first, C_max_level, lat_center.x, lat_center.y, lat_center.z);
    lat_center += D_vec(solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2);

    for (auto iFace = triFace_ptr->begin(); iFace != triFace_ptr->end() && fluid_inner->second.flag!=SOLID; ++iFace) 
    {
        double boxcenter[3] = {lat_center.x, lat_center.y, lat_center.z};
        double vertex_Coord[3][3] = {{iFace->vertex1.x, iFace->vertex1.y, iFace->vertex1.z},
                                      {iFace->vertex2.x, iFace->vertex2.y, iFace->vertex2.z},
                                       {iFace->vertex3.x, iFace->vertex3.y, iFace->vertex3.z}};
        D_vec face2boxcenter = D_vec(lat_center.x, lat_center.y, lat_center.z) 
                              - D_vec((iFace->vertex1.x + iFace->vertex2.x + iFace->vertex3.x)/3,
                                       (iFace->vertex1.y + iFace->vertex2.y + iFace->vertex3.y)/3,
                                        (iFace->vertex1.z + iFace->vertex2.z + iFace->vertex3.z)/3);

        if (Shape::triBoxOverlap(boxcenter, halfbox, vertex_Coord)) {
            if (Shape::judge_point_within_Shape(lat_center, *iFace))
                // Solid body immersing over lattice center : SOLID
                fluid_inner->second.flag = SOLID;
            else 
                // Solid body occuping a corner of lattice
                fluid_inner->second.flag = IB;
        }
        if (!Shape::triBoxOverlap(boxcenter, halfbox, vertex_Coord) && fluid_inner->second.flag!=SOLID && fluid_inner->second.flag!=IB)
            // outside the solid body
            fluid_inner->second.flag = FLUID;
    } // for iFace
}

int count_immersingSolid = 0, count_ib = 0, count_insideSolid = 0, count_fluid = 0;

for (auto fluid_inner : lat_f.at(C_max_level)) {
    switch (fluid_inner.second.flag)
    {
        case FLUID:
        {
            ++count_fluid;
            break;
        }
        
        case IB:
        {
            ++count_ib;
            break;
        }

        case SOLID:
        {
            ++count_immersingSolid;
            break;
        }
        
        default:
        {
            log_error("Not FLUID / IB / SOLID.", Log_function::logfile);
            break;
        }
    }
}
std::cout << "Fluid_inner" << std::endl;
std::cout << " - " << "count_immersingSolid: " << count_immersingSolid << std::endl;
std::cout << " - " << "count_ib: " << count_ib << std::endl;
std::cout << " - " << "count_insideSolid: " << count_insideSolid << std::endl;
std::cout << " - " << "count_fluid: " << count_fluid << std::endl;

for (auto srf_inner = lat_sf.begin(); srf_inner != lat_sf.end(); ++srf_inner) 
{
    D_vec lat_center;
    Morton_Assist::compute_coordinate(srf_inner->first, C_max_level, lat_center.x, lat_center.y, lat_center.z);
    lat_center += D_vec(solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2, solidgrid_ptr->get_dx()/2);

    for (auto iFace = triFace_ptr->begin(); iFace != triFace_ptr->end() && srf_inner->second.flag!=INSIDESOLID; ++iFace) 
    {
        double boxcenter[3] = {lat_center.x, lat_center.y, lat_center.z};
        double vertex_Coord[3][3] = {{iFace->vertex1.x, iFace->vertex1.y, iFace->vertex1.z},
                                      {iFace->vertex2.x, iFace->vertex2.y, iFace->vertex2.z},
                                       {iFace->vertex3.x, iFace->vertex3.y, iFace->vertex3.z}};
        D_vec face2boxcenter = D_vec(lat_center.x, lat_center.y, lat_center.z) 
                              - D_vec((iFace->vertex1.x + iFace->vertex2.x + iFace->vertex3.x)/3,
                                       (iFace->vertex1.y + iFace->vertex2.y + iFace->vertex3.y)/3,
                                        (iFace->vertex1.z + iFace->vertex2.z + iFace->vertex3.z)/3);

        if (Shape::triBoxOverlap(boxcenter, halfbox, vertex_Coord)) {
            bool is_solid = Shape::judge_point_within_Shape(lat_center, *iFace);

            if (is_solid)
                // Solid body immersing over lattice center : SOLID
                srf_inner->second.flag = SOLID;

            if (!is_solid && srf_inner->second.flag!=SOLID)
                // Solid body occuping a corner of lattice
                srf_inner->second.flag = IB;
        }
        else
            // within the solid body
            srf_inner->second.flag = INSIDESOLID;
    } // for iFace
}

for (auto srf_inner : lat_sf) {
    switch (srf_inner.second.flag)
    {
        case SOLID:
        {
            ++count_fluid;
            break;
        }
    
        case IB:
        {
            ++count_ib;
            break;
        }

        case INSIDESOLID:
        {
            ++count_immersingSolid;
            break;
        }
    
        default:
        {
            std::cout << srf_inner.second.flag << std::endl;
            log_error("Not SOLID / IB / INSIDESOLID.", Log_function::logfile);
            break;
        }
    }
}

std::cout << "Surface_inner" << std::endl;
std::cout << " - " << "count_immersingSolid: " << count_immersingSolid << std::endl;
std::cout << " - " << "count_ib: " << count_ib << std::endl;
std::cout << " - " << "count_insideSolid: " << count_insideSolid << std::endl;
std::cout << " - " << "count_fluid: " << count_fluid << std::endl;

*/

    // std::cout << "triFaceIdx_in_lat.size() " << triFaceIdx_in_lat.size() << std::endl;
#ifndef SPHERE_TEST
    initial_IBdistance(triFaceIdx_in_lat);
#endif

    for (D_int i_level = C_max_level; i_level >=0 ; --i_level)
    {
        std::cout << "ilevel = " << i_level << ", No. of lattices:" << lat_f.at(i_level).size() << std::endl;
    }
    std::cout << "dis.size() " << dis.size() << std::endl;


    // auto lat_f_inner = lat_f.at(C_max_level);
    // for (auto iLat : dis) {

    //     auto add_lat_sf = [&lat_f_inner, &iLat] (D_morton ngbr_code, D_mapLat& lat_sf)
    //     {
    //         if (lat_f_inner.find(ngbr_code) == lat_f_inner.end()) {
    //             Cell_type ngbr_lat_type = (lat_f_inner.at(iLat.first).flag == IB) ? SOLID : SURFACE;
    //             lat_sf.insert(std::make_pair(ngbr_code, Cell(Cell_Flag::ngbr_lat_type)));
    //         }
    //     };

    //     add_lat_sf(Morton_Assist::find_x1(iLat.first, C_max_level), lat_sf); // Q1 right
    //     add_lat_sf(Morton_Assist::find_x0(iLat.first, C_max_level), lat_sf); // Q2 left
    //     add_lat_sf(Morton_Assist::find_y1(iLat.first, C_max_level), lat_sf); // Q3 front
    //     add_lat_sf(Morton_Assist::find_y0(iLat.first, C_max_level), lat_sf); // Q4 back
    //     add_lat_sf(Morton_Assist::find_z1(iLat.first, C_max_level), lat_sf); // Q5 top
    //     add_lat_sf(Morton_Assist::find_z0(iLat.first, C_max_level), lat_sf); // Q6 bottom

    //     add_lat_sf(Morton_Assist::find_y1(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q7
    //     add_lat_sf(Morton_Assist::find_y0(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q8
    //     add_lat_sf(Morton_Assist::find_y0(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q9
    //     add_lat_sf(Morton_Assist::find_y1(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q10

    //     add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q11
    //     add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q12
    //     add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_x1(iLat.first, C_max_level),C_max_level), lat_sf); //Q13
    //     add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_x0(iLat.first, C_max_level),C_max_level), lat_sf); //Q14

    //     add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_y1(iLat.first, C_max_level),C_max_level), lat_sf); //Q15
    //     add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_y0(iLat.first, C_max_level),C_max_level), lat_sf); //Q16
    //     add_lat_sf(Morton_Assist::find_z0(Morton_Assist::find_y1(iLat.first, C_max_level),C_max_level), lat_sf); //Q17
    //     add_lat_sf(Morton_Assist::find_z1(Morton_Assist::find_y0(iLat.first, C_max_level),C_max_level), lat_sf); //Q18
    // }

    // std::cout << "lat_sf.size() " << lat_sf.size() << std::endl;

    // for (auto i_lat : lat_f.at(C_max_level)) {
    //     if (i_lat.second.flag == SOLID) {
    //         std::cout << "i_lat.second.flag " << i_lat.second.flag << std::endl;
    //     }
    // }

        

}
#endif

/**
 * @brief 
 * @date 2023.6.20
 */

#ifdef OLD_VERSION_LAT
/*
// void Lat_Manager::initial_triface_in_Lat(D_map_define<std::vector<std::set<uint>>> &triFaceIdx_in_lat_sf)
// {
//     D_real inner_dx = Grid_Manager::pointer_me->gr_inner.get_dx();

//     for (uint ishape = 0; ishape < Solid_Manager::pointer_me->numb_solids; ++ishape)
//     {
//         std::array<D_real, C_DIMS> xyz;  ///> The point coordinates of SOLID
// #if (C_CHECK_MORTON_BOUNDARY == 1)
// 		std::array<D_real, C_DIMS> temp_xyz; // the latest xyz which does not exceed the computational domain
// #endif
// 		for (D_uint ipoint = 0; ipoint < Solid_Manager::pointer_me->shape_solids.at(ishape).numb_nodes; ++ipoint)
// 		{
//             uint triFace_idx = Solid_Manager::pointer_me->shape_solids.at(ishape).triFaceIdx_of_node[ipoint];

// 			xyz.at(0) = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).x;
// 			xyz.at(1) = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).y;
// 			xyz.at(2) = Solid_Manager::pointer_me->shape_solids.at(ishape).node.at(ipoint).z;
// #if (C_CHECK_MORTON_BOUNDARY == 1)
// 			bool exceed_domain = false;
// 			if (xyz.at(0) > xb_domain + C_eps)
// 			{
// 				std::stringstream warning;
// 				warning << "coordinate x = " << xyz.at(0) << " of point " << ipoint << " in shape "<< ishape << " exceeds the compuational domain where xb_domain  = " << xb_domain << std::endl;
// 				log_warning(warning.str(), Log_function::logfile);
// 				exceed_domain = true;
// 			}
// 			if (xyz.at(1) > yb_domain + C_eps)
// 			{
// 				std::stringstream warning;
// 				warning << "coordinate y = " << xyz.at(1) << " of point " << ipoint << " in shape " << ishape << " exceeds the compuational domain where yb_domain  = " << yb_domain << std::endl;
// 				log_warning(warning.str(), Log_function::logfile);
// 				exceed_domain = true;
// 			}
// #if(C_DIMS == 2)
// 			if (exceed_domain)
// 			{
// 				xyz = temp_xyz;
// 			}
// 			else
// 			{
// 				temp_xyz = xyz;
// 			}
// #endif
// #if(C_DIMS == 3)
// 			if (xyz.at(2) > zb_domain  + C_eps)
// 			{
// 				std::stringstream warning;
// 				warning << "coordinate z = " << xyz.at(2) << "of point " << ipoint << " in shape " << ishape << " exceeds the compuational domain where zb_domain  = " << zb_domain  << std::endl;
// 				log_warning(warning.str(), Log_function::logfile);
// 				exceed_domain = true;
// 			}
// 			if (exceed_domain)
// 			{
// 				xyz = temp_xyz;
// 			}
// 			else
// 			{
// 				temp_xyz = xyz;
// 			}
// #endif
// #endif
//             D_real x_index = xyz.at(0) / inner_dx;
//             D_real y_index = xyz.at(1) / inner_dx;
//             D_real z_index = xyz.at(2) / inner_dx;
//             D_uint xint = static_cast<D_uint>(x_index);
//             D_uint yint = static_cast<D_uint>(y_index);
//             D_uint zint = static_cast<D_uint>(z_index);
// 			auto lat_morton = Morton_Assist::pointer_me->morton_encode(xint, yint, zint);

//             if (triFaceIdx_in_lat_sf.find(lat_morton) == triFaceIdx_in_lat_sf.end()) {
//                 std::vector<std::set<uint>> a(Solid_Manager::pointer_me->numb_solids);
//                 a.at(ishape).insert(triFace_idx);
//                 triFaceIdx_in_lat_sf.insert(make_pair(lat_morton, a));
//             }
//             else
//                 triFaceIdx_in_lat_sf.at(lat_morton).at(ishape).insert(triFace_idx);
// 		} // for ipoint
//     } // for ishape
    
// }

*/
#endif

#ifdef OLD_VERSION_LAT
void Lat_Manager::initial_IBdistance(const D_map_define<std::vector<std::set<uint>>> &triFaceIdx_in_lat_sf)
{    
    D_real inner_dx = Grid_Manager::pointer_me->gr_inner.get_dx();
    
    D_mapLat* latf_inner = &(lat_f.at(C_max_level));

    for (uint ishape = 0; ishape < Solid_Manager::pointer_me->numb_solids; ++ishape)
    {
        std::vector<Solid_Face>* triFace_ptr = &(Solid_Manager::pointer_me->shape_solids.at(ishape).triFace);
        for (auto surf_lat : triFaceIdx_in_lat_sf)
        {
            D_morton lat_sf_code = surf_lat.first;

            D_vec solid_lat_center;
            Morton_Assist::compute_coordinate(lat_sf_code, C_max_level, solid_lat_center.x, solid_lat_center.y, solid_lat_center.z);
            solid_lat_center += D_vec(inner_dx/2, inner_dx/2, inner_dx/2);

            // surf_lat is a "IB lattice"
            // surf_lat needs a IB_distance
            if (latf_inner->find(lat_sf_code) != latf_inner->end()) 
            {

                if (latf_inner->at(lat_sf_code).flag != IB) {
                    log_error("Lat_Manager error: Lattice containing solid points is signed as FLUID wrongly!",Log_function::logfile);
                }
                
                if (dis.find(lat_sf_code) == dis.end()) {
                    std::array<D_real,C_Q-1> init_dis;
                    init_dis.fill(NAN);
                    dis.insert(make_pair(lat_sf_code, init_dis));
                }

                auto& lat_sf_alias = lat_sf;
                std::cout << lat_sf_alias.size() << std::endl;
                /**
                 * @brief insert to dis
                 * @param[in] dir_no
                 * @param[in] x0_code
                 * @param[out] dis 
                 */
                auto add_dis = [&lat_sf_alias, &surf_lat, ishape, inner_dx, solid_lat_center, &triFace_ptr, lat_sf_code]
                 (int dir_no, D_morton x0_code, D_map_define<std::array<D_real, C_Q-1>> &dis)
                  {
                      if (lat_sf_alias.find(x0_code) != lat_sf_alias.end())
                      {
                          D_vec x0;
                          Morton_Assist::compute_coordinate(x0_code, C_max_level, x0.x, x0.y, x0.z);
                          x0 += D_vec(inner_dx/2, inner_dx/2, inner_dx/2);
  
                          Solid_Node projectPoint_on_triFace;
  
                          for (uint triFace_idx : surf_lat.second.at(ishape))
                          {
                              if (Shape::intersect_line_with_triangle(solid_lat_center, x0, triFace_ptr->at(triFace_idx), projectPoint_on_triFace)) {
                                  dis.at(lat_sf_code)[dir_no] = 
                                   Shape::two_points_length(D_vec(projectPoint_on_triFace.x, projectPoint_on_triFace.y, projectPoint_on_triFace.z), solid_lat_center) /
                                    Shape::two_points_length(x0, solid_lat_center);
                                  break;
                              }
                          }
                      }
                  };
                
                // neighbor 1
                D_morton x1_code = Morton_Assist::find_x1(lat_sf_code, C_max_level);
                add_dis(0, x1_code, dis);
                
                // neighbor 2
                D_morton x0_code = Morton_Assist::find_x0(lat_sf_code, C_max_level);
                add_dis(1, x0_code, dis);

                // neighbor 3
                D_morton y1_code = Morton_Assist::find_y1(lat_sf_code, C_max_level);
                add_dis(2, y1_code, dis);

                // neighbor 4
                D_morton y0_code = Morton_Assist::find_y0(lat_sf_code, C_max_level);
                add_dis(3, y0_code, dis);

                // neighbor 5
                D_morton z1_code = Morton_Assist::find_z1(lat_sf_code, C_max_level);
                add_dis(4, z1_code, dis);

                // neighbor 6
                D_morton z0_code = Morton_Assist::find_z0(lat_sf_code, C_max_level);
                add_dis(5, z0_code, dis);

                // neighbor 7
                D_morton x1y1_code = Morton_Assist::find_y1(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(6 ,x1y1_code, dis);

                // neighbor 8
                D_morton x0y0_code = Morton_Assist::find_y0(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(7, x0y0_code, dis);

                // neighbor 9
                D_morton x1y0_code = Morton_Assist::find_y0(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(8, x1y0_code, dis);

                // neighbor 10
                D_morton x0y1_code = Morton_Assist::find_y1(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(9, x0y1_code, dis);

                // neighbor 11
                D_morton x1z1_code = Morton_Assist::find_z1(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(10, x1z1_code, dis);

                // neighbor 12
                D_morton x0z0_code = Morton_Assist::find_z0(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(11, x0z0_code, dis);

                // neighbor 13
                D_morton x1z0_code = Morton_Assist::find_z0(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(12, x1z0_code, dis);

                // neighbor 14
                D_morton x0z1_code = Morton_Assist::find_z1(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(13, x0z1_code, dis);

                // neighbor 15
                D_morton y1z1_code = Morton_Assist::find_z1(Morton_Assist::find_y1(lat_sf_code, C_max_level), C_max_level);
                add_dis(14, y1z1_code, dis);

                // neighbor 16
                D_morton y0z0_code = Morton_Assist::find_z0(Morton_Assist::find_y0(lat_sf_code, C_max_level), C_max_level);
                add_dis(15, y0z0_code, dis);

                // neighbor 17
                D_morton y1z0_code = Morton_Assist::find_z0(Morton_Assist::find_y1(lat_sf_code, C_max_level), C_max_level);
                add_dis(16, y1z0_code, dis);

                // neighbor 18
                D_morton y0z1_code = Morton_Assist::find_z1(Morton_Assist::find_y0(lat_sf_code, C_max_level), C_max_level);
                add_dis(17, y0z1_code, dis);

            } // endif surf_lat is a IB_Lattice (i.e. fluid but containing solid body)

            // surf_lat is a "SOLID lattice"
            // surf_lat does not need IB_distance, but only provides IB_distance to neighbor FLUID/IB lattice
            if (lat_sf.find(lat_sf_code) != lat_sf.end()) 
            {
                auto add_dis = [&latf_inner, inner_dx, &triFaceIdx_in_lat_sf, ishape, &surf_lat, solid_lat_center, &triFace_ptr]
                 (int dir_no, D_morton nbrx1_code, D_map_define<std::array<D_real, C_Q-1>> &dis)
                  {
                    if (dis.find(nbrx1_code) == dis.end()) {
                        std::array<D_real,C_Q-1> init_dis;
                        init_dis.fill(NAN);
                        dis.insert(make_pair(nbrx1_code, init_dis));
                    }

                    if (latf_inner->find(nbrx1_code) != latf_inner->end())
                    {
                        D_vec x1;
                        Morton_Assist::compute_coordinate(nbrx1_code, C_max_level, x1.x, x1.y, x1.z);
                        x1 += D_vec(inner_dx/2, inner_dx/2, inner_dx/2);
                        if (latf_inner->at(nbrx1_code).flag == IB)
                        {
                            Solid_Node projectPoint_on_triFace;
                            bool intersect_within_nbrx1 = false;
                            for (uint triFace_inIB_idx : triFaceIdx_in_lat_sf.at(nbrx1_code).at(ishape)) {
                                if (Shape::intersect_line_with_triangle(solid_lat_center, x1, triFace_ptr->at(triFace_inIB_idx), projectPoint_on_triFace)) {
                                    intersect_within_nbrx1 = true;
                                    break;
                                }
                            }
                            
                            if (!intersect_within_nbrx1) {
                                for (uint triFace_idx : surf_lat.second.at(ishape)) {
                                    if (Shape::intersect_line_with_triangle(solid_lat_center, x1, triFace_ptr->at(triFace_idx), projectPoint_on_triFace)) {
                                        dis.at(nbrx1_code)[dir_no] = 
                                         Shape::two_points_length(D_vec(projectPoint_on_triFace.x, projectPoint_on_triFace.y, projectPoint_on_triFace.z), x1) / 
                                          Shape::two_points_length(x1, solid_lat_center);
                                        break;
                                    }
                                }
                            }
                        } // endif neighbor if IB
                        else if (latf_inner->at(nbrx1_code).flag == FLUID)
                        {
                            Solid_Node projectPoint_on_triFace;
                            for (uint triFace_idx : surf_lat.second.at(ishape)) {
                                if (Shape::intersect_line_with_triangle(solid_lat_center, x1, triFace_ptr->at(triFace_idx), projectPoint_on_triFace)) {
                                    dis.at(nbrx1_code)[dir_no] = 
                                     Shape::two_points_length(D_vec(projectPoint_on_triFace.x, projectPoint_on_triFace.y, projectPoint_on_triFace.z), x1) / 
                                      Shape::two_points_length(x1, solid_lat_center);
                                    break;
                                }
                            }
                        } // endif neighbor if FLUID

                    } // endif neighbor not Solid

                  };

                // neighbor 1
                D_morton nbrx1_code = Morton_Assist::find_x1(lat_sf_code, C_max_level);
                add_dis(1, nbrx1_code, dis);

                // neighbor 2
                D_morton nbrx0_code = Morton_Assist::find_x0(lat_sf_code, C_max_level);
                add_dis(0, nbrx0_code, dis);

                // neighbor 3
                D_morton nbry1_code = Morton_Assist::find_y1(lat_sf_code, C_max_level);
                add_dis(3, nbry1_code, dis);

                // neighbor 4
                D_morton nbry0_code = Morton_Assist::find_y0(lat_sf_code, C_max_level);
                add_dis(2, nbry0_code, dis);

                // neighbor 5
                D_morton nbrz1_code = Morton_Assist::find_z1(lat_sf_code, C_max_level);
                add_dis(5, nbrz1_code, dis);

                // neighbor 6
                D_morton nbrz0_code = Morton_Assist::find_z0(lat_sf_code, C_max_level);
                add_dis(4, nbrz0_code, dis);

                // neighbor 7
                D_morton nbrx1y1_code = Morton_Assist::find_y1(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(7, nbrx1y1_code, dis);

                // neighbor 8
                D_morton nbrx0y0_code = Morton_Assist::find_y0(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(6, nbrx0y0_code, dis);

                // neighbor 9
                D_morton nbrx1y0_code = Morton_Assist::find_y0(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(9, nbrx1y0_code, dis);

                // neighnor 10
                D_morton nbrx0y1_code = Morton_Assist::find_y1(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(8, nbrx0y1_code, dis);

                // neighbor 11
                D_morton nbrx1z1_code = Morton_Assist::find_z1(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(11, nbrx1z1_code, dis);

                // neighbor 12
                D_morton nbrx0z0_code = Morton_Assist::find_z0(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(10, nbrx0z0_code, dis);

                // neighbor 13
                D_morton nbrx1z0_code = Morton_Assist::find_z0(Morton_Assist::find_x1(lat_sf_code, C_max_level), C_max_level);
                add_dis(13, nbrx1z0_code, dis);

                // neighbor 14
                D_morton nbrx0z1_code = Morton_Assist::find_z1(Morton_Assist::find_x0(lat_sf_code, C_max_level), C_max_level);
                add_dis(12, nbrx0z1_code, dis);

                // neighbor 15
                D_morton nbry1z1_code = Morton_Assist::find_z1(Morton_Assist::find_y1(lat_sf_code, C_max_level), C_max_level);
                add_dis(15, nbry1z1_code, dis);

                // neighbor 16
                D_morton nbry0z0_code = Morton_Assist::find_z0(Morton_Assist::find_y0(lat_sf_code, C_max_level), C_max_level);
                add_dis(14, nbry0z0_code, dis);

                // neighbor 17
                D_morton nbry1z0_code = Morton_Assist::find_z0(Morton_Assist::find_y1(lat_sf_code, C_max_level), C_max_level);
                add_dis(17, nbry1z0_code, dis);

                // neighbor 18
                D_morton nbry0z1_code = Morton_Assist::find_z1(Morton_Assist::find_y0(lat_sf_code, C_max_level), C_max_level);
                add_dis(16, nbry0z1_code, dis);
                
            } // endif surf_lat is a Solid

            // ++surf_lat_index;
        } // for lattice containing triFace

    }

}
#endif

#ifdef SPHERE_TEST
bool Lat_Manager::get_intersectPoint_withSphere(const D_vec& pointA, const D_vec& pointB, std::vector<D_vec>& intersectPoint, const D_vec& sphere_center, double r)
{
    std::vector<D_vec>().swap(intersectPoint);
    D_vec d = pointB - pointA;

    double a = d.x * d.x + d.y * d.y + d.z * d.z;
    double b = 2 * d.x * (pointA.x - sphere_center.x) + 2 * d.y *(pointA.y - sphere_center.y) + 2 *d.z * (pointA.z - sphere_center.z);
    double c = (pointA.x - sphere_center.x)*(pointA.x - sphere_center.x) + (pointA.y - sphere_center.y)*(pointA.y - sphere_center.y) + (pointA.z - sphere_center.z)*(pointA.z - sphere_center.z) - r*r;

    std::vector<double> t;
    double delta = b*b - 4*a*c;
    if (delta < 0)
        return false; // no intersect with sphere
    if (abs(delta) < 1e-10) {
        t.push_back(-b / (2*a));
    }
    else {
        t.push_back((-b + sqrt(delta)) / (2*a));
        t.push_back((-b - sqrt(delta)) / (2*a));
    }
    for (auto it : t) {
        if (it >= 0 && it <= 1) {
            intersectPoint.push_back(pointA + D_vec(d.x*it, d.y*it, d.z*it));
        }
    }
    return true;

}
#endif


// void Lat_Manager::output_hdf5()
// {
    
// }