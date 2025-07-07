/**
* @file
* @brief This class used to manager IO classes
*.
*/
#ifndef IO_MANAGER_H
#define IO_MANAGER_H
#include "General.h"
#include "tecplot.h"
//#include "CGNS.h"
// #include "VTK.h"
// #include "io/HDF5H.h"
// #include "HDF5H.h"  // HDF5 mesh output disabled
class IO_Manager
{
public:
	static IO_Manager* pointer_me;    ///< static pointer point to the IO manager
	unsigned int method;              ///< format used to write flowfield data
	std::string outfile;              ///< name of the file to write flowfield data
	std::vector<unsigned int> vlevel; ///< refiement levels of blocks for flowfield output
public:
	void control();                   ///< the maind funtion of the IO manager
	void writeFlow(const D_int iter);
	void writeMesh();                 ///< function to write mesh in Tecplot format
private:
	IO_TECPLOTH io_tecplot;
//	CGNS io_cgns;
	// VTK io_vtk;
	// HDF5H io_hdf5;  // HDF5 mesh output disabled
};
#endif