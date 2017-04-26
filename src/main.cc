/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cc
 * Copyright (C) 2017 Lukas Gartmair <lukas@lgartmair-GA-H81M-D2V>
 * 
 * isosurface is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * isosurface is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "atomprobe/primitives/ionHit.h"
#include "atomprobe/io/dataFiles.h"
#include "atomprobe/io/ranges.h"
#include "../testsuite/src/CTF_functions.h"
#include "../testsuite/src/vdb_functions.h"
#include <string>
#include <openvdb/openvdb.h>
//#include <openvdb/math/Maps.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/io/Stream.h>
#include <openvdb/math/Coord.h>
#include <openvdb/metadata/MetaMap.h>

int main()
{
	std::cout << "Isosurface creation!" << std::endl;

	// load posfile
	std::string pos_filename = "../../data/distributed_ref_pos_atomic_density_25_precRad_8_precConc_10_excessL_25_excessR_25.pos";
	const char * pos_char_filename = pos_filename.c_str();
	std::vector<AtomProbe::IonHit> posIons;
	loadPosFile(posIons, pos_char_filename);

	// load range file
	std::string rng_filename = "../../data/sampled_ranges.rng";
	const char * rng_char_filename = rng_filename.c_str();

	AtomProbe::RangeFile rng;
	bool rng_check = rng.open(rng_char_filename);

	AtomProbe::Point3D position;
	position = posIons[0].getPos();

	// Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
	openvdb::initialize();

	const float background = 0.0;	

	// initialize a grid where the division result is stored
	openvdb::FloatGrid::Ptr calculation_result_grid = openvdb::FloatGrid::create(background);

	// initialize nominator and denominator grids
	openvdb::FloatGrid::Ptr denominator_grid = openvdb::FloatGrid::create(background);
	openvdb::FloatGrid::Accessor denominator_accessor = denominator_grid->getAccessor();

	openvdb::FloatGrid::Ptr numerator_grid = openvdb::FloatGrid::create(background);
	openvdb::FloatGrid::Accessor numerator_accessor = numerator_grid->getAccessor();

	const unsigned int xyzs = 3;

	const float voxelsize = 1.0; //nm

	const unsigned int numerator_id = 1;

	// positions size will have several mio - which for counter type to store that 
	for (double i=0;i<posIons.size();++i)
	{
	
		AtomProbe::Point3D current_position;
		current_position = posIons[i].getPos();

		AtomProbe::IonHit current_hit;
		current_hit = posIons[i];

		unsigned int current_ionID = rng.getIonID(current_hit);

		std::vector<double> atom_position(xyzs); 
		for (unsigned int j=0;j<xyzs;++j)
		{
			atom_position[j] = current_position[j];
		}

		// 1st step - project the current atom position to unit voxel i.e. from 0 to 1
		std::vector<double> position_in_unit_voxel;
		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxelsize);

		// 2nd step - determine each contribution to the adjecent 8 voxels outgoining from the position in the unit voxel
		std::vector<double> volumes_of_subcuboids;
		std::vector<double> contributions_to_adjacent_voxels;
		bool vertex_corner_coincidence = false;

		vertex_corner_coincidence = CTF::checkVertexCornerCoincidence(position_in_unit_voxel);

		// in case of coincidence of atom and voxel the contribution becomes 100 percent
		if (vertex_corner_coincidence == false)
		{
			volumes_of_subcuboids = CTF::calcSubvolumes(position_in_unit_voxel);
			contributions_to_adjacent_voxels = CTF::HellmanContributions(volumes_of_subcuboids);
		}
		else
		{
			contributions_to_adjacent_voxels = CTF::handleVertexCornerCoincidence(position_in_unit_voxel);
		}

		// 3rd step - determine the adjacent voxel indices in the actual grid
		std::vector<std::vector<double> > adjacent_voxel_vertices;
		adjacent_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxelsize);

		// 4th step - assign each of the 8 adjacent voxels the corresponding contribution that results from the atom position in the unit voxel

		const unsigned int number_of_adjacent_voxels = 8;
		std::vector<double> current_voxel_index;
		for (unsigned int k=0;k<number_of_adjacent_voxels;++k)
		{

			current_voxel_index = adjacent_voxel_vertices[k];
			// normalized voxel indices based on 00, 01, 02 etc. // very important otherwise there will be spacings
			openvdb::Coord ijk(current_voxel_index[0], current_voxel_index[1], current_voxel_index[2]);
	
			// write to all ions to the denominator grid
			denominator_accessor.setValue(ijk, contributions_to_adjacent_voxels[k] + denominator_accessor.getValue(ijk));

			// only write matching numerator ids to the numerator grid
			if (current_ionID == numerator_id) 
			{
				numerator_accessor.setValue(ijk, contributions_to_adjacent_voxels[k] + numerator_accessor.getValue(ijk));
			}
		}
	}

	std::cout << " active voxel count numerator_grid " << " = " << numerator_grid->activeVoxelCount() << std::endl;
	std::cout << " active voxel count denominator_grid" << " = " << denominator_grid->activeVoxelCount() << std::endl;

	// final voxelization calculation

	// composite.h operations modify the first grid and leave the second grid emtpy !!!!!!!!!!!!!!!!!!
	// compute a = a / b
	openvdb::tools::compDiv(*numerator_grid, *denominator_grid);

	calculation_result_grid = numerator_grid->deepCopy();

	//check for negative nans and infs introduced by the division
	//set them to zero in order not to obtain nan mesh coordinates

	for (openvdb::FloatGrid::ValueAllIter iter = calculation_result_grid->beginValueAll(); iter; ++iter)
	{   
		if (std::isfinite(iter.getValue()) == false)
		{
    			iter.setValue(0.0);
		}
	}

	// Associate a scaling transform with the grid that sets the voxel size
	// to voxelsize units in world space.

	openvdb::math::Transform::Ptr linearTransform = openvdb::math::Transform::createLinearTransform(voxelsize);
	calculation_result_grid->setTransform(linearTransform);

	std::vector<openvdb::Vec3s> points;
  	std::vector<openvdb::Vec3I> triangles;
  	std::vector<openvdb::Vec4I> quads;	

	float isovalue = 0.2;

	openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*calculation_result_grid, points, triangles, quads, isovalue);

	// check the mesh for nans and infinites and set them to zero
	for(unsigned int ui=0;ui<points.size();ui++)
	{
		for(unsigned int uj=0;uj<xyzs;uj++)
		{
			if (std::isfinite(points[ui][uj]) == false)
			{
				for(unsigned int uk=0;uk<xyzs;uk++)
				{
					points[ui][uk] = 0.0;
				}
			}
		}
	}	

	// create a triangular mesh
	double number_of_splitted_triangles = 2*quads.size();
	std::vector<std::vector<float> > triangles_from_splitted_quads(number_of_splitted_triangles, std::vector<float>(xyzs));

	triangles_from_splitted_quads = VDB::splitQuadsToTriangles(points, quads);

	std::vector<std::vector<float> > triangles_combined;
	triangles_combined = VDB::concatenateTriangleVectors(triangles, triangles_from_splitted_quads);

	std::cout << " active voxel count calculation_result_grid" << " = " << calculation_result_grid->activeVoxelCount() << std::endl;
	std::cout << "points size vdb isosurface" << " = " << points.size() << std::endl;
	std::cout << "triangles size vdb isosurface" << " = " << triangles.size() << std::endl;
	
	std::cout << " " << std::endl;
	std::cout << "points size isosurface" << " = " << points.size() << std::endl;
	std::cout << "triangles size triangulated isosurface" << " = " << triangles_combined.size() << std::endl;

	// export meshes
	VDB::exportTriangleMeshAsObj(points, triangles_combined);
	VDB::exportVDBMeshAsObj(points, triangles, quads);


	return 0;
}

