#ifndef CTF_H
#define CTF_H

#include <vector>
#include <math.h>       /* modf */

typedef struct coordinate {
double x;
double y;
double z;
} coordinate;

namespace CTF{

std::vector<std::vector<double> > initializeCubeVertices(double xmin=0, double ymin=0, double zmin=0);

std::vector<double> projectAtompositionToUnitvoxel(std::vector<double> atom_position, double voxel_size);

bool checkVertexCornerCoincidence(std::vector<double> atom_position);

std::vector<double> handleVertexCornerCoincidence(std::vector<double> atom_position);

std::vector<double> calcSubvolumes(std::vector<double> atom_position);

std::vector<double> calcVoxelContributions(std::vector<double> volumes_of_subcuboids);

std::vector<double> HellmanContributions(std::vector<double> volumes_of_subcuboids);

std::vector<std::vector<double> > determineAdjacentVoxelVertices(std::vector<double> atom_position, double voxel_size);

}

#endif
