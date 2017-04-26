#ifndef CTFTESTSUITE_H
#define CTFTESTSUITE_H

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>



#include "CTF_functions.h"
#include <vector>
#include <algorithm>



template <typename T>
bool compareVectors(std::vector<T> a, std::vector<T> b)
{

    const double EPSILON = 0.001;
    bool are_equal = true;
    for (int i = 0; i < a.size(); i++){

    	if ((fabs(a[i] - b[i])) > EPSILON){
        	are_equal = false;
                break;
               	}
    }
    return are_equal;
}

class TestCTF : public CppUnit::TestFixture {


public:
 
	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestCTF");
 
		suiteOfTests->addTest(new CppUnit::TestCaller<TestCTF>("Test0 - Test the Test itsself",
				&TestCTF::testCTF_TestTheTest ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestCTF>("Test1 - Calculate Subvolumes of a unit cubes eight subcuboids",
				&TestCTF::testCTF_CalculateSubvolumes));
	
		
		suiteOfTests->addTest(new CppUnit::TestCaller<TestCTF>("Test2 - Calculate Voxel Contributions",
				&TestCTF::testCTF_CalculateVoxelContributions));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestCTF>("Test3 - Project Atom Position to unit voxel",
				&TestCTF::testCTF_ProjectAtompositionToUnitvoxel));
/*
		suiteOfTests->addTest(new CppUnit::TestCaller<TestCTF>("Test4 - Determine Adjacent Voxel Vertices",
				&TestCTF::testCTF_DetermineAdjacentVoxelVertices));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestCTF>("Test6 - Hellman sawtooth contributions",
				&TestCTF::testCTF_HellmanSawtoothContributions));
*/


				
		return suiteOfTests;
	}
 
	/// Setup method
	void setUp() {}
 
	/// Teardown method
	void tearDown() {}
 
protected:
	void testCTF_TestTheTest() 
	{
		int x = 1;
		int z = 2;
		int u = x + z;
		CPPUNIT_ASSERT_EQUAL(3, u);
	}
	
	void testCTF_CalculateSubvolumes() 
	{	
		const unsigned int xyzs = 3;
		const unsigned int number_of_subvolumes = 8;
		// initial atom position
		std::vector<double> atom_position(xyzs);
		std::fill(atom_position.begin(), atom_position.end(),0.5);
		std::vector<double> volumes_of_subcuboids;
		volumes_of_subcuboids = CTF::calcSubvolumes(atom_position);
	
		std::vector<double> assertion_subvolumes(number_of_subvolumes);
		std::fill(assertion_subvolumes.begin(), assertion_subvolumes.end(),0.125);

		bool assertion_equal = false;
		assertion_equal = compareVectors(assertion_subvolumes, volumes_of_subcuboids);

		CPPUNIT_ASSERT(assertion_equal);
		
		// new atom position
		std::fill(atom_position.begin(), atom_position.end(),0.25);

		volumes_of_subcuboids = CTF::calcSubvolumes(atom_position);

		assertion_subvolumes = {0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625};
	
		assertion_equal = false;
		assertion_equal = compareVectors(assertion_subvolumes, volumes_of_subcuboids);

		CPPUNIT_ASSERT(assertion_equal);


	}
	
	void testCTF_CalculateVoxelContributions() 
	{	
		const unsigned int number_of_subvolumes = 8;
		std::vector<double> test_subvolumes(number_of_subvolumes); 
		std::fill(test_subvolumes.begin(), test_subvolumes.end(),0.125);
		
		std::vector<double> contributions_of_subcuboids;		
		contributions_of_subcuboids = CTF::calcVoxelContributions(test_subvolumes);
	
		std::vector<double> assertion_contributions(number_of_subvolumes);
		std::fill(assertion_contributions.begin(), assertion_contributions.end(),0.125);

		bool assertion_equal = false;
		assertion_equal = compareVectors(assertion_contributions, contributions_of_subcuboids);
		CPPUNIT_ASSERT(assertion_equal);
	
		test_subvolumes = {0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625};
		assertion_contributions = {0.421875,  0.140625,  0.046875,  0.140625,  0.140625,  0.046875, 0.015625,  0.046875};
		
		contributions_of_subcuboids = CTF::calcVoxelContributions(test_subvolumes);
		
		assertion_equal = false;
		assertion_equal = compareVectors(assertion_contributions, contributions_of_subcuboids);

/*
		for (auto &elem : contributions_of_subcuboids){

			std::cout << elem << std::endl;

		}
*/
		CPPUNIT_ASSERT(assertion_equal);

	}
	
	void testCTF_ProjectAtompositionToUnitvoxel()
	{
		const unsigned int xyzs = 3;
		std::vector<double> atom_position {0.5, 0.5, 0.5};
		std::vector<double> unit_position(xyzs);
		double voxel_size = 1;

		std::vector<double> position_in_unit_voxel(xyzs);
		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		std::vector<double> assert_position(xyzs);
		std::fill(assert_position.begin(), assert_position.end(),0.5);
		bool assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);

		atom_position = {2.5, 2.5, 2.5};
		voxel_size = 1;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		std::fill(assert_position.begin(), assert_position.end(),0.5);
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);

		atom_position = {-2.5, -2.5, -2.5};
		voxel_size = 1;

		std::fill(assert_position.begin(), assert_position.end(),0.5);
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);
		
		atom_position = {-2.5, -2.5, -2.5};
		voxel_size = 2;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		std::fill(assert_position.begin(), assert_position.end(),0.75);
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);
		

		atom_position = {-2.5, -2.5, -2.5};
		voxel_size = 3;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		std::fill(assert_position.begin(), assert_position.end(),0.166);
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);
		
		atom_position = {18.5, -2.3, 14.2};
		voxel_size = 1;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		assert_position = {0, 0, 0};
		assert_position[0] = 0.5;
		assert_position[1] = 0.7;
		assert_position[2] = 0.2;
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);

		// new case for a grid with much higher resolution with voxel sizes of let's say 0.1
		// how is this supposed to work?

		atom_position = {0.52, -1.31, -0.29};
		voxel_size = 0.1;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		assert_position = {0, 0, 0};
		assert_position[0] = 0.2;
		assert_position[1] = 0.9;
		assert_position[2] = 0.1;
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);

		// new case for a grid with much higher resolution with voxel sizes of let's say 0.1
		// positive
		atom_position = {0.05, 0.05, 0.05};
		voxel_size = 0.1;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		assert_position = {0, 0, 0};
		assert_position[0] = 0.5;
		assert_position[1] = 0.5;
		assert_position[2] = 0.5;
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);

		// new case for a grid with much higher resolution with voxel sizes of let's say 0.1

		atom_position = {0.09, 0.09, 0.09};
		unit_position = {0, 0, 0};
		voxel_size = 0.1;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		assert_position = {0, 0, 0};
		assert_position[0] = 0.9;
		assert_position[1] = 0.9;
		assert_position[2] = 0.9;
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);

		// new case for a grid with much higher resolution with voxel sizes of let's say 0.1

		atom_position = {0.9, 0.9, 0.9};
		unit_position = {0, 0, 0};
		voxel_size = 0.1;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		assert_position = {0, 0, 0};
		assert_position[0] = 1;
		assert_position[1] = 1;
		assert_position[2] = 1;
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);

		// new case for a grid with much higher resolution with voxel sizes of let's say 0.1

		atom_position = {0.25, 0.25, 0.25};
		unit_position = {0, 0, 0};
		voxel_size = 0.1;

		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);

		assert_position = {0, 0, 0};
		assert_position[0] = 0.5;
		assert_position[1] = 0.5;
		assert_position[2] = 0.5;
		assertion_equal = false;
		assertion_equal = compareVectors(assert_position, position_in_unit_voxel);
		CPPUNIT_ASSERT(assertion_equal);
		
	}
	
	void testCTF_DetermineAdjacentVoxelVertices()
	{
		// positive atom position
		std::vector<double>atom_position {2.5, 2.5, 2.5};
		double voxel_size = 1;
		const unsigned int number_of_vertices = 8;
		const unsigned int xyzs = 3;
		std::vector<std::vector<double> > assert_voxel_vertices = CTF::initializeCubeVertices(2,2,2);
		
		std::vector<std::vector<double> > surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}
		
		// bigger voxelsize
	
		atom_position = {2.5, 2.5, 2.5};
		voxel_size = 2;
		assert_voxel_vertices = CTF::initializeCubeVertices(1,1,1);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}

		// bigger voxelsize
	
		atom_position = {2.5, 2.5, 2.5};
		voxel_size = 3;
		assert_voxel_vertices = CTF::initializeCubeVertices(0,0,0);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}
		// negative atom position
	
		atom_position = {-2.5, -2.5, -2.5};
		voxel_size = 1;
		assert_voxel_vertices = CTF::initializeCubeVertices(-3,-3,-3);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}

		// negative atom position
	
		atom_position = {-1.7, -1.7, -1.7};
		voxel_size = 2;
		assert_voxel_vertices = CTF::initializeCubeVertices(-1,-1,-1);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}

		// negative atom position and double voxel size
	
		atom_position = {-1.7, -1.7, -1.7};
		voxel_size = 1.5;
		assert_voxel_vertices = CTF::initializeCubeVertices(-2,-2,-2);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}

		// negative atom position and double voxel size and arbitrary positions
	
		atom_position = {-1.9, -0.1, 0.7};
		voxel_size = 1.5;
		assert_voxel_vertices = CTF::initializeCubeVertices(-2,-1,0);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}
	
		// negative atom position and double voxel size and arbitrary positions
	
		atom_position = {-2.9, -0.1, 0.7};
		voxel_size = 2;
		
		assert_voxel_vertices = CTF::initializeCubeVertices(-2,-1,0);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}

		// new case for a grid with much higher resolution with voxel sizes of let's say 0.1
		// how is this supposed to work?

		atom_position = {0.5, -1.3, -0.2};
		voxel_size = 0.1;
		
		assert_voxel_vertices = CTF::initializeCubeVertices(5,-13,-2);
		
		surr_voxel_vertices = CTF::determineAdjacentVoxelVertices(atom_position, voxel_size);
		
		for (int i=0;i<assert_voxel_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_voxel_vertices[i][j], surr_voxel_vertices[i][j], 0.01);
			}
		}

	}
	
	void testCTF_HellmanSawtoothContributions()
	{

	// alternative ctf from http://nucapt.northwestern.edu/refbase/files/UM_95_199.pdf
	// the contribution is the volume of the opposing subcuboid divided by the total volume

		

		const unsigned int number_of_subcuboids = 8;
		const unsigned int xyzs = 3;
		double voxel_size =  1;
		std::vector<double> volumes_of_subcuboids(number_of_subcuboids);
		std::vector<double> vertex_contributions(number_of_subcuboids);
		std::vector<double> atom_position {0.5, 0.5, 0.5};
		std::vector<double> assertion_contributions(number_of_subcuboids);
		std::fill(assertion_contributions.begin(), assertion_contributions.end(),0.125);
	
		std::vector<double> position_in_unit_voxel(xyzs);
		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);
	
		bool vertex_corner_coincidence = CTF::checkVertexCornerCoincidence(position_in_unit_voxel);
	
		if (vertex_corner_coincidence == false)
		{
			volumes_of_subcuboids = CTF::calcSubvolumes(position_in_unit_voxel);
			vertex_contributions = CTF::HellmanContributions(volumes_of_subcuboids);
		}
		else
		{
			vertex_contributions = CTF::handleVertexCornerCoincidence(position_in_unit_voxel);
		}
	
		bool assertion_equal = false;
		assertion_equal = compareVectors(assertion_contributions, vertex_contributions);
		CPPUNIT_ASSERT(assertion_equal);
	
	
		voxel_size =  0.5;	
		atom_position = {0.25, 0.25, 0.25};
	
		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);
	
		vertex_corner_coincidence = CTF::checkVertexCornerCoincidence(position_in_unit_voxel);

		volumes_of_subcuboids = CTF::calcSubvolumes(position_in_unit_voxel);

		vertex_contributions = CTF::HellmanContributions(volumes_of_subcuboids);

		assertion_equal = false;
		assertion_equal = compareVectors(assertion_contributions, vertex_contributions);
		CPPUNIT_ASSERT(assertion_equal);


		voxel_size =  1;	
		atom_position = {0.99, 0.99, 0.99};

		assertion_contributions = {0, 0, 0, 0, 0, 0, 0.97, 0};
	
		position_in_unit_voxel = CTF::projectAtompositionToUnitvoxel(atom_position, voxel_size);
	
		vertex_corner_coincidence = CTF::checkVertexCornerCoincidence(position_in_unit_voxel);

		volumes_of_subcuboids = CTF::calcSubvolumes(position_in_unit_voxel);

		vertex_contributions = CTF::HellmanContributions(volumes_of_subcuboids);

		assertion_equal = false;
		assertion_equal = compareVectors(assertion_contributions, vertex_contributions);
		CPPUNIT_ASSERT(assertion_equal);

	}	


	
	

};

#endif