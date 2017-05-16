#ifndef ISOTESTSUITE_H
#define ISOTESTSUITE_H

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>


#include "atomprobe/primitives/ionHit.h"
#include "atomprobe/io/dataFiles.h"
#include "atomprobe/io/ranges.h"
#include <string>
#include <vector>
#include <typeinfo>

 
class TestIso : public CppUnit::TestFixture {


public:
 
	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Test iso");
 
		suiteOfTests->addTest(new CppUnit::TestCaller<TestIso>("Test - Test the Test itsself",
				&TestIso::testIso_TestTheTest ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestIso>("Test - Load Posfile",
				&TestIso::testIso_LoadPosFile));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestIso>("Test - Load Rangefile",
				&TestIso::testIso_LoadRangeFile));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestIso>("Test - Access Ion Positions",
				&TestIso::testIso_accessIonPositionsFromPos));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestIso>("Test - Access Ion IDs",
				&TestIso::testIso_accessIonIDs));

		return suiteOfTests;
	}
 
	/// Setup method
	void setUp() {}
 
	/// Teardown method
	void tearDown() {}
 
protected:
	void testIso_TestTheTest() 
	{
		int x = 1;
		int z = 2;
		int u = x + z;
		CPPUNIT_ASSERT_EQUAL(3, u);
	}	

	void testIso_LoadPosFile() 
	{
		std::string filename = "../../data/distributed_ref_pos_atomic_density_25_precRad_8_precConc_10_excessL_25_excessR_25.pos";
		const char * char_filename = filename.c_str();
		std::vector<AtomProbe::IonHit> posIons;
		loadPosFile(posIons, char_filename);

		double number_of_ions = posIons.size();
		double assert_ions = 1350750;
		CPPUNIT_ASSERT_EQUAL(assert_ions, number_of_ions);
	}

	void testIso_LoadRangeFile()
	{
		std::string filename = "../../data/sampled_ranges.rng";
		const char * char_filename = filename.c_str();

		AtomProbe::RangeFile rng;
		bool rng_check = rng.open(char_filename);

		CPPUNIT_ASSERT(rng_check);

		unsigned int number_of_ions {0};
		number_of_ions = rng.getNumIons();

		unsigned int noi_assert = 2;
		CPPUNIT_ASSERT_EQUAL(noi_assert,number_of_ions);

	}

	void testIso_accessIonPositionsFromPos()
	{
		std::string filename = "../../data/distributed_ref_pos_atomic_density_25_precRad_8_precConc_10_excessL_25_excessR_25.pos";
		const char * char_filename = filename.c_str();
		std::vector<AtomProbe::IonHit> posIons;
		loadPosFile(posIons, char_filename);

		AtomProbe::Point3D position;
		position = posIons[0].getPos();
		//std::cout << position[2] << std::endl;

		CPPUNIT_ASSERT_DOUBLES_EQUAL(1.48528,position[0],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(-7.99949407577515,position[1],0.01);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(-8.91878890991211,position[2],0.01);
	}

	void testIso_accessIonIDs()
	{
		std::string pos_filename = "../../data/distributed_ref_pos_atomic_density_25_precRad_8_precConc_10_excessL_25_excessR_25.pos";
		const char * pos_char_filename = pos_filename.c_str();
		std::vector<AtomProbe::IonHit> posIons;
		loadPosFile(posIons, pos_char_filename);

		std::string rng_filename = "../../data/sampled_ranges.rng";
		const char * rng_char_filename = rng_filename.c_str();

		AtomProbe::RangeFile rng;
		bool rng_check = rng.open(rng_char_filename);

		// short example
		AtomProbe::IonHit hit01;
		hit01 = posIons[0];
		unsigned int ionID = rng.getIonID(hit01);
		
		unsigned int zero_counter = 0;
		unsigned int one_counter = 1;
		unsigned int samplesize = 20000; // insert only up to ~ 4 mio cause of the data type -> overflow
		for (unsigned int i=0;i<samplesize;++i)
		{	
			AtomProbe::IonHit curr_hit;
			curr_hit = posIons[i];
		
			unsigned int curr_ionID = rng.getIonID(curr_hit);
			if (0 != curr_ionID)
			{
				one_counter += 1;
			}
			else 
			{
				zero_counter += 1;
			}
		}	
		// for the first 20 k ions
		unsigned int assert_zeros = 19733;
		unsigned int assert_ones = 268;

		CPPUNIT_ASSERT_EQUAL(assert_zeros,zero_counter);
		CPPUNIT_ASSERT_EQUAL(assert_ones,one_counter);



	}



};











#endif
