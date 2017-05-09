/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cc
 * Copyright (C) 2017 Lukas Gartmair <lukas@lgartmair-GA-H81M-D2V>
 * 
 * iso_testsuite is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * iso_testsuite is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
#include <iostream>

int main()
{
	std::cout << "Hello world!" << std::endl;
	return 0;
}
*/
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>
 
#include "iso_testsuite.h"
#include "CTF_testsuite.h"
#include "CTF_functions.h"
#include "vdb_testsuite.h"
#include "vdb_functions.h"

//sudo g++ TestRunner.cpp -lcppunit -o RunTests
//  sudo g++ TestRunner.cpp CTF_functions.cpp -lcppunit -std=c++11 -o RunTests


//sudo g++ TestRunner.cpp iso_functions.cpp /home/lukas/isosurface_executable/libatomprobe/src/io/dataFiles.cpp /home/lukas/isosurface_executable/libatomprobe/src/primitives/boundcube.cpp /home/lukas/isosurface_executable/libatomprobe/src/primitives/ionhit.cpp /home/lukas/isosurface_executable/libatomprobe/src/primitives/point3D.cpp -lcppunit -I/home/lukas/isosurface_executable/libatomprobe/include -I/home/lukas/isosurface_executable/libatomprobe/include/atomprobe/io -I/home/lukas/isosurface_executable/libatomprobe/src -L/home/lukas/isosurface_executable/libatomprobe/libatomprobe.so -std=c++11 -o RunTests


 
using namespace std;
 
int main() {
	CppUnit::TextUi::TestRunner runner;
 
	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestIso::suite());
	runner.addTest(TestCTF::suite());
	runner.addTest(TestOpenVDB::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();
 
	return 0;
}
