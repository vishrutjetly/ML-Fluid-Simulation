




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip/source/general.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2016 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Globally used macros and functions (e.g. time measurements),
 * and doxygen documentation collection.
 *
 ******************************************************************************/


/*! \mainpage Welcome to mantaflow!
 *
 *  Here you can find the auto-generated documentation of the mantaflow framework.
 *
 *  One of the most useful parts is probably the list of python functions, classes and the C++ kernels.
 *  Those can be found found in the ''Modules'' section. For python functions the parameters (e.g. Grids, Real or int values)
 *  are automatically transferred to and from python. Thus, this list is a good reference how to call the functions used in the
 *  example scenes.
 * 
 */

// Define plugin documentation group
// all kernels, plugin functions and classes will automatically be added to this group
//! @defgroup Plugins Functions callable from Python
//! @defgroup PyClasses Classes exposed to Python
//! @defgroup Kernels Computation Kernels

#include "general.h"

#include <iomanip>

#if defined(WIN32) || defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#  undef WIN32_LEAN_AND_MEAN
#  undef NOMINMAX
#else
#   include <sys/time.h>
#	include "gitinfo.h"
#endif

using namespace std;

namespace Manta {

int gDebugLevel = 1;
 
void MuTime::get() {    
#if defined(WIN32) || defined(_WIN32)
	LARGE_INTEGER liTimerFrequency;
	QueryPerformanceFrequency(&liTimerFrequency);
	LARGE_INTEGER liLastTime;
	QueryPerformanceCounter(&liLastTime);
	time = static_cast<unsigned long>(liLastTime.QuadPart*(1e6/liTimerFrequency.QuadPart));
#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	time = tv.tv_sec*1e6 + tv.tv_usec;
#endif    
}

MuTime MuTime::update() {
	MuTime o = *this;
	get();
	return *this - o;
}

string MuTime::toString() const {
	stringstream ss;
	ss << *this;
	return ss.str();
}

std::ostream& operator<<(std::ostream& os, const MuTime& t) {
	const unsigned long ms = static_cast<unsigned long>(static_cast<double>(t.time)/6e7); // minutes
	const double ss = static_cast<double>(t.time - ms*6e7)/1e6; // seconds

	if(ms>0) os << ms << "m";
	os << ss << "s";
	return os;
}

//! print info about this mantaflow build, used eg by printBuildInfo in fluidsolver.cpp
std::string buildInfoString() {
	std::ostringstream infoStr;
	infoStr << "mantaflow " << MANTAVERSION;

	// os
#if defined(WIN32) || defined(_WIN32)
		infoStr << " win";
#	endif
#	ifdef __APPLE__
		infoStr << " mac";
#	endif
#	ifdef LINUX
		infoStr << " linux";
#	endif

	// 32/64 bit
	if (sizeof(size_t) == 8) 
		infoStr << " 64bit";
	else
		infoStr << " 32bit";

	// fp precision
#	if FLOATINGPOINT_PRECISION==2
		infoStr << " fp2";
#	else
		infoStr << " fp1";
#	endif

	// other compile switches
#	ifdef DEBUG
		infoStr << " debug";
#	endif 
#	ifdef OPENMP
		infoStr << " omp";
#	endif

	// repository info (git commit id)
#	ifndef MANTA_GIT_VERSION
#	define MANTA_GIT_VERSION "<unknown-commit>"
#	endif
	infoStr << " "<< MANTA_GIT_VERSION;

	infoStr << " from "<< __DATE__<<", "<<__TIME__;
	return infoStr.str();
}

//! note - generic PYTHON helpers in fluidsolver.cpp , no python bindings here

} // namespace


