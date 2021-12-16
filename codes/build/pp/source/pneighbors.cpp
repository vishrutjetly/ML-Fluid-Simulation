




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip/source/pneighbors.cpp"
// ----------------------------------------------------------------------------
//
// MantaFlow fluid solver framework
// Copyright 2015 Kiwon Um, Nils Thuerey
//
// This program is free software, distributed under the terms of the
// GNU General Public License (GPL)
// http://www.gnu.org/licenses
//
// Particle-neighbors data
//
// ----------------------------------------------------------------------------

#include "pneighbors.h"

namespace Manta {





 struct knUpdateNeighbor : public KernelBase { knUpdateNeighbor( ParticleNeighbors::NeighborData &neighbor, const BasicParticleSystem &pts, const Grid<int> &index, const ParticleIndexSystem &indexSys, const Real radius, const ParticleDataImpl<int> *pT, const int exclude) :  KernelBase(neighbor.size()) ,neighbor(neighbor),pts(pts),index(index),indexSys(indexSys),radius(radius),pT(pT),exclude(exclude)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleNeighbors::NeighborData &neighbor, const BasicParticleSystem &pts, const Grid<int> &index, const ParticleIndexSystem &indexSys, const Real radius, const ParticleDataImpl<int> *pT, const int exclude )  {
	if(!pts.isActive(idx) || (pT && (*pT)[idx]&exclude)) return;
	const Vec3 &xi = pts[idx].pos;
	const Vec3i xidx = toVec3i(xi);
	if(!index.isInBounds(xidx)) return;

	// search cells around each particle, idx
	neighbor[idx].clear();
	const Real radiusSqr = square(radius);
	const int r  = static_cast<int>(radius) + 1;
	const int rZ = pts.getParent()->is3D() ? r : 0;
	for(int k=xidx.z-rZ; k<=xidx.z+rZ; ++k) {
		for(int j=xidx.y-r; j<=xidx.y+r; ++j) {
			for(int i=xidx.x-r; i<=xidx.x+r; ++i) {
				if(!index.isInBounds(Vec3i(i, j, k))) continue;

				// loop for particles in each cell
				const int isysIdxS = index.index(i, j, k);
				const int pStart = index(isysIdxS);
				const int pEnd = index.isInBounds(isysIdxS+1) ? index(isysIdxS+1) : indexSys.size();
				for(int p=pStart; p<pEnd; ++p) {
					const int psrc = indexSys[p].sourceIndex;
					const Vec3 &xj = pts[psrc].pos;
					const Real lensqr_xij = normSquare(xi - xj);
					if(lensqr_xij>radiusSqr) continue;
					neighbor[idx].push_back(std::make_pair(psrc, std::sqrt(lensqr_xij)));
				}
			}
		}
	}
}    inline ParticleNeighbors::NeighborData& getArg0() { return neighbor; } typedef ParticleNeighbors::NeighborData type0;inline const BasicParticleSystem& getArg1() { return pts; } typedef BasicParticleSystem type1;inline const Grid<int> & getArg2() { return index; } typedef Grid<int>  type2;inline const ParticleIndexSystem& getArg3() { return indexSys; } typedef ParticleIndexSystem type3;inline const Real& getArg4() { return radius; } typedef Real type4;inline const ParticleDataImpl<int> * getArg5() { return pT; } typedef ParticleDataImpl<int>  type5;inline const int& getArg6() { return exclude; } typedef int type6; void runMessage() { debMsg("Executing kernel knUpdateNeighbor ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,neighbor,pts,index,indexSys,radius,pT,exclude);  }   } ParticleNeighbors::NeighborData& neighbor; const BasicParticleSystem& pts; const Grid<int> & index; const ParticleIndexSystem& indexSys; const Real radius; const ParticleDataImpl<int> * pT; const int exclude;   };
#line 22 "pneighbors.cpp"


void ParticleNeighbors::update(
	const BasicParticleSystem &pts, const Grid<int> &index, const ParticleIndexSystem &indexSys, const Real radius,
	const ParticleDataImpl<int> *pT, const int exclude) {
	_neighbor.resize(pts.size());
	knUpdateNeighbor(_neighbor, pts, index, indexSys, radius, pT, exclude);
}


 struct knUpdateDistance : public KernelBase { knUpdateDistance(ParticleNeighbors::NeighborData &neighbor, const BasicParticleSystem &pts) :  KernelBase(neighbor.size()) ,neighbor(neighbor),pts(pts)   { runMessage(); run(); }   inline void op(IndexInt idx, ParticleNeighbors::NeighborData &neighbor, const BasicParticleSystem &pts )  {
	for(unsigned int jj=0; jj<neighbor[idx].size(); ++jj) {
		const int j = neighbor[idx][jj].first;
		const Real len_xij = norm(pts[idx].pos - pts[j].pos);
		neighbor[idx][jj].second = len_xij;
	}
}    inline ParticleNeighbors::NeighborData& getArg0() { return neighbor; } typedef ParticleNeighbors::NeighborData type0;inline const BasicParticleSystem& getArg1() { return pts; } typedef BasicParticleSystem type1; void runMessage() { debMsg("Executing kernel knUpdateDistance ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,neighbor,pts);  }   } ParticleNeighbors::NeighborData& neighbor; const BasicParticleSystem& pts;   };
#line 61 "pneighbors.cpp"


void ParticleNeighbors::updateDistance(const BasicParticleSystem &pts) {
	knUpdateDistance(_neighbor, pts);
}

} // namespace


