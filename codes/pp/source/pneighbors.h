




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip-mod/source/pneighbors.h"
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

#ifndef PNEIGHBORS_H
#define PNEIGHBORS_H

#include "manta.h"
#include "particle.h"

namespace Manta {

class ParticleNeighbors : public PbClass {public:
	typedef std::vector< std::pair<int,Real> > Neighbors;
	typedef std::vector<Neighbors> NeighborData;

	ParticleNeighbors(FluidSolver *parent) :PbClass(parent){} static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "ParticleNeighbors::ParticleNeighbors" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock);  obj = new ParticleNeighbors(parent); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleNeighbors::ParticleNeighbors" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("ParticleNeighbors::ParticleNeighbors",e.what()); return -1; } }

	
void update(const BasicParticleSystem &pts, const Grid<int> &index, const ParticleIndexSystem &indexSys, const Real radius, const ParticleDataImpl<int> *pT=NULL, const int exclude=0); static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleNeighbors* pbo = dynamic_cast<ParticleNeighbors*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "ParticleNeighbors::update" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const BasicParticleSystem& pts = *_args.getPtr<BasicParticleSystem >("pts",0,&_lock); const Grid<int> & index = *_args.getPtr<Grid<int>  >("index",1,&_lock); const ParticleIndexSystem& indexSys = *_args.getPtr<ParticleIndexSystem >("indexSys",2,&_lock); const Real radius = _args.get<Real >("radius",3,&_lock); const ParticleDataImpl<int> * pT = _args.getPtrOpt<ParticleDataImpl<int>  >("pT",4,NULL,&_lock); const int exclude = _args.getOpt<int >("exclude",5,0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->update(pts,index,indexSys,radius,pT,exclude);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleNeighbors::update" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("ParticleNeighbors::update",e.what()); return 0; } }
	void updateDistance(const BasicParticleSystem &pts); static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); ParticleNeighbors* pbo = dynamic_cast<ParticleNeighbors*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "ParticleNeighbors::updateDistance" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const BasicParticleSystem& pts = *_args.getPtr<BasicParticleSystem >("pts",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->updateDistance(pts);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"ParticleNeighbors::updateDistance" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("ParticleNeighbors::updateDistance",e.what()); return 0; } }

	int neighborIdx(const Neighbors::const_iterator &it) const { return it->first; }
	Real length(const Neighbors::const_iterator &it) const { return it->second; }
	Neighbors::const_iterator begin(const int idx) const { return _neighbor[idx].begin(); }
	Neighbors::const_iterator end(const int idx) const { return _neighbor[idx].end(); }
	int size(const int idx) const { return _neighbor[idx].size(); }

private: 	NeighborData _neighbor; public: PbArgs _args; }
#define _C_ParticleNeighbors
;

} // namespace

#endif	/* PNEIGHBORS_H */


