




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip/source/sphsolver.h"
// ----------------------------------------------------------------------------
//
// MantaFlow fluid solver framework
// Copyright 2015-2016 Kiwon Um, Nils Thuerey
//
// This program is free software, distributed under the terms of the
// GNU General Public License (GPL)
// http://www.gnu.org/licenses
//
// SPH simulation world
//
// ----------------------------------------------------------------------------

#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include "manta.h"

#include "particle.h"
#include "pneighbors.h"
#include "vectorbase.h"

namespace Manta {


class SphWorld : public PbClass {public:
	enum { FlagFluid=FlagGrid::TypeFluid, FlagObstacle=FlagGrid::TypeObstacle };

	static Real equationOfState(const Real d, const Real d0, const Real p0, const Real gamma, const Real chi=0.0) {
		return p0*(std::pow(d/d0, gamma) - 1.0) + chi;
	}

	


SphWorld( FluidSolver *parent, const Real alpha=0.08, const Real delta=0.5, const Real density=1e3, const Real eta=0.01, const Vec3 g=Vec3(0,-9.8,0), const Real gamma=7.0) :PbClass(parent),_alpha(alpha),_eta(eta),_g(g){
		updateDelta(delta);
		updateDensity(density);
		updateGamma(gamma);
		updateSoundSpeed(static_cast<Real>(parent->getGridSize().y)*std::fabs(_g.y)/_eta);
	} static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "SphWorld::SphWorld" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); const Real alpha = _args.getOpt<Real >("alpha",1,0.08,&_lock); const Real delta = _args.getOpt<Real >("delta",2,0.5,&_lock); const Real density = _args.getOpt<Real >("density",3,1e3,&_lock); const Real eta = _args.getOpt<Real >("eta",4,0.01,&_lock); const Vec3 g = _args.getOpt<Vec3 >("g",5,Vec3(0,-9.8,0),&_lock); const Real gamma = _args.getOpt<Real >("gamma",6,7.0,&_lock);  obj = new SphWorld(parent,alpha,delta,density,eta,g,gamma); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"SphWorld::SphWorld" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("SphWorld::SphWorld",e.what()); return -1; } }

	void bindParticleSystem(const BasicParticleSystem &p_system, const ParticleDataImpl<int> &p_type, const ParticleNeighbors &p_neighbor) {
		_pSystem = &p_system;
		_type = &p_type;
		_neighbor = &p_neighbor;
	} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::bindParticleSystem" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const BasicParticleSystem& p_system = *_args.getPtr<BasicParticleSystem >("p_system",0,&_lock); const ParticleDataImpl<int> & p_type = *_args.getPtr<ParticleDataImpl<int>  >("p_type",1,&_lock); const ParticleNeighbors& p_neighbor = *_args.getPtr<ParticleNeighbors >("p_neighbor",2,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->bindParticleSystem(p_system,p_type,p_neighbor);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::bindParticleSystem" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::bindParticleSystem",e.what()); return 0; } }

	void updateDelta(const Real d) { _delta = d; update_m0(); } static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::updateDelta" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real d = _args.get<Real >("d",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->updateDelta(d);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::updateDelta" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::updateDelta",e.what()); return 0; } }
	void updateDensity(const Real d) { _d0 = d; update_m0(); update_p0(); } static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::updateDensity" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real d = _args.get<Real >("d",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->updateDensity(d);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::updateDensity" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::updateDensity",e.what()); return 0; } }
	void updateSoundSpeed(const Real c) { _c = c; update_p0(); } static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::updateSoundSpeed" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real c = _args.get<Real >("c",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->updateSoundSpeed(c);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::updateSoundSpeed" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::updateSoundSpeed",e.what()); return 0; } }
	void updateGamma(const Real gamma) { _gamma = gamma; update_p0(); } static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::updateGamma" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real gamma = _args.get<Real >("gamma",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->updateGamma(gamma);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::updateGamma" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::updateGamma",e.what()); return 0; } }

	const BasicParticleSystem &particleSystem() const { return *_pSystem; }
	const ParticleDataImpl<int> &particleType() const { return *_type; }
	const ParticleNeighbors &neighborData() const { return *_neighbor; }

	Real limitDtByVmax(const Real dt, const Real h, const Real vmax, const Real a=0.4) const {
		return std::min(dt, a*h/(vmax+VECTOR_EPSILON));
	} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::limitDtByVmax" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real dt = _args.get<Real >("dt",0,&_lock); const Real h = _args.get<Real >("h",1,&_lock); const Real vmax = _args.get<Real >("vmax",2,&_lock); const Real a = _args.getOpt<Real >("a",3,0.4,&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->limitDtByVmax(dt,h,vmax,a));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::limitDtByVmax" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::limitDtByVmax",e.what()); return 0; } }
	Real limitDtByFmax(const Real dt, const Real h, const Real fmax, const Real a=0.25) const {
		return std::min(dt, a*std::sqrt(h/(fmax+VECTOR_EPSILON)));
	} static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::limitDtByFmax" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real dt = _args.get<Real >("dt",0,&_lock); const Real h = _args.get<Real >("h",1,&_lock); const Real fmax = _args.get<Real >("fmax",2,&_lock); const Real a = _args.getOpt<Real >("a",3,0.25,&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->limitDtByFmax(dt,h,fmax,a));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::limitDtByFmax" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::limitDtByFmax",e.what()); return 0; } }
	Real limitDtByGravity(const Real dt, const Real h, const Real a=0.25) const {
		const Real norm_g = norm(_g);
		return (norm_g<=VECTOR_EPSILON) ? dt : std::min(dt, a*std::sqrt(h/norm_g));
	} static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::limitDtByGravity" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real dt = _args.get<Real >("dt",0,&_lock); const Real h = _args.get<Real >("h",1,&_lock); const Real a = _args.getOpt<Real >("a",2,0.25,&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->limitDtByGravity(dt,h,a));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::limitDtByGravity" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::limitDtByGravity",e.what()); return 0; } }
	Real limitDtByViscous(const Real dt, const Real h, const Real a=0.4) const { // WCSPH [Monaghan 1992, Becker and Teschner 2007]
		return std::min(dt, a*h/(_c*(static_cast<Real>(1.0)+static_cast<Real>(0.6)*_alpha)));
	} static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::limitDtByViscous" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const Real dt = _args.get<Real >("dt",0,&_lock); const Real h = _args.get<Real >("h",1,&_lock); const Real a = _args.getOpt<Real >("a",2,0.4,&_lock);  pbo->_args.copy(_args);  _retval = toPy(pbo->limitDtByViscous(dt,h,a));  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::limitDtByViscous" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::limitDtByViscous",e.what()); return 0; } }

	void showParameters() const {
		std::cout << "SPH parameters:" << std::endl;
		std::cout << "\t" << "alpha = " << _alpha << std::endl;
		std::cout << "\t" << "c = " << _c << std::endl;
		std::cout << "\t" << "d0 = " << _d0 << std::endl;
		std::cout << "\t" << "delta = " << _delta << std::endl;
		std::cout << "\t" << "eta = " << _eta << std::endl;
		std::cout << "\t" << "g = " << _g << std::endl;
		std::cout << "\t" << "gamma = " << _gamma << std::endl;
		std::cout << "\t" << "m0 = " << _m0 << std::endl;
		std::cout << "\t" << "p0 = " << _p0 << std::endl;
	} static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "SphWorld::showParameters" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->showParameters();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"SphWorld::showParameters" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("SphWorld::showParameters",e.what()); return 0; } }

	Real _alpha;static PyObject* _GET__alpha(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_alpha); } static int _SET__alpha(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_alpha = fromPy<Real  >(val); return 0; }
	Real _c;static PyObject* _GET__c(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_c); } static int _SET__c(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_c = fromPy<Real  >(val); return 0; }
	Real _d0;static PyObject* _GET__d0(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_d0); } static int _SET__d0(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_d0 = fromPy<Real  >(val); return 0; }
	Real _delta;static PyObject* _GET__delta(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_delta); } static int _SET__delta(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_delta = fromPy<Real  >(val); return 0; } // particle spacing
	Real _eta;static PyObject* _GET__eta(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_eta); } static int _SET__eta(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_eta = fromPy<Real  >(val); return 0; }
	Vec3 _g;static PyObject* _GET__g(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_g); } static int _SET__g(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_g = fromPy<Vec3  >(val); return 0; }
	Real _gamma;static PyObject* _GET__gamma(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_gamma); } static int _SET__gamma(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_gamma = fromPy<Real  >(val); return 0; }
	Real _m0;static PyObject* _GET__m0(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_m0); } static int _SET__m0(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_m0 = fromPy<Real  >(val); return 0; }
	Real _p0;static PyObject* _GET__p0(PyObject* self, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); return toPy(pbo->_p0); } static int _SET__p0(PyObject* self, PyObject* val, void* cl) { SphWorld* pbo = dynamic_cast<SphWorld*>(Pb::objFromPy(self)); pbo->_p0 = fromPy<Real  >(val); return 0; }

private:
	void update_p0() { _p0 = _d0*_c*_c/_gamma; }
	void update_m0() { _m0 = _d0*std::pow(_delta, this->getParent()->is2D() ? 2.0 : 3.0); }

	// mandatory data
	const BasicParticleSystem *_pSystem;
	const ParticleDataImpl<int> *_type; 	const ParticleNeighbors *_neighbor; public: PbArgs _args; }
#define _C_SphWorld
;

// Functors
struct Pij {			// [Gingold and Monaghan, 1982, JCP]
	explicit Pij(const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const Real m, const ParticleDataImpl<int> &t) : _p(p), _d(d), _t(t), _mSqr(square(m)) {}
	Real operator()(const int i, int j) const { j=((_t[j]&SphWorld::FlagObstacle)?i:j); return _mSqr*(_p[i]/square(_d[i]) + _p[j]/square(_d[j])); }
	const ParticleDataImpl<Real> &_p, &_d;
	const ParticleDataImpl<int> &_t;
	const Real _mSqr;
};
struct PijMij {			// [Gingold and Monaghan, 1982, JCP] + [Akinci et al., 2012, SIG/TOG]
	explicit PijMij(const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &m, const ParticleDataImpl<int> &t) : _p(p), _d(d), _m(m), _t(t) {}
	Real operator()(const int i, int j) const { j=((_t[j]&SphWorld::FlagObstacle)?i:j); return _m[i]*_m[j]*(_p[i]/square(_d[i]) + _p[j]/square(_d[j])); }
	const ParticleDataImpl<Real> &_p, &_d, &_m;
	const ParticleDataImpl<int> &_t;
};
struct PijDfsph {		// [Bender and Koschier, 2015, SCA]
	explicit PijDfsph(const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const Real m, const ParticleDataImpl<int> &t) : _p(p), _d(d), _t(t), _mSqr(square(m)) {}
	Real operator()(const int i, int j) const { j=((_t[j]&SphWorld::FlagObstacle)?i:j); return _mSqr*(_p[i]/_d[i] + _p[j]/_d[j]); }
	const ParticleDataImpl<Real> &_p, &_d;
	const ParticleDataImpl<int> &_t;
	const Real _mSqr;
};
struct PijMijDfsph {		// [Bender and Koschier, 2015, SCA] + [Akinci et al., 2012, SIG/TOG]
	explicit PijMijDfsph(const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &m, const ParticleDataImpl<int> &t) : _p(p), _d(d), _m(m), _t(t) {}
	Real operator()(const int i, int j) const { j=((_t[j]&SphWorld::FlagObstacle)?i:j); return _m[i]*_m[j]*(_p[i]/_d[i] + _p[j]/_d[j]); }
	const ParticleDataImpl<Real> &_p, &_d, &_m;
	const ParticleDataImpl<int> &_t;
};
struct PijWeighted {		// [Hu and Adams, 2007, JCP]
	explicit PijWeighted(const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &Vsqr) : _p(p), _d(d), _Vsqr(Vsqr) {}
	Real operator()(const int i, const int j) const { return (_Vsqr[i]+_Vsqr[j])*(_p[i]*_d[j] + _p[j]*_d[i])/(_d[i]+_d[j]); }
	const ParticleDataImpl<Real> &_p, &_d, &_Vsqr;
};
struct PijBg {
	explicit PijBg(const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &Vsqr) : _p(p), _Vsqr(Vsqr) {}
	Real operator()(const int i, const int j) const { return 10e0*std::fabs(_p[i])*(_Vsqr[i]+_Vsqr[j]); }
	const ParticleDataImpl<Real> &_p, &_Vsqr;
};

} // namespace

#endif	/* SPHSOLVER_H */


