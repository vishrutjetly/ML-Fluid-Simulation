




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip-mod/source/sphkernel.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2015-2016 Kiwon Um, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL)
 * http://www.gnu.org/licenses
 *
 * SPH Kernels
 *
 ******************************************************************************/

#ifndef SPHKERNEL_H
#define SPHKERNEL_H

#include "manta.h"

namespace Manta {


class CubicSpline : public PbClass {public:
	CubicSpline(FluidSolver *parent, const Real h=1) :PbClass(parent){
		_dim = (parent->is2D()) ? 2 : 3;
		setRadius(h);
	} static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "CubicSpline::CubicSpline" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); const Real h = _args.getOpt<Real >("h",1,1,&_lock);  obj = new CubicSpline(parent,h); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"CubicSpline::CubicSpline" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("CubicSpline::CubicSpline",e.what()); return -1; } }
	void setRadius(const Real h) {
		const Real h2 = square(h), h3 = h2*h, h4 = h3*h, h5 = h4*h;
		_h = h;
		_sr = 2e0*h;
		_c[0]  = 2e0/(3e0*h);
		_c[1]  = 10e0/(7e0*M_PI*h2);
		_c[2]  = 1e0/(M_PI*h3);
		_gc[0] = 3e0/(2e0*h3);
		_gc[1] = 45e0/(14e0*M_PI*h4);
		_gc[2] = 9e0/(4e0*M_PI*h5);
	}
	Real radius() const { return _h; } static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); CubicSpline* pbo = dynamic_cast<CubicSpline*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "CubicSpline::radius" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->radius());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"CubicSpline::radius" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("CubicSpline::radius",e.what()); return 0; } }
	Real supportRadius() const { return _sr; } static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); CubicSpline* pbo = dynamic_cast<CubicSpline*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "CubicSpline::supportRadius" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->supportRadius());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"CubicSpline::supportRadius" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("CubicSpline::supportRadius",e.what()); return 0; } }

	Real f(const Real l) const {
		const Real q = l/_h;
		if(q<1e0) return _c[_dim-1]*(1e0 - 1.5*square(q) + 0.75*cubed(q));
		else if(q<2e0) return _c[_dim-1]*(0.25*cubed(2e0-q));
		return 0;
	}
	Real derivative_f(const Real l) const {
		const Real q = l/_h;
		if(q<=1e0) return _gc[_dim-1]*(q-4e0/3e0)*l;
		else if(q<2e0) return -_gc[_dim-1]*square(2e0-q)*_h/3e0;
		return 0;
	}

	Real w(const Vec3 &rij) const { return f(norm(rij)); }
	Vec3 grad_w(const Vec3 &rij) const { return grad_w(rij, norm(rij)); }
	Vec3 grad_w(const Vec3 &rij, const Real len) const { return derivative_f(len)*rij/len; }

private:
	unsigned int _dim; 	Real _h, _sr, _c[3], _gc[3]; public: PbArgs _args; }
#define _C_CubicSpline
;


class BndKernel : public PbClass {public:
	BndKernel(FluidSolver *parent, const Real c, const Real h=1) :PbClass(parent),_c(c){
		setRadius(h);
	} static int _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "BndKernel::BndKernel" , !noTiming ); { ArgLocker _lock; FluidSolver* parent = _args.getPtr<FluidSolver >("parent",0,&_lock); const Real c = _args.get<Real >("c",1,&_lock); const Real h = _args.getOpt<Real >("h",2,1,&_lock);  obj = new BndKernel(parent,c,h); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"BndKernel::BndKernel" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("BndKernel::BndKernel",e.what()); return -1; } }
	void setRadius(const Real h) {
		_h = h;
		_sr = 2e0*h;
	}
	Real radius() const { return _h; }
	Real supportRadius() const { return _sr; }

	Real f(const Real l) const {
		const Real q = l/_h;
		const Real s = 0.02*square(_c)/l;
		if(q<=2e0/3e0) return s*2e0/3e0;
		else if(q<1e0) return s*(2e0*q - square(q)*2e0/3e0);
		else if(q<2e0) return s*0.5*square(2e0 - q);
		return 0e0;
	}

private: 	Real _h, _sr, _c; public: PbArgs _args; }
#define _C_BndKernel
;

inline Real wSmooth(const Real &r2, const Real h2) { return (r2>h2) ? 0.0 : 1.0 - r2/h2; }

} // namespaces

#endif	/* SPHKERNEL_H */


