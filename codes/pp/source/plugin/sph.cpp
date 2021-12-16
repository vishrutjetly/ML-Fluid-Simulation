




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip-mod/source/plugin/sph.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2015-2016 Kiwon Um, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL)
 * http://www.gnu.org/licenses
 *
 * SPH (smoothed particle hydrodynamics) plugins
 *
 * TODO: Different kernels should be usable.
 *
 ******************************************************************************/

#include "matrixbase.h"
#include "particle.h"
#include "pneighbors.h"
#include "sphkernel.h"
#include "sphsolver.h"

namespace Manta {




 struct knSphUpdateVelocity : public KernelBase { knSphUpdateVelocity( ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &vn, const ParticleDataImpl<Vec3> &f, const Real dt, const SphWorld &sph, const int itype) :  KernelBase(v.size()) ,v(v),vn(vn),f(f),dt(dt),sph(sph),itype(itype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &vn, const ParticleDataImpl<Vec3> &f, const Real dt, const SphWorld &sph, const int itype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	v[idx] = vn[idx] + f[idx]*dt/sph._m0;
}    inline ParticleDataImpl<Vec3> & getArg0() { return v; } typedef ParticleDataImpl<Vec3>  type0;inline const ParticleDataImpl<Vec3> & getArg1() { return vn; } typedef ParticleDataImpl<Vec3>  type1;inline const ParticleDataImpl<Vec3> & getArg2() { return f; } typedef ParticleDataImpl<Vec3>  type2;inline const Real& getArg3() { return dt; } typedef Real type3;inline const SphWorld& getArg4() { return sph; } typedef SphWorld type4;inline const int& getArg5() { return itype; } typedef int type5; void runMessage() { debMsg("Executing kernel knSphUpdateVelocity ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,v,vn,f,dt,sph,itype);  }   } ParticleDataImpl<Vec3> & v; const ParticleDataImpl<Vec3> & vn; const ParticleDataImpl<Vec3> & f; const Real dt; const SphWorld& sph; const int itype;   };
#line 27 "plugin/sph.cpp"






void sphUpdateVelocity( ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &vn, const ParticleDataImpl<Vec3> &f, const Real dt, const SphWorld &sph, const int itype=SphWorld::FlagFluid) {
	knSphUpdateVelocity(v, vn, f, dt, sph, itype);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphUpdateVelocity" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",0,&_lock); const ParticleDataImpl<Vec3> & vn = *_args.getPtr<ParticleDataImpl<Vec3>  >("vn",1,&_lock); const ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",2,&_lock); const Real dt = _args.get<Real >("dt",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagFluid,&_lock);   _retval = getPyNone(); sphUpdateVelocity(v,vn,f,dt,sph,itype);  _args.check(); } pbFinalizePlugin(parent,"sphUpdateVelocity", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphUpdateVelocity",e.what()); return 0; } } static const Pb::Register _RP_sphUpdateVelocity ("","sphUpdateVelocity",_W_0);  extern "C" { void PbRegister_sphUpdateVelocity() { KEEP_UNUSED(_RP_sphUpdateVelocity); } } 



 struct knSphSetBndVelocity : public KernelBase { knSphSetBndVelocity( ParticleDataImpl<Vec3> &v, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(v.size()) ,v(v),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &v, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	v[idx] = Vec3(0.0);

	Real sum_w = 0.0;
	Vec3 sum_vw(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it); if(!(t[j]&jtype)) continue;
		const Real w = k.f(n.length(it));
		sum_w += w;
		sum_vw += v[j]*w;
	}
	if(sum_w>Real(0)) v[idx] = -sum_vw/sum_w;
}    inline ParticleDataImpl<Vec3> & getArg0() { return v; } typedef ParticleDataImpl<Vec3>  type0;inline const CubicSpline& getArg1() { return k; } typedef CubicSpline type1;inline const SphWorld& getArg2() { return sph; } typedef SphWorld type2;inline const int& getArg3() { return itype; } typedef int type3;inline const int& getArg4() { return jtype; } typedef int type4; void runMessage() { debMsg("Executing kernel knSphSetBndVelocity ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,v,k,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & v; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 44 "plugin/sph.cpp"







void sphSetBndVelocity( ParticleDataImpl<Vec3> &v, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagObstacle, const int jtype=SphWorld::FlagFluid) {
	knSphSetBndVelocity(v, k, sph, itype, jtype);
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphSetBndVelocity" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",0,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); const int itype = _args.getOpt<int >("itype",3,SphWorld::FlagObstacle,&_lock); const int jtype = _args.getOpt<int >("jtype",4,SphWorld::FlagFluid,&_lock);   _retval = getPyNone(); sphSetBndVelocity(v,k,sph,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphSetBndVelocity", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphSetBndVelocity",e.what()); return 0; } } static const Pb::Register _RP_sphSetBndVelocity ("","sphSetBndVelocity",_W_1);  extern "C" { void PbRegister_sphSetBndVelocity() { KEEP_UNUSED(_RP_sphSetBndVelocity); } } 



 struct knSphUpdatePosition : public KernelBase { knSphUpdatePosition( BasicParticleSystem &x, const ParticleDataImpl<Vec3> &v, const Real dt, const SphWorld &sph, const int itype) :  KernelBase(x.size()) ,x(x),v(v),dt(dt),sph(sph),itype(itype)   { runMessage(); run(); }   inline void op(IndexInt idx,  BasicParticleSystem &x, const ParticleDataImpl<Vec3> &v, const Real dt, const SphWorld &sph, const int itype )  {
	const ParticleDataImpl<int> &t = sph.particleType();

	if(!x.isActive(idx) || !(t[idx]&itype)) return;
	x[idx].pos += dt*v[idx];
}    inline BasicParticleSystem& getArg0() { return x; } typedef BasicParticleSystem type0;inline const ParticleDataImpl<Vec3> & getArg1() { return v; } typedef ParticleDataImpl<Vec3>  type1;inline const Real& getArg2() { return dt; } typedef Real type2;inline const SphWorld& getArg3() { return sph; } typedef SphWorld type3;inline const int& getArg4() { return itype; } typedef int type4; void runMessage() { debMsg("Executing kernel knSphUpdatePosition ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,x,v,dt,sph,itype);  }   } BasicParticleSystem& x; const ParticleDataImpl<Vec3> & v; const Real dt; const SphWorld& sph; const int itype;   };
#line 73 "plugin/sph.cpp"






void sphUpdatePosition( BasicParticleSystem &x, const ParticleDataImpl<Vec3> &v, const Real dt, const SphWorld &sph, const int itype=SphWorld::FlagFluid) {
	knSphUpdatePosition(x, v, dt, sph, itype);
} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphUpdatePosition" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& x = *_args.getPtr<BasicParticleSystem >("x",0,&_lock); const ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",1,&_lock); const Real dt = _args.get<Real >("dt",2,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",3,&_lock); const int itype = _args.getOpt<int >("itype",4,SphWorld::FlagFluid,&_lock);   _retval = getPyNone(); sphUpdatePosition(x,v,dt,sph,itype);  _args.check(); } pbFinalizePlugin(parent,"sphUpdatePosition", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphUpdatePosition",e.what()); return 0; } } static const Pb::Register _RP_sphUpdatePosition ("","sphUpdatePosition",_W_2);  extern "C" { void PbRegister_sphUpdatePosition() { KEEP_UNUSED(_RP_sphUpdatePosition); } } 



 struct knSphComputeDensity : public KernelBase { knSphComputeDensity( ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(d.size()) ,d(d),m(m),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	Real sum_w = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it); if(!(t[j]&jtype)) continue;
		const Real len_xij = n.length(it);
		sum_w += ((m) ? m->get(j) : sph._m0) * k.f(len_xij);
	}
	d[idx] = sum_w;
}    inline ParticleDataImpl<Real> & getArg0() { return d; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> * getArg1() { return m; } typedef ParticleDataImpl<Real>  type1;inline const CubicSpline& getArg2() { return k; } typedef CubicSpline type2;inline const SphWorld& getArg3() { return sph; } typedef SphWorld type3;inline const int& getArg4() { return itype; } typedef int type4;inline const int& getArg5() { return jtype; } typedef int type5; void runMessage() { debMsg("Executing kernel knSphComputeDensity ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,d,m,k,sph,itype,jtype);  }   } ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 89 "plugin/sph.cpp"







void sphComputeDensity( ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeDensity(d, m, k, sph, itype, jtype);
} static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeDensity" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",0,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",3,NULL,&_lock); const int itype = _args.getOpt<int >("itype",4,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",5,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeDensity(d,k,sph,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeDensity", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeDensity",e.what()); return 0; } } static const Pb::Register _RP_sphComputeDensity ("","sphComputeDensity",_W_3);  extern "C" { void PbRegister_sphComputeDensity() { KEEP_UNUSED(_RP_sphComputeDensity); } } 


 struct knSphComputeVolume : public KernelBase { knSphComputeVolume( ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(d.size()) ,d(d),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	Real sum_w = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it); if(!(t[j]&jtype)) continue;
		const Real len_xij = n.length(it);
		sum_w += k.f(len_xij);
	}
	d[idx] = 1.0/sum_w;
}    inline ParticleDataImpl<Real> & getArg0() { return d; } typedef ParticleDataImpl<Real>  type0;inline const CubicSpline& getArg1() { return k; } typedef CubicSpline type1;inline const SphWorld& getArg2() { return sph; } typedef SphWorld type2;inline const int& getArg3() { return itype; } typedef int type3;inline const int& getArg4() { return jtype; } typedef int type4; void runMessage() { debMsg("Executing kernel knSphComputeVolume ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,d,k,sph,itype,jtype);  }   } ParticleDataImpl<Real> & d; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 113 "plugin/sph.cpp"






void sphComputeVolume( ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid) {
	knSphComputeVolume(d, k, sph, itype, jtype);
} static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeVolume" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",0,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); const int itype = _args.getOpt<int >("itype",3,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",4,SphWorld::FlagFluid,&_lock);   _retval = getPyNone(); sphComputeVolume(d,k,sph,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeVolume", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeVolume",e.what()); return 0; } } static const Pb::Register _RP_sphComputeVolume ("","sphComputeVolume",_W_4);  extern "C" { void PbRegister_sphComputeVolume() { KEEP_UNUSED(_RP_sphComputeVolume); } } 



 struct knSphComputeDivergence : public KernelBase { knSphComputeDivergence( ParticleDataImpl<Real> &div, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(div.size()) ,div(div),d(d),v(v),m(m),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &div, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos, &vi = v[idx];
	Real sum_div = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos, &vj = v[j];
		Vec3 xij = xi-xj, vij = vi-vj; if(pts.getParent()->is2D()) xij.z = vij.z = 0.0;
		const Real Vj = ((m) ? m->get(j) : sph._m0)/d[j];
		sum_div += dot(vij, k.grad_w(xij, len_xij))*Vj;
	}
	div[idx] = -sum_div;
}    inline ParticleDataImpl<Real> & getArg0() { return div; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d; } typedef ParticleDataImpl<Real>  type1;inline const ParticleDataImpl<Vec3> & getArg2() { return v; } typedef ParticleDataImpl<Vec3>  type2;inline const ParticleDataImpl<Real> * getArg3() { return m; } typedef ParticleDataImpl<Real>  type3;inline const CubicSpline& getArg4() { return k; } typedef CubicSpline type4;inline const SphWorld& getArg5() { return sph; } typedef SphWorld type5;inline const int& getArg6() { return itype; } typedef int type6;inline const int& getArg7() { return jtype; } typedef int type7; void runMessage() { debMsg("Executing kernel knSphComputeDivergence ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,div,d,v,m,k,sph,itype,jtype);  }   } ParticleDataImpl<Real> & div; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Vec3> & v; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 137 "plugin/sph.cpp"







void sphComputeDivergence( ParticleDataImpl<Real> &div, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &v, const CubicSpline &k, const SphWorld &sph, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeDivergence(div, d, v, m, k, sph, itype, jtype);
} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeDivergence" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & div = *_args.getPtr<ParticleDataImpl<Real>  >("div",0,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",1,&_lock); const ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",5,NULL,&_lock); const int itype = _args.getOpt<int >("itype",6,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",7,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeDivergence(div,d,v,k,sph,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeDivergence", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeDivergence",e.what()); return 0; } } static const Pb::Register _RP_sphComputeDivergence ("","sphComputeDivergence",_W_5);  extern "C" { void PbRegister_sphComputeDivergence() { KEEP_UNUSED(_RP_sphComputeDivergence); } } 



 struct knSphComputeDivergenceSimple : public KernelBase { knSphComputeDivergenceSimple( ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(d.size()) ,d(d),v(v),m(m),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos, &vi = v[idx];
	Real sum_div = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos, &vj = v[j];
		Vec3 xij = xi - xj, vij = vi - vj; if(pts.getParent()->is2D()) xij.z = vij.z = 0.0;
		sum_div += dot(vij, k.grad_w(xij, len_xij))*((m) ? m->get(j) : sph._m0);
	}
	d[idx] = -sum_div;
}    inline ParticleDataImpl<Real> & getArg0() { return d; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Vec3> & getArg1() { return v; } typedef ParticleDataImpl<Vec3>  type1;inline const ParticleDataImpl<Real> * getArg2() { return m; } typedef ParticleDataImpl<Real>  type2;inline const CubicSpline& getArg3() { return k; } typedef CubicSpline type3;inline const SphWorld& getArg4() { return sph; } typedef SphWorld type4;inline const int& getArg5() { return itype; } typedef int type5;inline const int& getArg6() { return jtype; } typedef int type6; void runMessage() { debMsg("Executing kernel knSphComputeDivergenceSimple ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,d,v,m,k,sph,itype,jtype);  }   } ParticleDataImpl<Real> & d; const ParticleDataImpl<Vec3> & v; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 167 "plugin/sph.cpp"








void sphComputeDivergenceSimple( ParticleDataImpl<Real> &div, const ParticleDataImpl<Vec3> &v, const CubicSpline &k, const SphWorld &sph, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeDivergenceSimple(div, v, m, k, sph, itype, jtype);
} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeDivergenceSimple" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & div = *_args.getPtr<ParticleDataImpl<Real>  >("div",0,&_lock); const ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",1,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",2,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",3,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",4,NULL,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",6,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeDivergenceSimple(div,v,k,sph,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeDivergenceSimple", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeDivergenceSimple",e.what()); return 0; } } static const Pb::Register _RP_sphComputeDivergenceSimple ("","sphComputeDivergenceSimple",_W_6);  extern "C" { void PbRegister_sphComputeDivergenceSimple() { KEEP_UNUSED(_RP_sphComputeDivergenceSimple); } } 



 struct knSphReinitDensity : public KernelBase { knSphReinitDensity( ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &dn, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(d.size()) ,d(d),dn(dn),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &dn, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	Real sum_w = 0.0, sum_w_over_d = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it); if(!(t[j]&jtype)) continue;
		const Real len_xij = n.length(it);
		const Real w = k.f(len_xij);
		sum_w += w;
		sum_w_over_d += w/dn[j];
	}
	d[idx] = sum_w/sum_w_over_d; // NOTE: always sum_w_over_d > 0
}    inline ParticleDataImpl<Real> & getArg0() { return d; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return dn; } typedef ParticleDataImpl<Real>  type1;inline const CubicSpline& getArg2() { return k; } typedef CubicSpline type2;inline const SphWorld& getArg3() { return sph; } typedef SphWorld type3;inline const int& getArg4() { return itype; } typedef int type4;inline const int& getArg5() { return jtype; } typedef int type5; void runMessage() { debMsg("Executing kernel knSphReinitDensity ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,d,dn,k,sph,itype,jtype);  }   } ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> & dn; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 197 "plugin/sph.cpp"







void sphReinitDensity( ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &dn, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphReinitDensity(d, dn, k, sph, itype, jtype);
} static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphReinitDensity" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",0,&_lock); const ParticleDataImpl<Real> & dn = *_args.getPtr<ParticleDataImpl<Real>  >("dn",1,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",2,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",3,&_lock); const int itype = _args.getOpt<int >("itype",4,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",5,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphReinitDensity(d,dn,k,sph,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphReinitDensity", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphReinitDensity",e.what()); return 0; } } static const Pb::Register _RP_sphReinitDensity ("","sphReinitDensity",_W_7);  extern "C" { void PbRegister_sphReinitDensity() { KEEP_UNUSED(_RP_sphReinitDensity); } } 



 struct knSphComputePressure : public KernelBase { knSphComputePressure( ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const SphWorld &sph, const int itype) :  KernelBase(p.size()) ,p(p),d(d),sph(sph),itype(itype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const SphWorld &sph, const int itype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	p[idx] = sph.equationOfState(d[idx], sph._d0, sph._p0, sph._gamma);
}    inline ParticleDataImpl<Real> & getArg0() { return p; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d; } typedef ParticleDataImpl<Real>  type1;inline const SphWorld& getArg2() { return sph; } typedef SphWorld type2;inline const int& getArg3() { return itype; } typedef int type3; void runMessage() { debMsg("Executing kernel knSphComputePressure ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,p,d,sph,itype);  }   } ParticleDataImpl<Real> & p; const ParticleDataImpl<Real> & d; const SphWorld& sph; const int itype;   };
#line 225 "plugin/sph.cpp"






void sphComputePressure( ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const SphWorld &sph, const int itype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputePressure(p, d, sph, itype);
} static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputePressure" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",0,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); const int itype = _args.getOpt<int >("itype",3,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputePressure(p,d,sph,itype);  _args.check(); } pbFinalizePlugin(parent,"sphComputePressure", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputePressure",e.what()); return 0; } } static const Pb::Register _RP_sphComputePressure ("","sphComputePressure",_W_8);  extern "C" { void PbRegister_sphComputePressure() { KEEP_UNUSED(_RP_sphComputePressure); } } 



 struct knSphComputeCurl : public KernelBase { knSphComputeCurl( ParticleDataImpl<Vec3> &w, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(w.size()) ,w(w),v(v),d(d),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &w, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos, &vi = v[idx];
	Vec3 sum_v(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);
		const Real len_xij = n.length(it);
		if(!(t[j]&jtype) || idx==j || len_xij<=Real(0)) continue;

		const Vec3 &xj = pts[j].pos, &vj = v[j];
		Vec3 xij = xi - xj;
		Vec3 vij = vi - vj; if(pts.getParent()->is2D()) xij.z = vij.z = 0.0;
		const Real Vj = sph._m0/d[j];
		sum_v += cross(vij, k.grad_w(xij, len_xij))*Vj;
	}
	w[idx] = sum_v;
}    inline ParticleDataImpl<Vec3> & getArg0() { return w; } typedef ParticleDataImpl<Vec3>  type0;inline const ParticleDataImpl<Vec3> & getArg1() { return v; } typedef ParticleDataImpl<Vec3>  type1;inline const ParticleDataImpl<Real> & getArg2() { return d; } typedef ParticleDataImpl<Real>  type2;inline const CubicSpline& getArg3() { return k; } typedef CubicSpline type3;inline const SphWorld& getArg4() { return sph; } typedef SphWorld type4;inline const int& getArg5() { return itype; } typedef int type5;inline const int& getArg6() { return jtype; } typedef int type6; void runMessage() { debMsg("Executing kernel knSphComputeCurl ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,w,v,d,k,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & w; const ParticleDataImpl<Vec3> & v; const ParticleDataImpl<Real> & d; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 242 "plugin/sph.cpp"







void sphComputeCurl( ParticleDataImpl<Vec3> &w, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeCurl(w, v, d, k, sph, itype, jtype);
} static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeCurl" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & w = *_args.getPtr<ParticleDataImpl<Vec3>  >("w",0,&_lock); const ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",1,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",6,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeCurl(w,v,d,k,sph,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeCurl", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeCurl",e.what()); return 0; } } static const Pb::Register _RP_sphComputeCurl ("","sphComputeCurl",_W_9);  extern "C" { void PbRegister_sphComputeCurl() { KEEP_UNUSED(_RP_sphComputeCurl); } } 




 struct knSphAddJacobian : public KernelBase { knSphAddJacobian( ParticleDataImpl<Vec3> &diag, ParticleDataImpl<Vec3> &offu, ParticleDataImpl<Vec3> &offl, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(diag.size()) ,diag(diag),offu(offu),offl(offl),v(v),d(d),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &diag, ParticleDataImpl<Vec3> &offu, ParticleDataImpl<Vec3> &offl, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos, &vi = v[idx];
	Vec3 sum_diag(0.0), sum_offu(0.0), sum_offl(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);
		const Real len_xij = n.length(it);
		if(!(t[j]&jtype) || idx==j || len_xij<=Real(0)) continue;

		const Vec3 &xj = pts[j].pos, &vj = v[j];
		Vec3 xij = xi - xj;
		Vec3 vji = vj - vi; if(pts.getParent()->is2D()) xij.z = vji.z = 0.0;
		const Vec3 grad_w = k.grad_w(xij, len_xij);
		const Real Vj = sph._m0/d[j];
		sum_diag += vji*grad_w*Vj;
		sum_offu.x += vji.x*grad_w.y*Vj; // xy
		sum_offu.y += vji.y*grad_w.z*Vj; // yz
		sum_offu.z += vji.z*grad_w.x*Vj; // zx
		sum_offl.x += vji.y*grad_w.x*Vj; // yx
		sum_offl.y += vji.z*grad_w.y*Vj; // zy
		sum_offl.z += vji.x*grad_w.z*Vj; // xz
	}
	diag[idx] += sum_diag;
	offu[idx] += sum_offu;
	offl[idx] += sum_offl;
}    inline ParticleDataImpl<Vec3> & getArg0() { return diag; } typedef ParticleDataImpl<Vec3>  type0;inline ParticleDataImpl<Vec3> & getArg1() { return offu; } typedef ParticleDataImpl<Vec3>  type1;inline ParticleDataImpl<Vec3> & getArg2() { return offl; } typedef ParticleDataImpl<Vec3>  type2;inline const ParticleDataImpl<Vec3> & getArg3() { return v; } typedef ParticleDataImpl<Vec3>  type3;inline const ParticleDataImpl<Real> & getArg4() { return d; } typedef ParticleDataImpl<Real>  type4;inline const CubicSpline& getArg5() { return k; } typedef CubicSpline type5;inline const SphWorld& getArg6() { return sph; } typedef SphWorld type6;inline const int& getArg7() { return itype; } typedef int type7;inline const int& getArg8() { return jtype; } typedef int type8; void runMessage() { debMsg("Executing kernel knSphAddJacobian ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,diag,offu,offl,v,d,k,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & diag; ParticleDataImpl<Vec3> & offu; ParticleDataImpl<Vec3> & offl; const ParticleDataImpl<Vec3> & v; const ParticleDataImpl<Real> & d; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 276 "plugin/sph.cpp"


//	11 12 13
// S =	21 22 23
//	31 32 33
// diag=(11 22 33), offu = (12, 23, 13), offl = (21, 32, 31)






void sphAddJacobian( ParticleDataImpl<Vec3> &diag, ParticleDataImpl<Vec3> &offu, ParticleDataImpl<Vec3> &offl, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphAddJacobian(diag, offu, offl, v, d, k, sph, itype, jtype);
} static PyObject* _W_10 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphAddJacobian" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & diag = *_args.getPtr<ParticleDataImpl<Vec3>  >("diag",0,&_lock); ParticleDataImpl<Vec3> & offu = *_args.getPtr<ParticleDataImpl<Vec3>  >("offu",1,&_lock); ParticleDataImpl<Vec3> & offl = *_args.getPtr<ParticleDataImpl<Vec3>  >("offl",2,&_lock); const ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",3,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",4,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",5,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",6,&_lock); const int itype = _args.getOpt<int >("itype",7,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",8,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphAddJacobian(diag,offu,offl,v,d,k,sph,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphAddJacobian", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphAddJacobian",e.what()); return 0; } } static const Pb::Register _RP_sphAddJacobian ("","sphAddJacobian",_W_10);  extern "C" { void PbRegister_sphAddJacobian() { KEEP_UNUSED(_RP_sphAddJacobian); } } 



 struct knSphSetBndValues : public KernelBase { knSphSetBndValues( ParticleDataImpl<Real> &d, ParticleDataImpl<Real> &p, ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype, const bool clamp) :  KernelBase(d.size()) ,d(d),p(p),Vsqr(Vsqr),k(k),sph(sph),itype(itype),jtype(jtype),clamp(clamp)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &d, ParticleDataImpl<Real> &p, ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype, const bool clamp )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	const Vec3 &xi = pts[idx].pos;

	d[idx]	  = sph._d0;
	p[idx]	  = 0.0;
	Vsqr[idx] = square(sph._m0/sph._d0);

	Real sum_w = 0.0, sum_pw=0.0;
	Vec3 sum_dx_wf(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);
		if(!(t[j]&jtype)) continue;
		const Real w = k.f(n.length(it));
		const Vec3 &xj = pts[j].pos;
		sum_w += w;
		sum_pw += p[j]*w;
		sum_dx_wf += d[j]*(xi -xj)*w;
	}
	if(sum_w>Real(0)) {
		p[idx]	  = (sum_pw + dot(sph._g, sum_dx_wf))/sum_w;
		if(clamp) p[idx] = std::max(Real(0), p[idx]);
		d[idx]	  = sph._d0*std::pow(p[idx]/sph._p0 + Real(1), Real(1)/sph._gamma);
		Vsqr[idx] = square(sph._m0/d[idx]);
	}
}    inline ParticleDataImpl<Real> & getArg0() { return d; } typedef ParticleDataImpl<Real>  type0;inline ParticleDataImpl<Real> & getArg1() { return p; } typedef ParticleDataImpl<Real>  type1;inline ParticleDataImpl<Real> & getArg2() { return Vsqr; } typedef ParticleDataImpl<Real>  type2;inline const CubicSpline& getArg3() { return k; } typedef CubicSpline type3;inline const SphWorld& getArg4() { return sph; } typedef SphWorld type4;inline const int& getArg5() { return itype; } typedef int type5;inline const int& getArg6() { return jtype; } typedef int type6;inline const bool& getArg7() { return clamp; } typedef bool type7; void runMessage() { debMsg("Executing kernel knSphSetBndValues ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,d,p,Vsqr,k,sph,itype,jtype,clamp);  }   } ParticleDataImpl<Real> & d; ParticleDataImpl<Real> & p; ParticleDataImpl<Real> & Vsqr; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype; const bool clamp;   };
#line 323 "plugin/sph.cpp"








void sphSetBndValues( ParticleDataImpl<Real> &d, ParticleDataImpl<Real> &p, ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagObstacle, const int jtype=SphWorld::FlagFluid, const bool clamp=false) {
	knSphSetBndValues(d, p, Vsqr, k, sph, itype, jtype, clamp);
} static PyObject* _W_11 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphSetBndValues" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",0,&_lock); ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",1,&_lock); ParticleDataImpl<Real> & Vsqr = *_args.getPtr<ParticleDataImpl<Real>  >("Vsqr",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagObstacle,&_lock); const int jtype = _args.getOpt<int >("jtype",6,SphWorld::FlagFluid,&_lock); const bool clamp = _args.getOpt<bool >("clamp",7,false,&_lock);   _retval = getPyNone(); sphSetBndValues(d,p,Vsqr,k,sph,itype,jtype,clamp);  _args.check(); } pbFinalizePlugin(parent,"sphSetBndValues", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphSetBndValues",e.what()); return 0; } } static const Pb::Register _RP_sphSetBndValues ("","sphSetBndValues",_W_11);  extern "C" { void PbRegister_sphSetBndValues() { KEEP_UNUSED(_RP_sphSetBndValues); } } 



 struct knSphComputeIisphDii : public KernelBase { knSphComputeIisphDii( ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype) :  KernelBase(dii.size()) ,dii(dii),d(d),m(m),k(k),sph(sph),dt(dt),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	const Real overdi2 = safeDivide(Real(1.0), square(d[idx]));

	dii[idx].x = dii[idx].y = dii[idx].z = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		dii[idx] -= k.grad_w(xi - xj, len_xij)*((m) ? m->get(j) : sph._m0)*overdi2;
	}
}    inline ParticleDataImpl<Vec3> & getArg0() { return dii; } typedef ParticleDataImpl<Vec3>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d; } typedef ParticleDataImpl<Real>  type1;inline const ParticleDataImpl<Real> * getArg2() { return m; } typedef ParticleDataImpl<Real>  type2;inline const CubicSpline& getArg3() { return k; } typedef CubicSpline type3;inline const SphWorld& getArg4() { return sph; } typedef SphWorld type4;inline const Real& getArg5() { return dt; } typedef Real type5;inline const int& getArg6() { return itype; } typedef int type6;inline const int& getArg7() { return jtype; } typedef int type7; void runMessage() { debMsg("Executing kernel knSphComputeIisphDii ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,dii,d,m,k,sph,dt,itype,jtype);  }   } ParticleDataImpl<Vec3> & dii; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const Real dt; const int itype; const int jtype;   };
#line 365 "plugin/sph.cpp"







void sphComputeIisphDii( ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const Real dt, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeIisphDii(dii, d, m, k, sph, dt, itype, jtype);
} static PyObject* _W_12 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeIisphDii" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & dii = *_args.getPtr<ParticleDataImpl<Vec3>  >("dii",0,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",1,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",2,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",3,&_lock); const Real dt = _args.get<Real >("dt",4,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",5,NULL,&_lock); const int itype = _args.getOpt<int >("itype",6,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",7,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeIisphDii(dii,d,k,sph,dt,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeIisphDii", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeIisphDii",e.what()); return 0; } } static const Pb::Register _RP_sphComputeIisphDii ("","sphComputeIisphDii",_W_12);  extern "C" { void PbRegister_sphComputeIisphDii() { KEEP_UNUSED(_RP_sphComputeIisphDii); } } 




 struct knSphComputeIisphAii : public KernelBase { knSphComputeIisphAii( ParticleDataImpl<Real> &aii, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype) :  KernelBase(aii.size()) ,aii(aii),d(d),dii(dii),m(m),k(k),sph(sph),dt(dt),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &aii, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	const Real mi_di2 = safeDivide((m) ? m->get(idx) : sph._m0, square(d[idx]));

	aii[idx] = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		const Vec3 grad_w = k.grad_w(xi - xj, len_xij);
		const Vec3 dji = (t[j]&SphWorld::FlagFluid) ? grad_w*mi_di2 : Vec3(0.0);
		aii[idx] += dot(dii[idx] - dji, grad_w)*((m) ? m->get(j) : sph._m0);
	}
}    inline ParticleDataImpl<Real> & getArg0() { return aii; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d; } typedef ParticleDataImpl<Real>  type1;inline const ParticleDataImpl<Vec3> & getArg2() { return dii; } typedef ParticleDataImpl<Vec3>  type2;inline const ParticleDataImpl<Real> * getArg3() { return m; } typedef ParticleDataImpl<Real>  type3;inline const CubicSpline& getArg4() { return k; } typedef CubicSpline type4;inline const SphWorld& getArg5() { return sph; } typedef SphWorld type5;inline const Real& getArg6() { return dt; } typedef Real type6;inline const int& getArg7() { return itype; } typedef int type7;inline const int& getArg8() { return jtype; } typedef int type8; void runMessage() { debMsg("Executing kernel knSphComputeIisphAii ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,aii,d,dii,m,k,sph,dt,itype,jtype);  }   } ParticleDataImpl<Real> & aii; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Vec3> & dii; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const Real dt; const int itype; const int jtype;   };
#line 395 "plugin/sph.cpp"







void sphComputeIisphAii( ParticleDataImpl<Real> &aii, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Vec3> &dii, const CubicSpline &k, const SphWorld &sph, const Real dt, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeIisphAii(aii, d, dii, m, k, sph, dt, itype, jtype);
} static PyObject* _W_13 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeIisphAii" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & aii = *_args.getPtr<ParticleDataImpl<Real>  >("aii",0,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",1,&_lock); const ParticleDataImpl<Vec3> & dii = *_args.getPtr<ParticleDataImpl<Vec3>  >("dii",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const Real dt = _args.get<Real >("dt",5,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",6,NULL,&_lock); const int itype = _args.getOpt<int >("itype",7,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",8,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeIisphAii(aii,d,dii,k,sph,dt,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeIisphAii", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeIisphAii",e.what()); return 0; } } static const Pb::Register _RP_sphComputeIisphAii ("","sphComputeIisphAii",_W_13);  extern "C" { void PbRegister_sphComputeIisphAii() { KEEP_UNUSED(_RP_sphComputeIisphAii); } } 




 struct knSphComputeIisphDijPj : public KernelBase { knSphComputeIisphDijPj( ParticleDataImpl<Vec3> &dijpj, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype) :  KernelBase(dijpj.size()) ,dijpj(dijpj),d(d),p(p),m(m),k(k),sph(sph),dt(dt),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &dijpj, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;

	dijpj[idx].x = dijpj[idx].y = dijpj[idx].z = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		dijpj[idx] -= safeDivide(k.grad_w(xi - xj, len_xij)*((m) ? m->get(j) : sph._m0)*p[j], Vec3(square(d[j])));
	}
}    inline ParticleDataImpl<Vec3> & getArg0() { return dijpj; } typedef ParticleDataImpl<Vec3>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d; } typedef ParticleDataImpl<Real>  type1;inline const ParticleDataImpl<Real> & getArg2() { return p; } typedef ParticleDataImpl<Real>  type2;inline const ParticleDataImpl<Real> * getArg3() { return m; } typedef ParticleDataImpl<Real>  type3;inline const CubicSpline& getArg4() { return k; } typedef CubicSpline type4;inline const SphWorld& getArg5() { return sph; } typedef SphWorld type5;inline const Real& getArg6() { return dt; } typedef Real type6;inline const int& getArg7() { return itype; } typedef int type7;inline const int& getArg8() { return jtype; } typedef int type8; void runMessage() { debMsg("Executing kernel knSphComputeIisphDijPj ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,dijpj,d,p,m,k,sph,dt,itype,jtype);  }   } ParticleDataImpl<Vec3> & dijpj; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> & p; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const Real dt; const int itype; const int jtype;   };
#line 427 "plugin/sph.cpp"







void sphComputeIisphDijPj( ParticleDataImpl<Vec3> &dijpj, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &p, const CubicSpline &k, const SphWorld &sph, const Real dt, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeIisphDijPj(dijpj, d, p, m, k, sph, dt, itype, jtype);
} static PyObject* _W_14 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeIisphDijPj" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & dijpj = *_args.getPtr<ParticleDataImpl<Vec3>  >("dijpj",0,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",1,&_lock); const ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const Real dt = _args.get<Real >("dt",5,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",6,NULL,&_lock); const int itype = _args.getOpt<int >("itype",7,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",8,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeIisphDijPj(dijpj,d,p,k,sph,dt,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeIisphDijPj", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeIisphDijPj",e.what()); return 0; } } static const Pb::Register _RP_sphComputeIisphDijPj ("","sphComputeIisphDijPj",_W_14);  extern "C" { void PbRegister_sphComputeIisphDijPj() { KEEP_UNUSED(_RP_sphComputeIisphDijPj); } } 






 struct knSphComputeIisphP : public KernelBase { knSphComputeIisphP( ParticleDataImpl<Real> &p_next, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d_adv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &aii, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Vec3> &dijpj, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype) :  KernelBase(p_next.size()) ,p_next(p_next),p(p),d_adv(d_adv),d(d),aii(aii),dii(dii),dijpj(dijpj),m(m),k(k),sph(sph),dt(dt),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &p_next, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d_adv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &aii, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Vec3> &dijpj, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	if(aii[idx]==Real(0)) return; // isolated particle

	const Vec3 &xi = pts[idx].pos;
	const Real overdt2 = safeDivide(Real(1.0), square(dt));
	const Real mi_di2 = safeDivide((m) ? m->get(idx) : sph._m0, square(d[idx]));

	Real sumv=0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		const Vec3 grad_w = k.grad_w(xi-xj, len_xij);
		const Vec3 dji = grad_w*mi_di2;
		const Real mj = ((m) ? m->get(j) : sph._m0);
		const Vec3 dFi = dijpj[idx];
		const Vec3 dFj = (t[j]&SphWorld::FlagFluid) ? dii[j]*p[j] + dijpj[j] - dji*p[idx] : Vec3(0.0);
		sumv += dot(dFi - dFj, grad_w)*mj;
	}
	p_next[idx] = 0.5*p[idx] + 0.5*safeDivide((sph._d0 - d_adv[idx])*overdt2 - sumv, aii[idx]);
}    inline ParticleDataImpl<Real> & getArg0() { return p_next; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return p; } typedef ParticleDataImpl<Real>  type1;inline const ParticleDataImpl<Real> & getArg2() { return d_adv; } typedef ParticleDataImpl<Real>  type2;inline const ParticleDataImpl<Real> & getArg3() { return d; } typedef ParticleDataImpl<Real>  type3;inline const ParticleDataImpl<Real> & getArg4() { return aii; } typedef ParticleDataImpl<Real>  type4;inline const ParticleDataImpl<Vec3> & getArg5() { return dii; } typedef ParticleDataImpl<Vec3>  type5;inline const ParticleDataImpl<Vec3> & getArg6() { return dijpj; } typedef ParticleDataImpl<Vec3>  type6;inline const ParticleDataImpl<Real> * getArg7() { return m; } typedef ParticleDataImpl<Real>  type7;inline const CubicSpline& getArg8() { return k; } typedef CubicSpline type8;inline const SphWorld& getArg9() { return sph; } typedef SphWorld type9;inline const Real& getArg10() { return dt; } typedef Real type10;inline const int& getArg11() { return itype; } typedef int type11;inline const int& getArg12() { return jtype; } typedef int type12; void runMessage() { debMsg("Executing kernel knSphComputeIisphP ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,p_next,p,d_adv,d,aii,dii,dijpj,m,k,sph,dt,itype,jtype);  }   } ParticleDataImpl<Real> & p_next; const ParticleDataImpl<Real> & p; const ParticleDataImpl<Real> & d_adv; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> & aii; const ParticleDataImpl<Vec3> & dii; const ParticleDataImpl<Vec3> & dijpj; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const Real dt; const int itype; const int jtype;   };
#line 458 "plugin/sph.cpp"









void sphComputeIisphP( ParticleDataImpl<Real> &p_next, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d_adv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &aii, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Vec3> &dijpj, const CubicSpline &k, const SphWorld &sph, const Real dt, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeIisphP(p_next, p, d_adv, d, aii, dii, dijpj, m, k, sph, dt, itype, jtype);
} static PyObject* _W_15 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeIisphP" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & p_next = *_args.getPtr<ParticleDataImpl<Real>  >("p_next",0,&_lock); const ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",1,&_lock); const ParticleDataImpl<Real> & d_adv = *_args.getPtr<ParticleDataImpl<Real>  >("d_adv",2,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",3,&_lock); const ParticleDataImpl<Real> & aii = *_args.getPtr<ParticleDataImpl<Real>  >("aii",4,&_lock); const ParticleDataImpl<Vec3> & dii = *_args.getPtr<ParticleDataImpl<Vec3>  >("dii",5,&_lock); const ParticleDataImpl<Vec3> & dijpj = *_args.getPtr<ParticleDataImpl<Vec3>  >("dijpj",6,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",7,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",8,&_lock); const Real dt = _args.get<Real >("dt",9,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",10,NULL,&_lock); const int itype = _args.getOpt<int >("itype",11,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",12,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeIisphP(p_next,p,d_adv,d,aii,dii,dijpj,k,sph,dt,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeIisphP", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeIisphP",e.what()); return 0; } } static const Pb::Register _RP_sphComputeIisphP ("","sphComputeIisphP",_W_15);  extern "C" { void PbRegister_sphComputeIisphP() { KEEP_UNUSED(_RP_sphComputeIisphP); } } 





 struct knSphComputeIisphD : public KernelBase { knSphComputeIisphD( ParticleDataImpl<Real> &d_next, const ParticleDataImpl<Real> &d_adv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Vec3> &dijpj, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype) :  KernelBase(d_next.size()) ,d_next(d_next),d_adv(d_adv),d(d),p(p),dii(dii),dijpj(dijpj),m(m),k(k),sph(sph),dt(dt),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &d_next, const ParticleDataImpl<Real> &d_adv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Vec3> &dijpj, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real dt, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;

	Real sumv=0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		const Vec3 grad_w = k.grad_w(xi-xj, len_xij);
		const Real mj = ((m) ? m->get(j) : sph._m0);
		const Vec3 dFi = (t[j]&SphWorld::FlagFluid) ? dii[idx]*p[idx] + dijpj[idx] : -safeDivide(mj, square(d[idx]))*p[idx]*grad_w;
		const Vec3 dFj = (t[j]&SphWorld::FlagFluid) ? dii[j]*p[j] + dijpj[j] : Vec3(0.0);
		sumv += dot(dFi - dFj, grad_w)*mj;
	}
	d_next[idx] = d_adv[idx] + square(dt)*sumv;
}    inline ParticleDataImpl<Real> & getArg0() { return d_next; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d_adv; } typedef ParticleDataImpl<Real>  type1;inline const ParticleDataImpl<Real> & getArg2() { return d; } typedef ParticleDataImpl<Real>  type2;inline const ParticleDataImpl<Real> & getArg3() { return p; } typedef ParticleDataImpl<Real>  type3;inline const ParticleDataImpl<Vec3> & getArg4() { return dii; } typedef ParticleDataImpl<Vec3>  type4;inline const ParticleDataImpl<Vec3> & getArg5() { return dijpj; } typedef ParticleDataImpl<Vec3>  type5;inline const ParticleDataImpl<Real> * getArg6() { return m; } typedef ParticleDataImpl<Real>  type6;inline const CubicSpline& getArg7() { return k; } typedef CubicSpline type7;inline const SphWorld& getArg8() { return sph; } typedef SphWorld type8;inline const Real& getArg9() { return dt; } typedef Real type9;inline const int& getArg10() { return itype; } typedef int type10;inline const int& getArg11() { return jtype; } typedef int type11; void runMessage() { debMsg("Executing kernel knSphComputeIisphD ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,d_next,d_adv,d,p,dii,dijpj,m,k,sph,dt,itype,jtype);  }   } ParticleDataImpl<Real> & d_next; const ParticleDataImpl<Real> & d_adv; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> & p; const ParticleDataImpl<Vec3> & dii; const ParticleDataImpl<Vec3> & dijpj; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const Real dt; const int itype; const int jtype;   };
#line 499 "plugin/sph.cpp"








void sphComputeIisphD( ParticleDataImpl<Real> &d_next, const ParticleDataImpl<Real> &d_adv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Vec3> &dii, const ParticleDataImpl<Vec3> &dijpj, const CubicSpline &k, const SphWorld &sph, const Real dt, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeIisphD(d_next, d_adv, d, p, dii, dijpj, m, k, sph, dt, itype, jtype);
} static PyObject* _W_16 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeIisphD" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & d_next = *_args.getPtr<ParticleDataImpl<Real>  >("d_next",0,&_lock); const ParticleDataImpl<Real> & d_adv = *_args.getPtr<ParticleDataImpl<Real>  >("d_adv",1,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",2,&_lock); const ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",3,&_lock); const ParticleDataImpl<Vec3> & dii = *_args.getPtr<ParticleDataImpl<Vec3>  >("dii",4,&_lock); const ParticleDataImpl<Vec3> & dijpj = *_args.getPtr<ParticleDataImpl<Vec3>  >("dijpj",5,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",6,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",7,&_lock); const Real dt = _args.get<Real >("dt",8,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",9,NULL,&_lock); const int itype = _args.getOpt<int >("itype",10,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",11,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeIisphD(d_next,d_adv,d,p,dii,dijpj,k,sph,dt,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeIisphD", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeIisphD",e.what()); return 0; } } static const Pb::Register _RP_sphComputeIisphD ("","sphComputeIisphD",_W_16);  extern "C" { void PbRegister_sphComputeIisphD() { KEEP_UNUSED(_RP_sphComputeIisphD); } } 

Real sphComputePcisphDelta(const CubicSpline &k, const SphWorld &sph) {
	const int r  = static_cast<int>(k.supportRadius()/sph._delta) + 1;
	int	  k0 = -r, k1 = r;
	if(sph.particleSystem().getParent()->is2D()) { k0 = 0; k1 = 1; }
	const Real Rsqr = square(k.supportRadius());
	const Real over_beta = 0.5*square(sph._d0/sph._m0);

	Vec3 sum_gw(0.0);
	Real sum_gwgw=0.0;
	for(int di=-r; di<r; ++di) {
		for(int dj=-r; dj<r; ++dj) {
			for(int dk=k0; dk<k1; ++dk) {
				const Vec3 xj(sph._delta*static_cast<Real>(di), sph._delta*static_cast<Real>(dj), sph._delta*static_cast<Real>(dk));
				const Vec3 xij = -xj;
				const Real lensqr_xij = normSquare(xij);
				if(lensqr_xij<=Real(0) || Rsqr<lensqr_xij) continue;

				const Real len_xij = std::sqrt(lensqr_xij);
				const Vec3 grad_w = k.derivative_f(len_xij)*xij/len_xij;
				sum_gw += grad_w;
				sum_gwgw += dot(grad_w, grad_w);
			}
		}
	}

	return over_beta/(dot(sum_gw, sum_gw) + sum_gwgw);
} static PyObject* _W_17 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputePcisphDelta" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; const CubicSpline& k = *_args.getPtr<CubicSpline >("k",0,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",1,&_lock);   _retval = toPy(sphComputePcisphDelta(k,sph));  _args.check(); } pbFinalizePlugin(parent,"sphComputePcisphDelta", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputePcisphDelta",e.what()); return 0; } } static const Pb::Register _RP_sphComputePcisphDelta ("","sphComputePcisphDelta",_W_17);  extern "C" { void PbRegister_sphComputePcisphDelta() { KEEP_UNUSED(_RP_sphComputePcisphDelta); } } 




 struct knSphComputeDfsphAlpha : public KernelBase { knSphComputeDfsphAlpha( ParticleDataImpl<Real> &a, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(a.size()) ,a(a),d(d),m(m),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &a, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	Vec3 sum_mgw(0.0);
	Real sum_mgw_sqr = 0.0;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		const Vec3 gw_ij = k.grad_w(xi - xj, len_xij);
		const Real mj = ((m) ? m->get(j) : sph._m0);
		sum_mgw += mj*gw_ij;
		sum_mgw_sqr += (t[j]&SphWorld::FlagFluid) ? square(mj)*normSquare(gw_ij) : 0.0;
	}
	a[idx] = d[idx]/std::max(VECTOR_EPSILON, normSquare(sum_mgw) + sum_mgw_sqr);
}    inline ParticleDataImpl<Real> & getArg0() { return a; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d; } typedef ParticleDataImpl<Real>  type1;inline const ParticleDataImpl<Real> * getArg2() { return m; } typedef ParticleDataImpl<Real>  type2;inline const CubicSpline& getArg3() { return k; } typedef CubicSpline type3;inline const SphWorld& getArg4() { return sph; } typedef SphWorld type4;inline const int& getArg5() { return itype; } typedef int type5;inline const int& getArg6() { return jtype; } typedef int type6; void runMessage() { debMsg("Executing kernel knSphComputeDfsphAlpha ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,a,d,m,k,sph,itype,jtype);  }   } ParticleDataImpl<Real> & a; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 562 "plugin/sph.cpp"







void sphComputeDfsphAlpha( ParticleDataImpl<Real> &a, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphComputeDfsphAlpha(a, d, m, k, sph, itype, jtype);
} static PyObject* _W_18 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeDfsphAlpha" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & a = *_args.getPtr<ParticleDataImpl<Real>  >("a",0,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",1,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",2,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",3,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",4,NULL,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",6,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphComputeDfsphAlpha(a,d,k,sph,m,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphComputeDfsphAlpha", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeDfsphAlpha",e.what()); return 0; } } static const Pb::Register _RP_sphComputeDfsphAlpha ("","sphComputeDfsphAlpha",_W_18);  extern "C" { void PbRegister_sphComputeDfsphAlpha() { KEEP_UNUSED(_RP_sphComputeDfsphAlpha); } } 


 struct knSphStompIsolatedParticleValue : public KernelBase { knSphStompIsolatedParticleValue( ParticleDataImpl<Real> &v, const SphWorld &sph, const int itype) :  KernelBase(v.size()) ,v(v),sph(sph),itype(itype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &v, const SphWorld &sph, const int itype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	if(n.size(idx)<=1) v[idx] = 0.0;
}    inline ParticleDataImpl<Real> & getArg0() { return v; } typedef ParticleDataImpl<Real>  type0;inline const SphWorld& getArg1() { return sph; } typedef SphWorld type1;inline const int& getArg2() { return itype; } typedef int type2; void runMessage() { debMsg("Executing kernel knSphStompIsolatedParticleValue ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,v,sph,itype);  }   } ParticleDataImpl<Real> & v; const SphWorld& sph; const int itype;   };
#line 593 "plugin/sph.cpp"



void sphStompIsolatedParticleValue(ParticleDataImpl<Real> &v, const SphWorld &sph, const int itype=SphWorld::FlagFluid) {
	knSphStompIsolatedParticleValue(v, sph, itype);
} static PyObject* _W_19 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphStompIsolatedParticleValue" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & v = *_args.getPtr<ParticleDataImpl<Real>  >("v",0,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",1,&_lock); const int itype = _args.getOpt<int >("itype",2,SphWorld::FlagFluid,&_lock);   _retval = getPyNone(); sphStompIsolatedParticleValue(v,sph,itype);  _args.check(); } pbFinalizePlugin(parent,"sphStompIsolatedParticleValue", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphStompIsolatedParticleValue",e.what()); return 0; } } static const Pb::Register _RP_sphStompIsolatedParticleValue ("","sphStompIsolatedParticleValue",_W_19);  extern "C" { void PbRegister_sphStompIsolatedParticleValue() { KEEP_UNUSED(_RP_sphStompIsolatedParticleValue); } } 




 struct knSphComputeConstantForce : public KernelBase { knSphComputeConstantForce( ParticleDataImpl<Vec3> &f, const Vec3 &v, const SphWorld &sph, const int itype) :  KernelBase(f.size()) ,f(f),v(v),sph(sph),itype(itype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &f, const Vec3 &v, const SphWorld &sph, const int itype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	f[idx] += v;
}    inline ParticleDataImpl<Vec3> & getArg0() { return f; } typedef ParticleDataImpl<Vec3>  type0;inline const Vec3& getArg1() { return v; } typedef Vec3 type1;inline const SphWorld& getArg2() { return sph; } typedef SphWorld type2;inline const int& getArg3() { return itype; } typedef int type3; void runMessage() { debMsg("Executing kernel knSphComputeConstantForce ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,f,v,sph,itype);  }   } ParticleDataImpl<Vec3> & f; const Vec3& v; const SphWorld& sph; const int itype;   };
#line 609 "plugin/sph.cpp"







void sphComputeConstantForce( ParticleDataImpl<Vec3> &f, const Vec3 &v, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	knSphComputeConstantForce(f, v, sph, itype);
} static PyObject* _W_20 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeConstantForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const Vec3& v = _args.get<Vec3 >("v",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); const int itype = _args.getOpt<int >("itype",3,SphWorld::FlagFluid,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",4,true,&_lock);   _retval = getPyNone(); sphComputeConstantForce(f,v,sph,itype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeConstantForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeConstantForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputeConstantForce ("","sphComputeConstantForce",_W_20);  extern "C" { void PbRegister_sphComputeConstantForce() { KEEP_UNUSED(_RP_sphComputeConstantForce); } } 



template <typename T_Pij>  struct knSphComputePressureForce : public KernelBase { knSphComputePressureForce( ParticleDataImpl<Vec3> &f, const T_Pij &Pij, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(f.size()) ,f(f),Pij(Pij),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &f, const T_Pij &Pij, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	Vec3 fi(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		fi -= k.grad_w(xi - xj, len_xij)*Pij(idx, j);
	}
	f[idx] += fi;
}    inline ParticleDataImpl<Vec3> & getArg0() { return f; } typedef ParticleDataImpl<Vec3>  type0;inline const T_Pij& getArg1() { return Pij; } typedef T_Pij type1;inline const CubicSpline& getArg2() { return k; } typedef CubicSpline type2;inline const SphWorld& getArg3() { return sph; } typedef SphWorld type3;inline const int& getArg4() { return itype; } typedef int type4;inline const int& getArg5() { return jtype; } typedef int type5; void runMessage() { debMsg("Executing kernel knSphComputePressureForce ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,f,Pij,k,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & f; const T_Pij& Pij; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 628 "plugin/sph.cpp"








void sphComputePressureForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	if(m) knSphComputePressureForce<PijMij>(f, PijMij(p, d, *m, sph.particleType()), k, sph, itype, jtype);
	else knSphComputePressureForce<Pij>(f, Pij(p, d, sph._m0, sph.particleType()), k, sph, itype, jtype);
} static PyObject* _W_21 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputePressureForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",1,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",5,NULL,&_lock); const int itype = _args.getOpt<int >("itype",6,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",7,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",8,true,&_lock);   _retval = getPyNone(); sphComputePressureForce(f,p,d,k,sph,m,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputePressureForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputePressureForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputePressureForce ("","sphComputePressureForce",_W_21);  extern "C" { void PbRegister_sphComputePressureForce() { KEEP_UNUSED(_RP_sphComputePressureForce); } } 






void sphComputeDfsphPressureForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	if(m) knSphComputePressureForce<PijMijDfsph>(f, PijMijDfsph(p, d, *m, sph.particleType()), k, sph, itype, jtype);
	else knSphComputePressureForce<PijDfsph>(f, PijDfsph(p, d, sph._m0, sph.particleType()), k, sph, itype, jtype);
} static PyObject* _W_22 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeDfsphPressureForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",1,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",5,NULL,&_lock); const int itype = _args.getOpt<int >("itype",6,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",7,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",8,true,&_lock);   _retval = getPyNone(); sphComputeDfsphPressureForce(f,p,d,k,sph,m,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeDfsphPressureForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeDfsphPressureForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputeDfsphPressureForce ("","sphComputeDfsphPressureForce",_W_22);  extern "C" { void PbRegister_sphComputeDfsphPressureForce() { KEEP_UNUSED(_RP_sphComputeDfsphPressureForce); } } 







void sphComputeDensityWeightedPressureForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	knSphComputePressureForce<PijWeighted>(f, PijWeighted(p, d, Vsqr), k, sph, itype, jtype);
} static PyObject* _W_23 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeDensityWeightedPressureForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",1,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",2,&_lock); const ParticleDataImpl<Real> & Vsqr = *_args.getPtr<ParticleDataImpl<Real>  >("Vsqr",3,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",4,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",5,&_lock); const int itype = _args.getOpt<int >("itype",6,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",7,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",8,true,&_lock);   _retval = getPyNone(); sphComputeDensityWeightedPressureForce(f,p,d,Vsqr,k,sph,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeDensityWeightedPressureForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeDensityWeightedPressureForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputeDensityWeightedPressureForce ("","sphComputeDensityWeightedPressureForce",_W_23);  extern "C" { void PbRegister_sphComputeDensityWeightedPressureForce() { KEEP_UNUSED(_RP_sphComputeDensityWeightedPressureForce); } } 






void sphComputeBackgroundPressureForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	knSphComputePressureForce<PijBg>(f, PijBg(p, Vsqr), k, sph, itype, jtype);
} static PyObject* _W_24 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeBackgroundPressureForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",1,&_lock); const ParticleDataImpl<Real> & Vsqr = *_args.getPtr<ParticleDataImpl<Real>  >("Vsqr",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",6,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",7,true,&_lock);   _retval = getPyNone(); sphComputeBackgroundPressureForce(f,p,Vsqr,k,sph,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeBackgroundPressureForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeBackgroundPressureForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputeBackgroundPressureForce ("","sphComputeBackgroundPressureForce",_W_24);  extern "C" { void PbRegister_sphComputeBackgroundPressureForce() { KEEP_UNUSED(_RP_sphComputeBackgroundPressureForce); } } 



 struct knSphComputeBoundaryForce : public KernelBase { knSphComputeBoundaryForce( ParticleDataImpl<Vec3> &f, const BndKernel &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(f.size()) ,f(f),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &f, const BndKernel &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	Vec3 fi(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(k.supportRadius()<len_xij) continue;
		const Vec3 &xj = pts[j].pos;
		fi += (xi - xj)*k.f(len_xij)/len_xij; // NOTE: [Monaghan et al. 2005] used slightly different one.
	}
	f[idx] += fi*0.5;
}    inline ParticleDataImpl<Vec3> & getArg0() { return f; } typedef ParticleDataImpl<Vec3>  type0;inline const BndKernel& getArg1() { return k; } typedef BndKernel type1;inline const SphWorld& getArg2() { return sph; } typedef SphWorld type2;inline const int& getArg3() { return itype; } typedef int type3;inline const int& getArg4() { return jtype; } typedef int type4; void runMessage() { debMsg("Executing kernel knSphComputeBoundaryForce ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,f,k,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & f; const BndKernel& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 691 "plugin/sph.cpp"








void sphComputeBoundaryForce( ParticleDataImpl<Vec3> &f, const BndKernel &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	knSphComputeBoundaryForce(f, k, sph, itype, jtype);
} static PyObject* _W_25 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeBoundaryForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const BndKernel& k = *_args.getPtr<BndKernel >("k",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); const int itype = _args.getOpt<int >("itype",3,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",4,SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",5,true,&_lock);   _retval = getPyNone(); sphComputeBoundaryForce(f,k,sph,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeBoundaryForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeBoundaryForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputeBoundaryForce ("","sphComputeBoundaryForce",_W_25);  extern "C" { void PbRegister_sphComputeBoundaryForce() { KEEP_UNUSED(_RP_sphComputeBoundaryForce); } } 




 struct knSphSetBndPressure : public KernelBase { knSphSetBndPressure( ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const Real pb, const Real d_th, const SphWorld &sph, const int itype) :  KernelBase(p.size()) ,p(p),d(d),pb(pb),d_th(d_th),sph(sph),itype(itype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const Real pb, const Real d_th, const SphWorld &sph, const int itype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;
	if(d[idx]<d_th) p[idx] = pb;
}    inline ParticleDataImpl<Real> & getArg0() { return p; } typedef ParticleDataImpl<Real>  type0;inline const ParticleDataImpl<Real> & getArg1() { return d; } typedef ParticleDataImpl<Real>  type1;inline const Real& getArg2() { return pb; } typedef Real type2;inline const Real& getArg3() { return d_th; } typedef Real type3;inline const SphWorld& getArg4() { return sph; } typedef SphWorld type4;inline const int& getArg5() { return itype; } typedef int type5; void runMessage() { debMsg("Executing kernel knSphSetBndPressure ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,p,d,pb,d_th,sph,itype);  }   } ParticleDataImpl<Real> & p; const ParticleDataImpl<Real> & d; const Real pb; const Real d_th; const SphWorld& sph; const int itype;   };
#line 722 "plugin/sph.cpp"







void sphSetBndPressure( ParticleDataImpl<Real> &p, const ParticleDataImpl<Real> &d, const Real pb, const Real d_th, const SphWorld &sph, const int itype=SphWorld::FlagFluid) {
	knSphSetBndPressure(p, d, pb, d_th, sph, itype);
} static PyObject* _W_26 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphSetBndPressure" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Real> & p = *_args.getPtr<ParticleDataImpl<Real>  >("p",0,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",1,&_lock); const Real pb = _args.get<Real >("pb",2,&_lock); const Real d_th = _args.get<Real >("d_th",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagFluid,&_lock);   _retval = getPyNone(); sphSetBndPressure(p,d,pb,d_th,sph,itype);  _args.check(); } pbFinalizePlugin(parent,"sphSetBndPressure", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphSetBndPressure",e.what()); return 0; } } static const Pb::Register _RP_sphSetBndPressure ("","sphSetBndPressure",_W_26);  extern "C" { void PbRegister_sphSetBndPressure() { KEEP_UNUSED(_RP_sphSetBndPressure); } } 



 struct knSphComputeArtificialViscousForce : public KernelBase { knSphComputeArtificialViscousForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(f.size()) ,f(f),v(v),d(d),m(m),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos, &vi = v[idx];
	const Real di = d[idx], mi = (m) ? m->get(idx) : sph._m0;
	Vec3 fi(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos, &vj = v[j];
		const Vec3 xij = xi - xj, vij = vi - vj;
		const Real v_dot_x = dot(vij, xij) + ((pts.getParent()->is3D()) ? 0.0 : -vij.z*xij.z);
		if(v_dot_x<Real(0)) {
			const Real dj = ((t[j]&SphWorld::FlagObstacle) ? di : d[j]), mj = (m) ? m->get(j) : sph._m0;
			const Real nu = 2.0*sph._alpha*k.radius()*sph._c/(di+dj);
			fi += k.grad_w(xi - xj, len_xij)*mi*mj*nu*v_dot_x/(square(len_xij) + 0.01*square(k.radius()));
		}
	}
	f[idx] += fi;
}    inline ParticleDataImpl<Vec3> & getArg0() { return f; } typedef ParticleDataImpl<Vec3>  type0;inline const ParticleDataImpl<Vec3> & getArg1() { return v; } typedef ParticleDataImpl<Vec3>  type1;inline const ParticleDataImpl<Real> & getArg2() { return d; } typedef ParticleDataImpl<Real>  type2;inline const ParticleDataImpl<Real> * getArg3() { return m; } typedef ParticleDataImpl<Real>  type3;inline const CubicSpline& getArg4() { return k; } typedef CubicSpline type4;inline const SphWorld& getArg5() { return sph; } typedef SphWorld type5;inline const int& getArg6() { return itype; } typedef int type6;inline const int& getArg7() { return jtype; } typedef int type7; void runMessage() { debMsg("Executing kernel knSphComputeArtificialViscousForce ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,f,v,d,m,k,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & f; const ParticleDataImpl<Vec3> & v; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 740 "plugin/sph.cpp"








void sphComputeArtificialViscousForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Real> &d, const CubicSpline &k, const SphWorld &sph, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	knSphComputeArtificialViscousForce(f, v, d, m, k, sph, itype, jtype);
} static PyObject* _W_27 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeArtificialViscousForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",1,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",2,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",3,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",4,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",5,NULL,&_lock); const int itype = _args.getOpt<int >("itype",6,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",7,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",8,true,&_lock);   _retval = getPyNone(); sphComputeArtificialViscousForce(f,v,d,k,sph,m,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeArtificialViscousForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeArtificialViscousForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputeArtificialViscousForce ("","sphComputeArtificialViscousForce",_W_27);  extern "C" { void PbRegister_sphComputeArtificialViscousForce() { KEEP_UNUSED(_RP_sphComputeArtificialViscousForce); } } 
// [Becker and Teschner 2007]




 struct knSphComputeSurfTension : public KernelBase { knSphComputeSurfTension( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real kappa, const int itype, const int jtype) :  KernelBase(f.size()) ,f(f),m(m),k(k),sph(sph),kappa(kappa),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Real> *m, const CubicSpline &k, const SphWorld &sph, const Real kappa, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	Vec3 sum_f(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		const Vec3 xij = xi - xj;
		sum_f += xij*((m) ? m->get(j) : sph._m0)*k.f(len_xij);
	}
	f[idx] -= sum_f*kappa;
}    inline ParticleDataImpl<Vec3> & getArg0() { return f; } typedef ParticleDataImpl<Vec3>  type0;inline const ParticleDataImpl<Real> * getArg1() { return m; } typedef ParticleDataImpl<Real>  type1;inline const CubicSpline& getArg2() { return k; } typedef CubicSpline type2;inline const SphWorld& getArg3() { return sph; } typedef SphWorld type3;inline const Real& getArg4() { return kappa; } typedef Real type4;inline const int& getArg5() { return itype; } typedef int type5;inline const int& getArg6() { return jtype; } typedef int type6; void runMessage() { debMsg("Executing kernel knSphComputeSurfTension ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,f,m,k,sph,kappa,itype,jtype);  }   } ParticleDataImpl<Vec3> & f; const ParticleDataImpl<Real> * m; const CubicSpline& k; const SphWorld& sph; const Real kappa; const int itype; const int jtype;   };
#line 779 "plugin/sph.cpp"








void sphComputeSurfTension( ParticleDataImpl<Vec3> &f, const CubicSpline &k, const SphWorld &sph, const Real kappa, const ParticleDataImpl<Real> *m=NULL, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	knSphComputeSurfTension(f, m, k, sph, kappa, itype, jtype);
} static PyObject* _W_28 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeSurfTension" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); const Real kappa = _args.get<Real >("kappa",3,&_lock); const ParticleDataImpl<Real> * m = _args.getPtrOpt<ParticleDataImpl<Real>  >("m",4,NULL,&_lock); const int itype = _args.getOpt<int >("itype",5,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",6,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",7,true,&_lock);   _retval = getPyNone(); sphComputeSurfTension(f,k,sph,kappa,m,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeSurfTension", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeSurfTension",e.what()); return 0; } } static const Pb::Register _RP_sphComputeSurfTension ("","sphComputeSurfTension",_W_28);  extern "C" { void PbRegister_sphComputeSurfTension() { KEEP_UNUSED(_RP_sphComputeSurfTension); } } 





 struct knSphComputeTensorForce : public KernelBase { knSphComputeTensorForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &vv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(f.size()) ,f(f),v(v),vv(vv),d(d),Vsqr(Vsqr),k(k),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &vv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	const Real Vsqr_i = Vsqr[idx];
	const Matrix3x3f A_i = outerProduct(v[idx], vv[idx]-v[idx])*d[idx];

	Vec3 fi(0.0);
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j|| !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;

		const Vec3 &xj = pts[j].pos;
		const Real deriv_w = k.derivative_f(len_xij);
		const Real Vsqr_j = Vsqr[j];
		Vec3 sum_Vsqr_g_w = (xi - xj)*(deriv_w*(Vsqr_i+Vsqr_j)/len_xij);
		if(pts.getParent()->is2D()) sum_Vsqr_g_w.z = 0.0;

		const Matrix3x3f A_j = outerProduct(v[j], vv[j]-v[j])*d[j];

		fi += 0.5*(A_i+A_j)*sum_Vsqr_g_w;
	}
	if(pts.getParent()->is2D()) fi.z = 0.0;
	f[idx] += fi;
}    inline ParticleDataImpl<Vec3> & getArg0() { return f; } typedef ParticleDataImpl<Vec3>  type0;inline const ParticleDataImpl<Vec3> & getArg1() { return v; } typedef ParticleDataImpl<Vec3>  type1;inline const ParticleDataImpl<Vec3> & getArg2() { return vv; } typedef ParticleDataImpl<Vec3>  type2;inline const ParticleDataImpl<Real> & getArg3() { return d; } typedef ParticleDataImpl<Real>  type3;inline const ParticleDataImpl<Real> & getArg4() { return Vsqr; } typedef ParticleDataImpl<Real>  type4;inline const CubicSpline& getArg5() { return k; } typedef CubicSpline type5;inline const SphWorld& getArg6() { return sph; } typedef SphWorld type6;inline const int& getArg7() { return itype; } typedef int type7;inline const int& getArg8() { return jtype; } typedef int type8; void runMessage() { debMsg("Executing kernel knSphComputeTensorForce ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,f,v,vv,d,Vsqr,k,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & f; const ParticleDataImpl<Vec3> & v; const ParticleDataImpl<Vec3> & vv; const ParticleDataImpl<Real> & d; const ParticleDataImpl<Real> & Vsqr; const CubicSpline& k; const SphWorld& sph; const int itype; const int jtype;   };
#line 812 "plugin/sph.cpp"









void sphComputeTensorForce( ParticleDataImpl<Vec3> &f, const ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &vv, const ParticleDataImpl<Real> &d, const ParticleDataImpl<Real> &Vsqr, const CubicSpline &k, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle, const bool accumulate=true) {
	if(!accumulate) f.setConst(Vec3(0.0));
	knSphComputeTensorForce(f, v, vv, d, Vsqr, k, sph, itype, jtype);
} static PyObject* _W_29 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphComputeTensorForce" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & f = *_args.getPtr<ParticleDataImpl<Vec3>  >("f",0,&_lock); const ParticleDataImpl<Vec3> & v = *_args.getPtr<ParticleDataImpl<Vec3>  >("v",1,&_lock); const ParticleDataImpl<Vec3> & vv = *_args.getPtr<ParticleDataImpl<Vec3>  >("vv",2,&_lock); const ParticleDataImpl<Real> & d = *_args.getPtr<ParticleDataImpl<Real>  >("d",3,&_lock); const ParticleDataImpl<Real> & Vsqr = *_args.getPtr<ParticleDataImpl<Real>  >("Vsqr",4,&_lock); const CubicSpline& k = *_args.getPtr<CubicSpline >("k",5,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",6,&_lock); const int itype = _args.getOpt<int >("itype",7,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",8,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock); const bool accumulate = _args.getOpt<bool >("accumulate",9,true,&_lock);   _retval = getPyNone(); sphComputeTensorForce(f,v,vv,d,Vsqr,k,sph,itype,jtype,accumulate);  _args.check(); } pbFinalizePlugin(parent,"sphComputeTensorForce", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphComputeTensorForce",e.what()); return 0; } } static const Pb::Register _RP_sphComputeTensorForce ("","sphComputeTensorForce",_W_29);  extern "C" { void PbRegister_sphComputeTensorForce() { KEEP_UNUSED(_RP_sphComputeTensorForce); } } 




 struct knSphAddPositionCorrection : public KernelBase { knSphAddPositionCorrection( ParticleDataImpl<Vec3> &dp, const BasicParticleSystem &x, const SphWorld &sph, const int itype, const int jtype) :  KernelBase(dp.size()) ,dp(dp),x(x),sph(sph),itype(itype),jtype(jtype)   { runMessage(); run(); }   inline void op(IndexInt idx,  ParticleDataImpl<Vec3> &dp, const BasicParticleSystem &x, const SphWorld &sph, const int itype, const int jtype )  {
	const BasicParticleSystem	&pts = sph.particleSystem();
	const ParticleDataImpl<int>	&t   = sph.particleType();
	const ParticleNeighbors		&n   = sph.neighborData();

	if(!pts.isActive(idx) || !(t[idx]&itype)) return;

	const Vec3 &xi = pts[idx].pos;
	for(ParticleNeighbors::Neighbors::const_iterator it=n.begin(idx); it!=n.end(idx); ++it) {
		const int j = n.neighborIdx(it);   if(idx==j || !(t[j]&jtype)) continue;
		const Real len_xij = n.length(it); if(len_xij<=Real(0)) continue;
		const Vec3 &xj = pts[j].pos;
		const Vec3 xij = xi - xj;
		dp[idx] += xij*wSmooth(square(len_xij), square(sph._delta))/len_xij;
	}
}    inline ParticleDataImpl<Vec3> & getArg0() { return dp; } typedef ParticleDataImpl<Vec3>  type0;inline const BasicParticleSystem& getArg1() { return x; } typedef BasicParticleSystem type1;inline const SphWorld& getArg2() { return sph; } typedef SphWorld type2;inline const int& getArg3() { return itype; } typedef int type3;inline const int& getArg4() { return jtype; } typedef int type4; void runMessage() { debMsg("Executing kernel knSphAddPositionCorrection ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,dp,x,sph,itype,jtype);  }   } ParticleDataImpl<Vec3> & dp; const BasicParticleSystem& x; const SphWorld& sph; const int itype; const int jtype;   };
#line 856 "plugin/sph.cpp"






void sphAddPositionCorrection( ParticleDataImpl<Vec3> &dp, const BasicParticleSystem &x, const SphWorld &sph, const int itype=SphWorld::FlagFluid, const int jtype=SphWorld::FlagFluid|SphWorld::FlagObstacle) {
	knSphAddPositionCorrection(dp, x, sph, itype, jtype);
} static PyObject* _W_30 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphAddPositionCorrection" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & dp = *_args.getPtr<ParticleDataImpl<Vec3>  >("dp",0,&_lock); const BasicParticleSystem& x = *_args.getPtr<BasicParticleSystem >("x",1,&_lock); const SphWorld& sph = *_args.getPtr<SphWorld >("sph",2,&_lock); const int itype = _args.getOpt<int >("itype",3,SphWorld::FlagFluid,&_lock); const int jtype = _args.getOpt<int >("jtype",4,SphWorld::FlagFluid|SphWorld::FlagObstacle,&_lock);   _retval = getPyNone(); sphAddPositionCorrection(dp,x,sph,itype,jtype);  _args.check(); } pbFinalizePlugin(parent,"sphAddPositionCorrection", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphAddPositionCorrection",e.what()); return 0; } } static const Pb::Register _RP_sphAddPositionCorrection ("","sphAddPositionCorrection",_W_30);  extern "C" { void PbRegister_sphAddPositionCorrection() { KEEP_UNUSED(_RP_sphAddPositionCorrection); } } 


int sphCountParticles(const ParticleDataImpl<int> &t, const int type) {
	int cnt=0;
	for(int i=0; i<t.size(); ++i) if(t[i]&type) ++cnt;
	return cnt;
} static PyObject* _W_31 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "sphCountParticles" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; const ParticleDataImpl<int> & t = *_args.getPtr<ParticleDataImpl<int>  >("t",0,&_lock); const int type = _args.get<int >("type",1,&_lock);   _retval = toPy(sphCountParticles(t,type));  _args.check(); } pbFinalizePlugin(parent,"sphCountParticles", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("sphCountParticles",e.what()); return 0; } } static const Pb::Register _RP_sphCountParticles ("","sphCountParticles",_W_31);  extern "C" { void PbRegister_sphCountParticles() { KEEP_UNUSED(_RP_sphCountParticles); } } 

} // namespace


