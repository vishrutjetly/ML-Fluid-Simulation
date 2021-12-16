




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip/source/plugin/tfplugins.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2017 Kiwon Um, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * tensorflor/numpy plugins, mostly for MLFLIP
 * [https://arxiv.org/abs/1704.04456] for now (only compiled if NUMPY is
 * enabled)
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include <cmath>

namespace Manta {

//! simple test kernel and kernel with numpy array


 struct knSimpleNumpyTest : public KernelBase { knSimpleNumpyTest(Grid<Real>& grid, PyArrayContainer npAr, Real scale) :  KernelBase(&grid,0) ,grid(grid),npAr(npAr),scale(scale)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real>& grid, PyArrayContainer npAr, Real scale )  {
	const float* p = reinterpret_cast<float*>(npAr.pData);
	grid(i,j,k) += scale * (Real)p[j*grid.getSizeX()+i]; // calc access into numpy array, no size check here!
}   inline Grid<Real>& getArg0() { return grid; } typedef Grid<Real> type0;inline PyArrayContainer& getArg1() { return npAr; } typedef PyArrayContainer type1;inline Real& getArg2() { return scale; } typedef Real type2; void runMessage() { debMsg("Executing kernel knSimpleNumpyTest ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,grid,npAr,scale);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,grid,npAr,scale);  } }  } Grid<Real>& grid; PyArrayContainer npAr; Real scale;   };
#line 26 "plugin/tfplugins.cpp"



//! simple test function and kernel with numpy array

void simpleNumpyTest( Grid<Real>& grid, PyArrayContainer npAr, Real scale) {
	knSimpleNumpyTest(grid, npAr, scale);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "simpleNumpyTest" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& grid = *_args.getPtr<Grid<Real> >("grid",0,&_lock); PyArrayContainer npAr = _args.get<PyArrayContainer >("npAr",1,&_lock); Real scale = _args.get<Real >("scale",2,&_lock);   _retval = getPyNone(); simpleNumpyTest(grid,npAr,scale);  _args.check(); } pbFinalizePlugin(parent,"simpleNumpyTest", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("simpleNumpyTest",e.what()); return 0; } } static const Pb::Register _RP_simpleNumpyTest ("","simpleNumpyTest",_W_0);  extern "C" { void PbRegister_simpleNumpyTest() { KEEP_UNUSED(_RP_simpleNumpyTest); } } 

//! extract feature vectors





 struct knExtractFeatureVel : public KernelBase { knExtractFeatureVel( const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const MACGrid &vel, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h) :  KernelBase(p.size()) ,p(p),fv(fv),N_row(N_row),off_begin(off_begin),vel(vel),scale(scale),ptype(ptype),exclude(exclude),window(window),h(h)   { runMessage(); run(); }   inline void op(IndexInt idx,  const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const MACGrid &vel, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h )  {
	if(!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	const int _k = (vel.is3D()) ? -window : 0, K = (vel.is3D()) ? window : 0;
	const IndexInt D = (vel.is3D()) ? 3 : 2;
	const IndexInt off_idx = idx*N_row;

	IndexInt off_stencil = 0;
	for(int i=-window; i<=window; ++i) {
		for(int j=-window; j<=window; ++j) {
			for(int k=_k; k<=K; ++k) {
				const Vec3 off_pos(static_cast<Real>(i)*h, static_cast<Real>(j)*h, static_cast<Real>(k)*h);
				const Vec3 pos_s = p[idx].pos + off_pos;

				const Vec3 vel_s = vel.getInterpolated(pos_s)*scale;
				const IndexInt off_vel = off_idx + off_begin + off_stencil*D;
				fv[off_vel + 0] = vel_s[0];
				fv[off_vel + 1] = vel_s[1];
				if(vel.is3D()) fv[off_vel + 2] = vel_s[2];

				++off_stencil;
			}
		}
	}
}    inline const BasicParticleSystem& getArg0() { return p; } typedef BasicParticleSystem type0;inline Real* getArg1() { return fv; } typedef Real type1;inline const IndexInt& getArg2() { return N_row; } typedef IndexInt type2;inline const IndexInt& getArg3() { return off_begin; } typedef IndexInt type3;inline const MACGrid& getArg4() { return vel; } typedef MACGrid type4;inline const Real& getArg5() { return scale; } typedef Real type5;inline const ParticleDataImpl<int> * getArg6() { return ptype; } typedef ParticleDataImpl<int>  type6;inline const int& getArg7() { return exclude; } typedef int type7;inline const int& getArg8() { return window; } typedef int type8;inline const Real& getArg9() { return h; } typedef Real type9; void runMessage() { debMsg("Executing kernel knExtractFeatureVel ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,p,fv,N_row,off_begin,vel,scale,ptype,exclude,window,h);  }   } const BasicParticleSystem& p; Real* fv; const IndexInt N_row; const IndexInt off_begin; const MACGrid& vel; const Real scale; const ParticleDataImpl<int> * ptype; const int exclude; const int window; const Real h;   };
#line 43 "plugin/tfplugins.cpp"






 struct knExtractFeatureVelRel : public KernelBase { knExtractFeatureVelRel( const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const MACGrid &vel, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h) :  KernelBase(p.size()) ,p(p),fv(fv),N_row(N_row),off_begin(off_begin),vel(vel),scale(scale),ptype(ptype),exclude(exclude),window(window),h(h)   { runMessage(); run(); }   inline void op(IndexInt idx,  const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const MACGrid &vel, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h )  {
	if(!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	const int _k = (vel.is3D()) ? -window : 0, K = (vel.is3D()) ? window : 0;
	const IndexInt D = (vel.is3D()) ? 3 : 2;
	const IndexInt off_idx = idx*N_row;
	const Vec3 v0 = vel.getInterpolated(p[idx].pos);

	IndexInt off_stencil = 0;
	for(int i=-window; i<=window; ++i) {
		for(int j=-window; j<=window; ++j) {
			for(int k=_k; k<=K; ++k) {
				const Vec3 off_pos(static_cast<Real>(i)*h, static_cast<Real>(j)*h, static_cast<Real>(k)*h);
				const Vec3 pos_s = p[idx].pos + off_pos;

				const Vec3 vel_s = (vel.getInterpolated(pos_s) - v0)*scale;
				const IndexInt off_vel = off_idx + off_begin + off_stencil*D;
				fv[off_vel + 0] = vel_s[0];
				fv[off_vel + 1] = vel_s[1];
				if(vel.is3D()) fv[off_vel + 2] = vel_s[2];

				if(i==0 && j==0 && k==0) {
					fv[off_vel + 0] = scale;
					fv[off_vel + 1] = scale;
					if(vel.is3D()) fv[off_vel + 2] = scale;
				}

				++off_stencil;
			}
		}
	}
}    inline const BasicParticleSystem& getArg0() { return p; } typedef BasicParticleSystem type0;inline Real* getArg1() { return fv; } typedef Real type1;inline const IndexInt& getArg2() { return N_row; } typedef IndexInt type2;inline const IndexInt& getArg3() { return off_begin; } typedef IndexInt type3;inline const MACGrid& getArg4() { return vel; } typedef MACGrid type4;inline const Real& getArg5() { return scale; } typedef Real type5;inline const ParticleDataImpl<int> * getArg6() { return ptype; } typedef ParticleDataImpl<int>  type6;inline const int& getArg7() { return exclude; } typedef int type7;inline const int& getArg8() { return window; } typedef int type8;inline const Real& getArg9() { return h; } typedef Real type9; void runMessage() { debMsg("Executing kernel knExtractFeatureVelRel ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,p,fv,N_row,off_begin,vel,scale,ptype,exclude,window,h);  }   } const BasicParticleSystem& p; Real* fv; const IndexInt N_row; const IndexInt off_begin; const MACGrid& vel; const Real scale; const ParticleDataImpl<int> * ptype; const int exclude; const int window; const Real h;   };
#line 72 "plugin/tfplugins.cpp"







 struct knExtractFeaturePhi : public KernelBase { knExtractFeaturePhi( const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const Grid<Real> &phi, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h) :  KernelBase(p.size()) ,p(p),fv(fv),N_row(N_row),off_begin(off_begin),phi(phi),scale(scale),ptype(ptype),exclude(exclude),window(window),h(h)   { runMessage(); run(); }   inline void op(IndexInt idx,  const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const Grid<Real> &phi, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h )  {
	if(!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	const int _k = (phi.is3D()) ? -window : 0, K = (phi.is3D()) ? window : 0;
	const IndexInt off_idx = idx*N_row;

	IndexInt off_stencil = 0;
	for(int i=-window; i<=window; ++i) {
		for(int j=-window; j<=window; ++j) {
			for(int k=_k; k<=K; ++k) {
				const Vec3 off_pos(static_cast<Real>(i)*h, static_cast<Real>(j)*h, static_cast<Real>(k)*h);
				const Vec3 pos_s = p[idx].pos + off_pos;

				const Real phi_s = phi.getInterpolated(pos_s)*scale;
				const IndexInt off_phi = off_idx + off_begin + off_stencil;
				fv[off_phi] = phi_s;

				++off_stencil;
			}
		}
	}
}    inline const BasicParticleSystem& getArg0() { return p; } typedef BasicParticleSystem type0;inline Real* getArg1() { return fv; } typedef Real type1;inline const IndexInt& getArg2() { return N_row; } typedef IndexInt type2;inline const IndexInt& getArg3() { return off_begin; } typedef IndexInt type3;inline const Grid<Real> & getArg4() { return phi; } typedef Grid<Real>  type4;inline const Real& getArg5() { return scale; } typedef Real type5;inline const ParticleDataImpl<int> * getArg6() { return ptype; } typedef ParticleDataImpl<int>  type6;inline const int& getArg7() { return exclude; } typedef int type7;inline const int& getArg8() { return window; } typedef int type8;inline const Real& getArg9() { return h; } typedef Real type9; void runMessage() { debMsg("Executing kernel knExtractFeaturePhi ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,p,fv,N_row,off_begin,phi,scale,ptype,exclude,window,h);  }   } const BasicParticleSystem& p; Real* fv; const IndexInt N_row; const IndexInt off_begin; const Grid<Real> & phi; const Real scale; const ParticleDataImpl<int> * ptype; const int exclude; const int window; const Real h;   };
#line 109 "plugin/tfplugins.cpp"







 struct knExtractFeatureGeo : public KernelBase { knExtractFeatureGeo( const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const FlagGrid &geo, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h) :  KernelBase(p.size()) ,p(p),fv(fv),N_row(N_row),off_begin(off_begin),geo(geo),scale(scale),ptype(ptype),exclude(exclude),window(window),h(h)   { runMessage(); run(); }   inline void op(IndexInt idx,  const BasicParticleSystem &p, Real *fv, const IndexInt N_row, const IndexInt off_begin, const FlagGrid &geo, const Real scale, const ParticleDataImpl<int> *ptype, const int exclude, const int window, const Real h )  {
	if(!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;

	const int _k = (geo.is3D()) ? -window : 0, K = (geo.is3D()) ? window : 0;
	const IndexInt off_idx = idx*N_row;

	IndexInt off_stencil = 0;
	for(int i=-window; i<=window; ++i) {
		for(int j=-window; j<=window; ++j) {
			for(int k=_k; k<=K; ++k) {
				const Vec3 off_pos(static_cast<Real>(i)*h, static_cast<Real>(j)*h, static_cast<Real>(k)*h);
				const Vec3 pos_s = p[idx].pos + off_pos;

				const Real geo_s = static_cast<Real>(geo.getAt(pos_s))*scale;
				const IndexInt off_geo = off_idx + off_begin + off_stencil;
				fv[off_geo] = geo_s;

				++off_stencil;
			}
		}
	}
}    inline const BasicParticleSystem& getArg0() { return p; } typedef BasicParticleSystem type0;inline Real* getArg1() { return fv; } typedef Real type1;inline const IndexInt& getArg2() { return N_row; } typedef IndexInt type2;inline const IndexInt& getArg3() { return off_begin; } typedef IndexInt type3;inline const FlagGrid& getArg4() { return geo; } typedef FlagGrid type4;inline const Real& getArg5() { return scale; } typedef Real type5;inline const ParticleDataImpl<int> * getArg6() { return ptype; } typedef ParticleDataImpl<int>  type6;inline const int& getArg7() { return exclude; } typedef int type7;inline const int& getArg8() { return window; } typedef int type8;inline const Real& getArg9() { return h; } typedef Real type9; void runMessage() { debMsg("Executing kernel knExtractFeatureGeo ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,p,fv,N_row,off_begin,geo,scale,ptype,exclude,window,h);  }   } const BasicParticleSystem& p; Real* fv; const IndexInt N_row; const IndexInt off_begin; const FlagGrid& geo; const Real scale; const ParticleDataImpl<int> * ptype; const int exclude; const int window; const Real h;   };
#line 136 "plugin/tfplugins.cpp"








void extractFeatureVel( PyArrayContainer fv, const int N_row, const int off_begin, const BasicParticleSystem &p, const MACGrid &vel, const Real scale=1.0, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0, const int window=1, const Real h=1.0) {
	knExtractFeatureVel(p, reinterpret_cast<Real*>(fv.pData), N_row, off_begin, vel, scale, ptype, exclude, window, h);
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "extractFeatureVel" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; PyArrayContainer fv = _args.get<PyArrayContainer >("fv",0,&_lock); const int N_row = _args.get<int >("N_row",1,&_lock); const int off_begin = _args.get<int >("off_begin",2,&_lock); const BasicParticleSystem& p = *_args.getPtr<BasicParticleSystem >("p",3,&_lock); const MACGrid& vel = *_args.getPtr<MACGrid >("vel",4,&_lock); const Real scale = _args.getOpt<Real >("scale",5,1.0,&_lock); const ParticleDataImpl<int> * ptype = _args.getPtrOpt<ParticleDataImpl<int>  >("ptype",6,NULL,&_lock); const int exclude = _args.getOpt<int >("exclude",7,0,&_lock); const int window = _args.getOpt<int >("window",8,1,&_lock); const Real h = _args.getOpt<Real >("h",9,1.0,&_lock);   _retval = getPyNone(); extractFeatureVel(fv,N_row,off_begin,p,vel,scale,ptype,exclude,window,h);  _args.check(); } pbFinalizePlugin(parent,"extractFeatureVel", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("extractFeatureVel",e.what()); return 0; } } static const Pb::Register _RP_extractFeatureVel ("","extractFeatureVel",_W_1);  extern "C" { void PbRegister_extractFeatureVel() { KEEP_UNUSED(_RP_extractFeatureVel); } } 





void extractFeatureVelRel( PyArrayContainer fv, const int N_row, const int off_begin, const BasicParticleSystem &p, const MACGrid &vel, const Real scale=1.0, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0, const int window=1, const Real h=1.0) {
	knExtractFeatureVelRel(p, reinterpret_cast<Real*>(fv.pData), N_row, off_begin, vel, scale, ptype, exclude, window, h);
} static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "extractFeatureVelRel" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; PyArrayContainer fv = _args.get<PyArrayContainer >("fv",0,&_lock); const int N_row = _args.get<int >("N_row",1,&_lock); const int off_begin = _args.get<int >("off_begin",2,&_lock); const BasicParticleSystem& p = *_args.getPtr<BasicParticleSystem >("p",3,&_lock); const MACGrid& vel = *_args.getPtr<MACGrid >("vel",4,&_lock); const Real scale = _args.getOpt<Real >("scale",5,1.0,&_lock); const ParticleDataImpl<int> * ptype = _args.getPtrOpt<ParticleDataImpl<int>  >("ptype",6,NULL,&_lock); const int exclude = _args.getOpt<int >("exclude",7,0,&_lock); const int window = _args.getOpt<int >("window",8,1,&_lock); const Real h = _args.getOpt<Real >("h",9,1.0,&_lock);   _retval = getPyNone(); extractFeatureVelRel(fv,N_row,off_begin,p,vel,scale,ptype,exclude,window,h);  _args.check(); } pbFinalizePlugin(parent,"extractFeatureVelRel", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("extractFeatureVelRel",e.what()); return 0; } } static const Pb::Register _RP_extractFeatureVelRel ("","extractFeatureVelRel",_W_2);  extern "C" { void PbRegister_extractFeatureVelRel() { KEEP_UNUSED(_RP_extractFeatureVelRel); } } 





void extractFeaturePhi( PyArrayContainer fv, const int N_row, const int off_begin, const BasicParticleSystem &p, const Grid<Real> &phi, const Real scale=1.0, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0, const int window=1, const Real h=1.0) {
	knExtractFeaturePhi(p, reinterpret_cast<Real*>(fv.pData), N_row, off_begin, phi, scale, ptype, exclude, window, h);
} static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "extractFeaturePhi" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; PyArrayContainer fv = _args.get<PyArrayContainer >("fv",0,&_lock); const int N_row = _args.get<int >("N_row",1,&_lock); const int off_begin = _args.get<int >("off_begin",2,&_lock); const BasicParticleSystem& p = *_args.getPtr<BasicParticleSystem >("p",3,&_lock); const Grid<Real> & phi = *_args.getPtr<Grid<Real>  >("phi",4,&_lock); const Real scale = _args.getOpt<Real >("scale",5,1.0,&_lock); const ParticleDataImpl<int> * ptype = _args.getPtrOpt<ParticleDataImpl<int>  >("ptype",6,NULL,&_lock); const int exclude = _args.getOpt<int >("exclude",7,0,&_lock); const int window = _args.getOpt<int >("window",8,1,&_lock); const Real h = _args.getOpt<Real >("h",9,1.0,&_lock);   _retval = getPyNone(); extractFeaturePhi(fv,N_row,off_begin,p,phi,scale,ptype,exclude,window,h);  _args.check(); } pbFinalizePlugin(parent,"extractFeaturePhi", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("extractFeaturePhi",e.what()); return 0; } } static const Pb::Register _RP_extractFeaturePhi ("","extractFeaturePhi",_W_3);  extern "C" { void PbRegister_extractFeaturePhi() { KEEP_UNUSED(_RP_extractFeaturePhi); } } 





void extractFeatureGeo( PyArrayContainer fv, const int N_row, const int off_begin, const BasicParticleSystem &p, const FlagGrid &flag, const Real scale=1.0, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0, const int window=1, const Real h=1.0) {
	knExtractFeatureGeo(p, reinterpret_cast<Real*>(fv.pData), N_row, off_begin, flag, scale, ptype, exclude, window, h);
} static PyObject* _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "extractFeatureGeo" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; PyArrayContainer fv = _args.get<PyArrayContainer >("fv",0,&_lock); const int N_row = _args.get<int >("N_row",1,&_lock); const int off_begin = _args.get<int >("off_begin",2,&_lock); const BasicParticleSystem& p = *_args.getPtr<BasicParticleSystem >("p",3,&_lock); const FlagGrid& flag = *_args.getPtr<FlagGrid >("flag",4,&_lock); const Real scale = _args.getOpt<Real >("scale",5,1.0,&_lock); const ParticleDataImpl<int> * ptype = _args.getPtrOpt<ParticleDataImpl<int>  >("ptype",6,NULL,&_lock); const int exclude = _args.getOpt<int >("exclude",7,0,&_lock); const int window = _args.getOpt<int >("window",8,1,&_lock); const Real h = _args.getOpt<Real >("h",9,1.0,&_lock);   _retval = getPyNone(); extractFeatureGeo(fv,N_row,off_begin,p,flag,scale,ptype,exclude,window,h);  _args.check(); } pbFinalizePlugin(parent,"extractFeatureGeo", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("extractFeatureGeo",e.what()); return 0; } } static const Pb::Register _RP_extractFeatureGeo ("","extractFeatureGeo",_W_4);  extern "C" { void PbRegister_extractFeatureGeo() { KEEP_UNUSED(_RP_extractFeatureGeo); } } 

// non-numpy related helpers

//! region detection functions


void extendRegion(FlagGrid &flags, const int region, const int exclude, const int depth) {
	flags.extendRegion(region, exclude, depth);
} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "extendRegion" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); const int region = _args.get<int >("region",1,&_lock); const int exclude = _args.get<int >("exclude",2,&_lock); const int depth = _args.get<int >("depth",3,&_lock);   _retval = getPyNone(); extendRegion(flags,region,exclude,depth);  _args.check(); } pbFinalizePlugin(parent,"extendRegion", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("extendRegion",e.what()); return 0; } } static const Pb::Register _RP_extendRegion ("","extendRegion",_W_5);  extern "C" { void PbRegister_extendRegion() { KEEP_UNUSED(_RP_extendRegion); } } 

// particle helpers

// TODO: merge with KnUpdateVelocity (particle.h); duplicated

 struct KnAddForcePvel : public KernelBase { KnAddForcePvel(ParticleDataImpl<Vec3> &v, const Vec3 &da, const ParticleDataImpl<int> *ptype, const int exclude) :  KernelBase(v.size()) ,v(v),da(da),ptype(ptype),exclude(exclude)   { runMessage(); run(); }   inline void op(IndexInt idx, ParticleDataImpl<Vec3> &v, const Vec3 &da, const ParticleDataImpl<int> *ptype, const int exclude )  {
	if(ptype && ((*ptype)[idx] & exclude)) return;
	v[idx] += da;
}    inline ParticleDataImpl<Vec3> & getArg0() { return v; } typedef ParticleDataImpl<Vec3>  type0;inline const Vec3& getArg1() { return da; } typedef Vec3 type1;inline const ParticleDataImpl<int> * getArg2() { return ptype; } typedef ParticleDataImpl<int>  type2;inline const int& getArg3() { return exclude; } typedef int type3; void runMessage() { debMsg("Executing kernel KnAddForcePvel ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,v,da,ptype,exclude);  }   } ParticleDataImpl<Vec3> & v; const Vec3& da; const ParticleDataImpl<int> * ptype; const int exclude;   };
#line 205 "plugin/tfplugins.cpp"


//! add force to vec3 particle data (ie, a velocity)
// TODO: merge with ParticleSystem::updateVelocity (particle.h); duplicated

void addForcePvel(ParticleDataImpl<Vec3> &vel, const Vec3 &a, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude) {
	KnAddForcePvel(vel, a*dt, ptype, exclude);
} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "addForcePvel" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; ParticleDataImpl<Vec3> & vel = *_args.getPtr<ParticleDataImpl<Vec3>  >("vel",0,&_lock); const Vec3& a = _args.get<Vec3 >("a",1,&_lock); const Real dt = _args.get<Real >("dt",2,&_lock); const ParticleDataImpl<int> * ptype = _args.getPtr<ParticleDataImpl<int>  >("ptype",3,&_lock); const int exclude = _args.get<int >("exclude",4,&_lock);   _retval = getPyNone(); addForcePvel(vel,a,dt,ptype,exclude);  _args.check(); } pbFinalizePlugin(parent,"addForcePvel", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("addForcePvel",e.what()); return 0; } } static const Pb::Register _RP_addForcePvel ("","addForcePvel",_W_6);  extern "C" { void PbRegister_addForcePvel() { KEEP_UNUSED(_RP_addForcePvel); } } 

//! retrieve velocity from position change

void updateVelocityFromDeltaPos(BasicParticleSystem& parts, ParticleDataImpl<Vec3> &vel, const ParticleDataImpl<Vec3> &x_prev, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude) {
	parts.updateVelocityFromDeltaPos(vel, x_prev, dt, ptype, exclude);
} static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "updateVelocityFromDeltaPos" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<Vec3> & vel = *_args.getPtr<ParticleDataImpl<Vec3>  >("vel",1,&_lock); const ParticleDataImpl<Vec3> & x_prev = *_args.getPtr<ParticleDataImpl<Vec3>  >("x_prev",2,&_lock); const Real dt = _args.get<Real >("dt",3,&_lock); const ParticleDataImpl<int> * ptype = _args.getPtr<ParticleDataImpl<int>  >("ptype",4,&_lock); const int exclude = _args.get<int >("exclude",5,&_lock);   _retval = getPyNone(); updateVelocityFromDeltaPos(parts,vel,x_prev,dt,ptype,exclude);  _args.check(); } pbFinalizePlugin(parent,"updateVelocityFromDeltaPos", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("updateVelocityFromDeltaPos",e.what()); return 0; } } static const Pb::Register _RP_updateVelocityFromDeltaPos ("","updateVelocityFromDeltaPos",_W_7);  extern "C" { void PbRegister_updateVelocityFromDeltaPos() { KEEP_UNUSED(_RP_updateVelocityFromDeltaPos); } } 

//! simple foward Euler integration for particle system

void eulerStep(BasicParticleSystem& parts, const ParticleDataImpl<Vec3> &vel, const ParticleDataImpl<int> *ptype, const int exclude) {
	parts.advect(vel, ptype, exclude);
} static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "eulerStep" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); const ParticleDataImpl<Vec3> & vel = *_args.getPtr<ParticleDataImpl<Vec3>  >("vel",1,&_lock); const ParticleDataImpl<int> * ptype = _args.getPtr<ParticleDataImpl<int>  >("ptype",2,&_lock); const int exclude = _args.get<int >("exclude",3,&_lock);   _retval = getPyNone(); eulerStep(parts,vel,ptype,exclude);  _args.check(); } pbFinalizePlugin(parent,"eulerStep", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("eulerStep",e.what()); return 0; } } static const Pb::Register _RP_eulerStep ("","eulerStep",_W_8);  extern "C" { void PbRegister_eulerStep() { KEEP_UNUSED(_RP_eulerStep); } } 

// TODO: merge with KnSetType (particle.h); duplicated

 struct KnSetPartType : public KernelBase { KnSetPartType(ParticleDataImpl<int> &ptype, BasicParticleSystem &part, const int mark, const int stype, const FlagGrid &flags, const int cflag) :  KernelBase(ptype.size()) ,ptype(ptype),part(part),mark(mark),stype(stype),flags(flags),cflag(cflag)   { runMessage(); run(); }   inline void op(IndexInt idx, ParticleDataImpl<int> &ptype, BasicParticleSystem &part, const int mark, const int stype, const FlagGrid &flags, const int cflag )  {
	if(flags.isInBounds(part.getPos(idx), 0) && (flags.getAt(part.getPos(idx))&cflag) && (ptype[idx]&stype)) ptype[idx] = mark;
}    inline ParticleDataImpl<int> & getArg0() { return ptype; } typedef ParticleDataImpl<int>  type0;inline BasicParticleSystem& getArg1() { return part; } typedef BasicParticleSystem type1;inline const int& getArg2() { return mark; } typedef int type2;inline const int& getArg3() { return stype; } typedef int type3;inline const FlagGrid& getArg4() { return flags; } typedef FlagGrid type4;inline const int& getArg5() { return cflag; } typedef int type5; void runMessage() { debMsg("Executing kernel KnSetPartType ", 3); debMsg("Kernel range" <<  " size "<<  size  << " "   , 4); }; void run() {   const IndexInt _sz = size; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (IndexInt i = 0; i < _sz; i++) op(i,ptype,part,mark,stype,flags,cflag);  }   } ParticleDataImpl<int> & ptype; BasicParticleSystem& part; const int mark; const int stype; const FlagGrid& flags; const int cflag;   };
#line 230 "plugin/tfplugins.cpp"


// TODO: merge with ParticleSystem::setType (particle.h); duplicated

void setPartType(BasicParticleSystem &parts, ParticleDataImpl<int> &ptype, const int mark, const int stype, const FlagGrid &flags, const int cflag) {
	KnSetPartType(ptype, parts, mark, stype, flags, cflag);
} static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "setPartType" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; BasicParticleSystem& parts = *_args.getPtr<BasicParticleSystem >("parts",0,&_lock); ParticleDataImpl<int> & ptype = *_args.getPtr<ParticleDataImpl<int>  >("ptype",1,&_lock); const int mark = _args.get<int >("mark",2,&_lock); const int stype = _args.get<int >("stype",3,&_lock); const FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",4,&_lock); const int cflag = _args.get<int >("cflag",5,&_lock);   _retval = getPyNone(); setPartType(parts,ptype,mark,stype,flags,cflag);  _args.check(); } pbFinalizePlugin(parent,"setPartType", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("setPartType",e.what()); return 0; } } static const Pb::Register _RP_setPartType ("","setPartType",_W_9);  extern "C" { void PbRegister_setPartType() { KEEP_UNUSED(_RP_setPartType); } } 

} //namespace


