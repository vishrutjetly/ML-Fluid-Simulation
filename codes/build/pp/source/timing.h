




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/vishrut/Study/USC/Sem-1/CSCI 596: Scientific Computing and Visualization/project/mlflip/source/timing.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Plugin timing
 *
 ******************************************************************************/

#ifndef _TIMING_H
#define _TIMING_H

#include "manta.h"
#include <map>
namespace Manta {


class TimingData {
private:
	TimingData();
public:
	static TimingData& instance() { static TimingData a; return a; }

	void print();
	void saveMean(const std::string& filename);
	void start(FluidSolver* parent, const std::string& name);
	void stop(FluidSolver* parent, const std::string& name);
	void step();

protected:
	struct TimingSet {
		TimingSet() : num(0), updated(false) { cur.clear(); total.clear(); }
		MuTime cur, total;
		int num;
		bool updated;
		std::string solver;
	};
	bool updated;

	int num;
	MuTime mPluginTimer;
	std::string mLastPlugin;
	std::map< std::string, std::vector<TimingSet> > mData;
};

// Python interface
class Timings : public PbClass {public:
	Timings() :PbClass(NULL){} static int _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "Timings::Timings" , !noTiming ); { ArgLocker _lock;  obj = new Timings(); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Timings::Timings" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("Timings::Timings",e.what()); return -1; } }

	void display() { TimingData::instance().print(); } static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timings* pbo = dynamic_cast<Timings*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timings::display" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->display();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timings::display" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timings::display",e.what()); return 0; } }
	void step() { TimingData::instance().step(); } static PyObject* _W_2 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timings* pbo = dynamic_cast<Timings*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timings::step" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->step();  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timings::step" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timings::step",e.what()); return 0; } } 	void saveMean(const std::string& file) { TimingData::instance().saveMean(file); } static PyObject* _W_3 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timings* pbo = dynamic_cast<Timings*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timings::saveMean" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const std::string& file = _args.get<std::string >("file",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->saveMean(file);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timings::saveMean" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timings::saveMean",e.what()); return 0; } } public: PbArgs _args; }
#define _C_Timings
;

class Timer : public PbClass {public:
	Timer() :PbClass(NULL){} static int _W_4 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = Pb::objFromPy(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(0, "Timer::Timer" , !noTiming ); { ArgLocker _lock;  obj = new Timer(); obj->registerObject(_self, &_args); _args.check(); } pbFinalizePlugin(obj->getParent(),"Timer::Timer" , !noTiming ); return 0; } catch(std::exception& e) { pbSetError("Timer::Timer",e.what()); return -1; } }

	void push(const std::string &name) {
		mName.push_back(name);
		mTime.push_back(MuTime());
	} static PyObject* _W_5 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timer* pbo = dynamic_cast<Timer*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timer::push" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const std::string& name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->push(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timer::push" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timer::push",e.what()); return 0; } }
	void pop(const std::string &name) {
		mPopTiming = MuTime() - mTime.back();
		mPopName = mName.back();
		assertMsg(name==mPopName, "incorrect use of timer; you must pop " + mPopName);
		mName.pop_back();
		mTime.pop_back();
	} static PyObject* _W_6 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timer* pbo = dynamic_cast<Timer*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timer::pop" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock; const std::string& name = _args.get<std::string >("name",0,&_lock);  pbo->_args.copy(_args);  _retval = getPyNone(); pbo->pop(name);  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timer::pop" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timer::pop",e.what()); return 0; } }
	const std::string& popName() const { return mPopName; } static PyObject* _W_7 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timer* pbo = dynamic_cast<Timer*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timer::popName" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->popName());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timer::popName" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timer::popName",e.what()); return 0; } }
	int popTiming() const { return mPopTiming.time; } static PyObject* _W_8 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timer* pbo = dynamic_cast<Timer*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timer::popTiming" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->popTiming());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timer::popTiming" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timer::popTiming",e.what()); return 0; } }
	std::string popTimingStr() const { return mPopTiming.toString(); } static PyObject* _W_9 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); Timer* pbo = dynamic_cast<Timer*>(Pb::objFromPy(_self)); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(pbo->getParent(), "Timer::popTimingStr" , !noTiming); PyObject *_retval = 0; { ArgLocker _lock;  pbo->_args.copy(_args);  _retval = toPy(pbo->popTimingStr());  pbo->_args.check(); } pbFinalizePlugin(pbo->getParent(),"Timer::popTimingStr" , !noTiming); return _retval; } catch(std::exception& e) { pbSetError("Timer::popTimingStr",e.what()); return 0; } }

private:
	std::vector<std::string> mName;
	std::vector<MuTime> mTime;

	MuTime mPopTiming; 	std::string mPopName; public: PbArgs _args; }
#define _C_Timer
;

}

#endif


