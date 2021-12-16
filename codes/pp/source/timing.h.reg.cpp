




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep link).




#include "timing.h"
namespace Manta {
#ifdef _C_Timer
 static const Pb::Register _R_16 ("Timer","Timer","PbClass"); template<> const char* Namify<Timer >::S = "Timer"; 
 static const Pb::Register _R_17 ("Timer","Timer",Timer::_W_4); 
 static const Pb::Register _R_18 ("Timer","push",Timer::_W_5); 
 static const Pb::Register _R_19 ("Timer","pop",Timer::_W_6); 
 static const Pb::Register _R_20 ("Timer","popName",Timer::_W_7); 
 static const Pb::Register _R_21 ("Timer","popTiming",Timer::_W_8); 
 static const Pb::Register _R_22 ("Timer","popTimingStr",Timer::_W_9); 
#endif
#ifdef _C_Timings
 static const Pb::Register _R_23 ("Timings","Timings","PbClass"); template<> const char* Namify<Timings >::S = "Timings"; 
 static const Pb::Register _R_24 ("Timings","Timings",Timings::_W_0); 
 static const Pb::Register _R_25 ("Timings","display",Timings::_W_1); 
 static const Pb::Register _R_26 ("Timings","step",Timings::_W_2); 
 static const Pb::Register _R_27 ("Timings","saveMean",Timings::_W_3); 
#endif
extern "C" {
void PbRegister_file_16()
{
	KEEP_UNUSED(_R_16);
	KEEP_UNUSED(_R_17);
	KEEP_UNUSED(_R_18);
	KEEP_UNUSED(_R_19);
	KEEP_UNUSED(_R_20);
	KEEP_UNUSED(_R_21);
	KEEP_UNUSED(_R_22);
	KEEP_UNUSED(_R_23);
	KEEP_UNUSED(_R_24);
	KEEP_UNUSED(_R_25);
	KEEP_UNUSED(_R_26);
	KEEP_UNUSED(_R_27);
}
}}