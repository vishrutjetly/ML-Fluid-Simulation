




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep link).




#include "turbulencepart.h"
namespace Manta {
#ifdef _C_ParticleSystem
 static const Pb::Register _R_21 ("ParticleSystem<TurbulenceParticleData>","ParticleSystem<TurbulenceParticleData>","ParticleBase"); template<> const char* Namify<ParticleSystem<TurbulenceParticleData> >::S = "ParticleSystem<TurbulenceParticleData>"; 
 static const Pb::Register _R_22 ("ParticleSystem<TurbulenceParticleData>","ParticleSystem",ParticleSystem<TurbulenceParticleData>::_W_2); 
 static const Pb::Register _R_23 ("ParticleSystem<TurbulenceParticleData>","size",ParticleSystem<TurbulenceParticleData>::_W_3); 
 static const Pb::Register _R_24 ("ParticleSystem<TurbulenceParticleData>","pySize",ParticleSystem<TurbulenceParticleData>::_W_4); 
 static const Pb::Register _R_25 ("ParticleSystem<TurbulenceParticleData>","setPos",ParticleSystem<TurbulenceParticleData>::_W_5); 
 static const Pb::Register _R_26 ("ParticleSystem<TurbulenceParticleData>","getPos",ParticleSystem<TurbulenceParticleData>::_W_6); 
 static const Pb::Register _R_27 ("ParticleSystem<TurbulenceParticleData>","getPosPdata",ParticleSystem<TurbulenceParticleData>::_W_7); 
 static const Pb::Register _R_28 ("ParticleSystem<TurbulenceParticleData>","setPosPdata",ParticleSystem<TurbulenceParticleData>::_W_8); 
 static const Pb::Register _R_29 ("ParticleSystem<TurbulenceParticleData>","transformFrom",ParticleSystem<TurbulenceParticleData>::_W_9); 
 static const Pb::Register _R_30 ("ParticleSystem<TurbulenceParticleData>","clear",ParticleSystem<TurbulenceParticleData>::_W_10); 
 static const Pb::Register _R_31 ("ParticleSystem<TurbulenceParticleData>","killOutOfDomain",ParticleSystem<TurbulenceParticleData>::_W_11); 
 static const Pb::Register _R_32 ("ParticleSystem<TurbulenceParticleData>","killRegion",ParticleSystem<TurbulenceParticleData>::_W_12); 
 static const Pb::Register _R_33 ("ParticleSystem<TurbulenceParticleData>","killSelected",ParticleSystem<TurbulenceParticleData>::_W_13); 
 static const Pb::Register _R_34 ("ParticleSystem<TurbulenceParticleData>","sampleFrom",ParticleSystem<TurbulenceParticleData>::_W_14); 
 static const Pb::Register _R_35 ("ParticleSystem<TurbulenceParticleData>","setType",ParticleSystem<TurbulenceParticleData>::_W_15); 
 static const Pb::Register _R_36 ("ParticleSystem<TurbulenceParticleData>","advectInGrid",ParticleSystem<TurbulenceParticleData>::_W_16); 
 static const Pb::Register _R_37 ("ParticleSystem<TurbulenceParticleData>","advect",ParticleSystem<TurbulenceParticleData>::_W_17); 
 static const Pb::Register _R_38 ("ParticleSystem<TurbulenceParticleData>","updateVelocity",ParticleSystem<TurbulenceParticleData>::_W_18); 
 static const Pb::Register _R_39 ("ParticleSystem<TurbulenceParticleData>","updateVelocityFromDeltaPos",ParticleSystem<TurbulenceParticleData>::_W_19); 
 static const Pb::Register _R_40 ("ParticleSystem<TurbulenceParticleData>","getAtFromFlagGrid",ParticleSystem<TurbulenceParticleData>::_W_20); 
 static const Pb::Register _R_41 ("ParticleSystem<TurbulenceParticleData>","getInterpolatedFromMACGrid",ParticleSystem<TurbulenceParticleData>::_W_21); 
 static const Pb::Register _R_42 ("ParticleSystem<TurbulenceParticleData>","getInterpolatedFromGrid",ParticleSystem<TurbulenceParticleData>::_W_22); 
 static const Pb::Register _R_43 ("ParticleSystem<TurbulenceParticleData>","getInterpolatedFromVecGrid",ParticleSystem<TurbulenceParticleData>::_W_23); 
 static const Pb::Register _R_44 ("ParticleSystem<TurbulenceParticleData>","projectOutside",ParticleSystem<TurbulenceParticleData>::_W_24); 
 static const Pb::Register _R_45 ("ParticleSystem<TurbulenceParticleData>","projectOutOfBnd",ParticleSystem<TurbulenceParticleData>::_W_25); 
#endif
#ifdef _C_TurbulenceParticleSystem
 static const Pb::Register _R_46 ("TurbulenceParticleSystem","TurbulenceParticleSystem","ParticleSystem<TurbulenceParticleData>"); template<> const char* Namify<TurbulenceParticleSystem >::S = "TurbulenceParticleSystem"; 
 static const Pb::Register _R_47 ("TurbulenceParticleSystem","TurbulenceParticleSystem",TurbulenceParticleSystem::_W_0); 
 static const Pb::Register _R_48 ("TurbulenceParticleSystem","resetTexCoords",TurbulenceParticleSystem::_W_1); 
 static const Pb::Register _R_49 ("TurbulenceParticleSystem","seed",TurbulenceParticleSystem::_W_2); 
 static const Pb::Register _R_50 ("TurbulenceParticleSystem","synthesize",TurbulenceParticleSystem::_W_3); 
 static const Pb::Register _R_51 ("TurbulenceParticleSystem","deleteInObstacle",TurbulenceParticleSystem::_W_4); 
#endif
extern "C" {
void PbRegister_file_21()
{
	KEEP_UNUSED(_R_21);
	KEEP_UNUSED(_R_22);
	KEEP_UNUSED(_R_23);
	KEEP_UNUSED(_R_24);
	KEEP_UNUSED(_R_25);
	KEEP_UNUSED(_R_26);
	KEEP_UNUSED(_R_27);
	KEEP_UNUSED(_R_28);
	KEEP_UNUSED(_R_29);
	KEEP_UNUSED(_R_30);
	KEEP_UNUSED(_R_31);
	KEEP_UNUSED(_R_32);
	KEEP_UNUSED(_R_33);
	KEEP_UNUSED(_R_34);
	KEEP_UNUSED(_R_35);
	KEEP_UNUSED(_R_36);
	KEEP_UNUSED(_R_37);
	KEEP_UNUSED(_R_38);
	KEEP_UNUSED(_R_39);
	KEEP_UNUSED(_R_40);
	KEEP_UNUSED(_R_41);
	KEEP_UNUSED(_R_42);
	KEEP_UNUSED(_R_43);
	KEEP_UNUSED(_R_44);
	KEEP_UNUSED(_R_45);
	KEEP_UNUSED(_R_46);
	KEEP_UNUSED(_R_47);
	KEEP_UNUSED(_R_48);
	KEEP_UNUSED(_R_49);
	KEEP_UNUSED(_R_50);
	KEEP_UNUSED(_R_51);
}
}}