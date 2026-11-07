#include "include/p528.h"
#include "include/p676.h"
#include "include/p835.h"
#include <cstdio>
int main(){ setvbuf(stdout,NULL,_IONBF,0);
  double hs[]={0,1,5,11,15,20,30,47,60,80};
  for(double h: hs){ double T=GlobalTemperature(h),P=GlobalPressure(h),e=GlobalWetPressure(h);
    printf("atm %.17g %.17g %.17g %.17g\n",h,T,P,e); }
  double fs[]={0.1,1,8.16,22.235,30};
  for(double f: fs){ double T=GlobalTemperature(0),P=GlobalPressure(0),e=GlobalWetPressure(0);
    printf("gam %.17g %.17g\n",f,SpecificAttenuation(f,T,e,P)); }
  struct {double f,h1,h2,b;} sp[]={{8.16,0,0.12192,1.5707963267948966},{8.16,0.01,0.12192,1.4},{1,0,10,1.5707963267948966},{30,0,20,0.5},{8.16,0.12192,0.5,1.5737963267948966}};
  for(auto s: sp){ SlantPathAttenuationResult r; SlantPathAttenuation(s.f,s.h1,s.h2,s.b,&r);
    printf("slt %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",s.f,s.h1,s.h2,s.b,r.A_gas__db,r.bending__rad,r.a__km,r.angle__rad,r.delta_L__km); }
}
