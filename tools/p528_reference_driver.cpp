#include "include/p528.h"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){
  double d,h1,h2,f,p; int pol;
  while(scanf("%lf %lf %lf %lf %d %lf",&d,&h1,&h2,&f,&pol,&p)==6){
    Result r; int rtn=P528(d,h1,h2,f,pol,p,&r);
    printf("%.17g %.17g %.17g %.17g %d %.17g %d %.17g %.17g %.17g %.17g %d %d\n",d,h1,h2,f,pol,p,rtn,r.d__km,r.A__db,r.A_fs__db,r.A_a__db,r.propagation_mode,r.warnings);
  }
}
