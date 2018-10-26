/*######################  RANDOM NUMBER GENERATOR  #########################*/
#include <stdio.h>
#include "fluke.h"

#define NO 0
#define YES 1

int initialized = NO, i97, j97;
double urand[98], crand, cdrand, cmrand;

extern long seed1, seed2;

double Fluke(void)
{               /*   FlukeInit() must be called to initialize Fluke()       */
  double uni;

  if (initialized == NO) {
    FlukeInit();
    initialized=YES;
  }
  
  uni=urand[i97]-urand[j97];  if (uni < 0.0) uni=uni+1.0;  urand[i97]=uni;
  i97--;  if (i97==0) i97=97;  j97--;  if (j97==0) j97=97;
  crand=crand-cdrand;  if (crand<0.0) crand=crand+cmrand;
  uni=uni-crand;  if (uni<0.0) uni=uni+1.0;
  return (uni);
}
/*______________________ test random number generator ______________________

void TestFluke(void)
{
   int ii, jj;
   char xx;
   
   seed1 = 1802;
   seed2 = 9373;
   
   FlukeInit();
   for (ii = 0; ii < 50000000; ++ii) {
	 xx = (int)(256*Fluke());
	 //xx = rand()%256;
	 printf ("%c",xx);
   }
}
_________________ random number initialization ___________________________*/  

void FlukeInit(void)
{
  int ii, jj;  long int i,j,k,l,m;  float s, t;

  if ((seed1<=0)||(seed1>31328)||(seed2<=0)||(seed2>30081))
  {
    printf("\nThe first random number seed must be between 0 and 31328");
    printf("\nThe second seed must be between 0 and 30081\n");
  }
  i=((seed1/177)%177)+2; j=(seed1%177)+2; k=((seed2/169)%178)+1; l=seed2%169;
  for (ii=1;ii<=97;ii++)  { s=0.0;  t=0.5;
    for (jj=1;jj<=24;jj++)  {m=(((i*j)%179)*k)%179; i=j; j=k; k=m;
      l=(53*l+1)%169;  if (((l*m)%64)>=32) s=s+t;   t=0.5*t;       }
    urand[ii]=s;  }
  crand=362436.0/16777216.0; cdrand=7654321.0/16777216.0;
  cmrand=16777213.0/16777216.0;  i97=97; j97=33;
  initialized = YES;
}
