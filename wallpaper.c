#include "globals.h"
#include "groups.h"
#include "shapes.h"
#include "fluke.h"

extern double shape_var;

char directory[200];
FILE *fpout;
FILE *fp_isopointal_out;

int numoccsites;
/*int RANDSEED=778;  this may get changed by an argument */  
double max_step_size=0.01;
long int seed1=1802, seed2=9373;

void plot_shape (struct shapetype *s, char filename[MAXFILENAME]) {
  int i;
  double phi;
  
  FILE *fp;
  fp = fopen(filename, "w");
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    fprintf(fp, "%.12f %.12f\n", s->r[i]*cos(phi), s->r[i]*sin(phi));
  } 
  fclose(fp);
}

void print_group_definitions() {
  int wg, wk, mp;

  for (wg=0; wg<NUMWALLPAPERGROUPS; wg++) {
    printf("wallpaper group %d:\n", wg);
    for (wk=0; wk<wallpapergroups[wg].num_wyckoffs; wk++) { 
      printf("wyckoff site: %c\n", wallpapergroups[wg].wyckoffs[wk].letter);
      for (mp=0; mp<wallpapergroups[wg].wyckoffs[wk].multiplicity;mp++) {
	printf("x_new = %1.f*x_old + %1.f*y_old + %4.2f\n", 
	       wallpapergroups[wg].wyckoffs[wk].image[mp].coord_coeffs[0][0],
	       wallpapergroups[wg].wyckoffs[wk].image[mp].coord_coeffs[0][1],
	       wallpapergroups[wg].wyckoffs[wk].image[mp].coord_coeffs[0][2]
	       );
	printf("y_new = %1.f*x_old + %1.f*y_old + %4.2f\n",
	       wallpapergroups[wg].wyckoffs[wk].image[mp].coord_coeffs[1][0],
	       wallpapergroups[wg].wyckoffs[wk].image[mp].coord_coeffs[1][1],
	       wallpapergroups[wg].wyckoffs[wk].image[mp].coord_coeffs[1][2]
	       );
      }
      printf("\n");
    }
    printf("-------------------------------------------\n");
  }
}

int factorial(int n) {
  if((n==0)||(n==1))
    return(1);
  else
    return(n*factorial(n-1));
}

int power( int b , int p) {
  if( p > 1 )
    return b * power( b , p-1 );
  else
    return b;
}

int enumerate(int OS[][MAXNUMOCCSITES], int maxcombinations, int n_comb, int os, int bag_i, int combination_bag[], int bag_size) {
  int i, j;
  int count=0, new;

  if (os == numoccsites)
    return 1;
  else {
    for (i = bag_i; i < bag_size-(numoccsites-os-1); i++) {
      if ((i!=bag_i)&&(combination_bag[i]==combination_bag[i-1])) continue;
      new = enumerate(OS, maxcombinations, n_comb+count, os+1, i+1, combination_bag, bag_size);
      for (j=0; j<new; j++) { 
	OS[n_comb+j+count][os]=combination_bag[i];
	//printf("saving [%d][%d]=%d (new %d, count %d)\n", n_comb+j+count, os, i, new, count);
      }
      count+=new;
    }
  }
  return count;
}

double area(struct shapetype *shape) {
  /* a calculation of the area of the polygon */
  int i;
  double areasum=0.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    // old   areasum += 0.5*pow(0.5*(shape->r[i]+shape->r[(i+1)%SHAPE_RESOLUTION]),2.0)*(2.0*M_PI/SHAPE_RESOLUTION);
	// the following is the exact triangle area
	areasum += 0.5*(shape->r[i])*(shape->r[(i+1)%SHAPE_RESOLUTION])*sin(2.0*M_PI/SHAPE_RESOLUTION);
  }
  //printf("area %f\n", areasum);
  return areasum;
}

void fractcoords(double coord_coeffs[DIM][DIM+1], double *sitevariables[3], double newfcoords[DIM]) {  
  /* converts site variables and wyckoff site coefficients into the location of the actual wyckoff image in fractional coordinates */
  newfcoords[0] = coord_coeffs[0][0]*(*sitevariables[0])+ coord_coeffs[0][1]*(*sitevariables[1])+ coord_coeffs[0][2];
  while (newfcoords[0]<0.0) newfcoords[0] += 1.0;
  while (newfcoords[0]>1.0) newfcoords[0] -= 1.0;
  newfcoords[1] = coord_coeffs[1][0]*(*sitevariables[0])+ coord_coeffs[1][1]*(*sitevariables[1])+ coord_coeffs[1][2];
  while (newfcoords[1]<0.0) newfcoords[1] += 1.0;
  while (newfcoords[1]>1.0) newfcoords[1] -= 1.0;
  return;
}

void realcoords(double fcoords[DIM], double coords[DIM], double *cellsides[DIM], double *cellangles[DIM-1]) {
  //printf("cell %f %f %f\n", (*cellsides[0]), (*cellsides[1]), (*cellangles[0]));
  coords[0] = fcoords[0]*(*cellsides[0])+fcoords[1]*(*cellsides[1])*cos(*cellangles[0]);
  coords[1] = fcoords[1]*(*cellsides[1])*sin(*cellangles[0]);
  return;
}

bool segments_cross(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy) {
  //  Determines the intersection point of the line segment defined by points A and B
  //  with the line segment defined by points C and D.
  //
  //  Returns YES if the intersection point was found, and stores that point in X,Y.
  //  Returns NO if there is no determinable intersection point, in which case X,Y will
  //  be unmodified.
  
  double  distAB, theCos, theSin, newX, ABpos ;
  //double X, Y;
  
  //  Fail if either line segment is zero-length.
  //if (Ax==Bx && Ay==By || Cx==Dx && Cy==Dy) return NO;

  //  Count as a clash if the segments share an end-point.
  if (((Ax==Cx) && (Ay==Cy)) || ((Bx==Cx) && (By==Cy))
      ||  ((Ax==Dx) && (Ay==Dy)) || ((Bx==Dx) && (By==Dy))) {
    return 1; }

  //  (1) Translate the system so that point A is on the origin.
  Bx-=Ax; By-=Ay;
  Cx-=Ax; Cy-=Ay;
  Dx-=Ax; Dy-=Ay;

  //  Discover the length of segment A-B.
  distAB=sqrt(Bx*Bx+By*By);

  //  (2) Rotate the system so that point B is on the positive X axis.
  theCos=Bx/distAB;
  theSin=By/distAB;
  newX=Cx*theCos+Cy*theSin;
  Cy  =Cy*theCos-Cx*theSin; Cx=newX;
  newX=Dx*theCos+Dy*theSin;
  Dy  =Dy*theCos-Dx*theSin; Dx=newX;

  //  Fail if segment C-D doesn't cross line A-B.
  if (((Cy<0.) && (Dy<0.)) || ((Cy>=0.) && (Dy>=0.))) return 0;

  //  (3) Discover the position of the intersection point along line A-B.
  ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);

  //  Fail if segment C-D crosses line A-B outside of segment A-B.
  if (ABpos<0. || ABpos>distAB) return 0;

  //  (4) Apply the discovered position to line A-B in the original coordinate system.
  //X=Ax+ABpos*theCos;
  //Y=Ay+ABpos*theSin;

  //  Success (i.e. clash).
  return 1; 
}

void dump_clash(struct shapetype *shape_a, struct shapetype *shape_b, double coords_a[DIM], double coords_b[DIM], double angle_a, double angle_b, bool flip_a, bool flip_b,  double central_distsq){
  
  double theta_ia0=0.0, theta_ia1;
  double theta_tab[SHAPE_RESOLUTION/2];
  double xb_tab[SHAPE_RESOLUTION/2];
  double yb_tab[SHAPE_RESOLUTION/2];
  double xa0, xa1, ya0, ya1;
  double centraldist;
  int ia, ib;
 
  int angle_a_int, angle_b_int;
  
  int scan_sign=1;
  
  FILE *fp_test;

  centraldist=sqrt(central_distsq);

  angle_a_int = (int)(angle_a*SHAPE_RESOLUTION/(2.0*M_PI)+0.5);
  angle_b_int = (int)(angle_b*SHAPE_RESOLUTION/(2.0*M_PI)+0.5);
  
  fp_test = fopen("testclash.xy","w");

  for (scan_sign=1; scan_sign>=-1; scan_sign -= 2) {
    for (ib = -1; ib<(SHAPE_RESOLUTION/4)+1; ib++) {
      theta_tab[ib+1] = fabs((angle_b_int-(1-2*flip_b)*scan_sign*ib)*2.0*M_PI/SHAPE_RESOLUTION-angle_b);
      //sine_theta_ib_tab[ib] = sin(theta_tab[ib]);
      xb_tab[ib+1] = shape_b->r[(angle_b_int-(1-2*flip_b)*scan_sign*ib+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_tab[ib+1]);
      yb_tab[ib+1] = shape_b->r[(angle_b_int-(1-2*flip_b)*scan_sign*ib+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_tab[ib+1]);
      if (ib==-1) yb_tab[ib+1] *= -1.0; /* to correct an issue with the fabs of the very first point */ 
      
      fprintf(fp_test,"%f\t%f\n",  xb_tab[ib+1], yb_tab[ib+1]);

    }

    fprintf(fp_test,"%f\t%f\n",  0.0, 0.0);


    theta_ia1 = fabs((angle_a_int+(1-2*flip_a)*scan_sign*-1)*2.0*M_PI/SHAPE_RESOLUTION-angle_a);
    xa1 = sqrt(central_distsq)-shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*-1+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_ia0);
    ya1 = shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*-1+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_ia0);
    ya1 *= -1.0; /* to correct an issue with the fabs of the very first point */ 
    for (ia = -1; ia<(SHAPE_RESOLUTION/4)+1; ia++) {

      theta_ia0 = theta_ia1;
      theta_ia1 = fabs((angle_a_int+(1-2*flip_a)*scan_sign*(ia+1))*2.0*M_PI/SHAPE_RESOLUTION-angle_a);
      xa0 = xa1;
      ya0 = ya1;
      xa1 = centraldist-shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*(ia+1)+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_ia1);
      ya1 = shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*(ia+1)+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_ia1);

      fprintf(fp_test,"%f\t%f\n",  xa0, ya0);
    }
  }
  fclose(fp_test);
}

bool indirect_clash(struct shapetype *shape_a, struct shapetype *shape_b, double coords_a[DIM], double coords_b[DIM], double angle_a, double angle_b, bool flip_a, bool flip_b,  double central_distsq, int tests) {
  
  double theta_ia0=0.0, theta_ia1;
  double theta_tab[SHAPE_RESOLUTION/2];
  double xb_tab[SHAPE_RESOLUTION/2];
  double yb_tab[SHAPE_RESOLUTION/2];
  double xa0, xa1, ya0, ya1;
  double centraldist;
  int ia, ib;
 
  int angle_a_int, angle_b_int;
  
  int scan_sign=1;
  
  /* this is the most expensive routine of all - try to optimize (especially inside the double loops) */
  
  /* TSH THIS DOES NOT UTILIZE THE FLIP_A and FLIP_B parameters 
     .... it sort of does now, not sure if it's dealing with the round offs correctly */
  
  centraldist=sqrt(central_distsq);

  angle_a_int = (int)(angle_a*SHAPE_RESOLUTION/(2.0*M_PI)+0.5);
  angle_b_int = (int)(angle_b*SHAPE_RESOLUTION/(2.0*M_PI)+0.5);
    
  for (scan_sign=1; scan_sign>=-1; scan_sign -= 2) {

    if (tests==0) {
      //If we have just initialised (before doing any moves), and are checking for clashes, we want to check all the way around the edge of each shape, to make sure that by chance the two shapes weren't positioned on top of each other.
      for (ib = -1; ib<SHAPE_RESOLUTION/2; ib++) {
	theta_tab[ib+1] = fabs((angle_b_int-(1-2*flip_b)*scan_sign*ib)*2.0*M_PI/SHAPE_RESOLUTION-angle_b);
	//sine_theta_ib_tab[ib] = sin(theta_tab[ib]);
	xb_tab[ib+1] = shape_b->r[(angle_b_int-(1-2*flip_b)*scan_sign*ib+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_tab[ib+1]);
	yb_tab[ib+1] = shape_b->r[(angle_b_int-(1-2*flip_b)*scan_sign*ib+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_tab[ib+1]);
	if (ib==-1) yb_tab[ib+1] *= -1.0; /* to correct an issue with the fabs of the very first point */ 
      }
      
      theta_ia1 = fabs((angle_a_int+(1-2*flip_a)*scan_sign*-1)*2.0*M_PI/SHAPE_RESOLUTION-angle_a);
      xa1 = sqrt(central_distsq)-shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*-1+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_ia0);
      ya1 = shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*-1+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_ia0);
      ya1 *= -1.0; /* to correct an issue with the fabs of the very first point */ 
      for (ia = -1; ia<SHAPE_RESOLUTION/2; ia++) {
	theta_ia0 = theta_ia1;
	theta_ia1 = fabs((angle_a_int+(1-2*flip_a)*scan_sign*(ia+1))*2.0*M_PI/SHAPE_RESOLUTION-angle_a);
	xa0 = xa1;
	ya0 = ya1;
	xa1 = centraldist-shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*(ia+1)+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_ia1);
	ya1 = shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*(ia+1)+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_ia1);
	
	for (ib = -1; ib<SHAPE_RESOLUTION/2; ib++) {	  
	  if (segments_cross(xa0, ya0, xa1, ya1, xb_tab[ib+1], yb_tab[ib+1], xb_tab[ib+2], yb_tab[ib+2])) {
	    //printf("segments cross 1: ib %d\n",ib);
	    return 1;
	  }
	}
      }
    } else {
      for (ib = -1; ib<(SHAPE_RESOLUTION/4)+1; ib++) {
	theta_tab[ib+1] = fabs((angle_b_int-(1-2*flip_b)*scan_sign*ib)*2.0*M_PI/SHAPE_RESOLUTION-angle_b);
	//sine_theta_ib_tab[ib] = sin(theta_tab[ib]);
	xb_tab[ib+1] = shape_b->r[(angle_b_int-(1-2*flip_b)*scan_sign*ib+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_tab[ib+1]);
	yb_tab[ib+1] = shape_b->r[(angle_b_int-(1-2*flip_b)*scan_sign*ib+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_tab[ib+1]);
	if (ib==-1) yb_tab[ib+1] *= -1.0; /* to correct an issue with the fabs of the very first point */ 
      }
      
      theta_ia1 = fabs((angle_a_int+(1-2*flip_a)*scan_sign*-1)*2.0*M_PI/SHAPE_RESOLUTION-angle_a);
      xa1 = sqrt(central_distsq)-shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*-1+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_ia0);
      ya1 = shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*-1+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_ia0);
      ya1 *= -1.0; /* to correct an issue with the fabs of the very first point */ 
      for (ia = -1; ia<(SHAPE_RESOLUTION/4)+1; ia++) {
	theta_ia0 = theta_ia1;
	theta_ia1 = fabs((angle_a_int+(1-2*flip_a)*scan_sign*(ia+1))*2.0*M_PI/SHAPE_RESOLUTION-angle_a);
	xa0 = xa1;
	ya0 = ya1;
	xa1 = centraldist-shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*(ia+1)+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*cos(theta_ia1);
	ya1 = shape_a->r[(angle_a_int+(1-2*flip_a)*scan_sign*(ia+1)+5*SHAPE_RESOLUTION)%SHAPE_RESOLUTION]*sin(theta_ia1);
	
	for (ib = -1; ib<(SHAPE_RESOLUTION/4)+1; ib++) {
	  if (segments_cross(xa0, ya0, xa1, ya1, xb_tab[ib+1], yb_tab[ib+1], xb_tab[ib+2], yb_tab[ib+2])) {
	    //printf("segments cross 2: ib %d %f %f %f %f %f %f %f %f\n",ib, xa0, ya0, xa1, ya1, xb_tab[ib+1], yb_tab[ib+1], xb_tab[ib+2], yb_tab[ib+2]);

	    //dump_clash(shape_a, shape_b, coords_a, coords_b, angle_a, angle_b, flip_a, flip_b, central_distsq);
	    return 1;
	  }
	}
      }
    }
  }

  return 0;
}


bool pair_clash(struct shapetype *shape_a, struct shapetype *shape_b, double coords_a[DIM], double coords_b[DIM], double orient_a, double orient_b, int rotation_a, int rotation_b, bool flip_a, bool flip_b, int tests) {
  double central_distsq, maxshapecontact_distsq;
  double a_to_b_incline, b_to_a_incline;
  
  double mathtol = 1e-8;
 
  central_distsq = (coords_a[0]-coords_b[0])*(coords_a[0]-coords_b[0])+(coords_a[1]-coords_b[1])*(coords_a[1]-coords_b[1]);
  maxshapecontact_distsq = shape_a->maxr+shape_b->maxr;
  maxshapecontact_distsq *= maxshapecontact_distsq;
  if (central_distsq>maxshapecontact_distsq) {
    /* if they're further apart than even the maximum shape radii measures, no clash*/
    return (0);
  }
   
  /* we need a polygon overlap calculation */
  /* first calculate the direct overlap */
  a_to_b_incline = acos((coords_b[0]-coords_a[0])/sqrt(central_distsq));
  
  if (isnan(a_to_b_incline)) {
    //printf("nan atob incline 1: %f %f %f =>>> acos(%f)\n", coords_b[0], coords_a[0], central_distsq,(coords_b[0]-coords_a[0])/sqrt(central_distsq));
    if (fabs((coords_b[0]-coords_a[0])/sqrt(central_distsq) - 1.0) < mathtol) {
      a_to_b_incline = 0.0;
    } else if (fabs((coords_b[0]-coords_a[0])/sqrt(central_distsq) + 1.0) < mathtol) {
      a_to_b_incline = M_PI;
    }
  } 
  if (coords_b[1]<coords_a[1]) {
    a_to_b_incline = 2.0*M_PI-a_to_b_incline;
  }
  b_to_a_incline = a_to_b_incline+M_PI;
  if (flip_a) {
    a_to_b_incline = 2*M_PI - a_to_b_incline;
  }
  if (flip_b) {
    b_to_a_incline = 2*M_PI - b_to_a_incline;
  }

  /* now add in the rotation due to the orientation parameters */
  /* This is unaffected by the flip states */
  a_to_b_incline += orient_a;
  b_to_a_incline += orient_b;
  
  /* now add in the rotation due to the rotation of this image wrt the other images of the same wyckoff */
  a_to_b_incline += (1-2*flip_a)*rotation_a*2.0*M_PI/SHAPE_RESOLUTION;
  b_to_a_incline += (1-2*flip_b)*rotation_b*2.0*M_PI/SHAPE_RESOLUTION;

  while (a_to_b_incline<0.0) a_to_b_incline += 2.0*M_PI;
  while (b_to_a_incline<0.0) b_to_a_incline += 2.0*M_PI;
  while (a_to_b_incline >= 2.0*M_PI) a_to_b_incline -= 2.0*M_PI;
  while (b_to_a_incline >= 2.0*M_PI) b_to_a_incline -= 2.0*M_PI;

  return indirect_clash(shape_a, shape_b, coords_a, coords_b, a_to_b_incline, b_to_a_incline, flip_a, flip_b, central_distsq, tests);
}

bool clash_polygon(struct shapetype *shape_a, struct shapetype *shape_b, struct wallpaper_group_type group, double *cellsides[DIM], double *cellangles[DIM-1],  int occupiedsites[MAXNUMOCCSITES], double *sitevariables[MAXNUMOCCSITES][3], int rep_a[2], int rep_b[2], bool flips[numoccsites], int tests) {
  double fcoords_a[DIM], fcoords_b[DIM];
  double img_fcoords_b[DIM];
  double coords_a[DIM], coords_b[DIM];
  //double distsq; 
  int cellimgx, cellimgy;
  double orient_a, orient_b;
  int rotation_a, rotation_b;
  bool flip_a, flip_b;

  fractcoords(group.wyckoffs[occupiedsites[rep_a[0]]].image[rep_a[1]].coord_coeffs, sitevariables[rep_a[0]], fcoords_a);
  fractcoords(group.wyckoffs[occupiedsites[rep_b[0]]].image[rep_b[1]].coord_coeffs, sitevariables[rep_b[0]], fcoords_b);
  
  orient_a = *sitevariables[rep_a[0]][2];
  orient_b = *sitevariables[rep_b[0]][2];

  /* this means extra rotations due to the symmetries that generated this image from the site basis, done as an integer offset */
  rotation_a = group.wyckoffs[occupiedsites[rep_a[0]]].image[rep_a[1]].rotation_offset;
  rotation_b = group.wyckoffs[occupiedsites[rep_b[0]]].image[rep_b[1]].rotation_offset;

  flip_a = group.wyckoffs[occupiedsites[rep_a[0]]].image[rep_a[1]].flipped ^ flips[rep_a[0]];
  flip_b = group.wyckoffs[occupiedsites[rep_b[0]]].image[rep_b[1]].flipped ^ flips[rep_b[0]];

  /* this method of only looking at the 9 nearest cells fails for those with extreme angles, so I have limited to around 30 [45, as at 25/11/12] degrees */
  /* now expand to 5x5 instead of 3x3*/
  for (cellimgx=-2; cellimgx<=2; cellimgx++) { 

    /* a is fixed, copies of b are made to test for the clash */
    img_fcoords_b[0] = fcoords_b[0]+1.0*cellimgx;

    for (cellimgy=-2; cellimgy<=2; cellimgy++) { 
	  
      if ((rep_a[0]==rep_b[0]) && (rep_a[1]==rep_b[1]) && (cellimgx==0) && (cellimgy==0)) {
	//printf("discard checking for clashes of self with self\n");
	continue;
      }

      img_fcoords_b[1] = fcoords_b[1]+1.0*cellimgy;

      //printf("img_fcoords_a %f %f ", img_fcoords_a[0], img_fcoords_a[1] );

      realcoords(fcoords_a, coords_a, cellsides, cellangles);
      realcoords(img_fcoords_b, coords_b, cellsides, cellangles);

      if (pair_clash(shape_a, shape_b, coords_a, coords_b, orient_a, orient_b, rotation_a, rotation_b, flip_a, flip_b, tests)) {
	return 1;
      }

      //printf("clash ratios %f %f %d %d [%f %f] [%f %f]\n", min_clash_ratio, clash_ratio, cellimgx, cellimgy, fcoords_a[0], fcoords_a[1], img_fcoords_b[0], img_fcoords_b[1]);
    }
  }

  return 0;
}

double packing_fraction(struct shapetype *shape, struct wallpaper_group_type group, double *cellsides[DIM], double *cellangles[DIM-1],  int occupiedsites[MAXNUMOCCSITES], double *sitevariables[MAXNUMOCCSITES][3], bool flips[numoccsites], int tests, int *clash_rejection) {
  int os;
  int countreplicas=0;
  int replicaindex[MAXREPLICAS][2];
  int m, r, s;
  bool currclash = 1;
  //double fcoords[DIM];
  //double clashratio;
  //double min_clashratio=9999.999;
  int iter=0;
  //int double_count=0;
  *clash_rejection=0;

  for (os = 0; os<numoccsites; os++) {
    //printf("site variable %p\n", sitevariables[os][0]);
    //printf("site variable %f\n", *sitevariables[os][0]);
    for (m=0; m<group.wyckoffs[occupiedsites[os]].multiplicity; m++) {
      replicaindex[countreplicas+m][0]=os;
      replicaindex[countreplicas+m][1]=m;
    }
    countreplicas += group.wyckoffs[occupiedsites[os]].multiplicity;
  }
  
  //printf("count replicas %d\n", countreplicas);

  if (tests!=MAXSTEPS+1) {
    /* we skip the clash detection procedure entirely if tests=MAXSTEPS+1, which is the case when we are reloading the best structures in the cycle and across all cycles, because reloaded structures have already been assessed for clashes */
    
    while ((currclash)&&(*clash_rejection==0)&&(iter<MAXITER)) {
      iter++;
      currclash=0;    
      
      //if (iter>1) printf("iter %d min_clashratio %15.10f\n", iter, min_clashratio);
      
      //min_clashratio = 9999.999;
      
      for (r=0; r<countreplicas; r++) {
	// this line is new on 2012-12-03. It means we only check one replica of a shape against the other replicas (copies in the same cell, either from the multiplicity of the same occupied site or at other occupied sites and their multiplicities) when that replica is the 'first' (m=0) of its occupied site. If there is no clash from one replica to the others, then there won't be a clash between replicas on the same occupied site.
	if (replicaindex[r][1]==0) {
	  //if (1) {
	  for (s=0; s<countreplicas; s++) {
	    // s begins at double_count, which starts at 0. this is to avoid doubling up when checking pairs. the first time we take a shape, double_count is 0 so we check it for clashes against every replica. we then add 1 to double_count. we continue (without checking clashes) until we find a replica on a new occupied site, at which point m=0. we then check this replica against all replicas EXCEPT the original replica we checked when double_count was 0, as this pairing has already been checked for clashes.
	    if (tests==0) {
	      currclash = clash_polygon(shape,shape,group,cellsides,cellangles,occupiedsites,sitevariables,replicaindex[r],replicaindex[s],flips,tests);
	      if (currclash) {
		if (!EMPTY_OUTPUT) printf("found a clash when initialising\n");
		//exit(1);
		return -1.0;
	      }
	      /*clashratio = clash_spokes(shape,shape,group,cellsides,cellangles,occupiedsites,sitevariables,replicaindex[r],replicaindex[s],flips);
	      if (clashratio<1.0-EPS) {
		currclash=1;
		if ((WARNINGS)&&(clashratio<0.4)) {
		  printf("replica %d (%d, %d) and %d (%d, %d) clash ratio = %f cell %f %f\n", r, replicaindex[r][0], replicaindex[r][1], s, replicaindex[s][0], replicaindex[s][1], clashratio, *cellsides[0], *cellsides[1]);
		}
	      }
	      min_clashratio = MIN(clashratio, min_clashratio);*/



	    } else {
	      currclash = clash_polygon(shape,shape,group,cellsides,cellangles,occupiedsites,sitevariables,replicaindex[r],replicaindex[s],flips,tests);
	      if (currclash) {
		/* if we have found any clash, and we are not in the initial phases (tests>0), then there's no need to search through for any more. exit the two for loops and the while loop immediately. also, flag there has been a clash so we can reject the move that caused it */
		*clash_rejection=1;
		//printf("clash rejection\n");
	      }
	    }
	    if (*clash_rejection==1) break;
	  }
	}
	if (*clash_rejection==1) break;
      }
      
      //printf("minclash %f\n", min_clashratio);
      
      //if ((min_clashratio<1.0)&&(tests==0)) {
	/* if there is a clash, we rescale the cell to avoid the clash only if we have just begun, i.e. before the first step. otherwise, i.e. if there is a clash after any step, we don't rescale the cell */
	/**cellsides[0] /= min_clashratio;
	*cellsides[0] = MIN(*cellsides[0], 99.0);
	if (cellsides[1]!=cellsides[0]) { 
	  *cellsides[1] /= min_clashratio;
	  *cellsides[1] = MIN(*cellsides[1], 99.0);
	  }*/
	
	/* I HAVE NOT YET TAKEN CARE TO ENSURE THE CELL SIDES STAY WITHIN THEIR PROSCRIBED RANGE DURING THIS STEP 
	   FOR NOW, a VERY QUICK AND DIRTY VERSION OF THIS IS IMPLEMENTED IN THE ABOVE - hardcoded 99.0*/
      //} // else Don't compress at all... allow the MC to do that 
    }
  } /* else printf("no need to check for clash, just reloading best structure\n"); */
  
  //printf("iter %d area %f maxrad %f cell %f %f sintheta %f packing %f\n", iter, area(shape), shape->maxr, *cellsides[0], *cellsides[1], sin(*cellangles[0]), countreplicas*area(shape)/((*cellsides[0])*(*cellsides[1])*sin(*cellangles[0])));
  /////////////printf("cell %f %f angle %6.2f packing %f\n", *cellsides[0], *cellsides[1], *cellangles[0]*180.0/M_PI, countreplicas*area(shape)/((*cellsides[0])*(*cellsides[1])*sin(*cellangles[0])));
  
 
  if (isnan(countreplicas*area(shape)/((*cellsides[0])*(*cellsides[1])*fabs(sin(*cellangles[0]))))) {
    printf("nan encountered %f %f %f\n", (*cellsides[0]), (*cellsides[1]), fabs(sin(*cellangles[0])));
    exit(1);
  }
 
  return countreplicas*area(shape)/((*cellsides[0])*(*cellsides[1])*fabs(sin(*cellangles[0])));
}

int initialize_structure_in_group(struct shapetype *shape, struct wallpaper_group_type group, int occupiedsites[MAXNUMOCCSITES], double basis[MAXVAR], int basis_fixed[MAXVAR], double basis_range[MAXVAR][2], double *cellsides[DIM], double *cellangles[DIM-1],  double *sitevariables[MAXNUMOCCSITES][3], bool flips[MAXNUMOCCSITES], int checkfile, int *isfile) {

  double zero=0.0; /* keep constant, effectively a hidden part of the basis */

  int countvar=0;
  int countreplicas=0;
  int i;
  int os;
  char oslist[MAXFILENAME];
  char filename[MAXFILENAME];
  FILE *fp;
  int flips_int;

  for (os = 0; os<numoccsites; os++) {
    countreplicas += group.wyckoffs[occupiedsites[os]].multiplicity;
  }
  
  /* srand(RANDSEED);  always set up using the same random number seed */

  //cell angles.
  if (group.hexagonal) {
    if (!EMPTY_OUTPUT) printf("hexagonal group\n");
    basis[0] = M_PI/3.0;
    basis_range[0][0] = M_PI/3.0;
    basis_range[0][1] = M_PI/3.0;
    basis_fixed[0] = 1;
    countvar+=1; /* not actually variable, but always part of the basis, so must be included */
  } else if (group.rectangular) {
    basis[0] = M_PI/2.0;
    basis_range[0][0] = M_PI/2.0;
    basis_range[0][1] = M_PI/2.0;
    basis_fixed[0] = 1;
    countvar+=1; /* not actually variable, but always part of the basis, so must be included */
  } else {
    basis_range[0][0] = M_PI/4.0;
    basis_range[0][1] = 3.0*M_PI/4.0;
    basis[0] = basis_range[0][0] + Fluke()*(basis_range[0][1]-basis_range[0][0]); /*oblique. M_PI/2.5 angle is starting value only, this will be varied */
    basis_fixed[0] = 0;
    countvar+=1;
  }
  cellangles[0] = &basis[0];

  //cell sides.
  if (group.a_b_equal) {
    if (!EMPTY_OUTPUT) printf("cell sides equal\n");
    countvar+=1;
    basis[1] = 4.0*shape->maxr*countreplicas;     /* NEW: 2012-12-10: starting cell edge length, start large to avoid overlaps, and there will be an early scaling down to remove free space and get a reasonable cell size. determine initial size based on number of replicas per cell */ /* 2018-05-29 sqrt removed because some wyckoffs don't spread out in both dimensions */
    basis_range[1][0] = 0.1;
    basis_range[1][1] = basis[1];
    basis_fixed[1] = 0;
    cellsides[0] = &basis[1];
    cellsides[1] = &basis[1];
  } else {
    countvar+=2;
    basis[1] = 4.0*shape->maxr*countreplicas;    
    basis[2] = 4.0*shape->maxr*countreplicas;     /* NEW 2012-12-10: starting cell edge length */ /* 2018-05-29 sqrt removed because some wyckoffs don't spread out in both dimensions */ 
    basis_range[1][0] = 0.1;
    basis_range[1][1] = basis[1];
    basis_range[2][0] = 0.1;
    basis_range[2][1] = basis[2];
    basis_fixed[1] = 0;
    basis_fixed[2] = 0;
    cellsides[0] = &basis[1];
    cellsides[1] = &basis[2];
  } 

  //now position the particles.
  for (os = 0; os<numoccsites; os++) {
    countreplicas += group.wyckoffs[occupiedsites[os]].multiplicity;
    if (!EMPTY_OUTPUT) printf("%c ", group.wyckoffs[occupiedsites[os]].letter);
    sitevariables[os][0] = &zero; //initialize, for some sites there are no site variables
    sitevariables[os][1] = &zero;
    sitevariables[os][2] = &zero;
    if (fabs(group.wyckoffs[occupiedsites[os]].image[0].coord_coeffs[0][0])>0.1) {
      /* then x is variable*/
      basis[countvar]=Fluke();
      basis_range[countvar][0] = 0.0;
      basis_range[countvar][1] = 1.0;
      basis_fixed[countvar] = 0;
      sitevariables[os][0]=&basis[countvar];
      //printf("site variable %p\n", sitevariables[os][0]);
      if (!EMPTY_OUTPUT) printf("site x variable[%d] %f\n", os, *sitevariables[os][0]);
      countvar++;
    }
    if (fabs(group.wyckoffs[occupiedsites[os]].image[0].coord_coeffs[1][1])>0.1) {
      /* then y is variable*/
      basis[countvar]=Fluke();
      basis_range[countvar][0] = 0.0;
      basis_range[countvar][1] = 1.0;
      basis_fixed[countvar] = 0;
      sitevariables[os][1]=&basis[countvar];
      if (!EMPTY_OUTPUT) printf("site y variable[%d] %f\n", os, *sitevariables[os][1]);
      countvar++;
    } 
      
	/*choose the orientation of the zeroth image of this particle */ /* This could also be done within the Wyckoff position SITEROTATION parameter?... */
	if (group.wyckoffs[occupiedsites[os]].sitemirrors) {
	  if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror0) {
		basis[countvar]=0.0;
	  } else if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror90) {
		basis[countvar]= M_PI/2.0;
	  } else if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror45) {
		basis[countvar]= 3.0*M_PI/4.0; /*not sure why I need this and the next one the wrong way around... but it works */
	  } else if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror135) {
		basis[countvar]= M_PI/4.0;
	  } else if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror30) {
		basis[countvar]= M_PI/6.0;
	  } else if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror60) {
		basis[countvar]= M_PI/3.0;
	  } else if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror330) {
		basis[countvar]= 11.0*M_PI/6.0;
	  } else if (group.wyckoffs[occupiedsites[os]].image[0].sitemirror300) {
		basis[countvar]= 5.0*M_PI/3.0;
	  }
	  basis_range[countvar][0] = basis[countvar];
	  basis_range[countvar][1] = basis[countvar];
	  basis_fixed[countvar] = 2; /*this is a different kind of fixed ... we can rotate by certain angles, as long as we match the sitemirrors */
	  sitevariables[os][2]=&basis[countvar]; 
	} else {
	  basis[countvar]=Fluke()*2.0*M_PI;
	  basis_range[countvar][0] = 0.0;
	  basis_range[countvar][1] = 2.0*M_PI;
	  basis_fixed[countvar] = 0; /* TO DO: new type of fixing - periodic in 2Pi */
	  sitevariables[os][2]=&basis[countvar];
	  if (!EMPTY_OUTPUT) printf("site offset-angle is variable[%d] %f\n", os, *sitevariables[os][2]);
	}
	countvar++;
	/* not always variable, but part of the basis, so should be included?? */
  }
  
  if (!EMPTY_OUTPUT) printf("replicas %d variables %d os ", countreplicas, countvar);
  fflush(NULL);

  for (i=0; i<numoccsites; i++) {
    oslist[i]=group.wyckoffs[occupiedsites[i]].letter;
    flips[i] = 0;
  }
  oslist[numoccsites] = '\0';

  sprintf(filename, "%s/solution_%s_%s_%s_%.3f.svg", directory, shape->name, group.label, oslist, shape_var);

  if (checkfile==1) {
    /* try to load existing solution */
    if ((fp=fopen(filename, "r"))) {
      if (!EMPTY_OUTPUT) printf("try to load existing solution\n");
      *isfile=1;
      char line [ 1024 ]; /* doesn't matter that the really long lines don't fit in one array?*/
      
      while ( fgets ( line, sizeof line, fp ) != NULL ) /* read a line */
	{
	  if ((strncmp(line, "<!--Copyright: Toby Hudson",26)==0)) {
	    if (!EMPTY_OUTPUT) printf("found copyright line, next line is data\n");
	    for (i=0; i<countvar; i++) {
	      if ((fscanf(fp, "%lf ", &basis[i]))==EOF) {
		fprintf(stderr, "ERROR: reload format incorrect\n");
		exit(1);
	      }
	    }
	    for (i=0; i<numoccsites; i++) {
	      if ((fscanf(fp, "%d ", &flips_int))==EOF) {
		fprintf(stderr, "ERROR: reload format incorrect\n");
		exit(1);		      
	      }
	      flips[i] = flips_int;
	    }
	    if (!EMPTY_OUTPUT) printf("existing structure loaded - cold quench only\n");
	  }
	}
      fclose(fp);
    } else {
      *isfile=0;
      if (!EMPTY_OUTPUT) printf("No existing solution -  File: %s cannot be loaded.\n", filename);
    }
  }
  return(countvar);
}

void uniform_best_packing_in_isopointal_group(struct shapetype *shape, struct wallpaper_group_type group, int occupiedsites[MAXNUMOCCSITES]) {
  double phi;
  double kT_start = 0.1;
  double kT_finish = 0.0005;
  double kT_ratio;
  double kT;
  int vary_index; /* which of our basis values are we currently trying to vary */
  double old_value;
  double old_phi; /* phi == packing fraction */
  double max_phi=0.0;
  double max_phi_in_cycle;
  //double old_angle;
  double old_cell0, old_cell1;
  int i, os, m;

  bool flips[MAXNUMOCCSITES];
  
  double best_basis[MAXVAR];
  double best_basis_in_cycle[MAXVAR][CYCLES];
  bool best_flips[numoccsites];
  bool best_flips_in_cycle[numoccsites][CYCLES];
 
  double fcoords[2]; /*fractional coordinates (scaled between 0 and 1) */
  double coords[2];

  double image_scale, shape_scale;
  int cell_img_x, cell_img_y; /* indexes of x and y periodic cells used to check for overlap */
  double image_width, image_height;
  double origin_translate_x, origin_translate_y;

  long int tests=0, rejections=0;
  int cycle=0, isfile=0, checkfile=1;

  int flip_site = -1;

  char oslist[MAXFILENAME];
  char filename[MAXFILENAME];
  FILE *fp;
  
  int basis_size=0;
  double basis[MAXVAR];
  int    basis_fixed[MAXVAR];
  double basis_range[MAXVAR][2];
  double *cellsides[DIM];
  double *cellangles[DIM-1];
  double *sitevariables[MAXNUMOCCSITES][3]; /* x_center, y_center, theta */

  int clash_rejection=0;

  int countreplicas=0;

  for (os = 0; os<numoccsites; os++) {
    countreplicas += group.wyckoffs[occupiedsites[os]].multiplicity;
  }
  
  /* construct the filename */
  for (i=0; i<numoccsites; i++) {
    oslist[i]=group.wyckoffs[occupiedsites[i]].letter;
  }
  oslist[numoccsites] = '\0';
  sprintf(filename, "%s/solution_%s_%s_%s_%.3f.svg", directory, shape->name, group.label, oslist, shape_var);


  for (cycle=0; cycle<CYCLES; cycle++) {
    //if ((cycle!=0)&&(reload=1)) reload=0;
    tests = 0;

    phi=-1.0;
    while (phi<0.0) {
      if ((cycle==0) || (isfile==0)) basis_size = initialize_structure_in_group(shape, group, occupiedsites, basis, basis_fixed, basis_range, cellsides, cellangles, sitevariables, flips, checkfile, &isfile);  

      checkfile=0;

      if (isfile==1) {
	kT_start=kT_finish*10.0;
      }

      phi = packing_fraction(shape, group, cellsides, cellangles, occupiedsites, sitevariables, flips, tests, &clash_rejection);
    }
    if (!EMPTY_OUTPUT) printf("initial packing fraction = %f\n", phi);
    
    kT_ratio = pow(kT_finish/kT_start, 1.0/MAXSTEPS); 
      
    kT = kT_start;
    rejections=0;
    max_phi_in_cycle=0.0;
    
    while (tests<MAXSTEPS) {
      kT *= kT_ratio;
      clash_rejection=0;

      vary_index=-1; 
      while (vary_index==-1) vary_index = rand()%basis_size; 
      //if (tests%100==0) vary_index=0; /* only occasionally allow the cell angle to change */
      
      while (basis_fixed[vary_index]==1) {
	/* search until a variable basis is found */
	vary_index = rand()%basis_size; 
      }
      
      if ((tests%100)&&(ALLOWFLIPS)) {
	flip_site = rand()%numoccsites;
	flips[flip_site] ^= 1;
      } else {
	flip_site = -1;
      }
      
      old_value = basis[vary_index];
      old_cell0 = *cellsides[0];
      old_cell1 = *cellsides[1];
      
      old_phi = phi;
      
      if (basis_fixed[vary_index]==0) { 
	if (basis_range[vary_index][1]<3.0) {
	  /* This is a regular structural parameter (non cell-side) which can be varied within its given range */
	  //basis[vary_index] = old_value *(1.0 + max_step_size*(Fluke()-0.5));
	  basis[vary_index] = old_value + max_step_size*(basis_range[vary_index][1]-basis_range[vary_index][0])*(Fluke()-0.5);
	  //printf("A kT = %g, phi=%f %d %f, (%f, %f) theta=%f\n", kT, phi, vary_index, old_value, basis_range[vary_index][0], basis_range[vary_index][1], *cellangles[0]);
	} else {
	  /* This is a cell side.  Since these have a direct effect on the packing fraction, make their step size in proportion to kT */
	  basis[vary_index] = old_value *(1.0 + MIN(3.0*kT,0.1)*(Fluke()-0.5));		
	}
	
	if (basis[vary_index]<basis_range[vary_index][0]) {
	  basis[vary_index]=basis_range[vary_index][0];
	} else if (basis[vary_index]>basis_range[vary_index][1]) {
	  basis[vary_index]=basis_range[vary_index][1];
	}
	
	if (vary_index==0) {
	  /* this is the cell angle. stretch cell sides so that the total area is preserved (to prevent runaway angle collapse) */
	  *cellsides[0] = old_cell0 * sqrt(sin(old_value)/sin(basis[0]));
	  *cellsides[1] = old_cell1 * sqrt(sin(old_value)/sin(basis[0]));
	}
	
      } else {
	/* basis_fixed[vary_index]==2 which means it's a particle orientation which must vary in a quantized manner to ensure site mirrors still lie on symmetry mirror planes */
	
	if ((shape->mirrors%2==0) && (Fluke()<0.5)) {
	  /* turn it 90 degrees to switch the x and y mirror planes */
	  if (basis[vary_index]<M_PI*3.0/4.0) 
	    basis[vary_index] += M_PI/(double)shape->mirrors;
	  else
	    basis[vary_index] -= M_PI/(double)shape->mirrors;
	} else {
	  /* turn it 180 degrees so that all mirror planes are preserved (just reversed) */
	  if (basis[vary_index]<M_PI) 
	    basis[vary_index] += M_PI;
	  else
	    basis[vary_index] -= M_PI;
	} 
      }
      
      tests++; /* increase tests here so that when calculate packing fraction on first test, tests=1 */
      
      phi = packing_fraction(shape, group, cellsides, cellangles, occupiedsites, sitevariables, flips, tests, &clash_rejection);
      
      if ((clash_rejection==1) || (Fluke()>exp(((1.0/old_phi-1.0/phi)/kT) + countreplicas*log(old_phi/phi)))) {
	// if (phi<old_phi) {
	/* some probability of rejecting move if packing fraction decreased, or definitely reject it if it caused a clash */
	rejections++;
	basis[vary_index] = old_value;
	*cellsides[0] = old_cell0;
	*cellsides[1] = old_cell1;
	if (flip_site != -1) flips[flip_site]^=1;
	phi = old_phi;
	/*if (clash_rejection!=1) {
	  printf("######not a clash rejection, return to old values\n");
	}
	if (clash_rejection==1) {
	  printf("############################a clash rejection, return to old values\n");
	  }*/
      }

      if (phi > max_phi) {
	/* best packing seen yet ... save data */
	for (i=0; i<basis_size; i++) {
	  best_basis[i] = basis[i];
	}
	for (i=0; i<numoccsites; i++) {
	  best_flips[i] = flips[i];
	} 
	max_phi = phi;
      } 
      
      if (phi > max_phi_in_cycle) {
	/* best packing seen yet in current cycle ... save data */
	for (i=0; i<basis_size; i++) {
	  best_basis_in_cycle[i][cycle] = basis[i];
	}
	for (i=0; i<numoccsites; i++) {
	  best_flips_in_cycle[i][cycle] = flips[i];
	} 
	max_phi_in_cycle = phi;
      }

      if (!EMPTY_OUTPUT) {
	if (tests%500==0) printf("cycle %d of %d, step %ld of %d, kT=%g, packing %f, angle %f, b/a=%f, rejection %f percent\n", cycle+1, CYCLES, tests, MAXSTEPS, kT, phi, *cellangles[0]*180.0/M_PI, (*cellsides[1])/(*cellsides[0]), (100.0*rejections)/tests);
      }
    }

    tests=MAXSTEPS+1;
    /* add 1 to tests so that tests=MAXSTEPS+1, to signify to packing_fraction function not to check for clashes when reloading best structures */

    /*now reload the best structure ever seen in the current cycle*/
    for (i=0; i<basis_size; i++) {
      basis[i] = best_basis_in_cycle[i][cycle];
    }
    for (i=0; i<numoccsites; i++) {
      flips[i] = best_flips_in_cycle[i][cycle];
    }
    
    phi = packing_fraction(shape, group, cellsides, cellangles, occupiedsites, sitevariables, flips, tests, &clash_rejection);
    if (!EMPTY_OUTPUT) {
      printf("BEST: cell %f %f angle %6.2f packing %f rejection (%f\%%) ", *cellsides[0], *cellsides[1], *cellangles[0]*180.0/M_PI, max_phi_in_cycle, (100.0*rejections)/tests);
      for (i=0; i<basis_size; i++) {
	printf("%f ", basis[i]);
      }
      for (i=0; i<numoccsites; i++) {
	printf("%d ", flips[i]);
      }
      printf("\n");
    }
    #define btoa(x) ((x)?"true":"false")
    char chiral_state;
    int chiralsum=0;
    int totalsum=0;
    for ( os=0; os<numoccsites; os++ ) {
        chiralsum+=(best_flips[os]*2-1)*group.wyckoffs[occupiedsites[os]].multiplicity;
        totalsum+=abs((best_flips[os]*2-1)*group.wyckoffs[occupiedsites[os]].multiplicity);
    }
    if ( chiralsum!=0 ) {
        chiral_state = 'c';
        if (chiralsum!=totalsum) {
            chiral_state='s';
        }
    }
    else {
        chiral_state = 'a';
    }
    
    fprintf(fpout,"%8s\t%s\t%d\t%.6f\t%.7f\t%.3f\t%s\t%d\t%.5f\t%d\t%.3f\t%d\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%c", group.label, oslist, cycle+1, max_phi, area(shape), (100.0*rejections)/tests, shape->name, MAXSTEPS, max_step_size, SHAPE_RESOLUTION, shape_var, numoccsites, shape->minr, shape->maxr, shape->minr/shape->maxr, max_phi_in_cycle, *cellsides[0], *cellsides[1], *cellangles[0]*180.0/M_PI, chiral_state);
    
    for ( os=0; os<numoccsites; os++ ) {
            fprintf(fpout, "\t%s\t%d", btoa(best_flips[os]), group.wyckoffs[occupiedsites[os]].multiplicity);
    }
    fprintf(fpout,"\n");
    fprintf(fp_isopointal_out,"%.3f\t%.6f\t%.7f\t%.3f\t%s\t%8s\t%s\t%d\t%d\t%.5f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\n", shape_var, max_phi, area(shape), (100.0*rejections)/tests, shape->name, group.label, oslist, cycle+1, MAXSTEPS, max_step_size, SHAPE_RESOLUTION, numoccsites, shape->minr, shape->maxr, shape->minr/shape->maxr, max_phi_in_cycle, *cellsides[0], *cellsides[1], *cellangles[0]*180.0/M_PI);

    fflush(NULL);
    
    
    /*now reload the best structure ever seen, for the image*/
    for (i=0; i<basis_size; i++) {
      basis[i] = best_basis[i];
    }
    for (i=0; i<numoccsites; i++) {
      flips[i] = best_flips[i];
    }
    phi = packing_fraction(shape, group, cellsides, cellangles, occupiedsites, sitevariables, flips, tests, &clash_rejection);
    
    
    image_scale=(MAXGRAPHICSIZE-2*GRAPHICFRAME)/(3.0*MAX(*cellsides[0]+fabs(*cellsides[1]*cos(*cellangles[0])),*cellsides[1]*sin(*cellangles[0])));
    shape_scale= 1.0;  /* (image_scale-1.0)/image_scale;  this was originally an attempt to correct for the 1px stroke width, but since that applies to all directions additively, it is better to remove it in the definition of the shape itself rather than try to scale it */
    origin_translate_x = GRAPHICFRAME+MAX(0.0,-image_scale**cellsides[1]*cos(*cellangles[0]));
    origin_translate_y = GRAPHICFRAME;
    image_width  = 2*GRAPHICFRAME+(MAXGRAPHICSIZE-2*GRAPHICFRAME)*(*cellsides[0]+fabs(*cellsides[1]*cos(*cellangles[0])))/MAX(*cellsides[0]+fabs(*cellsides[1]*cos(*cellangles[0])),*cellsides[1]*sin(*cellangles[0]));
    image_height = 2*GRAPHICFRAME+(MAXGRAPHICSIZE-2*GRAPHICFRAME)*(*cellsides[1]*sin(*cellangles[0]))/MAX(*cellsides[0]+fabs(*cellsides[1]*cos(*cellangles[0])),*cellsides[1]*sin(*cellangles[0]));
    
    fp = fopen(filename, "w");
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
    fprintf(fp, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
    fprintf(fp, "<svg width=\"%f\" height=\"%f\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", image_width, image_height);
    fprintf(fp, "<defs>\n");
    fprintf(fp, "<g id=\"cell\" opacity=\"1.0\" >\n");
    fprintf(fp, "<polygon points=\"0,0 %f,0 %f,%f %f,%f\" style=\"stroke:#000000;stroke-width:2\"/>\n", image_scale**cellsides[0], image_scale**cellsides[0]+image_scale**cellsides[1]*cos(*cellangles[0]), image_scale**cellsides[1]*sin(*cellangles[0]), image_scale**cellsides[1]*cos(*cellangles[0]),image_scale**cellsides[1]*sin(*cellangles[0]) );
    fprintf(fp, "</g>\n");
    fprintf(fp, "<g id=\"shape\" opacity=\"0.7\" >\n");
    if (POLYGON_REP==1) {
      fprintf(fp, "<polygon points=\"");
      for (i=0; i<SHAPE_RESOLUTION; i++) {
	fprintf(fp, "%f,%f ", (image_scale*shape->r[i]-0.5)*cos(i*2.0*M_PI/SHAPE_RESOLUTION),  (image_scale*shape->r[i]-0.5)*sin(i*2.0*M_PI/SHAPE_RESOLUTION));
      }
      fprintf(fp, "\" style=\"stroke:#000000;stroke-width:1\"/>\n");
    } else if (POLYGON_REP==2) {
      fprintf(fp, "<polygon points=\"");
      for (i=0; i<SHAPE_RESOLUTION; i++) {
	fprintf(fp, "%f,%f ", image_scale*shape->r[i]*cos(i*2.0*M_PI/SHAPE_RESOLUTION),  image_scale*shape->r[i]*sin(i*2.0*M_PI/SHAPE_RESOLUTION));
      }
      fprintf(fp, "\" style=\"stroke:#000000;stroke-width:0\"/>\n");
    } else {
      /* draw radial lines */
      for (i=0; i<SHAPE_RESOLUTION; i++) {
	fprintf(fp, "<line x1=\"0\" y1=\"0\" x2=\"%f\" y2=\"%f\" stroke=\"purple\" stroke-width=\"1\" />\n", image_scale*shape->r[i]*cos(i*2.0*M_PI/SHAPE_RESOLUTION),  image_scale*shape->r[i]*sin(i*2.0*M_PI/SHAPE_RESOLUTION));
      }    
    }
    if (AXES_SHOWN) { 
      fprintf(fp, "<line x1=\"0\" y1=\"0\" x2=\"20\" y2=\"0\"  stroke=\"black\" stroke-width=\"1\"  />\n");
      fprintf(fp, "<line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"20\"  stroke=\"red\" stroke-width=\"1\"  />\n");
    }
    if (CENTRE_SPOT) {
      fprintf(fp, "<circle cx=\"0.0\" cy=\"0.0\" r=\"3.0\" fill=\"black\"  stroke=\"black\" stroke-width=\"0\" />\n");
    }
    fprintf(fp, "</g>\n");
    fprintf(fp, "</defs>\n");
    for (cell_img_x=0; cell_img_x<3; cell_img_x++) {
      for (cell_img_y=0; cell_img_y<3; cell_img_y++) {
	fcoords[0] = 1.0*cell_img_x;
	fcoords[1] = 1.0*cell_img_y;
	realcoords(fcoords, coords, cellsides, cellangles);
	if (cell_img_x==1 && cell_img_y==1) 
	  fprintf(fp, "<use xlink:href=\"#cell\" transform=\"translate(%f,%f) scale(1)\" style=\"fill:grey\"/>\n", origin_translate_x+image_scale*coords[0], origin_translate_y+image_scale*coords[1]);
	else
	  fprintf(fp, "<use xlink:href=\"#cell\" transform=\"translate(%f,%f) scale(1)\" style=\"fill:white\"/>\n", origin_translate_x+image_scale*coords[0], origin_translate_y+image_scale*coords[1]);
      }
    }
    for (os=0; os<numoccsites; os++) {
      for (m=0; m<group.wyckoffs[occupiedsites[os]].multiplicity; m++) {
	fractcoords(group.wyckoffs[occupiedsites[os]].image[m].coord_coeffs, sitevariables[os], fcoords);
	for (cell_img_x=0; cell_img_x<3; cell_img_x++) {
	  for (cell_img_y=0; cell_img_y<3; cell_img_y++) {
	    fcoords[0] += 1.0*cell_img_x;
	    fcoords[1] += 1.0*cell_img_y;
	    realcoords(fcoords, coords, cellsides, cellangles);
	    fcoords[0] -= 1.0*cell_img_x;
	    fcoords[1] -= 1.0*cell_img_y;
	    fprintf(fp, "<use xlink:href=\"#shape\" transform=\"translate(%f,%f) rotate(%f, 0, 0) scale(%f %f)\" style=\"fill:", 
		    origin_translate_x+image_scale*coords[0], 
		    origin_translate_y+image_scale*coords[1], 
		    -((group.wyckoffs[occupiedsites[os]].image[m].flipped^flips[os])*(-2)+1)*(*sitevariables[os][2]*180.0/M_PI)-group.wyckoffs[occupiedsites[os]].image[m].rotation_offset*360.0/SHAPE_RESOLUTION, 
		    shape_scale, 
		    ((group.wyckoffs[occupiedsites[os]].image[m].flipped^flips[os])*(-2)+1)*shape_scale
		    );	
	    if (cell_img_x==1 && cell_img_y==1) {
	      fprintf(fp,"rgb(%d,%d,%d)", (int)(0.1*256/(os+1)),  (int)(0.3*256/(os+1)),  (int)(0.6*256/(os+1)));
	    } else {
	      fprintf(fp,"rgb(%d,%d,%d)", (int)(0.1*256/(os+1)),  (int)(0.6*256/(os+1)),  (int)(0.3*256/(os+1)));
	    }
	    fprintf(fp, "\"/>\n");
	    
	  }
	}
      }
    }
    fprintf(fp, "</svg>\n");
    fprintf(fp, "<!--Copyright: Toby Hudson\n");
    for (i=0; i<basis_size; i++) {
      fprintf(fp, "%f ", basis[i]);
    }
    for (i=0; i<numoccsites; i++) {
      fprintf(fp, "%d ", flips[i]);
    }
    fprintf(fp, "\ncell %f %f angle %6.2f packing %f (max %f) rejection (%f\%%) ", *cellsides[0], *cellsides[1], *cellangles[0]*180.0/M_PI, phi, max_phi, (100.0*rejections)/tests);
    fprintf(fp, "-->\n");
    fclose(fp);
  }
}


void uniform_best_packing_in_group(struct shapetype *shape, struct wallpaper_group_type group) {
  bool validsite[group.num_wyckoffs];
  int numvalid=0, numvalid_and_variable=0;
  int wk, os;

  int n, a, v, maxcombinations, numcombinations;

  int combination_bag[MAXNUMOCCSITES*group.num_wyckoffs]; /* temporary variable only */
  int curr_bag=0;
  
  int i;

  // first test if the shape has the required symmetries for various sites
  // at the moment only tests rotations & mirrors... that's all?
  for (wk=0; wk<group.num_wyckoffs; wk++) {
    if ((shape->rotsym%group.wyckoffs[wk].siterotations==0)&&
	((group.wyckoffs[wk].sitemirrors==0)||((shape->mirrors%group.wyckoffs[wk].sitemirrors==0)&&(shape->mirrors!=0)))) {
      if (group.wyckoffs[wk].somevariability) {
	printf("site %c is valid and variable\n", group.wyckoffs[wk].letter);
      } else {
	printf("site %c is valid\n", group.wyckoffs[wk].letter);
      }
      validsite[wk] = 1;
    } else {
      printf("site %c does not have the required rotational symmetry\n", group.wyckoffs[wk].letter);
      validsite[wk] = 0;
    }
  }

  /* Count up how many were valid in total,
     and amongst those, how many had some variability (which means they can be mulitply occupied). */
  /* The last wyckoff site should always be the general one,
     which permits any particle symmetry and is always variable. */
  for (wk=0; wk<group.num_wyckoffs; wk++) {
    if (validsite[wk]) {
      numvalid++;
      if (group.wyckoffs[wk].somevariability) {
	for (os=0; os<numoccsites; os++) {
	  combination_bag[curr_bag] = wk;
	  curr_bag++;
	}
	numvalid_and_variable++;
      } else {
	/* since not variable, only add one option to the bag for this */
	  combination_bag[curr_bag] = wk;
	  curr_bag++;	
      }
    }
  }

  v = numvalid - numvalid_and_variable;
  a = numvalid_and_variable;
  n = numoccsites;

  maxcombinations = factorial(n*a+v)/(factorial(n)*factorial(n*a+v-n));

  if (n*a+v > 12) {
	// factorial fails
	maxcombinations = 2000;
  }

  printf("estimate: there should be less than %d ways of filling %d sites with choices from %d variable wyckoffs and %d non-variable wyckoffs \n", maxcombinations, n, a, v);
  
  int occupiedsites[maxcombinations][MAXNUMOCCSITES];

  numcombinations = enumerate(occupiedsites, maxcombinations, 0, 0, 0, combination_bag, curr_bag);
  
  printf("enumeration complete: in fact there were %d ways\n", numcombinations);


  for (i=0; i<numcombinations; i++) {
    if (!EMPTY_OUTPUT) {
      printf("occupiedsites[%d]", i);
      for (os=0; os<numoccsites; os++) printf(" %c", group.wyckoffs[occupiedsites[i][os]].letter);
      printf("\n");
    }

    char format[] = "%s/output_%s_%s_%c_%d.dat"; //DAN this paragraph
    char filename[sizeof format+200];
    sprintf(filename,format,directory,shape->name,group.label,group.wyckoffs[occupiedsites[i][0]].letter,numoccsites);

    fp_isopointal_out = fopen(filename, "a"); //DAN
 
    uniform_best_packing_in_isopointal_group(shape, group, occupiedsites[i]);
    fclose(fp_isopointal_out);
  }
}


int main (int argc, char *argv[]) {

  struct shapetype fourier;
  struct shapetype circle;
  struct shapetype square;
  struct shapetype wobblycircle;
  struct shapetype pointyegg;
  struct shapetype ellipse;
  struct shapetype wobbly3;
  struct shapetype wobbly4;
  struct shapetype wobbly5;
  struct shapetype wobblyfan3;
  struct shapetype ninjastar;
  struct shapetype asymcomma;
  struct shapetype pentagon;
  struct shapetype pentagoncircle;
  struct shapetype heptagoncircle;
  struct shapetype bowtie;
  struct shapetype zigzagbone;
  struct shapetype bells;
  struct shapetype bone;
  struct shapetype triangle;
  struct shapetype smoothedpentagon;
  struct shapetype hexagon;
  struct shapetype heptagon;
  struct shapetype nonagon;
  struct shapetype octagon;
  //struct shapetype modifiedoctagon;
  struct shapetype randomsmoothedoctagon;
  struct shapetype smoothedoctagon;
  struct shapetype dodecagon;
  struct shapetype fan2;  
  struct shapetype fan3;
  struct shapetype fan4;
  struct shapetype flower5;
  struct shapetype fan6;
  struct shapetype fan12;
  struct shapetype kite;
  struct shapetype RJpentagon;
  struct shapetype RJ2pentagon;
  struct shapetype MRpentagon;

  struct shapetype *chosen_shape;

  char datfilename[200];
  char command[200];
  double fcoeff[2][FOURIER_TERMS];
  int dd, ff;

  //fprintf(stderr, "%d\n", argv[1]);

  int do_all;
  int study_wallpaper_group=0;

  numoccsites = atoi(argv[1]);

  if (strcmp(argv[2], "fourier")==0) {
    if (argc!=22) {
      fprintf(stderr, "wrong number of arguments. wallpaper.exe num_occ_sites shape_name wallpaper_group_num/all/best5/best8 coefficients(a1->a9,b1->b9)\n");
      exit(1);
    }
  } else {
    if (argc!=4) {
      fprintf(stderr, "wrong number of arguments. wallpaper.exe num_occ_sites shape_name wallpaper_group_num/all/best5/best8\n");
      exit(1);
    }
  }

  if (strcmp(argv[2], "fourier")==0) {
    for (dd=0; dd<2; dd++) {
      for (ff=0; ff<FOURIER_TERMS; ff++) {
	fcoeff[dd][ff] = atof(argv[4+FOURIER_TERMS*dd+ff]);
      }
    }
    define_fourier(&fourier, fcoeff);
    chosen_shape = &fourier;
    //fouriers will make their own subdirectory, to avoid collisions
    sprintf(command, "mkdir %4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f", fcoeff[0][0], fcoeff[0][1], fcoeff[0][2], fcoeff[0][3], fcoeff[0][4], fcoeff[0][5], fcoeff[0][6], fcoeff[0][7], fcoeff[0][8], fcoeff[1][0], fcoeff[1][1], fcoeff[1][2], fcoeff[1][3], fcoeff[1][4], fcoeff[1][5], fcoeff[1][6], fcoeff[1][7], fcoeff[1][8]);
    sprintf(directory, "%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f", fcoeff[0][0], fcoeff[0][1], fcoeff[0][2], fcoeff[0][3], fcoeff[0][4], fcoeff[0][5], fcoeff[0][6], fcoeff[0][7], fcoeff[0][8], fcoeff[1][0], fcoeff[1][1], fcoeff[1][2], fcoeff[1][3], fcoeff[1][4], fcoeff[1][5], fcoeff[1][6], fcoeff[1][7], fcoeff[1][8]);
    system(command);
  } else if (strcmp(argv[2], "circle")==0) {
    define_circle(&circle);
    chosen_shape = &circle;
  } else if (strcmp(argv[2], "square")==0) {
    define_square(&square);
    chosen_shape = &square;
  } else if (strcmp(argv[2], "wobblycircle")==0) {
    define_wobblycircle(&wobblycircle);
    chosen_shape = &wobblycircle;
  } else if (strcmp(argv[2], "pointyegg")==0) {
    define_pointyegg(&pointyegg);
    chosen_shape = &pointyegg;
  } else if (strcmp(argv[2], "ellipse")==0) {
    define_ellipse(&ellipse);
    chosen_shape = &ellipse;
  } else if (strcmp(argv[2], "wobbly3")==0) {
    define_wobbly3(&wobbly3);
    chosen_shape = &wobbly3;
  } else if (strcmp(argv[2], "wobbly4")==0) {
    define_wobbly4(&wobbly4);
    chosen_shape = &wobbly4;
  } else if (strcmp(argv[2], "wobbly5")==0) {
    define_wobbly5(&wobbly5);
    chosen_shape = &wobbly5;
  } else if (strcmp(argv[2], "wobblyfan3")==0) {
    define_wobblyfan3(&wobblyfan3);
    chosen_shape = &wobblyfan3;
  } else if (strcmp(argv[2], "ninjastar")==0) {
    define_ninjastar(&ninjastar);
    chosen_shape = &ninjastar;
  } else if (strcmp(argv[2], "asymcomma")==0) {
    define_asymcomma(&asymcomma);
    chosen_shape = &asymcomma;
  } else if (strcmp(argv[2], "pentagon")==0) {
    define_pentagon(&pentagon);
    chosen_shape = &pentagon;
  } else if (strcmp(argv[2], "pentagoncircle")==0) {
    define_pentagoncircle(&pentagoncircle);
    chosen_shape = &pentagoncircle;
  } else if (strcmp(argv[2], "heptagoncircle")==0) {
    define_heptagoncircle(&heptagoncircle);
    chosen_shape = &heptagoncircle;
  } else if (strcmp(argv[2], "bowtie")==0) {
    define_bowtie(&bowtie); 
    chosen_shape = &bowtie;
  } else if (strcmp(argv[2], "zigzagbone")==0) {
    define_zigzagbone(&zigzagbone);
    chosen_shape = &zigzagbone;
  } else if (strcmp(argv[2], "bells")==0) {
    define_bells(&bells);
    chosen_shape = &bells;
  } else if (strcmp(argv[2], "bone")==0) {
    define_bone(&bone);
    chosen_shape = &bone;
  } else if (strcmp(argv[2], "triangle")==0) {
    define_triangle(&triangle);
    chosen_shape = &triangle;
  } else if (strcmp(argv[2], "smoothedpentagon")==0) {
    define_smoothedpentagon(&smoothedpentagon);
    chosen_shape = &smoothedpentagon;
  } else if (strcmp(argv[2], "hexagon")==0) {
    define_hexagon(&hexagon);
    chosen_shape = &hexagon;
  } else if (strcmp(argv[2], "heptagon")==0) {
    define_heptagon(&heptagon);
    chosen_shape = &heptagon;
  } else if (strcmp(argv[2], "nonagon")==0) {
    define_nonagon(&nonagon);
    chosen_shape = &nonagon;
  } else if (strcmp(argv[2], "octagon")==0) {
    define_octagon(&octagon);
    chosen_shape = &octagon;
  } else if (strcmp(argv[2], "randomsmoothedoctagon")==0) {
    // define_randomsmoothedoctagon(&randomsmoothedoctagon);
    chosen_shape = &randomsmoothedoctagon;
  } else if (strcmp(argv[2], "smoothedoctagon")==0) {
    define_smoothedoctagon(&smoothedoctagon);
    chosen_shape = &smoothedoctagon;
  } else if (strcmp(argv[2], "dodecagon")==0) {
    define_dodecagon(&dodecagon);
    chosen_shape = &dodecagon;
  } else if (strcmp(argv[2], "fan2")==0) {
    define_fan2(&fan2);
    chosen_shape = &fan2;
  } else if (strcmp(argv[2], "fan3")==0) {
    define_fan3(&fan3);
    chosen_shape = &fan3;
  } else if (strcmp(argv[2], "fan4")==0) {
    define_fan4(&fan4);
    chosen_shape = &fan4;
  } else if (strcmp(argv[2], "flower5")==0) {
    define_flower5(&flower5);
    chosen_shape = &flower5;
  } else if (strcmp(argv[2], "fan6")==0) {
    define_fan6(&fan6);
    chosen_shape = &fan6;
  } else if (strcmp(argv[2], "fan12")==0) {
    define_fan12(&fan12);
    chosen_shape = &fan12;
  } else if (strcmp(argv[2], "kite")==0) {
    define_kite(&kite);
    chosen_shape = &kite;
  } else if (strcmp(argv[2], "RJpentagon")==0) {
    define_RichardJames_pentagon(&RJpentagon, 1.4); /* a bit of an arbitrary one */
    chosen_shape = &RJpentagon;
  } else if (strcmp(argv[2], "RJ2pentagon")==0) {
    define_RichardJames_pentagon(&RJ2pentagon, 1.2); /* a bit of an arbitrary one */
    chosen_shape = &RJ2pentagon;
  } else if (strcmp(argv[2], "MRpentagon")==0) {
    define_MarjorieRice_pentagon(&MRpentagon, 1.9); /* angle must be greater than pi/2 */
    chosen_shape = &MRpentagon;
  } else {
    fprintf(stderr, "shape not recognized\n");
    exit(1);
  }
  sprintf(datfilename, "%s/%s_%d_%.3f.dat", directory, argv[2], SHAPE_RESOLUTION, shape_var);
  plot_shape(chosen_shape, datfilename);
  
  
  if (strcmp(argv[3], "all")==0) {
    do_all=1;
  } else if (strcmp(argv[3], "best5")==0) {     //DAN
    do_all=2;
  } else if (strcmp(argv[3], "best8")==0) {     //DAN
    do_all=3;
  } else if (strcmp(argv[3], "chiral5")==0) {     //DAN
    do_all=4;
  } else if (strcmp(argv[3], "worst9")==0) {     //DAN
    do_all=5;
  } else if (strcmp(argv[3], "chiral10")==0) {     //DAN
    do_all=6;
  } else if (strcmp(argv[3], "worst7")==0) {     //DAN
    do_all=7;
  } else if (strcmp(argv[3], "chiraltest")==0) {   //CLARE
    do_all=8;
  } else {
    do_all=0;
    study_wallpaper_group = atoi(argv[3]);
  }

  /*
  plot_shape(&unitcircle, "unitcircle.dat");
  plot_shape(&square, "square.dat");
  plot_shape(&wobblycircle, "wobblycircle.dat");
  plot_shape(&pointyegg, "pointyegg.dat");
  plot_shape(&ellipse, "ellipse.dat");
  plot_shape(&wobbly3, "wobbly3.dat");
  plot_shape(&wobbly4, "wobbly4.dat");
  plot_shape(&wobbly5, "wobbly5.dat");
  plot_shape(&wobblyfan3, "wobblyfan3.dat");
  plot_shape(&ninjastar, "ninjastar.dat");
  plot_shape(&asymcomma, "asymcomma.dat");
  plot_shape(&pentagon, "pentagon.dat");
  plot_shape(&bowtie, "bowtie.dat");
  plot_shape(&zigzagbone, "zigzagbone.dat");
  plot_shape(&bells, "bells.dat");
  plot_shape(&bone, "bone.dat");
  plot_shape(&triangle, "triangle.dat");
  plot_shape(&smoothedpentagon, "smoothedpentagon.dat");
  plot_shape(&hexagon, "hexagon.dat");
  plot_shape(&heptagon, "heptagon.dat");
  plot_shape(&octagon, "octagon.dat");
  plot_shape(&randomsmoothedoctagon, "randomsmoothedoctagon.dat");
  plot_shape(&smoothedoctagon, "smoothedoctagon.dat");
  plot_shape(&dodecagon, "dodecagon.dat");
  plot_shape(&fan2, "fan2.dat");
  plot_shape(&fan3, "fan3.dat");
  plot_shape(&fan4, "fan4.dat");
  plot_shape(&flower5, "flower5.dat");
  plot_shape(&fan6, "fan6.dat");
  plot_shape(&fan12, "fan12.dat");
  plot_shape(&kite, "kite.dat");
  */
  //plot_shape(&RJpentagon, "RJpentagon.dat");

  //print_group_definitions();

  char format[] = "%s/output_%s_%.3f_%d.dat";  //DAN this paragraph
  char filename[sizeof format+200];
  sprintf(filename,format,directory,chosen_shape->name,shape_var,numoccsites);

  fpout = fopen(filename, "a");


  if (do_all == 1) { //all
	for (study_wallpaper_group=0; study_wallpaper_group<17; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
  } else if (do_all == 2) {  //DAN: best5 - p1, p2, pg, p2mg, p2gg
	for (study_wallpaper_group=0; study_wallpaper_group<2; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=3; study_wallpaper_group<4; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=6; study_wallpaper_group<8; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
  } else if (do_all == 3) {  //DAN: best8 - p1, p2, pg, cm, p2mg, p2gg, p3, p31m
	for (study_wallpaper_group=0; study_wallpaper_group<2; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=3; study_wallpaper_group<5; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=6; study_wallpaper_group<8; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=12; study_wallpaper_group<13; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=14; study_wallpaper_group<15; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
  } else if (do_all == 4) {  //DAN: chiral5 - p1, p2, p3, p4, p6
	for (study_wallpaper_group=0; study_wallpaper_group<2; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=9; study_wallpaper_group<10; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=12; study_wallpaper_group<13; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=15; study_wallpaper_group<16; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
  } else if (do_all == 5) {  //DAN: worst9 - pm, p2mm, c2mm, p4, p4mm, p4gm, p3m1, p6, p6mm
	for (study_wallpaper_group=2; study_wallpaper_group<3; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=5; study_wallpaper_group<6; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=8; study_wallpaper_group<12; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=13; study_wallpaper_group<14; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=15; study_wallpaper_group<17; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
  } else if (do_all == 6) {  //DAN: chiral10 - p1, p2, pg, cm, p2mg, p2gg, p4, p3, p31m, p6
	for (study_wallpaper_group=0; study_wallpaper_group<2; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=3; study_wallpaper_group<5; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=6; study_wallpaper_group<8; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=9; study_wallpaper_group<10; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=12; study_wallpaper_group<13; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=14; study_wallpaper_group<16; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
  } else if (do_all == 7) {  //DAN: worst7 - pm, p2mm, c2mm, p4mm, p4gm, p3m1, p6mm
	for (study_wallpaper_group=2; study_wallpaper_group<3; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=5; study_wallpaper_group<6; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=8; study_wallpaper_group<9; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=10; study_wallpaper_group<12; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=13; study_wallpaper_group<14; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
	for (study_wallpaper_group=16; study_wallpaper_group<17; study_wallpaper_group++) {
	  printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	  uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
	}
  } else if (do_all == 8) { //CLARE: just the one s group
    for (study_wallpaper_group=17; study_wallpaper_group<18; study_wallpaper_group++) {
      printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
      uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]);
    }
  } else {
    printf("find the best packing of %s in group %s with %d sites occupied\n", chosen_shape->name, wallpapergroups[study_wallpaper_group].label, numoccsites);
	uniform_best_packing_in_group(chosen_shape, wallpapergroups[study_wallpaper_group]); 
  }

  //printf("number of mirrors = %d\n", chosen_shape->mirrors);
  fclose(fpout);

  //printf("radius of unit circle at 180 = %f\n", unitcircle.r[180]);
  return 0;
}

