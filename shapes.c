#include "globals.h"
double shape_var=1;

int hcf(int u, int v) {
  int t;
  if (v>u) {
    t=u;
    u=v;
    v=t;
  }
  while (v!=0) {
    t=v;
    v=u%v;
    u=t;
  }
  return u;
}

/* DAN */
void define_fourier(struct shapetype *shape, double fcoeff[2][FOURIER_TERMS]) {
  int i, j, p, q;
  double phi;
  sprintf(shape->name, "%s", "fourier");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;   
    shape->r[i] = 1.0;
    for (j=0; j<FOURIER_TERMS; j++) {
      shape->r[i] += fcoeff[0][j]*cos((j+1)*phi) + fcoeff[1][j]*sin((j+1)*phi);
    }
    if (shape->r[i] <= 0) {
      fprintf(stderr, "shape has 0 or negative radius at phi=%f\n", phi);
      exit(1);
    }
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }

  //find rotsym. if all coeffs are zero, rotsym is number of points. if any 1 coeffs are nonzero, rotsym is 1. otherwise, proceed through all other coeffs, find hcf.
  if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][2]==0) && (fcoeff[0][3]==0) && (fcoeff[0][4]==0) && (fcoeff[0][5]==0) && (fcoeff[0][6]==0) && (fcoeff[0][7]==0) && (fcoeff[0][8]==0) && (fcoeff[1][0]==0) && (fcoeff[1][1]==0) && (fcoeff[1][2]==0) && (fcoeff[1][3]==0) && (fcoeff[1][4]==0) && (fcoeff[1][5]==0) && (fcoeff[1][6]==0) && (fcoeff[1][7]==0) && (fcoeff[1][8]==0) ) {
    shape->rotsym=SHAPE_RESOLUTION;
  } else if ((fcoeff[0][0]!=0) || (fcoeff[1][0]!=0)) {
    shape->rotsym=1;
  } else {
    q = 0;
    for (i=0; i<2; i++) {
      for (j=1; j<FOURIER_TERMS; j++) {
	if (fcoeff[i][j]!=0) {
	  p=j+1;
	  if (q == 0) {
	    q = p;
	  }
	  q=hcf(p,q);
	}
      }
    }
    shape->rotsym=q;
  }

  //if there are any sine terms, then there is not a mirror on the x-axis. otherwise, there is.
  //BUT on 2013-04-02 - we realised that although this statement is strictly true, (i) shapes with sine terms can have mirrors elsewhere which we must consider, and (ii) shapes with cos terms can have more than one mirror. So we temporarily patched this part of code up, so now it works (a) for non-zero coefficients only from 2 to 6, and (b) when you input the coefficients as all cos terms if there was originally a sine term and the shape had one or more mirrors. (i.e. input the same shape using different coefficient values, by expressing shape with only cos's, effectively rotating the shape so it's mirror is on the x axis). mirrors after the first mirror on the x axis must be in certain positions (e.g. perpendicular if there are 2 mirrors).
  //so first: if any (non-zero) sine terms, no mirrors.
  if ( (fcoeff[1][0]!=0) || (fcoeff[1][1]!=0) || (fcoeff[1][2]!=0) || (fcoeff[1][3]!=0) || (fcoeff[1][4]!=0) || (fcoeff[1][5]!=0) || (fcoeff[1][6]!=0) || (fcoeff[1][7]!=0) || (fcoeff[1][8]!=0) ) {
    shape->mirrors = 0;
  } else {
    //so now no sine terms. next, if shape is a single (i.e. only one term, a cos), then same number of mirrors as order of cos term. NB i have ignored the possibility of the single cos term being first order, since i never use these.
    if ( (fcoeff[0][0]==0) && (fcoeff[0][2]==0) && (fcoeff[0][3]==0) && (fcoeff[0][4]==0) && (fcoeff[0][5]==0) && (fcoeff[0][6]==0) && (fcoeff[0][7]==0) && (fcoeff[0][8]==0) ) {
      shape->mirrors = 2;
    } else if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][3]==0) && (fcoeff[0][4]==0) && (fcoeff[0][5]==0) && (fcoeff[0][6]==0) && (fcoeff[0][7]==0) && (fcoeff[0][8]==0) ) {
      shape->mirrors = 3;
    } else if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][2]==0) && (fcoeff[0][4]==0) && (fcoeff[0][5]==0) && (fcoeff[0][6]==0) && (fcoeff[0][7]==0) && (fcoeff[0][8]==0) ) {
      shape->mirrors = 4;
    } else if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][2]==0) && (fcoeff[0][3]==0) && (fcoeff[0][5]==0) && (fcoeff[0][6]==0) && (fcoeff[0][7]==0) && (fcoeff[0][8]==0) ) {
      shape->mirrors = 5;
    } else if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][2]==0) && (fcoeff[0][3]==0) && (fcoeff[0][4]==0) && (fcoeff[0][6]==0) && (fcoeff[0][7]==0) && (fcoeff[0][8]==0) ) {
      shape->mirrors = 6;
    } else if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][2]==0) && (fcoeff[0][3]==0) && (fcoeff[0][4]==0) && (fcoeff[0][5]==0) && (fcoeff[0][7]==0) && (fcoeff[0][8]==0) ) {
      shape->mirrors = 7;
    } else if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][2]==0) && (fcoeff[0][3]==0) && (fcoeff[0][4]==0) && (fcoeff[0][5]==0) && (fcoeff[0][6]==0) && (fcoeff[0][8]==0) ) {
      shape->mirrors = 8;
    } else if ( (fcoeff[0][0]==0) && (fcoeff[0][1]==0) && (fcoeff[0][2]==0) && (fcoeff[0][3]==0) && (fcoeff[0][4]==0) && (fcoeff[0][5]==0) && (fcoeff[0][6]==0) && (fcoeff[0][7]==0) ) {
      shape->mirrors = 9;
    } else if ( (fcoeff[0][0]!=0) || (fcoeff[0][4]!=0) || (fcoeff[0][6]!=0) || (fcoeff[0][8]!=0) ) {
      //so now no sine terms and not a single. next, if any (non-zero) cos order 1, 5, 7 or 9 (i.e. all odds except 3), there is one mirror, on the x axis.
      shape->mirrors = 1;
    } else if ( (fcoeff[0][2]!=0) ) {
      //so now, no sine terms and not a single and no cos terms of order 1, 5, 7, or 9. next, if cos order 3 is non-zero, there are three mirrors if only other term is cos 6, otherwise, 1 mirror.
      if ( (fcoeff[0][5]!=0) && (fcoeff[0][1]==0) && (fcoeff[0][3]==0) && (fcoeff[0][7]==0) ) {
	shape->mirrors = 3;
      } else {
	shape->mirrors = 1;
      }
    } else if ( (fcoeff[0][0]==0) || (fcoeff[0][2]==0) || (fcoeff[0][4]==0) || (fcoeff[0][6]==0) || (fcoeff[0][8]==0) ) {
      //so now, no sine terms and not a single and no cos terms of odd order. then there are 2 mirrors.
      shape->mirrors = 2;
    }
  }
}

void define_circle(struct shapetype *c) {
  int i;
  sprintf(c->name, "%s", "circle");
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    c->r[i] = 1.0;
  }
  c->maxr = 1.0;
  c->minr = 1.0;
  c->rotsym = 360;  /*shouldn't that be in terms of SHAPE_RESOLUTION? */
  c->mirrors = 360;
}

/* DAN shape_var = 0 gives a circle, shape_var 1 gives a smoothed gear. */
void define_wobblycircle(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "wobblycircle");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    shape->r[i] = 0.8+0.2*shape_var*sin(10*phi);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i];
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i];
  }
  shape->rotsym = 10;
  shape->mirrors = 10;
}

void define_pointyegg(struct shapetype *egg) {
  int i;
  double phi;
  sprintf(egg->name, "%s", "pointyegg");
  egg->maxr=0.0;
  egg->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;   
    phi = phi/4;
    egg->r[i] = fabs(sin(M_PI/4.0)*1.0/(sin(3.0*M_PI/4.0-phi)));
    egg->r[i] *= egg->r[i];
    if (egg->r[i]>egg->maxr) egg->maxr = egg->r[i]; 
    if (egg->r[i]<egg->minr) egg->minr = egg->r[i]; 
  }
  egg->rotsym = 1;
  egg->mirrors = 1;
}

/* DAN shape_var = 1 gives a circle, shape_var = 0 gives a straight line. eccentricity is sqrt((0.25-pow((0.5*shape_var),2)/0.25) */
void define_ellipse(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "ellipse");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;   
    shape->r[i] = (0.5*(0.5*shape_var))/sqrt(pow((0.5*cos(phi)),2)+pow(((0.5*shape_var)*sin(phi)),2));
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 2;
  shape->mirrors = 2;
}

/* DAN */
void define_wobbly4(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "wobbly4");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;   
    shape->r[i] = 0.6+0.4*cos(4*phi);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 4;
  shape->mirrors = 4;
}

void define_wobbly3(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "wobbly3");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;   
    shape->r[i] = 0.6+0.4*cos(3*phi);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 3;
  shape->mirrors = 3;
}

void define_wobbly5(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "wobbly5");
  shape->maxr=0.0;
  shape->maxr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;   
    shape->r[i] = 0.6+0.4*cos(5*phi);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 5;
  shape->mirrors = 5;
}

void define_wobblyfan3(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "wobblyfan3");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;   
    shape->r[i] = 0.2+0.6*fabs(sin(3*phi/2.0))*fmod(3.0*phi/(2.0*M_PI),1.0);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 3;
  shape->mirrors = 0;
}

void define_ninjastar(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "ninjastar");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    while (phi>2.0*M_PI/5.0) phi -= 2.0*M_PI/5.0;
    shape->r[i] = sin(2.0*M_PI/10.0)/sin(M_PI-2.0*M_PI/10.0-phi);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 5;
  shape->mirrors = 0;
}

void define_pentagon(struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "pentagon");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    while (phi>2.0*M_PI/5.0) phi -= 2.0*M_PI/5.0;
    shape->r[i] = sin(0.5*M_PI*(5-2)/5.0)/sin(M_PI-0.5*M_PI*(5-2)/5.0-phi);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 5;
  shape->mirrors = 5;
}

/* DAN - see notebook*/
void define_pentagoncircle(struct shapetype *shape) {
  int i;
  double phi;
  double beta;                     //angle of one segment of pentagon
  double alpha;
  double theta;
  double gamma;
  double d_centres;
  double theta_dash;
  sprintf(shape->name, "%s", "pentagoncircle");
  shape->maxr=0.0;
  shape->minr=100.0;
  beta = 2*M_PI/5;
  alpha = M_PI-beta/2;
  theta = 2*asin(sin(alpha)/shape_var);
  gamma = M_PI-alpha-theta/2;
  d_centres = sqrt(pow(shape_var,2)+1-2*shape_var*cos(gamma));
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    while (phi>=beta) phi -= beta;
    theta_dash = M_PI-phi-alpha-asin(d_centres*sin(phi+alpha)/shape_var);
    shape->r[i] = sqrt(pow(shape_var,2)+pow(d_centres,2)-2*shape_var*d_centres*cos(theta_dash));
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 5;
  shape->mirrors = 5;
}

/* DAN - see notebook*/
void define_heptagoncircle(struct shapetype *shape) {
  int i;
  double phi;
  double beta;                     //angle of one segment of heptagon
  double alpha;
  double theta;
  double gamma;
  double d_centres;
  double theta_dash;
  sprintf(shape->name, "%s", "heptagoncircle");
  shape->maxr=0.0;
  shape->minr=100.0;
  beta = 2*M_PI/7;
  alpha = M_PI-beta/2;
  theta = 2*asin(sin(alpha)/shape_var);
  gamma = M_PI-alpha-theta/2;
  d_centres = sqrt(pow(shape_var,2)+1-2*shape_var*cos(gamma));
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    while (phi>=beta) phi -= beta;
    theta_dash = M_PI-phi-alpha-asin(d_centres*sin(phi+alpha)/shape_var);
    shape->r[i] = sqrt(pow(shape_var,2)+pow(d_centres,2)-2*shape_var*d_centres*cos(theta_dash));
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 7;
  shape->mirrors = 7;
}

void define_RichardJames_pentagon(struct shapetype *shape, double angle_D) {
  int i;
  double phi;
  double theta_1, theta_2;
  double xb, yb, xc, yc;
  double x_int, y_int;
  double incline;
  double x1;
  double side_d;


  side_d = (1.0 + cos(angle_D)+(tan(angle_D/2.0))*(sin(angle_D)-1))/(sin(angle_D)+cos(angle_D)+tan(angle_D/2.0)*(sin(angle_D)-cos(angle_D)));


  sprintf(shape->name, "%s_%.3f", "RJpentagon", angle_D);
  shape->maxr=0.0;
  shape->minr=100.0;


  xc = side_d*cos(angle_D-(M_PI/2.0));
  yc = 1.0 + side_d*sin(angle_D-(M_PI/2.0));

  xb = 1.0+(1.0-side_d)*cos(angle_D);
  yb = (1.0-side_d)*sin(angle_D);

  theta_1 = atan((yb-0.5)/(xb-0.5));
  theta_2 = atan((xc-0.5)/(yc-0.5));

  printf("theta %f theta %f\n", theta_1, theta_2);

  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;

	// using the origin at (0.5,0.5) from E
	if ((phi>5.0*M_PI/4.0)&&(phi<7.0*M_PI/4.0)) {
	  // this is segment AE
	  shape->r[i] = 0.5/cos(phi-(3.0*M_PI/2.0)); 
	} else if ((phi>3.0*M_PI/4.0)&&(phi<5.0*M_PI/4.0)) {
	  // this is segment DE
	  shape->r[i] = 0.5/cos(phi-M_PI); 
	} else if (((phi>7.0*M_PI/4.0)&&(phi<2.0*M_PI+theta_1))||(phi<theta_1)) {
	  // this is segment AB
	  x1 = 1.0;
	  incline = angle_D;
	  x_int = (x1*tan(incline)+0.5-0.5*tan(phi))/(tan(incline)-tan(phi));
	  y_int = (x_int-x1)*tan(incline);
	  shape->r[i] = sqrt(pow(x_int-0.5,2.0)+pow(y_int-0.5,2.0));
	} else if ((phi>(M_PI/2.0-theta_2))&&(phi<3.0*M_PI/4.0)) {
	  // this is segment CD
	  incline = angle_D-(M_PI/2.0);
	  x_int = (1.0+tan(phi))/(2.0*(tan(phi)-tan(incline)));
	  y_int = 1+x_int*tan(incline);
	  shape->r[i] = sqrt(pow(x_int-0.5,2.0)+pow(y_int-0.5,2.0));	  
	} else {
	  // this is segment BC
	  incline = (angle_D-M_PI)/2.0;
	  x1 = side_d*cos(angle_D-(M_PI/2.0))-((1.0+side_d*sin(angle_D-(M_PI/2.0)))/tan((angle_D-M_PI)/2.0));
	  x_int = (x1*tan(incline)+0.5-0.5*tan(phi))/(tan(incline)-tan(phi));
	  y_int = (x_int-x1)*tan(incline);
	  shape->r[i] = sqrt(pow(x_int-0.5,2.0)+pow(y_int-0.5,2.0));	  
	}

	if (SHAPE_RESOLUTION%8==0) {
	  // fix exact corners
	  shape->r[SHAPE_RESOLUTION*3/8] = sqrt(0.5);
	  shape->r[SHAPE_RESOLUTION*5/8] = sqrt(0.5);
	  shape->r[SHAPE_RESOLUTION*7/8] = sqrt(0.5);
	}

    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 1;
  shape->mirrors = 0;
}

void define_MarjorieRice_pentagon(struct shapetype *shape, double angle_B) {
  int i;
  double phi;
  double theta_b, theta_c, theta_d;
  double xd, yd, xc, yc;
  double x_int, y_int;
  double incline;
  double x1;
  double side_b, side_c;


  side_c = 1.0/(sqrt(2.0)*sin(angle_B-M_PI/4.0)+2.0*cos(angle_B));
  side_b = 2.0*side_c*sin(angle_B)+sqrt(2.0)*side_c*cos(angle_B-M_PI/4.0);

  sprintf(shape->name, "%s_%.3f", "MRpentagon", angle_B);
  shape->maxr=0.0;
  shape->minr=100.0;


  xc = side_b-side_c*cos(angle_B);
  yc = side_c*sin(angle_B);

  xd = 2.0*side_c*cos(angle_B-M_PI/2.0);
  yd = 1.0+2.0*side_c*sin(angle_B-M_PI/2.0);

  // using the origin at (0.5,0.5) from A
  theta_b = atan((-0.5)/(side_b-0.5));
  theta_c = atan((yc-0.5)/(xc-0.5));
  theta_d = atan((yd-0.5)/(xd-0.5));
 
  if (theta_b<0.0) theta_b += 2.0*M_PI;
  if (theta_c<0.0) theta_c += 2.0*M_PI;
  if (theta_d<0.0) theta_d += 2.0*M_PI;

  printf("MR pentagon theta_b %f theta_c %f theta_d %f side_b %f side_c %f xc %f yc %f xd %f yd %f\n", theta_b, theta_c, theta_d, side_b, side_c, xc, yc, xd, yd);

  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;

	// using the origin at (0.5,0.5) from A
	if ((phi>5.0*M_PI/4.0)&&(phi<theta_b)) {
	  // this is segment AB
	  shape->r[i] = 0.5/cos(phi-(3.0*M_PI/2.0)); 
	} else if ((phi>3.0*M_PI/4.0)&&(phi<5.0*M_PI/4.0)) {
	  // this is segment EA
	  shape->r[i] = 0.5/cos(phi-M_PI); 
	} else if ((phi<3.0*M_PI/4.0)&&(phi>theta_d)) {
	  // this is segment DE
	  incline = angle_B-M_PI/2.0;
	  x1 = -1.0/tan(angle_B-M_PI/2.0);
	  x_int = (x1*tan(incline)+0.5-0.5*tan(phi))/(tan(incline)-tan(phi));
	  y_int = (x_int-x1)*tan(incline);
	  shape->r[i] = sqrt(pow(x_int-0.5,2.0)+pow(y_int-0.5,2.0));
	} else if ((phi>theta_c)&&(phi<theta_d)) {
	  // this is segment CD
	  incline = 3.0*M_PI/2.0 - angle_B;
	  x1=(xc-(yc/yd)*xd)/(1.0-yc/yd);
	  //printf("x1 CD %f\n", x1);
	  x_int = (x1*tan(incline)+0.5-0.5*tan(phi))/(tan(incline)-tan(phi));
	  y_int = (x_int-x1)*tan(incline);
	  shape->r[i] = sqrt(pow(x_int-0.5,2.0)+pow(y_int-0.5,2.0));	  
	} else {
	  // this is segment BC
	  incline = M_PI - angle_B;
	  x1 = side_b;
	  x_int = (x1*tan(incline)+0.5-0.5*tan(phi))/(tan(incline)-tan(phi));
	  y_int = (x_int-x1)*tan(incline);
	  shape->r[i] = sqrt(pow(x_int-0.5,2.0)+pow(y_int-0.5,2.0));	  
	}
  }
  
  if (SHAPE_RESOLUTION%8==0) {
	// fix exact corners
	shape->r[SHAPE_RESOLUTION*3/8] = sqrt(0.5);
	shape->r[SHAPE_RESOLUTION*5/8] = sqrt(0.5);
  }
    
  for (i=0; i<SHAPE_RESOLUTION; i++) {
	shape->r[i] *= 0.25; /* to make it a similar size to the unit radius shapes */
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  

  shape->rotsym = 1;
  shape->mirrors = 0;
}

void define_asymcomma(struct shapetype *shape) {
  int i;
  double phi;
  double phi_rescaled;
  sprintf(shape->name, "%s", "asymmetriccomma");
  shape->maxr=0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    phi_rescaled = sqrt(phi/(2.0*M_PI));
    shape->r[i] = 0.4+0.6*pow(2.0*(phi_rescaled-0.5), 4.0);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 1;
  shape->mirrors = 0;
}

void define_square(struct shapetype *sq) {
  int i;
  double phi;
  int rot_symmetry = 4;
  sprintf(sq->name, "%s", "square");
  sq->maxr=0.0;
  sq->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    } 
    sq->r[i] = fabs(sin(M_PI/4.0)*1.0/(sin(3.0*M_PI/4.0-phi)));
    if (sq->r[i]>sq->maxr) sq->maxr = sq->r[i]; 
    if (sq->r[i]<sq->minr) sq->minr = sq->r[i]; 
  }
  sq->rotsym = rot_symmetry;
  sq->mirrors = 4;
}

void define_bowtie(struct shapetype *bowtie) {
  int i;
  double phi;
  int rot_symmetry = 2;
  sprintf(bowtie->name, "%s", "bowtie");
  bowtie->maxr=1.0;                             
  bowtie->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    /* to implement the rotational symmetry */
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    /* to implement the y-mirror symmetry */
    if (phi>M_PI/2.0) {
      phi = M_PI - phi;
    }

    /* now specify all values of the radius for phi between 0 and 90 degrees */
    
    if (phi==0*M_PI) {
      bowtie->r[i] = 0.22;
    } else if (phi<=(M_PI)/6) {
      bowtie->r[i] = (((12*phi)/M_PI)+13)/60;
    } else if (phi<=(5*M_PI)/18) {
      bowtie->r[i] = 1/3-sqrt(((10*M_PI)/36-phi)/320);
    } else if (phi>(5*M_PI)/18 && phi<=(7*M_PI)/18) {
      bowtie->r[i] = 1-(pow((phi-(14*M_PI)/36),2.0))/24;
    } else if (phi>(7*M_PI)/18 && phi<=(M_PI)/2) {
      bowtie->r[i] = 1-(pow((phi-(14*M_PI)/36),2.0))/80;
    } else if (bowtie->r[i]>bowtie->maxr) {
      bowtie->maxr = bowtie->r[i];
    } else if (bowtie->r[i]<bowtie->minr) bowtie->minr = bowtie->r[i]; 
  }
  bowtie->rotsym = rot_symmetry;
  bowtie->mirrors = 2;
}

void define_zigzagbone(struct shapetype *b) {
  int i;
  double phi;
  int rot_symmetry = 2;
  sprintf(b->name, "%s", "zigzagbone");
  b->maxr=1.5;                             
  b->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    /* to implement the rotational symmetry */
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    /* to implement the y-mirror symmetry */
    if (phi>M_PI/2.0) {
      phi = M_PI - phi;
    }

    /* now specify all values of the radius for phi between 0 and 90 degrees */
    
    if (phi>=0*M_PI && phi<M_PI/4) {
      b->r[i] = fabs(sin(M_PI/4)*1.0/(sin(3.0*M_PI/4-phi)));
    }
    if (phi>=M_PI/4 && phi<=M_PI/2) {
      b->r[i] = 1.5;
    }
    if (b->r[i]<b->minr) b->minr = b->r[i]; 
  } 
  b->rotsym = rot_symmetry;
  b->mirrors = 2;
}

void define_bells(struct shapetype *b3) {
  int i;
  double phi;
  int rot_symmetry = 2;
  sprintf(b3->name, "%s", "bells");
  b3->maxr=1.5;                             
  b3->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    /* to implement the rotational symmetry */
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    /* to implement the y-mirror symmetry */
    if (phi>M_PI/2.0) {
      phi = M_PI - phi;
    }
    
    /* now specify all values of the radius for phi between 0 and 90 degrees */
    
    if (phi>=0*M_PI && phi<M_PI/4) {
      b3->r[i] = fabs(sin(2*phi));
    }
    if (phi>=M_PI/4 && phi<=M_PI/2) {
      b3->r[i] = 1.5;
    }
    if (b3->r[i]<b3->minr) b3->minr = b3->r[i]; 
  } 
  b3->rotsym = rot_symmetry;
  b3->mirrors = 2;
}

void define_bone(struct shapetype *b4) {
  int i;
  double phi;
  int rot_symmetry = 2;
  sprintf(b4->name, "%s", "bone");
  b4->maxr=1.5;                             
  b4->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    /* to implement the rotational symmetry */
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    /* to implement the y-mirror symmetry */
    if (phi>M_PI/2.0) {
      phi = M_PI - phi;
    }

    /* now specify all values of the radius for phi between 0 and 90 degrees */
    
    if (phi>=0*M_PI && phi<M_PI/4) {
      b4->r[i] = 0.5/sin(M_PI/2-phi);
    }
    if (phi>=M_PI/4 && phi<=M_PI/2) {
      b4->r[i] = 1.5;
    }
    if (b4->r[i]<b4->minr) b4->minr = b4->r[i]; 
  } 
  b4->rotsym = rot_symmetry;
  b4->mirrors = 2;
}

void define_triangle (struct shapetype *triangle) {
  int i;
  double phi;
  int rot_symmetry = 3; 
  sprintf(triangle->name, "%s", "triangle");
  triangle->maxr=sqrt(3);
  triangle->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    
    /* to implement the rotational symmetry */
    while (phi>2.0*M_PI/rot_symmetry) {  /* && phi != 2*M_PI/3 */
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //to implement the x-mirror symmetry
    while (phi>M_PI) {  /* && phi != 2*M_PI/3 && phi != 5*M_PI/3 && phi != 7*M_PI/6 */
      phi = 2*M_PI - phi;
    }
    
    //specifying all values of the radius for phi between -90 and 90 degrees
    if (phi >= 0*M_PI && phi < M_PI/2 && phi != 5*M_PI/6) {
      triangle->r[i] = fabs(sin(M_PI/6)*1.0/(sin(5*M_PI/6-phi)));
    }
    else if (phi >= M_PI/2 && phi < 2*M_PI/3 && phi != M_PI/6) {
      triangle->r[i] = fabs(sin(M_PI/6)*1.0/(sin(M_PI/6+phi)));
    }
    else if (phi >= 2*M_PI/3 && phi < 2*M_PI && phi != 3*M_PI/2) {
      triangle->r[i] = fabs(sin(M_PI/6)*1.0/(sin(3*M_PI/2-phi)));
    }
    if (triangle->r[i]<triangle->minr) triangle->minr = triangle->r[i]; 
  }
  triangle->rotsym = rot_symmetry;
  triangle->mirrors = 3;
} 
 
void define_smoothedpentagon (struct shapetype *shape) { /* NOT WORKING- needs to be fixed*/
  int i;
  double phi;
  double angle2 = M_PI/5;
  double A = pow(tan(angle2),2);
  //double B = pow(tan(3*M_PI/10),2);
  //double C = (2.6*B+sqrt(6.76*pow(B,2)+(A-B)*(6.76-5.6*A)))/(2*B-2*A);
  //double D = tan(3*M_PI/10)*(C-1.3);
  //double E = sqrt(pow(3.3-C,2)+pow(D,2));
  //double angle1 = asin(D/E);
  double angle1 = 0.05; 
  int rot_symmetry = 5;
  sprintf(shape->name, "%s", "smoothedpentagon");
  shape->maxr=2.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    
    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry) { 
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //to implement the x-mirror symmetry
    while (phi>M_PI && phi != M_PI/2) { 
      phi = 2*M_PI - phi;
    }
    
    //specifying all values of the radius for phi between 0 and 72 degrees
    if (phi >= 0*M_PI && phi <= 0*M_PI+angle1) {
      double a = A*pow(cos(phi),2)-pow(sin(phi),2);
      double b = (-6.6)*A*cos(phi);
      double c = 8.93*A;
      shape->r[i] = (-b-sqrt(pow(b,2)-4*a*c))/(2*a);
    }
    else if (phi >= 0*M_PI+angle1 && phi <= 2*M_PI/5-angle1) {
      shape->r[i] = 2.0*sin(3*M_PI/10)/sin(7*M_PI/10-phi);
    }
    else if (phi >= 2*M_PI-angle1 && phi <= 2*M_PI/5) {
      double a = A*pow(cos(2*M_PI/5-phi),2)-pow(sin(2*M_PI/5-phi),2);
      double b = (-6.6)*A*cos(2*M_PI/5-phi);
      double c = 8.93*A;
      shape->r[i] = (-b-sqrt(pow(b,2)-4*a*c))/(2*a);
    }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 5;
}

//DAN
void define_hexagon (struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "hexagon");
  shape->maxr=1;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;

    // to implement the rotational symmetry
    while (phi>=2.0*M_PI/6) {
      phi -= 2.0*M_PI/6;
    }
    
    shape->r[i]=sqrt(3)/(sin(phi)+sqrt(3)*cos(phi));

    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 6;
  shape->mirrors = 6;
}

/* before DAN - incorrect!
void define_hexagon (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 6;
  sprintf(shape->name, "%s", "hexagon");
  shape->maxr=2;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;

    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry && phi!=2.0*M_PI/3 && phi!=M_PI) {   //is that correct?
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //to implement y-mirror symmetry
    if (phi>M_PI/2.0 && phi!=2.0*M_PI/3 && phi!=M_PI) {
      phi = M_PI - phi;
    }

    //specifying all values of the radius for phi between 0 and 90 degrees
     if (phi >= 0*M_PI && phi < M_PI/3 && phi!=2.0*M_PI/3) {
       shape->r[i] = fabs(sin(M_PI/3)*2.0/(sin(2*M_PI/3-phi)));
     }
     else if (phi >= M_PI/3 && phi < M_PI/2 && phi!=M_PI) {
       shape->r[i] = fabs(sin(M_PI/3)*2.0/(sin(M_PI-phi)));
     }
     else if (phi == 2.0*M_PI/3) {
       shape->r[i] = 2;
	 }
     else {
       printf("error phi is %f\n", phi);
     }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 6;
}
*/

void define_heptagon (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 7;
  sprintf(shape->name, "%s", "heptagon");
  shape->maxr=1.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    
    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry) { 
      phi -= 2.0*M_PI/rot_symmetry;
    };
    
    //specifying all values of the radius for phi between 0 and 2pi/7
    if (phi >= 0*M_PI && phi <= 2*M_PI/7) {
      shape->r[i] = sin(5*M_PI/14)/(sin(9*M_PI/14-phi));
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
    }
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 7;
}

void define_nonagon (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 9;
  sprintf(shape->name, "%s", "nonagon");
  shape->maxr=1.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    
    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry) { 
      phi -= 2.0*M_PI/rot_symmetry;
    };
    
    //specifying all values of the radius for phi between 0 and 2pi/9
    //if (phi >= 0*M_PI && phi <= 2*M_PI/9) {
      shape->r[i] = sin(7*M_PI/18)/(sin(11*M_PI/18-phi));
	  //}
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 9;
}

void define_octagon (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 8;
  sprintf(shape->name, "%s", "octagon");
  shape->maxr=1.5;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    
    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry) { 
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //to implement the x-mirror symmetry
    while (phi>M_PI && phi != M_PI/2) { 
      phi = 2*M_PI - phi;
    }
    
    //specifying all values of the radius for phi between 0 and 45 degrees
    if (phi >= 0*M_PI && phi <= M_PI/8) {
      shape->r[i] = sin(M_PI/2)/(sin(M_PI/2-phi));
    }
    else if (phi >= M_PI/8 && phi <= M_PI/4) {
      shape->r[i] = sin(M_PI/2)/(sin(M_PI/4+phi));
    }
    /*  else {
      printf("errorrrr phi is %f\n", phi);
      }*/
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 8;
}

/*void define_randomsmoothedoctagon (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 8;
  sprintf(shape->name, "%s", "randomsmoothedoctagon");
  shape->maxr=1.5;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    
    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry) { 
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //to implement the x-mirror symmetry
    while (phi>M_PI && phi != M_PI/2) { 
      phi = 2*M_PI - phi;
    }
    
    //specifying all values of the radius for phi between 0 and 45 degrees
    if (phi >= 0*M_PI && phi <= 3*M_PI/32) {
      shape->r[i] = sin(M_PI/2)/(sin(M_PI/2-phi));
    }
    else if (phi >= 3*M_PI/32 && phi <= 5*M_PI/32) {
      shape->r[i] =  (sin(M_PI/2)/sin(13*M_PI/32))+(0.005*sin(16*phi-3*M_PI/2));
      }*/
//     else if (phi >= M_PI/16 && phi <= 3*M_PI/16) {
//       shape->r[i] =  (sin(M_PI/2)/sin(7*M_PI/16)) ;
//       }
/*    else if (phi >= 5*M_PI/32 && phi <= M_PI/4) {
      shape->r[i] = sin(M_PI/2)/(sin(M_PI/4+phi));
      }*/
//  else {
//      printf("errorrrr phi is %f\n", phi);
//      }
/*    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 8;
  }*/

void define_smoothedoctagon (struct shapetype *shape) {
  int i;
  double phi;
  //see notes to understand where these constants/values come from
  double m = sqrt(6*sqrt(2)-8)*(sqrt(2)+1)/2;
  //double m2 = pow(m,2);
  double l = sqrt(2)-1;
  //double l2 = pow(l,2);
  double d1 = sqrt(pow(((2*sqrt(2)+2)/(sqrt(2)+1-l))-2,2)+pow((2*l*(sqrt(2)+1))/(sqrt(2)+1-l),2));
  //double d1 = 1.0823922;
  double d2 = (d1*sin(3*M_PI/8))/sin(M_PI/4);
  //double d2 = 1.414213562;
  double n = ((-4)*pow(tan(3*M_PI/8),2)+sqrt(4.7*pow(10,-10)))/(6-4*sqrt(2)-2*pow(tan(3*M_PI/8),2));
  //double n = 0.107771425;
  double d3 = sqrt(pow(n-(d2+2.0),2)+pow((n-2.0)*(sqrt(2)+1),2));
  double d4 = sqrt(pow((n-2.0),2)+pow((n-2.0)*(sqrt(2)+1),2));
  double angle = asin(d4*sin(3*M_PI/8)/d3);
  //double angle = 0.107771425; 
  int rot_symmetry = 8;
  sprintf(shape->name, "%s", "smoothedoctagon");
  shape->maxr = d2;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    
    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry) { 
      phi -= 2.0*M_PI/rot_symmetry;
    };
    
    //to implement the x-mirror symmetry
    while (phi>M_PI) { 
      phi = M_PI - phi;
    }
    
    //specifying all values of the radius for phi between 0 and 45 degrees
    if (phi >= 0*M_PI+angle && phi <= M_PI/4-angle) {
      shape->r[i] = (sin(3*M_PI/8)*d2)/(sin(5*M_PI/8-phi));
    }
    else if (phi >= 0*M_PI && phi <= 0*M_PI+angle) {
      double a = pow(l*cos(phi),2)-pow(sin(phi),2);
      double b = (-2)*(d2+2)*cos(phi)*pow(l,2);
      double c = pow(l*(d2+2),2)-pow(m,2); 
      shape->r[i] = fabs((-b-sqrt(pow(b,2)-4*a*c))/(2*a));
    }
    else if (phi >= M_PI/4-angle && phi <= M_PI/4) {
      double a = pow(l*cos(M_PI/4-phi),2)-pow(sin(M_PI/4-phi),2); 
      double b = (-2)*(d2+2)*cos(M_PI/4-phi)*pow(l,2);
      double c = pow(l*(d2+2),2)-pow(m,2); 
      shape->r[i] = fabs((-b-sqrt(pow(b,2)-4*a*c))/(2*a));
    }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 8;
}  

void define_dodecagon (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 12;
  sprintf(shape->name, "%s", "dodecagon");
  shape->maxr=1.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    
    // to implement the rotational symmetry
    while (phi>2.0*M_PI/rot_symmetry) { 
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //to implement y-mirror symmetry
    if (phi>M_PI/2.0) {
      phi = M_PI - phi;
    }
    
    //specifying all values of the radius for phi between 0 and 30 degrees
     if (phi >= 0*M_PI && phi <= M_PI/6) {
       shape->r[i] = sin(5*M_PI/12)/sin(7*M_PI/12-phi);
     }
     /*  else if (phi >= M_PI/6 && phi < M_PI/3) {
	 shape->r[i] = sin(5*M_PI/12)/sin(3*M_PI/4-phi);
	 }
	 else if (phi >= M_PI/3 && phi < M_PI/2) {
	 shape->r[i] = sin(5*M_PI/12)/sin(11*M_PI/12-phi);
	 }
	 else {
	 printf("errorrrr phi is %f\n", phi);
	 }*/
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 12;
}

void define_fan2 (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 2;
  sprintf(shape->name, "%s", "fan2");
  shape->maxr = 0.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    //to implement rotation symmetry
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //specifying all the values of the radius
    /*if (phi>0*M_PI && phi<M_PI/2) {
      phi = phi/4;
      shape->r[i] = 1.0+(sin(M_PI/4.0)*1.0/(sin(3.0*M_PI/4.0-phi)));
      shape->r[i] *= shape->r[i];
      if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
      }*/
    if (phi>=0*M_PI && phi<=M_PI) {
      phi = phi/6;
      shape->r[i] = fabs(sin(M_PI/6.0)*1.0/(sin(5.0*M_PI/6.0-phi)));
      shape->r[i] *= shape->r[i];
      if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    }   
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 0;
}

void define_fan3 (struct shapetype *shape) {
 int i;
  double phi;
  int rot_symmetry = 3;
  sprintf(shape->name, "%s", "fan3");
  shape->maxr = 2.89;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    //to implement rotation symmetry
    while (phi>=2*M_PI/rot_symmetry) {
      phi -= 2*M_PI/rot_symmetry;
    }
    //specifying all the values of the radius
    if ((phi>=0*M_PI) && (phi<M_PI/3)) {
      shape->r[i] = 0.7+(fabs((sin(M_PI/6))/(sin(5*M_PI/6-phi/6))));
      shape->r[i] *= shape->r[i];
    } else if ((phi>=M_PI/3) && (phi<2*M_PI/3)) {
      shape->r[i] = (sin(2*phi-M_PI/6))+1.15;
    } 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 0;
}

void define_fan4 (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 4;
  sprintf(shape->name, "%s", "fan4");
  shape->maxr = 2.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    //to implement rotation symmetry
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //specifying all the values of the radius
    if (phi>0*M_PI && phi<=M_PI/4) {
      phi = phi/8;
      shape->r[i] = fabs((sin(M_PI/8))/(sin(7*M_PI/8-phi)));
    }
    else if (phi>M_PI/4 && phi<=M_PI/2) {
      shape->r[i] = (sin(M_PI/2-phi))+0.1;
    }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 0;
}

void define_flower5 (struct shapetype *shape) {
  int i;
  double phi; 
  int rot_symmetry = 5;
  sprintf(shape->name, "%s", "flower5");
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    //to implement rotation symmetry
    while (phi>=2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //specifying all the values of the radius
    /*   if (phi>0*M_PI && phi<M_PI/5) {
	 phi_rescaled = fabs(1/(sqrt(M_PI/2-(5*phi/2))));
	 shape->r[i] = fabs(sin(phi_rescaled));
	 }*/
    if (phi>=0*M_PI && phi<2*M_PI/15) {
      shape->r[i] = 0.7+fabs(sin((15*phi)/4));
    } else if ((phi>=2*M_PI/15) && (phi<4*M_PI/15)) {
      shape->r[i] = 0.7+(fabs((sin(M_PI/15))/(sin(14*M_PI/15-phi/15))));
    } else if (phi>=4*M_PI/15 && phi<2*M_PI/5) {
      shape->r[i] = (sin(2*phi-M_PI/10))+0.8;
    }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->maxr = 1.8; //=fabs(sin(fabs(1/(sqrt(M_PI/4-phi)))));
  shape->rotsym = rot_symmetry;
  shape->mirrors = 0;
}

void define_fan6 (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 6;
  sprintf(shape->name, "%s", "fan6");
  shape->maxr = 2.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    //to implement rotation symmetry
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //specifying all the values of the radius
    if (phi>=0*M_PI && phi<=M_PI/6) {
      phi = phi/12;
      shape->r[i] = fabs((sin(M_PI/12))/(sin(11*M_PI/12-phi)))-0.1;
      shape->r[i] *= shape->r[i];
    }
    else if (phi>M_PI/6 && phi<=M_PI/3) {
      shape->r[i] = fabs(sin(M_PI/3-phi))+0.1;
    }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 0;
}

void define_fan12 (struct shapetype *shape) {
  int i;
  double phi;
  int rot_symmetry = 12;
  sprintf(shape->name, "%s", "fan12");
  shape->maxr = 2.0;
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    //to implement rotation symmetry
    while (phi>2.0*M_PI/rot_symmetry) {
      phi -= 2.0*M_PI/rot_symmetry;
    };
    //specifying all the values of the radius
    if (phi>=0*M_PI && phi<=M_PI/6) {
      shape->r[i] = ((sqrt(3)/2)*sin(M_PI/24)/(sin(M_PI/2-phi)))+(0.5*cos(3*phi));
    }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = rot_symmetry;
  shape->mirrors = 0;
}

void define_kite (struct shapetype *shape) {
  int i;
  double phi;
  sprintf(shape->name, "%s", "kite");
  shape->maxr = sqrt(3.0);
  shape->minr=100.0;
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = (2.0*M_PI*i)/SHAPE_RESOLUTION;
    //to implement x-mirror symmetry
    if (phi>M_PI) {
      phi = 2*M_PI - phi;
    }
    //specifying all the values of the radius
    if (phi>=0*M_PI && phi<=M_PI/2) {
      shape->r[i] = sin(M_PI/4)/(sin(3*M_PI/4-phi));
    }
    else if (phi>=M_PI/2 && phi<=M_PI) {
      shape->r[i] = (sqrt(3)*sin(M_PI/6))/(sin(phi-M_PI/6));
    }
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  shape->rotsym = 1;
  shape->mirrors = 2;
}

//ALI
void define_fileshape(struct shapetype *shape) {
  int i;
  double phi;
  double x,y;
  char filename[MAXFILENAME];
  FILE *fp;
  double modcheck;
  sprintf(filename, "fileshape.dat");
    
  sprintf(shape->name, "%s", "fileshape");
  shape->maxr=0.0;
  shape->minr=100.0;
  fp = fopen(filename, "r");
  //check that number of lines in the file is equal to the SHAPE_RESOLUTION
  
  for (i=0; i<SHAPE_RESOLUTION; i++) {
    phi = 2.0*M_PI*i/SHAPE_RESOLUTION;
    //once we get the x,y from the file, we could check that it is also at this phi value

    fscanf(fp,"%lf %lf", &x, &y);

    if ( fabs(modf( 0.5+(phi-atan(y/x))/M_PI, &modcheck)-0.5) > 0.001 ) {
      printf("ERROR: fileshape.dat may have a different number of points to SHAPE_RESOLUTION\n These should be equal: %lf %lf %lf\n " , atan(y/x), phi, modf( (phi-atan(y/x))/M_PI, &modcheck));
      exit(1);
    }
    
    shape->r[i] = pow((x*x+y*y),0.5);
    if (shape->r[i]>shape->maxr) shape->maxr = shape->r[i]; 
    if (shape->r[i]<shape->minr) shape->minr = shape->r[i]; 
  }
  //could go around the shape and check for symmetries, but for now assume no symmetry
  shape->rotsym = 1;
  shape->mirrors = 0;
}
