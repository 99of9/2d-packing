#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a > b ? b : a)

#define MAXFILENAME 200
#define MAXGROUPNAME 50
#define MAXMULTIPLICITY 12
#define MAXWYCKOFFS 9
#define DIM 2
#define NUMWALLPAPERGROUPS 18
#define MAXNUMOCCSITES 3
#define MAXVAR 30
#define MAXREPLICAS 500
#define MAXITER 20
#define MAXGRAPHICSIZE 1000
#define GRAPHICFRAME 100
#define CENTRE_SPOT 1
#define AXES_SHOWN 1

#define EPS 1e-9

#define SHAPE_RESOLUTION 72
#define ALLOWFLIPS 1
#define MAXSTEPS 360000
#define CYCLES 32
#define POLYGON_REP 1

#define EMPTY_OUTPUT 1
#define WARNINGS 0
#define FOURIER_TERMS 9

struct shapetype {
  char name[20];
  double r[SHAPE_RESOLUTION];
  double maxr;
  double minr;
  int rotsym; /* number of rotational symmetries about central point*/
  int mirrors;
};

struct shape {
  struct shapetype *type;
  double x, y;
  double theta; /* rotation from original orientation in radians */
};


struct image_type {
  double coord_coeffs[DIM][DIM+1];
  /* first index is the output coordinate,
   second index is input coordinate (+1 is for the constant that is added) 
   
   for example, to work out the new values:
   x_new = coord_coeffs[0][0]*x_old + coord_coeffs[0][1]*y_old + coord_coeffs[0][2];
   y_new = coord_coeffs[1][0]*x_old + coord_coeffs[1][1]*y_old + coord_coeffs[1][2];
  */

  int rotation_offset; /*angles move clockwise*/
  bool flipped;  /* =1, when flipped using the x-axis as mirror, rotation_offset is then employed if mirror is at y-axis for example */
  int sitemirror0;
  int sitemirror90; /* sitemirror is on the cartesian y-axis */
  int sitemirror45; /*sitemirror is on forward diagonal: x,x */
  int sitemirror135; /*sitemirror is on backward diagonal: x,-x */
  int sitemirror30; /*sitemirror is on diagonal: x,2x */
  int sitemirror60; /*sitemirror is on diagonal: 2x,x */
  int sitemirror330; /*sitemirror is on */
  int sitemirror300; /*sitemirror is on */


};

struct wyckoff_type {
  int multiplicity;
  char letter;
  int somevariability;
  int siterotations;
  int sitemirrors;
  struct image_type image[MAXMULTIPLICITY];
};

struct wallpaper_group_type {
  char label[MAXGROUPNAME];
  bool a_b_equal;
  bool hexagonal;
  bool rectangular; /* includes square, which is also specified with a_b_equal=TRUE */
  int num_symmetries;
  int num_wyckoffs;
  struct wyckoff_type wyckoffs[MAXWYCKOFFS];
};
