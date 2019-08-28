/*****************************************************************************
 *
 *  capillary.c
 *
 *  Compile and link with -lm for the maths library. 
 *
 *  This utility produces an output file suitable for initialising
 *  a capillary structure in Ludwig. 
 *
 *  It is assumed the periodic dimension is z, and that the system
 *  is square in the x-y directions. No 'leakage' in the x-y
 *  directions is ensured by making the last site in each direction
 *  solid.
 *
 *  The various system parameters should be set at comiple time,
 *  and are described below. The output file is always in BINARY
 *  format.
 *
 *  1. Output capaillary structure
 *  Set the required parameters and invoke with no argument
 *      ./a.out
 *   
 *  2. Profiles
 *  If the program is invoked with a single phi output file
 *  argument, e.g.,
 *      ./a.out phi-001000.001-001
 *  a scatter plot of the profile of the interface in the
 *  wetting half of the capillary will be produced. That
 *  is height vs r, the radial distance from the centre.
 *  The output file should match the capillary structure!
 *
 *  $Id: capillary.c,v 1.3 2008-11-26 19:34:32 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistcal Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2008 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This is a copy from ../src/map.h; it would be better to include
 * directly, but that incurs additional dependencies on targetDP.h */

enum map_status {MAP_FLUID, MAP_BOUNDARY, MAP_COLLOID, MAP_STATUS_MAX};

/* SYSTEM SIZE */
/* Set the system size as desired. Clearly, this must match the system
 * set in the main input file for Ludwig. */

const int xmax = 32;
const int ymax = 32;
const int zmax = 1;

/* CROSS SECTION */
/* You can choose a square or circular cross section */

enum {CIRCLE, SQUARE, XWALL, YWALL, ZWALL, XWALL_OBSTACLES, XWALL_BOTTOM};
const int xsection = YWALL;

/*Modify the local geometry of the wall*/

int obstacle_number = 1; /* number of obstacles per wall */
int obstacle_length = 6; /* along the wall direction */
int obstacle_height = 10; /* perpendicular from wall */
int obstacle_depth  = 6; /* perpendicular to length and height */
			 /* NOTE: obstacle_depth == xmax/ymax/zmax 
				  means obstacles don't have a z-boundary */

/* SURFACE CHARGE */

const double sigma = 0.125;

/* FREE ENERGY PARAMETERS */
/* Set the fluid and solid free energy parameters. The fluid parameters
 * must match those used in the main calculation. See Desplat et al.
 * Comp. Phys. Comm. (2001) for details. */

const double kappa1 = 0.01;
const double kappa2 = 0.01;
const double kappa3 = 0.01;
const double alpha = 1.0;
const double H1 = 0.002;
const double H2 = -0.002;
const double H3 = 0.00;
const double C = 0.000;	// Following Desplat et al.

/* WETTING */
/* A section of capillary between z1 and z2 (inclusive) will have
 * wetting property H = H, the remainder H = 0 */

const int z1 = 1;
const int z2 = 36;

/* OUTPUT */
/* You can generate a file with solid/fluid status information only,
 * or one which includes the wetting parameter H or charge Q. */

enum {STATUS_ONLY, STATUS_WITH_H, STATUS_WITH_C_H, STATUS_WITH_SIGMA};
const int output_type = STATUS_WITH_C_H;

/* OUTPUT FILENAME */

const char * filename = "capillary.001-001";

static void profile(const char *);

/*****************************************************************************
 *
 *  main program
 *
 *****************************************************************************/

int main(int argc, char ** argv) {


  char * map_in;
  FILE * fp_orig;
  int i, j, k, n;
  int nsolid = 0;

  double * map_h1; // for wetting coefficient H1
  double * map_h2; // for wetting coefficient H2
  double * map_h3; // for wetting coefficient H3
  double * map_c; // for additional wetting coefficient C

  double * map_sig; // for (surface) charge

  double rc = 0.5*(xmax-2);
  double x0 = 0.5*xmax + 0.5;
  double y0 = 0.5*ymax + 0.5;
  double x, y, r;
  double h, h1, h2, theta,theta1, theta2;
    double aa,ab,ac,ad,ae,af;
  int iobst;
  int obst_start[2*obstacle_number][3];
  int obst_stop[2*obstacle_number][3];
  int gap_length;

  FILE  * WriteFile;	
  char  file[800];

  if (argc == 2) profile(argv[1]);

  if (output_type == STATUS_WITH_C_H) {

    printf("Free energy parameters:\n");
    printf("free energy parameter kappa1 = %f\n", kappa1);
    printf("free energy parameter kappa2 = %f\n", kappa2);
    printf("free energy parameter kappa3 = %f\n", kappa3);
    printf("free energy parameter alpha = %f\n", alpha);
    printf("surface free energy   H1     = %f\n", H1);
    printf("surface free energy   H2     = %f\n", H2);
    printf("surface free energy   H3     = %f\n", H3);
    
      aa = alpha*kappa1 + 4*H1;
      ab = alpha*kappa1 - 4*H1;
      ac = alpha*kappa2 + 4*H2;
      ad = alpha*kappa2 - 4*H2;
      ae = alpha*kappa3 + 4*H3;
      af = alpha*kappa3 - 4*H3;
      
      h = (pow(aa, 1.5) - pow(ab, 1.5))/(2*(kappa1 + kappa2)*pow(alpha*kappa2,0.5)) - (pow(ac, 1.5) - pow(ad, 1.5))/(2*(kappa1 + kappa2)*pow(alpha*kappa1,0.5));
      h1 = (pow(ac, 1.5) - pow(ad, 1.5))/(2*(kappa2 + kappa3)*pow(alpha*kappa3,0.5)) - (pow(ae, 1.5) - pow(af, 1.5))/(2*(kappa2 + kappa3)*pow(alpha*kappa2,0.5));
      h2 = (pow(ae, 1.5) - pow(af, 1.5))/(2*(kappa3 + kappa1)*pow(alpha*kappa1,0.5)) - (pow(aa, 1.5) - pow(ab, 1.5))/(2*(kappa3 + kappa1)*pow(alpha*kappa3,0.5));
      
    printf("dimensionless parameter h=cos(theta)   = %f\n", h);
    printf("dimensionless parameter h1=cos(theta1)   = %f\n", h1);
    printf("dimensionless parameter h2=cos(theta2)   = %f\n", h2);
    theta = acos(h);
    theta1 = acos(h1);
    theta2 = acos(h2);
    printf("contact angle theta         = %f radians\n", theta);
    printf("contact angle theta1         = %f radians\n", theta1);
    printf("contact angle theta2         = %f radians\n", theta2);
    theta = theta*180.0/(4.0*atan(1.0));
    theta1 = theta1*180.0/(4.0*atan(1.0));
    theta2 = theta2*180.0/(4.0*atan(1.0));
    printf("                            = %f degrees\n", theta);
    printf("                            = %f degrees\n", theta1);
    printf("                            = %f degrees\n", theta2);

  }
  
    map_in = (char *) malloc(xmax*ymax*zmax*sizeof(char));
    if (map_in == NULL) exit(-1);
    
    map_h1 = (double *) malloc(xmax*ymax*zmax*sizeof(double));
    if (map_h1 == NULL) exit(-1);
    
    map_h2 = (double *) malloc(xmax*ymax*zmax*sizeof(double));
    if (map_h2 == NULL) exit(-1);
    
    map_h3 = (double *) malloc(xmax*ymax*zmax*sizeof(double));
    if (map_h3 == NULL) exit(-1);
    
    map_c = (double *) malloc(xmax*ymax*zmax*sizeof(double));
    if (map_c == NULL) exit(-1);

   switch (xsection) {


  case YWALL:

    for (i = 0; i < xmax; i++) {
      for (j = 0; j < ymax; j++) {
	for (k = 0; k < zmax; k++) {
	  n = ymax*zmax*i + zmax*j + k;
	  map_in[n] = MAP_FLUID;

	  if (j == 0 || j == ymax-1 || i == 0 || i == xmax-1 ) {
	    map_in[n] = MAP_BOUNDARY;
            if (output_type == STATUS_WITH_C_H) {  map_h1[n] = H1; map_h2[n] = H2; map_h3[n] = H3; map_c[n] = C; }
	    ++nsolid;
	  }
	}
      }
    }
    break;
  default:
    printf("No cross-section!\n");
    /* End switch */
  }

  /* picture */

  printf("\nCross section (%d = fluid, %d = solid)\n", MAP_FLUID, MAP_BOUNDARY);

  k = 0;
  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
	n = ymax*zmax*i + zmax*j + k;
      
	if (map_in[n] == MAP_BOUNDARY) printf(" %d", MAP_BOUNDARY);
	if (map_in[n] == MAP_FLUID)    printf(" %d", MAP_FLUID);
    }
    printf("\n");
  }





  if (output_type == STATUS_WITH_C_H)  {
    sprintf(file,"Configuration_capillary.dat");
    WriteFile=fopen(file,"w");
    fprintf(WriteFile,"#x y z n map H1 H2 H3 C\n");

    for (i = 0; i < xmax; i++) {
      for (j = 0; j < ymax; j++) {
	for (k = 0; k < zmax; k++) {

	n = ymax*zmax*i + zmax*j + k;
      
	if (map_in[n] == MAP_BOUNDARY) { fprintf(WriteFile,"%i %i %i %i %d %f %f %f %f \n", i, j, k, n, MAP_BOUNDARY, map_h1[n], map_h2[n], map_h3[n], map_c[n]); }
	if (map_in[n] == MAP_FLUID)    { fprintf(WriteFile,"%i %i %i %i %d %f %f %f %f \n", i, j, k, n, MAP_FLUID, map_h1[n], map_h2[n], map_h3[n], map_c[n]); }

	}  
      }  
    }
    fclose(WriteFile);
  }


  printf("n = %d nsolid = %d nfluid = %d\n", xmax*ymax*zmax, nsolid,
	 xmax*ymax*zmax - nsolid);

  /* Write new data as char */

  fp_orig = fopen(filename, "w");
  if (fp_orig == NULL) {
    printf("Cant open output\n");
    exit(-1);
  }

  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
      for (k = 0; k < zmax; k++) {
	n = ymax*zmax*i + zmax*j + k;

	fputc(map_in[n], fp_orig);

	if (output_type == STATUS_WITH_C_H) {
	  fwrite(map_c + n, sizeof(double), 1, fp_orig);
	  fwrite(map_h1 + n, sizeof(double), 1, fp_orig);
      fwrite(map_h2 + n, sizeof(double), 1, fp_orig);
      fwrite(map_h3 + n, sizeof(double), 1, fp_orig);
	}


      }
    }
  }

  fclose(fp_orig);

  free(map_in);
  free(map_c);
  free(map_h1);
  free(map_h2);
  free(map_h3);


  return 0;
}

/*****************************************************************************
 *
 *  profile
 *
 *  This attempts to give a profile of the interface as a function
 *  of the radial distance from the centre of the capillary defined
 *  above.
 *
 *  For each position (i, j) we examine the 1-d profile phi(z) and
 *  locate the position of zero by linear interpolation. The results
 *  go to standard output.
 *
 *  Note that the test for the interface assumes the phi = -1 block
 *  is in the middle of the system (see block initialisation in
 *  src/phi_stats.c).
 *
 *****************************************************************************/

static void profile(const char * filename) {

  int ic, jc, kc, index;
  int inside;
  double rc = 0.5*(xmax-2);
  double x0 = 0.5*xmax + 0.5;
  double y0 = 0.5*ymax + 0.5;
  double r, x, y;
  double * phi;
  FILE * fp;

  phi = (double *) malloc(xmax*ymax*zmax*sizeof(double));
  if (phi == NULL) {
    printf("malloc(phi) failed\n");
    exit(-1);
  }

  /* Read the data */

  printf("Reading phi data from %s...\n", filename);

  fp = fopen(filename, "r");
  if (fp == NULL) {
    printf("Failed to open %s\n", filename);
    exit(-1);
  }

  for (ic = 0; ic < xmax; ic++) {
    for (jc = 0; jc < ymax; jc++) {
      for (kc = 0; kc < zmax; kc++) {
	index = ic*zmax*ymax + jc*zmax + kc;
	fread(phi + index, 1, sizeof(double), fp);
      }
    }
  }

  fclose(fp);

  /* Find the interface for each solid location */

  for (ic = 0; ic < xmax; ic++) {
    x = 1.0 + ic - x0;
    for (jc = 0; jc < ymax; jc++) {
      y = 1.0 + jc - y0;

      r = sqrt(x*x + y*y);

      /* Work out whether solid or fluid */
      inside = 0;
      if (xsection == SQUARE) {
	if (ic > 0 && ic < xmax-1 && jc > 0 && jc < ymax-1) inside = 1;
      }
      if (xsection == CIRCLE) {
	if (r <= rc) inside = 1;
      }

      if (inside) {
	/* Examine the profile */
	double h, dh;
	h = -1.0;

	for (kc = z1; kc <= z2; kc++) {
	  index = ic*zmax*ymax + jc*zmax + kc;
	  if (phi[index] > 0.0 && phi[index+1] < 0.0) {
	    /* Linear interpolation to get surface position */
	    dh = phi[index] / (phi[index] - phi[index+1]);
	    h = 1.0 + kc + dh;
	  }
	}
	printf("%f %f\n", r, h);
      }
    }
  }

  free(phi);

  /* Do not return! */
  exit(0);

  return;
}
