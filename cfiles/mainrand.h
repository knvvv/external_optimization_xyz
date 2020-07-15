#ifndef MAIN_RAND_H
#define MAIN_RAND_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <sys/time.h>
#include <math.h>

#define N 41
#define MAXE 0.5
#define MINDSQ 16  // SQUARED !!!
#define STERCOEFF_E (double)MAXE/MINDSQ/MINDSQ
#define STERCOEFF_G (double)4*MAXE/MINDSQ/MINDSQ
#define NPAIRS 702

using namespace std;
class xyz_arr {
    public:
	// See main_rand.cpp
	xyz_arr();
	~xyz_arr();
	bool gd(double,double);
	bool makeamove(gsl_vector*,double);
	void print_vector(string,gsl_vector*);
	void writeout(char *sym[]);
	void loadgeom(double*);
	void getgeom(double*);
    double getStericEnergy();
    double getStericEnergy_clean();
    void updateStericGrad();

	// cycldef.cpp
	void updgrad();
	void initSteric();
	double geterror();

    private:
        int ** distatoms;
        int confcount;
        double maxSE, minD,minDSQ,stercoeff_e,stercoeff_g;

        double * coord; // x0 y0 z0 x1 y1 z1 ...
        double * mygrad;
        double * coord_prev;
        double * grad_prev;
        double * hessp;
        string * atsymbols;

        gsl_eigen_symmv_workspace * ews;
        gsl_rng * r;

        gsl_matrix *mhess, *evec;
        gsl_vector *vgrad, *vcoord, *vgrad_prev, *vcoord_prev, *eval;
};

#endif