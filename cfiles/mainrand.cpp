#include "mainrand.h"

using namespace std;

xyz_arr::xyz_arr(){
	gsl_rng_env_setup();
	r = gsl_rng_alloc (gsl_rng_ranlxd2);
	struct timeval tv;
	gettimeofday(&tv,0);
	unsigned long mySeed = tv.tv_sec + tv.tv_usec;
	gsl_rng_set(r, mySeed);
	ews = gsl_eigen_symmv_alloc (3*N);
	atsymbols = new string[N];

    confcount = 0;

	vcoord = gsl_vector_alloc(3*N);
	vgrad = gsl_vector_alloc(3*N);
	vcoord_prev = gsl_vector_alloc(3*N);
	vgrad_prev = gsl_vector_alloc(3*N);
	eval = gsl_vector_alloc(3*N);
	mhess = gsl_matrix_alloc (3*N, 3*N);
	evec = gsl_matrix_alloc (3*N, 3*N);

	coord = vcoord->data;
	mygrad = vgrad->data;
	coord_prev = vcoord_prev->data;
	grad_prev = vgrad_prev->data;
	hessp = mhess->data;

    distatoms = new int*[2];
    distatoms[0] = new int[NPAIRS];
    distatoms[1] = new int[NPAIRS];
	initSteric();
}

xyz_arr::~xyz_arr(){
    delete [] atsymbols;
    delete [] distatoms[0];
    delete [] distatoms[1];
    delete [] distatoms;
	gsl_rng_free(r);
	gsl_eigen_symmv_free (ews);
	
	gsl_vector_free(vgrad);
	gsl_vector_free(vcoord);
	gsl_vector_free(vgrad_prev);
	gsl_vector_free(vcoord_prev);
	gsl_vector_free(eval);

	gsl_matrix_free(evec);
	gsl_matrix_free(mhess);
}

bool xyz_arr::gd(double e,double d) {
    maxSE = e;
    minD = d;
    stercoeff_e = maxSE/minD/minD;
    stercoeff_g = 2*maxSE/minD/minD;
    minDSQ = minD*minD;
	int step = 0;
	int sign;
	double stepsize,mynorm,mydot;
	while((step<404)&&(geterror() > 1e-5)){
		updgrad();
		if(step == 0){
			gsl_blas_ddot(vgrad,vgrad,&mynorm);
			stepsize = 1/sqrt(mynorm);
		} else {
			gsl_blas_dscal(-1,vcoord_prev);
			gsl_blas_daxpy(1,vcoord,vcoord_prev);
			gsl_blas_dscal(-1,vgrad_prev);
			gsl_blas_daxpy(1,vgrad,vgrad_prev);
			gsl_blas_ddot(vgrad_prev,vcoord_prev,&mydot);
			gsl_blas_ddot(vgrad_prev,vgrad_prev,&mynorm);
			if (abs(mynorm) < 1e-10 ) {
			    cout << "The norm is " << mynorm << " so gd is finished" << endl;
			    break;
			} else {
			    stepsize = abs(mydot)/mynorm;
			}
		}
		print_vector("Gradient is ",vgrad);
        cout << "Gradnorm is "<<mynorm << endl;
        cout << "Error is "<< geterror() << endl;
		gsl_blas_dcopy(vcoord,vcoord_prev);
		gsl_blas_dcopy(vgrad,vgrad_prev);
		makeamove(vgrad,-stepsize); // MakeMove changes stepsize, not really honest for GD method.
		step++;
	}          

	cout << "Error = " << geterror() << endl;
	cout << "Done in " << step << " steps" << endl;
	cout << "Here alive 1"<< endl;
}

void xyz_arr::writeout(char *sym[]) {
    confcount += 1;
    char filename[100];
    for(int j =0;j<100;j++)
			filename[j]='\0';
	sprintf(filename,"geometries/relaxedconf%d.xyz",confcount);
	FILE *pFile = fopen (filename,"w");
	char newline[200];
	string newstring;
	newstring.clear();

	for(int i=0;i<N;++i){
		for(int j =0;j<200;j++)
			newline[j]='\0';
		sprintf(newline,"%2s %f %f %f",sym[i],coord[3*i],coord[3*i+1],coord[3*i+2]);
		newstring += newline;
		newstring += "\n";
	}
	fprintf(pFile,"%s",newstring.c_str());
	fclose (pFile);
}

void xyz_arr::loadgeom(double *arr){
    memcpy(coord,arr,sizeof(double)*3*N);
}

void xyz_arr::getgeom(double *arr){
    memcpy(arr,coord,sizeof(double)*3*N);
}

bool xyz_arr::makeamove(gsl_vector * vec,double size){
    cout << "Entering makemove" << endl;
	double sterr = geterror();
	double myerr;
	gsl_blas_daxpy(size,vec,vcoord);
	int count = 0;
	double preverr = sterr;
	bool mycase = false;
	myerr = geterror();
	while (myerr < preverr){
		mycase = true;
		size *= 2;
		gsl_blas_daxpy(size/2,vec,vcoord);
		preverr = myerr;
		myerr = geterror();
	}
	
	count = 0;
	while ((myerr>sterr)&&(count < 30)){
		count++;
		gsl_blas_daxpy(-size/pow(2,count),vec,vcoord);
		myerr = geterror();
	}
	
	if(count == 30){
		gsl_blas_daxpy(-size/pow(2,count),vec,vcoord);
	}
	
	if(geterror() - sterr>1e-10) {
		cout << "MakeMove is broken!!" << endl;
		cout << "Increase in error = " << geterror() - sterr << endl;
	}
	return (count<30);
}

void xyz_arr::print_vector(string desc, gsl_vector *vec){
	cout << "\n" << desc << endl;
	cout << "[";
	for(int i = 0; i < 3*N; ++i)
		cout << "  "<<vec->data[i];
	cout << "  ]" << endl;
}

double xyz_arr::getStericEnergy(){
    double res = 0;
    double dssq; // DISTANCE SQUARED
    for(int i = 0; i < NPAIRS; ++i){
        dssq = pow(coord[3*distatoms[0][i]] - coord[3*distatoms[1][i]],2) + pow(coord[3*distatoms[0][i]+1] - coord[3*distatoms[1][i]+1],2)+pow(coord[3*distatoms[0][i]+2] - coord[3*distatoms[1][i]+2],2);
        if (dssq<minDSQ)
            res += stercoeff_e*pow(sqrt(dssq)-minD,2);
    }
    return res;
}

double xyz_arr::getStericEnergy_clean(){
    double res = 0;
    double dssq; // DISTANCE SQUARED
    for(int i = 0; i < NPAIRS; ++i){
        dssq = pow(coord[3*distatoms[0][i]] - coord[3*distatoms[1][i]],2) + pow(coord[3*distatoms[0][i]+1] - coord[3*distatoms[1][i]+1],2)+pow(coord[3*distatoms[0][i]+2] - coord[3*distatoms[1][i]+2],2);
        if (dssq<minDSQ)
            res += pow(sqrt(dssq)-minD,2);
    }
    return res;
}

void xyz_arr::updateStericGrad(){
    //cout << "Entering StericGrad" << endl;
    double dssq,ds,mynorm;
    gsl_blas_ddot(vgrad,vgrad,&mynorm);
    //cout << "Gradiend norm before = " << mynorm;
    gsl_vector* dvec = gsl_vector_alloc(3);
    double * darr = dvec->data;
    for(int i = 0; i < NPAIRS; ++i){
        dssq = pow(coord[3*distatoms[0][i]] - coord[3*distatoms[1][i]],2) + pow(coord[3*distatoms[0][i]+1] - coord[3*distatoms[1][i]+1],2)+pow(coord[3*distatoms[0][i]+2] - coord[3*distatoms[1][i]+2],2);
        if (dssq<MINDSQ){
            ds = sqrt(dssq);
            darr[0] = coord[3*distatoms[0][i]] - coord[3*distatoms[1][i]];
            darr[1] = coord[3*distatoms[0][i]+1] - coord[3*distatoms[1][i]+1];
            darr[2] = coord[3*distatoms[0][i]+2] - coord[3*distatoms[1][i]+2];

            gsl_blas_ddot(dvec,dvec,&mynorm);
            //if(abs(mynorm)>1e-10)
            //    cout << "dvec norm^2 =" << mynorm << endl;

            gsl_vector_scale (dvec, 1/ds);
            gsl_vector_scale (dvec, stercoeff_g*(ds-minD));
            //cout << "Subj1 = " << 1/ds << endl;
            //cout << "Subj2_1 = " << STERCOEFF_G << endl;
            //cout << "Subj2_2 = " << dssq-MINDSQ << endl;
            mygrad[3*distatoms[0][i]] += darr[0];
            mygrad[3*distatoms[0][i]+1] += darr[1];
            mygrad[3*distatoms[0][i]+2] += darr[2];
            mygrad[3*distatoms[1][i]] -= darr[0];
            mygrad[3*distatoms[1][i]+1] -= darr[1];
            mygrad[3*distatoms[1][i]+2] -= darr[2];
        }
    }
    gsl_blas_ddot(vgrad,vgrad,&mynorm);
    cout << "Gradiend norm after = " << mynorm;
    gsl_vector_free(dvec);
}