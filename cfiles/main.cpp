#include <iostream>
#include <stdio.h>
#include "mainrand.h"

using namespace std;

extern "C" struct molgeometry /* Nx3 */
{
    int nconf;
    double arr[3*N];
    char *symbols[N];
};

extern "C" void relaxgeom(molgeometry * mg)
{
    xyz_arr *conf = new xyz_arr();
    conf->loadgeom(mg->arr);
    cout << "The error is " << conf->geterror() << endl;
    double vmax = 0;
    double lmin = 2.5;
    conf->gd(vmax,lmin);
    while ((conf->getStericEnergy_clean()>10)&&(vmax<10)){
        conf->gd(vmax,lmin);
        cout << "[!!!!!!] My steric energy = " << conf->getStericEnergy_clean() << endl;
        conf->writeout(mg->symbols);
        vmax += 0.1;
        lmin += 0.01;
    }
    vmax = 5;
    lmin = 2.5;
    for (int i = 0 ; i < N; ++i)
        cout << i << " = " << mg->symbols[i] << endl;
    //conf->gd(1,2.5);
    //conf->gd(10,4);
    //conf->gd(1,2.5);
    //conf->gd();
    cout << "Here alive 2"<< endl;
    conf->getgeom(mg->arr);
    cout << "Here alive 3"<< endl;
}
