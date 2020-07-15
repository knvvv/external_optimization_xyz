#include <iostream>
#include <stdio.h>
#include <math.h>
#include "mainrand.h"

void xyz_arr::initSteric(){
    $ATCONST$
}

void xyz_arr::updgrad() {
    $GRAD$
    updateStericGrad();
}

double xyz_arr::geterror() {
    $ERRF$
}