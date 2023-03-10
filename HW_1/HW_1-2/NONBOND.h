//
//  NONBOND.h
//  HW_1
//
//  Created by 陳啟瑋 on 2023/3/10.
//

#ifndef NONBOND_h
#define NONBOND_h

#include <iostream>
#include <cmath>
#include "ATOM.h"

using namespace std;

class NONBOND{
    public:
        NONBOND(){
            r=r2=vdw=cou=eng=0;
        };
        double r; //separation distance
        double r2; //r*r
        double vdw; //van der waals energy in kcal/mol
        double cou; //coulomb energy in kcal/mol
        double eng; //nonbond energy in kcal/mol
        double cal_r(ATOM &,ATOM &); //calculate r and r2
        double cal_eng(ATOM &,ATOM &); //calculate nonbond energy
};

double NONBOND::cal_r(ATOM &a,ATOM &b){
    r2 = (a.x[0]-b.x[0])*(a.x[0]-b.x[0])
        +(a.x[1]-b.x[1])*(a.x[1]-b.x[1])
        +(a.x[2]-b.x[2])*(a.x[2]-b.x[2]);
    r = sqrt(r2);
    return r;
}

double NONBOND::cal_eng(ATOM &a,ATOM &b){
    
    cal_r(a, b); // calc r2 and r
    cou = 332.0637 * a.chg * b.chg/r;
    
    double D0, R0 = 0.0;
    
    D0 = sqrt(a.Do * b.Do);
    R0 = (a.Ro + b.Ro)*0.5;
    vdw = D0*(pow(R0/r, 12) - 2.0*pow(R0/r, 6));
    eng = vdw + cou;
    return eng;
}

#endif /* NONBOND_h */
