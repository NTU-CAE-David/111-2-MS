//
//  BOND.h
//  HW_1_v1.0
//
//  Created by 陳啟瑋 on 2023/3/8.
//

#ifndef BOND_h
#define BOND_h

#include <iostream>
#include <cmath>
#include "ATOM.h"

using namespace std;

class BOND{
    public:
        BOND(){
            len=len0=Kb=eng=0;
        };
        double len; //in angstroms
        double len0; //equilib bond length in angstroms
        double Kb; //force constant in kcal/mol A2
        double eng; //bond energy in kcal/mol
        double cal_len(ATOM &a,ATOM &b); //calculate bond length from two vectors
        double cal_eng(ATOM &a,ATOM &b); //calculate bond energy
};

double BOND::cal_len(ATOM &a,ATOM &b){
    
    double r2;
    
    r2 = (a.x[0]-b.x[0])*(a.x[0]-b.x[0])
        +(a.x[1]-b.x[1])*(a.x[1]-b.x[1])
        +(a.x[2]-b.x[2])*(a.x[2]-b.x[2]);
    len = sqrt(r2);
    return len;
}

double BOND::cal_eng(ATOM &a,ATOM &b){
    
    cal_len(a, b);
    eng = 0.5*Kb*(len-len0)*(len-len0);
    return eng;
}

#endif /* BOND_h */
