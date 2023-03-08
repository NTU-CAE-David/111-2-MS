//
//  ATOM.h
//  HW_1_v1.0
//
//  Created by 陳啟瑋 on 2023/3/8.
//

#ifndef ATOM_h
#define ATOM_h

#include <iostream>
#include <cmath>

using namespace std;

class ATOM{
    public:
        ATOM(){
            x[0]=x[1]=x[2]=chg=Ro=Do=0;
        };
    public:
        double x[3]; //position vector of atom in angstroms
        double chg; //partial charges in electrons
        double Ro; //LJ Ro parameter in Angstroms
        double Do; //LJ Do parameter in kcal/mol
};

#endif /* ATOM_h */
