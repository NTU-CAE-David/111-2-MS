//
//  main.cpp
//  W2
//
//  Created by 陳啟瑋 on 2023/3/2.
//

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

class NONBOND{
    public:
        NONBOND(){
            r=r2=vdw=cou=0;
        };
        double r; //separation distance
        double r2; //r*r
        double vdw; //van der waals energy in kcal/mol
        double cou; //coulomb energy in kcal/mol
        double cal_r(ATOM &,ATOM &); //calculate r and r2
        double cal_eng(ATOM &,ATOM &); //calculate nonbond energy
};

class BOND{
    public:
        BOND(){
            len=len0=Kb=eng=0;
        };
//        ATOM *atm[2];
        double len; //in angstroms
        double len0; //equilib bond length in angstroms
        double Kb; //force constant in kcal/mol A2
        double eng; //bond energy in kcal/mol
        double cal_len(ATOM &a,ATOM &b); //calculate bond length from two vectors
        double cal_eng(ATOM &a,ATOM &b); //calculate bond energy
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
    
    return vdw + cou;
}

double BOND::cal_len(ATOM &a,ATOM &b){
    
    double r2;
    
//    r2 = (atm[0]->x[0] - atm[1]->x[0]) * (atm[0]->x[0] - atm[1]->x[0])
//        +(atm[0]->x[1] - atm[1]->x[1]) * (atm[0]->x[1] - atm[1]->x[1])
//        +(atm[0]->x[2] - atm[1]->x[2]) * (atm[0]->x[2] - atm[1]->x[2]);
    
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


int main() {
    
    ATOM a[2];
    
    //data for Oxygen
    a[0].x[0]=0;
    a[0].x[1]=0;
    a[0].x[2]=0;
    a[0].chg=-0.82;
    a[0].Ro=3.5532;
    a[0].Do=0.1848;
    
    //data for Hydrogen
    a[1].x[0]=0;
    a[1].x[1]=0.7;
    a[1].x[2]=-0.7;
    a[1].chg=0.41;
    a[1].Ro=0.9;
    a[1].Do=0.01;
    
    // Nonbond
    NONBOND nb;
    nb.cal_eng(a[0], a[1]);
    
    cout << "vdw: " << nb.vdw
         << " kcal/mol, cou: "
         << nb.cou << " kcal/mol" << endl;
    
    // bond
    BOND bond;
//    bond.atm[0] = &a[0];
//    bond.atm[1] = &a[1];
    bond.Kb = 500; // kcal/mole
    bond.len0 = 1;
    bond.cal_eng(a[0], a[1]);
    cout<<"bond length is "<<bond.len<<" A, energy is "<<bond.eng<<" kcal/mol"<<endl;

    
    return 0;

}
