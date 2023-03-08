//
//  main.cpp
//  HW_1
//
//  Created by 陳啟瑋 on 2023/3/3.
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

class ANGLE {
    public:
        ANGLE(){
            k0=theta0=eng=0;
        }
        double k0;      // Force constant in kcal/mol
        double theta0;  // Equilibrium angle in degrees
        double eng;     // angle energy in kcal/mol
    
        double energy(double theta) {
            double rad_theta = theta * M_PI / 180.0;
            double rad_theta0 = theta0 * M_PI / 180.0;
            double cos_sin_delta = cos(rad_theta - rad_theta0) / sin(theta0);
            eng = 0.5 * k0 * cos_sin_delta * cos_sin_delta;
            return eng;
        }
        
        double potential_energy(double theta) const {
            double delta_theta = theta - theta0;
            return 0.5 * k0 * delta_theta * delta_theta;
        }

        double force(double theta) const {
            double delta_theta = theta - theta0;
            return -k0 * delta_theta;
        }
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
    
    return eng = vdw + cou;
}

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



int main() {
    
    ATOM a[3];
    
    // H-O-H
    
    //data for Oxygen
    a[0].x[0]=0;
    a[0].x[1]=0;
    a[0].x[2]=0;
    a[0].chg=-0.82;
    a[0].Ro=3.5532;
    a[0].Do=0.1848;
    
    //data for Hydrogen-1
    a[1].x[0]=0;
    a[1].x[1]=0;
    a[1].x[2]=1;
    a[1].chg=0.41;
    a[1].Ro=0.9;
    a[1].Do=0.01;
    
    //data for Hydrogen-2
    a[2].x[0]=0;
    a[2].x[1]=0.9;
    a[2].x[2]=-0.35;
    a[2].chg=0.41;
    a[2].Ro=0.9;
    a[2].Do=0.01;
    
    double total_eng = 0;
    
    // Nonbond
    NONBOND nbs[3]; // 創建一個大小為 3 的 NONBOND 物件陣列
    nbs[0].cal_eng(a[0], a[1]); // 計算第一對原子之間的非鍵能
    nbs[1].cal_eng(a[0], a[2]); // 計算第二對原子之間的非鍵能
    nbs[2].cal_eng(a[1], a[2]); // 計算第三對原子之間的非鍵能
    
    for (int i = 0; i < 3; i++) {
        cout << "NONMOND: " << i << endl;
        cout << "vdw: " << nbs[i].vdw
             << " kcal/mol, cou: " << nbs[i].cou << " kcal/mol\n"
             << "eng: " << nbs[i].eng << " kcal/mol\n" << endl;
        total_eng += nbs[i].eng;
    }
    
    // bond
    BOND bonds[2];
    bonds[0].Kb = 500; // kcal/mole
    bonds[0].len0 = 1;
    bonds[0].cal_eng(a[0], a[1]);

    bonds[1].Kb = 500; // kcal/mole
    bonds[1].len0 = 1;
    bonds[1].cal_eng(a[0], a[2]);

    for (int i = 0; i < 2; i++) {
        cout << "BOND: " << i << endl;
        cout<<"bond length is "<<bonds[i].len
            <<" A, energy is "<<bonds[i].eng<<" kcal/mol\n"<<endl;
        total_eng += bonds[i].eng;
    }
    
    // angle
    ANGLE angle;
    double theta = 111.251;  // degrees
    angle.k0 = 120;  // kcal/mol
    angle.theta0 = 109.47;  // degrees
    angle.energy(theta);
    
    double pe = angle.potential_energy(theta);
    double f = angle.force(theta);

    cout <<  "Energy of HOH angle at " << theta
        << " degrees: " << angle.eng << " kcal/mol." << endl;
    cout << "Potential energy: " << pe << " kcal/mol\n";
    cout << "Force: " << f << " kcal/(mol.deg)\n";
    cout << "Equilibrium angle: " << angle.theta0 << " degrees" << endl;
    
    
    total_eng += angle.eng;
    cout << "Total Energy = " << total_eng << " kcal/mol." <<endl;
    

    
    return 0;

}
