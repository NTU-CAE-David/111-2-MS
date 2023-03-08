//
//  main.cpp
//  HW_1_v1.0
//
//  Created by 陳啟瑋 on 2023/3/8.
//

#include <iostream>
#include <cmath>
#include "ATOM.h"
#include "BOND.h"
#include "NONBOND.h"
#include "ANGLE.h"

using namespace std;

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
