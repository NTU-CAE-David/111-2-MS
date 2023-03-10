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
    
    ATOM a[2][3]; // a["分子編號"]["分子內原子編號"]
    
    // H-O-H (No.1)
    //data for Oxygen
    a[0][0].x[0]=0;
    a[0][0].x[1]=0;
    a[0][0].x[2]=3.36;
    a[0][0].chg=-0.82;
    a[0][0].Ro=3.5532;
    a[0][0].Do=0.1848;
    
    //data for Hydrogen-1
    a[0][1].x[0]=0.8;
    a[0][1].x[1]=0;
    a[0][1].x[2]=3.95;
    a[0][1].chg=0.41;
    a[0][1].Ro=0.9;
    a[0][1].Do=0.01;
    
    //data for Hydrogen-2
    a[0][2].x[0]=-0.8;
    a[0][2].x[1]=0;
    a[0][2].x[2]=-0.95;
    a[0][2].chg=0.41;
    a[0][2].Ro=0.9;
    a[0][2].Do=0.01;
    
    // H-O-H (No.2)
    //data for Oxygen
    a[1][0].x[0]=0;
    a[1][0].x[1]=0;
    a[1][0].x[2]=0;
    a[1][0].chg=-0.82;
    a[1][0].Ro=3.5532;
    a[1][0].Do=0.1848;
    
    //data for Hydrogen-1
    a[1][1].x[0]=0;
    a[1][1].x[1]=0;
    a[1][1].x[2]=1;
    a[1][1].chg=0.41;
    a[1][1].Ro=0.9;
    a[1][1].Do=0.01;
    
    //data for Hydrogen-2
    a[1][2].x[0]=0;
    a[1][2].x[1]=0.9;
    a[1][2].x[2]=-0.4;
    a[1][2].chg=0.41;
    a[1][2].Ro=0.9;
    a[1][2].Do=0.01;
    
    double total_eng = 0;
    
    // Nonbond
    NONBOND nbs[2][3]; // 創建一個大小為 3 的 NONBOND 物件陣列
    nbs[0][0].cal_eng(a[0][0], a[0][1]); // 計算第一對原子之間的非鍵能
    nbs[0][1].cal_eng(a[0][0], a[0][2]); // 計算第二對原子之間的非鍵能
    nbs[0][2].cal_eng(a[0][1], a[0][2]); // 計算第三對原子之間的非鍵能
    
    // TODO 找出有多少個 non-bond 並將所有cal_eng列出來並計算。
    
    for (int i = 0; i < 3; i++) {
        cout << "NONMOND: " << i << endl;
        cout << "vdw: " << nbs[0][i].vdw
             << " kcal/mol, cou: " << nbs[0][i].cou << " kcal/mol\n"
             << "eng: " << nbs[0][i].eng << " kcal/mol\n" << endl;
        total_eng += nbs[0][i].eng;
    }
    
    // bond
    BOND bonds[2][2]; // bonds["分子編號"]["分子內鍵結編號"]
    bonds[0][0].Kb = 500; // kcal/mole
    bonds[0][0].len0 = 1;
    bonds[0][0].cal_eng(a[0][0], a[0][1]);

    bonds[0][1].Kb = 500; // kcal/mole
    bonds[0][1].len0 = 1;
    bonds[0][1].cal_eng(a[0][0], a[0][2]);
    
    bonds[1][0].Kb = 500; // kcal/mole
    bonds[1][0].len0 = 1;
    bonds[1][0].cal_eng(a[1][0], a[1][1]);

    bonds[1][1].Kb = 500; // kcal/mole
    bonds[1][1].len0 = 1;
    bonds[1][1].cal_eng(a[1][0], a[1][2]);

    for (int j = 0; j < 2; j++){
        for (int i = 0; i < 2; i++) {
            cout << "分子編號: " << j
                << " BOND: " << i << endl;
            cout << "bond length is " << bonds[j][i].len
                <<" A, energy is "<< bonds[j][i].eng << " kcal/mol\n"<<endl;
            total_eng += bonds[j][i].eng;
        }
    }
    
    // angle
    ANGLE angle[2];
    
    double theta[2];  // degrees
    theta[0] = angle[0].angle(a[0][0], a[0][1], a[0][2]);
    theta[1] = angle[1].angle(a[1][0], a[1][1], a[1][2]);
    
    angle[0].k0 = 120;  // kcal/mol
    angle[0].theta0 = 109.47;  // degrees
    angle[0].energy(theta[0]);
    
    angle[1].k0 = 120;  // kcal/mol
    angle[1].theta0 = 109.47;  // degrees
    angle[1].energy(theta[1]);
    
    // TODO 求出兩個分子各自的 angle
    
    double pe = angle[0].potential_energy(theta[0]);
    double f = angle[0].force(theta[0]);

    cout <<  "Energy of HOH angle[0] at " << theta[0]
        << " degrees: " << angle[0].eng << " kcal/mol." << endl;
    cout << "Potential energy: " << pe << " kcal/mol\n";
    cout << "Force: " << f << " kcal/(mol.deg)\n";
    cout << "Equilibrium angle[0]: " << angle[0].theta0 << " degrees" << endl;
    
    
    total_eng += angle[0].eng;
    cout << "Total Energy = " << total_eng << " kcal/mol." <<endl;
    
    return 0;

}
