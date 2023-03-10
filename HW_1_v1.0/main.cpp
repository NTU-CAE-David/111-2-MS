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
    
    double internal_eng[2] = {}; // Internal Energy
    double total_eng = 0.0; // Total Energy
    
    // Nonbond
    NONBOND nbs[3][3]; // 創建一個大小為 3*3 的 NONBOND 物件陣列
    // nbs[No.1->原子編號][No.2->原子編號]
    nbs[0][0].cal_eng(a[0][0], a[0][1]); // 計算第一對原子之間的非鍵能
    nbs[0][1].cal_eng(a[0][0], a[0][2]); // 計算第二對原子之間的非鍵能
    nbs[0][2].cal_eng(a[0][1], a[0][2]); // 計算第三對原子之間的非鍵能
    
    // 找出有多少個 non-bond 並將所有cal_eng列出來並計算。
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            nbs[i][j].cal_eng(a[0][i], a[1][j]); // 計算第一對原子之間的非鍵能
            
            cout << "NONMOND: "
                << "0" << i << "->"
                << "1" << j << endl;
            cout << "vdw: " << nbs[i][j].vdw
                 << " kcal/mol, cou: " << nbs[i][j].cou << " kcal/mol\n"
                 << "eng: " << nbs[i][j].eng << " kcal/mol\n" << endl;
            
            // calculate Energy
//            internal_eng[i] += nbs[i][j].eng;
            total_eng += nbs[i][j].eng;
        }
    }
    
    // bond
    BOND bonds[2][2];
    // bonds["分子編號"]["分子內鍵結編號"]
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
            
            // calculate Energy
            internal_eng[j] += bonds[j][i].eng;
            total_eng += bonds[j][i].eng;
        }
    }
    
    // angle
    // angle, theta["分子編號"]
    ANGLE angle[2];
    double theta[2];  // degrees
    
    // Initi
    for (int i = 0; i < 2; i++) {
        
        theta[i] = angle[i].angle(a[i][0], a[i][1], a[i][2]);
        
        angle[i].k0 = 120;  // kcal/mol
        angle[i].theta0 = 109.47;  // degrees
        angle[i].energy(theta[i]);
    }
    
    // 求出兩個分子各自的 angle
    for (int i = 0; i < 2; i++) {
        angle[i].potential_energy(theta[i]);
        angle[i].force(theta[i]);
        
        cout <<  "Energy of HOH angle[" << i << "] at " << theta[i]
        << " degrees: " << angle[i].eng << " kcal/mol." << endl;
        cout << "Potential energy: " << angle[i].pe << " kcal/mol\n";
        cout << "Force: " << angle[i].f << " kcal/(mol.deg)\n";
        cout << "Equilibrium angle[" << i << "]: " << angle[i].theta0 << " degrees" << endl;
        cout << endl;
        
        // calculate Energy
        internal_eng[i] += angle[i].eng;
        total_eng += angle[i].eng;
    }
    
    
    cout << "=================== HW 問題回答 ===================\n" << endl;
    
    cout << "Q.2 (a): " << endl;
    cout << "Internal Energy [0] = " << internal_eng[0] << " kcal/mol." <<endl;
    cout << "Internal Energy [1] = " << internal_eng[1] << " kcal/mol." <<endl;
    
    cout << "Q.2 (b): " << endl;
    cout << "Total Energy = " << total_eng << " kcal/mol." <<endl;
    
    cout << "Q.2 (c): " << endl;
    double binding_eng =  total_eng - internal_eng[0] - internal_eng[1];
    cout << "Binding Energy = " << binding_eng << " kcal/mol." <<endl;
    
    return 0;

}
