//
//  ANGLE.h
//  HW_1_v1.0
//
//  Created by 陳啟瑋 on 2023/3/8.
//

#ifndef ANGLE_h
#define ANGLE_h

#include <iostream>
#include <cmath>
#include "ATOM.h"

using namespace std;

class ANGLE{
    public:
        ANGLE(){
            k0=theta0=eng=0;
        }
        double k0;      // Force constant in kcal/mol
        double theta0;  // Equilibrium angle in degrees
        double eng;     // angle energy in kcal/mol
        double energy(double theta);
        double potential_energy(double theta) const;
        double force(double theta) const;
    
};

double ANGLE::energy(double theta) {
    double rad_theta = theta * M_PI / 180.0;
    double rad_theta0 = theta0 * M_PI / 180.0;
    double cos_sin_delta = cos(rad_theta - rad_theta0) / sin(theta0);
    eng = 0.5 * k0 * cos_sin_delta * cos_sin_delta;
    return eng;
}

double ANGLE::potential_energy(double theta) const {
    double delta_theta = theta - theta0;
    return 0.5 * k0 * delta_theta * delta_theta;
}

double ANGLE::force(double theta) const {
    double delta_theta = theta - theta0;
    return -k0 * delta_theta;
}

#endif /* ANGLE_h */
