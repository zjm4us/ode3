///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
///
/// DO NOT CHANGE the required print block.
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <vector>

using namespace std;

struct Params {
    double g;
    double m;
    double d;
    double b;
    double c;
};

//
// Derivative functions for RK4StepN
//
// y = { x, y_dummy, z, vx, vy_dummy, vz }
//

double fx(double t, const vector<double> &y, void *p) {
    return y[3];       // dx/dt = vx
}

double fy(double t, const vector<double> &y, void *p) {
    return 0.0;        // no y motion
}

double fz(double t, const vector<double> &y, void *p) {
    return y[5];       // dz/dt = vz
}

double fvx(double t, const vector<double> &y, void *p) {
    Params *pars = (Params*) p;

    double vx = y[3];
    double vz = y[5];
    double v = sqrt(vx*vx + vz*vz);
    if (v < 1e-12) v = 1e-12;

    double D = pars->b * v + pars->c * v * v;
    double Fdx = -D * (vx / v);
    return Fdx / pars->m;
}

double fvy(double t, const vector<double> &y, void *p) {
    return 0.0;        // no y motion
}

double fvz(double t, const vector<double> &y, void *p) {
    Params *pars = (Params*) p;

    double vx = y[3];
    double vz = y[5];
    double v = sqrt(vx*vx + vz*vz);
    if (v < 1e-12) v = 1e-12;

    double D = pars->b * v + pars->c * v * v;
    double Fdz = -D * (vz / v);
    return Fdz / pars->m - pars->g;
}

int main(int argc, char **argv){

    Params pars;
    pars.g = 9.81;
    pars.m = 0.145;
    pars.d = 0.075;                         // correct diameter
    pars.b = 1.6e-4 * pars.d;
    pars.c = 0.25 * pars.d * pars.d;

    void *p_par = (void*) &pars;

    double xend = 18.5;
    double z0 = 1.4;
    double theta0 = 1;
    bool showPlot = false;

    int c;
    while ((c = getopt(argc, argv, "x:z:t:p")) != -1)
        switch (c) {
        case 'x': xend = atof(optarg); break;
        case 'z': z0 = atof(optarg); break;
        case 't': theta0 = atof(optarg); break;
        case 'p': showPlot = true; break;
        }

    TApplication theApp("App", &argc, argv);

    double theta = theta0 * M_PI / 180.0;

    // List of derivative functions
    vector<pfunc_t> fnlist = {fx, fy, fz, fvx, fvy, fvz};

    //
    // function to simulate trajectory for a given v0
    //
    auto simulate = [&](double v0) {

        vector<double> y(6);

        y[0] = 0.0;
        y[1] = 0.0;
        y[2] = z0;

        y[3] = v0 * cos(theta);
        y[4] = 0.0;
        y[5] = v0 * sin(theta);

        double t = 0.0;
        double dt = 0.0005;

        while (y[0] < xend && y[2] > 0.0) {
            y = RK4StepN(fnlist, y, t, dt, p_par);
            t += dt;
        }

        return y[2];   // return z at plate
    };

    //
    // Bisection to solve for vPitch
    //
    double targetZ = 0.9;
    double low = 25.0;
    double high = 60.0;

    for (int i = 0; i < 50; i++) {
        double mid = 0.5 * (low + high);
        double zfinal = simulate(mid);

        if (zfinal > targetZ)
            high = mid;
        else
            low = mid;
    }

    double vPitch = 0.5 * (low + high);

    //
    // REQUIRED OUTPUT (DO NOT CHANGE)
    //
    printf("********************************\n");
    printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n", xend, z0, theta0);
    printf("v_pitch = %lf m/s\n", vPitch);
    printf("********************************\n");

    if (showPlot) {
        cout << "Press ^c to exit" << endl;
        theApp.SetIdleTimer(30,".q");
        theApp.Run();
    }

    return 0;
}

