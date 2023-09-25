#ifndef SPHERE_H
#define SPHERE_H
#include "helper.h"


class sphere{
public:
    sphere(string filename);
    ~sphere(){};
    void run_simulation();

private:
    double m_state[13];
    int m_i_max=100000;
    double m_time_step,m_rho_sphere,m_mass,m_r,m_W,m_Ixx,m_Iyy,m_Izz;
    int m_counter=0;

    void get_state_array_delta(double* y0,double* ans);
    void aerodynamics(double* y0, double* ans);
    void rk4();
    

};

#endif