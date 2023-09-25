#include <iostream>
#include <fstream>
#include <cmath>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;


//declare constants
const double pi = 3.141592653589793238462643383279502884197;

// declare structs
struct Atmosphere
{
    const double R_earth=6356766, g0=9.806645, R=287.0528, gamma=1.4;
    double geopotential_altitude;
    double temperature;
    double pressure;
    double density;
    double speed_of_sound;

};

// declare functions
// helpfull additional functions
void print_array(double* values, int size);
void array_copy(double* A, double* B,int size);
// void print_square_array(double* in,int n);
// void matMul(double m1[][3],double m2[][3],double result[][3],int rows1,int cols1, int rows2, int cols2);
// functions from atmosphere
void get_atmospheric_properties_si(double altitude, Atmosphere&atm);
void get_atmospheric_properties_english(double altitude, Atmosphere&atm);
void print_atmosphere();
double gravity_si(double altitude);
double gravity_english(double altitude);

// functions from RK4
void testRK4();
double f_scalar(double t,double y);
double rk4_scalar(double t0, double y0, double dt);
void f_array(double t, double* x,double* ans);
void rk4_array(double t0, double y0[], double dt, int size, double* ans);

// quaternion functions
void quat_mult(double* A, double* B, double* ans);
void quat_norm(double* quat);
void quat_norm_quick(double* quat);
void euler_to_quat(double* eul, double* quat);
void quat_to_euler(double* quat, double* eul);
void body_to_fixed(double* vec, double* quat, double* ans);
void fixed_to_body(double* vec, double* quat, double* ans);
