#include "helper.h"

void get_atmospheric_properties_si(double altitude, Atmosphere&atm)
{
    // variables
    int i;
    bool check=false;
    double Z,T,temp,rho,a,mu;
    double p=101325;
    double geopotential_altitude;
    double Zi[9] = {0,11000,20000,32000,47000,52000,61000,79000,90000};
    double Ti[8] = {288.15,216.65,216.65,228.65,270.65,270.65,252.65,180.65};
    double dTi[8] = {-0.0065,0,0.001,0.0028,0,-0.002,-0.004,0};
    double pi[8];
    // find geopotential altitude
    Z = atm.R_earth * altitude / (atm.R_earth + altitude);
    
    // find temperature
    for(i = 0; i<9;i++){
        if (Z < Zi[i+1]){
            T = Ti[i] + dTi[i]*(Z-Zi[i]);
            break;
        } else{
            continue;
        }
    }

    // find pressure
    int j = 0;
    while(!check){
        if(Z<Zi[j+1]){
            temp=Z;
            check = true;
        }else{
            temp = Zi[j+1];
        }
        if (dTi[j]==0){
            p *= exp(-atm.g0*(temp-Zi[j])/(atm.R*Ti[j]));
        }else{
            p *= pow(((Ti[j]+dTi[j]*(temp-Zi[j]))/Ti[j]),(-atm.g0/(atm.R*dTi[j])));
        }
        j++;
    }

    //find density
    rho = p / atm.R / T;
    
    // find speed of sound
    a = sqrt(atm.gamma * atm.R * T);
 
    // find viscocity
    double T0 = 273.15; // Kelvin
    double mu0 = 0.00001716; //kg.m-s
    double C = 110.4; // Southerland constant
    mu = mu0*(T0+C)/(T+C)*pow(T/T0,1.5); 

    // save answer
    atm.geopotential_altitude = Z;
    atm.temperature = T;
    atm.pressure = p;
    atm.density = rho;
    atm.speed_of_sound = a;
    atm.viscosity = mu;
    

}

void get_atmospheric_properties_english(double altitude, Atmosphere& atm){
    // convert altitude from ft to m
    altitude *= 0.3048;
    
    // get SI properties
    get_atmospheric_properties_si(altitude, atm);

    // convert SI properties to english units
    atm.geopotential_altitude /= 0.3048; // from https://www.grc.nasa.gov/www/winddocs/cff/factors.html
    atm.temperature *= 1.8; // from https://www.grc.nasa.gov/www/winddocs/cff/factors.html
    atm.pressure *= 0.020885434304801722; // from email
    atm.density *= 0.00194032032363104; // from email
    atm.speed_of_sound /= 0.3048; // from https://www.grc.nasa.gov/www/winddocs/cff/factors.html
    atm.viscosity /= 47.88025898; // slugs/(ft-s) from Dynamic Viscosity https://www.cfd-online.com/Wiki/Sutherland%27s_law
}

void print_atmosphere()
{
    printf("Printing Atmospheric Tables\n");

    Atmosphere my_atm;
    FILE* si_file = fopen("stdatmos_si.txt","w");

    fprintf(si_file,"  Geometric_Alt[m]     Geopotential_Alt[m]  Temperature[K]       Pressure[N/m^2]      Density[kg/m^3]      Speed_of_Sound[m/s]\n");
    for (double H = 0.0; H < 72000.0; H += 2000.0) {
        get_atmospheric_properties_si(H,my_atm);
        fprintf(si_file,"%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n",H,my_atm.geopotential_altitude,my_atm.temperature, my_atm.pressure, my_atm.density, my_atm.speed_of_sound);
    }
    fclose(si_file);

    FILE* en_file = fopen("stdatmos_english.txt", "w");
    fprintf(en_file,"  Geometric_Alt[ft]    Geopotential_Alt[ft] Temperature[R]       Pressure[lbf/ft^2]   Density[slugs/ft^3]  Speed_of_Sound[ft/s]\n");
    for (double H = 0.0; H < 180000.0; H += 5000.0) {
        get_atmospheric_properties_english(H,my_atm);
        fprintf(en_file,"%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n",H,my_atm.geopotential_altitude,my_atm.temperature, my_atm.pressure, my_atm.density, my_atm.speed_of_sound);
    }
    fclose(en_file);
}

double f_scalar(double t,double y){
    double ans;
    ans = 1 + tan(y);
    return ans;
}
double rk4_scalar(double t0, double y0, double dt){
    double y,k1,k2,k3,k4;

    k1 = f_scalar(t0,y0);
    k2 = f_scalar(t0 + dt/2, y0 + k1*dt/2);
    k3 = f_scalar(t0 + dt/2,y0 + k2*dt/2);
    k4 = f_scalar(t0 + dt, y0 + k3*dt);
    y = y0 + (dt/6)*(k1+2*k2+2*k3+k4);
    return y;
}

void f_array(double t, double* x,double* ans){

    ans[0] = 1+pow(x[1],2)*sin(x[0]);
    ans[1] = 1 + x[0]*cos(x[1]);
}

void rk4_array(double t0, double y0[], double dt, int size, double* ans){
    double* k1 = new double[size];
    double* k2 = new double[size];
    double* k3 = new double[size];
    double* k4 = new double[size];
    double* y_in = new double[size];
    // calculate k1
    f_array(t0,y0,k1);
    // calculate k2
    for (int i = 0; i<size;i++){
        y_in[i] = y0[i] + k1[i] * dt/2;
    }
    f_array(t0 + dt/2, y_in ,k2);
    //calculate k3
    for (int i = 0; i<size;i++){
        y_in[i] = y0[i] + k2[i] * dt/2;
    }
    f_array(t0 + dt/2, y_in,k3);
    // calculate k4
    for (int i = 0; i<size;i++){
        y_in[i] = y0[i] + k3[i] * dt;
    }
    f_array(t0 + dt, y_in,k4);
    //calculate y next 
    for (int i = 0; i < size; i++) {
        ans[i] = y0[i] + (dt/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }

    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
}


void testRK4(){
    double dt = 0.025;
    double t = 0;
    double y0[2] = {0,0};
    double y[2] = {0,0};

    FILE* rk_out = fopen("rk4_out.txt","w");
    fprintf(rk_out,"Time(s), X, Z\n");
    while (t <= 0.25) {
        fprintf(rk_out,"%20.12e,%20.12e,%20.12e\n",t,y[0],y[1]);
        rk4_array(t, y0, dt,2,y);
        memcpy(y0,y,sizeof(y));
        t += dt;
    }
    fclose(rk_out);

}

double gravity_si(double altitude){
    return 9.806645 * pow(6356766.0/(6356766.0+altitude),2);
}
double gravity_english(double altitude){
    double ans;
    altitude *= 0.3048;
    ans = gravity_si(altitude);
    ans /= 0.3048;
    return ans; 
}

void print_array(double* values, int size)
{
    for(int i=0; i<size; i++)
    {
        printf("%20.12e",values[i]);
    }
    printf("\n\n");
}

void array_copy(double* A, double* B, int size){
    for(int i=0; i<size; i++)
    {
        B[i]=A[i];
    }
}


void quat_mult(double* A, double* B, double* ans){
// this function assumes A and B are quaternions and of size 4
    // finding scalar component
    ans[0] = A[0] * B[0] - A[1] * B[1] - A[2] * B[2] - A[3] * B[3];
    // finding x component
    ans[1] = A[0] * B[1] + A[1] * B[0] + A[2] * B[3] - A[3] * B[2];
    // finding y component
    ans[2] = A[0] * B[2] - A[1] * B[3] + A[2] * B[0] + A[3] * B[1];
    // finding z component
    ans[3] = A[0] * B[3] + A[1] * B[2] - A[2] * B[1] + A[3] * B[0];
}

void quat_norm(double* quat){
    // this function takes a size 4 array quaternion as an input and normalizes it with its length.  It updates the input variable with the result.
    double mag = (sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]));
    quat[0] /= mag;
    quat[1] /= mag;
    quat[2] /= mag;
    quat[3] /= mag;
}
void quat_norm_quick(double* quat){
    double eps = 1 - (quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);
    quat[0] *= (1+0.5*eps);
    quat[1] *= (1+0.5*eps);
    quat[2] *= (1+0.5*eps);
    quat[3] *= (1+0.5*eps);
}
void euler_to_quat(double* eul, double* quat){
    // this function convers the input eul (a size 3 radian array) to quaternions quat (a size 4 array)
    double phi = eul[0]*0.5;
    double theta = eul[1]*0.5;
    double psi = eul[2]*0.5;
    double s_phi2 = sin(phi);
    double c_phi2 = cos(phi);
    double s_theta2 = sin(theta);
    double c_theta2 = cos(theta);
    double s_psi2 = sin(psi);
    double c_psi2 = cos(psi);

    quat[0] = c_phi2*c_theta2*c_psi2 + s_phi2*s_theta2*s_psi2; //e_0
    quat[1] = s_phi2*c_theta2*c_psi2 - c_phi2*s_theta2*s_psi2; //e_x
    quat[2] = c_phi2*s_theta2*c_psi2 + s_phi2*c_theta2*s_psi2; //e_y
    quat[3] = c_phi2*c_theta2*s_psi2 - s_phi2*s_theta2*c_psi2; //e_z

    // quat[0] = cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2); //e_0
    // quat[1] = sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2); //e_x
    // quat[2] = cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2); //e_y
    // quat[3] = cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2); //e_z
}

void quat_to_euler(double* quat, double* eul){
    // this function converts the input quat (a normalized size 4 quaternian array) into euler angles eul(a size 3 radian array)
    double e0 = quat[0];
    double ex = quat[1];
    double ey = quat[2];
    double ez = quat[3];
    double test = e0 * ey - ex * ez;

    if (test == 0.5){
        eul[0] = 2. * asin(ex/cos(pi*0.25)); //phi
        eul[1] = pi/2.; //theta
        eul[2] = 0.; // psi
    } else if (test == -0.5){
        eul[0] = 2. * asin(ex/cos(pi*0.25)); //phi
        eul[1] = -pi/2.; //theta
        eul[2] = 0.; // psi
    } else{
        eul[0] = atan2((2. * (e0 * ex + ey * ez)),(e0 * e0 - ex * ex - ey * ey + ez * ez  )); //phi
        eul[1] = asin((2. * (e0 * ey - ex * ez))); //theta
        eul[2] = atan2((2. * (e0 * ez + ex * ey)),(e0 * e0 + ex * ex - ey * ey - ez * ez)); // psi
    }



}
void body_to_fixed(double* vec, double* quat, double* ans){
    // this function converts the inputs vec (a size 3 array in body-fixed coordinates) and quat (a size 4 quaternion array) into ans(a size 3 array in earth fixed coordinates)
    double e0 = quat[0];
    double ex = quat[1];
    double ey = quat[2];
    double ez = quat[3];

    double x = vec[0];
    double y = vec[1];
    double z = vec[2];

    ans[0] = (ex * ex + e0 * e0 - ey * ey - ez * ez) * x + 2 * (ex * ey - ez * e0) * y + 2 * (ex * ez + ey * e0) * z; 
    ans[1] = 2 * (ex * ey + ez * e0) * x + (ey * ey + e0 * e0 - ex * ex - ez * ez) * y + 2 * (ey * ez - ex * e0) * z;
    ans[2] = 2 * (ex * ez - ey * e0) * x + 2 * (ey * ez + ex * e0) * y + (ez * ez + e0 * e0 - ex * ex - ey * ey) * z;
}

void fixed_to_body(double* vec, double* quat, double* ans){
    // this function converts the inputs vec (a size 3 array in earth fixed coordinates) and quat (a size 4 quaternion array) into ans(a size 3 array in body-fixed coordinates) 
    double e0 = quat[0];
    double ex = quat[1];
    double ey = quat[2];
    double ez = quat[3];

    double x = vec[0];
    double y = vec[1];
    double z = vec[2];

    ans[0] = (ex * ex + e0 * e0 - ey * ey - ez * ez) * x + 2 * (ex * ey + ez * e0) * y + 2 * (ex * ez - ey * e0) * z; 
    ans[1] = 2 * (ex * ey - ez * e0) * x + (ey * ey + e0 * e0 - ex * ex - ez * ez) * y + 2 * (ey * ez + ex * e0) * z;
    ans[2] = 2 * (ex * ez + ey * e0) * x + 2 * (ey * ez - ex * e0) * y + (ez * ez + e0 * e0 - ex * ex - ey * ey) * z;
}

// void print_square_array(double* in,int n){
    
//      for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
            
//             cout << *(in + i * n + j) << " ";
//         }
//         cout << endl;
//     }
// }

// void matMul(double m1[][3],double m2[][3],double result[][3],int rows1,int cols1, int rows2, int cols2){


//     if (cols1 != rows2) {
//             cerr << "Error: The number of columns in the first matrix must be equal to the number of rows in the second matrix." << endl;
//             return;
//         }

//     for (int i = 0; i < rows1; ++i) {
//         for (int j = 0; j < cols2; ++j) {
//             result[i][j] = 0;
//             for (int k = 0; k < cols1; ++k) {
//                 result[i][j] += m1[i][k] * m2[k][j];
//             }
//         }
//     }
// }
