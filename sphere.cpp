#include "sphere.h"


sphere::sphere(string filename){
    std::ifstream f(filename);
    json data = json::parse(f);
    double airspeed,altitude,elevation_angle,temp[3],mass;

    // simulation data
    m_time_step = data["simulation"]["time_step[s]"];

    // initial state
    airspeed = data["initial"]["airspeed[ft/s]"];
    elevation_angle = data["initial"]["elevation_angle[deg]"];
    altitude = data["initial"]["altitude[ft]"];
    temp[0] = 0;
    temp[1] = pi/180*elevation_angle;
    temp[2] = 0;

    // update this still
    m_state[0] = airspeed;
    m_state[1] = 0;
    m_state[2] = 0;
    m_state[3] = data["initial"]["initial_p"];
    m_state[4] = data["initial"]["initial_q"];
    m_state[5] = data["initial"]["initial_r"];
    m_state[6] = 0;
    m_state[7] = 0;
    m_state[8] = -altitude;
    euler_to_quat(temp,&m_state[9]);

    //mass props
    
    m_rho_sphere = data["geometry"]["density[slug/ft^3]"];
    m_r = data["geometry"]["radius[in]"];
    m_r /= 12;
    m_mass = 4/3*pi*m_r*m_r*m_r*m_rho_sphere;
    m_Ixx = 0.4*m_mass*m_r*m_r;
    m_Iyy = m_Ixx;
    m_Izz = m_Ixx;
    m_Sref = 2.0 *pi * m_r*m_r;
    m_lref = 2. * m_r;


}



void sphere::aerodynamics(double* y0,double* ans){
    // get inputs
    double u = y0[0];
    double v = y0[1];
    double w = y0[2];
    double p = y0[3];
    double q = y0[4];
    double r = y0[5];
    double xf = y0[6];
    double yf = y0[7];
    double zf = y0[8];
    double e0 = y0[9];
    double ex = y0[10];
    double ey = y0[11];
    double ez = y0[12];

    // init variables
    double V,CL,CS,CD,Cl,Cm,Cn,alpha,beta,Fx,Fy,Fz,Mx,My,Mz,Fx_spin,Fy_spin,Fz_spin;
    double salpha,sbeta,calpha,cbeta,coef1,coef2;
    double RE;
    double mu;
    double g;


    // get atmospheric properties
    Atmosphere atm;
    get_atmospheric_properties_english(-zf,atm);

    // get total velocity 
    V = sqrt(u * u + v * v + w * w);

    // get alpha and beta
    alpha = atan2(w,u);
    beta = atan2(v,V);

    // get forces and moments from rotation
    coef2 = 0.15;
    coef1 = 16/3* coef2 * pi*pi * m_r*m_r*m_r *atm.density;
    Fx_spin = coef1 * (q*w - r*v);
    Fy_spin = coef1 * (r*u - p*w);
    Fz_spin = coef1 * (p*v - q*u);

     // get gravity for current state
    g = gravity_english(-zf);

    // Calculate the viscosity of the air
    //mu = ((m_rho_sphere-atm.density)/VT)*(2.0/9.0)*g*m_r*m_r;

    // Calculate the reynolds number
    RE = atm.density*V*2.0*m_r/atm.viscosity;

    // Calculate CD for a sphere
    if (RE>0.0 && RE <=450000) {
        CD = 24.0/RE + 6/(1+sqrt(RE)) + 0.4;
    } else if (RE>450000 && RE <=560000) {
        CD = 1.0*(pow(10.0,29.0))*pow(RE,-5.211);
    } else if (RE>560000 && RE <=14000000) {
        CD = -2.0*(pow(10.0,-23.0))*RE*RE*RE - 1.0*(pow(10.0,-16.0))*RE*RE + 9.0*(pow(10.0,-9.0))*RE + 0.069;
    } else {
        CD = 0.12;
    }


    // precompute sin and cos for alpha and beta
    calpha = cos(alpha);
    cbeta = cos(beta);
    salpha = sin(alpha);
    sbeta = sin(beta);

    // get forces and moments 
    Fx = -0.5 * atm.density * V * V * m_Sref * CD * calpha * cbeta + Fx_spin;
    Fy = -0.5 * atm.density * V * V * m_Sref * CD * sbeta + Fy_spin;
    Fz = -0.5 * atm.density * V * V * m_Sref * CD * salpha * cbeta + Fz_spin;
    Mx = 0;
    My = 0;
    Mz = 0;
    // Mx =  0.5 * atm.density * V * V * m_Sref * m_lref * Cl;
    // My =  0.5 * atm.density * V * V * m_Sref * m_lref * Cm;
    // Mz =  0.5 * atm.density * V * V * m_Sref * m_lref * Cn;

    // set outputs
    ans[0] = Fx;
    ans[1] = Fy;
    ans[2] = Fz;
    ans[3] = Mx;
    ans[4] = My;
    ans[5] = Mz;
}

void sphere::get_state_array_delta(double* y0,double* ans){
    // declare input variables
    double u = y0[0];
    double v = y0[1];
    double w = y0[2];
    double p = y0[3];
    double q = y0[4];
    double r = y0[5];
    double xf = y0[6];
    double yf = y0[7];
    double zf = y0[8];
    double e0 = y0[9];
    double ex = y0[10];
    double ey = y0[11];
    double ez = y0[12];
    // declare local variables
    double FM[6];
    double g,Fx,Fy,Fz,Mx,My,Mz;
    double temp1[4], temp2[4],temp3[4];
    

    // get gravity and weight for current state
    g = gravity_english(-zf);
    m_W = g * m_mass;

    // get psudo aerodynamic forces for current states
    this->aerodynamics(y0,FM);
    
    Fx = FM[0];
    Fy = FM[1];
    Fz = FM[2];
    Mx = FM[3];
    My = FM[4];
    Mz = FM[5];
   
    // u v w dots
    ans[0] = g * Fx / m_W + g * 2. * (ex*ez - ey*e0) + (r * v - q * w);
    ans[1] = g * Fy / m_W + g * 2. * (ey*ez + ex*e0) + (p * w - r * u);
    ans[2] = g * Fz / m_W + g * (ez*ez + e0*e0 - ex*ex - ey*ey) + (q * u - p * v);
    // p q r dots
    ans[3] = Mx / m_Ixx;
    ans[4] = My  / m_Iyy;
    ans[5] = Mz  / m_Izz;
    // xf yf zf dots
    // get temps 
    temp1[0] = 0;
    temp1[1] = u;
    temp1[2] = v;
    temp1[3] = w;

    temp2[0] = e0;
    temp2[1] = -ex;
    temp2[2] = -ey;
    temp2[3] = -ez;

    quat_mult(temp1,temp2,temp3);
    quat_mult(&y0[9],temp3,temp1);
    ans[6] = temp1[1];
    ans[7] = temp1[2];
    ans[8] = temp1[3];

    // phi theta psi dots
    ans[9] = 0.5 * (-ex*p - ey*q -ez*r);
    ans[10] = 0.5 * (e0*p - ez*q + ey*r);
    ans[11] = 0.5 * (ez*p + e0*q - ex*r);
    ans[12] = 0.5 * (-ey*p + ex*q + e0*r);
    
 
    
}

void sphere::run_simulation(){
    double t=0;
    int i=0;

    // check to see if file is already open
    const char * file_name = "Sphere_out.txt";
    ifstream f(file_name);
    if (!f.is_open()) {
        cerr << "Error opening file:" << file_name << endl;
    }
    // open file 
    FILE* sphere_out = fopen(file_name,"w");
    fprintf(sphere_out,"Time[s],u[ft/s],v[ft/s],w[ft/s],p[rad/s],q[rad/s],r[rad/s],x[ft],y[ft],z[ft],e0[rad],ex[rad],ey[rad],ez[rad]\n");
   

    quat_norm(&m_state[9]);
    
    // start simulation
    do{
        fprintf(sphere_out,"%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e\n",t,m_state[0],m_state[1],m_state[2],m_state[3],m_state[4],m_state[5],m_state[6],m_state[7],m_state[8],m_state[9],m_state[10],m_state[11],m_state[12]);
        this->rk4();
        if (i>m_i_max){
            cout << "Simulation Stoped at Max Iteration Count:" << i <<endl;
            break;
        }
        quat_norm(&m_state[9]);
        t += m_time_step;
        i++;
    }while(m_state[8]<=0);
    fprintf(sphere_out,"%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e,%20.12e\n",t,m_state[0],m_state[1],m_state[2],m_state[3],m_state[4],m_state[5],m_state[6],m_state[7],m_state[8],m_state[9],m_state[10],m_state[11],m_state[12]);

    fclose(sphere_out);
}

void sphere::rk4(){
    double* k1 = new double[13];
    double* k2 = new double[13];
    double* k3 = new double[13];
    double* k4 = new double[13];
    double y_in[13];
    // calculate k1
    this->get_state_array_delta(m_state,k1);

    // calculate k2
    for (int i = 0; i<13;i++){
        y_in[i] = m_state[i] + k1[i] * m_time_step/2;
    }
    this->get_state_array_delta(y_in,k2);

    //calculate k3
    for (int i = 0; i<13;i++){
        y_in[i] = m_state[i] + k2[i] * m_time_step/2;
    }
    this->get_state_array_delta(y_in,k3);

    // calculate k4
    for (int i = 0; i<13;i++){
        y_in[i] = m_state[i] + k3[i] * m_time_step;
    }
    this->get_state_array_delta(y_in,k4);

    //calculate y next 
    for (int i = 0; i<13; i++) {
        m_state[i] = m_state[i] + (m_time_step/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }
    
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
}
