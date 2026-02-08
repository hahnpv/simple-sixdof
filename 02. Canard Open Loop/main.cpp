//
//	Homework 3
//	Full 6DOF
//

#include <iostream>
#include <vector>
#include <cmath>

#include "tbl.h"
#include "vec3.h"
#include "mat3.h"

const double PI = 3.14159;
using namespace std;
using namespace dsf::util;

	//
	//	Inputs
	double m;		// lookup table
	double rho  = 1.225;
	Vec3    g( 0, 0, 9.81);
	double D    = 0.0706;		
	double phi, theta, psi;
	//
	//

double dt = 0.001;
bool safe_to_sample;	// don't let guidance iterate during RK integration
void do_aeros();	// aerodynamic equations
void do_rockets();	// rocket power	
void update();		// Update derivitives
void rk4();		// RK4 integrator declaration
void addToIntegrator(double &x, double &dx); 	// put diffeq pairs into our derivitive struct

struct derivitive {			// struct to hold derivitives and their intermediate RK4 values
	std::vector<double*>x;		// value
	std::vector<double*>xd;		// derivitive
	std::vector<double>x0;		// initial value
	std::vector<double>xd0;		// initial derivitive value
	std::vector<double>xd1;		// first  pass value
	std::vector<double>xd2;		// second pass value
	std::vector<double>xd3;		// third  pass value
};
double t;				// sim time
derivitive deriv;			// derivitive struct

	Vec3 dphidthetadpsi;
	Vec3 xyz;
	Vec3 dxdydz;
	Vec3 uvw;
	Vec3 dudvdw;
	Vec3 pqr;
	Vec3 dpdqdr;

	Vec3 XYZ;
	Vec3 LMN;

	Vec3 XYZ_c;	// canard
	Vec3 LMN_c;

	Vec3 XYZ_r;	// rocket 
	Vec3 LMN_r;

	Mat3 iframe;

	Table  *Cx0tab;
	Table  *Cx2tab;
	Table  *Cy0tab;
	Table  *Cz0tab;
	Table  *CnAtab; 
	Table  *CypAtab;
	Table  *Clptab;
	Table  *Cmqtab;
	Table  *Clddtab;
	Table  *SLDEL;
	Table  *Ixx;
	Table  *Iyy;
	Table  *slcg;		// control system
	Table  *thrust;
	Table  *mass;

	// canard transform matrix initialization

	vector<Mat3>T_c;
	vector<Vec3>sl;

	// control system
	double d[4] = {0,0,0,0};	// 	Canard Deflection

	double integral_error;		//	Roll Authority
	double error;			//

	double angle_integral_error;	//	Angular Authority	
	double angle_error;		//
	
	double position_integral_error;	//	Position Authority
	double position_error;		//
	double position_derivitive_error;
	double y_old = 0;

int main () {

	addToIntegrator(phi,   dphidthetadpsi.x);
	addToIntegrator(theta, dphidthetadpsi.y);
	addToIntegrator(psi,   dphidthetadpsi.z);

	addToIntegrator(xyz.x,dxdydz.x);
	addToIntegrator(xyz.y,dxdydz.y);
	addToIntegrator(xyz.z,dxdydz.z);

	addToIntegrator(uvw.x,dudvdw.x);
	addToIntegrator(uvw.y,dudvdw.y);
	addToIntegrator(uvw.z,dudvdw.z);

	addToIntegrator(pqr.x,dpdqdr.x);
	addToIntegrator(pqr.y,dpdqdr.y);
	addToIntegrator(pqr.z,dpdqdr.z);	

	char *fname = "6dofTables.dat";
	Cx0tab = new Table(fname,"[Cx0]");
	Cx2tab = new Table(fname,"[Cx2]");
	Cy0tab = new Table(fname,"[Cy0]");
	Cz0tab = new Table(fname,"[Cz0]");
	CnAtab = new Table(fname,"[CnA]"); 
	CypAtab = new Table(fname,"[CypA]");
	Clptab = new Table(fname,"[Clp]");
	Cmqtab = new Table(fname,"[Cmq]");
	Clddtab   = new Table(fname,"[Cldd]");
	SLDEL  = new Table(fname,"[SLDEL]");
	Ixx    = new Table(fname,"[Ixx]");
	Iyy    = new Table(fname,"[Iyy]");
	slcg   = new Table(fname,"[SLCG]");
	thrust = new Table(fname,"[THRUST]");
	mass   = new Table(fname,"[MASS]");

	// Canard generation
	double phi_c = 0;			// rotation about x-axis of body
	double gamma_c = 0;			// tilt forward/back (null for now)
		Mat3 y(  cos(gamma_c),-sin(gamma_c),0, 
			cos(phi_c)*sin(gamma_c),cos(phi_c)*cos(gamma_c),-sin(phi_c), 
			sin(phi_c)*sin(gamma_c),sin(phi_c)*cos(gamma_c),cos(phi_c));
		T_c.push_back( y );	// [0]
		T_c.push_back( y );	// placeholder
		phi_c += PI;
		gamma_c += 0;
		y = Mat3(  cos(gamma_c),-sin(gamma_c),0, 
			cos(phi_c)*sin(gamma_c),cos(phi_c)*cos(gamma_c),-sin(phi_c), 
			sin(phi_c)*sin(gamma_c),sin(phi_c)*cos(gamma_c),cos(phi_c));

		T_c.push_back( y );	// [2]
		T_c.push_back( y );	// placeholder


	// hard code sl values
	Vec3 x;
	x = Vec3(0.05,0.02,0);
	sl.push_back( x );
	x = Vec3(0,0,0);
	sl.push_back( x );		// not using 1 for now
	x = Vec3(0.05,-0.02,0);
	sl.push_back( x );
	x = Vec3(0,0,0);			// not using 3 for now
	sl.push_back( x );

	// Set IC's

	xyz(0,0,0);				// position
	uvw(38,0,0);				// velocity
	pqr(0,0,0);				// rotation rates

	phi   = 0;				// orientation
	theta = 0.1309;				//
	psi   = 0;				//


	// Simulation Loop 
	t = 0;
	cout << " time x y z negz range u v w V p q r phi theta psi" << endl;
	while (t < 5.0) {
		update();
		rk4();
		cout 	<< t << " " 
			<< xyz.x << " "
			<< xyz.y << " "
			<< xyz.z << " "
//			<< -xyz.z << " "
//			<< xyz.mag() << " "
//			<< uvw.x << " " 
//			<< uvw.y << " "
//			<< uvw.z << " "
			<< "\t"
			<< uvw.mag() << " "
			<< "\t\t" 
			<< pqr.x << " "
//			<< pqr.y << " "
//			<< pqr.z << " "
			<< phi << " "
//			<< theta << " "
//			<< psi 
//			<< "\t"
//			<< XYZ_c << " "
//			<< LMN_c 
			<< endl;

	}

	return 0;
}

void update() {

	Mat3 transframe(   0,  pqr.z, -pqr.y,
			 -pqr.z,  0,  pqr.x,
			  pqr.y, -pqr.x,  0);

	Mat3 eulerframe( 1, sin(phi)*tan(theta), cos(phi)*tan(theta),
			0, cos(phi), -sin(phi),
			0, sin(phi)/cos(theta), cos(phi)/cos(theta));

	Mat3 T_b_e( 	cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi),
			cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi),
			-sin(theta),					  sin(phi)*cos(theta),				   cos(phi)*cos(theta));


	iframe( Ixx->interp(t), 0, 0,		
		0, Iyy->interp(t), 0,
		0, 0, Iyy->interp(t));

	m = mass->interp(t);

	do_aeros();
	do_rockets();

	dxdydz = T_b_e * uvw;

	dphidthetadpsi = eulerframe * pqr;

	dudvdw = transframe * uvw + XYZ + XYZ_c/m + XYZ_r/m + T_b_e.inv() * g;

	dpdqdr = iframe.inv()*transframe*iframe*pqr + iframe.inv()*( LMN + LMN_c + LMN_r );

}

void do_aeros(void) {

	double M = uvw.mag() / sqrt( 1.4 * 287 * 300 );
	double Cx0, Cx2, Cy0, Cz0, CnA, CypA, Clp, Cmq, sl, Cldd;

	Cx0  = Cx0tab->interp(M);
	Cx2  = Cx2tab->interp(M);
	Cy0  = Cy0tab->interp(M);
	Cz0  = Cz0tab->interp(M);
	CnA  = CnAtab->interp(M); 
	CypA = CypAtab->interp(M);
	Clp  = Clptab->interp(M);  
	Cmq  = Cmqtab->interp(M);
	Cldd = Clddtab->interp(M);
	sl   = SLDEL->interp(M);

	double V = uvw.mag();

	double v = uvw.y;
	double w = uvw.z;

	double p = pqr.x;
	double q = pqr.y;
	double r = pqr.z;

	double qS  = (PI/8)*rho*pow(V,2)*pow(D,2);
	double qSd = (PI/8)*rho*pow(V,2)*pow(D,3);

		// Normal Aero Forces
	double Xa, Ya, Za;
	Xa = -qS*(Cx0 + Cx2*( ( pow(v,2) + pow(w,2) ) /pow(V,2) ) );
	Ya = -qS*(Cy0 + CnA*(v/V) - ((p * D) / (2 * V))*CypA*(w/V));
	Za = -qS*(Cz0 + CnA*(w/V) + ((p * D) / (2 * V))*CypA*(v/V));
	
	double La, Ma, Na;
	La = qSd*( Cldd + (p*D*Clp)/(2*V) );
	Ma = qSd*( (q*D*Cmq)/(2*V) );
	Na = qSd*( (r*D*Cmq)/(2*V) );	

	LMN = Vec3( La, Ma - sl*Za, Na + sl*Ya);

	XYZ   ( Xa/m, Ya/m, Za/m);
}

void do_rockets(void) {

	double theta   = 1.97/57.294;	// not 100% sure which is which...
	double psi = 0;		// 
	Mat3 T_r( cos(psi)*cos(theta), -sin(psi)*cos(theta), sin(theta),
		sin(psi), cos(psi), 0,
		-cos(psi)*sin(theta), sin(psi)*sin(theta), cos(theta)		// typo in presentation
	);

	XYZ_r = T_r * Vec3( thrust->interp(t), 0, 0);
	LMN_r = Vec3( 1.997 - slcg->interp(t), 0.0254, 0).cross( XYZ_r );

//	cout << "Rocket: " << XYZ_r << " " << LMN_r <<  endl;
}

void rk4() {
	safe_to_sample = false;
	derivitive *d;

	for (unsigned int pass=0; pass<=3; pass++) {

		d = &deriv;

		// for each derivitive
		for(unsigned int i=0; i < d->x.size(); i++) {
			if ( pass == 0 ) {
				d->x0[i]  = *d->x[i];
				d->xd0[i] = *d->xd[i];
				*d->x[i]  =  d->x0[i] + dt / 2 * d->xd0[i];
			} else if ( pass == 1 ) {
				d->xd1[i] = *d->xd[i];
				*d->x[i]  =  d->x0[i] + dt / 2 * d->xd1[i];
			} else if ( pass == 2 ) {
				d->xd2[i] = *d->xd[i];
				*d->x[i]  =  d->x0[i] + dt * d->xd2[i];
			} else {
				d->xd3[i] = *d->xd[i];
				*d->x[i]  =  d->x0[i] + dt / 6 * ( d->xd0[i] + 2 * d->xd1[i] + 2 * d->xd2[i] + d->xd3[i] );
			}
		}
	if (pass == 0 || pass == 2)
		t += dt/2;

	update();
	}
	safe_to_sample = true;
}

void addToIntegrator(double &x, double &dx) {
	double *x_ptr, *dx_ptr;
	 x_ptr =  &x;
	dx_ptr = &dx;
		
 	 deriv.x.push_back(  x_ptr );
	deriv.xd.push_back( dx_ptr );

	deriv.x0.push_back ( 0 );
	deriv.xd0.push_back( 0 );
	deriv.xd1.push_back( 0 );
	deriv.xd2.push_back( 0 );
	deriv.xd3.push_back( 0 );
}


