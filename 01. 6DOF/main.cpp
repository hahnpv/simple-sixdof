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
	double m    = 5.029;
	double rho  = 1.225;
	Vec3    g( 0, 0, 9.81);
	double D    = 0.0268;
	double Ixx  = 0.00033;
	double Iyy  = 0.24021;
	double Izz  = 0.24021;
	double phi, theta, psi;
	//
	//

void do_aeros();	// aerodynamic equations	
void update();		// Update derivitives
void rk4();			// RK4 integrator declaration
void addToIntegrator(double &x, double &dx); 	// put diffeq pairs into our derivitive struct

struct derivitive {				// struct to hold derivitives and their intermediate RK4 values
	std::vector<double*>x;		// value
	std::vector<double*>xd;		// derivitive
	std::vector<double>x0;		// initial value
	std::vector<double>xd0;		// initial derivitive value
	std::vector<double>xd1;		// first  pass value
	std::vector<double>xd2;		// second pass value
	std::vector<double>xd3;		// third  pass value
};
double t;
derivitive deriv;			// derivitive struct

	Vec3 dphidthetadpsi;		// Vector declatartions of states and derivitives
	Vec3 xyz;
	Vec3 dxdydz;
	Vec3 uvw;
	Vec3 dudvdw;
	Vec3 pqr;
	Vec3 dpdqdr;

	Vec3 XYZ;
	Vec3 LMN;

	Mat3 iframe;

	Table *Cx0tab;			// Table object declarations
	Table *Cx2tab;
	Table *Cy0tab;
	Table *Cz0tab;
	Table *CnAtab; 
	Table *CypAtab;
	Table *Clptab;
	Table *Cmqtab;
	Table *Clddtab;
	Table *SLDEL;

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
	Clddtab = new Table(fname, "[Cldd]");
	SLDEL = new Table(fname, "[SLDEL]");

	// Set IC's

	xyz(0,0,0);
	uvw(1676,0,0);
	pqr(0,0,0);

	phi   = 0;
	theta = 0.01745;
	psi   = 0;

	iframe( Ixx, 0, 0,
		0, Iyy, 0,
		0, 0, Izz);

	// Simulation Loop 
	t = 0;
	cout << " time x y z negz range u v w V p q r phi theta psi" << endl;
	while (t < 2.0) {
		update();
		rk4();
		cout 	<< t << " " 
			<< xyz.x << " "
			<< xyz.y << " "
			<< xyz.z << " "
			<< -xyz.z << " "
			<< xyz.mag() << " "
			<< uvw.x << " " 
			<< uvw.y << " "
			<< uvw.z << " "
			<< uvw.mag() << " "
			<< pqr.x << " "
			<< pqr.y << " "
			<< pqr.z << " "
			<< phi << " "
			<< theta << " "
			<< psi 
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

	do_aeros();

	dxdydz = T_b_e * uvw;

	dphidthetadpsi = eulerframe * pqr;

	dudvdw = transframe * uvw + XYZ + T_b_e.inv() * g;

	dpdqdr = iframe.inv()*transframe*iframe*pqr + iframe.inv()*LMN;

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

void rk4() {
	double dt = 0.001;				// CLOCK RATE

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


