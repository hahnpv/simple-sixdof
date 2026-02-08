//
//	Homework 3
//	Full 6DOF
//

#include <iostream>
#include <iomanip>
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
void do_gnc();		// canard GNC
void do_pulsejets();	// pulseget GNC
void update();		// Update derivitives
void rk4();		// RK4 integrator declaration
void addToIntegrator(double &x, double &dx); 	// put diffeq pairs into our derivitive struct

double getUniform(double min, double max) {
	return min + (max - min) * ( (double) rand() / RAND_MAX);
}

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

	Table *Cx0tab;
	Table *Cx2tab;
	Table *Cy0tab;
	Table *Cz0tab;
	Table *CnAtab; 
	Table *CypAtab;
	Table *Clptab;
	Table *Cmqtab;
	Table *Clddtab;
	Table *SLDEL;
	Table *Ixx;
	Table *Iyy;
	Table *slcg;		// control system
	Table *thrust;
	Table *mass;

	// canard transform matrix initialization
	vector<Mat3>T_c;
	vector<Vec3>sl;

	// control system
	double d[4] = {0,0,0,0};	// 	Canard Deflection

	double integral_error;		//	Roll Authority

	double angle_integral_error;	//	Angular Authority	

	double position_integral_error;	//	Position Authority

	Vec3 IC;
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

	// set the random seed
	srand( (unsigned)time(NULL) );

	// Canard generation
	double phi_c = 0;			// rotation about x-axis of body
	double gamma_c = 0;			// tilt forward/back (null for now)
		Mat3 y(  cos(gamma_c),-sin(gamma_c),0, 
			cos(phi_c)*sin(gamma_c),cos(phi_c)*cos(gamma_c),-sin(phi_c), 
			sin(phi_c)*sin(gamma_c),sin(phi_c)*cos(gamma_c),cos(phi_c));
		T_c.push_back( y );	// [0]
		T_c.push_back( y );	//	junk
		phi_c += PI;
		gamma_c += 0;
		y = Mat3(  cos(gamma_c),-sin(gamma_c),0, 
			cos(phi_c)*sin(gamma_c),cos(phi_c)*cos(gamma_c),-sin(phi_c), 
			sin(phi_c)*sin(gamma_c),sin(phi_c)*cos(gamma_c),cos(phi_c));

		T_c.push_back( y );	// [2] phi_c = pi, gamma_c = 0
		T_c.push_back( y );	// 	junk


	// hard code sl values		-> NEEDS A REWRITE ... cg = 1.217m after t=1.10s. Max distance from tail for control=1.3m
	//				-> so del-x_max = 0.083
	//				-> del-y = 0.0706/2 = 0.0353

	// try giving x component negative sense ...
	Vec3 x;
	x = Vec3(0.08,0.0453,0);	// was: 0.05, 0.02, should be 0.083, 0.0353, 0 ... value mods severity but doesn't negate
	sl.push_back( x );	// second value should be midpoint of the airfoil ( 0.0353 + span/2 )
	x = Vec3(0,0,0);	
	sl.push_back( x );
	x = Vec3(0.08,-0.0453,0);
	sl.push_back( x );
	x = Vec3(0,0,0);	
	sl.push_back( x );

	// Set IC's

	xyz(0,0,0);				// position
	uvw(38,0,0);				// velocity
//	pqr(0,0,0);				// rotation rates
	pqr(0,getUniform(-0.15,0.15),getUniform(-0.15,0.15));		// DISTURBED case
	// REQUIREMENTS CHANGED:
	// STDEV = 0.15 rad/s for p,q

	cout << "pqr: " << pqr << endl;

	IC = pqr;

	phi   = 0;				// orientation
	theta = 0.1309;			//				// SHOULD BE 0.1309
	psi   = 0;				//

	// Simulation Loop 
	t = 0;
	int i = 1;
	cout << "t x y z V p q r phi theta psi  d[0] d[2] AoA" << endl;
	while (t < 10) {
		update();
		rk4();
		if ( i == 10 ) {
		cout 	<< setw(5) << setprecision(4) << t << " " 
			<< xyz << "\t"
			<< uvw.mag() << "\t"
			<< pqr << "\t"
			<< phi   << " "
			<< theta << " "
			<< psi   << "\t"
			<< d[0] << ","
			<< d[2] << "\t"
			<< atan2( uvw.z, uvw.x ) * 57.295
			<< endl;
			i = 0;
		}
		i++;
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

	Mat3 T_e_b( 	cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi),
			cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi),
			-sin(theta),					  sin(phi)*cos(theta),				   cos(phi)*cos(theta));

	iframe( Ixx->interp(t), 0, 0,		
		0, Iyy->interp(t), 0,
		0, 0, Iyy->interp(t));

	m = mass->interp(t);

	do_aeros();

	if (t <=1.1) {
		do_rockets();
	} else {
		XYZ_r = Vec3(0,0,0);
		LMN_r = Vec3(0,0,0);
	}

	if (t >= 1.25)
		do_gnc();

	//
	//	Undisturbed graphs match almost perfectly.
	//	(theta doesn't jump like his graph ... )
	//

	dxdydz = T_e_b * uvw;

	dphidthetadpsi = eulerframe * pqr;

	dudvdw = transframe * uvw + T_e_b.inv() * g + XYZ/m + XYZ_c/m + XYZ_r/m;

	dpdqdr = iframe.inv()*transframe*iframe*pqr + iframe.inv()*( LMN  + LMN_c + LMN_r );
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
	Ya = -qS*(Cy0 + CnA*(v/V) - ((p * D) / (2 * V))*CypA*(w/V));	// CypA = 0
	Za = -qS*(Cz0 + CnA*(w/V) + ((p * D) / (2 * V))*CypA*(v/V));	// CypA = 0

	double La, Ma, Na;
	La = qSd*( Cldd + (p*D*Clp)/(2*V) );
	Ma = qSd*( (q*D*Cmq)/(2*V) );
	Na = qSd*( (r*D*Cmq)/(2*V) );	

	XYZ(Xa, Ya, Za);
	LMN = Vec3(La, Ma, Na) + Vec3(-sl,0,0).cross(XYZ);	// MINUS sl not +sl!!!! was screwing this up ... 
								// pitch is correct now
}
void do_rockets(void) {

	double theta = -1.97/57.294;
	double psi   = 0;
	Mat3 T_r_1( cos(psi)*cos(theta), -sin(psi)*cos(theta), sin(theta),
		 sin(psi), cos(psi), 0,
		-cos(psi)*sin(theta), sin(psi)*sin(theta), cos(theta)
	);
	theta = +1.97/57.294;
	Mat3 T_r_2( cos(psi)*cos(theta), -sin(psi)*cos(theta), sin(theta),
		 sin(psi), cos(psi), 0,
		-cos(psi)*sin(theta), sin(psi)*sin(theta), cos(theta)
	);

	Vec3 XYZ_r1, XYZ_r2, LMN_r1, LMN_r2;
	XYZ_r1 = T_r_1 * Vec3( thrust->interp(t)/2, 0, 0);
	XYZ_r2 = T_r_2 * Vec3( thrust->interp(t)/2, 0, 0);
	LMN_r1 = Vec3( 1.997 - slcg->interp(t),  0.0254, 0).cross( XYZ_r1 ); 
	LMN_r2 = Vec3( 1.997 - slcg->interp(t), -0.0254, 0).cross( XYZ_r2 );

	XYZ_r = XYZ_r1 + XYZ_r2;
	LMN_r = LMN_r1 + LMN_r2;
}

void do_gnc(void) {
	double max_canard_angle = 0.174;	// 8 degrees (0.1396)=146
						// Baseline Case WORKS!!!

	if ( t >= 1.25 && safe_to_sample) {				// Despin Control (working)
		double error = pqr.x;

		if ( integral_error > 0.5 )				// integral error limiters
			integral_error = 0.5;
		else if (integral_error < -0.5 )
			integral_error = -0.5;

		integral_error += dt*pqr.x;

		double asym_control = 0.07 * error + 0.1 * integral_error; 
		//		      0.07	     0.1	baseline case (CL=1.83)

		if ( asym_control >= max_canard_angle ) {		// physical canard limiters
			d[0] = max_canard_angle;
			d[2] = max_canard_angle;
		} else if ( asym_control < -max_canard_angle ) {
			d[0] = -max_canard_angle;
			d[2] = -max_canard_angle;
		} else {
			d[0] =  asym_control;
			d[2] =  asym_control;
		}
	}

	if ( t >= 2.1 && safe_to_sample ) {				// Roll Control (mod into range control)
		double angle_error = phi - 2*PI*floor(phi/(2*PI));

		if ( angle_integral_error > 0.5 )			// integral error limiters
			angle_integral_error = 0.5;
		else if (angle_integral_error < -0.5 )
			angle_integral_error = -0.5;

		angle_integral_error +=  dt * angle_error ;

		// 0.15, 0.325 w/CL = 1.5
		double control = 0.05*angle_error + 0.145*angle_integral_error;
		//		 0.05		    0.22	baseline case (CL=1.83) (old drag)

		//		 0.05		    0.145 	baseline case, new drag
		
		d[0] += control;
		d[2] += control;

		if ( d[0] >= max_canard_angle) {			// physical canard limiters
			d[0] = max_canard_angle;
			d[2] = max_canard_angle;
		} else if ( d[0] < -max_canard_angle ) {
			d[0] = -max_canard_angle;
			d[2] = -max_canard_angle;
		}
	}

	if ( t >= 3.3 && safe_to_sample ) {		// Divert Control 
		double position_error = xyz.y; 

		if ( position_integral_error > 0.5 )
			position_integral_error = 0.5;
		else if ( position_integral_error < -0.5)
			position_integral_error = -0.5;

		position_integral_error += dt * xyz.y;

		
		// works for positive values
		double slope, y;
		if (abs(IC.z) < 0.05){
			slope = 3.44;
			y = slope * abs(IC.z) + 0.0105; 
		} else {
			slope = -1.25;
			y = 0.1825 + slope*(abs(IC.z)-0.05);
		}

		// adjust the period on this one till t=10 and call it good.
		double divert = 0.005*position_error + y *position_integral_error;
		//		0.005		       0.115	baseline case (CL=1.83) (old drag)

		// 0,0,0	0.005			0.0105	first  swerve
		// 0,0,0	0.005			0.23	second swerve
		// 0,0,0.025	0.005			0.0885	first swerve
		// 0,0,0.05	0.005			0.1825	first swerve	
		// 0,0,0.075	0.005			0.152	first swerve
		// 0,0,0.1	0.005			0.12	first swerve
		// 0,0,0.125	0.005			0.0875	first swerve
		// 0,0,0.15	0.005			0.0575	first swerve

		d[0] -= divert;
		d[2] += divert;

		if ( d[0] >= max_canard_angle )
			d[0] = max_canard_angle;
		else if ( d[0] < -max_canard_angle )
			d[0] = -max_canard_angle;
		
		if ( d[2] >= max_canard_angle )
			d[2] = max_canard_angle;
		else if ( d[2] < -max_canard_angle )
			d[2] = -max_canard_angle;
		}

	double a[4];		// angle of attack
	double V[4];		// velocity wrt canards
	double M[4];		// Mach, wrt canards
	double qS[4];		// dynamic pressure
	double S = 0.003;	// surface area

	double CLa = 1.83;
	double CL[4];
	// Need to find velocities wrt each fin to calculate alpha
	for (int i=0; i<4; i++)	{
		Mat3 transframe (     0, pqr.z,  pqr.y,		// sleger's slides. 
				-pqr.z,     0, -pqr.z,		// implement WITHOUT  " * -1 "
				-pqr.y, pqr.x,     0);		// might fix something ...

		Vec3 uvw_c = T_c[i].inv() * (uvw + transframe*sl[i] );
		a[i] = d[i] + atan2(uvw_c.z, uvw_c.x);
		CL[i] = CLa * a[i];
		V[i] = uvw_c.mag();
		M[i] = V[i] / sqrt( 1.4 * 287 * 300 );
		qS[i] = 0.5 * rho * pow ( V[i] , 2 )*S;
	}

	double AR = 0.83;	// Aspect ratio
	double e  = 0.99;	// Oswald efficiency factor,  1/(1+delta) 
	double k = 1/(PI * AR * e);
	double CD[4] = { 0.02 + k * pow(CL[0],2) , 0.02, 0.02 + k * pow(CL[0],2), 0.02 };

	XYZ_c = Vec3(0,0,0);
	LMN_c = Vec3(0,0,0);
	for (int i=0; i<4; i+=2) {
		Vec3 xyz = T_c[i] * Vec3( CL[i]*sin(a[i]-d[i]) - CD[i]*cos(a[i]-d[i]), 0, -CL[i]*cos(a[i]-d[i]) - CD[i]*sin(a[i]-d[i])) * qS[i];
		XYZ_c += xyz;
		LMN_c += sl[i].cross( xyz );
	}
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
