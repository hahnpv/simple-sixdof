//
//	Homework 2
//	Propogate 4 differential equations forward in time to model 
//	a ballistic trajectory
//

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// struct to hold derivitives and their intermediate RK4 values
struct derivitive {			
	std::vector<double*>x;		// value
	std::vector<double*>xd;		// derivitive
	std::vector<double>x0;		// initial value
	std::vector<double>xd0;		// initial derivitive value
	std::vector<double>xd1;		// first  pass value
	std::vector<double>xd2;		// second pass value
	std::vector<double>xd3;		// third  pass value
};

	//
	//	Inputs
	double V = 5589.0;		// velocity, ft/s
	double theta = 0.006;		// elevation, rad
	//
	//

//
//	The following variables need global scope since we are using functions in stead of classes.
//

double t;					// time
double x, z;				// position
double dx, dz;				
double u, v;				// velocity
double du, dv;
derivitive deriv;			// derivitive struct

void update();		// Update derivitives
void rk4();			// RK4 integrator declaration

// this creates vectors of variables to be integrated by RK4. 
// it also creates the vectors to hold the intermediate RK4 values.
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
};

int main () {


	addToIntegrator(x,dx);		// position
	addToIntegrator(z,dz);

	addToIntegrator(u,du);		// velocity
	addToIntegrator(v,dv);

	// Set IC's
	
	x = 0;
	z = 0;

	dx = V * cos ( theta );
	dz = V * sin ( theta );

	u = dx;
	v = dz;

	t = 0;
	while (t <= 2.0) {

		update();
		rk4();

		cout 	<< t << "\t" 
			<< x << "\t"
			<< z << "\t"
			<< u << "\t"
			<< v << "\t"
			<< endl;
	}

	return 0;
}

void update() {					// update derivitive values

	double rho  = 0.00237847 * pow((1-0.00000688*z),4.256);
	double a    = 49.0124 * sqrt( 518.4 - 0.003566 * z );

	double mass = 7.0/32.2;

	double V    = sqrt( dx*dx + dz*dz ); 
	
	double M    = V / a;
	
	double CD;
	if ( M > 6 ) 
		CD = 1.2;
	else if ( M > 4)
		CD = 1.2;
	else if ( M > 3.5)
		CD = 1.31;
	else if ( M > 3.0)
		CD = 1.42;
	else if ( M > 2.5)
		CD = 1.68;
	else if ( M > 2.0)
		CD = 1.90;
	else if ( M > 1.75)
		CD = 2.12;
	else if ( M > 1.5)
		CD = 2.64;
	else if ( M > 1.35)
		CD = 3.16;

	double beta = (3.14159/8) * rho * V * pow(0.1243,2) * CD / mass;

	dx = u;
	du = -beta * dx;
	dz = v;
	dv = -32.2 - beta*dz;
}

// runge-Kutta 4th order integrator.
void rk4() {
	double dt = 0.001;				// CLOCK RATE

	derivitive *d;

	for (unsigned int pass=0; pass<=3; pass++) {

		d = &deriv;

		for(unsigned int i=0; i < d->x.size(); i++) {
			if ( pass == 0 ) {
				d->x0[i]  = *d->x[i];
				d->xd0[i] = *d->xd[i];
				*d->x[i]  =  d->x0[i] + dt / 2.0 * d->xd0[i];
			} else if ( pass == 1 ) {
				d->xd1[i] = *d->xd[i];
				*d->x[i]  =  d->x0[i] + dt / 2.0 * d->xd1[i];
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
