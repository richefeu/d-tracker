#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

struct data_t
{
	double x,y,dx,dy;
	double dxmean,dymean;
};

vector <data_t> data;

int equals( double a, double b, double tolerance )
{
	return ( a == b ) ||
		( ( a <= ( b + tolerance ) ) &&
		( a >= ( b - tolerance ) ) );
}

double cross2( double x0, double y0, double x1, double y1 )
{
	return x0*y1 - y0*x1;
}


int in_range( double val, double range_min, double range_max, double tol )
{
	return ((val+tol) >= range_min) && ((val-tol) <= range_max);
}


// Returns number of solutions found.  If there is one valid solution, it will be put in s and t 
//   p2 --- p3
//   |      |
// t |  p   |
//   |      |
//   p0 --- p1
//      s
int inverseBilinear( double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double* sout, double* tout, double* s2out, double* t2out )
{
	int t_valid, t2_valid;

	double a  = cross2( x0-x, y0-y, x0-x2, y0-y2 );
	double b1 = cross2( x0-x, y0-y, x1-x3, y1-y3 );
	double b2 = cross2( x1-x, y1-y, x0-x2, y0-y2 );
	double c  = cross2( x1-x, y1-y, x1-x3, y1-y3 );
	double b  = 0.5 * (b1 + b2);

	double s, s2, t, t2;

	double am2bpc = a-2*b+c;
	// this is how many valid s values we have
	int num_valid_s = 0;

	if ( equals( am2bpc, 0, 1e-10 ) ) {
		if ( equals( a-c, 0, 1e-10 ) ) {
		// Looks like the input is a line
		// You could set s=0.5 and solve for t if you wanted to
			return 0;
		}
		s = a / (a-c);
		if ( in_range( s, 0, 1, 1e-10 ) )
			num_valid_s = 1;
	}
	else {
		double sqrtbsqmac = sqrt( b*b - a*c );
		s  = ((a-b) - sqrtbsqmac) / am2bpc;
		s2 = ((a-b) + sqrtbsqmac) / am2bpc;
		num_valid_s = 0;
		if ( in_range( s, 0, 1, 1e-10 ) ) {
			num_valid_s++;
			if ( in_range( s2, 0, 1, 1e-10 ) )
				num_valid_s++;
		}
		else {
			if ( in_range( s2, 0, 1, 1e-10 ) ) {
				num_valid_s++;
				s = s2;
			}
		}
	}

	if ( num_valid_s == 0 )
		return 0;

	t_valid = 0;
	if ( num_valid_s >= 1 ) {
		double tdenom_x = (1-s)*(x0-x2) + s*(x1-x3);
		double tdenom_y = (1-s)*(y0-y2) + s*(y1-y3);
		t_valid = 1;
		if ( equals( tdenom_x, 0, 1e-10 ) && equals( tdenom_y, 0, 1e-10 ) ) {
			t_valid = 0;
		}
		else {
			// Choose the more robust denominator
			if ( fabs( tdenom_x ) > fabs( tdenom_y ) ) {
				t = ( (1-s)*(x0-x) + s*(x1-x) ) / ( tdenom_x );
			}
			else {
				t = ( (1-s)*(y0-y) + s*(y1-y) ) / ( tdenom_y );
			}
			if ( !in_range( t, 0, 1, 1e-10 ) )
				t_valid = 0;
		}
	}

	// Same thing for s2 and t2
	t2_valid = 0;
	if ( num_valid_s == 2 ) {
		double tdenom_x = (1-s2)*(x0-x2) + s2*(x1-x3);
		double tdenom_y = (1-s2)*(y0-y2) + s2*(y1-y3);
		t2_valid = 1;
		if ( equals( tdenom_x, 0, 1e-10 ) && equals( tdenom_y, 0, 1e-10 ) ) {
			t2_valid = 0;
		}
		else {
			// Choose the more robust denominator
			if ( fabs( tdenom_x ) > fabs( tdenom_y ) ) {
				t2 = ( (1-s2)*(x0-x) + s2*(x1-x) ) / ( tdenom_x );
			}
			else {
				t2 = ( (1-s2)*(y0-y) + s2*(y1-y) ) / ( tdenom_y );
			}
			if ( !in_range( t2, 0, 1, 1e-10 ) )
				t2_valid = 0;
		}
	}

	// Final cleanup
	if ( t2_valid && !t_valid ) {
		s = s2;
		t = t2;
		t_valid = t2_valid;
		t2_valid = 0;
	}

	// Output
	if ( t_valid ) {
		*sout = s;
		*tout = t;
	}

	if ( t2_valid ) {
		*s2out = s2;
		*t2out = t2;
	}

	return t_valid + t2_valid;
}

void bilinear( double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double s, double t, double* x, double* y )
{
	*x = t*(s*x3+(1-s)*x2) + (1-t)*(s*x1+(1-s)*x0);
	*y = t*(s*y3+(1-s)*y2) + (1-t)*(s*y1+(1-s)*y0);
}

void read_field (string filename)
{
	if (!data.empty()) data.clear();
	ifstream file(filename.c_str());

	data_t D;
	while (file) {
		file >> D.x >> D.y >> D.dx >> D.dy; // Format a voir avec Gael
		data.push_back(D);
	}
}


// A revoir...
int main(int argc, char ** argv)
{	
	string ext = ".nod"; // default extension

	if (argc != 5) { cout << "usage: " << argv[0] << " basename num_ini num_fin overwrite" << endl; return 0; }
	string basename = argv[1];
	unsigned int num_ini = atoi(argv[2]);
	unsigned int num_fin = atoi(argv[3]);
	//unsigned int step = atoi(argv[4]);
	unsigned int overwrite = atoi(argv[4]);
	// if overwrite = 1 we will overwrite the file on the previous one, else we will create a new one
	cout << "the chosen initial number is" << num_ini << endl;
	cout << "the chosen final number is" << num_fin << endl;

	//unsigned int num_ini = 2955581,num_fin = 2955607;
	//unsigned int num=2955584;
	int nok;
	double s,t,s2,t2;

	for (unsigned int num = num_ini; num <= num_fin ;num++)
	{
		stringstream name_in;
		name_in << basename << num << ext; 
		read_field (name_in.str());

		stringstream name_out;
		name_out << basename << num; 
		if (!overwrite) name_out << "_mod";
		name_out << ext;

		cout << name_in.str() << " --> " << name_out.str() << endl;

		ofstream out(name_out.str().c_str());

		for (int i = 0 ; i < 4 ; i++) {
			data[i].dxmean = data[i].dx;
			data[i].dymean = data[i].dy;
		}

		for (int i = 4 ; i < data.size() ; i++) {
			nok = inverseBilinear( data[0].x, data[0].y, data[3].x, data[3].y, 
				data[1].x, data[1].y, data[2].x, data[2].y, data[i].x, data[i].y, &s,&t,&s2,&t2);
			bilinear( data[0].dx, data[0].dy, data[3].dx, data[3].dy,
				data[1].dx, data[1].dy, data[2].dx, data[2].dy, s,t, &(data[i].dxmean), &(data[i].dymean) );
			out << data[i].x << "\t" << data[i].y << "\t" 
				<< data[i].dx << "\t" << data[i].dy << "\t" 
				<< data[i].dxmean << "\t" << data[i].dymean << endl;	
		}
	}

}





