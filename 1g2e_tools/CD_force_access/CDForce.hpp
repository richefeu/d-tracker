#pragma once
// Version avec uniquement fn et ft
// author: Vincent Richefeu
// year: 2011

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

typedef unsigned int uint;
typedef double real;

// STRUCTURES
struct grain
{
	real x,y,rot;    // position
	real vx,vy,vrot; // velocity
	real fx,fy,M;    // resultant force/moment
	real mass,mom;   // mass and inertial moment
	real rad;        // radius (for display)
	
	grain():fx(0.0),fy(0.0),M(0.0) { }
};

struct control
{
	bool vxImposed;
	bool vyImposed;
	bool vrotImposed;
	real xvalue,yvalue,rotvalue;
};


struct contact
{
	uint idi,idj;
	grain *i,*j;
	real x,y;
	real nx,ny;
	real vn,vt;
	real fn,ft;
	
	real invmi,invmj,invmomi,invmomj;
	real fac_en,fac_et;
	
	real Wnn,Wtt,Wnt;
	real invWnn,invWtt,invWnt; // Effective mass = inverse of W matrix

	void Res() // transfer the contact forces (and moment) to the bodies 
	{	
		real tx = -ny, ty = nx;
		real fx = fn * nx + ft * tx;
		real fy = fn * ny + ft * ty;
		real cix, ciy, cjx, cjy;

		cix = x - i->x;
		ciy = y - i->y; 

		cjx = x - j->x;
		cjy = y - j->y;

		i->fx   += fx;
		i->fy   += fy;
		i->M += (cix * fy - ciy * fx);  // +fs 

		j->fx   -= fx;
		j->fy   -= fy;
		j->M += (-cjx * fy + cjy * fx); // +fs
	}
	
	void Res(const real dfn, const real dft/*, const real dfs*/) // Increment de resultant force of the bodies in contact
	{	
		real tx = -ny, ty = nx;
		real dfx = dfn * nx + dft * tx;
		real dfy = dfn * ny + dft * ty;
		real cix, ciy, cjx, cjy;

		cix = x - i->x;
		ciy = y - i->y; 

		cjx = x - j->x;
		cjy = y - j->y;

		i->fx   += dfx;
		i->fy   += dfy;
		i->M += (cix * dfy - ciy * dfx);  // +dfs 

		j->fx   -= dfx;
		j->fy   -= dfy;
		j->M += (-cjx * dfy + cjy * dfx); // +dfs
	}
		
	void Kin()
	{
		real tx = -ny, ty = nx;
		real cix = x - i->x;
		real ciy = y - i->y; 

		real cjx = x - j->x;
		real cjy = y - j->y;

		real vx = i->vx - j->vx - ciy * i->vrot + cjy * j->vrot;
		real vy = i->vy - j->vy + cix * i->vrot - cjx * j->vrot;

		vn = vx * nx + vy * ny;
		vt = vx * tx + vy * ty;
	}
	
	void CDcoeff(real en, real et)
	{	
		real tx = -ny, ty = nx;

		invmi   = 1.0 / i->mass;
		invmomi = 1.0 / i->mom;
		real cix = x - i->x;
		real ciy = y - i->y;
		real cin = cix * nx + ciy * ny;
		real cit = cix * tx + ciy * ty; 		

		invmj   = 1.0 / j->mass;
		invmomj = 1.0 / j->mom;
		real cjx = x - j->x;
		real cjy = y - j->y;
		real cjn = cjx * nx + cjy * ny;
		real cjt = cjx * tx + cjy * ty;

		fac_en = (1.0 + en);
		fac_et = (1.0 + et);

		Wnn = invmi + invmj + (cit*cit)*invmomi + (cjt*cjt)*invmomj;
		Wtt = invmi + invmj + (cin*cin)*invmomi + (cjn*cjn)*invmomj;
		Wnt = (cin*cjt)*invmomi + (cin*cjt)*invmomj;
		
		real invdet = 1.0 / (Wnn*Wtt - Wnt*Wnt); // it can be shown that det(W) > 0
		invWnn =  Wtt * invdet;
		invWtt =  Wnn * invdet;
		invWnt = -Wnt * invdet;
	}
	
	real An(real invdt)
	{
		return (
			- fac_en * vn*invdt + Wnn*fn + Wnt*ft
			- invmi  * (i->fx*nx + i->fy*ny)
			+ invmj  * (j->fx*nx + j->fy*ny)
			);
	}
	
	real At(real invdt)
	{
		real tx = -ny, ty = nx;
		return (
			- fac_et * vt*invdt + Wnt*fn + Wtt*ft
			- invmi  * (i->fx*tx + i->fy*ty)
			+ invmj  * (j->fx*tx + j->fy*ty)
			);
	}
};


// GLOBAL DATA

vector <contact>     contacts;
vector <grain>       grains;
vector <control>     controls; // imposed controls on the first bodies

uint nitermn,nitermx,niterconv;
real epsf;
real dt;
real en,et;
real mu;
real xgrav,ygrav;

bool vel0 = false;

// FUNCTIONS

void read_data(const char* filename)
{
	ifstream file(filename);
	string token;
	
	file >> token;
	while(file)
	{	
		if      (token == "xgrav")            file >> xgrav;
		else if (token == "ygrav")            file >> ygrav;
		else if (token == "epsf")             file >> epsf;
		else if (token == "niterconv")        file >> niterconv;
		else if (token == "nitermx")          file >> nitermx;
		else if (token == "nitermn")          file >> nitermn;
		else if (token == "mu")               file >> mu;
		else if (token == "en")               file >> en;
		else if (token == "et")               file >> et;
		else if (token == "dt")               file >> dt;
		else if (token == "vel0")             vel0 = true;
		else if (token == "END_PARAMETERS")   break;
		else cerr << "!! Parameter inconnu : " << token << endl;

		file >> token;
	}
	
	cout << "       xgrav = " << xgrav << endl;
	cout << "       ygrav = " << ygrav << endl;
	cout << "        epsf = " << epsf << endl;
	cout << "   niterconv = " << niterconv << endl;
	cout << "     nitermx = " << nitermx << endl;
	cout << "     nitermn = " << nitermn << endl;
	cout << "          mu = " << mu << endl;
	cout << "          en = " << en << endl;
	cout << "          et = " << et << endl;
	cout << "          dt = " << dt << endl;
	
	uint ngrains;
	file >> ngrains;
	cout << "Total number of bodies: " << ngrains << endl;
	grain G;
	for (uint i = 0 ; i < ngrains ; ++i) {
		file >> G.mass >> G.mom >> G.rad >> G.x >> G.y >> G.rot >> G.vx >> G.vy >> G.vrot;
		grains.push_back(G);
	}
	
	uint ncontacts;
	file >> ncontacts;
	cout << "Total number of contacts: " << ncontacts << endl;
	contact C;
	for (uint c = 0 ; c < ncontacts ; ++c) {
		file >> C.idi >> C.idj; // Start at 1
		C.idi -= 1; C.idj -= 1; // Now it starts at 0
		C.i = &grains[C.idi];
		C.j = &grains[C.idj];
		file >> C.x >> C.y >> C.nx >> C.ny;
		contacts.push_back(C);
	}
	
	uint ncontrols;
	file >> ncontrols;
	cout << "Total number of controls: " << ncontrols << endl;
	control ctr;
	for (uint c = 0 ; c < ncontrols ; ++c) {
		file >> ctr.vxImposed >> ctr.vyImposed >> ctr.vrotImposed;
		file >> ctr.xvalue >> ctr.yvalue >> ctr.rotvalue;
		controls.push_back(ctr);
	}
}


void read_data_next(const char* filename)
{
	ifstream file(filename);
	string token;
	
	// Only the parameters in the first file is accounted for 
	file >> token;
	while(file)
	{	
		if (token == "END_PARAMETERS")   break;
		file >> token;
	}
	
	uint ngrains;
	file >> ngrains;
	if (ngrains != grains.size()) cout << "Bad number of grains" << endl;
	for (uint i = 0 ; i < ngrains ; ++i) {
		file >> grains[i].mass >> grains[i].mom >> grains[i].rad 
		>> grains[i].x >> grains[i].y >> grains[i].rot 
		>> grains[i].vx >> grains[i].vy >> grains[i].vrot;
	}
	
	// save and clear current contacts
	vector <contact> contacts_old;
	for (uint c = 0 ; c < contacts.size() ; ++c) {
		contacts_old.push_back( contacts[c] );
	}
	contacts.clear();
	
	// new contacts
	uint ncontacts;
	file >> ncontacts;
	contact C;
	for (uint c = 0 ; c < ncontacts ; ++c) {
		file >> C.idi >> C.idj; // Start at 1
		C.idi -= 1; C.idj -= 1; // Now it starts at 0
		C.i = &grains[C.idi];
		C.j = &grains[C.idj];
		file >> C.x >> C.y >> C.nx >> C.ny;
		contacts.push_back(C);
	}
	
	// restore known forces
	size_t k, kold = 0;
	for (k = 0 ; k < contacts.size() ; ++k) {	
		while (kold < contacts_old.size() && contacts_old[kold].idi < contacts[k].idi) ++kold;
		if (kold == contacts_old.size()) break;

		while (kold < contacts_old.size() && contacts_old[kold].idi == contacts[k].idi && contacts_old[kold].idj < contacts[k].idj) ++kold;
		if (kold == contacts_old.size()) break;

		if (contacts_old[kold].idi == contacts[k].idi && contacts_old[kold].idj == contacts[k].idj) {
			contacts[k].fn = contacts_old[kold].fn;
			contacts[k].ft = contacts_old[kold].ft;
			
			++kold;
		}
	}
	
	// controls are ignored for not-first reading
}

void fres()
{
	cout << "Take care of forces apply on boudaries and volumic forces." << endl;
	// External forces
	// Force-controled bodies are NOT subjected to acceleration field
	for (uint i=0 ; i < controls.size() ; ++i) {	
		if (controls[i].vxImposed) grains[i].fx = 0.0;
		else grains[i].fx = controls[i].xvalue;
		
		if (controls[i].vyImposed) grains[i].fy = 0.0;
		else grains[i].fy = controls[i].yvalue;
		
		if (controls[i].vrotImposed) grains[i].M = 0.0;
		else grains[i].M = controls[i].rotvalue;
	} 

	for (uint i = controls.size() ; i < grains.size() ; ++i) {
		grains[i].fx =  xgrav * grains[i].mass;
		grains[i].fy =  ygrav * grains[i].mass;
		grains[i].M  =  0.0;
	}
}

int gs_iter()
{
	uint c;

	real fn0,an;
	real ft0,at;
	real gn,gt;
	real dfn = 0.0, dft = 0.0;
	real fz1 = 0.0, fz2 = 0.0, dfz = 0.0;
	uint nstop = 0;
	real invdt = 1.0/dt;
	contact* oxo;

	// Main CD loop
	for (uint kiter = 1 ; kiter <= nitermx ; ++kiter)
	{
		c = 0;
		while (c < contacts.size()) {
			oxo = &(contacts[c]);
			
			// BIDOUILLAGE***********
			mu = 0.3;
			if (oxo->idi<=3) mu =0.5; 
			// **********************

			// .......... Normal forces
			fn0 = oxo->fn;
			ft0 = oxo->ft;
			
			an = oxo->An(invdt);
			at = oxo->At(invdt);

			gn = oxo->invWnn * an + oxo->invWnt * at;
			
			// Signorini condition
			if (gn > 0.0) oxo->fn = gn;
			else          oxo->fn = 0.0;

			dfn = oxo->fn - fn0;

			// .......... Tangential forces
			real mufn = mu * (oxo->fn);
			gt = oxo->invWnt * an + oxo->invWtt * at;

			// Coulomb friction law
			if      (gt >=  mufn) oxo->ft =  mufn;
			else if (gt <= -mufn) oxo->ft = -mufn;
			else                  oxo->ft =  gt;

			dft =  oxo->ft - ft0;

			// .......... Normal force sum
			fz2 += fabs(oxo->fn);

			// .......... Write new values of force resultants
			oxo->Res(dfn, dft);

			for (uint i = 0 ; i < controls.size() ; ++i) {
				if (controls[i].vxImposed)   grains[i].fx = 0.0;
				if (controls[i].vyImposed)   grains[i].fy = 0.0;
				if (controls[i].vrotImposed) grains[i].M  = 0.0;
			}

			++c;
		} // while-loop on contacts

		// .......... Check convergence
		if( kiter >= nitermn ) {
			if (fz2 != 0.0) dfz = fabs(fz2 - fz1) / fz2;
			else            dfz = 0.0;

			if (dfz < epsf) {
				++nstop;
				if (nstop == niterconv) return kiter;
			} else nstop = 0;
		}

		fz1 = fz2;
		fz2 = 0.0;
	} // for-loop on iter

	return nitermx;
}

// _________________  SAVING FUNCTIONS

void savePS()
{
	//ofstream file("forces.ps");
	// Todo	 
}

void saveRAW()
{
	ofstream fg("sample");
	for (uint i=0;i<grains.size();++i) {
		fg << grains[i].x << ' ' << grains[i].y << ' ' << grains[i].rot << ' '
			<< grains[i].vx << ' ' << grains[i].vy << ' ' << grains[i].vrot << ' '
			<< grains[i].fx << ' ' << grains[i].fy << ' ' << grains[i].M << ' '
			<< grains[i].mass << ' ' << grains[i].mom << ' ' << grains[i].rad << endl; 	
	}

	ofstream fc("network");
	for (uint c=0;c<contacts.size();++c) {
		fc << contacts[c].x << ' ' << contacts[c].y << ' ' << contacts[c].fn << ' ' << contacts[c].ft << ' '
			<< contacts[c].nx << ' ' << contacts[c].ny << ' ' 
			<< contacts[c].idi + 1 << ' ' <<  contacts[c].idj + 1 << ' '
			<< contacts[c].vn << ' ' << contacts[c].vt << endl;
	}

	ofstream fctr("control");
	for (uint c=0;c<controls.size();++c) {
		fctr << controls[c].vxImposed << ' ' << controls[c].vyImposed << ' ' << controls[c].vrotImposed << ' ' 
			<< controls[c].xvalue << ' ' << controls[c].yvalue << ' ' << controls[c].rotvalue  << endl;
	}	 
}

void saveRAW(int n)
{
	char name[256];
	printf(name,"sample%d", n);
	ofstream fg(name);
	for (uint i=0;i<grains.size();++i) {
		fg << grains[i].x << ' ' << grains[i].y << ' ' << grains[i].rot << ' '
			<< grains[i].vx << ' ' << grains[i].vy << ' ' << grains[i].vrot << ' '
			<< grains[i].fx << ' ' << grains[i].fy << ' ' << grains[i].M << ' '
			<< grains[i].mass << ' ' << grains[i].mom << ' ' << grains[i].rad << endl; 	
	}

	sprintf(name,"network%d", n);
	ofstream fc(name);
	for (uint c=0;c<contacts.size();++c) {
		fc << contacts[c].x << ' ' << contacts[c].y << ' ' << contacts[c].fn << ' ' << contacts[c].ft << ' '
			<< contacts[c].nx << ' ' << contacts[c].ny << ' ' 
			<< contacts[c].idi + 1 << ' ' <<  contacts[c].idj + 1 << ' '
			<< contacts[c].vn << ' ' << contacts[c].vt << endl;
	}
}

void saveForces()
{
	ofstream file("forces.dat");
	if (!file) {
		cerr << "cannot open file forces.dat" << endl;
		return;
	}
	
	for (uint c=0;c<contacts.size();++c) {
		file << contacts[c].fn << ' ' << contacts[c].ft << ' '
			<< contacts[c].nx << ' ' << contacts[c].ny << endl;
	}
	
}

void saveDEMbox_disks(/*uint N = 0*/)
{
	ofstream hisfile("data_0.dem");
	if (!hisfile) {
		cerr << "cannot open file data_0.dem" << endl;
		return;
	}

	// Basic header with fake data

	hisfile <<"Simulation{"<< endl;
	hisfile <<"dt " << dt << endl;
	hisfile <<"t " << 0 << endl;
	hisfile <<"ihis "<< 0 << endl;
	hisfile <<"iter "<< 0 << endl;
	hisfile <<"}"<<endl;

	hisfile <<"Physic{"<<endl;
	hisfile <<"}"<<endl;

	hisfile <<"Data_table{"<<endl;
	hisfile <<"ngroup 1"<<endl;
	hisfile <<"}"<<endl;

	hisfile <<"Properties{"<<endl;
	hisfile <<"add density"<<endl;
	hisfile <<"set density 0 1000"<<endl;
	hisfile <<"}"<<endl;

	//hisfile << "System{" << endl;
	//hisfile << "period " << hxx << " " << hyy << endl; // pour LBSE prog
	//hisfile << "}" << endl << endl;

	hisfile.setf(ios_base::scientific);

	hisfile << "Sample{" << endl;
	for (uint i = 0 ; i < grains.size() ; ++i) {
		hisfile << "disk 0 " << grains[i].rad << " " 
			<< grains[i].x << " " << grains[i].y << " " << grains[i].rot << "  "
			<< grains[i].vx << " " << grains[i].vy << " " << grains[i].vrot << endl;
	}
	hisfile << "}" << endl << endl;

	hisfile << "Network{" << endl;
	for (uint c=0;c<contacts.size();++c) {

		hisfile << "dkdk " << contacts[c].idi << ' ' << contacts[c].idj << ' '
			<< contacts[c].x << ' ' << contacts[c].y << ' '
			<< contacts[c].nx << ' ' << contacts[c].ny << ' '
			<< contacts[c].fn << ' ' << contacts[c].ft << ' ' << 0.0
			<< endl;
	}      

	hisfile << "}" << endl << endl;

	hisfile << flush;
}
