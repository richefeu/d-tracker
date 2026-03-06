// Programme pour la determination des forces de contacts
// a partir des positions des grains, leurs vitesses, des positions des points de contact
// et les reperes locaux.
//
// author: Vincent Richefeu
// year: 2011

#include "CDForce.hpp"


void resolve()
{
	fres();
	for (uint k = 0 ; k < contacts.size() ; ++k) {
		contacts[k].Kin();
		contacts[k].CDcoeff(en,et);
	}
	
	cout << "Begin iterations." << endl;
	int nbi = gs_iter();
	cout << "done." << endl;
	cout << "Number of iterations: " << nbi << endl;
}

int main()
{
	int nbeg = 2;
	int nend = 4;
	
	char name[256];
	
	sprintf(name,"1g2e_data_%d.txt", nbeg);
	read_data(name);
	resolve();
	saveRAW();
	saveRAW(nbeg);
	
	for (int n = nbeg+1 ; n < nend; ++n) {
		sprintf(name, "1g2e_data_%d.txt", n);
		read_data_next(name);
		resolve();
		saveRAW(n);
	}
		
	return 0;
}

/*
int main()
{
	// Read data from matlab script
	read_data("1g2e_data.txt");
	
	// Eventually set all velocities to zero
	if (vel0) {
		cout << "Velocities set to ZERO." << endl;
		for (uint i = 0 ; i < grains.size() ; ++i) {
			grains[i].vx = grains[i].vy = grains[i].vrot = 0.0;
		}
	} 
	
	
	fres();
	for (uint k = 0 ; k < contacts.size() ; ++k) {
		contacts[k].Kin();
		contacts[k].CDcoeff(en,et);
	}
	
	cout << "Begin iterations." << endl;
	int nbi = gs_iter();
	cout << "done." << endl;
	cout << "Number of iterations: " << nbi << endl;

	// Eventually update the velocities (?)
	//savePS();
	
	saveRAW();
	saveDEMbox_disks();
	saveForces();
	
	return 0;
}

*/




