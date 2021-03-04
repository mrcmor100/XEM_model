#ifndef _XEM_H
#define _XEM_H

namespace XEM 
{

	//input
	//       Z,N: proton and neutron number of the nucleus.	;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//Tb and Ta will be used for radiated XS only if they are both positive
	//Tb: material thickness in unit of radiation length before scattering;
	//Ta: material thickness in unit of radiation length after scattering;

double GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb=-0.001, double Ta=-0.001);
}
#endif

