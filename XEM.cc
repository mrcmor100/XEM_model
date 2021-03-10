
namespace XEM 
{
	extern "C" 
	{
	  //(e,ep,theta,a,z,sigdis_new,sig_qe_new)
	  void xem_model_(double* Ei, double* Ep, double* ang, double* A, double* Z, double* dis_xs, double *qe_xs);
	}

	//input
	//       Z,N: proton and neutron number of the nucleus.	;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//Tb and Ta will be used for radiated XS only if they are both positive
	//Tb: material thickness in unit of radiation length before scattering;
	//Ta: material thickness in unit of radiation length after scattering;
	double GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb, double Ta)
	{
		double NZ, NA;
		NZ = Z;
		NA = Z+N;
		theta = theta * 180. / 3.14159;  //My model takes theta in Degrees!
		double XS, dis_XS, qe_XS;
#ifdef WIN32
		//this fortran routine does not work in windows, do not know why
		return 1.0;
#else
		xem_model_(&Ei, &Ef, &theta, &NA, &NZ, &dis_XS, &qe_XS);
		XS = dis_XS+qe_XS;
#endif
		return XS / 1000.;
	}
}

