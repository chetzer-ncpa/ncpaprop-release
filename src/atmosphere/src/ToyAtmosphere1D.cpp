#include "ToyAtmosphere1D.h"



NCPA::ToyAtmosphere1D::ToyAtmosphere1D() : NCPA::Atmosphere1D() {
	size_t n = 901;
	double *alt = new double [ n ];
	double *zw  = new double [ n ];
	double *mw  = new double [ n ];
	double *T   = new double [ n ];
	double *rho = new double [ n ];
	double *pr  = new double [ n ];

	// Make temperature, density profiles, based on form from Lingevitch et al.,1999
	double T_o    = 288.2;
	double rho_o  = 1.225;
	double A[8] = { -3.9082017E-02, -1.1526465E-03,  3.2891937E-05, -2.0494958E-07, 
	           -4.7087295E-02,  1.2506387E-03, -1.5194498E-05,  6.518877E-08 };
	double B[8] = { -4.9244637E-03, -1.2984142E-06, -1.5701595E-06,  1.5535974E-08,
	           -2.7221769E-02,  4.247473E-04, -3.958318E-06,  1.7295795E-08 }; 
	double T_nm    = 1.0;
	double T_dnm   = 1.0;
	double rho_nm  = 0.0;
	double rho_dnm = 1.0;
	double dz      = 150.0 / ((double)(n-1));

	for (size_t i=0; i<n; i++) {
	  alt[i] = i*dz;
	  for (size_t j=0; j<8; j++) {
	      if (j<4) {
	          rho_nm  = rho_nm  + A[j]*pow((alt[i]),j+1);
	          rho_dnm = rho_dnm + B[j]*pow((alt[i]),j+1);
	      }
	      else {
	          T_nm  = T_nm  + A[j]*pow((alt[i]),j-3);
	          T_dnm = T_dnm + B[j]*pow((alt[i]),j-3);
	      }
	  }
	  T[i]   = T_o*(T_nm/T_dnm);
	  rho[i] = rho_o*pow(10,rho_nm/rho_dnm);
	  pr[i]  = rho[i]*GASCONSTANT*T[i];
	  zw[i]  = 0.0;
	  mw[i]  = 0.0;
	  T_nm    = 1.0;
	  T_dnm   = 1.0;
	  rho_nm  = 0.0;
	  rho_dnm = 1.0;
	}

	// add parameters
	z_ = NCPA::VectorWithUnits( n, alt, NCPA::Units::fromString( "km" ) );
	add_property( "Z0", 0.0, NCPA::Units::fromString( "km" ) );
	add_property( "T", n, T, NCPA::Units::fromString( "K" ) );
	add_property( "V", n, mw, NCPA::Units::fromString( "m/s" ) );
	add_property( "P", n, pr, NCPA::Units::fromString( "Pa" ) );
	add_property( "RHO", n, rho, NCPA::Units::fromString( "kg/m3" ) );
	make_gaussian_parameter_( "U", 50.0, 60000.0, 12500.0 );

	delete [] alt;
	delete [] T;
	delete [] mw;
	delete [] pr;
	delete [] rho;
	delete [] zw;
}



NCPA::ToyAtmosphere1D::~ToyAtmosphere1D() { }



void NCPA::ToyAtmosphere1D::make_gaussian_parameter_( const std::string &new_key,
	double amplitude, double height, double width ) {

	double gaussarg;
	double *newprop = new double[ nz() ];
	for (size_t i=0; i < nz(); i++) {
		double z = z_[ i ].get();
		double alt_m = NCPA::Units::convert( z, get_altitude_units(), NCPA::Units::fromString( "m" ) );
	    gaussarg   = -( alt_m - height ) * ( alt_m - height ) / ( 2 * width * width );
	    newprop[i] = amplitude * std::exp( gaussarg );
    }

    add_property( new_key, nz(), newprop, NCPA::Units::fromString( "m/s" ) );
    delete [] newprop;
}
