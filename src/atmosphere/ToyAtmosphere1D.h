#ifndef NCPAPROP_TOYATMOSPHERE1D_H_INCLUDED
#define NCPAPROP_TOYATMOSPHERE1D_H_INCLUDED

#include "Atmosphere1D.h"
#include <string>

#ifndef GASCONSTANT
#define GASCONSTANT 287.058
#endif

namespace NCPA {

	class ToyAtmosphere1D : public Atmosphere1D {

	public:
		ToyAtmosphere1D();
		~ToyAtmosphere1D();

	private:
		void make_gaussian_parameter_( const std::string &new_key, double amplitude, double center,
			double width );

	};
}



#endif