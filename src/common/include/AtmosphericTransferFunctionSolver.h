#ifndef NCPAPROP_ATMOSPHERICTRANSFERFUNCTIONSOLVER_H_INCLUDED
#define NCPAPROP_ATMOSPHERICTRANSFERFUNCTIONSOLVER_H_INCLUDED

namespace NCPA {

	class AtmosphericTransferFunctionSolver {

	public:
		virtual int solve() = 0;
		virtual ~AtmosphericTransferFunctionSolver() { }

	};

}






#endif
