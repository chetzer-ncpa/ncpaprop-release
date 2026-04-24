#pragma once

namespace NCPA {

	class AtmosphericTransferFunctionSolver {

	public:
		virtual int solve() = 0;
		virtual ~AtmosphericTransferFunctionSolver() { }

	};

}

