#ifndef NCPA__UNIT_TESTER_H_INCLUDED
#define NCPA__UNIT_TESTER_H_INCLUDED

#include <iostream>

namespace NCPA {

	class UnitTester {
	public:
		UnitTester(bool b = false) {
			verbose = b;
		}

		bool test( bool b ) {
			total++;
			b ? passed++ : failed++;
			return b;
		}

		bool test( std::string a, const char *b ) {
			return this->test( a, std::string(b) );
		}

		template<typename T>
		bool test( T a, T b ) {
			total++;
			if (verbose) {
				std::cout << "Test: " << a << " == " << b << ": ";
			}
			if (a == b) {
				passed++;
				if (verbose) {
					std::cout << "passed" << std::endl;
				}
				return true;
			} else {
				failed++;
				if (verbose) {
					std::cout << "failed" << std::endl;
				}
				return false;
			}
		}

		void print_results() {
			std::cout << passed << " out of " << total << " tests passed." << std::endl;
		}

		bool all() { return (passed == total); }

		void set_verbose( bool b ) {
			verbose = b;
		}

		size_t passed = 0;
		size_t failed = 0;
		size_t total = 0;
		bool verbose = false;

	};
}

#endif
