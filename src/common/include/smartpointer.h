#ifndef NCPA__SMARTPOINTER_H_INCLUDED
#define NCPA__SMARTPOINTER_H_INCLUDED

#include <cstddef>

namespace NCPA {

	template<class C>
	class SmartPointer {
		private:
			C* ptr;

		public:
			typedef C element_type;

			explicit SmartPointer( C* tgt = 0 ) throw() { ptr = tgt; }

			SmartPointer(const SmartPointer<C>& ptrCopy) throw() {
				ptr = ptrCopy.get();
//				if (ptrCopy) {
//					ptr = ptrCopy;
//				} else {
//					ptr = NULL;
//				}
			}
			~SmartPointer() { if (ptr) { delete (ptr); } }

			SmartPointer<C>& operator=(SmartPointer<C> &ptrCopy) throw() {
				if (ptrCopy.ptr) {
					ptr = ptrCopy.ptr;
				} else {
					ptr = 0;
				}
				return *this;
			}

			C& operator*() const throw() { return *ptr; }
			C *operator->() const throw() { return ptr; }
			C *get() const throw() { return ptr; }
	};

}

#endif
