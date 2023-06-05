#ifndef BINARYREADER_H
#define BINARYREADER_H

#include <iostream>
#include <cstring>

typedef long long int64;
typedef unsigned int uint;
typedef unsigned char byte; 

namespace NCPA {

/**
 * A class for reading binary data in a portable way.  The class determines the native byte
 * order of the system and does any byte shifting necessary to read data of a specified
 * endian-ness.
 */
class BinaryReader {

protected:
        static bool bigEndianNative;
        static bool endianChecked;
        float swap( float ) const;
        int swap( int ) const;
        double swap( double ) const;
        short swap( short ) const;
        int64 swap( int64 ) const;
	uint swap( uint ) const;
	void swap( byte *buffer, int nBytes, int nSamples ) const;

public:

/**
Reads 4-byte big-endian IEEE floating-point data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readBigFloatArray( std::istream *is, int nsamp, float *target ) const;
        void readBigFloatArray( std::istream is, int nsamp, float *target ) const;

/**
Reads 4-byte little-endian IEEE floating-point data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readLittleFloatArray( std::istream *is, int nsamp, float *target ) const;
        void readLittleFloatArray( std::istream is, int nsamp, float *target ) const;
	void readLittleFloatArray( std::istream *is, int nsamp, double *target ) const;
	void readLittleFloatArray( std::istream is, int nsamp, double *target ) const;

/**
Reads 4-byte IEEE floating-point data from an input stream of a specified endian-ness.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
@param isBig TRUE for big-endian sources, FALSE for little-endian sources.
*/
        void readFloatArray( std::istream *is, int nsamp, float *target, bool isBig ) const;
        void readFloatArray( std::istream is, int nsamp, float *target, bool isBig ) const;

/**
Reads 4-byte big-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readBigIntArray( std::istream *is, int nsamp, int *target ) const;
        void readBigIntArray( std::istream is, int nsamp, int *target ) const;

/**
Reads 4-byte little-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readLittleIntArray( std::istream *is, int nsamp, int *target ) const;
        void readLittleIntArray( std::istream is, int nsamp, int *target ) const;

/**
Reads 4-byte integer data from an input stream of a specified endian-ness.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
@param isBig TRUE for big-endian sources, FALSE for little-endian sources.
*/
        void readIntArray( std::istream *is, int nsamp, int *target, bool isBig ) const;
        void readIntArray( std::istream is, int nsamp, int *target, bool isBig ) const;

/**
Reads 8-byte big-endian IEEE floating-point data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readBigDoubleArray( std::istream *is, int nsamp, double *target ) const;
        void readBigDoubleArray( std::istream is, int nsamp, double *target ) const;

/**
Reads 8-byte little-endian IEEE floating-point data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readLittleDoubleArray( std::istream *is, int nsamp, double *target ) const;
        void readLittleDoubleArray( std::istream is, int nsamp, double *target ) const;

/**
Reads 8-byte IEEE floating-point data from an input stream of specified endian-ness.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
@param isBig TRUE for big-endian sources, FALSE for little-endian sources.
*/
        void readDoubleArray( std::istream *is, int nsamp, double *target, bool isBig ) const;
        void readDoubleArray( std::istream is, int nsamp, double *target, bool isBig ) const;

/**
Reads 2-byte big-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readBigShortArray( std::istream *is, int nsamp, short *target ) const;
        void readBigShortArray( std::istream is, int nsamp, short *target ) const;

/**
Reads 2-byte big-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readLittleShortArray( std::istream *is, int nsamp, short *target ) const;
        void readLittleShortArray( std::istream is, int nsamp, short *target ) const;

/**
Reads 2-byte integer data from an input stream of specified endian-ness.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
@param isBig TRUE for big-endian sources, FALSE for little-endian sources.
*/
        void readShortArray( std::istream *is, int nsamp, short *target, bool isBig ) const;
        void readShortArray( std::istream is, int nsamp, short *target, bool isBig ) const;

/**
Reads 8-byte big-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readBigLongArray( std::istream *is, int nsamp, int64 *target ) const;
        void readBigLongArray( std::istream is, int nsamp, int64 *target ) const;

/**
Reads 8-byte little-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readLittleLongArray( std::istream *is, int nsamp, int64 *target ) const;
        void readLittleLongArray( std::istream is, int nsamp, int64 *target ) const;

/**
Reads 8-byte integer data from an input stream of specified endian-ness.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
@param isBig TRUE for big-endian sources, FALSE for little-endian sources.
*/
        void readLongArray( std::istream *is, int nsamp, int64 *target, bool isBig ) const;
        void readLongArray( std::istream is, int nsamp, int64 *target, bool isBig ) const;

/**
Reads 3-byte big-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readBig3ByteIntArray( std::istream *is, int nsamp, int *target ) const;
        void readBig3ByteIntArray( std::istream is, int nsamp, int *target ) const;

/**
Reads 3-byte little-endian integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readLittle3ByteIntArray( std::istream *is, int nsamp, int *target ) const;
        void readLittle3ByteIntArray( std::istream is, int nsamp, int *target ) const;

/**
Reads 4-byte big-endian unsigned integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
	void readBigUnsignedIntArray( std::istream *is, int nsamp, uint *target ) const;
	void readBigUnsignedIntArray( std::istream is, int nsamp, uint *target ) const;

/**
Reads 4-byte little-endian unsigned integer data from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
	void readLittleUnsignedIntArray( std::istream *is, int nsamp, uint *target ) const;
	void readLittleUnsignedIntArray( std::istream is, int nsamp, uint *target ) const;

/**
Reads 4-byte unsigned integer data from an input stream of specified endian-ness.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
@param isBig TRUE if the source is big-endian, FALSE if it is little-endian.
*/
	void readUnsignedIntArray( std::istream *is, int nsamp, uint *target, bool isBig ) const;
	void readUnsignedIntArray( std::istream is, int nsamp, uint *target, bool isBig ) const;


/**
Reads binary data in a specified format from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param datatype The 2-character CSS datatype code indicating the source format.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
        void readFromCSSCode( std::istream *is, int nsamp, const char *datatype, 
				double *target ) const;
        void readFromCSSCode( std::istream is, int nsamp, const char *datatype, 
				double *target ) const;

/**
Reads binary data in a specified format from an input stream.
@param is The stream from which to read the binary data.
@param nsamp The number of values to read.
@param datatype The 2-character CSS datatype code indicating the source format.
@param target The array in which to place the values.  It is assumed that the array has
		already been allocated.
*/
	void readFromCSSCode( std::istream *is, int nsamp, const std::string datatype,
				double *target ) const;
	void readFromCSSCode( std::istream is, int nsamp, const std::string datatype,
				double *target ) const;

/**
Returns the sample size associated with a CSS datatype code.
@param datatype a 2-character CSS datatype code.
@return The number of bytes per sample associated with the supplied code.
*/
	int getSampleSize( const std::string datatype ) const;

/**
Returns the sample size associated with a CSS datatype code.
@param datatype a 2-character CSS datatype code.
@return The number of bytes per sample associated with the supplied code.
*/
	int getSampleSize( const char *datatype ) const;

/**
Default constructor.
*/
        BinaryReader();

/**
Returns the endian-ness of the native system.
@return TRUE if the system is big-endian, FALSE if it is little-endian.
*/
	bool nativeIsBigEndian() const;

};
}
#endif // BINARYREADER_H
