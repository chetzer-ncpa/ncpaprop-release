#include "binaryreader.h"
#include <iostream>
#include <fstream>
#include <vector>

/*
    BinaryReader.cpp: Reads binary data from files and returns arrays of
    numbers.
*/

bool NCPA::BinaryReader::bigEndianNative = false;
bool NCPA::BinaryReader::endianChecked = false;

typedef unsigned char byte;

/**
 * BinaryReader(): Default constructor
 * Checks the endian-ness of the native system and stores it for later reference.
 */
NCPA::BinaryReader::BinaryReader()
{
        if (!endianChecked) {
                // Cute test for the endian-ness of the system.  Big thanks to a program
                // by Promit Roy, who credits the source code for Quake 2.
                byte SwapTest[2] = { 1, 0 };

                if( *(short *) SwapTest == 1 )
                {
                        bigEndianNative = false;
                } else {
                        bigEndianNative = true;
                }
                endianChecked = true;
        }
        // InitEndian();
}

/**
 * swap( byte*, int, int ) const: Swap bytes within an array of raw bytes.
 */
void NCPA::BinaryReader::swap( byte *buffer, int nBytes, int nSamples ) const {
	byte temp;
	int startsamp;
	for (int j = 0; j < nSamples; j++) {
		startsamp = nBytes * j;
		for (int i = 0; i < nBytes/2; i++) {
			temp = buffer[ startsamp + nBytes - 1 - i ];
			buffer[ startsamp + nBytes - 1 - i ] = buffer[ startsamp + i ];
			buffer[ startsamp + i ] = temp;
		}
	}
}


/*
Read a 4-byte floating-point array, with endianness indicated by a boolean parameter.
*/
void NCPA::BinaryReader::readFloatArray( std::istream* infile, int samples, float* buffer, bool bigEndian ) const {
	if (bigEndian) {
		readBigFloatArray( infile, samples, buffer );
	} else {
		readLittleFloatArray( infile, samples, buffer );
	}
}

void NCPA::BinaryReader::readFloatArray( std::istream infile, int samples, float* buffer, bool bigEndian ) const {
	readFloatArray( &infile, samples, buffer, bigEndian );
}

void NCPA::BinaryReader::readBigFloatArray( std::istream* infile, int samples, float* buffer ) const {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

void NCPA::BinaryReader::readLittleFloatArray( std::istream* infile, int samples, float* buffer ) const {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

float NCPA::BinaryReader::swap( float f ) const {
        union {
                byte b[4];
                float f;
        } u1,u2;

        u1.f = f;
        u2.b[0] = u1.b[3];
        u2.b[1] = u1.b[2];
        u2.b[2] = u1.b[1];
        u2.b[3] = u1.b[0];
        return u2.f;
}



void NCPA::BinaryReader::readIntArray( std::istream* infile, int samples, int* buffer, bool bigEndian ) const {
        if (bigEndian) {
                readBigIntArray( infile, samples, buffer );
        } else {
                readLittleIntArray( infile, samples, buffer );
        }
}


void NCPA::BinaryReader::readBigIntArray( std::istream* infile, int samples, int* buffer ) const {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void NCPA::BinaryReader::readLittleIntArray( std::istream* infile, int samples, int* buffer ) const {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

int NCPA::BinaryReader::swap( int i ) const
{
        union {
                byte b[4];
                int i;
        } u1,u2;

        u1.i = i;
        u2.b[0] = u1.b[3];
        u2.b[1] = u1.b[2];
        u2.b[2] = u1.b[1];
        u2.b[3] = u1.b[0];

        return u2.i;
}

uint NCPA::BinaryReader::swap( uint i ) const
{
        union {
                byte b[4];
                uint i;
        } u1,u2;

        u1.i = i;
        u2.b[0] = u1.b[3];
        u2.b[1] = u1.b[2];
        u2.b[2] = u1.b[1];
        u2.b[3] = u1.b[0];

        return u2.i;
}



void NCPA::BinaryReader::readUnsignedIntArray( std::istream* infile, int samples, uint* buffer, bool bigEndian ) const {
	if (bigEndian) {
		readBigUnsignedIntArray( infile, samples, buffer );
	} else {
		readLittleUnsignedIntArray( infile, samples, buffer );
	}
}

void NCPA::BinaryReader::readBigUnsignedIntArray( std::istream *infile, int samples, uint *buffer ) const {
	int nBytes = samples * 4;
	infile->read((char*)buffer, nBytes);
	if (!bigEndianNative) {
		for (int i = 0; i < samples; i++) {
			buffer[i] = swap( buffer[i] );
		}
	}
}

void NCPA::BinaryReader::readLittleUnsignedIntArray( std::istream *infile, int samples, uint *buffer ) const {
        int nBytes = samples * 4;
        infile->read((char*)buffer, nBytes);
        if (bigEndianNative) {
                for (int i = 0; i < samples; i++) {
                        buffer[i] = swap( buffer[i] );
                }
        }
}

void NCPA::BinaryReader::readDoubleArray( std::istream* infile, int samples, double* buffer, bool bigEndian ) const {
        if (bigEndian) {
                readBigDoubleArray( infile, samples, buffer );
        } else {
                readLittleDoubleArray( infile, samples, buffer );
        }
}


 void NCPA::BinaryReader::readBigDoubleArray( std::istream* infile, int samples, double* buffer ) const {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void NCPA::BinaryReader::readLittleDoubleArray( std::istream* infile, int samples, double* buffer ) const {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

double NCPA::BinaryReader::swap( double d ) const {
        union {
                byte b[8];
                double d;
        } u1,u2;
        u1.d = d;

        for (int k = 0; k < 8; k++) {
            u2.b[k] = u1.b[8-1-k];
        }

        return u2.d;
}

  void NCPA::BinaryReader::readShortArray( std::istream* infile, int samples, short* buffer, bool bigEndian ) const {
        if (bigEndian) {
                readBigShortArray( infile, samples, buffer );
        } else {
                readLittleShortArray( infile, samples, buffer );
        }
}


 void NCPA::BinaryReader::readBigShortArray( std::istream* infile, int samples, short* buffer ) const {
        int nBytes = samples * 2;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void NCPA::BinaryReader::readLittleShortArray( std::istream* infile, int samples, short* buffer ) const {
        int nBytes = samples * 2;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

short NCPA::BinaryReader::swap( short s ) const
{
        union {
                byte b[2];
                short s;
        } u1,u2;
        u1.s = s;

        for (int k = 0; k < 2; k++) {
            u2.b[k] = u1.b[1-k];
        }

        return u2.s;
}

void NCPA::BinaryReader::readLongArray( std::istream* infile, int samples, int64* buffer, bool bigEndian ) const {
        if (bigEndian) {
                readBigLongArray( infile, samples, buffer );
        } else {
                readLittleLongArray( infile, samples, buffer );
        }
}


 void NCPA::BinaryReader::readBigLongArray( std::istream* infile, int samples, int64* buffer ) const {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void NCPA::BinaryReader::readLittleLongArray( std::istream* infile, int samples, int64* buffer ) const {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

int64 NCPA::BinaryReader::swap( int64 s ) const {
        union {
                byte b[8];
                int64 s;
        } u1,u2;
        u1.s = s;

        for (int k = 0; k < 8; k++) {
            u2.b[k] = u1.b[7-k];
        }

        return u2.s;
}

void NCPA::BinaryReader::readFromCSSCode( std::istream *in, int samples, const std::string datatype,
					double *data ) const {
	this->readFromCSSCode( in, samples, datatype.c_str(), data );
}

void NCPA::BinaryReader::readFromCSSCode( std::istream *in, int samples, 
					const char *datatype, double *data ) const {
	
	int *idata;
	short *sdata;
	float *fdata;

	if (std::strncmp(datatype,"s4",2) == 0 || std::strncmp(datatype,"S4",2) == 0) {
		idata = new int[ samples ];
		this->readBigIntArray( in, samples, idata );
		for (int i = 0; i < samples; i++) {
			data[ i ] = (double)idata[ i ];
		}
		delete [] idata;
	} else if (std::strncmp(datatype,"i4",2) == 0 || std::strncmp(datatype,"I4",2) == 0) {
		idata = new int[ samples ];
                this->readLittleIntArray( in, samples, idata );
                for (int i = 0; i < samples; i++) {
                        data[ i ] = (double)idata[ i ];
                }
                delete [] idata;
	} else if (std::strncmp(datatype,"s2",2) == 0 || std::strncmp(datatype,"S2",2) == 0) {
                sdata = new short[ samples ];
                this->readBigShortArray( in, samples, sdata );
                for (int i = 0; i < samples; i++) {
                        data[ i ] = (double)sdata[ i ];
                }
                delete [] sdata;
        } else if (std::strncmp(datatype,"i2",2) == 0 || std::strncmp(datatype,"I2",2) == 0) {
                sdata = new short[ samples ];
                this->readLittleShortArray( in, samples, sdata );
                for (int i = 0; i < samples; i++) {
                        data[ i ] = (double)sdata[ i ];
                }
                delete [] sdata;
        } else if (std::strncmp(datatype,"s3",2) == 0 || std::strncmp(datatype,"S3",2) == 0) {
                idata = new int[ samples ];
                this->readBig3ByteIntArray( in, samples, idata );
                for (int i = 0; i < samples; i++) {
                        data[ i ] = (double)idata[ i ];
                }
                delete [] idata;
        } else if (std::strncmp(datatype,"i3",2) == 0 || std::strncmp(datatype,"I3",2) == 0) {
                idata = new int[ samples ];
                this->readLittle3ByteIntArray( in, samples, idata );
                for (int i = 0; i < samples; i++) {
                        data[ i ] = (double)idata[ i ];
                }
                delete [] idata;
        } else if (std::strncmp(datatype,"t4",2) == 0 || std::strncmp(datatype,"T4",2) == 0) {
                fdata = new float[ samples ];
                this->readBigFloatArray( in, samples, fdata );
                for (int i = 0; i < samples; i++) {
                        data[ i ] = (double)fdata[ i ];
                }
                delete [] fdata;
        } else if (std::strncmp(datatype,"f4",2) == 0 || std::strncmp(datatype,"F4",2) == 0) {
                fdata = new float[ samples ];
                this->readLittleFloatArray( in, samples, fdata );
                for (int i = 0; i < samples; i++) {
                        data[ i ] = (double)fdata[ i ];
                }
                delete [] fdata;
        } else if (std::strncmp(datatype,"t8",2) == 0 || std::strncmp(datatype,"T8",2) == 0) {
		this->readBigDoubleArray( in, samples, data );
        } else if (std::strncmp(datatype,"f8",2) == 0 || std::strncmp(datatype,"F8",2) == 0) {
		this->readLittleDoubleArray( in, samples, data );
	} else {
		std::cerr << "Unrecognized CSS datatype!" << std::endl;
	}
}

void NCPA::BinaryReader::readBig3ByteIntArray( std::istream *in, int samples, int *data ) const {
	
        int nBytes = samples * 3;
	int temp = 0;
	byte *buffer = new byte[ nBytes ];
        in->read((char*)buffer,nBytes);
	for (int i = 0; i < samples; i++) {
		temp = 0;
		temp |= ((unsigned int)buffer[3*i]) << 24 
			| ((unsigned int)buffer[3*i+1]) << 16
			| ((unsigned int)buffer[3*i+2]) << 8;
		data[ i ] = temp >> 8;		
	}
	delete [] buffer;
}

void NCPA::BinaryReader::readLittle3ByteIntArray( std::istream *in, int samples, int *data ) const  {
	int nBytes = samples * 3;
	int temp = 0;
	byte *buffer = new byte[ nBytes ];
	in->read((char*)buffer,nBytes);
	for (int i = 0; i < samples; i++) {
		temp = 0;
		temp |= ((unsigned int)buffer[3*i+2]) << 24
			| ((unsigned int)buffer[3*i+1]) << 16
			| ((unsigned int)buffer[3*i]) << 8;
		data[ i ] = temp >> 8;
	}
	delete [] buffer;
}

int NCPA::BinaryReader::getSampleSize( const std::string datatype ) const {
	return this->getSampleSize( datatype.c_str() );
}

int NCPA::BinaryReader::getSampleSize( const char *datatype ) const {
	
	if (std::strncmp(datatype,"s4",2) == 0 || std::strncmp(datatype,"S4",2) == 0) {
		return 4;
	} else if (std::strncmp(datatype,"i4",2) == 0 || std::strncmp(datatype,"I4",2) == 0) {
		return 4;
	} else if (std::strncmp(datatype,"s2",2) == 0 || std::strncmp(datatype,"S2",2) == 0) {
		return 2;
        } else if (std::strncmp(datatype,"i2",2) == 0 || std::strncmp(datatype,"I2",2) == 0) {
		return 2;
        } else if (std::strncmp(datatype,"s3",2) == 0 || std::strncmp(datatype,"S3",2) == 0) {
		return 3;
        } else if (std::strncmp(datatype,"i3",2) == 0 || std::strncmp(datatype,"I3",2) == 0) {
		return 3;
        } else if (std::strncmp(datatype,"t4",2) == 0 || std::strncmp(datatype,"T4",2) == 0) {
		return 4;
        } else if (std::strncmp(datatype,"f4",2) == 0 || std::strncmp(datatype,"F4",2) == 0) {
		return 4;
        } else if (std::strncmp(datatype,"t8",2) == 0 || std::strncmp(datatype,"T8",2) == 0) {
		return 8;
        } else if (std::strncmp(datatype,"f8",2) == 0 || std::strncmp(datatype,"F8",2) == 0) {
		return 8;
	} else {
		return 0;
	}
	
}

bool NCPA::BinaryReader::nativeIsBigEndian() const {
	return bigEndianNative;
}

void NCPA::BinaryReader::readLittleFloatArray( std::istream *in, int nsamp, double *target ) const {

	float *temp = new float[ nsamp ];
	readLittleFloatArray( in, nsamp, temp );
	for (int i = 0; i < nsamp; i++) {
		target[ i ] = (double)temp[ i ];
	}
	delete [] temp;
}



void NCPA::BinaryReader::readBigFloatArray( std::istream infile, int samples, float* buffer ) const {
        readBigFloatArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readLittleFloatArray( std::istream infile, int samples, float* buffer ) const {
        readLittleFloatArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readIntArray( std::istream infile, int samples, int* buffer, bool bigEndian ) const {
        readIntArray( &infile, samples, buffer, bigEndian );
}


void NCPA::BinaryReader::readBigIntArray( std::istream infile, int samples, int* buffer ) const {
        readBigIntArray( &infile, samples, buffer );
}

 void NCPA::BinaryReader::readLittleIntArray( std::istream infile, int samples, int* buffer ) const {
        readLittleIntArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readUnsignedIntArray( std::istream infile, int samples, uint* buffer, bool bigEndian ) const {
	readUnsignedIntArray( &infile, samples, buffer, bigEndian );
}

void NCPA::BinaryReader::readBigUnsignedIntArray( std::istream infile, int samples, uint *buffer ) const {
	readBigUnsignedIntArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readLittleUnsignedIntArray( std::istream infile, int samples, uint *buffer ) const {
        readLittleUnsignedIntArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readDoubleArray( std::istream infile, int samples, double* buffer, bool bigEndian ) const {
        readDoubleArray( &infile, samples, buffer, bigEndian );
}


void NCPA::BinaryReader::readBigDoubleArray( std::istream infile, int samples, double* buffer ) const {
        readBigDoubleArray( &infile, samples, buffer );
}

 void NCPA::BinaryReader::readLittleDoubleArray( std::istream infile, int samples, double* buffer ) const {
        readLittleDoubleArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readShortArray( std::istream infile, int samples, short* buffer, bool bigEndian ) const {
        readShortArray( &infile, samples, buffer, bigEndian );
}


 void NCPA::BinaryReader::readBigShortArray( std::istream infile, int samples, short* buffer ) const {
        readBigShortArray( &infile, samples, buffer );
}

 void NCPA::BinaryReader::readLittleShortArray( std::istream infile, int samples, short* buffer ) const {
        readLittleShortArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readLongArray( std::istream infile, int samples, int64* buffer, bool bigEndian ) const {
        readLongArray( &infile, samples, buffer, bigEndian );
}


 void NCPA::BinaryReader::readBigLongArray( std::istream infile, int samples, int64* buffer ) const {
        readBigLongArray( &infile, samples, buffer );
}

 void NCPA::BinaryReader::readLittleLongArray( std::istream infile, int samples, int64* buffer ) const {
        readLittleLongArray( &infile, samples, buffer );
}

void NCPA::BinaryReader::readFromCSSCode( std::istream in, int samples, const std::string datatype,
					double *data ) const {
	this->readFromCSSCode( &in, samples, datatype.c_str(), data );
}

void NCPA::BinaryReader::readFromCSSCode( std::istream in, int samples, 
					const char *datatype, double *data ) const {
	
	readFromCSSCode( &in, samples, datatype, data );
}

void NCPA::BinaryReader::readBig3ByteIntArray( std::istream in, int samples, int *data ) const {
	
        readBig3ByteIntArray( &in, samples, data );
}

void NCPA::BinaryReader::readLittle3ByteIntArray( std::istream in, int samples, int *data ) const  {
	readLittle3ByteIntArray( &in, samples, data );
}

void NCPA::BinaryReader::readLittleFloatArray( std::istream in, int nsamp, double *target ) const {
	readLittleFloatArray( &in, nsamp, target );
}




