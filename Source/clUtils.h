/*******************************************************************************
  A collection of general utility functions
  Author: Rolf Verberg
*******************************************************************************/

#ifndef CLUTILS_H__INCLUDED_
#define CLUTILS_H__INCLUDED_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FALSE
#define FALSE	0
#endif
#ifndef TRUE
#define TRUE	1
#endif
#ifndef NULL
#define NULL	0
#endif

#ifdef _DEBUG
#ifndef DEBUG
#define DEBUG	1
#endif
#endif


namespace clUtils_namespace
{
//	typedef __int32 my_int;
	typedef int my_int;
	typedef float my_float;
	typedef double my_double;

	extern int _little_endian_io;
	extern int _n_warn;
	extern int _size_my_int;
	extern int _size_my_float;
	extern int _size_my_double;

	const int MAX_WARN=10;
	const char* const _name0="0123456789";
	const double PI=acos(-1.0);
	const double TWO_PI=2.0*acos(-1.0);
	const int LINE_LENGTH=256;

	struct vector3d
	{
		double x,y,z;
	};

	struct function3d
	{
		double x,y,z,vol,*f;
	};

	struct function2d
	{
		double x,y,area,*f;
	};

	void warning(const char* const);
	void warning(const char* const,const char* const);
	void fatal_err(const char* const);
	void fatal_err(const char* const,const int);
	void fatal_err(int,const char* const);

	char* lowercase(char* const);
	char* uppercase(char* const);

	double degtorad(const double x);
	double radtodeg(const double x);

	int remove_fileextension(char*);
	int fileextension(char*,const char* const,const int);
	int filename3(char*,int);
	int filename3(char*,int,const int);
	int filename3(char*,const char* const,int,const int);
	int filename(char*,const char* const,int,const int,const int adddot=0);

	char** chararray(char**,const int,const int);
	char** delete_chararray(char**,const int);
	double** matrix_low_row(double**,const int);
	double** matrix_low_row_mpi(double**,const int,const int* const);
	double** matrix_low_col(double**,const int);
	double*** matrix_low_row(double***,const int[]);
	double*** matrix_low_row(double***,const int,const int);
	double*** matrix_low_row_mpi(double***,const int,const int* const,const int);
	void print_matrix_low_col(const char* const,const double* const* const,int);
	vector3d*** vectorfield(vector3d***,const int[]);
	vector3d*** vectorfield(vector3d***,const int,const int,const int);
	function2d** functionfield(function2d**,const int[],const int);
	function2d** functionfield(function2d**,const int,const int,const int);
	function2d** delete_functionfield(function2d**);
	function3d*** functionfield(function3d***,const int[],const int);
	function3d*** functionfield(function3d***,const int,const int,const int,const int);
	function3d*** delete_functionfield(function3d***);

	double distance(const double* const,const double* const,const int num_dim=3);
	double distance(const double* const,const double* const,double&,const int num_dim=3);
	double distance(const vector3d* const,const vector3d* const);
	double distance(const vector3d* const,const vector3d* const,double&);
	double innerproduct(const double* const,const double* const);
	void outerproduct(const double* const,const double* const,double* const);
	int get_barycentric_coords(const double* const,const double* const,const double* const,const double* const,double&,double&,double&);

	char* new_cstring(char**,const int);
	char* new_cstring(char**,const char* const);
	void delete_cstring(char*);
	int readdqstring(FILE*,char*,const int);
	int splitstring(const char* const,char,char* const,const int);
	int splitstring(const char* const,char,char* const,const int,int&);
    int remove_sur_white(const char* const, char* const);
	int readcsvstring(FILE*,char*,const int);
	int readcsvint(FILE*,int&);
	int readcsvint(FILE*,int*);
	int readcsvdouble(FILE*,double&);
	int readcsvdouble(FILE*,double*);
	int readletterstring(FILE*,char*,const int);
	int readletternumberstring(FILE*,char*,const int);
	int readuntilstring(FILE*,const char* const);
	int readuntilstring(FILE*,const char* const,const char);
	int readuntilstring_nocase(FILE*,const char* const);
	int readuntilstring_nocase(FILE*,const char* const,const char);

	int check_file(const char* const);
	int get_num_line(const char* const,const char* const);
	int get_num_line(const char* const,const bool comm_skip=FALSE);

	int read_int(FILE*);
	int read_int(FILE*,int&);
	int read_int(FILE*,int* const,const int num=1);
	double read_double(FILE*);
	int read_double(FILE*,double&);
	int read_double(FILE*,double* const,const int num=1);
	int read_char(FILE*,char*);
	int skip_read(FILE*,const int numskip=1);

	int read_binentry_int(FILE*,int&);
	int read_binentry_int(FILE*,int*,int);
	int read_binentry_int(FILE*,int&,int);
	void big_to_little_endian_int(my_int&);

	int is_power_two(int);
	int nearest_power_two(int n);
	
	/*  Return the maximum value of a and b  */
	template <class T> T max(T a,T b)
	{
		return (a>b ? a : b);
	};

	/*  Return the minimum value of a and b  */
	template <class T> T min(T a,T b)
	{
		return (a<b ? a : b);
	};

	/*  Return the maximum value of an array  */
	template <class T> T max(const T* const a,const int num)
	{
		const T *aptr=a;
		T amax=(*aptr++);
		for (int n=1; n<num; n++,aptr++) if (amax<(*aptr)) amax=(*aptr);
		return amax;
	}

	/*  Return the maximum value of an array and its index */
	template <class T> void max(const T* const a,const int num,T &amax,int &index)
	{
		amax=max(a,num,index);
	}
	template <class T> T max(const T* const a,const int num,int &index)
	{
		const T *aptr=a;
		T amax=(*aptr++);
		index=0;
		for (int n=1; n<num; n++,aptr++) 
			if (amax<(*aptr)) {
				amax=(*aptr);
				index=n;
			}
		return amax;
	}

	/*  Return the minimum value of an array  */
	template <class T> T min(const T* const a,const int num)
	{
		const T *aptr=a;
		T amin=(*aptr++);
		for (int n=1; n<num; n++,aptr++) if (amin>(*aptr)) amin=(*aptr);
		return amin;
	}

	/*  Return the minimum value of an array and its index */
	template <class T> void min(const T* const a,const int num,T &amin,int &index)
	{
		amin=min(a,num,index);
	}
	template <class T> T min(const T* const a,const int num,int &index)
	{
		const T *aptr=a;
		T amin=(*aptr++);
		index=0;
		for (int n=1; n<num; n++,aptr++) 
			if (amin>(*aptr)) {
				amin=(*aptr);
				index=n;
			}
		return amin;
	}

	/*  Return the absolute value of a  */
	template <class T> T abs(T x)
	{
		return (x<0.0 ? -x : x);
	};

	/*  Return the sign of a  */
	template <class T> T sign(T a)
	{
		return (a<0.0 ? -1.0 : 1.0);
	};

	/*  Return the absolute value of a times the sign of b  */
	template <class T> T sign(T a,T b)
	{
		return (b<0.0 ? -abs(a) : abs(a));
	};

	/*  Swap two variables */
	template <class T> void swap(T &a,T &b)
	{
		T temp;
		temp=a;
		a=b;
		b=temp;
	};

	/*  Return the square a  */
	template <class T> T sqr(T a)
	{
		return (a==0.0 ? 0.0 : a*a);
	};

	/*  Return base^power  */
	template <class T> T power(T base,int power)
	{
		if (!base && power<0) fatal_err("!base && power<0 in clUtils_namespace::power");
		T result=base;
		if (!power) return (T)(1.0);
		for (int n=1; n<abs(power); n++) result*=base;
		if (power>0) return result;
		else return (T)(1.0/result);
	};

	/*  Is within range of eps from targe  */
	template <class T> int within_eps(const T value,const T target,const double eps=1.e-6)
	{
		return (abs(value-target)<eps ? 1 : 0);
	};
	template <class T> int rel_within_eps(const T value,const T target,const double eps=1.e-6)
	{
		T error=abs(value-target);
		relerror_less_than_eps(error,target,eps);
	};

	/*  Is error smaller than epsilon? */
	template <class T> int abserror_less_than_eps(const T error,const double eps=1.e-6)
	{
		return (abs(error)<eps ? 1 : 0);
	};
	template <class T> int relerror_less_than_eps(const T error,const T norm,const double eps=1.e-6)
	{
		if (norm) return (abs(error/norm)<eps ? 1 : 0);
		else return (abs(error)<eps ? 1 : 0);
	};

	/*  Allocate a 2D matrix with subscript range m[0..idim-1][0..jdim-1]  */
	template <class T> T** matrix(T **m,const int idim,const int jdim)
	{
		int i;
		m=new T* [idim];
		if (!m) return NULL;
		m[0]=new T [idim*jdim];
		if (!m[0]) return NULL;
		for (i=1; i<idim; i++) m[i]=m[i-1]+jdim;
#ifdef ALLOCDEBUG
		fprintf(stdout,"utils_matrix %p %dx%d allocated\n",m,idim,jdim);
#endif
		return m;
	};

	/*  Allocate a 3D matrix with subscript range m[0..idim-1][0..jdim-1][0..kdim-1]  */
	template <class T> T*** matrix(T ***m,const int idim,const int jdim,const int kdim)
	{
		int i,j;
		m=new T** [idim];
		if (!m) return NULL;
		m[0]=new T* [idim*jdim];
		if (!m[0]) return NULL;
		m[0][0]=new T [idim*jdim*kdim];
		if (!m[0][0]) return NULL;
		for (j=1; j<jdim; j++) m[0][j]=m[0][j-1]+kdim;
		for (i=1; i<idim; i++) {
			m[i]=m[i-1]+jdim;
			m[i][0]=m[i-1][0]+jdim*kdim;
			for (j=1; j<jdim; j++) m[i][j]=m[i][j-1]+kdim;
		}
#ifdef ALLOCDEBUG
		fprintf(stdout,"utils_matrix3D %p %dx%dx%d allocated\n",m,idim,jdim,kdim);
#endif
		return m;
	};

	/*  Allocate a 4D matrix with subscript range m[0..dim[0]-1][0..dim[1]-1][0..dim[2]-1][0..dim[3]-1]  */
	template <class T> T**** matrix(T ****m,const int idim,const int jdim,const int kdim,const int ldim)
	{
		int i,j,k;
		m=new T*** [idim];
		if (!m) return NULL;
		m[0]=new T** [idim*jdim];
		if (!m[0]) return NULL;
		m[0][0]=new T* [idim*jdim*kdim];
		if (!m[0][0]) return NULL;
		m[0][0][0]=new T [idim*jdim*kdim*ldim];
		if (!m[0][0][0]) return NULL;
		for (k=1; k<kdim; k++) m[0][0][k]=m[0][0][k-1]+ldim;
		for (j=1; j<jdim; j++) {
			m[0][j]=m[0][j-1]+kdim;
			m[0][j][0]=m[0][j-1][0]+kdim*ldim;
			for (k=1; k<kdim; k++) m[0][j][k]=m[0][j][k-1]+ldim;
		}
		for (i=1; i<idim; i++) {
			m[i]=m[i-1]+jdim;
			m[i][0]=m[i-1][0]+jdim*kdim;
			m[i][0][0]=m[i-1][0][0]+jdim*kdim*ldim;
			for (j=1; j<jdim; j++) {
				m[i][j]=m[i][j-1]+kdim;
				m[i][j][0]=m[i][j-1][0]+kdim*ldim;
				for (k=1; k<kdim; k++) m[i][j][k]=m[i][j][k-1]+ldim;
			}
		}
#ifdef ALLOCDEBUG
		fprintf(stdout,"utils_matrix4D %p %dx%dx%dx%d allocated\n",m,idim,jdim,kdim,ldim);
#endif
		return m;
	};

	/*  Return the transpose of a symmetric matrix  */
	template <class T> void matrix_transpose(T **m,int dim[2])
	{
		if (dim[0]!=dim[1]) fatal_err("asymmetric matrix in matrix_transpose");
		int i,j;
		T swap;
		for (i=0; i<dim[0]; i++) for (j=i+1; j<dim[1]; j++) {
			swap=m[i][j];
			m[i][j]=m[j][i];
			m[j][i]=swap;
		}
	};
	template <class T> void matrix_transpose(T **m,int ddim)
	{
		int dim[2];
		dim[0]=ddim;
		dim[1]=ddim;
		matrix_transpose(m,dim);
	};

	/*  Print an array to file  */
	template <class T> void print_array(const char* const filename,T *a,int dim)
	{
		FILE *file_ptr=fopen(filename,"w");
		if (file_ptr==0) fatal_err("Could not open filename in print_array");
//		for (int i=0; i<dim; i++) fprintf(file_ptr," % .16g",(double)a[i]);
		for (int i=0; i<dim; i++) fprintf(file_ptr," % 24.16le",(double)a[i]);
		fclose(file_ptr);
	};
	/*  Print the transpose of an array to file  */
	template <class T> void print_array_transpose(const char* const filename,T *a,int dim)
	{
		FILE *file_ptr=fopen(filename,"w");
		if (file_ptr==0) fatal_err("Could not open filename in print_array");
//		for (int i=0; i<dim; i++) fprintf(file_ptr," % .16g\n",(double)a[i]);
		for (int i=0; i<dim; i++) fprintf(file_ptr," % 24.16le",(double)a[i]);
		fclose(file_ptr);
	};
	/*  Print a matrix to file  */
	template <class T> void print_matrix(const char* const filename,T **m,int dim[2])
	{
		FILE *file_ptr=fopen(filename,"w");
		if (file_ptr==0) fatal_err("Could not open filename in print_matrix");
		int i,j;
		for (i=0; i<dim[0]; i++) {
//			for (j=0; j<dim[1]; j++) fprintf(file_ptr," % .16g",(double)m[i][j]);
//			for (j=0; j<dim[1]; j++) fprintf(file_ptr," % 10.3le",(double)m[i][j]);
			for (j=0; j<dim[1]; j++) fprintf(file_ptr," % 24.16le",(double)m[i][j]);
			fprintf(file_ptr,"\n");
		}
		fclose(file_ptr);
	};
	/*  Print a matrix to file  */
	template <class T> void print_matrix(const char* const filename,T **m,int idim,int jdim)
	{
		int dim[2];
		dim[0]=idim;
		dim[1]=jdim;
		print_matrix(filename,m,dim);
	};
	/*  Print a matrix to file  */
	template <class T> void print_matrix(const char* const filename,const T* const *m,int dim[2])
	{
		FILE *file_ptr=fopen(filename,"w");
		if (file_ptr==0) fatal_err("Could not open filename in print_matrix");
		int i,j;
		for (i=0; i<dim[0]; i++) {
//			for (j=0; j<dim[1]; j++) fprintf(file_ptr," % .16g",(double)m[i][j]);
			for (j=0; j<dim[1]; j++) fprintf(file_ptr," % 24.16le",(double)m[i][j]);
			fprintf(file_ptr,"\n");
		}
		fclose(file_ptr);
	};
	/*  Print a matrix to file  */
	template <class T> void print_matrix(const char* const filename,const T* const *m,int idim,int jdim)
	{
		int dim[2];
		dim[0]=idim;
		dim[1]=jdim;
		print_matrix(filename,m,dim);
	};

	/*  Print the transpose of a matrix to file  */
	template <class T> void print_matrix_transpose(const char* const filename,T **m,int dim[2])
	{
		FILE *file_ptr=fopen(filename,"w");
		if (file_ptr==0) fatal_err("Could not open filename in print_matrix_transpose");
		int i,j;
		for (j=0; j<dim[1]; j++) {
//			for (i=0; i<dim[0]; i++) fprintf(file_ptr," % .16g",(double)m[i][j]);
			for (i=0; i<dim[0]; i++) fprintf(file_ptr," % 24.16le",(double)m[i][j]);
			fprintf(file_ptr,"\n");
		}
		fclose(file_ptr);
	};
	/*  Print the transpose of a matrix to file  */
	template <class T> void print_matrix_transpose(const char* const filename,T **m,int idim,int jdim)
	{
		int dim[2];
		dim[0]=idim;
		dim[1]=jdim;
		print_matrix_transpose(filename,m,dim);
	};

	/*  Read an entry from a binary file  */
	template <class T> int read(FILE *file_ptr,T &readvalue)
	{
		void *io_ptr=(void *)(&readvalue);
		if (fread(io_ptr,sizeof(T),1,file_ptr)!=1) return 0;
		return 1;
	};
	/*  Read num entries from a binary file  */
	template <class T> int read(FILE *file_ptr,T *readvalue,const int num)
	{
		unsigned int size=num*sizeof(T);
		void *io_ptr=(void *)(readvalue);
		if (fread(io_ptr,1,size,file_ptr)!=size) return 0;
		return 1;
	};

	/*  Check an entry from a binary file  */
	template <class T> int readcheck(FILE *file_ptr,T checkvalue)
	{
		T readvalue;
		void *io_ptr=(void *)(&readvalue);
		if (fread(io_ptr,sizeof(T),1,file_ptr)!=1) return 0;
		if (readvalue!=checkvalue) return 0;
		return 1;
	};

	/*  Write an entry to a binary file  */
	template <class T> int write(FILE *file_ptr,T writevalue)
	{
		void *io_ptr=(void *)(&writevalue);
		if (fwrite(io_ptr,sizeof(T),1,file_ptr)!=1) return 0;
		return 1;
	};
	/*  Write num entries to a binary file  */
	template <class T> int write(FILE *file_ptr,T *writevalue,const int num)
	{
		unsigned int size=num*sizeof(T);
		void *io_ptr=(void *)(writevalue);
		if (fwrite(io_ptr,1,size,file_ptr)!=size) return 0;
		return 1;
	};

	/*  Read an entry from a binary file while checking endian  */
	template <class T> int read_binentry(FILE *file_ptr,T &readvalue)
	{
		unsigned size=sizeof(readvalue);
		if (_little_endian_io) {
			void *io_ptr=(void *)(&readvalue);
			if (fread(io_ptr,size,1,file_ptr)!=1) return 0;
		} else {
			int readchar;
			unsigned char *bytes=new unsigned char [size];
			for (unsigned int i=0; i<size; i++) {
				if ((readchar=fgetc(file_ptr))==EOF) return 0;
				bytes[size-1-i]=readchar & 0xFF;
			}
			memcpy((void *)(&readvalue),(void *)bytes,size);
			delete [] bytes;
		}
		return 1;
	};
	/*  Read num entries from a binary file while checking endian  */
	template <class T> int read_binarray(FILE *file_ptr,T *readvalue,int num)
	{
		for (int n=0; n<num; n++)
			if (!read_binentry(file_ptr,readvalue[n])) return 0;
		return 1;
	};

}; // namespace clUtils_namespace

#endif // CLUTILS_H__INCLUDED_
