/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/

#ifndef CLCOMPLEX_H__INCLUDED_
#define CLCOMPLEX_H__INCLUDED_

#include <stdio.h>
#include <stdlib.h>
#include "clVector.h"

#ifdef DIST_MPI
#include "clMpi.h"
#endif

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
#ifndef ALLOCDEBUG
#define ALLOCDEBUG	1
#endif
#endif

namespace clComplex_namespace
{

	class clComplex
	{
		friend class clVector_complex;
		friend class clMatrix_complex;

		private:

			double _real;	// Real part
			double _imag;	// Imaginary part

		public:

			/*  Constructors  */
			clComplex(void) :
				_real(0.0),
				_imag(0.0)
			{ };
			explicit clComplex(double real,double imag=0.0) :
				_real(real),
				_imag(imag)
			{ };

			/*  Destructor  */
			virtual ~clComplex(void) { };

			/*  Default Copy Constructor  */
			clComplex(const clComplex &arg) :
				_real(0.0),
				_imag(0.0)
			{
				_real=arg._real;
				_imag=arg._imag;
			};

			/*  Default Assignment Constructor  */
			clComplex& operator= (const clComplex &arg)
			{
				if (this==&arg) return *this;

				_real=arg._real;
				_imag=arg._imag;

				return *this;
			};

			/*  Operator overloading functions  */

			clComplex& operator()(const clComplex &arg)
			{
				_real=arg._real;
				_imag=arg._imag;
				return *this;
			};
			clComplex& operator()(const double real,const double imag=0.0)
			{
				_real=real;
				_imag=imag;
				return *this;
			};

			clComplex operator+ (const double arg) const
			{
				clComplex sum;
				sum._real=_real+arg;
				sum._imag=_imag;
				return sum;
			};
			clComplex operator+ (const clComplex &arg) const
			{
				clComplex sum;
				sum._real=_real+arg._real;
				sum._imag=_imag+arg._imag;
				return sum;
			};
			friend clComplex operator+ (const double arg1,const clComplex& arg2)
			{
				clComplex sum;
				sum._real=arg1+arg2._real;
				sum._imag=arg2._imag;
				return sum;
			};

			clComplex operator- (const double arg) const
			{
				clComplex sum;
				sum._real=_real-arg;
				sum._imag=_imag;
				return sum;
			};
			clComplex operator- (const clComplex &arg) const
			{
				clComplex sum;
				sum._real=_real-arg._real;
				sum._imag=_imag-arg._imag;
				return sum;
			};
			friend clComplex operator- (const double arg1,const clComplex& arg2)
			{
				clComplex sum;
				sum._real=arg1-arg2._real;
				sum._imag=-arg2._imag;
				return sum;
			};

			clComplex operator* (const double arg) const
			{
				clComplex sum;
				sum._real=arg*_real;
				sum._imag=arg*_imag;
				return sum;
			};
			clComplex operator* (const clComplex &arg) const
			{
				clComplex sum;
				sum._real=_real*arg._real-_imag*arg._imag;
				sum._imag=_real*arg._imag+_imag*arg._real;
				return sum;
			};
			friend clComplex operator* (const double arg1,const clComplex& arg2)
			{
				clComplex sum;
				sum._real=arg1*arg2._real;
				sum._imag=arg1*arg2._imag;
				return sum;
			};

			clComplex operator/ (const double &arg) const
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clComplex::operator/\n");
					exit(1);
				}
				clComplex sum;
				sum._real=_real/arg;
				sum._imag=_imag/arg;
				return sum;
			};
			clComplex operator/ (const clComplex &arg) const
			{
				double den=arg._real*arg._real+arg._imag*arg._imag;
				if (!den) {
					fprintf(stdout,"FATAL ERROR: den=0 in clComplex::operator/\n");
					exit(1);
				}
				clComplex sum;
				sum._real=(_real*arg._real+_imag*arg._imag)/den;
				sum._imag=(_imag*arg._real-_real*arg._imag)/den;
				return sum;
			};
			friend clComplex operator/ (const double arg1,const clComplex& arg2)
			{
				double den=arg2._real*arg2._real+arg2._imag*arg2._imag;
				if (!den) {
					fprintf(stdout,"FATAL ERROR: den=0 in clComplex::operator/\n");
					exit(1);
				}
				clComplex sum;
				sum._real=(arg1*arg2._real)/den;
				sum._imag=-(arg1*arg2._imag)/den;
				return sum;
			};

			clComplex& operator+= (const double arg)
			{
				_real+=arg;
				return *this;
			};
			clComplex& operator+= (const clComplex &arg)
			{
				_real+=arg._real;
				_imag+=arg._imag;
				return *this;
			};

			clComplex& operator-= (const double arg)
			{
				_real-=arg;
				return *this;
			};
			clComplex& operator-= (const clComplex &arg)
			{
				_real-=arg._real;
				_imag-=arg._imag;
				return *this;
			};

			clComplex& operator*= (const double arg)
			{
				_real*=arg;
				_imag*=arg;
				return *this;
			};
			clComplex& operator*= (const clComplex &arg)
			{
				double real=_real*arg._real-_imag*arg._imag;
				_imag=_real*arg._imag+_imag*arg._real;
				_real=real;
				return *this;
			};

			clComplex& operator/= (const double arg)
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clComplex::operator/=\n");
					exit(1);
				}
				_real/=arg;
				_imag/=arg;
				return *this;
			};
			clComplex& operator/= (const clComplex &arg)
			{
				double den=arg._real*arg._real+arg._imag*arg._imag;
				if (!den) {
					fprintf(stdout,"FATAL ERROR: den=0 in clComplex::operator/=\n");
					exit(1);
				}
				double real=(_real*arg._real+_imag*arg._imag)/den;
				_imag=(_imag*arg._real-_real*arg._imag)/den;
				_real=real;
				return *this;
			};

			clComplex& operator++ (void)
			{
				_real++;
				return *this;
			};

			clComplex& operator-- (void)
			{
				_real--;
				return *this;
			};

			double real(void) const {return _real;};
			void set_real(const double real) {_real=real;};
			double imag(void) const {return _imag;};
			void set_imag(const double imag) {_imag=imag;};
//			double norm(void) const {return (_real*_real+_imag*_imag);};

			void write_ascii(const char *filename) const
			{
				FILE *fileptr=fopen(filename,"w");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clVector_complex::write_ascii\n");
					exit(1);
				}
				write_ascii(fileptr);
				fclose(fileptr);
			};
			void write_ascii(FILE *fileptr) const
			{
				fprintf(fileptr," %12.4le %12.4le\n",_real,_imag);
			};
	};

	/*  clVector_complex class  */
	class clVector_complex
	{
		private:

			int _dim;			// Number of complex elements in the vector
			int _dimvec;		// Actual storage space of _vector (2*_dim)
			double *_vector;	// Vector
#ifdef ALLOCDEBUG
			static int _num_vector_complex;
#endif

		public:

			/*  Constructors  */
			clVector_complex(void) :
				_dim(0),
				_dimvec(0),
				_vector(NULL)
			{ };
			explicit clVector_complex(const int dim) :
				_dim(0),
				_dimvec(0),
				_vector(NULL)
			{
				allocate(dim);
			};
			explicit clVector_complex(const int dim,const double real,const double imag=0.0) :
				_dim(0),
				_dimvec(0),
				_vector(NULL)
			{
				allocate(dim,real,imag);
			};
			explicit clVector_complex(const int dim,const clComplex &arg) :
				_dim(0),
				_dimvec(0),
				_vector(NULL)
			{
				allocate(dim,arg);
			};

			/*  Destructor  */
			virtual ~clVector_complex(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clVector_complex(const clVector_complex &arg) :
				_dim(0),
				_dimvec(0),
				_vector(NULL)
			{
				allocate(arg._dim);
				double *vectorptr=_vector;
				const double *argvectorptr=arg._vector;
				for (int n=0; n<_dimvec; n++) (*vectorptr++)=(*argvectorptr++);
			};

			/*  Default Assignment Constructor  */
			clVector_complex& operator= (const clVector_complex &arg)
			{
				if (this==&arg) return *this;

				if (!_vector) {
					if (_dim) {
						fprintf(stdout,"FATAL ERROR: illegal _dim in clVector_complex::operator=\n");
						exit(1);
					} else if (!arg._dim) {
#ifdef DEBUG
						fprintf(stdout,"WARNING: empty argument vector in clVector_complex::operator=\n");
#endif
						return *this;
					}
					allocate(arg._dim);
				}
#ifdef DEBUG
				test_indices(arg);
#endif
				const int dimvec=(_dimvec <= arg._dimvec ? _dimvec : arg._dimvec);
				double *vectorptr=_vector;
				const double *argvectorptr=arg._vector;
				for (int n=0; n<dimvec; n++) (*vectorptr++)=(*argvectorptr++);
				return *this;
			};

			/*  Operator overloading functions  */

			void operator()(const int n,const clComplex &arg)
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector_complex::operator()\n");
					exit(1);
				}
				test_index(n);
#endif
				double *vectorptr=_vector+n+n;
				(*vectorptr++)=arg._real;
				(*vectorptr)=arg._imag;
			};
			void operator()(const int n,const double real,const double imag=0.0)
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector_complex::operator()\n");
					exit(1);
				}
				test_index(n);
#endif
				double *vectorptr=_vector+n+n;
				(*vectorptr++)=real;
				(*vectorptr)=imag;
			};

			clVector_complex& operator= (const double arg)
			{
				if (!_vector) {
					fprintf(stdout,"WARNING: _vector not allocated in clVector_complex::operator=\n");
					return *this;
				}
				double *vectorptr=_vector;
				for (int n=0; n<_dim; n++,vectorptr+=2) (*vectorptr)=arg;
				return *this;
			};
			clVector_complex& operator= (const clComplex &arg)
			{
				if (!_vector) {
					fprintf(stdout,"WARNING: _vector not allocated in clVector_complex::operator=\n");
					return *this;
				}
				const double argreal=arg._real;
				const double argimag=arg._imag;
				double *vectorptr=_vector;
				for (int n=0; n<_dim; n++) {
					(*vectorptr++)=argreal;
					(*vectorptr++)=argimag;
				}
				return *this;
			};
			/*  Converts a regular array to a vector, does not allow bounce checking  */
			clVector_complex& operator= (const double *arg)
			{
				if (!arg) {
					fprintf(stdout,"WARNING: arg not allocated in clVector_complex::operator=\n");
					return *this;
				}
				if (!_vector) {
					fprintf(stdout,"WARNING: _vector not allocated in clVector_complex::operator=\n");
					return *this;
				}
				double *vectorptr=_vector;
				const double *argptr=arg;
				for (int n=0; n<_dim; n++,vectorptr+=2) (*vectorptr)=(*argptr++);
				return *this;
			};

			clVector_complex operator+ (const clVector_complex &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				const int dim=(_dim <= arg._dim ? _dim : arg._dim);
				clVector_complex sum(dim);
				double *sumptr=sum._vector;
				const double *vectorptr=_vector;
				const double *argvectorptr=arg._vector;
				for (int n=0; n<2*dim; n++) (*sumptr++)=(*vectorptr++)+(*argvectorptr++);
				return sum;
			};
			clVector_complex operator+ (const double arg) const
			{
				clVector_complex sum(_dim);
				double *sumptr=sum._vector;
				const double *vectorptr=_vector;
				for (int n=0; n<_dim; n++,sumptr+=2,vectorptr+=2) (*sumptr)=(*vectorptr)+arg;
				return sum;
			};
			clVector_complex operator+ (const clComplex &arg) const
			{
				clVector_complex sum(_dim);
				const double argreal=arg._real;
				const double argimag=arg._imag;
				double *sumptr=sum._vector;
				const double *vectorptr=_vector;
				for (int n=0; n<_dim; n++) {
					(*sumptr++)=(*vectorptr++)+argreal;
					(*sumptr++)=(*vectorptr++)+argimag;
				}
				return sum;
			};
			friend clVector_complex operator+ (const double arg1,const clVector_complex& arg2)
			{
				clVector_complex sum(arg2._dim);
				double *sumptr=sum._vector;
				const double *arg2ptr=arg2._vector;
				for (int n=0; n<arg2._dim; n++,sumptr+=2,arg2ptr+=2) (*sumptr)=arg1+(*arg2ptr);
				return sum;
			};
			friend clVector_complex operator+ (const clComplex &arg1,const clVector_complex& arg2)
			{
				clVector_complex sum(arg2._dim);
				const double argreal=arg1.real();
				const double argimag=arg1.imag();
				double *sumptr=sum._vector;
				const double *arg2ptr=arg2._vector;
				for (int n=0; n<arg2._dim; n++) {
					(*sumptr++)=argreal+(*arg2ptr++);
					(*sumptr++)=argimag+(*arg2ptr++);
				}
				return sum;
			};

			clVector_complex& operator+= (const clVector_complex &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				double *vectorptr=_vector;
				const double *argvectorptr=arg._vector;
				for (int n=0; n<_dimvec; n++) (*vectorptr++)+=(*argvectorptr++);
				return *this;
			};
			clVector_complex& operator+= (const double arg)
			{
				double *vectorptr=_vector;
				for (int n=0; n<_dim; n++,vectorptr+=2) (*vectorptr)+=arg;
				return *this;
			};
			clVector_complex& operator+= (const clComplex &arg)
			{
				const double argreal=arg._real;
				const double argimag=arg._imag;
				double *vectorptr=_vector;
				for (int n=0; n<_dim; n++) {
					(*vectorptr++)+=argreal;
					(*vectorptr++)+=argimag;
				}
				return *this;
			};

			clVector_complex& operator++ (void)
			{
				double *vectorptr=_vector;
				for (int n=0; n<_dim; n++,vectorptr+=2) (*vectorptr)++;
				return *this;
			};

			/*  Remaining member functions  */

			int get_dim(void) const {return _dim;};
			double *get_startptr(void) {return _vector;};
			const double *get_startptr(void) const {return _vector;};

			double real(const int n) const
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector_complex::real\n");
					exit(1);
				}
				test_index(n);
#endif
				return *(_vector+n+n);
			};
			void set_real(const double real) 
			{
				if (!_vector) {
					fprintf(stdout,"WARNING: _vector not allocated in clVector_complex::set_real\n");
					return;
				}
				double *vectorptr=_vector;
				for (int n=0; n<_dim; n++,vectorptr+=2) (*vectorptr)=real;
			};
			void set_real(const int n,const double real) 
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector_complex::set_real\n");
					exit(1);
				}
				test_index(n);
#endif
				*(_vector+n+n)=real;
			};
			double imag(const int n) const
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector_complex::imag\n");
					exit(1);
				}
				test_index(n);
#endif
				return *(_vector+n+n+1);
			};
			void set_imag(const double imag) 
			{
				if (!_vector) {
					fprintf(stdout,"WARNING: _vector not allocated in clVector_complex::set_imag\n");
					return;
				}
				double *vectorptr=_vector+1;
				for (int n=0; n<_dim; n++,vectorptr+=2) (*vectorptr)=imag;
			};
			void set_imag(const int n,const double imag) 
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector_complex::set_imag\n");
					exit(1);
				}
				test_index(n);
#endif
				*(_vector+n+n+1)=imag;
			};

			/*  Converts a vector to a regular array, does not allow bounce checking  */
			void copy(double *arg) const
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg not allocated in clVector_complex::copy\n");
					exit(1);
				}
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector_complex::copy\n");
					exit(1);
				}
				const double *vectorptr=_vector;
				double *argptr=arg;
				for (int n=0; n<_dimvec; n++) (*argptr++)=(*vectorptr++);
			};

			/*  Allocate vector of dimension dim  */
			double* allocate(const int dim)
			{
				if (_dim!=dim) {
					if (_vector) deallocate();
					_dim=dim;
					_dimvec=dim+dim;
					if (!_dim) return _vector;
					_vector=new double [_dimvec];
					if (!_vector) {
						fprintf(stdout,"FATAL ERROR: Could not allocate _vector in clVector_complex::allocate\n");
						exit(1);
					}
#ifdef ALLOCDEBUG
					_num_vector_complex++;
					fprintf(stdout,"vector %p allocated %d size %d\n",_vector,_num_vector_complex,_dimvec);
#endif
				}
				return _vector;
			};
			/*  Allocate vector of dimension dim and initialize real part to initreal  */
			double* allocate(const int dim,const double initreal)
			{
				allocate(dim);
				if (_vector) {
					double *vectorptr=_vector;
					for (int n=0; n<_dim; n++,vectorptr+=2) (*vectorptr)=initreal;
				}
				return _vector;
			};
			/*  Allocate vector of dimension dim and initialize to (initreal,initimag)  */
			double* allocate(const int dim,const double initreal,const double initimag)
			{
				allocate(dim);
				if (_vector) {
					double *vectorptr=_vector;
					for (int n=0; n<_dim; n++) {
						(*vectorptr++)=initreal;
						(*vectorptr++)=initimag;
					}
				}
				return _vector;
			};
			/*  Allocate vector of dimension dim and initialize to (initreal,initimag)  */
			double* allocate(const int dim,const clComplex &initval)
			{
				allocate(dim);
				if (_vector) {
					double initreal=initval._real;
					double initimag=initval._imag;
					double *vectorptr=_vector;
					for (int n=0; n<_dim; n++) {
						(*vectorptr++)=initreal;
						(*vectorptr++)=initimag;
					}
				}
				return _vector;
			};

			/*  Resize a vector to a new dimension dim  */
			/*  Data that does not fit is lost  */
			/*  Zero out cells that are not used yet */
			double* redimension(const int new_dim)
			{
				if (!new_dim) {
					deallocate();
					return _vector;
				}
				if (!_vector) {
					allocate(new_dim);
					return _vector;
				}
				if (!_dim) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector_complex::redimension\n");
					exit(1);
				}
				const int new_dimvec=new_dim+new_dim;
				double *newvector=new double [new_dimvec];
				if (!newvector) {
					fprintf(stdout,"FATAL ERROR: Could not allocate newvector in clVector_complex::redimension\n");
					exit(1);
				}
				const int dimvec=(_dimvec <= new_dimvec ? _dimvec : new_dimvec);
				double *newvectorptr=newvector;
				const double *vectorptr=_vector;
				int n=0;
				for (; n<dimvec; n++) (*newvectorptr++)=(*vectorptr++);
				for (; n<new_dimvec; n++) (*newvectorptr++)=(double)0.0;
				double *swap=_vector;
				_vector=newvector;
				_dim=new_dim;
				_dimvec=new_dimvec;
				delete [] swap;
				return _vector;
			};

			/*  Deallocate vector  */
			void deallocate(void)
			{
				if (_vector) {
#ifdef ALLOCDEBUG
					_num_vector_complex--;
					fprintf(stdout,"vector %p deallocated %d size %d\n",_vector,_num_vector_complex,_dim);
#endif
					delete [] _vector;
					_vector=NULL;
				}
				_dim=0;
				_dimvec=0;
			};

			void psd(clVector_namespace::clVector<double> &psd_vector) const
			{
				if (!_dim) return;
				if (_dim!=psd_vector.get_dim()) {
					fprintf(stdout,"FATAL ERROR: Inconsistent dimensions in clVector_complex::psd\n");
					exit(1);
				}
				double real,imag;
				const double norm=1.0/(_dim*_dim);
				const double *vectorptr=_vector;
				double *psd_vectorptr=psd_vector.get_startptr();
				for (int n=0; n<_dim; n++) {
					real=(*vectorptr++);
					imag=(*vectorptr++);
					(*psd_vectorptr++)=norm*(real*real+imag*imag);;
				}
			};

			void write_ascii(const char *filename) const
			{
				FILE *fileptr=fopen(filename,"w");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clVector_complex::write_ascii\n");
					exit(1);
				}
				write_ascii(fileptr);
				fclose(fileptr);
			};
			void write_ascii(FILE *fileptr) const
			{
				double *vectorptr=_vector;
				fprintf(fileptr,"%d\n",_dim);
				for (int n=0; n<_dim; n++) {
					fprintf(fileptr," %12.4le",(*vectorptr++));
					fprintf(fileptr," %12.4le\n",(*vectorptr++));
				}
				fprintf(fileptr,"\n");
			};

		private:

#ifdef DEBUG
			void test_index(const int n) const
			{
				if (n<0 || n>=_dim) {
					fprintf(stdout,"FATAL ERROR: index out of range in clVector_complex::test_index %d %d\n",n,_dim);
					exit(1);
				}
			};
			void test_indices(const clVector_complex &arg) const
			{
				if (!_vector || !arg._vector) {
					fprintf(stdout,"FATAL ERROR: _vector or arg._vector not allocated in clVector_complex::test_indices\n");
					exit(1);
				}
				if (_dim!=arg._dim) {
					fprintf(stdout,"WARNING: different source and target dimension in clVector_complex::test_indices\n");
					fprintf(stdout,"         use smaller dimension of %d and %d\n",_dim,arg._dim);
				}
			};
#endif

	};

	/*  clMatrix_complex class  */
	class clMatrix_complex
	{
		private:

			int _dim1,_dim2,_dim;		// Dimensions of complex elements in the matrix
			int _dimmat;				// Actual storage space of _vector (2*_dim)
			double *_startptr;			// Matrix start pointer
			clVector_complex _matrix;	// Matrix stored as a single vector
#ifdef ALLOCDEBUG
			static int _num_matrix_complex;
#endif

		public:

			/*  Constructors  */
			clMatrix_complex(void) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_dimmat(0),
				_startptr(NULL)
			{ };
			explicit clMatrix_complex(const int dim1,const int dim2) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_dimmat(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2);
			};
			explicit clMatrix_complex(const int dim1,const int dim2,const double real,const double imag=0.0) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_dimmat(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2,real,imag);
			};
			explicit clMatrix_complex(const int dim1,const int dim2,const clComplex &arg) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_dimmat(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2,arg);
			};

			/*  Destructor  */
			virtual ~clMatrix_complex(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clMatrix_complex(const clMatrix_complex &mm) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_dimmat(0),
				_startptr(NULL)
			{
				_dim1=mm._dim1;
				_dim2=mm._dim2;
				_dim=mm._dim;
				_dimmat=mm._dimmat;
				_matrix=mm._matrix;
				_startptr=_matrix.get_startptr();
			};

			/*  Default Assignment Constructor  */
			clMatrix_complex& operator= (const clMatrix_complex &mm)
			{
				if (this==&mm) return *this;

				if (!_startptr) {
					if (_dim1 || _dim2 || _dim || _dimmat) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clMatrix_complex::operator=\n");
						exit(1);
					}
					allocate(mm._dim1,mm._dim2);
				} else {
					if (_dim1!=mm._dim1 || _dim2!=mm._dim2 || _dim!=mm._dim || _dimmat!=mm._dimmat) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clMatrix_complex::operator=\n");
						exit(1);
					}
				}
				_matrix=mm._matrix;
				return *this;
			};

			/*  Operator overloading functions  */

			double* operator()(const int n1)
			/* This returns a pointer to the start of the n1-th row */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator()\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+2*_dim2*n1);
			};
			double* operator()(const int n1) const
			/* This returns a pointer to the start of the n1-th row */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator()\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+2*_dim2*n1);
			};

			double* operator()(const int n1,const int n2)
			/* This returns a pointer to n1,n2 real element */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index1(n2);
#endif
				return (_startptr+2*(n2+_dim2*n1));
			};
			double* operator()(const int n1,const int n2) const
			/* This returns a pointer to n1,n2 real element */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index1(n2);
#endif
				return (_startptr+2*(n2+_dim2*n1));
			};

			void operator()(const int n1,const int n2,const clComplex &arg)
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator[]\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				double *vectorptr=_startptr+2*(n2+_dim2*n1);
				(*vectorptr++)=arg._real;
				(*vectorptr)=arg._imag;
			};
			void operator()(const int n1,const int n2,const double real,const double imag=0.0)
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator[]\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				double *vectorptr=_startptr+2*(n2+_dim2*n1);
				(*vectorptr++)=real;
				(*vectorptr)=imag;
			};

			clMatrix_complex& operator= (const double arg)
			{
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator=\n");
					exit(1);
				}
				double *matrixptr=_startptr;
				for (int n=0; n<_dim; n++,matrixptr+=2) (*matrixptr)=arg;
				return *this;
			};
			clMatrix_complex& operator= (const clComplex &arg)
			{
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix_complex::operator=\n");
					exit(1);
				}
				const double real=arg._real;
				const double imag=arg._imag;
				double *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) {
					(*matrixptr++)=real;
					(*matrixptr++)=imag;
				}
				return *this;
			};
			/*  Converts a regular 2D array to a matrix, does not allow bounce checking  */
			clMatrix_complex& operator= (const double* const *mm)
			{
				if (!mm || !mm[0]) {
					fprintf(stdout,"FATAL ERROR: mm not allocated in clMatrix_complex::operator=\n");
					exit(1);
				}
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: _matrix not allocated in clMatrix_complex::operator=\n");
					exit(1);
				}
				double *matrixptr=_startptr;
				const double *mmptr=mm[0];
				for (int n=0; n<_dimmat; n++) (*matrixptr++)=(*mmptr++);
				return *this;
			};

			clMatrix_complex operator+ (const clMatrix_complex &mm) const
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				clMatrix_complex sum(_dim1,_dim2);
				double *sumptr=sum._startptr;
				const double *matrixptr=_startptr;
				const double *mmptr=mm._startptr;
				for (int n=0; n<_dimmat; n++) (*sumptr++)=(*matrixptr++)+(*mmptr++);
				return sum;
			};
			clMatrix_complex operator+ (const double arg) const
			{
				clMatrix_complex sum(_dim1,_dim2);
				double *sumptr=sum._startptr;
				const double *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) {
					(*sumptr++)=(*matrixptr++)+arg;
					(*sumptr++)=(*matrixptr++);
				}
				return sum;
			};
			clMatrix_complex operator+ (const clComplex &arg) const
			{
				const double real=arg._real;
				const double imag=arg._imag;
				clMatrix_complex sum(_dim1,_dim2);
				double *sumptr=sum._startptr;
				const double *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) {
					(*sumptr++)=(*matrixptr++)+real;
					(*sumptr++)=(*matrixptr++)+imag;
				}
				return sum;
			};
			friend clMatrix_complex operator+ (const double arg,const clMatrix_complex &mm)
			{
				clMatrix_complex sum(mm._dim1,mm._dim2);
				double *sumptr=sum._startptr;
				const double *matrixptr=mm._startptr;
				for (int n=0; n<mm._dim; n++) {
					(*sumptr++)=(*matrixptr++)+arg;
					(*sumptr++)=(*matrixptr++);
				}
				return sum;
			};
			friend clMatrix_complex operator+ (const clComplex &arg,const clMatrix_complex &mm)
			{
				const double real=arg.real();
				const double imag=arg.imag();
				clMatrix_complex sum(mm._dim1,mm._dim2);
				double *sumptr=sum._startptr;
				const double *matrixptr=mm._startptr;
				for (int n=0; n<mm._dim; n++) {
					(*sumptr++)=(*matrixptr++)+real;
					(*sumptr++)=(*matrixptr++)+imag;
				}
				return sum;
			};

			/*  Remaining member functions  */

			int get_dim1(void) {return _dim1;};
			int get_dim1(void) const {return _dim1;};
			int get_dim2(void) {return _dim2;};
			int get_dim2(void) const {return _dim2;};
			int get_dim(void) {return _dim;};
			int get_dim(void) const {return _dim;};
			double* get_startptr(const int n1=0) {return _startptr+2*_dim2*n1;};
			const double* get_startptr(const int n1=0) const {return _startptr+2*_dim2*n1;};
			double* get_startptr(const int n1,const int n2) {return _startptr+2*(n2+_dim2*n1);};
			const double* get_startptr(const int n1,const int n2) const {return _startptr+2*(n2+_dim2*n1);};

			/*  Converts a matrix to a regular 2D array, does not allow bounce checking  */
			void copy(double **mm) const
			{
				if (!mm || !mm[0]) {
					fprintf(stdout,"FATAL ERROR: mm not allocated in clMatrix_complex::copy\n");
					exit(1);
				}
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: _startptr not allocated in clMatrix_complex::copy\n");
					exit(1);
				}
				const double *matrixptr=_startptr;
				double *mmptr=mm[0];
				for (int n=0; n<_dim; n++) (*mmptr++)=(*matrixptr++);
			};

			/*  Allocate matrix of dimension (dim1,dim2)  */
			double* allocate(const int dim1,const int dim2)
			{
				if (_dim==dim1*dim2) {
					_dim1=dim1;
					_dim2=dim2;
				} else {
					if (_startptr) deallocate();
					_dim=dim1*dim2;
					if (!_dim) return _startptr;
					_dim1=dim1;
					_dim2=dim2;
					_dimmat=2*_dim;
					_matrix.allocate(_dim);
					_startptr=_matrix.get_startptr();
#ifdef ALLOCDEBUG
					_num_matrix_complex++;
					fprintf(stdout,"complex matrix %p allocated %d size %dx%d\n",_startptr,_num_matrix_complex,dim1,dim2);
#endif
				}
				return _startptr;
			};
			/*  Allocate matrix of dimension (dim1,dim2) and initialize to initreal  */
			double* allocate(const int dim1,const int dim2,const double initreal,const double initimag=0.0)
			{
				clComplex initval(initreal,initimag);
				allocate(dim1,dim2,initval);
				return _startptr;
			};
			double* allocate(const int dim1,const int dim2,const clComplex &initval)
			{
				allocate(dim1,dim2);
				if (_startptr) _matrix=initval;
				return _startptr;
			};

			/*  Deallocate matrix  */
			void deallocate(void)
			{
				if (_startptr) {
#ifdef ALLOCDEBUG
					_num_matrix_complex--;
					fprintf(stdout,"complex matrix %p deallocated %d size %dx%d\n",_startptr,_num_matrix_complex,0.5*_dim1,0.5*_dim2);
#endif
					_startptr=NULL;
				}
				_matrix.deallocate();
				_dim1=0;
				_dim2=0;
				_dim=0;
			};

			void write_ascii_transpose(const char *filename) const
			{
				FILE *fileptr=fopen(filename,"w");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clMatrix_complex::write_ascii_transpose\n");
					exit(1);
				}
				write_ascii_transpose(fileptr);
				fclose(fileptr);
			};
			void write_ascii_transpose(FILE *fileptr) const
			{
				int n1,n2;
				const double *matrixptr=NULL;
				for (n2=0; n2<_dim2; n2++) {
					matrixptr=get_startptr(n1,n2);
					for (n1=0; n1<_dim1; n1++) {
						fprintf(fileptr," %12.4le",(*matrixptr++));
						fprintf(fileptr," %12.4le",(*matrixptr));
					}
					fprintf(fileptr,"\n");
				}
			};

			void write_ascii(const char *filename) const
			{
				FILE *fileptr=fopen(filename,"w");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clMatrix_complex::write_ascii\n");
					exit(1);
				}
				write_ascii(fileptr);
				fclose(fileptr);
			};
			void write_ascii(FILE *fileptr) const
			{
				int n1,n2;
				double *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++) {
					for (n2=0; n2<_dim2; n2++) {
						fprintf(fileptr," %12.4le",(*matrixptr++));
						fprintf(fileptr," %12.4le",(*matrixptr++));
					}
					fprintf(fileptr,"\n");
				}
			};

		private:

#ifdef DEBUG
			void test_index1(const int n1) const
			{
				if (n1<0 || n1>=_dim1) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix_complex::test_index1: %d %d\n",n1,_dim1);
					exit(1);
				}
			};
			void test_index2(const int n2) const
			{
				if (n2<0 || n2>=_dim2) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix_complex::test_index2: %d %d\n",n2,_dim2);
					exit(1);
				}
			};
			void test_indices(const clMatrix_complex &mm) const
			{
				if (!_startptr || !mm._startptr) {
					fprintf(stdout,"FATAL ERROR: matrix or mm.matrix not allocated in clMatrix_complex::test_indices\n");
					exit(1);
				}
				if (_dim1!=mm._dim1 || _dim2!=mm._dim2) {
					fprintf(stdout,"FATAL ERROR: different source and target dimension in clMatrix_complex::test_indices\n");
				}
			};
#endif

	};

	/* clVector_complex functions */

	void bcast(clVector_complex&,const int ip=0);

}; // namespace clComplex_namespace

#endif // CLCOMPLEX_H__INCLUDED_
