/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/

#ifndef CLVECTOR_H__INCLUDED_
#define CLVECTOR_H__INCLUDED_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

namespace clVector_namespace
{

	/*  clVector class  */
	template <class T> class clVector
	{
		private:

			int _dim;	// Dimension
			T *_vector;	// Vector
#ifdef ALLOCDEBUG
			static int _num_vector;
#endif

		public:

			/*  Constructors  */
			clVector(void) :
				_dim(0),
				_vector(NULL)
			{ };
			explicit clVector(int dim) :
				_dim(0),
				_vector(NULL)
			{
				allocate(dim);
			};
			clVector(int dim,const T initval) :
				_dim(0),
				_vector(NULL)
			{
				allocate(dim,initval);
			};

			/*  Destructor  */
			virtual ~clVector(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clVector(const clVector &arg) :
				_dim(0),
				_vector(NULL)
			{
				allocate(arg._dim);
				T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)=(*argvectorptr++);
			};

			/*  Default Assignment Constructor  */
			clVector& operator= (const clVector &arg)
			{
				if (this==&arg) return *this;

				if (!_vector) {
					if (_dim) {
						fprintf(stdout,"FATAL ERROR: illegal _dim in clVector::operator=\n");
						exit(1);
					} else if (!arg._dim) {
#ifdef DEBUG
						fprintf(stdout,"WARNING: empty argument vector in clVector::operator=\n");
#endif
						return *this;
					}
					allocate(arg._dim);
				}
#ifdef DEBUG
				test_indices(arg);
#endif
				int dim=(_dim <= arg._dim ? _dim : arg._dim);
				T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<dim; n++) (*vectorptr++)=(*argvectorptr++);
				return *this;
			};

			/*  Operator overloading functions  */

			T& operator[](const int n)
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector::operator[]\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};

			T& operator[](const int n) const
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector::operator[]\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};

			T& operator()(const int n)
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector::operator()\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};
			T& operator()(const int n) const
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector::operator()\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};

			clVector& operator= (const T arg)
			{
				if (!_vector) {
					fprintf(stdout,"WARNING: _vector not allocated in clVector::operator=\n");
					return *this;
				}
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)=arg;
				return *this;
			};
			/*  Converts a regular array to a vector, does not allow bounce checking  */
			clVector& operator= (const T *arg)
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg not allocated in clVector::operator=\n");
					exit(1);
				}
				if (!_vector) {
					fprintf(stdout,"WARNING: _vector not allocated in clVector::operator=\n");
					return *this;
				}
				T *vectorptr=_vector;
				const T *argptr=arg;
				for (int n=0; n<_dim; n++) (*vectorptr++)=(*argptr++);
				return *this;
			};

			clVector operator+ (const clVector &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				int dim=(_dim <= arg._dim ? _dim : arg._dim);
				clVector<T> sum(dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<dim; n++) (*sumptr++)=(*vectorptr++)+(*argvectorptr++);
				return sum;
			};
			clVector operator+ (const T arg) const
			{
				if (!arg) return *this;
				clVector<T> sum(_dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*vectorptr++)+arg;
				return sum;
			};
			friend clVector<T> operator+ (const T arg1,const clVector<T>& arg2)
			{
				if (!arg1) return arg2;
				clVector<T> sum(arg2._dim);
				T *sumptr=sum._vector;
				const T *arg2ptr=arg2._vector;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1+(*arg2ptr++);
				return sum;
			};

			clVector operator- (const clVector &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				int dim=(_dim <= arg._dim ? _dim : arg._dim);
				clVector<T> sum(dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<dim; n++) (*sumptr++)=(*vectorptr++)-(*argvectorptr++);
				return sum;
			};
			clVector operator- (const T arg) const
			{
				if (!arg) return *this;
				clVector<T> sum(_dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*vectorptr++)-arg;
				return sum;
			};
			friend clVector<T> operator- (const T arg1,const clVector<T>& arg2)
			{
				if (!arg1) return arg2;
				clVector<T> sum(arg2._dim);
				T *sumptr=sum._vector;
				const T *arg2ptr=arg2._vector;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1-(*arg2ptr++);
				return sum;
			};

			clVector operator* (const clVector &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				int dim=(_dim <= arg._dim ? _dim : arg._dim);
				clVector<T> sum(dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<dim; n++) (*sumptr++)=(*vectorptr++)*(*argvectorptr++);
				return sum;
			};
			clVector operator* (const T arg) const
			{
				if (arg==1.0) return *this;
				clVector<T> sum(_dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*vectorptr++)*arg;
				return sum;
			};
			friend clVector<T> operator* (const T arg1,const clVector<T>& arg2)
			{
				if (arg1==1.0) return arg2;
				clVector<T> sum(arg2._dim);
				T *sumptr=sum._vector;
				const T *arg2ptr=arg2._vector;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1*(*arg2ptr++);
				return sum;
			};

			clVector operator/ (const clVector &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				int dim=(_dim <= arg._dim ? _dim : arg._dim);
				clVector<T> sum(dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<dim; n++) {
					if (!(*argvectorptr)) {
						fprintf(stdout,"FATAL ERROR: arg=0 in clVector::operator/\n");
						exit(1);
					}
					(*sumptr++)=(*vectorptr++)/(*argvectorptr++);
				}
				return sum;
			};
			clVector operator/ (const T arg) const
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clVector::operator/\n");
					exit(1);
				}
				clVector<T> sum(_dim);
				T *sumptr=sum._vector;
				const T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*vectorptr++)/arg;
				return sum;
			};
			friend clVector<T> operator/ (const T arg1,const clVector<T>& arg2)
			{
				clVector<T> sum(arg2._dim);
				T *sumptr=sum._vector;
				const T *arg2ptr=arg2._vector;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1/(*arg2ptr++);
				return sum;
			};

			clVector& operator+= (const clVector &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)+=(*argvectorptr++);
				return *this;
			};
			clVector& operator+= (const T arg)
			{
				if (!arg) return *this;
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)+=arg;
				return *this;
			};

			clVector& operator-= (const clVector &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)-=(*argvectorptr++);
				return *this;
			};
			clVector& operator-= (const T arg)
			{
				if (!arg) return *this;
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)-=arg;
				return *this;
			};

			clVector& operator*= (const clVector &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)*=(*argvectorptr++);
				return *this;
			};
			clVector& operator*= (const T arg)
			{
				if (arg==1.0) return *this;
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)*=arg;
				return *this;
			};

			clVector& operator/= (const clVector &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *vectorptr=_vector;
				const T *argvectorptr=arg._vector;
				for (int n=0; n<_dim; n++) {
					if (!(*argvectorptr)) {
						fprintf(stdout,"FATAL ERROR: arg=0 in clVector::operator/=\n");
						exit(1);
					}
					(*vectorptr++)/=(*argvectorptr++);
				}
				return *this;
			};
			clVector& operator/= (const T arg)
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clVector::operator/=\n");
					exit(1);
				}
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)/=arg;
				return *this;
			};

			clVector& operator++ (void)
			{
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)++;
				return *this;
			};

			clVector& operator-- (void)
			{
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)--;
				return *this;
			};

			/*  Remaining member functions  */

			int get_dim(void) {return _dim;};
			int get_dim(void) const {return _dim;};
			T *get_startptr(void) {return _vector;};
			const T *get_startptr(void) const {return _vector;};

			/*  Move data from arg to current */
			void move(clVector &arg)
			{
				if (this==&arg) return;

				if (_vector) deallocate();

				_dim=arg._dim;
				_vector=arg._vector;

				arg._dim=0;
				arg._vector=NULL;
			};

			/*  Swap data content between current and arg, dimensions must match */
			void swap(clVector &arg)
			{
				if (this==&arg) return;

				if (_dim!=arg._dim) {
					fprintf(stdout,"FATAL ERROR: inconsistent dimension in clVector::swap\n");
					exit(1);
				}
				T *swap=_vector;
				_vector=arg._vector;
				arg._vector=swap;
			};

			/*  Converts a vector to a regular array, does not allow bounce checking  */
			void copy(T *arg) const
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg not allocated in clVector::copy\n");
					exit(1);
				}
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clVector::copy\n");
					exit(1);
				}
				const T *vectorptr=_vector;
				T *argptr=arg;
				for (int n=0; n<_dim; n++) (*argptr++)=(*vectorptr++);
			};

			/*  Allocate vector of dimension dim  */
			T* allocate(const int dim)
			{
				if (_dim!=dim) {
					if (_vector) deallocate();
					_dim=dim;
					if (!_dim) return _vector;
					_vector=new T [_dim];
					if (!_vector) {
						fprintf(stdout,"FATAL ERROR: Could not allocate _vector in clVector::allocate\n");
						exit(1);
					}
#ifdef ALLOCDEBUG
					_num_vector++;
					fprintf(stdout,"vector %p size %d allocated %d\n",_vector,_dim,_num_vector);
#endif
				}
				return _vector;
			};
			/*  Allocate vector of dimension dim and initialize to initval  */
			T* allocate(const int dim,const T initval)
			{
				allocate(dim);
				if (_vector) {
					T *vectorptr=_vector;
					for (int n=0; n<_dim; n++) (*vectorptr++)=initval;
				}
				return _vector;
			};

			/*  Resize a vector to a new dimension new_dim  */
			/*  Data that does not fit is lost  */
			/*  Zero out cells that are not used yet */
			T* redimension(const int new_dim)
			{
				if (!new_dim) {
					deallocate();
					return _vector;
				}
				if (!_vector) {
					allocate(new_dim,(T)0.0);
					return _vector;
				}
				if (!_dim) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::redimension\n");
					exit(1);
				}
				T *newvector=new T [new_dim];
				if (!newvector) {
					fprintf(stdout,"FATAL ERROR: Could not allocate newvector in clVector::redimension\n");
					exit(1);
				}
				const int dim=(_dim <= new_dim ? _dim : new_dim);
				T *newvectorptr=newvector;
				const T *vectorptr=_vector;
				int n=0;
				for (; n<dim; n++) (*newvectorptr++)=(*vectorptr++);
				for (; n<new_dim; n++) (*newvectorptr++)=(T)0.0;
				T *swap=_vector;
				_vector=newvector;
				_dim=new_dim;
				delete [] swap;
				return _vector;
			};

			/*  Deallocate vector  */
			void deallocate(void)
			{
				if (_vector) {
#ifdef ALLOCDEBUG
					_num_vector--;
					fprintf(stdout,"vector %p size %d deallocated %d\n",_vector,_dim,_num_vector);
#endif
					delete [] _vector;
					_vector=NULL;
				}
				_dim=0;
			};

			/*  Allocate vector of dimension dim and initialize to (0,1,2,...,dim-1)  */
			clVector& index_offsetzero(const int dim)
			{
				allocate(dim);
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)=(T)(n);
				return *this;
			};
			/*  Initialize existing vector to (0,1,2,...,dim-1)  */
			clVector& index_offsetzero(void)
			{
				T *vectorptr=_vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)=(T)(n);
				return *this;
			};
			/*  Allocate vector of dimension dim and initialize to (1,2,3,...,dim)  */
			clVector& index(const int dim)
			{
				allocate(dim);
				T *vectorptr=_vector;
				for (int n=1; n<=_dim; n++) (*vectorptr++)=(T)(n);
				return *this;
			};
			/*  Initialize existing vector to (1,2,3,...,dim)  */
			clVector& index(void)
			{
				T *vectorptr=_vector;
				for (int n=1; n<=_dim; n++) (*vectorptr++)=(T)(n);
				return *this;
			};

			/*  Normalize an existing vector  */
			void normalize(void)
			{
				if (!_dim) return;
				T *vectorptr=_vector;
				double nrm=norm();
				if (nrm) {
					nrm=(1.0)/sqrt(nrm);
					for (int n=0; n<_dim; n++) (*vectorptr++)*=(T)(nrm);
				}
			};

			/*  Return L1-norm of vector  */
			T sum(void) const
			{
				int dim=0;
				return sum(dim);
			}
			T sum(int dim,const int offset=0) const
			{
				if (!_dim) return 0;
				if (dim<0) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::sum\n");
					exit(1);
				}
				if (offset<0) {
					fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::sum\n");
					exit(1);
				}
				if (!dim) {
					if (offset) {
						fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::sum\n");
						exit(1);
					}
					dim=_dim;
				}
				dim+=offset;
				if (dim>_dim) dim=_dim;
				const T *vectorptr=_vector+offset;
				T ssum=0;
				for (int n=offset; n<dim; n++) ssum+=(*vectorptr++);
				return ssum;
			};

			/*  Return L2-norm of vector  */
			T norm(void) const
			{
				int dim=0;
				return norm(dim);
			}
			T norm(int dim,const int offset=0) const
			{
				if (!_dim) return 0;
				if (dim<0) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::norm\n");
					exit(1);
				}
				if (offset<0) {
					fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::norm\n");
					exit(1);
				}
				if (!dim) {
					if (offset) {
						fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::norm\n");
						exit(1);
					}
					dim=_dim;
				}
				dim+=offset;
				if (dim>_dim) dim=_dim;
				const T *vectorptr=_vector+offset;
				T nnorm=0;
				for (int n=offset; n<dim; n++,vectorptr++) nnorm+=(*vectorptr)*(*vectorptr);
				return nnorm;
			};

			/*  Return RMS of vector  */
			double rms(int dim=0,const int offset=0) const
			{
				if (!_dim) return 0;
				if (dim<0) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::rms\n");
					exit(1);
				}
				if (offset<0) {
					fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms\n");
					exit(1);
				}
				if (!dim) {
					if (offset) {
						fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms\n");
						exit(1);
					}
					dim=_dim;
				}
				dim+=offset;
				if (dim>_dim) dim=_dim;
				const T *vectorptr=_vector+offset;
				double sumsq=0;
				for (int n=offset; n<dim; n++,vectorptr++) sumsq+=(double)((*vectorptr)*(*vectorptr));
				return sqrt(sumsq/((double)(dim-offset)));
			};

			/*  Return RMS of vector, ignore zero entries  */
			double rms_skipzero(int dim=0,const int offset=0) const
			{
				if (!_dim) return 0;
				if (dim<0) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::rms_skipzero\n");
					exit(1);
				}
				if (offset<0) {
					fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms_skipzero\n");
					exit(1);
				}
				if (!dim) {
					if (offset) {
						fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms_skipzero\n");
						exit(1);
					}
					dim=_dim;
				}
				dim+=offset;
				if (dim>_dim) dim=_dim;
				const T *vectorptr=_vector+offset;
				int nn=0;
				double sumsq=0;
				for (int n=offset; n<dim; n++,vectorptr++) {
					if (*vectorptr) {
						nn++;
						sumsq+=(double)((*vectorptr)*(*vectorptr));
					}
				}
				if (!nn) return 0;
				else return sqrt(sumsq/((double)(nn)));
			};

			/*  Return RMS error between the current vector and the one passed to the routine  */
			double rms_error(const clVector<T> &arg,int dim=0,const int offset=0) const
			{
				if (!_dim) return 0;
				if (dim<0) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::rms_error\n");
					exit(1);
				}
				if (offset<0) {
					fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms_error\n");
					exit(1);
				}
				if (!dim) {
					if (offset) {
						fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms_error\n");
						exit(1);
					}
					dim=_dim;
				}
				if (!arg.get_dim()) return 0.0;
				dim+=offset;
				if (dim>_dim) dim=_dim;
				if (dim>arg.get_dim()) dim=arg.get_dim();
				const T *vectorptr=_vector+offset;
				const T *argvectorptr=arg.get_startptr()+offset;
				double dummy,sumsq=0;
				for (int n=offset; n<dim; n++) {
					dummy=(double)((*vectorptr++)-(*argvectorptr++));
					sumsq+=dummy*dummy;
				}
				return sqrt(sumsq/(double)((dim-offset)));
			};

			/*  Return RMS error between the current vector and the one passed to the routine */
			/*    skip entries where the current vector is zero  */
			double rms_error_skipzero(const clVector<T> &arg,int dim=0,const int offset=0) const
			{
				if (!_dim) return 0;
				if (dim<0) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::rms_error_skipzero\n");
					exit(1);
				}
				if (offset<0) {
					fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms_error_skipzero\n");
					exit(1);
				}
				if (!dim) {
					if (offset) {
						fprintf(stdout,"FATAL ERROR: Illegal offset in clVector::rms_error_skipzero\n");
						exit(1);
					}
					dim=_dim;
				}
				if (!arg.get_dim()) return 0.0;
				dim+=offset;
				if (dim>_dim) dim=_dim;
				if (dim>arg.get_dim()) dim=arg.get_dim();
				const T *vectorptr=_vector+offset;
				const T *argvectorptr=arg.get_startptr()+offset;
				int nn=0;
				double dummy,sumsq=0;
				for (int n=offset; n<dim; n++,vectorptr++,argvectorptr++) {
					if (*vectorptr) {
						nn++;
						dummy=(double)((*vectorptr)-(*argvectorptr));
						sumsq+=dummy*dummy;
					}
				}
				if (!nn) return 0;
				else return sqrt(sumsq/((double)(nn)));
			};

			/*  Return max of vector  */
			T max(void) const
			{
				if (!_dim) return 0;
				const T *vectorptr=_vector;
				T mmax=(*vectorptr++);
				for (int n=1; n<_dim; n++,vectorptr++) if ((*vectorptr)>mmax) mmax=(*vectorptr);
				return mmax;
			};
			/*  Return max of vector in subrange  */
			T max(const int dimlow,const int dimupp) const
			{
				if (!_dim || dimlow>dimupp) return 0;
#ifdef DEBUG
				test_index(dimlow);
				test_index(dimupp);
#endif
				int num=dimupp-dimlow+1;
				const T *vectorptr=_vector+dimlow;
				T mmax=(*vectorptr++);
				for (int n=1; n<num; n++,vectorptr++) if ((*vectorptr)>mmax) mmax=(*vectorptr);
				return mmax;
			};
			/*  Return index of max of vector  */
			int maxindex(void) const
			{
				if (!_dim) return -1;
				const T *vectorptr=_vector;
				int mmaxindex=0;
				T mmax=(*vectorptr++);
				for (int n=1; n<_dim; n++,vectorptr++) if ((*vectorptr)>mmax) {
					mmax=(*vectorptr);
					mmaxindex=n;
				}
				return mmaxindex;
			};
			/*  Return max of vector and its index  */
			void max(T &mmax,int &mmaxindex) const
			{
				if (!_dim) {
					mmax=0;
					mmaxindex=-1;
					return;
				}
				const T *vectorptr=_vector;
				mmaxindex=0;
				mmax=(*vectorptr++);
				for (int n=1; n<_dim; n++,vectorptr++) if ((*vectorptr)>mmax) {
					mmax=(*vectorptr);
					mmaxindex=n;
				}
			};

			/*  Return min of vector  */
			T min(void) const
			{
				if (!_dim) return 0;
				const T *vectorptr=_vector;
				T mmin=(*vectorptr++);
				for (int n=1; n<_dim; n++,vectorptr++) if ((*vectorptr)<mmin) mmin=(*vectorptr);
				return mmin;
			};
			/*  Return min of vector in subrange  */
			T min(const int dimlow,const int dimupp) const
			{
				if (!_dim || dimlow>dimupp) return 0;
#ifdef DEBUG
				test_index(dimlow);
				test_index(dimupp);
#endif
				int num=dimupp-dimlow+1;
				const T *vectorptr=_vector+dimlow;
				T mmin=(*vectorptr++);
				for (int n=1; n<num; n++,vectorptr++) if ((*vectorptr)<mmin) mmin=(*vectorptr);
				return mmin;
			};
			/*  Return index of min of vector  */
			int minindex(void) const
			{
				if (!_dim) return -1;
				const T *vectorptr=_vector;
				int mminindex=0;
				T mmin=(*vectorptr++);
				for (int n=1; n<_dim; n++,vectorptr++) if ((*vectorptr)<mmin) {
					mmin=(*vectorptr);
					mminindex=n;
				}
				return mminindex;
			};
			/*  Return min of vector and its index  */
			void min(T &mmin,int &mminindex) const
			{
				if (!_dim) {
					mmin=0;
					mminindex=-1;
					return;
				}
				const T *vectorptr=_vector;
				mminindex=0;
				mmin=(*vectorptr++);
				for (int n=1; n<_dim; n++,vectorptr++) if ((*vectorptr)<mmin) {
					mmin=(*vectorptr);
					mminindex=n;
				}
			};

			/*  Return max of the absolute value of vector  */
			T maxabs(void) const
			{
				if (!_dim) return 0;
				const T *vectorptr=_vector;
				T mmax=(*vectorptr++);
				if (mmax<0.0) mmax=-mmax;
				for (int n=1; n<_dim; n++,vectorptr++) {
					if ((*vectorptr)>mmax) mmax=(*vectorptr);
					else if ((*vectorptr)<-mmax) mmax=-(*vectorptr);
				}
				return mmax;
			};
			/*  Return index of max of the absolute value of vector  */
			int maxabsindex(void) const
			{
				if (!_dim) return -1;
				const T *vectorptr=_vector;
				int mmaxindex=0;
				T mmax=(*vectorptr++);
				if (mmax<0.0) mmax=-mmax;
				for (int n=1; n<_dim; n++,vectorptr++) {
					if ((*vectorptr)>mmax) {
						mmax=(*vectorptr);
						mmaxindex=n;
					} else if ((*vectorptr)<-mmax) {
						mmax=-(*vectorptr);
						mmaxindex=n;
					}
				}
				return mmaxindex;
			};

			/*  Return mean value of vector  */
			double mean(void) const
			{
				if (!_dim) return 0;
				const T *vectorptr=_vector;
				double mmean=0;
				for (int n=0; n<_dim; n++) mmean+=(*vectorptr++);
				return (mmean/_dim);
			};

			/*  Return standard deviation of vector  */
			double stdv(void) const
			{
				if (_dim<2) return 0;
				double mmean=mean();
				const T *vectorptr=_vector;
				double dummy,sstdv=0;
				for (int n=0; n<_dim; n++) {
					dummy=((*vectorptr++)-mmean);
					sstdv+=dummy*dummy;
				}
				return (sqrt(sstdv/_dim));
			};

			/*  Get the outerproduct of two 3D vectors and store in this */
			void outerproduct(const clVector<T> &arg1,const clVector<T> &arg2)
			{
				if (_dim!=3 || arg1.get_dim()!=3 || arg2.get_dim()!=3) {
					fprintf(stdout,"FATAL ERROR: illegal dimension in clVector::outerproduct\n");
					exit(1);
				}
				_vector[0]=arg1(1)*arg2(2)-arg1(2)*arg2(1);
				_vector[1]=arg1(2)*arg2(0)-arg1(0)*arg2(2);
				_vector[2]=arg1(0)*arg2(1)-arg1(1)*arg2(0);
			};

			/*  Read a vector from an ascii file  */
			int read_ascii(const char *filename)
			{
				FILE *fileptr=fopen(filename,"r");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clVector::read_ascii\n");
					exit(1);
				}
				int error_flag=read_ascii(fileptr);
				fclose(fileptr);
				return error_flag;
			};
			/*  Read a vector from an ascii stream  */
			int read_ascii(FILE *fileptr)
			{
				int dim;
				if (fscanf(fileptr,"%d\n",&dim)!=1) return 1;
				if (dim) {
					allocate(dim);
					T *vectorptr=_vector;
					for (int n=0; n<_dim; n++) if (fscanf(fileptr,"%le",(vectorptr++))!=1) return 1;
				}
				return 0;
			};

			/*  Write the transpose of a vector to an ascii file  */
			/*    append_flag=0:  Create a new file (default)  */
			/*    append_flag!=0: Append to file if existing  */
			void write_ascii_transpose(const char *filename,const int append_flag=0) const
			{
				FILE *fileptr=NULL;
				if (!append_flag) fileptr=fopen(filename,"w");
				else fileptr=fopen(filename,"a");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clVector::write_ascii_transpose\n");
					exit(1);
				}
				write_ascii_transpose(fileptr);
				fclose(fileptr);
			};
			/*  Write the transpose of a vector to an ascii stream  */
			void write_ascii_transpose(FILE *fileptr) const
			{
				T *vectorptr=_vector;
				fprintf(fileptr,"%d\n",_dim);
//				for (int n=0; n<_dim; n++) fprintf(fileptr," %12.4le\n",(double)(*vectorptr++));
				for (int n=0; n<_dim; n++) fprintf(fileptr," %24.16le\n",(double)(*vectorptr++));
			};

			/*  Write a vector to an ascii file  */
			/*    append_flag=0:  Create a new file (default)  */
			/*    append_flag!=0: Append to file if existing  */
			/*    full_precision!=0: Write in %16.8le format (default)  */
			/*    full_precision=0:  write in full precision (default)  */
			void write_ascii(const char *filename,const int append_flag=0,const int full_precision=0) const
			{
				FILE *fileptr=NULL;
				if (append_flag) fileptr=fopen(filename,"a");
				else fileptr=fopen(filename,"w");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clVector::write_ascii\n");
					exit(1);
				}
				write_ascii(fileptr,full_precision);
				fclose(fileptr);
			};
			/*  Write a vector to an ascii stream  */
			/*    full_precision!=0: Write in %16.8le format (default)  */
			/*    full_precision=0:  write in full precision (default)  */
			void write_ascii(FILE *fileptr,const int full_precision=0) const
			{
				T *vectorptr=_vector;
				fprintf(fileptr,"%d\n",_dim);
				if (full_precision) {
					for (int n=0; n<_dim; n++) fprintf(fileptr," %24.16le",(double)(*vectorptr++));
				} else {
					for (int n=0; n<_dim; n++) fprintf(fileptr," %16.8le",(double)(*vectorptr++));
				}
				fprintf(fileptr,"\n");
			};

			/*  Read a vector from a binary file  */
			int read(const char *filename)
			{
				FILE *fileptr=fopen(filename,"rb");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clVector::read\n");
					exit(1);
				}
				const int error_flag=read(fileptr);
				fclose(fileptr);
				return error_flag;
			};
			/*  Read a vector from a binary stream  */
			int read(FILE *fileptr)
			{
				int dim;
				void *io_ptr;
				io_ptr=(void *)&dim;
				if (fread(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
				if (dim) {
					allocate(dim);
					io_ptr=(void *)_vector;
					if (fread(io_ptr,sizeof(T),_dim,fileptr)!=(unsigned)_dim) return 0;
				}
				return 1;
			};

			/*  Write a vector to a binary file  */
			int write(const char *filename) const
			{
				FILE *fileptr=fopen(filename,"wb");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clVector::write\n");
					exit(1);
				}
				const int error_flag=write(fileptr);
				fclose(fileptr);
				return error_flag;
			};
			/*  Write a vector to a binary stream  */
			int write(FILE *fileptr) const
			{
				void *io_ptr;
				io_ptr=(void *)&_dim;
				if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
				if (_dim) {
					io_ptr=(void *)_vector;
					if (fwrite(io_ptr,sizeof(T),_dim,fileptr)!=(unsigned)_dim) return 0;
				}
				return 1;
			};

		private:

#ifdef DEBUG
			void test_index(const int n) const
			{
				if (n<0 || n>=_dim) {
					fprintf(stdout,"FATAL ERROR: index out of range in clVector::test_index %d %d\n",n,_dim);
					exit(1);
				}
			};
			void test_indices(const clVector &arg) const
			{
				if (!_vector || !arg._vector) {
					fprintf(stdout,"FATAL ERROR: _vector or arg._vector not allocated in clVector::test_indices\n");
					exit(1);
				}
				if (_dim!=arg._dim) {
					fprintf(stdout,"WARNING: different source and target dimension in clVector::test_indices\n");
					fprintf(stdout,"  use smaller dimension of %d and %d\n",_dim,arg._dim);
				}
			};
#endif

	};

	/* clVector functions */

	/*  Return the inner product of two vectors (only op to the smallest dimension)  */
	template <class T> T contract(const clVector<T> &arg1,const clVector<T> &arg2)
	{
		int dim=(arg1.get_dim() <= arg2.get_dim() ? arg1.get_dim() : arg2.get_dim());
		const T *arg1vectorptr=arg1.get_startptr();
		const T *arg2vectorptr=arg2.get_startptr();
		T sum=0;
		for (int n=0; n<dim; n++) sum+=(*arg1vectorptr++)*(*arg2vectorptr++);
		return sum;
	};

	/*  Return the RMS error between two vectors (must have the same dimension)  */
	template <class T> double rms_error(const clVector<T> &arg1,const clVector<T> &arg2,int dim=0)
	{
		if (dim<0) {
			fprintf(stdout,"FATAL ERROR: Illegal dimension in clVector::rms_error\n");
			exit(1);
		}
		if (!dim) {
			dim=(arg1.get_dim() <= arg2.get_dim() ? arg1.get_dim() : arg2.get_dim());
		} else {
			if (arg1.get_dim()<dim) dim=arg1.get_dim();
			if (arg2.get_dim()<dim) dim=arg2.get_dim();
		}
		const T *arg1vectorptr=arg1.get_startptr();
		const T *arg2vectorptr=arg2.get_startptr();
		T dummy,rms_error=0;
		for (int n=0; n<dim; n++) {
			dummy=(*arg1vectorptr++)-(*arg2vectorptr++);
			rms_error+=dummy*dummy;
		}
		return sqrt(rms_error/dim);
	};

	/*  Broadcast a vector from processor ip  */
	template <class T> void bcast(clVector<T> &vptr,const int ip=0)
	{
#ifdef DIST_MPI
		int num_proc=clMpi_namespace::get_num_proc();
		int proc_id=clMpi_namespace::get_proc_id();
		if (num_proc<2) return;
		int dim=-1;
		T *vector;
		if (proc_id==ip) dim=vptr.get_dim();
		clMpi_namespace::bcast(&dim,1,ip);
		if (!dim) {
			vptr.deallocate();
			return;
		}
		if (proc_id!=ip) vptr.allocate(dim);
		vector=vptr.get_startptr();
		clMpi_namespace::bcast(vector,dim,ip);
#endif
	};

}; // namespace clVector_namespace

#endif // CLVECTOR_H__INCLUDED_
