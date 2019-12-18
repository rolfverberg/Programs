/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/

#ifndef CLMATRIX_H__INCLUDED_
#define CLMATRIX_H__INCLUDED_

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

namespace clVector_namespace
{
	template <class T> class clMatrix
	{
		private:

			int _dim1,_dim2,_dim;	// Dimensions
			T *_startptr;			// Matrix start pointer
			clVector<T> _matrix;	// Matrix stored as a single vector
#ifdef ALLOCDEBUG
			static int _num_matrix;
#endif

		public:

			/*  Constructors  */
			clMatrix(void) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_startptr(NULL)
			{ };
			clMatrix(const int dim1,const int dim2) : 
				_dim1(0),
				_dim2(0),
				_dim(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2);
			};
			clMatrix(const int dim1,const int dim2,const T initval) : 
				_dim1(0),
				_dim2(0),
				_dim(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2,initval);
			};

			/*  Destructor  */
			virtual ~clMatrix(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clMatrix(const clMatrix &mm) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_startptr(NULL)
			{
				_dim1=mm._dim1;
				_dim2=mm._dim2;
				_dim=mm._dim;
				_matrix=mm._matrix;
				_startptr=_matrix.get_startptr();
			};

			/*  Default Assignment Constructor  */
			clMatrix& operator= (const clMatrix &mm)
			{
				if (this==&mm) return *this;

				if (!_startptr) {
					if (_dim1 || _dim2 || _dim) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clMatrix::operator=\n");
						exit(1);
					}
					allocate(mm._dim1,mm._dim2);
				} else {
					if (_dim1!=mm._dim1 || _dim2!=mm._dim2 || _dim!=mm._dim) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clMatrix::operator=\n");
						exit(1);
					}
				}
				_matrix=mm._matrix;
				return *this;
			};

			/*  Operator overloading functions  */

			/* This returns a pointer to the n1-th row
			   It allows you to use value=matrix[n1][n2], however it prevents column bounds checking */
			T* operator[](const int n1)
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::operator[]\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim2*n1);
			};
			T* operator[](const int n1) const
			/* This allows you to use value=matrix[n1][n2], however it prevents column bounds checking */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::operator[]\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim2*n1);
			};

			T* operator()(const int n1)
			/* This returns a pointer to the n1-th row, same as using [] */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::operator()\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim2*n1);
			};
			T* operator()(const int n1) const
			/* This returns a pointer to the n1-th row, same as using [] */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::operator()\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim2*n1);
			};

			T& operator()(const int n1,const int n2)
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				return *(_startptr+n2+_dim2*n1);
			};
			T& operator()(const int n1,const int n2) const
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				return *(_startptr+n2+_dim2*n1);
			};

			clMatrix& operator= (const T arg)
			{
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::operator=\n");
					exit(1);
				}
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)=arg;
				return *this;
			};
			/*  Converts a regular 2D array to a matrix, does not allow bounce checking  */
			clMatrix& operator= (const T* const *mm)
			{
				if (!mm || !mm[0]) {
					fprintf(stdout,"FATAL ERROR: mm not allocated in clMatrix::operator=\n");
					exit(1);
				}
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: _matrix not allocated in clMatrix::operator=\n");
					exit(1);
				}
				T *matrixptr=_startptr;
				const T *mmptr=mm[0];
				for (int n=0; n<_dim; n++) (*matrixptr++)=(*mmptr++);
				return *this;
			};

			clMatrix operator+ (const clMatrix &mm) const
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)+(*mmptr++);
				return sum;
			};
			clMatrix operator+ (const T arg) const
			{
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)+arg;
				return sum;
			};
			friend clMatrix<T> operator+ (const T arg,const clMatrix<T>& mm)
			{
				clMatrix<T> sum(mm._dim1,mm._dim2);
				T *sumptr=sum._startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<mm._dim; n++) (*sumptr++)=arg+(*mmptr++);
				return sum;
			};

			clMatrix operator- (const clMatrix &mm) const
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)-(*mmptr++);
				return sum;
			};
			clMatrix operator- (const T arg) const
			{
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)-arg;
				return sum;
			};
			friend clMatrix<T> operator- (const T arg,const clMatrix<T>& mm)
			{
				clMatrix<T> sum(mm._dim1,mm._dim2);
				T *sumptr=sum._startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<mm._dim; n++) (*sumptr++)=arg-(*mmptr++);
				return sum;
			};

			clMatrix operator* (const clMatrix &mm) const
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)*(*mmptr++);
				return sum;
			};
			clMatrix operator* (const T arg) const
			{
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)*arg;
				return sum;
			};
			friend clMatrix<T> operator* (const T arg,const clMatrix<T>& mm)
			{
				clMatrix<T> sum(mm._dim1,mm._dim2);
				T *sumptr=sum._startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<mm._dim; n++) (*sumptr++)=arg*(*mmptr++);
				return sum;
			};

			clMatrix operator/ (const clMatrix &mm) const
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)/(*mmptr++);
				return sum;
			};
			clMatrix operator/ (const T arg) const
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clMatrix::operator/\n");
					exit(1);
				}
				clMatrix<T> sum(_dim1,_dim2);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)/arg;
				return sum;
			};
			friend clMatrix<T> operator/ (const T arg,const clMatrix<T>& mm)
			{
				clMatrix<T> sum(mm._dim1,mm._dim2);
				T *sumptr=sum._startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<mm._dim; n++) (*sumptr++)=arg/(*mmptr++);
				return sum;
			};

			clMatrix& operator+= (const clMatrix &mm)
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)+=(*mmptr++);
				return *this;
			};
			clMatrix& operator+= (const T arg)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)+=arg;
				return *this;
			};

			clMatrix& operator-= (const clMatrix &mm)
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)-=(*mmptr++);
				return *this;
			};
			clMatrix& operator-= (const T arg)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)-=arg;
				return *this;
			};

			clMatrix& operator*= (const clMatrix &mm)
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)*=(*mmptr++);
				return *this;
			};
			clMatrix& operator*= (const T arg)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)*=arg;
				return *this;
			};

			clMatrix& operator/= (const clMatrix &mm)
			{
#ifdef DEBUG
				test_indices(mm);
#endif
				T *matrixptr=_startptr;
				const T *mmptr=mm._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)/=(*mmptr++);
				return *this;
			};
			clMatrix& operator/= (const T arg)
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clMatrix::operator/=\n");
					exit(1);
				}
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)/=arg;
				return *this;
			};

			clMatrix& operator++ (void)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)++;
				return *this;
			};

			clMatrix& operator-- (void)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)--;
				return *this;
			};

			/*  Remaining member functions  */

			int get_dim1(void) {return _dim1;};
			int get_dim1(void) const {return _dim1;};
			int get_dim2(void) {return _dim2;};
			int get_dim2(void) const {return _dim2;};
			int get_dim(void) {return _dim;};
			int get_dim(void) const {return _dim;};
			T* get_startptr(const int n1=0) {return _startptr+_dim2*n1;};
			const T* get_startptr(const int n1=0) const {return _startptr+_dim2*n1;};
			T* get_startptr(const int n1,const int n2) {return _startptr+n2+_dim2*n1;};
			const T* get_startptr(const int n1,const int n2) const {return _startptr+n2+_dim2*n1;};

			/*  Move data from arg to current */
			void move(clMatrix &arg)
			{
				if (this==&arg) return;

				if (_startptr) deallocate();

				_dim1=arg._dim1;
				_dim2=arg._dim2;
				_dim=arg._dim;
				_startptr=arg._startptr;
				_matrix.move(arg._matrix);

				arg._dim1=0;
				arg._dim2=0;
				arg._dim=0;
				arg._startptr=NULL;
			};

			/*  Swap data content from arg to current, dimensions must match */
			void swap(clMatrix &arg)
			{
				if (this==&arg) return;

				if (_dim1!=arg._dim1 || _dim2!=arg._dim2) {
					fprintf(stdout,"FATAL ERROR: inconsistent dimension in clVector::swap\n");
					exit(1);
				}
				_matrix.swap(arg._matrix);
				_startptr=_matrix.get_startptr();
				arg._startptr=arg._matrix.get_startptr();
			};

			/*  Converts a matrix to a regular 2D array, does not allow bounce checking  */
			void copy(T **mm) const
			{
				if (!mm || !mm[0]) {
					fprintf(stdout,"FATAL ERROR: mm not allocated in clMatrix::copy\n");
					exit(1);
				}
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: _startptr not allocated in clMatrix::copy\n");
					exit(1);
				}
				const T *matrixptr=_startptr;
				T *mmptr=mm[0];
				for (int n=0; n<_dim; n++) (*mmptr++)=(*matrixptr++);
			};

			/*  Allocate matrix of dimension (dim1,dim2)  */
			T* allocate(const int dim1,const int dim2)
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
					_matrix.allocate(_dim);
					_startptr=_matrix.get_startptr();
#ifdef ALLOCDEBUG
					_num_matrix++;
					fprintf(stdout,"matrix %p size %dx%d allocated %d\n",_startptr,_dim1,_dim2,_num_matrix);
#endif
				}
				return _startptr;
			};
			/*  Allocate matrix of dimension (dim1,dim2) and initialize to initval  */
			T* allocate(const int dim1,const int dim2,const T initval)
			{
				allocate(dim1,dim2);
				if (_startptr) {
					T *matrixptr=_startptr;
					for (int n=0; n<_dim; n++) (*matrixptr++)=initval;
				}
				return _startptr;
			};

			/*  Resize a matrix to new row dimensions (new_dim1,_dim2) */
			/*  Data that does not fit is lost  */
			/*  Zero out cells that are not used yet */
			T* redimension(const int new_dim1)
			{
				if (!new_dim1) {
					deallocate();
					return _startptr;
				}
				if (!_startptr) {
					allocate(new_dim1,_dim2);
					return _startptr;
				}
				if (!_dim) {
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clMatrix::redimension\n");
					exit(1);
				}
				const int new_dim=new_dim1*_dim2;
				_matrix.redimension(new_dim);
				_dim1=new_dim1;
				_dim=_dim1*_dim2;
				_startptr=_matrix.get_startptr();
				return _startptr;
			};

			/*  Deallocate matrix  */
			void deallocate(void)
			{
				if (_startptr) {
#ifdef ALLOCDEBUG
					_num_matrix--;
					fprintf(stdout,"matrix %p size %dx%d deallocated %d\n",_startptr,_dim1,_dim2,_num_matrix);
#endif
					_startptr=NULL;
				}
				_matrix.deallocate();
				_dim1=0;
				_dim2=0;
				_dim=0;
			};

			/*  Allocate or set an existing matrix of dimensions (dim1,dim2) to the identity matrix  */
			clMatrix& identity(int dim=0)
			{
				if (!dim) {
					if (!_dim) return *this;
					dim=(_dim1<_dim2 ? _dim1 : _dim2);
					_matrix=(T)(0.0);
				} else {
					if (_dim) {
						fprintf(stdout,"WARNING: matrix already allocated in clMatrix::identity, overwrite existing one\n");
					}
					allocate(dim,dim,(T)0.0);
				}
				const int num_skip=_dim2+1;
				T *matrixptr=_startptr;
				for (int n=0; n<dim; n++,matrixptr+=num_skip) (*matrixptr)=(T)(1.0);
				return *this;
			};

			/*  Allocate matrix of dimensions (dim1,dim2) and initialize to (1,2,3,...,dim1*dim2)  */
			clMatrix& index(const int dim1,const int dim2)
			{
				allocate(dim1,dim2);
				T *matrixptr=_startptr;
				for (int n=1; n<=_dim; n++) (*matrixptr++)=(T)(n);
				return *this;
			};
			/*  Initialize existing matrix to (1,2,3,...,dim1*dim2)  */
			clMatrix& index(void)
			{
				T *matrixptr=_startptr;
				for (int n=1; n<=_dim; n++) (*matrixptr++)=(T)(n);
				return *this;
			};

			/*  Return L2-norm of matrix  */
			T norm(void) const { return (T)_matrix.norm(); };

			/*  Return vector with L2-norm of matrix over rows (dir=1) or columns (dir=2) direction */
			/*    return vector dimension is dim2 and dim1, respectively  */
			/*    (requires copy out)  */
			clVector<T> norm(const int dir) const
			{
				int n1,n2;
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::norm %d\n",dir);
					exit(1);
				}
				clVector<T> vnorm;
				if (dir==1) {
					vnorm.allocate(_dim2,0);
					T *vnormptr=vnorm.get_startptr();
					const T *matrixptr;
					for (n2=0; n2<_dim2; n2++,vnormptr++) {
						matrixptr=_startptr+n2;
						for (n1=0; n1<_dim1; n1++,matrixptr+=_dim2) (*vnormptr)+=(*matrixptr)*(*matrixptr);
					}
				} else if (dir==2) {
					vnorm.allocate(_dim1,0);
					T *vnormptr=vnorm.get_startptr();
					const T *matrixptr=_startptr;
					for (n1=0; n1<_dim1; n1++,vnormptr++) {
						for (n2=0; n2<_dim2; n2++,matrixptr++) (*vnormptr)+=(*matrixptr)*(*matrixptr);
					}
				} else {
					fprintf(stdout,"FATAL ERROR: illegal value of dir clMatrix::norm %d\n",dir);
					exit(1);
				}
				return vnorm;
			};

			/*  Return RMS of matrix  */
			double rms(void) const { return _matrix.rms(); };

			/*  Return vector with RMS of matrix over rows (dir=1) or columns (dir=2) direction */
			/*    return vector dimension is dim2 and dim1, respectively  */
			/*    (requires copy out)  */
			clVector<double> rms(const int dir) const
			{
				int n1,n2;
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::rms %d\n",dir);
					exit(1);
				}
				clVector<double> vrms;
				if (dir==1) {
					vrms.allocate(_dim2,0);
					double *vrmsptr=vrms.get_startptr();
					const T *matrixptr;
					for (n2=0; n2<_dim2; n2++,vrmsptr++) {
						matrixptr=_startptr+n2;
						for (n1=0; n1<_dim1; n1++,matrixptr+=_dim2) (*vrmsptr)+=(*matrixptr)*(*matrixptr);
						(*vrmsptr)=sqrt((*vrmsptr)/(double)(_dim1));
					}
				} else if (dir==2) {
					vrms.allocate(_dim1,0);
					double *vrmsptr=vrms.get_startptr();
					const T *matrixptr=_startptr;
					for (n1=0; n1<_dim1; n1++,vrmsptr++) {
						for (n2=0; n2<_dim2; n2++,matrixptr++) (*vrmsptr)+=(*matrixptr)*(*matrixptr);
						(*vrmsptr)=sqrt((*vrmsptr)/(double)(_dim2));
					}
				} else {
					fprintf(stdout,"FATAL ERROR: illegal value of dir clMatrix::rms %d\n",dir);
					exit(1);
				}
				return vrms;
			};

			/*  Return RMS of the off-diagonal elements of square matrix  */
			double rms_offdiagonal(void) const
			{
				if (!_dim) return 0.0;
				if (_dim1!=_dim2) {
					fprintf(stdout,"FATAL ERROR: non-square matrix in clMatrix::rms_offdiagonal\n");
					exit(1);
				}
				int n1,n2;
				double sumsq=0;
				const T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++) {
					for (n2=0; n2<_dim2; n2++,matrixptr++) {
						if (n1!=n2) sumsq+=(double)((*matrixptr)*(*matrixptr));
					}
				}
				return sqrt(sumsq/((double)(_dim1*(_dim2-1))));
			};

			/*  Return max of matrix  */
			T max(void) const { return _matrix.max(); };
			/*  Return max of matrix in subrange  */
			T max(const int dim1low,const int dim1upp,const int dim2low,const int dim2upp) const
			{
				if (!_dim || dim1low>dim1upp || dim2low>dim2upp) return 0;
#ifdef DEBUG
				test_index1(dim1low);
				test_index1(dim1upp);
				test_index2(dim2low);
				test_index2(dim2upp);
#endif
				int n1,n2;
				const T *matrixptr=_matrix.get_startptr()+dim2low+_dim2*dim1low;
				T mmax=(*matrixptr);
				for (n1=dim1low; n1<=dim1upp; n1++) {
					matrixptr=_matrix.get_startptr()+dim2low+_dim2*n1;
					for (n2=dim2low; n2<=dim2upp; n2++,matrixptr++)	if ((*matrixptr)>mmax) mmax=(*matrixptr);
				}
				return mmax;
			};
			/*  Return global index of max of matrix  */
			int maxindex(void) const { return _matrix.maxindex(); };
			/*  Return directional indices of max of matrix  */
			void maxindex(int &mmaxindex1,int &mmaxindex2) const
			{
				if (!_startptr) {
					mmaxindex1=-1;
					mmaxindex2=-1;
					return;
				}
				const T *matrixptr=_startptr;
				int n1,n2;
				mmaxindex1=0;
				mmaxindex2=0;
				T mmax=(*matrixptr);
				for (n1=0; n1<_dim1; n1++) for (n2=0; n2<_dim2; n2++,matrixptr++) 
					if ((*matrixptr)>mmax) {
						mmax=(*matrixptr);
						mmaxindex1=n1;
						mmaxindex2=n2;
					}
			};
			/*  Return max of matrix and its global index  */
			void max(T &mmax,int &mmaxindex) const
			{
				if (!_startptr) {
					mmax=0;
					mmaxindex=-1;
					return;
				}
				const T *matrixptr=_startptr;
				mmaxindex=0;
				mmax=(*matrixptr++);
				for (int n=1; n<_dim; n++,matrixptr++) 
					if ((*matrixptr)>mmax) {
						mmax=(*matrixptr);
						mmaxindex=n;
					}
			};
			/*  Return max of matrix and its directional indices  */
			void max(T &mmax,int &mmaxindex1,int &mmaxindex2) const
			{
				if (!_startptr) {
					mmax=0;
					mmaxindex1=-1;
					mmaxindex2=-1;
					return;
				}
				const T *matrixptr=_startptr;
				int n1,n2,mmaxindex=0;
				mmax=(*matrixptr);
				mmaxindex1=0;
				mmaxindex2=0;
				for (n1=0; n1<_dim1; n1++) for (n2=0; n2<_dim2; n2++,matrixptr++) 
					if ((*matrixptr)>mmax) {
						mmax=(*matrixptr);
						mmaxindex1=n1;
						mmaxindex2=n2;
					}
			};
			/*  Return vector with max of matrix over rows (dir=1) or columns (dir=2) direction */
			/*    return vector dimension is dim2 and dim1, respectively  */
			/*    (requires copy out)  */
			clVector<T> max(const int dir) const
			{
				int n1,n2;
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::max %d\n",dir);
					exit(1);
				}
				clVector<T> vmax;
				if (dir==1) {
					vmax.allocate(_dim2);
					T *vmaxptr=vmax.get_startptr();
					const T *matrixptr;
					for (n2=0; n2<_dim2; n2++,vmaxptr++) {
						matrixptr=_startptr+n2;
						(*vmaxptr)=(*matrixptr);
						matrixptr+=_dim2;
						for (n1=1; n1<_dim1; n1++,matrixptr+=_dim2) 
							if ((*matrixptr)>(*vmaxptr)) (*vmaxptr)=(*matrixptr);
					}
				} else if (dir==2) {
					vmax.allocate(_dim1);
					T *vmaxptr=vmax.get_startptr();
					const T *matrixptr=_startptr;
					for (n1=0; n1<_dim1; n1++,vmaxptr++) {
						(*vmaxptr)=(*matrixptr++);
						for (n2=1; n2<_dim2; n2++,matrixptr++) 
							if ((*matrixptr)>(*vmaxptr)) (*vmaxptr)=(*matrixptr);
					}
				} else {
					fprintf(stdout,"FATAL ERROR: illegal value of dir clMatrix::max %d\n",dir);
					exit(1);
				}
				return vmax;
			};

			/*  Return min of matrix  */
			T min(void) const { return _matrix.min(); };
			/*  Return min of matrix in subrange  */
			T min(const int dim1low,const int dim1upp,const int dim2low,const int dim2upp) const
			{
				if (!_dim || dim1low>dim1upp || dim2low>dim2upp) return 0;
#ifdef DEBUG
				test_index1(dim1low);
				test_index1(dim1upp);
				test_index2(dim2low);
				test_index2(dim2upp);
#endif
				int n1,n2;
				const T *matrixptr=_matrix.get_startptr()+dim2low+_dim2*dim1low;
				T mmin=(*matrixptr);
				for (n1=dim1low; n1<=dim1upp; n1++) {
					matrixptr=_matrix.get_startptr()+dim2low+_dim2*n1;
					for (n2=dim2low; n2<=dim2upp; n2++,matrixptr++)	if ((*matrixptr)<mmin) mmin=(*matrixptr);
				}
				return mmin;
			};
			/*  Return global index of min of matrix  */
			int minindex(void) const { return _matrix.minindex(); };
			/*  Return directional indices of min of matrix  */
			void minindex(int &mminindex1,int &mminindex2) const
			{
				if (!_startptr) {
					mminindex1=-1;
					mminindex2=-1;
					return;
				}
				const T *matrixptr=_startptr;
				int n1,n2;
				mminindex1=0;
				mminindex2=0;
				T mmin=(*matrixptr);
				for (n1=0; n1<_dim1; n1++) for (n2=0; n2<_dim2; n2++,matrixptr++) 
					if ((*matrixptr)<mmin) {
						mmin=(*matrixptr);
						mminindex1=n1;
						mminindex2=n2;
					}
			};
			/*  Return min of matrix and its global index  */
			void min(T &mmin,int &mminindex) const
			{
				if (!_startptr) {
					mmin=0;
					mminindex=-1;
					return;
				}
				const T *matrixptr=_startptr;
				mminindex=0;
				mmin=(*matrixptr++);
				for (int n=1; n<_dim; n++,matrixptr++) 
					if ((*matrixptr)<mmin) {
						mmin=(*matrixptr);
						mminindex=n;
					}
			};
			/*  Return min of matrix and its directional indices  */
			void min(T &mmin,int &mminindex1,int &mminindex2) const
			{
				if (!_startptr) {
					mmin=0;
					mminindex1=-1;
					mminindex2=-1;
					return;
				}
				const T *matrixptr=_startptr;
				int n1,n2,mminindex=0;
				mmin=(*matrixptr);
				mminindex1=0;
				mminindex2=0;
				for (n1=0; n1<_dim1; n1++) for (n2=0; n2<_dim2; n2++,matrixptr++) 
					if ((*matrixptr)<mmin) {
						mmin=(*matrixptr);
						mminindex1=n1;
						mminindex2=n2;
					}
			};
			/*  Return vector with min of matrix over rows (dir=1) or columns (dir=2) */
			/*    return vector dimension is dim2 and dim1, respectively  */
			/*    (requires copy out)  */
			clVector<T> min(const int dir) const
			{
				int n1,n2;
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::min %d\n",dir);
					exit(1);
				}
				clVector<T> vmin;
				if (dir==1) {
					vmin.allocate(_dim2);
					T *vminptr=vmin.get_startptr();
					const T *matrixptr;
					for (n2=0; n2<_dim2; n2++,vminptr++) {
						matrixptr=_startptr+n2;
						(*vminptr)=(*matrixptr);
						matrixptr+=_dim2;
						for (n1=1; n1<_dim1; n1++,matrixptr+=_dim2) 
							if ((*matrixptr)<(*vminptr)) (*vminptr)=(*matrixptr);
					}
				} else if (dir==2) {
					vmin.allocate(_dim1);
					T *vminptr=vmin.get_startptr();
					const T *matrixptr=_startptr;
					for (n1=0; n1<_dim1; n1++,vminptr++) {
						(*vminptr)=(*matrixptr++);
						for (n2=1; n2<_dim2; n2++,matrixptr++) 
							if ((*matrixptr)<(*vminptr)) (*vminptr)=(*matrixptr);
					}
				} else {
					fprintf(stdout,"FATAL ERROR: illegal value of dir clMatrix::min %d\n",dir);
					exit(1);
				}
				return vmin;
			};

			/*  Return mean of matrix  */
			double mean(void) const { return _matrix.mean(); };
			/*  Return vector with mean of matrix over rows (dir=1) or columns (dir=2) */
			/*    return vector dimension is dim2 and dim1, respectively  */
			/*    (requires copy out)  */
			clVector<double> mean(const int dir) const
			{
				int n1,n2;
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::mean %d\n",dir);
					exit(1);
				}
				clVector<double> vmean;
				if (dir==1) {
					vmean.allocate(_dim2,0);
					double *vmeanptr=vmean.get_startptr();
					const T *matrixptr;
					for (n2=0; n2<_dim2; n2++,vmeanptr++) {
						matrixptr=_startptr+n2;
						for (n1=0; n1<_dim1; n1++,matrixptr+=_dim2) (*vmeanptr)+=(*matrixptr);
						(*vmeanptr)/=_dim1;
					}
				} else if (dir==2) {
					vmean.allocate(_dim1,0);
					double *vmeanptr=vmean.get_startptr();
					const T *matrixptr=_startptr;
					for (n1=0; n1<_dim1; n1++,vmeanptr++) {
						for (n2=0; n2<_dim2; n2++) (*vmeanptr)+=(*matrixptr++);
						(*vmeanptr)/=_dim2;
					}
				} else {
					fprintf(stdout,"FATAL ERROR: illegal value of dir clMatrix::mean %d\n",dir);
					exit(1);
				}
				return vmean;
			};

			/*  Return stdv of matrix  */
			double stdv(void) const { return _matrix.stdv(); };
			/*  Return vector with stdv of matrix over rows (dir=1) or columns (dir=2) */
			/*    return vector dimension is dim2 and dim1, respectively  */
			/*    (requires copy out)  */
			clVector<double> stdv(const int dir) const
			{
				int n1,n2;
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix::stdv %d\n",dir);
					exit(1);
				}
				double dummy,*vmeanptr=NULL;
				clVector<double> vmean=mean(dir);
				vmeanptr=vmean.get_startptr();
				clVector<double> vstdv;
				if (dir==1) {
					vstdv.allocate(_dim2,0.0);
					if (_dim1<2) return vstdv;
					double *vstdvptr=vstdv.get_startptr();
					const T *matrixptr;
					for (n2=0; n2<_dim2; n2++,vstdvptr++,vmeanptr++) {
						matrixptr=_startptr+n2;
						for (n1=0; n1<_dim1; n1++,matrixptr+=_dim2) {
							dummy=(*matrixptr)-(*vmeanptr);
							(*vstdvptr)+=dummy*dummy;
						}
						(*vstdvptr)=sqrt((*vstdvptr)/_dim1);
					}
				} else if (dir==2) {
					vstdv.allocate(_dim1,0.0);
					if (_dim2<2) return vstdv;
					double *vstdvptr=vstdv.get_startptr();
					const T *matrixptr=_startptr;
					for (n1=0; n1<_dim1; n1++,vstdvptr++,vmeanptr++) {
						for (n2=0; n2<_dim2; n2++) {
							dummy=(*matrixptr++)-(*vmeanptr);
							(*vstdvptr)+=dummy*dummy;
						}
						(*vstdvptr)=sqrt((*vstdvptr)/_dim2);
					}
				} else {
					fprintf(stdout,"FATAL ERROR: illegal value of dir clMatrix::stdv %d\n",dir);
					exit(1);
				}
				return vstdv;
			};

			/* set matrix to the row correlation coefficient */
			clMatrix& rowcorr(const clMatrix<T> &mm)
			{
				if (!_startptr) allocate(mm.get_dim1(),mm.get_dim1());
				if (_dim1!=_dim2 || _dim1!=mm.get_dim1()) {
					fprintf(stdout,"FATAL ERROR: illegal matrix dimensions in clMatrix::rowcorr\n");
					exit(1);
				}
				int n1,n2,n_sum;
				const T *mm1ptr=NULL,*mm2ptr=NULL;
				const double *mean1ptr=NULL,*mean2ptr=NULL;
				const int dim_sum=mm.get_dim2();
				clVector<double> mean=mm.mean(2);
				mean1ptr=mean.get_startptr();
				_matrix=(T)0;
				T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++,mean1ptr++) {
					for (n2=0; n2<n1; n2++) (*matrixptr++)=*(_startptr+n1+_dim2*n2);
					mean2ptr=mean.get_startptr()+n1;
					for (n2=n1; n2<_dim2; n2++,mean2ptr++,matrixptr++) {
						mm1ptr=mm(n1);
						mm2ptr=mm(n2);
						for (n_sum=0; n_sum<dim_sum; n_sum++) 
							(*matrixptr)+=((*mm1ptr++)-(T)(*mean1ptr))*((*mm2ptr++)-(T)(*mean2ptr));
					}
				}
				const int skip=_dim1+1;
				clVector<double> norm(_dim1);
				double *norm1ptr=norm.get_startptr();
				matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++,matrixptr+=skip) {
					if (!(*matrixptr)) (*norm1ptr++)=(T)1;
					else (*norm1ptr++)=sqrt(*matrixptr);
				}
				norm1ptr=norm.get_startptr();
				double *norm2ptr=NULL;
				matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++,norm1ptr++) {
					norm2ptr=norm.get_startptr();
					for (n2=0; n2<_dim2; n2++) (*matrixptr++)/=(T)((*norm1ptr)*(*norm2ptr++));
				}
				return *this;
			};

			/* set matrix to the column correlation coefficient */
			clMatrix& colcorr(const clMatrix<T> &mm)
			{
				if (!_startptr) allocate(mm.get_dim2(),mm.get_dim2());
				else if (_dim1!=_dim2 || _dim1!=mm.get_dim2()) {
					fprintf(stdout,"FATAL ERROR: illegal matrix dimensions in clMatrix::colcorr\n");
					exit(1);
				}
				int n1,n2,n_sum;
				const T *mm1ptr=NULL,*mm2ptr=NULL;
				const double *mean1ptr=NULL,*mean2ptr=NULL;
				const int dim_sum=mm.get_dim1();
				clVector<double> mean=mm.mean(1);
				mean1ptr=mean.get_startptr();
				_matrix=(T)0;
				T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++,mean1ptr++) {
					for (n2=0; n2<n1; n2++) (*matrixptr++)=*(_startptr+n1+_dim2*n2);
					mean2ptr=mean.get_startptr()+n1;
					for (n2=n1; n2<_dim2; n2++,mean2ptr++,matrixptr++) {
						mm1ptr=mm.get_startptr()+n1;
						mm2ptr=mm.get_startptr()+n2;
						for (n_sum=0; n_sum<dim_sum; n_sum++,mm1ptr+=_dim1,mm2ptr+=_dim2) 
							(*matrixptr)+=((*mm1ptr)-(T)(*mean1ptr))*((*mm2ptr)-(T)(*mean2ptr));
					}
				}
				const int skip=_dim1+1;
				clVector<double> norm(_dim1);
				double *norm1ptr=norm.get_startptr();
				matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++,matrixptr+=skip) {
					if (!(*matrixptr)) (*norm1ptr++)=(T)1;
					else (*norm1ptr++)=sqrt(*matrixptr);
				}
				norm1ptr=norm.get_startptr();
				double *norm2ptr=NULL;
				matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++,norm1ptr++) {
					norm2ptr=norm.get_startptr();
					for (n2=0; n2<_dim2; n2++) (*matrixptr++)/=(T)((*norm1ptr)*(*norm2ptr++));
				}
				return *this;
			};

			/* set matrix to aa*bb */
			clMatrix& mmult_ab(const clMatrix<T> &aa,const clMatrix<T> &bb)
			{
				if (!_startptr) allocate(aa.get_dim1(),bb.get_dim2());
				else if (_dim1!=aa.get_dim1() || _dim2!=bb.get_dim2() || aa.get_dim2()!=bb.get_dim1()) {
					fprintf(stdout,"FATAL ERROR: illegal matrix dimensions in clMatrix::mmult_ab\n");
					exit(1);
				}
				int n1,n2,n_sum;
				const T *aaptr=NULL,*bbptr=NULL;
				const int dim_sum=aa.get_dim2();
				_matrix=(T)0;
				T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++) {
					for (n2=0; n2<_dim2; n2++,matrixptr++) {
						aaptr=aa.get_startptr()+n1;
						bbptr=bb.get_startptr()+n2;
						for (n_sum=0; n_sum<dim_sum; n_sum++,bbptr+=_dim2) (*matrixptr)+=(*aaptr++)*(*bbptr);
					}
				}
				return *this;
			};

			/* set matrix to mm*mm_transpose */
			clMatrix& mmult_aat(const clMatrix<T> &mm)
			{
				if (!_startptr) allocate(mm.get_dim1(),mm.get_dim1());
				if (_dim1!=_dim2 || _dim1!=mm.get_dim1()) {
					fprintf(stdout,"FATAL ERROR: illegal matrix dimensions in clMatrix::mmult_aat\n");
					exit(1);
				}
				int n1,n2,n_sum;
				const T *mm1ptr=NULL,*mm2ptr=NULL;
				const int dim_sum=mm.get_dim2();
				_matrix=(T)0;
				T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++) {
					for (n2=0; n2<n1; n2++) (*matrixptr++)=*(_startptr+n1+_dim2*n2);
					for (n2=n1; n2<_dim2; n2++,matrixptr++) {
						mm1ptr=mm(n1);
						mm2ptr=mm(n2);
						for (n_sum=0; n_sum<dim_sum; n_sum++) {
							(*matrixptr)+=(*mm1ptr++)*(*mm2ptr++);
						}
					}
				}
				return *this;
			};

			/* set matrix to aa*bb_transpose */
			clMatrix& mmult_abt(const clMatrix<T> &aa,const clMatrix<T> &bb)
			{
				if (!_startptr) allocate(aa.get_dim1(),bb.get_dim1());
				else if (_dim1!=aa.get_dim1() || _dim2!=bb.get_dim1() || aa.get_dim2()!=bb.get_dim2()) {
					fprintf(stdout,"FATAL ERROR: illegal matrix dimensions in clMatrix::mmult_abt\n");
					exit(1);
				}
				int n1,n2,n_sum;
				const T *aaptr=NULL,*bbptr=NULL;
				const int dim_sum=aa.get_dim2();
				_matrix=(T)0;
				T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++) {
					for (n2=0; n2<_dim2; n2++,matrixptr++) {
						aaptr=aa(n1);
						bbptr=bb(n2);
						for (n_sum=0; n_sum<dim_sum; n_sum++) (*matrixptr)+=(*aaptr++)*(*bbptr++);
					}
				}
				return *this;
			};

			/* set matrix to mm_transpose*mm */
			clMatrix& mmult_ata(const clMatrix<T> &mm)
			{
				if (!_startptr) allocate(mm.get_dim2(),mm.get_dim2());
				else if (_dim1!=_dim2 || _dim1!=mm.get_dim2()) {
					fprintf(stdout,"FATAL ERROR: illegal matrix dimensions in clMatrix::mmult_ata\n");
					exit(1);
				}
				int n1,n2,n_sum;
				const T *mm1ptr=NULL,*mm2ptr=NULL;
				const int dim_sum=mm.get_dim1();
				_matrix=(T)0;
				T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++) {
					for (n2=0; n2<n1; n2++) (*matrixptr++)=*(_startptr+n1+_dim2*n2);
					for (n2=n1; n2<_dim2; n2++,matrixptr++) {
						mm1ptr=mm.get_startptr()+n1;
						mm2ptr=mm.get_startptr()+n2;
						for (n_sum=0; n_sum<dim_sum; n_sum++,mm1ptr+=_dim1,mm2ptr+=_dim2) (*matrixptr)+=(*mm1ptr)*(*mm2ptr);
					}
				}
				return *this;
			}

			/* set matrix to aa_transpose*bb */
			clMatrix& mmult_atb(const clMatrix<T> &aa,const clMatrix<T> &bb)
			{
				if (!_startptr) allocate(aa.get_dim2(),bb.get_dim2());
				else if (_dim1!=aa.get_dim2() || _dim2!=bb.get_dim2() || aa.get_dim1()!=bb.get_dim1()) {
					fprintf(stdout,"FATAL ERROR: illegal matrix dimensions in clMatrix::mmult_atb\n");
					exit(1);
				}
				int n1,n2,n_sum;
				const T *aaptr=NULL,*bbptr=NULL;
				const int dim_sum=aa.get_dim1();
				_matrix=(T)0;
				T *matrixptr=_startptr;
				for (n1=0; n1<_dim1; n1++) {
					for (n2=0; n2<n1; n2++) (*matrixptr++)=*(_startptr+n1+_dim2*n2);
					for (n2=n1; n2<_dim2; n2++,matrixptr++) {
						aaptr=aa.get_startptr()+n1;
						bbptr=bb.get_startptr()+n2;
						for (n_sum=0; n_sum<dim_sum; n_sum++,aaptr+=_dim1,bbptr+=_dim2) (*matrixptr)+=(*aaptr)*(*bbptr);
					}
				}
				return *this;
			};

			/*  Read a matrix from an ascii file  */
			void read_ascii(const char *filename)
			{
				FILE *fileptr=fopen(filename,"r");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clMatrix::read_ascii\n");
					exit(1);
				}
				read_ascii(fileptr);
				fclose(fileptr);
			};
			/*  Read a matrix from an ascii stream  */
			void read_ascii(FILE *fileptr)
			{
				int dim1,dim2;
				if (fscanf(fileptr,"%d %d\n",&dim1,&dim2)!=2) {
					fprintf(stdout,"FATAL ERROR: error reading matrix dimensions in clMatrix::read_ascii\n");
					exit(1);
				}
				if (dim1 && dim2) {
					allocate(dim1,dim2);
					double readvalue;
					T *matrixptr=_startptr;
					for (int n=0; n<_dim; n++) {
						if (fscanf(fileptr,"%le",&readvalue)!=1) {
							fprintf(stdout,"FATAL ERROR: error reading matrix in clMatrix::read_ascii\n");
							exit(1);
						}
						(*matrixptr++)=(T)readvalue;
					}
				} else deallocate();
			};

			/*  Write the transpose of a matrix to an ascii file  */
			void write_ascii_transpose(const char *filename) const
			{
				FILE *fileptr=fopen(filename,"w");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clMatrix::write_ascii_transpose\n");
					exit(1);
				}
				write_ascii_transpose(fileptr);
				fclose(fileptr);
			};
			/*  Write the transpose of a matrix to an ascii stream  */
			void write_ascii_transpose(FILE *fileptr) const
			{
				int n1,n2;
				fprintf(fileptr,"%d %d\n",_dim1,_dim2);
				for (n2=0; n2<_dim2; n2++) {
					for (n1=0; n1<_dim1; n1++) fprintf(fileptr," %12.4le",(double)(*this)(n1,n2));
					fprintf(fileptr,"\n");
				}
			};

			/*  Write a vector to an ascii file  */
			/*    full_precision!=0: Write in %16.8le format (default)  */
			/*    full_precision=0:  write in full precision (default)  */
			void write_ascii(const char *filename,const int full_precision=0) const
			{
				FILE *fileptr=fopen(filename,"w");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clMatrix::write_ascii\n");
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
				int n1,n2;
				T *matrixptr=_startptr;
				fprintf(fileptr,"%d %d\n",_dim1,_dim2);
				for (n1=0; n1<_dim1; n1++) {
					if (full_precision) {
						for (n2=0; n2<_dim2; n2++) fprintf(fileptr," %24.16le",(double)(*matrixptr++));
					} else {
						for (n2=0; n2<_dim2; n2++) fprintf(fileptr," %10.3le",(double)(*matrixptr++));
					}
					fprintf(fileptr,"\n");
				}
			};

			/*  Read a matrix from a binary file  */
			void read(const char *filename)
			{
				FILE *fileptr=fopen(filename,"rb");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clMatrix::read\n");
					exit(1);
				}
				read(fileptr);
				fclose(fileptr);
			};
			/*  Read a matrix from a binary stream  */
			int read(FILE *fileptr)
			{
				int dim1,dim2;
				void *io_ptr;
				io_ptr=(void *)&dim1;
				if (fread(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
				io_ptr=(void *)&dim2;
				if (fread(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
				if (dim1 && dim2) {
					allocate(dim1,dim2);
					io_ptr=(void *)_matrix.get_startptr();
					if (fread(io_ptr,sizeof(T),_dim,fileptr)!=(unsigned)_dim) return 0;
				}
				return 1;
			};

			/*  Write a vector to a binary file  */
			void write(const char *filename) const
			{
				FILE *fileptr=fopen(filename,"wb");
				if (!fileptr) {
					fprintf(stdout,"FATAL ERROR: unable to open file in clMatrix::write\n");
					exit(1);
				}
				write(fileptr);
				fclose(fileptr);
			};
			/*  Write a vector to a binary stream  */
			int write(FILE *fileptr) const
			{
				void *io_ptr;
				io_ptr=(void *)&_dim1;
				if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
				io_ptr=(void *)&_dim2;
				if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
				if (_dim) {
					io_ptr=(void *)_matrix.get_startptr();
					if (fwrite(io_ptr,sizeof(T),_dim,fileptr)!=(unsigned)_dim) return 0;
				}
				return 1;
			};

		private:

#ifdef DEBUG
			void test_index1(const int n1) const
			{
				if (n1<0 || n1>=_dim1) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix::test_index1: %d %d\n",n1,_dim1);
					exit(1);
				}
			};
			void test_index2(const int n2) const
			{
				if (n2<0 || n2>=_dim2) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix::test_index2: %d %d\n",n2,_dim2);
					exit(1);
				}
			};
			void test_indices(const clMatrix &mm) const
			{
				if (!_startptr || !mm._startptr) {
					fprintf(stdout,"FATAL ERROR: matrix or mm.matrix not allocated in clMatrix::test_indices\n");
					exit(1);
				}
				if (_dim1!=mm._dim1 || _dim2!=mm._dim2) {
					fprintf(stdout,"FATAL ERROR: different source and target dimension in clMatrix::test_indices\n");
				}
			};
#endif

	};

	/* clMatrix functions */

	/* Contract the last index of the input matrix with the input vector and return the result vector */
	/*   This routine needs a copy operation upon leaving, which may be less efficient for large arrays */
	template <class T> clVector<T> contract(const clMatrix<T> &arg1,const clVector<T> &arg2)
	{
		const int dim1=arg1.get_dim1();
		const int dim2=arg1.get_dim2();
		if (dim2!=arg2.get_dim()) {
			fprintf(stdout,"FATAL ERROR: inconsistent dimensions in clVector_namespace::contract\n");
			fprintf(stdout,"%d %d\n",dim2,arg2.get_dim());
			exit(1);
		}
		clVector<T> sum(dim1,0.0);
		T *sumvectorptr=sum.get_startptr();
		const T *arg1matrixptr=arg1.get_startptr();
		const T *arg2vectorptr=NULL;
		int n1,n2;
		for (n1=0; n1<dim1; n1++,sumvectorptr++) {
			arg2vectorptr=arg2.get_startptr();
			for (n2=0; n2<dim2; n2++) {
				(*sumvectorptr)+=(*arg1matrixptr++)*(*arg2vectorptr++);
			}
		}
		return sum;
	};

	/* Contract the last index of the input matrix with the input vector and return the result vector */
	/*   This routine avoids the copy operation upon leaving */
	template <class T> void contract(const clMatrix<T> &arg1,const clVector<T> &arg2,clVector<T> &sum)
	{
		const int dim1=arg1.get_dim1();
		const int dim2=arg1.get_dim2();
		if (dim1!=sum.get_dim()) {
			fprintf(stdout,"FATAL ERROR: inconsistent row dimensions in clVector_namespace::contract\n");
			fprintf(stdout,"%d %d\n",dim1,sum.get_dim());
			exit(1);
		}
		if (dim2!=arg2.get_dim()) {
			fprintf(stdout,"FATAL ERROR: inconsistent column dimensions in clVector_namespace::contract\n");
			fprintf(stdout,"%d %d\n",dim2,arg2.get_dim());
			exit(1);
		}
		sum=0.0;
		T *sumvectorptr=sum.get_startptr();
		const T *arg1matrixptr=arg1.get_startptr();
		const T *arg2vectorptr=NULL;
		int n1,n2;
		for (n1=0; n1<dim1; n1++,sumvectorptr++) {
			arg2vectorptr=arg2.get_startptr();
			for (n2=0; n2<dim2; n2++) {
				(*sumvectorptr)+=(*arg1matrixptr++)*(*arg2vectorptr++);
			}
		}
	};

	/* Contract the last index of the two matrices and store the result in vector */
/*	template <class T> void contract_lastindex(clVector<T> &vector,const clMatrix<T> &arg1,const clMatrix<T> &arg2)
	{
		const int dim1=vector.get_dim();
		const int dim2=arg1.get_dim2();
		if (dim1!=arg1.get_dim1() || dim1!=arg2.get_dim1() || dim2!=arg2.get_dim2()) {
			fprintf(stdout,"FATAL ERROR: inconsistent dimensions in clVector_namespace::contract_columns\n");
			exit(1);
		}
		int n1,n2;
		T *vectorptr=vector.get_startptr();		
		const T *arg1vectorptr=arg1.get_startptr();
		const T *arg2vectorptr=arg2.get_startptr();
		vector=0.0;
		for (n1=0; n1<dim1; n1++,vectorptr++) {
			for (n2=0; n2<dim2; n2++) (*vectorptr)+=(*arg1vectorptr++)*(*arg2vectorptr++);
		}
	};*/

	/* take the inverse of a matrix up to dimension 3 */
	/*    return 0 for zero determinant or illegal dimensions  */
	/*    return 1 for successful completion  */
	template <class T> int inverse(const clMatrix<T> &mm,clMatrix<T> &mminv)
	{
		const int dim=mm.get_dim1();
		if (dim!=mm.get_dim2() || dim!=mminv.get_dim1() || dim!=mminv.get_dim2()) return 0;
		double det;
		if (dim==1) {
			det=mm(0,0);
			if (!det) return 0;
			mminv(0,0)=1.0/det;
		} else if (dim==2) {
			det=mm(0,0)*mm(1,1)-mm(0,1)*mm(1,0);
			if (!det) return 0;
			det=1.0/det;
			mminv(0,0)=mm(1,1)*det;
			mminv(0,1)=-mm(0,1)*det;
			mminv(1,0)=-mm(1,0)*det;
			mminv(1,1)=mm(0,0)*det;
		} else if (dim==3) {
			double det=mm(0,0)*(mm(1,1)*mm(2,2)-mm(2,1)*mm(1,2))-
				mm(0,1)*(mm(1,0)*mm(2,2)-mm(1,2)*mm(2,0))+
				mm(0,2)*(mm(1,0)*mm(2,1)-mm(1,1)*mm(2,0));
			if (!det) return 0;
			det=1.0/det;
			mminv(0,0)=(mm(1,1)*mm(2,2)-mm(2,1)*mm(1,2))*det;
			mminv(0,1)=(mm(0,2)*mm(2,1)-mm(0,1)*mm(2,2))*det;
			mminv(0,2)=(mm(0,1)*mm(1,2)-mm(0,2)*mm(1,1))*det;
			mminv(1,0)=(mm(1,2)*mm(2,0)-mm(1,0)*mm(2,2))*det;
			mminv(1,1)=(mm(0,0)*mm(2,2)-mm(0,2)*mm(2,0))*det;
			mminv(1,2)=(mm(1,0)*mm(0,2)-mm(0,0)*mm(1,2))*det;
			mminv(2,0)=(mm(1,0)*mm(2,1)-mm(2,0)*mm(1,1))*det;
			mminv(2,1)=(mm(2,0)*mm(0,1)-mm(0,0)*mm(2,1))*det;
			mminv(2,2)=(mm(0,0)*mm(1,1)-mm(1,0)*mm(0,1))*det;
		} else return 0;
		return 1;
	};

	/*  Broadcast a matrix from processor ip  */
	template <class T> void bcast(clMatrix<T> &mptr,const int ip=0)
	{
#ifdef DIST_MPI
		int num_proc=clMpi_namespace::get_num_proc();
		int proc_id=clMpi_namespace::get_proc_id();
		if (num_proc<2) return;
		int dim1,dim2;
		T *matrix;
		if (proc_id==ip) {
			dim1=mptr.get_dim1();
			dim2=mptr.get_dim2();
		}
		clMpi_namespace::bcast(&dim1,1,ip);
		clMpi_namespace::bcast(&dim2,1,ip);
		if (!dim1 && !dim2) {
			mptr.deallocate();
			return;
		}
		if (proc_id!=ip) mptr.allocate(dim1,dim2);
		matrix=mptr.get_startptr();
		clMpi_namespace::bcast(matrix,dim1*dim2,ip);
#endif
	};

}; // namespace clVector_namespace

#endif // CLMATRIX_H__INCLUDED_
