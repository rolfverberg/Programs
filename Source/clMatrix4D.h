/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/

#ifndef CLMATRIX4D_H__INCLUDED_
#define CLMATRIX4D_H__INCLUDED_

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
	template <class T> class clMatrix4D
	{
		private:

			int _dim1,_dim2,_dim3,_dim4,_dim;	// Dimensions
			T *_startptr;				// Matrix start pointer
			clVector<T> _matrix;		// Matrix stored as a single vector
#ifdef ALLOCDEBUG
			static int _num_matrix;
#endif

		public:

			/*  Constructors  */
			clMatrix4D(void) :
				_dim1(0),
				_dim2(0),
				_dim3(0),
				_dim4(0),
				_dim(0),
				_startptr(NULL)
			{ };
			clMatrix4D(const int dim1,const int dim2,const int dim3,const int dim4) : 
				_dim1(0),
				_dim2(0),
				_dim3(0),
				_dim4(0),
				_dim(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2,dim3,dim4);
			};
			clMatrix4D(const int dim1,const int dim2,const int dim3,const int dim4,const T initval) : 
				_dim1(0),
				_dim2(0),
				_dim3(0),
				_dim4(0),
				_dim(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2,dim3,dim4,initval);
			};

			/*  Destructor  */
			virtual ~clMatrix4D(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clMatrix4D(const clMatrix4D &arg) :
				_dim1(0),
				_dim2(0),
				_dim3(0),
				_dim4(0),
				_dim(0),
				_startptr(NULL)
			{
				_dim1=arg._dim1;
				_dim2=arg._dim2;
				_dim3=arg._dim3;
				_dim4=arg._dim4;
				_dim=arg._dim;
				_matrix=arg._matrix;
				_startptr=_matrix.get_startptr();
			};

			/*  Default Assignment Constructor  */
			clMatrix4D& operator= (const clMatrix4D &arg)
			{
				if (this==&arg) return *this;

				if (!_startptr) {
					if (_dim1 || _dim2 || _dim3 || _dim4 || _dim) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clMatrix4D::operator=\n");
						exit(1);
					}
					allocate(arg._dim1,arg._dim2,arg._dim3,arg._dim4);
				} else {
					if (_dim1!=arg._dim1 || _dim2!=arg._dim2 || _dim3!=arg._dim3 || _dim4!=arg._dim4 || _dim!=arg._dim) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clMatrix4D::operator=\n");
						exit(1);
					}
				}
				_matrix=arg._matrix;
				return *this;
			};

			/*  Operator overloading functions  */

			T* operator()(const int n1)
			/* This returns a pointer to the n1-th n2/n3/n4-cube */
			/* However it prevents bounds checking in the n2, n3 or n4 direction */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator()\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim4*_dim3*_dim2*n1);
			};
			const T* operator()(const int n1) const
			/* This returns a pointer to the n1-th n2/n3/n4-cube */
			/* However it prevents bounds checking in the n2, n3 or n4 direction */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator()\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim4*_dim3*_dim2*n1);
			};

			T* operator()(const int n1,const int n2)
			/* This returns a pointer to the (n1,n2)-th plane */
			/* However it prevents bounds checking in the n3 or n4 direction */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				return (_startptr+_dim4*_dim3*(n2+_dim2*n1));
			};
			const T* operator()(const int n1,const int n2) const 
			/* This returns a pointer to the (n1,n2)-th plane */
			/* However it prevents bounds checking in the n3 or n4 direction */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				return (_startptr+_dim4*_dim3*(n2+_dim2*n1));
			};

			T* operator()(const int n1,const int n2,const int n3)
			/* This returns a pointer to the (n1,n2,n3)-th row */
			/* However it prevents bounds checking in the n4 direction */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				return (_startptr+_dim4*(n3+_dim3*(n2+_dim2*n1)));
			};
			const T* operator()(const int n1,const int n2,const int n3) const 
			/* This returns a pointer to the (n1,n2,n3)-th row */
			/* However it prevents bounds checking in the n4 direction */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator()\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
				test_index2(n3);
#endif
				return (_startptr+_dim4*(n3+_dim3*(n2+_dim2*n1)));
			};

			T& operator()(const int n1,const int n2,const int n3,const int n4)
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator[]\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
				test_index3(n3);
				test_index3(n4);
#endif
				return *(_startptr+n4+_dim4*(n3+_dim3*(n2+_dim2*n1)));
			};
			T& operator()(const int n1,const int n2,const int n3,const int n4) const
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator[]\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
				test_index3(n3);
				test_index3(n4);
#endif
				return *(_startptr+n4+_dim4*(n3+_dim3*(n2+_dim2*n1)));
			};

			clMatrix4D& operator= (const T arg)
			{
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clMatrix4D::operator=\n");
					exit(1);
				}
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)=arg;
				return *this;
			};
			/*  Converts a regular 4D array to a matrix, does not allow bounce checking  */
			clMatrix4D& operator= (const T* const* const* const *arg)
			{
				if (!arg || !arg[0] || !arg[0][0] || !arg[0][0][0]) {
					fprintf(stdout,"FATAL ERROR: arg not allocated in clMatrix4D::operator=\n");
					exit(1);
				}
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: _matrix not allocated in clMatrix4D::operator=\n");
					exit(1);
				}
				T *matrixptr=_startptr;
				const T *argptr=arg[0][0][0];
				for (int n=0; n<_dim; n++) (*matrixptr++)=(*argptr++);
				return *this;
			};

			clMatrix4D operator+ (const clMatrix4D &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)+(*argmatrixptr++);
				return sum;
			};
			clMatrix4D operator+ (const T arg) const
			{
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)+arg;
				return sum;
			};
			friend clMatrix4D<T> operator+ (const T arg1,const clMatrix4D<T>& arg2)
			{
				clMatrix4D<T> sum(arg2._dim1,arg2._dim2,arg2._dim3,arg2._dim4);
				T *sumptr=sum._startptr;
				const T *arg2ptr=arg2._startptr;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1+(*arg2ptr++);
				return sum;
			};

			clMatrix4D operator- (const clMatrix4D &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)-(*argmatrixptr++);
				return sum;
			};
			clMatrix4D operator- (const T arg) const
			{
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)-arg;
				return sum;
			};
			friend clMatrix4D<T> operator- (const T arg1,const clMatrix4D<T>& arg2)
			{
				clMatrix4D<T> sum(arg2._dim1,arg2._dim2,arg2._dim3,arg2._dim4);
				T *sumptr=sum._startptr;
				const T *arg2ptr=arg2._startptr;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1-(*arg2ptr++);
				return sum;
			};

			clMatrix4D operator* (const clMatrix4D &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)*(*argmatrixptr++);
				return sum;
			};
			clMatrix4D operator* (const T arg) const
			{
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)*arg;
				return sum;
			};
			friend clMatrix4D<T> operator* (const T arg1,const clMatrix4D<T>& arg2)
			{
				clMatrix4D<T> sum(arg2._dim1,arg2._dim2,arg2._dim3,arg2._dim4);
				T *sumptr=sum._startptr;
				const T *arg2ptr=arg2._startptr;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1*(*arg2ptr++);
				return sum;
			};

			clMatrix4D operator/ (const clMatrix4D &arg) const
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)/(*argmatrixptr++);
				return sum;
			};
			clMatrix4D operator/ (const T arg) const
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clMatrix4D::operator/\n");
					exit(1);
				}
				T arginv=1.0/arg;
				clMatrix4D<T> sum(_dim1,_dim2,_dim3,_dim4);
				T *sumptr=sum._startptr;
				const T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*sumptr++)=(*matrixptr++)*arginv;
				return sum;
			};
			friend clMatrix4D<T> operator/ (const T arg1,const clMatrix4D<T>& arg2)
			{
				clMatrix4D<T> sum(arg2._dim1,arg2._dim2,arg2._dim3,arg2._dim4);
				T *sumptr=sum._startptr;
				const T *arg2ptr=arg2._startptr;
				for (int n=0; n<arg2._dim; n++) (*sumptr++)=arg1/(*arg2ptr++);
				return sum;
			};

			clMatrix4D& operator+= (const clMatrix4D &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)+=(*argmatrixptr++);
				return *this;
			};
			clMatrix4D& operator+= (const T arg)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)+=arg;
				return *this;
			};

			clMatrix4D& operator-= (const clMatrix4D &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)-=(*argmatrixptr++);
				return *this;
			};
			clMatrix4D& operator-= (const T arg)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)-=arg;
				return *this;
			};

			clMatrix4D& operator*= (const clMatrix4D &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)*=(*argmatrixptr++);
				return *this;
			};
			clMatrix4D& operator*= (const T arg)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)*=arg;
				return *this;
			};

			clMatrix4D& operator/= (const clMatrix4D &arg)
			{
#ifdef DEBUG
				test_indices(arg);
#endif
				T *matrixptr=_startptr;
				const T *argmatrixptr=arg._startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)/=(*argmatrixptr++);
				return *this;
			};
			clMatrix4D& operator/= (const T arg)
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg=0 in clMatrix4D::operator/=\n");
					exit(1);
				}
				T arginv=1.0/arg;
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)*=arginv;
				return *this;
			};

			clMatrix4D& operator++ (void)
			{
				T *matrixptr=_startptr;
				for (int n=0; n<_dim; n++) (*matrixptr++)++;
				return *this;
			};

			clMatrix4D& operator-- (void)
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
			int get_dim3(void) {return _dim3;};
			int get_dim3(void) const {return _dim3;};
			int get_dim4(void) {return _dim4;};
			int get_dim4(void) const {return _dim4;};
			int get_dim(void) {return _dim;};
			int get_dim(void) const {return _dim;};
			T* get_startptr(void) {return _startptr;};
			const T* get_startptr(void) const {return _startptr;};

			/*  Move data from arg to current */
			void move(clMatrix4D &arg)
			{
				if (this==&arg) return;

				if (_startptr) deallocate();

				_dim1=arg._dim1;
				_dim2=arg._dim2;
				_dim3=arg._dim3;
				_dim4=arg._dim4;
				_dim=arg._dim;
				_startptr=arg._startptr;
				_matrix.move(arg._matrix);

				arg._dim1=0;
				arg._dim2=0;
				arg._dim3=0;
				arg._dim4=0;
				arg._dim=0;
				arg._startptr=NULL;
			};

			/*  Converts a matrix to a regular 4D array, does not allow bounce checking  */
			void copy(T ****arg) const
			{
				if (!arg || !arg[0] || !arg[0][0] || arg[0][0][0]) {
					fprintf(stdout,"FATAL ERROR: arg not allocated in clMatrix4D::copy\n");
					exit(1);
				}
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: _startptr not allocated in clMatrix4D::copy\n");
					exit(1);
				}
				const T *matrixptr=_startptr;
				T *argptr=arg[0][0][0];
				for (int n=0; n<_dim; n++) (*argptr++)=(*matrixptr++);
			};

			/*  Allocate matrix of dimension (dim1,dim2,dim3,dim4)  */
			T* allocate(const int dim1,const int dim2,const int dim3,const int dim4)
			{
				if (_dim==dim1*dim2*dim3*dim4) {
					_dim1=dim1;
					_dim2=dim2;
					_dim3=dim3;
					_dim4=dim4;
				} else {
					if (_startptr) deallocate();
					_dim=dim1*dim2*dim3*dim4;
					if (!_dim) return _startptr;
					_dim1=dim1;
					_dim2=dim2;
					_dim3=dim3;
					_dim4=dim4;
					_matrix.allocate(_dim);
					_startptr=_matrix.get_startptr();
#ifdef 	ALLOCDEBUG
					_num_matrix++;
					fprintf(stdout,"matrix4d %p allocated %d size %dx%dx%dx%d\n",_startptr,_num_matrix,_dim1,_dim2,_dim3,_dim4);
#endif
				}
				return _startptr;
			};
			/*  Allocate matrix of dimension (dim1,dim2) and initialize to initval  */
			T* allocate(const int dim1,const int dim2,const int dim3,const int dim4,const T initval)
			{
				allocate(dim1,dim2,dim3,dim4);
				if (_startptr) {
					T *matrixptr=_startptr;
					for (int n=0; n<_dim; n++) (*matrixptr++)=initval;
				}
				return _startptr;
			};

			/*  Deallocate matrix  */
			void deallocate(void)
			{
				if (_startptr) {
#ifdef ALLOCDEBUG
					_num_matrix--;
					fprintf(stdout,"matrix4d %p deallocated %d size %dx%dx%dx%d\n",_startptr,_num_matrix,_dim1,_dim2,_dim3,_dim4);
#endif
					_startptr=NULL;
				}
				_matrix.deallocate();
				_dim1=0;
				_dim2=0;
				_dim3=0;
				_dim4=0;
				_dim=0;
			};

			/*  Allocate matrix of dimensions (dim1,dim2,dim3,dim4) and initialize to (1,2,3,...,dim1*dim2*dim3*dim4)  */
			clMatrix4D& index(const int dim1,const int dim2,const int dim3,const int dim4)
			{
				allocate(dim1,dim2,dim3,dim4);
				T *matrixptr=_startptr;
				for (int n=1; n<=_dim; n++) (*matrixptr++)=(T)(n);
				return *this;
			};
			/*  Initialize existing matrix to (1,2,3,...,dim1*dim2*dim3*dim4)  */
			clMatrix4D& index(void)
			{
				T *matrixptr=_startptr;
				for (int n=1; n<=_dim; n++) (*matrixptr++)=(T)(n);
				return *this;
			};

			/*  Return L2-norm of matrix  */
			T norm(void) const { return (T)_matrix.norm(); };

			/*  Return RMS of matrix  */
			double rms(void) const { return _matrix.rms(); };

			/*  Return max of matrix  */
			T max(void) const { return _matrix.max(); };
			/*  Return global index of max of matrix  */
			int maxindex(void) const { return _matrix.maxindex(); };
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

			/*  Return min of matrix  */
			T min(void) const { return _matrix.min(); };
			/*  Return global index of min of matrix  */
			int minindex(void) const { return _matrix.minindex(); };
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

			/*  Return mean of matrix  */
			double mean(void) const { return _matrix.mean(); };

			/*  Return stdv of matrix  */
			double stdv(void) const { return _matrix.stdv(); };

		private:

#ifdef DEBUG
			void test_index1(const int n1) const
			{
				if (n1<0 || n1>=_dim1) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix4D::test_index1: %d %d\n",n1,_dim1);
					exit(1);
				}
			};
			void test_index2(const int n2) const
			{
				if (n2<0 || n2>=_dim2) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix4D::test_index2: %d %d\n",n2,_dim2);
					exit(1);
				}
			};
			void test_index3(const int n3) const
			{
				if (n3<0 || n3>=_dim3) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix4D::test_index3: %d %d\n",n3,_dim3);
					exit(1);
				}
			};
			void test_index4(const int n4) const
			{
				if (n4<0 || n4>=_dim4) {
					fprintf(stdout,"FATAL ERROR: index out of range in clMatrix4D::test_index4: %d %d\n",n4,_dim4);
					exit(1);
				}
			};
			void test_indices(const clMatrix4D &arg) const
			{
				if (!_startptr || !arg._startptr) {
					fprintf(stdout,"FATAL ERROR: matrix or arg.matrix not allocated in clMatrix4D::test_indices\n");
					exit(1);
				}
				if (_dim1!=arg._dim1 || _dim2!=arg._dim2 || _dim3!=arg._dim3 || _dim4!=arg._dim4) {
					fprintf(stdout,"FATAL ERROR: different source and target dimension in clMatrix4D::test_indices\n");
				}
			};
#endif

	};

	/* clMatrix4d functions */

	template <class T> void bcast(clMatrix4D<T> &mptr,const int ip=0)
	{
#ifdef DIST_MPI
		int num_proc=clMpi_namespace::get_num_proc();
		int proc_id=clMpi_namespace::get_proc_id();
		if (num_proc<2) return;
		int dim1,dim2,dim3,dim4;
		T *matrix;
		if (proc_id==ip) {
			dim1=mptr.get_dim1();
			dim2=mptr.get_dim2();
			dim3=mptr.get_dim3();
			dim4=mptr.get_dim4();
		}
		clMpi_namespace::bcast(&dim1,1,ip);
		clMpi_namespace::bcast(&dim2,1,ip);
		clMpi_namespace::bcast(&dim3,1,ip);
		clMpi_namespace::bcast(&dim4,1,ip);
		if (!dim1 && !dim2 && !dim3 && !dim4) {
			mptr.deallocate();
			return;
		}
		if (proc_id!=ip) mptr.allocate(dim1,dim2,dim3,dim4);
		matrix=mptr.get_startptr();
		clMpi_namespace::bcast(matrix,dim1*dim2*dim3,dim4,ip);
#endif
	};

}; // namespace clVector_namespace

#endif // CLMATRIX4D_H__INCLUDED_
