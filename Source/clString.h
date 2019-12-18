/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/

#ifndef CLSTRING_H__INCLUDED_
#define CLSTRING_H__INCLUDED_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

namespace clString_namespace
{
	enum constants {MAX_LENGTH=256};

	class clString
	{
		private:

			int _len;		// String length (excluding the terminating '\0')
			char *_string;	// Character string
#ifdef ALLOCDEBUG
			static int _num_string;
#endif

		public:

			/*  Constructors  */
			clString(void) : _len(0), _string(NULL)
			{
				allocate(0);
			};
			explicit clString(const char* const init) : _len(0), _string(NULL)
			{
				allocate(init);
			};
			clString(const int len,const char* const init="") : _len(0), _string(NULL)
			{
				allocate(len,init);
			};

			/*  Destructor  */
			virtual ~clString(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clString(const clString &string) : _len(0), _string(NULL)
			{
				allocate(string._len,string._string);
			};

			/*  Default Assignment Constructor  */
			clString& operator= (const clString &string)
			{
				if (this==&string) return *this;

				allocate(string._len,string._string);
				return *this;
			};

			/*  Operator overloading functions  */

			/*  Converts a regular string to a clString string  */
			clString& operator= (const char* const chararray)
			{
#ifdef DEBUG
				if (!chararray || chararray=="") {
					fprintf(stdout,"FATAL ERROR: chararray not allocated in clString::operator=\n");
					exit(1);
				}
#endif
				allocate(chararray);
				return *this;
			};

			/* Compare to string */
			int operator== (const clString &string) const
			{
				if (!strcmp(_string,string._string)) return 1;
				else return 0;
			};
			int operator== (const char* const chararray) const
			{
				if (!strcmp(_string,chararray)) return 1;
				else return 0;
			};
			int operator!= (const clString &string) const
			{
				if (strcmp(_string,string._string)) return 1;
				else return 0;
			};
			int operator!= (const char* const chararray) const
			{
				if (strcmp(_string,chararray)) return 1;
				else return 0;
			};

			/* Catenate two strings */
			clString operator+ (const clString &string) const
			{
				int len=_len+string._len;
				clString sum(len,_string);
				strcat(sum._string,string._string);
				return sum;
			};
			clString operator+ (const char* const chararray) const
			{
#ifdef DEBUG
				if (!chararray || chararray=="") {
					fprintf(stdout,"FATAL ERROR: chararray not allocated in clString::operator+\n");
					exit(1);
				}
#endif
				int len=_len+strlen(chararray);
				clString sum(len,_string);
				strcat(sum._string,chararray);
				return sum;
			};
			friend clString operator+ (const char* const chararray,const clString& string)
			{
#ifdef DEBUG
				if (!chararray || chararray=="") {
					fprintf(stdout,"FATAL ERROR: chararray not allocated in clString::operator+\n");
					exit(1);
				}
#endif
				int len=string._len+strlen(chararray);
				clString sum(len,chararray);
				strcat(sum._string,string._string);
				return sum;
			};

			clString& operator+= (const clString &string)
			{
				int len=_len+string._len;
				char *sum=new char [len+1];
				strcpy(sum,_string);
				strcat(sum,string._string);
				allocate(sum);
				delete [] sum;
				return *this;
			};
			clString& operator+= (const char* const chararray)
			{
#ifdef DEBUG
				if (!chararray || chararray=="") {
					fprintf(stdout,"FATAL ERROR: chararray not allocated in clString::operator+=\n");
					exit(1);
				}
#endif
				int len=_len+strlen(chararray);
				char *sum=new char [len+1];
				strcpy(sum,_string);
				strcat(sum,chararray);
				allocate(sum);
				delete [] sum;
				return *this;
			};

			/*  convert to uppercase  */
			clString& operator++ (void)
			{
				unsigned char c;
				for (int i=0; _string[i]; i++) {
					c=_string[i];
					if (c>96 && c<123) {
						c-=32;
						_string[i]=(char)c;
					}
				}
				return *this;
			};

			/*  convert to lowercase  */
			clString& operator-- (void)
			{
				unsigned char c;
				for (int i=0; _string[i]; i++) {
					c=_string[i];
					if (c>64 && c<91) {
						c+=32;
						_string[i]=(char)c;
					}
				}
				return *this;
			};

			/* Deallocate a string */
			void deallocate(void)
			{
				if (_string) {
#ifdef ALLOCDEBUG
					_num_string--;
					fprintf(stdout,"string %s %p deallocated %d\n",_string,_string,_num_string);
#endif
					delete [] _string;
					_string=NULL;
				}
				_len=0;
			};

			/*  Remaining member functions  */

			int get_len(void) {return _len;};
			int get_len(void) const {return _len;};
			char* get_string(void) {return _string;};
			const char* get_string(void) const {return _string;};

			void copy_string(char* const,const int);
			int read_ascii(FILE*);
			int read_ascii(const char* const);
			void write_ascii(FILE*,const int nlflag=0) const;
			void write_ascii(const char* const,const int iotype=0,const int nlflag=0) const;
			int read(FILE*);
			int write(FILE*) const;

		private:

			char* allocate(const char* const init="")
			{
				if (!init || init=="") _len=0;
				else _len=strlen(init);
				_string=allocate(_len,init);
				return _string;
			};
			char* allocate(const int len,const char* const init="")
			{
				if (_string) deallocate();
				if (len<0) {
					fprintf(stdout,"FATAL ERROR: Illegal value of len in clString::allocate\n");
					exit(1);
				}
				_len=len;
				_string=new char [_len+1];
				if (!_string) {
					fprintf(stdout,"FATAL ERROR: Could not allocate _string in clString::allocate\n");
					exit(1);
				}
				if (!init || init=="") {
					_string[0]='\0';
				} else {
					strcpy(_string,init);
				}
#ifdef ALLOCDEBUG
				_num_string++;
				fprintf(stdout,"string %s %p len %d allocated %d\n",_string,_string,_num_string,_len);
#endif
				return _string;
			};

	};

	class clString_vector
	{
		private:

			int _dim;			// Dimension
			clString *_vector;	// Vector
#ifdef ALLOCDEBUG
			static int _num_string_vector;
#endif

		public:

			/*  Constructors  */
			clString_vector(void) :
				_dim(0),
				_vector(NULL)
			{ };
			explicit clString_vector(const int dim) :
				_dim(0),
				_vector(NULL)
			{
				allocate(dim);
			};

			/*  Destructor  */
			virtual ~clString_vector(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clString_vector(const clString_vector &arg) :
				_dim(0),
				_vector(NULL)
			{
				allocate(arg._dim);
				clString *vectorptr=_vector;
				const clString *argvectorptr=arg._vector;
				for (int n=0; n<_dim; n++) (*vectorptr++)=(*argvectorptr++);
			};

			/*  Default Assignment Constructor  */
			clString_vector& operator= (const clString_vector &arg)
			{
				if (this==&arg) return *this;

				if (!_vector) {
					if (_dim) {
						fprintf(stdout,"FATAL ERROR: illegal _dim in clString_vector::operator=\n");
						exit(1);
					}
					allocate(arg._dim);
				}
#ifdef DEBUG
				test_indices(arg);
#endif
				int dim=(_dim <= arg._dim ? _dim : arg._dim);
				clString *vectorptr=_vector;
				const clString *argvectorptr=arg._vector;
				for (int n=0; n<dim; n++) (*vectorptr++)=(*argvectorptr++);
				return *this;
			};

			/*  Set a string_vector to a regular array of strings (does not allow bounds checking)  */
			clString_vector& operator= (const clString *arg)
			{
				if (!_vector || !_dim) {
					fprintf(stdout,"FATAL ERROR: illegal _vector or _dim in clString_vector::operator=\n");
					exit(1);
				}
				clString *vectorptr=_vector;
				const clString *argstringptr=arg;
				for (int n=0; n<_dim; n++) (*vectorptr++)=(*argstringptr++);
				return *this;
			};

			/*  Operator overloading functions  */

			clString& operator[](const int n)
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clString_vector::operator[]\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};
			clString& operator[](const int n) const
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clString_vector::operator[]\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};

			clString& operator()(const int n)
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clString_vector::operator[]\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};
			clString& operator()(const int n) const
			{
#ifdef DEBUG
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clString_vector::operator[]\n");
					exit(1);
				}
				test_index(n);
#endif
				return _vector[n];
			};

			/*  Remaining member functions  */

			int get_dim(void) {return _dim;};
			int get_dim(void) const {return _dim;};
			clString *get_startptr(void) {return _vector;};
			const clString *get_startptr(void) const {return _vector;};
			int read(FILE*);
			int write(FILE*) const;
			void read_ascii(FILE*);
			void read_ascii(const char* const);
			void write_ascii(FILE*,const int nlflag=0) const;
			void write_ascii(const char* const,const int iotype=0,const int nlflag=0) const;

			/*  Converts a vector to a regular array, does not allow bounce checking  */
			void copy(clString* const arg)
			{
				if (!arg) {
					fprintf(stdout,"FATAL ERROR: arg not allocated in clString_vector::copy\n");
					exit(1);
				}
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: _vector not allocated in clString_vector::copy\n");
					exit(1);
				}
				const clString *vectorptr=_vector;
				clString *argptr=arg;
				for (int n=0; n<_dim; n++) (*argptr++)=(*vectorptr++);
			};

			/*  Allocate vector of dimension dim  */
			clString* allocate(const int dim)
			{
				if (_vector) deallocate();
				_dim=dim;
				if (!_dim) return _vector;
				_vector=new clString [_dim];
				if (!_vector) {
					fprintf(stdout,"FATAL ERROR: Could not allocate _vector in clString_vector::allocate\n");
					exit(1);
				}
#ifdef ALLOCDEBUG
				_num_string_vector++;
				fprintf(stdout,"string vector %p size %d allocated %d\n",_vector,_dim,_num_string_vector);
#endif
				return _vector;
			};

			/*  Reallocate vector of dimension dim  */
			/*  Data that does not fit is lost  */
			/*  Zero out cells that are not used yet */
			clString* redimension(const int new_dim)
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
					fprintf(stdout,"FATAL ERROR: Illegal dimension in clString_vector::redimension\n");
					exit(1);
				}
				clString *newvector=new clString [new_dim];
				if (!newvector) {
					fprintf(stdout,"FATAL ERROR: Could not allocate newvector in clString_vector::redimension\n");
					exit(1);
				}
				const int dim=(_dim <= new_dim ? _dim : new_dim);
				clString *newvectorptr=newvector;
				const clString *vectorptr=_vector;
				int n=0;
				for (; n<dim; n++) (*newvectorptr++)=(*vectorptr++);
				clString *swap=_vector;
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
					_num_string_vector--;
					fprintf(stdout,"string vector %p size %d deallocated %d\n",_vector,_dim,_num_string_vector);
#endif
					delete [] _vector;
					_vector=NULL;
				}
				_dim=0;
			};

#ifdef DEBUG
			void test_index(const int n) const
			{
				if (n<0 || n>=_dim) {
					fprintf(stdout,"FATAL ERROR: index out of range in clString_vector::test_index %d %d\n",n,_dim);
					exit(1);
				}
			};
			void test_indices(const clString_vector &arg) const
			{
				if (!_vector || !arg._vector) {
					fprintf(stdout,"FATAL ERROR: _vector or arg._vector not allocated in clString_vector::test_indices\n");
					exit(1);
				}
				if (_dim!=arg._dim) {
					fprintf(stdout,"WARNING: different source and target dimension in clString_vector::test_indices\n");
					fprintf(stdout,"  use smaller dimension of %d and %d\n",_dim,arg._dim);
				}
			};
#endif

	};

	class clString_matrix
	{
		private:

			int _dim1,_dim2,_dim;		// Dimension
			clString *_startptr;		// Matrix start pointer
			clString_vector _matrix;	// Matrix stored as a single string_vector
#ifdef ALLOCDEBUG
			static int _num_string_matrix;
#endif

		public:

			/*  Constructors  */
			clString_matrix(void) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_startptr(NULL)
			{ };
			explicit clString_matrix(const int dim1,const int dim2) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_startptr(NULL)
			{
				allocate(dim1,dim2);
			};

			/*  Destructor  */
			virtual ~clString_matrix(void)
			{
				deallocate();
			};

			/*  Default Copy Constructor  */
			clString_matrix(const clString_matrix &arg) :
				_dim1(0),
				_dim2(0),
				_dim(0),
				_startptr(NULL)
			{
				_dim1=arg._dim1;
				_dim2=arg._dim2;
				_dim=arg._dim;
				_matrix=arg._matrix;
				_startptr=_matrix.get_startptr();
			};

			/*  Default Assignment Constructor  */
			clString_matrix& operator= (const clString_matrix &arg)
			{
				if (this==&arg) return *this;

				if (!_startptr) {
					if (_dim1 || _dim2 || _dim) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clString_matrix::operator=\n");
						exit(1);
					}
					allocate(arg._dim1,arg._dim2);
				} else {
					if (_dim1!=arg._dim1 || _dim2!=arg._dim2 || _dim!=arg._dim) {
						fprintf(stdout,"FATAL ERROR: illegal dimensions in clString_matrix::operator=\n");
						exit(1);
					}
				}
				_matrix=arg._matrix;
				return *this;
			};

			/*  Operator overloading functions  */

			clString* operator[](const int n1)
			/* This allows you to use value=matrix[n1][n2], however it prevents column bounds checking */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clString_matrix::operator[]\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim2*n1);
			};

			clString* operator()(const int n1)
			/* Same as using [] */
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clString_matrix::operator()\n");
					exit(1);
				}
				test_index1(n1);
#endif
				return (_startptr+_dim2*n1);
			};
			clString& operator()(const int n1,const int n2)
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clString_matrix::operator[]\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				return *(_startptr+n2+_dim2*n1);
			};
			clString& operator()(const int n1,const int n2) const
			{
#ifdef DEBUG
				if (!_startptr) {
					fprintf(stdout,"FATAL ERROR: matrix not allocated in clString_matrix::operator[]\n");
					exit(1);
				}
				test_index1(n1);
				test_index2(n2);
#endif
				return *(_startptr+n2+_dim2*n1);
			};

			/*  Remaining member functions  */

			int get_dim1(void) {return _dim1;};
			int get_dim1(void) const {return _dim1;};
			int get_dim2(void) {return _dim2;};
			int get_dim2(void) const {return _dim2;};
			int get_dim(void) {return _dim;};
			int get_dim(void) const {return _dim;};
			clString *get_startptr(void) {return _startptr;};
			const clString *get_startptr(void) const {return _startptr;};

			/*  Allocate string_matrix of dimension (dim1,dim2)  */
			clString* allocate(const int dim1,const int dim2)
			{
				if (_startptr) deallocate();
				_dim=dim1*dim2;
				if (!_dim) return _startptr;
				_dim1=dim1;
				_dim2=dim2;
				_matrix.allocate(_dim);
				_startptr=_matrix.get_startptr();
#ifdef ALLOCDEBUG
				_num_string_matrix++;
				fprintf(stdout,"string_matrix %p size %dx%d allocated %d\n",_startptr,_dim1,_dim2,_num_string_matrix);
#endif
				return _startptr;
			};

			/*  Deallocate string_matrix  */
			void deallocate(void)
			{
				if (_startptr) {
#ifdef ALLOCDEBUG
					_num_string_matrix--;
					fprintf(stdout,"string_matrix %p size %dx%d deallocated %d\n",_startptr,_dim1,_dim2,_num_string_matrix);
#endif
					_startptr=NULL;
				}
				_matrix.deallocate();
				_dim1=0;
				_dim2=0;
				_dim=0;
			};

#ifdef DEBUG
			void test_index1(const int n1) const
			{
				if (n1<0 || n1>=_dim1) {
					fprintf(stdout,"FATAL ERROR: index out of range in clString_matrix::test_index1 %d %d\n",n1,_dim1);
					exit(1);
				}
			};
			void test_index2(const int n2) const
			{
				if (n2<0 || n2>=_dim2) {
					fprintf(stdout,"FATAL ERROR: index out of range in clString_matrix::test_index2 %d %d\n",n2,_dim2);
					exit(1);
				}
			};
			void test_indices(const clString_matrix &arg) const
			{
				if (!_startptr || !arg._startptr) {
					fprintf(stdout,"FATAL ERROR: matrix or arg.matrix not allocated in clString_matrix::test_indices\n");
					exit(1);
				}
				if (_dim1!=arg._dim1 || _dim2!=arg._dim2) {
					fprintf(stdout,"FATAL ERROR: different source and target dimension in clString_matrix::test_indices\n");
				}
			};
#endif

	};

	void bcast(clString *sptr,const int n=1,const int ip=0);
	void bcast(clString_vector&,const int ip=0);
	void bcast(clString_matrix&,const int ip=0);

}; // namespace clString_namespace

#endif // CLSTRING_H__INCLUDED_
