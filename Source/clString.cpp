/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/
#include "clString.h"

/*  Copy a string into a character array  */
void clString_namespace::clString::copy_string(char* const charstring,const int maxlen)
{
	if (_len<maxlen) { 
		strcpy(charstring,_string);
	} else {
		fprintf(stdout,"FATAL ERROR: charstring exceeds max length in clString::copy_string");
		exit(1);
	}
};

/*  Read a string from an ascii stream  */
int clString_namespace::clString::read_ascii(FILE* fileptr)
{
	if (fileptr) {
		char charstring[MAX_LENGTH];
		if (fscanf(fileptr,"%s",charstring)!=1) {
			fprintf(stdout,"WARNING: Could not read from stream in clString::read_ascii\n");
			return 0;
		}
		allocate(charstring);
	} else {
		fprintf(stdout,"FATAL ERROR: Could not open stream in clString::read_ascii");
		exit(1);
	}
	return _len;
};
/*  Read a string from an ascii file  */
int clString_namespace::clString::read_ascii(const char* const filename)
{
	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) {
		fprintf(stdout,"FATAL ERROR: Could not open %s\n",filename);
		exit(1);
	}
	read_ascii(fileptr);
	fclose(fileptr);
	return _len;
};

/*  Write a string to an ascii stream  */
/*    nlflag=0:  No newline after writing the string (default)  */
/*    nlflag!=0: Write a newline after writing the string  */
void clString_namespace::clString::write_ascii(FILE* fileptr,const int nlflag) const
{
	if (fileptr) {
		fprintf(fileptr,"%s",_string);
		if (nlflag) fprintf(fileptr,"\n");
	} else {
		fprintf(stdout,"FATAL ERROR: Could not open stream in clString::write_ascii");
		exit(1);
	}
};
/*  Write a string to an ascii file  */
/*    append_flag=0:  Create a new file (default)  */
/*    append_flag!=0: Append to file if existing  */
/*    nlflag=0:  No newline after writing the string (default)  */
/*    nlflag!=0: Write a newline after writing the string  */
void clString_namespace::clString::write_ascii(const char* const filename,const int append_flag,const int nlflag) const
{
	FILE *fileptr=NULL;
	if (!append_flag) fileptr=fopen(filename,"w");
	else fileptr=fopen(filename,"a");
	if (!fileptr) {
		fprintf(stdout,"FATAL ERROR: Could not open %s\n",filename);
		exit(1);
	}
	fprintf(fileptr,"%s",_string);
	if (nlflag) fprintf(fileptr,"\n");
	fclose(fileptr);
};

/*  Read a string from a binary stream  */
/*  Return 0 for failure or 1 for success  */
int clString_namespace::clString::read(FILE* fileptr)
{
	int len;
	void *io_ptr;
	io_ptr=(void *)&len;
	if (fread(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
	allocate(len);
	if (_len) {
		io_ptr=(void *)_string;
		if (fread(io_ptr,sizeof(char),_len,fileptr)!=(unsigned)_len) return 0;
		_string[_len]='\0';
	}
	return 1;
};

/*  Write a string to a binary stream  */
/*  Return 0 for failure or 1 for success  */
int clString_namespace::clString::write(FILE* fileptr) const
{
	void *io_ptr;
	io_ptr=(void *)&_len;
	if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
	if (_len) {
		io_ptr=(void *)_string;
		if (fwrite(io_ptr,sizeof(char),_len,fileptr)!=(unsigned)_len) return 0;
	}
	return 1;
};

/*  Read a vector of strings from an ascii stream  */
/*  The first entry in the file must be the number of strings, followed by the strings  */
void clString_namespace::clString_vector::read_ascii(FILE *fileptr)
{
	if (fileptr) {
		int dim;
		if (fscanf(fileptr,"%d",&dim)!=1) {
			fprintf(stdout,"FATAL ERROR: Could not read from stream in clString_vector::read_ascii");
			exit(1);
		}
		if (!dim) return;
		allocate(dim);
		clString *vectorptr=_vector;
		for (int n=0; n<_dim; n++,vectorptr++) vectorptr->read_ascii(fileptr);
	} else {
		fprintf(stdout,"FATAL ERROR: Could not open stream in clString_vector::read_ascii");
		exit(1);
	}
};
/*  Read a vector of strings from an ascii file  */
/*  The first entry in the file must be the number of strings, followed by the strings  */
void clString_namespace::clString_vector::read_ascii(const char* const filename)
{
	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) {
		fprintf(stdout,"FATAL ERROR: Could not open %s\n",filename);
		exit(1);
	}
	read_ascii(fileptr);
	fclose(fileptr);
};

/*  Write a vector of strings to an ascii file  */
/*  The first entry written to file is the number of strings, followed by the strings  */
/*    nlflag=0:  Write a space between each string (default)  */
/*    nlflag!=0: Write a newline after writing each string  */
void clString_namespace::clString_vector::write_ascii(FILE *fileptr,const int nlflag) const
{
	fprintf(fileptr,"%d\n",_dim);
	if (!_dim) return;
	clString *vectorptr=_vector;
	for (int n=0; n<_dim; n++,vectorptr++) {
		vectorptr->write_ascii(fileptr,nlflag);
		if (!nlflag) fprintf(fileptr," ");
	}
};

/*  Write a vector of strings to an ascii file  */
/*  The first entry written to file is the number of strings, followed by the strings  */
/*    append_flag=0:  Create a new file (default)  */
/*    append_flag!=0: Append to file if existing  */
/*    nlflag=0:  Write a space between each string (default)  */
/*    nlflag!=0: Write a newline after writing each string  */
void clString_namespace::clString_vector::write_ascii(const char* const filename,const int append_flag,const int nlflag) const
{
	FILE *fileptr=NULL;
	if (!append_flag) fileptr=fopen(filename,"w");
	else fileptr=fopen(filename,"a");
	if (!fileptr) {
		fprintf(stdout,"FATAL ERROR: Could not open %s\n",filename);
		exit(1);
	}
	fprintf(fileptr,"%d\n",_dim);
	if (_dim) {
		clString *vectorptr=_vector;
		for (int n=0; n<_dim; n++,vectorptr++) {
			vectorptr->write_ascii(fileptr,nlflag);
			if (!nlflag) fprintf(fileptr," ");
		}
	}
	fclose(fileptr);
};

/*  Read a vector of strings from a binary stream  */
/*  Return 0 for failure or 1 for success  */
int clString_namespace::clString_vector::read(FILE* fileptr)
{
	int dim;
	void *io_ptr;
	io_ptr=(void *)&dim;
	if (fread(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
	allocate(dim);
	if (_dim) {
		clString *stringptr=get_startptr();
		for (int n=0; n<_dim; n++,stringptr++) if (!(stringptr->read(fileptr))) return 0;
	}
	return 1;
};

/*  Write a vector of strings to a binary stream  */
/*  Return 0 for failure or 1 for success  */
int clString_namespace::clString_vector::write(FILE* fileptr) const
{
	void *io_ptr;
	io_ptr=(void *)&_dim;
	if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) return 0;
	if (_dim) {
		const clString *stringptr=get_startptr();
		for (int n=0; n<_dim; n++,stringptr++) if (!(stringptr->write(fileptr))) return 0;
	}
	return 1;
};

/*  Broadcast an array of n strings from processor ip  */
void clString_namespace::bcast(clString *sptr,const int n,const int ip)
{
#ifdef DIST_MPI
	int num_proc=clMpi_namespace::get_num_proc();
	int proc_id=clMpi_namespace::get_proc_id();
	if (num_proc<2) return;
	int len;
	char *string=NULL;
	clString *ssptr=sptr;
	for (int i=0; i<n; i++,ssptr++) {
		if (proc_id==ip) len=ssptr->get_len()+1;
		clMpi_namespace::bcast(&len,1,ip);
		string=new char [len];
		if (!string) {
			fprintf(stdout,"FATAL ERROR: Could not allocate string in clString_namespace::bcast");
			exit(1);
		}
		if (proc_id==ip) ssptr->copy_string(string,len);
		clMpi_namespace::bcast(string,len,ip);
		if (proc_id!=ip) (*ssptr)=string;
		delete [] string;
	}
#endif
};

/*  Broadcast an vector of strings from processor ip  */
void clString_namespace::bcast(clString_vector &sptr,const int ip)
{
#ifdef DIST_MPI
	int num_proc=clMpi_namespace::get_num_proc();
	int proc_id=clMpi_namespace::get_proc_id();
	if (num_proc<2) return;
	int dim;
	clString *vector;
	if (proc_id==ip) dim=sptr.get_dim();
	clMpi_namespace::bcast(&dim,1,ip);
	if (!dim) {
		sptr.deallocate();
		return;
	}
	if (proc_id!=ip) sptr.allocate(dim);
	vector=sptr.get_startptr();
	bcast(vector,dim,ip);
#endif
};

/*  Broadcast an matrix of strings from processor ip  */
void clString_namespace::bcast(clString_matrix &mptr,const int ip)
{
#ifdef DIST_MPI
	int num_proc=clMpi_namespace::get_num_proc();
	int proc_id=clMpi_namespace::get_proc_id();
	if (num_proc<2) return;
	int dim1,dim2;
	clString *matrix;
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
	bcast(matrix,dim1*dim2,ip);
#endif
};
