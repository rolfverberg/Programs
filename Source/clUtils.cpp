/*******************************************************************************
  A collection of general utility functions
  Author: Rolf Verberg
*******************************************************************************/
#include "clUtils.h"

namespace clUtils_namespace
{
	int _little_endian_io=1;
	int _n_warn=0;
	int _size_my_int=sizeof(my_int);
	int _size_my_float=sizeof(my_float);
	int _size_my_double=sizeof(my_double);
}

/*  Convert degrees to radians  */
double clUtils_namespace::degtorad(const double x)
{
	return (x*PI/180.0);
};
/*  Convert radians to degrees  */
double clUtils_namespace::radtodeg(const double x)
{
	return (x*180.0/PI);
};

/*  Remove the file extension including the dot */
/*  Return 0 for successful completion  */
/*  Return 1 on failure  */
int clUtils_namespace::remove_fileextension(char *name)
{
	char *index=strrchr(name,'.');
	if (!index) {
		return 1;
	} else {
		int n=index-name;
		name[n]='\0';
	}
	return 0;
};

/*  Overwrite the file extension  */
/*  Return 0 for successful completion  */
/*  Return 1 on failure  */
int clUtils_namespace::fileextension(char *name,const char* const ext,const int lenmax)
{
	int namelength=strlen(name);
	int extlength=strlen(ext);
	char *newname=NULL;
	new_cstring(&newname,lenmax);
	char *index=strrchr(name,'.');
	if (!index) {
		if (namelength+extlength+2>lenmax) return 1;
		strcpy(newname,name);
		strcat(newname,".");
	} else {
		int n=index-name+1;
		if (n+extlength+1>lenmax) return 1;
		strncpy(newname,name,n);
		newname[n]='\0';
	}
	strcat(newname,ext);
	strcpy(name,newname);
	delete [] newname;
	return 0;
};
/*  Overwrite a three digit extension with a three digit file number label  */
/*  Return 0 for successful completion  */
/*  Return 1 if the string length is ill-defined or if the file number is out of bounds  */
int clUtils_namespace::filename3(char *name,int file_number)
{
	if (file_number<0 || file_number>999) return 0;
	int index=(int)strlen(name)-3;
	if (index<0) return 0;
	name[index++]=_name0[file_number/100];
	file_number%=100;
	name[index++]=_name0[file_number/10];
	file_number%=10;
	name[index++]=_name0[file_number];
	name[index]='\0';
	return 0;
};
/*  Construct a file name labeled with a three digit file number  */
/*  Return 0 for successful completion  */
/*  Return 1 if the maximum string length is exceeded or if the file number is out of bounds  */
int clUtils_namespace::filename3(char *name,int file_number,const int lenmax)
{
	if (file_number<0 || file_number>999) return 1;
	int index=(int)strlen(name);
	name[index++]=_name0[file_number/100];
	if (index==lenmax) return 1;
	file_number%=100;
	name[index++]=_name0[file_number/10];
	if (index==lenmax) return 1;
	file_number%=10;
	name[index++]=_name0[file_number];
	if (index==lenmax) return 1;
	name[index]='\0';
	return 0;
};
/*  Add a three digit extension with file number to a base file name  */
/*  Leave base file name unmodified  */
/*  Return 0 for successful completion  */
/*  Return 1 if the maximum string length is exceeded or if the file number is out of bounds  */
int clUtils_namespace::filename3(char *name,const char* const basename,int file_number,const int lenmax)
{
	if (file_number<0 || file_number>999) return 1;
	strcpy(name,basename);
	int index=(int)strlen(name);
	if (index>=lenmax) return 1;
	if (name[index-1]!='.') {
		name[index++]='.';
		if (index>=lenmax) return 1;
	}
	name[index++]=_name0[file_number/100];
	if (index>=lenmax) return 1;
	file_number%=100;
	name[index++]=_name0[file_number/10];
	if (index>=lenmax) return 1;
	file_number%=10;
	name[index++]=_name0[file_number];
	if (index>=lenmax) return 1;
	name[index]='\0';
	return 0;
};
/*  Add a file number to a base file name, but do not fix the number of digits  */
/*  Leave base file name unmodified  */
/*  Return 0 for successful completion  */
/*  Return 1 if the maximum string length is exceeded or if the file number is out of bounds  */
int clUtils_namespace::filename(char *name,const char* const basename,int file_number,const int lenmax,const int adddot)
{
	if (file_number<0) return 1;
	int len=strlen(basename);
	if (!file_number) len++;
	else len+=int(1+log10(file_number));
	if (adddot) len++;
	if (len>=lenmax) return 1;
	if (adddot) sprintf(name,"%s.%d",basename,file_number);
	else sprintf(name,"%s%d",basename,file_number);
	return 0;
};

/*  Print a warning message */
void clUtils_namespace::warning(const char* const warning_msg)
{
	fprintf(stderr,"WARNING: %s\n",warning_msg);
	_n_warn++;
	if (_n_warn>=MAX_WARN) fatal_err("maximum number of warnings exceeded");
};
void clUtils_namespace::warning(const char* const warning_msg1,const char* const warning_msg2)
{
	fprintf(stderr,"WARNING: %s %s\n",warning_msg1,warning_msg2);
	_n_warn++;
	if (_n_warn>=MAX_WARN) fatal_err("maximum number of warnings exceeded");
};

/*  Fatal error handler  */
void clUtils_namespace::fatal_err(const char* const error_msg)
{
	fprintf(stderr,"\nFATAL ERROR: %s\n",error_msg);
	fflush(stderr);
	fflush(stdout);
	exit(1);
};
void clUtils_namespace::fatal_err(const char* const error_msg,const int ierror)
{
	fprintf(stderr,"\nFATAL ERROR: %s %d\n",error_msg,ierror);
	fflush(stderr);
	fflush(stdout);
	exit(1);
};
void clUtils_namespace::fatal_err(const int flag,const char* const error_msg)
{
	if (flag==0)
		fprintf(stderr,"\nFATAL ERROR: Could not allocate %s\n",error_msg);
	else if (flag==1)
		fprintf(stderr,"\nFATAL ERROR: Could not open %s\n",error_msg);
	else if (flag==2)
		fprintf(stderr,"\nFATAL ERROR: Read error from %s\n",error_msg);
	else if (flag==3)
		fprintf(stderr,"\nFATAL ERROR: Write error to %s\n",error_msg);
	else
		fprintf(stderr,"\nFATAL ERROR: %s\n",error_msg);
	fflush(stderr);
	fflush(stdout);
	exit(1);
};

/*  Allocate an array of dim strings of length size  */
char** clUtils_namespace::chararray(char **m,const int dim,const int size)
{
	int i;
	m=new char* [dim];
	if (!m) return NULL;
	for (i=0; i<dim; i++) {
		m[i]=new char [size];
		if (!m[i]) return NULL;
		strcpy(m[i],"");
	}
	return m;
};
/*  Deallocate an array of dim strings of length size  */
char** clUtils_namespace::delete_chararray(char **m,const int dim)
{
	int i;
	for (i=0; i<dim; i++) delete [] m[i];
	delete [] m;
	m=NULL;
	return m;
};

/*  Allocate a lower triangular matrix with subscript range m[0..dim-1][0..dim-1]  */
/*  without storing the elements above the diagonal  */
/*  Store in row based indexing  */
double** clUtils_namespace::matrix_low_row(double **m,const int dim)
{
	int i,n=(dim*(dim+1))/2;
	m=new double* [dim];
	if (!m) return NULL;
	m[0]=new double [n];
	if (!m[0]) return NULL;
	for (i=1; i<dim; i++) m[i]=m[i-1]+i;
	return m;
};

/*  Allocate a lower triangular matrix with subscript range m[0..dimi-1][0..dimj-1]  */
/*  without storing the elements above the diagonal  */
/*  Store in row based mpi row-wrapped indexing  */
double** clUtils_namespace::matrix_low_row_mpi(double **m,const int dimi,const int* const dimj)
{
	int i,n=0;
	for (i=0; i<dimi; i++) n+=dimj[i];
	m=new double* [dimi];
	if (!m) return NULL;
	m[0]=new double [n];
	if (!m[0]) return NULL;
	for (i=1; i<dimi; i++) m[i]=m[i-1]+dimj[i-1];
	return m;
};

/*  Allocate a lower triangular matrix with subscript range m[0..dim-1][0..dim-1]  */
/*  without storing the elements above the diagonal  */
/*  Store in column based indexing  */
double** clUtils_namespace::matrix_low_col(double **m,const int dim)
{
	int i,j,n=(dim*(dim+1))/2;
	m=new double* [dim];
	if (!m) return NULL;
	m[0]=new double [n];
	if (!m[0]) return NULL;
	for (i=1,j=dim; i<dim; i++,j--) m[i]=m[i-1]+j;
	return m;
};

/*  Allocate a 3D matrix with subscript range m[0..dim[0]-1][0..dim[1]-1][0..dim[2]-1] */
/*  that is lower triangular matrix in the last two deminsions  */
/*  without storing the elements above the diagonal  */
/*  Store in row based indexing  */
double*** clUtils_namespace::matrix_low_row(double ***m,const int dim[3])
{
	int i,j,n=(dim[1]*(dim[1]+1))/2;
	m=new double** [dim[0]];
	if (!m) return NULL;
	m[0]=new double* [dim[0]*dim[1]];
	if (!m[0]) return NULL;
	m[0][0]=new double [dim[0]*n];
	if (!m[0][0]) return NULL;
	for (j=1; j<dim[1]; j++) m[0][j]=m[0][j-1]+j;
	for (i=1; i<dim[0]; i++) {
		m[i]=m[i-1]+dim[1];
		m[i][0]=m[i-1][0]+n;
		for (j=1; j<dim[1]; j++) m[i][j]=m[i][j-1]+j;
	}
	return m;
};
double*** clUtils_namespace::matrix_low_row(double ***m,const int idim,const int jkdim)
{
	int dim[3];
	dim[0]=idim;
	dim[1]=jkdim;
	dim[2]=jkdim;
	m=matrix_low_row(m,dim);
	return m;
};

/*  Allocate a 3D matrix with subscript range m[0..dimi-1][0..dimj-1][0..dimk-1]  */
/*  that is symmetric in the first two indices without storing the elements above that diagonal  */
/*  Store in row based mpi row-wrapped indexing for the first two indices **/
double*** clUtils_namespace::matrix_low_row_mpi(double ***m,const int dimi,const int* const dimj,const int dimk)
{
	int i,j,n=0;
	for (i=0; i<dimi; i++) n+=dimj[i];
	m=new double** [dimi];
	if (!m) return NULL;
	m[0]=new double* [n];
	if (!m[0]) return NULL;
	m[0][0]=new double [n*dimk];
	if (!m[0][0]) return NULL;
	for (i=1; i<dimi; i++) m[i]=m[i-1]+dimj[i-1];
	for (j=1; j<n; j++) m[0][j]=m[0][j-1]+dimk;
	return m;
};

/*  Print arrays or matrices to file  */
void clUtils_namespace::print_matrix_low_col(const char* const filename,const double* const* const m,int dim)
{
	FILE *fileptr=fopen(filename,"w");
	if (!fileptr) fatal_err("Could not open filename print_matrix_low_col");
	int i,j;
	for (i=0; i<dim; i++) {
		for (j=0; j<=i; j++) 
			fprintf(fileptr," % 24.16e",m[j][i-j]);
//		for (j=i; j<dim; j++) 
//			fprintf(fileptr," % 24.16e",m[i][j-i]);
		fprintf(fileptr,"\n");
	}
	fclose(fileptr);
};

/*  Allocate a 3D vector field with subscript range m[0..dim[0]-1][0..dim[1]-1][0..dim[2]-1]  */
clUtils_namespace::vector3d*** clUtils_namespace::vectorfield(vector3d ***m,const int dim[3])
{
	int i,j;
	m=new vector3d** [dim[0]];
	if (!m) return NULL;
	m[0]=new vector3d* [dim[0]*dim[1]];
	if (!m[0]) return NULL;
	m[0][0]=new vector3d [dim[0]*dim[1]*dim[2]];
	if (!m[0][0]) return NULL;
	for (j=1; j<dim[1]; j++) m[0][j]=m[0][j-1]+dim[2];
	for (i=1; i<dim[0]; i++) {
		m[i]=m[i-1]+dim[1];
		m[i][0]=m[i-1][0]+dim[1]*dim[2];
		for (j=1; j<dim[1]; j++) m[i][j]=m[i][j-1]+dim[2];
	}
	return m;
};
/*  Allocate a 3D vector field with subscript range m[0..idim-1][0..jdim-1][0..dimk-1]  */
clUtils_namespace::vector3d*** clUtils_namespace::vectorfield(vector3d ***m,const int idim,const int jdim,const int kdim)
{
	int dim[3];
	dim[0]=idim;
	dim[1]=jdim;
	dim[2]=kdim;
	m=vectorfield(m,dim);
	return m;
};

/*  Allocate a 2D function field with subscript range m[0..dim[0]-1][0..dim[1]-1],  */
/*  where each entry is an array of nfunc function values  */
clUtils_namespace::function2d** clUtils_namespace::functionfield(function2d **m,const int dim[2],const int nfunc)
{
	int i,n,vol;
	function2d *fptr;
	vol=dim[0]*dim[1];
	m=new function2d* [dim[0]];
	if (!m) return NULL;
	m[0]=new function2d [vol];
	if (!m[0]) return NULL;
	for (i=1; i<dim[0]; i++) m[i]=m[i-1]+dim[1];
	for (fptr=m[0],n=0; n<vol; n++,fptr++) {
		fptr->x=0.0;
		fptr->y=0.0;
		fptr->area=0.0;
		if (nfunc>0) {
			fptr->f=new double [nfunc];
			if (!fptr->f) return NULL;
		} else fptr->f=NULL;
	}
	return m;
};
/*  Allocate a 2D function field with subscript range m[0..idim-1][0..jdim-1],  */
/*  where each entry is an array of nfunc function values  */
clUtils_namespace::function2d** clUtils_namespace::functionfield(function2d **m,const int idim,const int jdim,const int nfunc)
{
	int dim[2];
	dim[0]=idim;
	dim[1]=jdim;
	m=functionfield(m,dim,nfunc);
	return m;
};
/*  Deallocate a 2D function field  */
clUtils_namespace::function2d** clUtils_namespace::delete_functionfield(function2d **m)
{
	delete [] m[0];
	delete [] m;
	m=NULL;
	return m;
};

/*  Allocate a 3D function field with subscript range m[0..dim[0]-1][0..dim[1]-1][0..dim[2]-1],  */
/*  where each entry is an array of nfunc function values  */
clUtils_namespace::function3d*** clUtils_namespace::functionfield(function3d ***m,const int dim[3],const int nfunc)
{
	int i,j,n,vol;
	function3d *fptr;
	vol=dim[0]*dim[1]*dim[2];
	m=new function3d** [dim[0]];
	if (!m) return NULL;
	m[0]=new function3d* [dim[0]*dim[1]];
	if (!m[0]) return NULL;
	m[0][0]=new function3d [vol];
	if (!m[0][0]) return NULL;
	for (j=1; j<dim[1]; j++) m[0][j]=m[0][j-1]+dim[2];
	for (i=1; i<dim[0]; i++) {
		m[i]=m[i-1]+dim[1];
		m[i][0]=m[i-1][0]+dim[1]*dim[2];
		for (j=1; j<dim[1]; j++) m[i][j]=m[i][j-1]+dim[2];
	}
	for (fptr=m[0][0],n=0; n<vol; n++,fptr++) {
		fptr->x=0.0;
		fptr->y=0.0;
		fptr->z=0.0;
		fptr->vol=0.0;
		if (nfunc>0) {
			fptr->f=new double [nfunc];
			if (!fptr->f) return NULL;
		} else fptr->f=NULL;
	}
	return m;
};
/*  Allocate a 3D function field with subscript range m[0..idim-1][0..jdim-1][0..dimk-1],  */
/*  where each entry is an array of nfunc function values  */
clUtils_namespace::function3d*** clUtils_namespace::functionfield(function3d ***m,const int idim,const int jdim,const int kdim,const int nfunc)
{
	int dim[3];
	dim[0]=idim;
	dim[1]=jdim;
	dim[2]=kdim;
	m=functionfield(m,dim,nfunc);
	return m;
};
/*  Deallocate a 3D function field  */
clUtils_namespace::function3d*** clUtils_namespace::delete_functionfield(function3d ***m)
{
	delete [] m[0][0];
	delete [] m[0];
	delete [] m;
	m=NULL;
	return m;
};

/*  Get the distance between two points in num_dim dimensions (default 3)  */
double clUtils_namespace::distance(const double* const r1,const double* const r2,const int num_dim)
{
	double rsq=0.0,dummy;
	for (int n_dim=0; n_dim<num_dim; n_dim++) {
		dummy=r1[n_dim]-r2[n_dim];
		rsq+=dummy*dummy;
	}
	return (rsq==0 ? 0 : sqrt(rsq));
};
double clUtils_namespace::distance(const double* const r1,const double* const r2,double &rsq,const int num_dim)
{
	double dummy;
	rsq=0.0;
	for (int n_dim=0; n_dim<num_dim; n_dim++) {
		dummy=r1[n_dim]-r2[n_dim];
		rsq+=dummy*dummy;
	}
	return (rsq==0 ? 0 : sqrt(rsq));
};
double clUtils_namespace::distance(const vector3d* const r1,const vector3d* const r2)
{
	double rsq=(r1->x-r2->x)*(r1->x-r2->x)+(r1->y-r2->y)*(r1->y-r2->y)+(r1->z-r2->z)*(r1->z-r2->z);
	return (rsq==0 ? 0 : sqrt(rsq));
};
double clUtils_namespace::distance(const vector3d* const r1,const vector3d* const r2,double &rsq)
{
	rsq=(r1->x-r2->x)*(r1->x-r2->x)+(r1->y-r2->y)*(r1->y-r2->y)+(r1->z-r2->z)*(r1->z-r2->z);
	return (!rsq ? 0 : sqrt(rsq));
};

/*  Get the innerproduct of two 3D vectors */
double clUtils_namespace::innerproduct(const double* const a,const double* const b)
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
};

/*  Get the outerproduct of two 3D vectors */
void clUtils_namespace::outerproduct(const double* const a,const double* const b,double* const axb)
{
	axb[0]=a[1]*b[2]-a[2]*b[1];
	axb[1]=a[2]*b[0]-a[0]*b[2];
	axb[2]=a[0]*b[1]-a[1]*b[0];
};

/*  Get the barycentric coordinates for a point on a triangle  */
/*    return 1/0 if p lies inside/outside triangle spanned by a,b,c  */
int clUtils_namespace::get_barycentric_coords(const double* const a,const double* const b,const double* const c,
	const double* const p,double &alpha,double &beta,double &gamma)
{
	const double eps=1.e-14;
	const double upp=1.0+eps;

	double ab[3],ac[3],abxac[3];
	ab[0]=b[0]-a[0];
	ab[1]=b[1]-a[1];
	ab[2]=b[2]-a[2];
	ac[0]=c[0]-a[0];
	ac[1]=c[1]-a[1];
	ac[2]=c[2]-a[2];
	clUtils_namespace::outerproduct(ab,ac,abxac);
	const double norm_abxac=abxac[0]*abxac[0]+abxac[1]*abxac[1]+abxac[2]*abxac[2];
	if (!norm_abxac) {
		alpha=1.0;
		beta=0.0;
		gamma=0.0;
		return 1;
	} else {
		double pb[3],pc[3],pbxpc[3];
		pb[0]=b[0]-p[0];
		pb[1]=b[1]-p[1];
		pb[2]=b[2]-p[2];
		pc[0]=c[0]-p[0];
		pc[1]=c[1]-p[1];
		pc[2]=c[2]-p[2];
		clUtils_namespace::outerproduct(pb,pc,pbxpc);
		alpha=(pbxpc[0]*abxac[0]+pbxpc[1]*abxac[1]+pbxpc[2]*abxac[2])/norm_abxac;
		double ap[3],apxac[3];
		ap[0]=p[0]-a[0];
		ap[1]=p[1]-a[1];
		ap[2]=p[2]-a[2];
		clUtils_namespace::outerproduct(ap,ac,apxac);
		beta=(apxac[0]*abxac[0]+apxac[1]*abxac[1]+apxac[2]*abxac[2])/norm_abxac;
		gamma=1.0-alpha-beta;
#ifdef DEBUG
		double abxap[3],gamma_check;
		clUtils_namespace::outerproduct(ab,ap,abxap);
		gamma_check=(abxap[0]*abxac[0]+abxap[1]*abxac[1]+abxap[2]*abxac[2])/norm_abxac;
		if (!within_eps(gamma,gamma_check,1.e-6))
			if (gamma>1.e-6) fatal_err("Inconsistent gamma in clUtils_namespace::get_barycentric_coords");
#endif
	}
	if (alpha<-eps || alpha>upp) return 0;
	if (beta<-eps || beta>upp) return 0;
	if (gamma<-eps || gamma>upp) return 0;
	return 1;
};

/*  Convert a character string to all lowercase  */
char* clUtils_namespace::lowercase(char* const s)
{
	unsigned char c;
	for (unsigned int i=0; i<strlen(s); i++) {
		c=s[i];
		if (c>64 && c<91) {
			c+=32;
			s[i]=(char)c;
		}
	}
	return s;
};

/*  Convert a character string to all uppercase  */
char* clUtils_namespace::uppercase(char* const s)
{
	unsigned char c;
	for (unsigned int i=0; i<strlen(s); i++) {
		c=s[i];
		if (c>96 && c<123) {
			c-=32;
			s[i]=(char)c;
		}
	}
	return s;
};

/*  Allocates a new character string of size len */
char* clUtils_namespace::new_cstring(char **s,const int len)
{
	*s=new char [len+1];
	if (!(*s)) fatal_err(0,"*s in clUtils_namespace::new_cstring");
	(*s)[0]='\0';
	return *s;
};
/*  Allocates and initializes a character string  */
char* clUtils_namespace::new_cstring(char **s,const char* const init)
{
	*s=new char [strlen(init)+1];
	if (!(*s)) fatal_err(0,"*s in clUtils_namespace::new_cstring");
	strcpy(*s,init);
	return *s;
};
/*  Deallocates a character string  */
void clUtils_namespace::delete_cstring(char *s)
{
	if (strlen(s)>0) delete [] s;
};

/*  Read a string enclosed in double quotes from an input stream  */
/*  Start reading after the first double quote and keep reading until the second one  */
/*  Interpret \" as an embedded double quote and include in string  */
/*  Return EOF for a maximum string length error or an embedded new line  */
/*  Return EOF for a read error or an embedded EOF  */
/*  Return string length after successful completion  */
int clUtils_namespace::readdqstring(FILE *fileptr,char *s,const int lenmax)
{
	int c,index;
	do {
		if ((c=getc(fileptr))==EOF) return c;
	} while (c!='"');
	for (index=0,c=getc(fileptr); c!='"'; c=getc(fileptr)) {
		if (c==EOF) return EOF;
		if (c=='\n') return EOF;
		if (c=='\\') {
			c=getc(fileptr);
			if (c==EOF) return EOF;
			if (c=='\n') return EOF;
			if (c!='"') {
				s[index++]='\\';
				if (index==lenmax) return EOF;
			}
		}
		s[index++]=char(c);
		if (index==lenmax) return EOF;
	}
	s[index]='\0';
	return index;
};

/*  Read a "separator" separated string from a char string  */
/*  Return -1 for reading past end of string */
/*  Return -1 for a maximum string length error  */
/*  Return substring length after successful completion  */
int clUtils_namespace::splitstring(const char* const charstring,const char separator,char* const s,const int lenmax)
{
	int startindex=0;
	return splitstring(charstring,separator,s,lenmax,startindex);
}
int clUtils_namespace::splitstring(const char* const charstring,const char separator,char* const s,const int lenmax,int &startindex)
{
	int len=strlen(charstring)-startindex;
	if (len<=0) {
		s[0]='\0';
		return -1;
	}

	const char *ptr=strchr(charstring+startindex,separator);
	if (ptr) len=ptr-charstring-startindex;
	if (len>=lenmax) return -1;
	strncpy(s,charstring+startindex,len);
	s[len]='\0';
	startindex+=len+1;

	return len;
};


/*  Read a comma or newline limited string from an input stream  */
/*  Return -1 for a maximum string length error  */
/*  Return -1 for an embedded EOF  */
/*  Return string length after successful completion  */
int clUtils_namespace::readcsvstring(FILE *fileptr,char *s,const int lenmax)
{
	int index;
	int c=getc(fileptr);
	for (index=0; c!=','; c=getc(fileptr)) {
		if (c==EOF) return -1;
		if (c=='\n') break;
		s[index++]=char(c);
		if (index==lenmax) return -1;
	}
	s[index]='\0';
	return index;
};

/*  Read a comma or newline limited integer from an input stream  */
/*  Return -1 for a maximum string length error  */
/*  Return -1 for an embedded EOF or an illegal character  */
/*  Return 0 after successful completion  */
int clUtils_namespace::readcsvint(FILE *fileptr,int &ivalue)
{
	return readcsvint(fileptr,&ivalue);
};
int clUtils_namespace::readcsvint(FILE *fileptr,int *ivalue)
{
	char s[LINE_LENGTH];
	int index;
	int c=getc(fileptr);
	for (index=0; c!=','; c=getc(fileptr)) {
		if (c==EOF) return -1;
		if (c=='\n') break;
		if (c<'0' || c>'9') {
			if (index) return -1;
			else if (!(c=='+' || c=='-')) return -1;
		}
		s[index++]=char(c);
		if (index==LINE_LENGTH) return -1;
	}
	s[index]='\0';
	*ivalue=atoi(s);
	return 0;
};

/*  Read a comma or newline limited double from an input stream  */
/*  Return -1 for a maximum string length error  */
/*  Return -1 for an embedded EOF or an illegal character  */
/*  Return 0 after successful completion  */
int clUtils_namespace::readcsvdouble(FILE *fileptr,double &dvalue)
{
	return readcsvdouble(fileptr,&dvalue);
};
int clUtils_namespace::readcsvdouble(FILE *fileptr,double *dvalue)
{
	char s[LINE_LENGTH];
	int index;
	int c=getc(fileptr);
	for (index=0; c!=','; c=getc(fileptr)) {
		if (c==EOF) return -1;
		if (c=='\n') break;
		if (c<'0' || c>'9') {
			if (!(c=='+' || c=='-' || c=='.' || c=='e' || c=='E')) return -1;
		}
		s[index++]=char(c);
		if (index==LINE_LENGTH) return -1;
	}
	s[index]='\0';
	*dvalue=atof(s);
	return 0;
}

/*  Read an uninterupted string of letters from an input stream  */
/*  Start reading after the first letter and keep reading  */
/*    letters until any non-letter character  */
/*  Return EOF for an EOF or read error or 0 for a maximum string length error  */
/*  Return string length after successful completion  */
int clUtils_namespace::readletterstring(FILE *fileptr,char *s,const int lenmax)
{
	int c;
	for (;;) {
		if ((c=getc(fileptr))==EOF) return EOF;
		if ((c>='a' && c<='z') || (c>='A' && c<='Z')) {
			ungetc(c,fileptr);
			break;
		}
	}
	int index=0;
	for (;;) {
		if ((c=getc(fileptr))==EOF) return EOF;
		if ((c>='a' && c<='z') || (c>='A' && c<='Z')) {
			s[index++]=char(c);
			if (index==lenmax) return 0;
		} else {
			ungetc(c,fileptr);
			break;
		}
	}
	s[index]='\0';
	return index;
};

/*  Read an un-interupted string of letters or numbers from an input stream  */
/*  Start reading after the first letter or number and keep reading  */
/*    until any non-letter/number character  */
/*  Return EOF for an EOF or read error or 0 for a maximum string length error  */
/*  Return string length after successful completion  */
int clUtils_namespace::readletternumberstring(FILE *fileptr,char *s,const int lenmax)
{
	int c;
	for (;;) {
		if ((c=getc(fileptr))==EOF) return EOF;
		if ((c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='0' && c<='9')) {
			ungetc(c,fileptr);
			break;
		}
	}
	int index=0;
	for (;;) {
		if ((c=getc(fileptr))==EOF) return EOF;
		if ((c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='0' && c<='9')) {
			s[index++]=char(c);
			if (index==lenmax) return 0;
		} else {
			ungetc(c,fileptr);
			break;
		}
	}
	s[index]='\0';
	return index;
};

/*  Read an input stream until a match with a certain string is found  */
/*  Return 0 for a read error or EOF or 1 for successful completion  */
/*  NOTE that this routine is case sensitive  */
int clUtils_namespace::readuntilstring(FILE *fileptr,const char* const matchstring)
{
	char readstring[LINE_LENGTH];
	for (;;) {
		if (fscanf(fileptr,"%s",readstring)!=1) return 0;
		if (!strcmp(readstring,matchstring)) return 1;
	}
};

/*  Read an input stream until a match with a certain string is found  */
/*  Skip the remainder of any line after the skip_char character  */
/*  Return 0 for a read error or EOF or 1 for successful completion  */
/*  NOTE that this routine is case sensitive  */
int clUtils_namespace::readuntilstring(FILE *fileptr,const char* const matchstring,const char skip_char)
{
	char readstring[LINE_LENGTH];
	for (;;) {
		if (fscanf(fileptr,"%s",readstring)!=1) return 0;
		if (readstring[0]==skip_char) fgets(readstring,LINE_LENGTH,fileptr);
		else if (!strcmp(readstring,matchstring)) return 1;
	}
};

/*  Read an input stream until a match with a certain string is found  */
/*  Return 0 for a read error or EOF or 1 for successful completion  */
/*  NOTE that this routine is case insensitive  */
int clUtils_namespace::readuntilstring_nocase(FILE *fileptr,const char* const matchstring)
{
	char readstring[LINE_LENGTH],mmatchstring[LINE_LENGTH];
	strcpy(mmatchstring,matchstring);
	for (;;) {
		if (fscanf(fileptr,"%s",readstring)!=1) return 0;
		if (!strcmp(uppercase(readstring),uppercase(mmatchstring))) return 1;
	}
};

/*  Read an input stream until a match with a certain string is found  */
/*  Skip the remainder of any line after the skip_char character  */
/*  Return 0 for a read error or EOF or 1 for successful completion  */
/*  NOTE that this routine is case insensitive  */
int clUtils_namespace::readuntilstring_nocase(FILE *fileptr,const char* const matchstring,const char skip_char)
{
	char readstring[LINE_LENGTH],mmatchstring[LINE_LENGTH];
	strcpy(mmatchstring,matchstring);
	for (;;) {
		if (fscanf(fileptr,"%s",readstring)!=1) return 0;
		if (readstring[0]==skip_char) fgets(readstring,LINE_LENGTH,fileptr);
		else if (!strcmp(uppercase(readstring),uppercase(mmatchstring))) return 1;
	}
};

/*  Check if a file exists  */
int clUtils_namespace::check_file(const char* const filename)
{
	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) {
		return 0;
	} else {
		fclose(fileptr);
		return 1;
	}
}

/*  Get the first line number in a file that contains a string  */
/*    return -1 if the file does not exists  */
/*    return -2 if the string is not found  */
int clUtils_namespace::get_num_line(const char* const filename,const char* const findstring)
{
	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) {
		return -1;
	} else {
		int c,num_line=0;
		char line[LINE_LENGTH];
		for (;;) {
			c=getc(fileptr);
			if (c==EOF) break;
			ungetc(c,fileptr);
			fgets(line,LINE_LENGTH,fileptr);
			num_line++;
			if (strstr(line,findstring)) {
				fclose(fileptr);
				return num_line;
			}
		}
		fclose(fileptr);
		return num_line;
	}
}

/*  Get the number of lines in a file (return -1 if the file does not exists  */
int clUtils_namespace::get_num_line(const char* const filename,const bool comm_skip)
{
	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) {
		return -1;
	} else {
		int c,num_line=0;
		char line[LINE_LENGTH];
		for (;;) {
			c=getc(fileptr);
			if (c==EOF) break;
			ungetc(c,fileptr);
			fgets(line,LINE_LENGTH,fileptr);
			if (comm_skip && line[0]=='#') continue;
			num_line++;
		}
		fclose(fileptr);
		return num_line;
	}
}

/*  Read an integer from an ascii input stream  */
int clUtils_namespace::read_int(FILE *fileptr)
{
	int readvalue;
	if (fscanf(fileptr,"%d",&readvalue)!=1) fatal_err("read error on clUtils_namespace::read_int");
	return readvalue;
};
int clUtils_namespace::read_int(FILE *fileptr,int &readvalue)
{
	if (fscanf(fileptr,"%d",&readvalue)!=1) return 0;
	return 1;
};
int clUtils_namespace::read_int(FILE *fileptr,int* const readvalue,const int num)
{
	int *readvalueptr=readvalue;
	for (int n=0; n<num; n++) 
		if (fscanf(fileptr,"%d",readvalueptr++)!=1) return 0;
	return num;
};

/*  Read a double from an ascii input stream  */
double clUtils_namespace::read_double(FILE *fileptr)
{
	double readvalue;
	if (fscanf(fileptr,"%le",&readvalue)!=1) fatal_err("read error on clUtils_namespace::read_double");
	return readvalue;
};
int clUtils_namespace::read_double(FILE *fileptr,double &readvalue)
{
	if (fscanf(fileptr,"%le",&readvalue)!=1) return 0;
	return 1;
};
int clUtils_namespace::read_double(FILE *fileptr,double* const readvalue,const int num)
{
	double *readvalueptr=readvalue;
	for (int n=0; n<num; n++) 
		if (fscanf(fileptr,"%le",readvalueptr++)!=1) return 0;
	return num;
};

/*  Read a char string from an ascii input stream  */
int clUtils_namespace::read_char(FILE *fileptr,char *readvalue)
{
	if (fscanf(fileptr,"%s",readvalue)!=1) return 0;
	return 1;
};

/*  Skip a given number of entries from an ascii input stream  */
int clUtils_namespace::skip_read(FILE *fileptr,const int num_skip)
{
	char readstring[LINE_LENGTH];
	for (int n_skip=0; n_skip<num_skip; n_skip++)
		if (fscanf(fileptr,"%s",readstring)!=1) return 0;
	return 1;
};

/*  Read an integer from a binary input stream  */
int clUtils_namespace::read_binentry_int(FILE *fileptr,int &readvalue)
{
	my_int readint;
	unsigned size=sizeof(readint);
	if (_little_endian_io) {
		void *ioptr=(void *)(&readint);
		if (fread(ioptr,size,1,fileptr)!=1) return 0;
	} else {
		readint=0;
		int readchar;
		for (; size; --size) {
			if ((readchar=fgetc(fileptr))==EOF) return 0;
			readint=(readint<<8)|readchar;
		}
// If needed at a later stage: to read little endian for internal big endian use:
//				for (unsigned size=0; size<sizeof(readvalue); ++size) {
//					if ((readchar=fgetc(fileptr))==EOF) return 1;
//					readvalue|=readchar<<(8*size);			
//				}
	}
	readvalue=(int)(readint);
	return 1;
};
/*  Read an integer array from a binary input stream  */
int clUtils_namespace::read_binentry_int(FILE *fileptr,int *readvalue,int num)
{
	for (int n=0; n<num; n++) {
		if (!read_binentry_int(fileptr,readvalue[n])) return 0;
	}
	return 1;
};
/*  Read an integer from a binary input stream and check endian  */
int clUtils_namespace::read_binentry_int(FILE *fileptr,int &readvalue,int checkvalue)
{
	_little_endian_io=1;

	my_int readint;
	unsigned size=sizeof(readint);
	if (_little_endian_io) {
		void *ioptr=(void *)(&readint);
		if (fread(ioptr,size,1,fileptr)!=1) return 0;
	} else {
		readint=0;
		int readchar;
		for (; size; --size) {
			if ((readchar=fgetc(fileptr))==EOF) return 0;
			readint=(readint<<8)|readchar;
		}
// If needed at a later stage: to read little endian for internal big endian use:
//				for (unsigned size=0; size<sizeof(readvalue); ++size) {
//					if ((readchar=fgetc(fileptr))==EOF) return 1;
//					readvalue|=readchar<<(8*size);			
//				}
	}
	readvalue=(int)(readint);
	if (readvalue!=checkvalue) {
		big_to_little_endian_int(readint);
		readvalue=(int)(readint);
		if (readvalue==checkvalue) _little_endian_io=0;
		else return -1;
	}

	return 1;
};

/*  Convert an integer from big to little endian  */
void clUtils_namespace::big_to_little_endian_int(my_int &readint)
{
	my_int value=0;
	int size=0,maxsize=8*sizeof(my_int);
	for (; size<maxsize; size+=8) {
		value=((readint>>size)&0xFF)|(value<<8);
	}
	readint=value;
};

/*  Check if an integer is a power of 2  */
/*    return the power if it is or 0 if not  */
int clUtils_namespace::is_power_two(int n)
{
	int pow=1;
	while (n>2 && !(n%2)) {
		n>>=1;
		pow++;
	}
	if (n!=2) {
//		fprintf(stdout,"not a power of two\n");
		return 0;
	} else {
//		fprintf(stdout,"%d to the power %d = %d\n",n,pow,power(2,pow));
		return pow;
	}
};

/*  Get the nearest smaller integer that is a power of 2  */
/*    return the largest value nn such that nn=2**pow<=n  */
int clUtils_namespace::nearest_power_two(int n)
{
	int nn=1;
	while (nn<=n) {
		nn<<=1;
	}
	if (nn>n) nn>>=1;
	return nn;
};

/*  Remove leading and trailing edge white space from char  */
/*  Also removes trailing \n  */
int clUtils_namespace::remove_sur_white(const char* const charstring, char* const s)
{
	int iskip1=0;
	int iskip2=strlen(charstring)-1; // -1 removes \n
	while (charstring[iskip1]<33) iskip1++;  // Assumes ASCII table values
	while (charstring[iskip2-1]<33) iskip2--;  // Assumes ASCII table values
	int len=iskip2-iskip1;
	strncpy(s,charstring+iskip1,len);
	s[len]='\0';
	return len;
};
