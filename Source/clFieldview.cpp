/*******************************************************************************
  A collection of Fieldview usage functions
  Author: Rolf Verberg
*******************************************************************************/

#include "clFieldview.h"

using clUtils_namespace::warning;
using clUtils_namespace::fatal_err;
using clUtils_namespace::readdqstring;
using clUtils_namespace::readletterstring;
using clUtils_namespace::uppercase;
using clUtils_namespace::matrix;

const int clFieldview::IO_LINE_LENGTH=80;

const int clFieldview::FV_MAGIC=0x00010203;	/* decimal 66051 */

const int clFieldview::FV_GRIDS_FILE=1;
const int clFieldview::FV_RESULTS_FILE=2;
const int clFieldview::FV_COMBINED_FILE=3;

const int clFieldview::FV_NODES=1001;
const int clFieldview::FV_FACES=1002;
const int clFieldview::FV_ELEMENTS=1003;
const int clFieldview::FV_VARIABLES=1004;
const int clFieldview::FV_BNDRY_VARS=1006;
const int clFieldview::FV_ARB_POLY_FACES=1007;
const int clFieldview::FV_ARB_POLY_ELEMENTS=1008;
const int clFieldview::FV_ARB_POLY_BNDRY_VARS=1009;

const int clFieldview::FV_TET_ELEM_ID=1;
const int clFieldview::FV_HEX_ELEM_ID=2;
const int clFieldview::FV_PRISM_ELEM_ID=3;
const int clFieldview::FV_PYRA_ELEM_ID=4;
const int clFieldview::FV_ARB_POLY_ELEM_ID=5;

static int read_bin_debug=0;

/* fread_str80: read in a string padded to 80 characters */
int clFieldview::fread_str80(FILE *file_ptr,char *string)
{
	strcpy(string,"");
	void *io_ptr=(void *)string;
	if (fread(io_ptr,80*sizeof(char),1,file_ptr)!=1) return 0;
	/* Remove trailing spaces */
	for (int i=IO_LINE_LENGTH-1; i>=0; i--) {
		if (string[i]<'!' || string[i]>'~') string[i]=0;
	}
	return 1;
}

void clFieldview::read(const char *filename)
{
	using clUtils_namespace::read_binentry;
	using clUtils_namespace::read_binarray;
	using clUtils_namespace::read_binentry_int;

	if (clUtils_namespace::_size_my_int!=4) fatal_err("Incompatible int32 type");
	if (clUtils_namespace::_size_my_float!=4) fatal_err("Incompatible float type");

	int readint;
	char *readstring=new char [IO_LINE_LENGTH];

	FILE *file_ptr=fopen(filename,"rb");
	if (!file_ptr) fatal_err(1,filename);

	/*  Magic number  */
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"Magic number=%d %d\n",readint,FV_MAGIC);
	if (readint!=FV_MAGIC) fatal_err("Illegal magic number");

	if (!fread_str80(file_ptr,readstring)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"%s\n",readstring);

	/* Version */
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	_version=(double)(readint);
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	_version+=0.1*(double)(readint);
	if (read_bin_debug) fprintf(stdout,"%.1f\n",_version);

	/* File type */
	if (!read_binentry_int(file_ptr,_filetype)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"File type code=%d\n",_filetype);
	if (_filetype!=FV_GRIDS_FILE && _filetype!=FV_RESULTS_FILE && _filetype!=FV_COMBINED_FILE) 
		fatal_err("Input file not a grid, a result or combined file");

	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);

	/* Read header */
	read_header(file_ptr,filename);

	/* Read grids */
	int n_grid;
	clFieldview_grid *currentgrid=NULL,*newgrid=NULL;
	for (n_grid=0; n_grid<_num_grid; n_grid++) {
		newgrid=new clFieldview_grid(_filetype,_num_facetype,_num_var,_num_bvar);
		newgrid->read_grid(file_ptr,filename);
		if (!n_grid) _firstgrid=newgrid;
		else currentgrid->_next=newgrid;
		currentgrid=newgrid;
	}

	/* Consolidate all grids */
	_grids=new clFieldview_grid* [_num_grid];
	if (!_grids) fatal_err(0,"_grids in clFieldview::read");
	currentgrid=_firstgrid;
	for (n_grid=0; n_grid<_num_grid; n_grid++) {
		if (!currentgrid) fatal_err("currentgrid not allocated in clFieldview::read");
		_grids[n_grid]=currentgrid;
		currentgrid=currentgrid->_next;
	}

	fclose(file_ptr);

	delete [] readstring;
};

void clFieldview::read_header(FILE *file_ptr,const char *filename)
{
	using clUtils_namespace::read_binentry;
	using clUtils_namespace::read_binentry_int;

	int i;
	clUtils_namespace::my_float readfloat;
	char *readstring=new char [IO_LINE_LENGTH];

	if (_filetype!=FV_GRIDS_FILE) {
		if (!read_binentry(file_ptr,readfloat)) fatal_err(2,filename);
		_time=(double)(readfloat);
		if (read_bin_debug) fprintf(stdout,"TIME=%le\n",_time);
		if (!read_binentry(file_ptr,readfloat)) fatal_err(2,filename);
		_mach=(double)(readfloat);
		if (read_bin_debug) fprintf(stdout,"FSMACH=%le\n",_mach);
		if (!read_binentry(file_ptr,readfloat)) fatal_err(2,filename);
		_alpha=(double)(readfloat);
		if (read_bin_debug) fprintf(stdout,"ALPHA=%le\n",_alpha);
		if (!read_binentry(file_ptr,readfloat)) fatal_err(2,filename);
		_re=(double)(readfloat);
		if (read_bin_debug) fprintf(stdout,"RE=%le\n",_re);
	}

	if (!read_binentry_int(file_ptr,_num_grid)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"Number of grids=%d\n",_num_grid);

	if (_filetype!=FV_RESULTS_FILE) {
		if (!read_binentry_int(file_ptr,_num_facetype)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Number of boundary types=%d\n",_num_facetype);
		_faceresultflag=new int [_num_facetype];
		if (!_faceresultflag) fatal_err(0,"_faceresultflag in clFieldview::read_header");
		_facenormalflag=new int [_num_facetype];
		if (!_facenormalflag) fatal_err(0,"_facenormalflag in clFieldview::read_header");
		_facetypename=new char* [_num_facetype];
		if (!_facetypename) fatal_err(0,"_facetypename in clFieldview::read_header");
		for (i=0; i<_num_facetype; i++) {
			if (!read_binentry_int(file_ptr,_faceresultflag[i])) fatal_err(2,filename);
			if (read_bin_debug) fprintf(stdout,"  %d",_faceresultflag[i]);
			if (!read_binentry_int(file_ptr,_facenormalflag[i])) fatal_err(2,filename);
			if (read_bin_debug) fprintf(stdout,"  %d",_facenormalflag[i]);
			if (!fread_str80(file_ptr,readstring)) fatal_err(2,filename);
			_facetypename[i]=new char [strlen(readstring)+1];
			if (!_facetypename[i]) fatal_err(0,"_facetypename[i] in clFieldview::read_header");
			strcpy(_facetypename[i],readstring);
			if (read_bin_debug) fprintf(stdout,"  %s\n",_facetypename[i]);
		}
	}

	if (_filetype!=FV_GRIDS_FILE) {
		if (!read_binentry_int(file_ptr,_num_var)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Number of function variables=%d\n",_num_var);
		_varname=new char* [_num_var];
		if (!_varname) fatal_err(0,"_varname in clFieldview::read_header");
		for (i=0; i<_num_var; i++) {
			if (!fread_str80(file_ptr,readstring)) fatal_err(2,filename);
			_varname[i]=new char [strlen(readstring)+1];
			if (!_varname[i]) fatal_err(0,"_varname[i] in clFieldview::read_header");
			strcpy(_varname[i],readstring);
			if (read_bin_debug) fprintf(stdout,"  %s\n",_varname[i]);
		}
		if (!read_binentry_int(file_ptr,_num_bvar)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Number of boundary function variables=%d\n",_num_bvar);
		_bvarname=new char* [_num_bvar];
		if (!_bvarname) fatal_err(0,"_bvarname in clFieldview::read_header");
		for (i=0; i<_num_bvar; i++) {
			if (!fread_str80(file_ptr,readstring)) fatal_err(2,filename);
			_bvarname[i]=new char [strlen(readstring)+1];
			if (!_bvarname[i]) fatal_err(0,"_bvarname[i] in clFieldview::read_header");
			strcpy(_bvarname[i],readstring);
			if (read_bin_debug) fprintf(stdout,"  %s\n",_bvarname[i]);
		}
	}

	delete [] readstring;
};

void clFieldview_grid::read_grid(FILE *file_ptr,const char *filename)
{
	using clUtils_namespace::read_binarray;
	using clUtils_namespace::read_binentry_int;

	int i,readint,sectionheader;

	/* Read header */
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	if (readint!=clFieldview::FV_NODES) fatal_err("Illegal nodes header");
	if (!read_binentry_int(file_ptr,_num_node)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"Number of nodes=%d\n",_num_node);

	if (_filetype==clFieldview::FV_RESULTS_FILE) {

		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		sectionheader=readint;

	} else {

		/* Read coordinates */
		_xcoord=new clUtils_namespace::my_float [_num_node];
		if (!_xcoord) fatal_err(0,"_xcoord in clFieldview_grid::read_grid");
		_ycoord=new clUtils_namespace::my_float [_num_node];
		if (!_ycoord) fatal_err(0,"_ycoord in clFieldview_grid::read_grid");
		_zcoord=new clUtils_namespace::my_float [_num_node];
		if (!_zcoord) fatal_err(0,"_zcoord in clFieldview_grid::read_grid");
		if (!read_binarray(file_ptr,_xcoord,_num_node)) fatal_err(2,filename);
		if (!read_binarray(file_ptr,_ycoord,_num_node)) fatal_err(2,filename);
		if (!read_binarray(file_ptr,_zcoord,_num_node)) fatal_err(2,filename);

		/* Read boundary faces */
		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		sectionheader=readint;
		for (;;) {
			if (sectionheader==clFieldview::FV_FACES) 
				read_standard_face(file_ptr,filename);
			else if (sectionheader==clFieldview::FV_ARB_POLY_FACES) 
				read_arbitrary_face(file_ptr,filename);
			else break;
			if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
			sectionheader=readint;
		}

		/* Read volume elements */
		for (;;) {
			if (sectionheader==clFieldview::FV_ELEMENTS)
				read_standard_element(file_ptr,filename);
			else if (sectionheader==clFieldview::FV_ARB_POLY_ELEMENTS)
				fatal_err("Arbitrary polyhedron elements not implemented");
			else break;
			if (!read_binentry_int(file_ptr,readint)) {
				if (_filetype==clFieldview::FV_GRIDS_FILE) break;
				else fatal_err(2,filename);
			}
			sectionheader=readint;
		}
	}

	if (_filetype!=clFieldview::FV_GRIDS_FILE) {

		/* Read function values */
		if (sectionheader!=clFieldview::FV_VARIABLES) fatal_err("Illegal variables header");
		_func=matrix(_func,_num_var,_num_node);
		if (!_func) fatal_err(0,"_func in clFieldview_grid::read_grid");
		for (i=0; i<_num_var; i++)
			if (!read_binarray(file_ptr,_func[i],_num_node)) fatal_err(2,filename);

		/* Read boundary values */
		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		sectionheader=readint;
		if (sectionheader!=clFieldview::FV_BNDRY_VARS) fatal_err("Illegal boundary variables header");
		if (_num_bvar) fatal_err("boundary variables not implemented yeat");
/*		_bfunc=matrix(_bfunc,_num_bvar,_num_node);
		if (!_func) fatal_err(0,"_func in clFieldview_grid::read_grid");
		for (i=0; i<_num_bvar; i++) {
			fatal_err("Boundary variables not implemented");
			if (!read_binarray(file_ptr,_bfunc[i],_num_node)) fatal_err(2,filename);
		}

		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		sectionheader=readint;
		if (sectionheader!=clFieldview::FV_ARB_POLY_BNDRY_VARS) fatal_err("Illegal boundary variables header");
		for (i=0; i<_num_bvar; i++) {
			fatal_err("Boundary variables not implemented");
			if (!read_binarray(file_ptr,func[i],num_node)) fatal_err(2,filename);
		}*/
	}
};

void clFieldview_grid::read_standard_face(FILE *file_ptr,const char *filename)
{
	using clUtils_namespace::read_binentry_int;

	int n_facetype,num_face;
	if (!read_binentry_int(file_ptr,n_facetype)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"Boundary type=%d\n",n_facetype);
	n_facetype--;
	if (n_facetype<0 || n_facetype>=_num_facetype) fatal_err("Illegal boundary type");
	if (!read_binentry_int(file_ptr,num_face)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"  number of faces=%d\n",num_face);
	int **conn=NULL,*connptr=NULL;
	conn=matrix(conn,num_face,4);
	if (!conn) fatal_err(0,"conn in clFieldview_grid::read_standard_face");
	connptr=conn[0];
	if (!read_binentry_int(file_ptr,connptr,4*num_face)) fatal_err(2,filename);
	standard_face *faces=new standard_face;
	faces->_num_face=num_face;
	faces->_conn=conn;
	faces->_next=NULL;
	if (!_standard_face[n_facetype]) {
		_standard_face[n_facetype]=faces;
	} else {
		standard_face *standard_face_ptr=_standard_face[n_facetype];
		while (standard_face_ptr->_next) standard_face_ptr=standard_face_ptr->_next;
		standard_face_ptr->_next=faces;
	}
};

void clFieldview_grid::read_arbitrary_face(FILE *file_ptr,const char *filename)
{
	using clUtils_namespace::read_binentry_int;

	fatal_err("clFieldview_grid::read_arbitrary_face not implemented");
	int n_facetype,n_face,num_face,num_facenode=0;
	if (!read_binentry_int(file_ptr,n_facetype)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"Boundary type=%d\n",n_facetype);
	n_facetype--;
	if (n_facetype<0 || n_facetype>=_num_facetype) fatal_err("Illegal boundary type");
	if (!read_binentry_int(file_ptr,num_face)) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"  number of faces=%d\n",num_face);
	for (n_face=0; n_face<num_face; n_face++) {
		int *face=new int [num_facenode];
		if (!read_binentry_int(file_ptr,face,num_facenode)) fatal_err(2,filename);
	}
};

/* RV: This isn't implemented yet: elements won't be stored */
void clFieldview_grid::read_standard_element(FILE *file_ptr,const char *filename)
{
	using clUtils_namespace::read_binentry_int;

	int i,j,readint;

	int num_tet,num_hex,num_prism,num_pyramid;
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	num_tet=readint;
	if (read_bin_debug) fprintf(stdout,"Number of tetrahedrons=%d\n",num_tet);
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	num_hex=readint;
	if (read_bin_debug) fprintf(stdout,"Number of hexahedrons=%d\n",num_hex);
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	num_prism=readint;
	if (read_bin_debug) fprintf(stdout,"Number of prisms=%d\n",num_prism);
	if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
	num_pyramid=readint;
	if (read_bin_debug) fprintf(stdout,"Number of pyramids=%d\n",num_pyramid);
	for (i=0; i<num_tet; i++) {
		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Element header=%d\n",readint);
		if (read_bin_debug) fprintf(stdout,"  Connectivities:");
		for (j=0; j<4; j++) {
			if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
			if (read_bin_debug) fprintf(stdout," %d",readint);
		}
		if (read_bin_debug) fprintf(stdout,"\n");
	}
	for (i=0; i<num_hex; i++) {
		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Element header=%d\n",readint);
		if (read_bin_debug) fprintf(stdout,"  Connectivities:");
		for (j=0; j<8; j++) {
			if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
			if (read_bin_debug) fprintf(stdout," %d",readint);
		}
		if (read_bin_debug) fprintf(stdout,"\n");
	}
	for (i=0; i<num_prism; i++) {
		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Element header=%d\n",readint);
		if (read_bin_debug) fprintf(stdout,"  Connectivities:");
		for (j=0; j<6; j++) {
			if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
			if (read_bin_debug) fprintf(stdout," %d",readint);
		}
		if (read_bin_debug) fprintf(stdout,"\n");
	}
	for (i=0; i<num_pyramid; i++) {
		if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Element header=%d\n",readint);
		if (read_bin_debug) fprintf(stdout,"  Connectivities:");
		for (j=0; j<5; j++) {
			if (!read_binentry_int(file_ptr,readint)) fatal_err(2,filename);
			if (read_bin_debug) fprintf(stdout," %d",readint);
		}
		if (read_bin_debug) fprintf(stdout,"\n");
	}
};

void clFieldview::convert_to_tecplot_surface(clTecplot *tecplot)
{
	int num_conn,n_data,num_data,n_elem,num_elem,n_funcvar,n_facetype,n_grid,n_node,num_node;
	int nshare,n_var,num_var,num_funcvar,n_zone,nn_zone,num_zone;
	int *connptr=NULL,*fieldviewconnptr=NULL;
	double **data=NULL,*dataptr=NULL;
	clVector_namespace::clMatrix<int> conn;
	clTecplot_zone *tecplotzone=NULL;
	clUtils_namespace::my_float *fieldviewdataptr=NULL;
	clFieldview_grid *currentgrid=NULL;
	clFieldview_grid::standard_face *standard_face_ptr=NULL;

	if (_filetype==clFieldview::FV_RESULTS_FILE) {
		num_zone=_num_grid;
	} else {
		num_zone=0;
		for (n_grid=0; n_grid<_num_grid; n_grid++) {
			currentgrid=_grids[n_grid];
			if (!currentgrid) 
				fatal_err("currentgrid not allocated in clFieldview::convert_to_tecplot_surface");
			for (n_facetype=0; n_facetype<_num_facetype; n_facetype++) {
				standard_face_ptr=currentgrid->_standard_face[n_facetype];
				for (num_elem=0; ; standard_face_ptr=standard_face_ptr->_next) {
					num_elem+=standard_face_ptr->_num_face;
					if (!standard_face_ptr->_next) break;
				}
				if (num_elem) num_zone++;
			}
		}
	}
	if (read_bin_debug) fprintf(stdout,"Number of Tecplot zones=%d\n",num_zone);

	num_funcvar=_num_var;
	num_var=num_funcvar;
	if (_filetype!=clFieldview::FV_RESULTS_FILE) num_var+=3;

	tecplot->alloc_zones(num_zone);
	tecplot->set_num_var(num_var);
	tecplot->alloc_varname();
	n_var=0;
	if (_filetype!=clFieldview::FV_RESULTS_FILE) {
		tecplot->set_varname(n_var++,"X");
		tecplot->set_varname(n_var++,"Y");
		tecplot->set_varname(n_var++,"Z");
	}
	for (n_funcvar=0; n_var<num_var; n_var++,n_funcvar++) {
		tecplot->set_varname(n_var,_varname[n_funcvar]);
	}

	n_zone=0;
	for (n_grid=0; n_grid<_num_grid; n_grid++) {

		currentgrid=_grids[n_grid];
		if (!currentgrid) 
			fatal_err("currentgrid not allocated in clFieldview::convert_to_tecplot_surface");

		num_node=currentgrid->_num_node;

		if (_filetype==clFieldview::FV_RESULTS_FILE) {

			tecplotzone=new clTecplot_zone_unstructured(num_var);
			if (!tecplotzone) fatal_err(0,"tecplotzone in clFieldview::convert_to_tecplot_surface");
			tecplotzone->set_zonetype_num_conn("FEQUADRILATERAL");
			tecplotzone->set_datatype("BLOCK");
			tecplotzone->set_num_node(num_node);
			tecplotzone->allocate_data();
			if (num_funcvar) {
				fieldviewdataptr=currentgrid->_func[0];
				for (n_funcvar=0; n_funcvar<num_funcvar; n_funcvar++) {
					dataptr=tecplotzone->get_datastartptr(n_funcvar);
					for (n_node=0; n_node<num_node; n_node++) (*dataptr++)=(double)(*fieldviewdataptr++);
				}
			}

			tecplot->set_zone(n_zone,tecplotzone);

			n_zone++;

		} else {

			nn_zone=0;
			for (n_facetype=0; n_facetype<_num_facetype; n_facetype++) {

				standard_face_ptr=currentgrid->_standard_face[n_facetype];
				for (num_elem=0; ; standard_face_ptr=standard_face_ptr->_next) {
					num_elem+=standard_face_ptr->_num_face;
					if (!standard_face_ptr->_next) break;
				}
				if (!num_elem) continue;

				tecplotzone=new clTecplot_zone_unstructured(num_var);
				if (!tecplotzone) fatal_err(0,"tecplotzone in clFieldview::convert_to_tecplot_surface");
				tecplotzone->set_zonetype_num_conn("FEQUADRILATERAL");
				tecplotzone->set_datatype("BLOCK");
				if (nn_zone) {
					tecplotzone->set_num_varsharing(num_var);
					tecplotzone->set_varsharingzone(nshare);
				}
				tecplotzone->set_title(_facetypename[n_facetype]);
				tecplotzone->set_num_node(num_node);
				tecplotzone->set_num_elem(num_elem);

				if (!nn_zone) {
					nshare=n_zone;
					tecplotzone->allocate_data();
					dataptr=tecplotzone->get_datastartptr(0);
					fieldviewdataptr=currentgrid->_xcoord;
					for (n_node=0; n_node<num_node; n_node++) (*dataptr++)=(double)(*fieldviewdataptr++);
					dataptr=tecplotzone->get_datastartptr(1);
					fieldviewdataptr=currentgrid->_ycoord;
					for (n_node=0; n_node<num_node; n_node++) (*dataptr++)=(double)(*fieldviewdataptr++);
					dataptr=tecplotzone->get_datastartptr(2);
					fieldviewdataptr=currentgrid->_zcoord;
					for (n_node=0; n_node<num_node; n_node++) (*dataptr++)=(double)(*fieldviewdataptr++);
					if (num_funcvar) {
						fieldviewdataptr=currentgrid->_func[0];
						for (n_funcvar=0; n_funcvar<num_funcvar; n_funcvar++) {
							dataptr=tecplotzone->get_datastartptr(3+n_funcvar);
							for (n_node=0; n_node<num_node; n_node++) (*dataptr++)=(double)(*fieldviewdataptr++);
						}
					}
				}

				tecplotzone->set_num_elem(num_elem);
				num_conn=tecplotzone->get_num_conn();
				if (num_conn!=4) fatal_err("num_conn!=4 in clFieldview::convert_to_tecplot_surface");
				conn.allocate(num_elem,num_conn);
				connptr=conn.get_startptr();
				standard_face_ptr=currentgrid->_standard_face[n_facetype];
				for (; ; standard_face_ptr=standard_face_ptr->_next) {
					num_data=4*standard_face_ptr->_num_face;
					fieldviewconnptr=standard_face_ptr->_conn[0];
					for (n_data=0; n_data<num_data; n_data++) (*connptr++)=(*fieldviewconnptr++)-1;
					if (!standard_face_ptr->_next) break;
				}
				connptr=conn.get_startptr();
				for (n_elem=0; n_elem<num_elem; n_elem++,connptr+=num_conn) if (connptr[3]<0) connptr[3]=connptr[2];
				tecplotzone->set_conn(conn);

				tecplot->set_zone(n_zone,tecplotzone);

				n_zone++;
				nn_zone++;
			}
		}
	}

	if (n_zone!=num_zone) fatal_err("n_zone!=num_zone in clFieldview::convert_to_tecplot_surface");

	if (_filetype!=clFieldview::FV_RESULTS_FILE) tecplot->set_zonesharing();
};

void clFieldview::convert_to_tecplot_volume(clTecplot *tecplot)
{
	if (_num_grid!=1) fatal_err("convert_to_tecplot not implemented for _num_grid!=1");
	fatal_err("clFieldview::convert_to_tecplot_volume not implemented");
};

void clFieldview::add_result_to_tecplot_surface(clTecplot *tecplot)
{
	int n_grid,n_node,num_node,n_var,n_funcvar,n_zone,num_zonetecplot,sharingzone;
	double *dataptr=NULL;
	clTecplot_zone *tecplotzone=NULL;
	clUtils_namespace::my_float *fieldviewdataptr=NULL;
	clFieldview_grid *currentgrid=NULL;

	if (!_num_var) fatal_err("No functions variables in result file");

	if (_filetype!=clFieldview::FV_RESULTS_FILE) 
		fatal_err("_filetype!=FV_RESULTS_FILE in clFieldview::add_result_to_tecplot_surface");

	if (read_bin_debug) fprintf(stdout,"Number of result zones=%d\n",_num_grid);
	num_zonetecplot=tecplot->get_num_zone();
	if (_num_grid>num_zonetecplot) fatal_err("Inconsistent number of zones");
	if (read_bin_debug) fprintf(stdout,"Number of grid zones=%d\n",num_zonetecplot);

	const int num_funcvar=_num_var;
	const int num_var=num_funcvar+3;
	const int num_vartecplot=tecplot->get_num_var();

	if (num_vartecplot<3) fatal_err("Inconsistent number of variables in tecplot file");
	for (n_var=3; n_var<num_vartecplot; n_var++) tecplot->remove_variable(3);
	tecplot->add_variables(num_funcvar);

	tecplot->set_num_var(num_var);
	tecplot->alloc_varname();
	tecplot->set_varname(0,"X");
	tecplot->set_varname(1,"Y");
	tecplot->set_varname(2,"Z");
	for (n_var=3; n_var<num_var; n_var++) tecplot->set_varname(n_var,_varname[n_var-3]);

	n_grid=0;
	for (n_zone=0; n_zone<num_zonetecplot; n_zone++) {

		tecplotzone=tecplot->get_zoneptr(n_zone);
		if (!tecplotzone) 
			fatal_err("tecplotzone not allocated in clFieldview::add_result_to_tecplot_surface");
		if (num_vartecplot==tecplotzone->get_num_varsharing()) {
			tecplotzone->set_num_varsharing(num_var);
			sharingzone=tecplotzone->get_varsharingzone(0);
			tecplotzone->dealloc_varsharingzone();
			tecplotzone->set_num_var(num_var);
			tecplotzone->set_varsharingzone(sharingzone);
			continue;
		}
		if (tecplotzone->get_num_varsharing()) fatal_err("Illegal number of shared variables");
		if (tecplotzone->datatype_point()) fatal_err("Illegal datatype in clFieldview::add_result_to_tecplot_surface");;

		currentgrid=_grids[n_grid];
		if (!currentgrid) 
			fatal_err("currentgrid not allocated in clFieldview::add_result_to_tecplot_surface");

		num_node=currentgrid->_num_node;
		if (num_node!=tecplotzone->get_num_node()) fatal_err("Inconsistent number of nodes");

		fieldviewdataptr=currentgrid->_func[0];
		for (n_funcvar=0; n_funcvar<num_funcvar; n_funcvar++) {
			dataptr=tecplotzone->get_datastartptr(3+n_funcvar);
			for (n_node=0; n_node<num_node; n_node++) (*dataptr++)=(double)(*fieldviewdataptr++);
		}

		tecplotzone->dealloc_varsharingzone();
		tecplotzone->set_num_var(num_var);
		tecplotzone->alloc_varsharingzone();

		n_grid++;

	}

	if (n_grid!=_num_grid) fatal_err("n_grid!=_num_grid in clFieldview::add_result_to_tecplot_surface");
};
