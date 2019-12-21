/*******************************************************************************
  A collection of Tecplot usage functions
  Author: Rolf Verberg
*******************************************************************************/
#include "clTecplot.h"

using clUtils_namespace::warning;
using clUtils_namespace::fatal_err;
using clUtils_namespace::readdqstring;
using clUtils_namespace::readletterstring;
using clUtils_namespace::readletternumberstring;
using clUtils_namespace::uppercase;
using clUtils_namespace::matrix;

const int clTecplot::LINE_LENGTH=256;
const int clTecplot::NAME_LENGTH=256;
const int clTecplot::IOBUFFERSIZE=10000;

static int read_bin_debug=0;

/* Constructor for a node-to-node table with num_node nodes */
clNton_table::clNton_table(const int num_node)
{
	_num_node=num_node;
	_nton=new nton_entry* [_num_node];
	if (!_nton) fatal_err(0,"_nton in clNton_table::clNton_table");
	nton_entry **ntonptr=_nton;
	for (int n_node=0; n_node<num_node; n_node++) (*ntonptr++)=NULL;
	_ntondist=NULL;
};
/* Constructor for a node-to-node table */
/*   alloc_flag: 0: Allocate memory for both the node-to-node table and the node-to-node distance table */
/*               1: Allocate memory for just the node-to-node table */
/*               2: Allocate memory for both the node-to-node distance table */
clNton_table::clNton_table(const int num_node,const int alloc_flag)
{
	_num_node=num_node;
	_nton=NULL;
	_ntondist=NULL;
	if (alloc_flag==0 || alloc_flag==1) {
		_nton=new nton_entry* [_num_node];
		if (!_nton) fatal_err(0,"_nton in clNton_table::clNton_table");
		nton_entry **ntonptr=_nton;
		for (int n_node=0; n_node<num_node; n_node++) (*ntonptr++)=NULL;
	}
	if (alloc_flag==0 || alloc_flag==2) {
		_ntondist=new ntondist_entry* [_num_node];
		if (!_ntondist) fatal_err(0,"_ntondist in clNton_table::clNton_table");
		ntondist_entry **ntondistptr=_ntondist;
		for (int n_node=0; n_node<num_node; n_node++) (*ntondistptr++)=NULL;
	}
};

/* Allocate memory for an node-to-node entry */
void clNton_table::nton_entry_init(nton_entry* nton,const int num_node)
{
	nton->_num_node=num_node;
	if (!num_node) {
		nton->_nodes=NULL;
		return;
	}
	nton->_nodes=new int [num_node];
	if (!nton->_nodes) fatal_err(0,"nton->_nodes in clNton_table::nton_entry_init");
	int n_node,*nodesptr=nton->_nodes;
	for (n_node=0; n_node<num_node; n_node++) (*nodesptr++)=-1;
};

/* Allocate memory for an node-to-node and node-to-node distance entry */
void clNton_table::ntondist_entry_init(ntondist_entry* ntondist,const int num_node)
{
	ntondist->_num_node=num_node;
	ntondist->_sum_r2inv_inv=0.0;
	if (!num_node) {
		ntondist->_nodes=NULL;
		ntondist->_r2inv=NULL;
		return;
	}
	ntondist->_nodes=new int [num_node];
	if (!ntondist->_nodes) fatal_err(0,"ntondist->_nodes in clNton_table::ntondist_entry_init");
	int n_node,*nodesptr=ntondist->_nodes;
	for (n_node=0; n_node<num_node; n_node++) (*nodesptr++)=-1;
	ntondist->_r2inv=new double [num_node];
	if (!ntondist->_r2inv) fatal_err(0,"ntondist->_r2inv in clNton_table::ntondist_entry_init");
	double *r2invptr=ntondist->_r2inv;
	for (n_node=0; n_node<num_node; n_node++) (*r2invptr++)=0.0;
};

/* Create the node to node nearest neighbor lookup table for a structured grid or an unstructured grid without connectivities */
/* RV slow approach, think about splitting to speed up later */
/*   num_node (>0): number of nodes */
/*   num_var (>=3): number of data variables */
/*   dataptr: array of pointers to the start of each data variables */
/*            (This routine assumes that x,y,z are the first three data variables) */
/*   range (>=0.0): 0.0:  Only create the node-to-node table */
/*                  >0.0: Also create the node-to-node distance table, only search for neighbors within the given range */
/*   maxnode (>0): The maximum number of neighbors to keep track of in the node-to-node distance */
void clNton_table::construct_nton_table(const int num_node,const int num_var,
	const double* const* const dataptr,const double range,const int maxnode)
{
	int n_index,indexmax,n_node,nnum_node,nn_node;
	int *indexptr=NULL,*nodesptr=NULL,*ntonptr=NULL;
	double x,y,z,delx,dely,delz,rsq,rsqmax;
	const double rangesq=range*range;
	const double* const xstartptr=dataptr[0];
	const double* const ystartptr=dataptr[1];
	const double* const zstartptr=dataptr[2];
	nton_entry *nton_entryptr=NULL;
	clVector_namespace::clVector<int> index;
	clVector_namespace::clVector<double> ntondist;
	clVector_namespace::clMatrix<int> nton;

	if (num_node<1) fatal_err("num_node<1 in clNton_table::construct_nton_table");
	if (num_var<3) fatal_err("num_var<3 in clNton_table::construct_nton_table");
	if (range<0.0) fatal_err("range<=0.0 in clNton_table::construct_nton_table");
	if (maxnode<1) fatal_err("maxnode<1 in clNton_table::construct_nton_table");

	if (_nton && _num_node!=num_node) {
		delete [] _nton;
		_nton=new nton_entry* [num_node];
		if (!_nton) fatal_err(0,"_nton in clNton_table::construct_nton_table");
		nton_entry **ntonptr=_nton;
		for (n_node=0; n_node<num_node; n_node++) (*ntonptr++)=NULL;
	}
	_num_node=num_node;

	/* Store node to node table if range=0.0 and return */
	if (!range) {
		for (n_node=0; n_node<num_node; n_node++) {
			nton_entryptr=_nton[n_node]=new nton_entry;
			if (!nton_entryptr) fatal_err(0,"nton_entryptr in clNton_table::construct_nton_table");
			nton_entry_init(nton_entryptr,0);
		}
		return;
	}

	/* Allocate and initialize index table */
	index.allocate(num_node,-1);
	ntondist.allocate(maxnode);

	/* Allocate and create node to node table */
	nton.allocate(num_node,maxnode,-1);
	ntonptr=nton.get_startptr();
	indexptr=index.get_startptr();
	for (n_node=0; n_node<num_node; n_node++,indexptr++) {
		x=(*(xstartptr+n_node));
		y=(*(ystartptr+n_node));
		z=(*(zstartptr+n_node));
		ntonptr=nton(n_node);
		nnum_node=0;
		rsqmax=0.0;
		for (nn_node=0; nn_node<num_node; nn_node++) {
			if (n_node==nn_node) continue;
			delx=x-(*(xstartptr+nn_node));
			dely=y-(*(ystartptr+nn_node));
			delz=z-(*(zstartptr+nn_node));
			rsq=delx*delx+dely*dely+delz*delz;
			if (range) if (rsq>rangesq) continue;
			if (nnum_node==maxnode) {
				if (rsq<rsqmax) {
					ntondist(indexmax)=rsq;
					ntonptr[indexmax]=nn_node;
					rsqmax=rsq;
					for (n_index=0; n_index<maxnode; n_index++) {
						if (ntondist(n_index)>rsqmax) {
							indexmax=n_index;
							rsqmax=ntondist(n_index);
						}
					}
				}
			} else {
				if (rsq>rsqmax) {
					indexmax=nnum_node;
					rsqmax=rsq;
				}
				ntondist(nnum_node)=rsq;
				ntonptr[nnum_node++]=nn_node;
			}
		}
		(*indexptr)=nnum_node;
	}

	/* Store node to node table */
	indexptr=index.get_startptr();
	for (n_node=0; n_node<num_node; n_node++,indexptr++) {
		if (*indexptr==-1) fatal_err("failed to get node to node index");
		nton_entryptr=_nton[n_node]=new nton_entry;
		if (!nton_entryptr) fatal_err(0,"nton_entryptr in clNton_table::construct_nton_table");
		if (!(*indexptr)) continue;
		nton_entry_init(nton_entryptr,*indexptr);
		ntonptr=nton(n_node);
		nodesptr=nton_entryptr->_nodes;
		for (n_index=0; n_index<(*indexptr); n_index++) (*nodesptr++)=(*ntonptr++);
	}
};

/* Create the node to node connection lookup table for an unstructured grid */
/*   num_node (>0): number of nodes */
/*   num_elem (>0): number of elements */
/*   num_conn (>0): number of connectivities */
/*   conn: connectivities */
void clNton_table::construct_nton_table(const int num_node,const int num_elem,
	const int num_conn,const int* const conn)
{
	int maxindex=0,n_conn,n_elem,n_index,n_node,nn_node,nnum_node;
	int *indexptr=NULL,*nodesptr=NULL,*ntoeptr=NULL,*ntonptr=NULL;
	const int *connptr=NULL;
	nton_entry *nton_entryptr=NULL;
	clVector_namespace::clVector<int> index;
	clVector_namespace::clMatrix<int> ntoe;
	clVector_namespace::clMatrix<int> nton;

	if (num_node<1) fatal_err("num_node<1 in clNton_table::construct_nton_table");
	if (num_elem<1) fatal_err("num_elem<1 in clNton_table::construct_nton_table");
	if (num_conn<1) fatal_err("num_conn<1 in clNton_table::construct_nton_table");
	
	if (_nton && _num_node!=num_node) {
		delete [] _nton;
		_nton=new nton_entry* [num_node];
		if (!_nton) fatal_err(0,"_nton in clNton_table::construct_nton_table");
		nton_entry **ntonptr=_nton;
		for (n_node=0; n_node<num_node; n_node++) (*ntonptr++)=NULL;
	}
	_num_node=num_node;

	/* Allocate and initialize index table */
	index.allocate(num_node,0);
	const int num_data=num_elem*num_conn;
	connptr=conn;
	for (int n_data=0; n_data<num_data; n_data++) index(*connptr++)++;

	/* Allocate and create node to element table */
	indexptr=index.get_startptr();
	for (n_node=0,maxindex=0; n_node<num_node; n_node++,indexptr++) 
		if ((*indexptr)>maxindex) maxindex=(*indexptr);
	ntoe.allocate(num_node,maxindex,-1);
	index=0;
	connptr=conn;
	for (n_elem=0; n_elem<num_elem; n_elem++) for (n_conn=0; n_conn<num_conn; n_conn++) {
		n_node=(*connptr++);
		ntoe(n_node,index(n_node))=n_elem;
		index(n_node)++;
	}

	/* Allocate and create node to node table */
	const int maxnode=(num_conn-1)*maxindex;
	nton.allocate(num_node,maxnode,-1);
	indexptr=index.get_startptr();
	for (n_node=0; n_node<num_node; n_node++,indexptr++) {
		if (!(*indexptr)) fatal_err("failed to get node to element index");
		ntoeptr=ntoe(n_node);
		ntonptr=nton(n_node);
		nnum_node=0;
		for (n_index=0; n_index<(*indexptr); n_index++) {
			connptr=conn+(*ntoeptr++)*num_conn;
			for (n_conn=0; n_conn<num_conn; n_conn++,connptr++) if ((*connptr)!=n_node) {
				for (nn_node=0; nn_node<nnum_node; nn_node++) 
					if ((*connptr)==ntonptr[nn_node]) break;
				if (nn_node==nnum_node) ntonptr[nnum_node++]=(*connptr);
			}
		}
		(*indexptr)=nnum_node;
	}

	/* Store node to node table */
	indexptr=index.get_startptr();
	for (n_node=0; n_node<num_node; n_node++,indexptr++) {
		if (!(*indexptr)) fatal_err("failed to get node to node index");
		nton_entryptr=_nton[n_node]=new nton_entry;
		if (!nton_entryptr) fatal_err(0,"nton_entryptr in clNton_table::construct_nton_table");
		nton_entry_init(nton_entryptr,*indexptr);
		ntonptr=nton(n_node);
		nodesptr=nton_entryptr->_nodes;
		for (n_index=0; n_index<(*indexptr); n_index++) (*nodesptr++)=(*ntonptr++);
	}
};

/* Set the zone type and number of connectivities */
int clTecplot_zone::set_zonetype_num_conn(const char* const zonetype)
{
	if (!strcmp(zonetype,"ORDERED")) {
		_zonetype=clTecplot::ORDERED;
	} else if (!strcmp(zonetype,"FELINESEG")) {
		_zonetype=clTecplot::FELINESEG;
		set_num_conn(2);
	} else if (!strcmp(zonetype,"FETRIANGLE")) {
		_zonetype=clTecplot::FETRIANGLE;
		set_num_conn(3);
	} else if (!strcmp(zonetype,"FEQUADRILATERAL")) {
		_zonetype=clTecplot::FEQUADRILATERAL;
		set_num_conn(4);
	} else if (!strcmp(zonetype,"FETETRAHEDRON")) {
		_zonetype=clTecplot::FETETRAHEDRON;
		set_num_conn(4);
	} else if (!strcmp(zonetype,"FEBRICK")) {
		_zonetype=clTecplot::FEBRICK;
		set_num_conn(8);
	} else if (!strcmp(zonetype,"FEPOLYGON")) {
		_zonetype=clTecplot::FEPOLYGON;
	} else if (!strcmp(zonetype,"FEPOLYHEDRON")) {
		_zonetype=clTecplot::FEPOLYHEDRON;
	} else {
		return 0;
	}
	return 1;
};

/* Set the zone type and number of connectivities */
/* Use only to read zonetype from binary Tecplot files */
int clTecplot_zone::set_zonetype_num_conn(const int zonetype)
{
	if (zonetype==0) {
		_zonetype=clTecplot::ORDERED;
	} else if (zonetype==1) {
		_zonetype=clTecplot::FELINESEG;
		set_num_conn(2);
	} else if (zonetype==2) {
		_zonetype=clTecplot::FETRIANGLE;
		set_num_conn(3);
	} else if (zonetype==3) {
		_zonetype=clTecplot::FEQUADRILATERAL;
		set_num_conn(4);
	} else if (zonetype==4) {
		_zonetype=clTecplot::FETETRAHEDRON;
		set_num_conn(4);
	} else if (zonetype==5) {
		_zonetype=clTecplot::FEBRICK;
		set_num_conn(8);
	} else if (zonetype==6) {
		_zonetype=clTecplot::FEPOLYGON;
	} else if (zonetype==7) {
		_zonetype=clTecplot::FEPOLYHEDRON;
	} else {
		return 0;
	}
	return 1;
};

bool clTecplot_zone::zonetype_ordered(void) const {return (_zonetype==clTecplot::ORDERED ? TRUE : FALSE);};
//bool clTecplot_zone::zonetype_fequadrilateral(void) const {return (_zonetype==clTecplot::FEQUADRILATERAL ? TRUE : FALSE);};
//bool clTecplot_zone::zonetype_fepolygon(void) const {return (_zonetype==clTecplot::FEPOLYGON ? TRUE : FALSE);};
//bool clTecplot_zone::zonetype_fepolyhedron(void) const {return (_zonetype==clTecplot::FEPOLYHEDRON ? TRUE : FALSE);};

/* Write the zone type and number of connectivities to a stream */
void clTecplot_zone::print_zonetype(FILE *fileptr) const
{
	if (_zonetype==clTecplot::FELINESEG) 
		fprintf(fileptr,"ZONETYPE=FELINESEG\n");
	else if (_zonetype==clTecplot::FETRIANGLE) 
		fprintf(fileptr,"ZONETYPE=FETRIANGLE\n");
	else if (_zonetype==clTecplot::FEQUADRILATERAL) 
		fprintf(fileptr,"ZONETYPE=FEQUADRILATERAL\n");
	else if (_zonetype==clTecplot::FETETRAHEDRON) 
		fprintf(fileptr,"ZONETYPE=FETETRAHEDRON\n");
	else if (_zonetype==clTecplot::FEBRICK) 
		fprintf(fileptr,"ZONETYPE=FEBRICK\n");
	else if (_zonetype==clTecplot::FEPOLYGON) 
		fprintf(fileptr,"ZONETYPE=FEPOLYGON\n");
	else if (_zonetype==clTecplot::FEPOLYGON) 
		fprintf(fileptr,"ZONETYPE=FEPOLYGON\n");
	else if (_zonetype==clTecplot::FEPOLYHEDRON) 
		fprintf(fileptr,"ZONETYPE=FEPOLYHEDRON\n");
};

/* Set the data type */
int clTecplot_zone::set_datatype(const char* const datatype)
{
	if (!strcmp(datatype,"POINT")) {
		_datatype=clTecplot::POINT;
	} else if (!strcmp(datatype,"BLOCK")) {
		_datatype=clTecplot::BLOCK;
	} else {
		return 0;
	}
	return 1;
};

/* Check if the data type is point (return TRUE) or not (return FALSE) */
bool clTecplot_zone::datatype_point(void) const {return (_datatype==clTecplot::POINT ? TRUE : FALSE);};

/* Write the data type to a stream */
void clTecplot_zone::write_datatype(FILE *fileptr) const
{
	if (zonetype_ordered()) fprintf(fileptr,"F=");
	else fprintf(fileptr,"DATAPACKING=");
	if (datatype_point()) fprintf(fileptr,"%s\n","POINT");
	else fprintf(fileptr,"%s\n","BLOCK");
};

/* Add data from the sourcezone zone to the current zone */
void clTecplot_zone::add_data(const clTecplot_zone* const sourcezone)
{
	clVector_namespace::clVector<int> varindex(_num_var,1);
	add_data(sourcezone,varindex);
};
/* Add data from the sourcezone zone to the current zone for each variable that has varindex set */
void clTecplot_zone::add_data(const clTecplot_zone* const sourcezone,const clVector_namespace::clVector<int> &varindex)
{
	int n_data,num_data,n_var,num_var;
	const double *sourcedataptr;
	double *dataptr=NULL;

	if (_zonetype!=sourcezone->get_zonetype()) fatal_err("Incompatible zonetype in clTecplot_zone::add_data");
	if (_datatype!=sourcezone->get_datatype()) fatal_err("Incompatible datatype in clTecplot_zone::add_data");
	if (_num_var!=sourcezone->get_num_var()) fatal_err("Incompatible num_var in clTecplot_zone::add_data");
	if (_num_varsharing!=sourcezone->_num_varsharing) fatal_err("Incompatible num_varsharing in clTecplot_zone::add_data");
	if (_num_varsharing) {
		for (n_var=0; n_var<_num_var; n_var++) 
			if (_varsharingzone(n_var)!=sourcezone->get_varsharingzone(n_var)) fatal_err("Incompatible varsharingzone in clTecplot_zone::add_data");
	}
	if (zonetype_ordered()) {
		if (get_idim()!=sourcezone->get_idim()) fatal_err("Incompatible idim in clTecplot_zone::add_data");
		if (get_jdim()!=sourcezone->get_jdim()) fatal_err("Incompatible jdim in clTecplot_zone::add_data");
		if (get_kdim()!=sourcezone->get_kdim()) fatal_err("Incompatible kdim in clTecplot_zone::add_data");
	} else {
		if (get_num_node()!=sourcezone->get_num_node()) fatal_err("Incompatible num_node in clTecplot_zone::add_data");
	}
	_varmin.deallocate();
	_varmax.deallocate();
	num_var=_num_var-_num_varsharing;
	if (num_var>0) {
		if (zonetype_ordered()) num_data=get_idim()*get_jdim()*get_kdim();
		for (n_var=0; n_var<_num_var; n_var++) if (varindex(n_var)) {
			if (!_num_varsharing || get_varsharingzone(n_var)==-1) {
				if (!zonetype_ordered()) {
					if (get_varlocation(n_var)) num_data=get_num_elem();
					else num_data=get_num_node();
				}
				dataptr=get_datastartptr(n_var);
				sourcedataptr=sourcezone->get_datastartptr(n_var);
				for (n_data=0; n_data<num_data; n_data++) (*dataptr++)+=(*sourcedataptr++);
			}
		}
	}
};

/* Add delta to the data in the current zone */
void clTecplot_zone::translate_data(const double delta)
{
	if (!delta) return;

	clVector_namespace::clVector<int> varindex(_num_var,1);
	translate_data(delta,varindex);
};
/* Add delta to variable n_var in the current zone */
void clTecplot_zone::translate_data(const double delta,const int n_var)
{
	if (!delta) return;

	clVector_namespace::clVector<int> varindex(_num_var,0);
	varindex(n_var)=1;
	translate_data(delta,varindex);
};
/* Add delta to the data in the current zone for each variable that has varindex set */
void clTecplot_zone::translate_data(const double delta,const clVector_namespace::clVector<int> &varindex)
{
	if (!delta) return;

	int n_data,num_data,n_var,num_var;
	double *dataptr=NULL;

	num_var=_num_var-_num_varsharing;
	if (num_var>0) {
		if (zonetype_ordered()) num_data=get_idim()*get_jdim()*get_kdim();
		for (n_var=0; n_var<_num_var; n_var++) if (varindex(n_var)) {
			if (!_num_varsharing || get_varsharingzone(n_var)==-1) {
				if (!zonetype_ordered()) {
					if (get_varlocation(n_var)) num_data=get_num_elem();
					else num_data=get_num_node();
				}
				dataptr=get_datastartptr(n_var);
				for (n_data=0; n_data<num_data; n_data++) (*dataptr++)+=delta;
			}
		}
	}
};

/* Multiply the data in the current zone by factor */
void clTecplot_zone::multiply_data(const double factor)
{
	if (factor==1.0) return;

	clVector_namespace::clVector<int> varindex(_num_var,1);
	multiply_data(factor,varindex);
};
/* Multiply variable n_var in the current zone by factor */
void clTecplot_zone::multiply_data(const double factor,const int n_var)
{
	if (factor==1.0) return;

	clVector_namespace::clVector<int> varindex(_num_var,0);
	varindex(n_var)=1;
	multiply_data(factor,varindex);
};
/* Multiply the data in the current zone by factor for each variable that has varindex set */
void clTecplot_zone::multiply_data(const double factor,const clVector_namespace::clVector<int> &varindex)
{
	if (factor==1.0) return;

	int n_data,num_data,n_var,num_var;
	double *dataptr=NULL;

	num_var=_num_var-_num_varsharing;
	if (num_var>0) {
		if (zonetype_ordered()) num_data=get_idim()*get_jdim()*get_kdim();
		for (n_var=0; n_var<_num_var; n_var++) if (varindex(n_var)) {
			if (!_num_varsharing || get_varsharingzone(n_var)==-1) {
				if (!zonetype_ordered()) {
					if (get_varlocation(n_var)) num_data=get_num_elem();
					else num_data=get_num_node();
				}
				dataptr=get_datastartptr(n_var);
				for (n_data=0; n_data<num_data; n_data++) (*dataptr++)*=factor;
			}
		}
	}
};

/* Find the smallest and/or largest vertex */
double clTecplot_zone::find_min_vertex(const double mincutoff) const
{
	double minvertex=1.e30,maxvertex=0.0;
	find_minmax_vertex(minvertex,maxvertex,mincutoff);
	return minvertex;
};
double clTecplot_zone::find_max_vertex(const double mincutoff) const
{
	double minvertex=1.e30,maxvertex=0.0;
	find_minmax_vertex(minvertex,maxvertex,mincutoff);
	return maxvertex;
};
void clTecplot_zone::find_minmax_vertex(double &minvertex,double &maxvertex,const double mincutoff) const
{
	if (minvertex<0.0) minvertex=1.e30;
	if (maxvertex<0.0) maxvertex=0.0;
	if (_num_varsharing)
		fatal_err("Shared variables in _zones[n_zone] not yet implemeted in clTecplot::find_minmax_vertex");
	double delx,dely,delz,r;
	if (zonetype_ordered()) {
		int n_data,nn_data;
		const int num_data=get_idim()*get_jdim()*get_kdim();
		const double *xxptr=NULL,*yyptr=NULL,*zzptr=NULL;
		const double *xptr=get_datastartptr(0);
		const double *yptr=get_datastartptr(1);
		const double *zptr=get_datastartptr(2);
		for (nn_data=0; nn_data<num_data; nn_data++,xptr++,yptr++,zptr++) {
			xxptr=xptr+1;
			yyptr=yptr+1;
			zzptr=zptr+1;
			for (n_data=nn_data+1; n_data<num_data; n_data++) {
				delx=(*xxptr++)-(*xptr);
				dely=(*yyptr++)-(*yptr);
				delz=(*zzptr++)-(*zptr);
				r=sqrt(delx*delx+dely*dely+delz*delz);
				if (r>=mincutoff) {
					if (r<minvertex) minvertex=r;
					if (r>maxvertex) maxvertex=r;
				}
			}
		}
	} else {
		const int num_conn=get_num_conn();
		if (num_conn!=3) fatal_err("clTecplot_zone::find_minmax_vertex only implemented yet for num_conn=3");
		int n_elem,conn1,conn2,conn3;
		const int num_elem=get_num_elem();
		const int *connptr=get_connstartptr();
		const double* const xstartptr=get_datastartptr(0);
		const double* const ystartptr=get_datastartptr(1);
		const double* const zstartptr=get_datastartptr(2);
		for (n_elem=0; n_elem<num_elem; n_elem++) {
			conn1=(*connptr++);
			conn2=(*connptr++);
			conn3=(*connptr++);
			delx=xstartptr[conn2]-xstartptr[conn1];
			dely=ystartptr[conn2]-ystartptr[conn1];
			delz=zstartptr[conn2]-zstartptr[conn1];
			r=sqrt(delx*delx+dely*dely+delz*delz);
			if (r>=mincutoff) {
				if (r<minvertex) minvertex=r;
				if (r>maxvertex) maxvertex=r;
			}
			delx=xstartptr[conn3]-xstartptr[conn1];
			dely=ystartptr[conn3]-ystartptr[conn1];
			delz=zstartptr[conn3]-zstartptr[conn1];
			r=sqrt(delx*delx+dely*dely+delz*delz);
			if (r>=mincutoff) {
				if (r<minvertex) minvertex=r;
				if (r>maxvertex) maxvertex=r;
			}
			delx=xstartptr[conn3]-xstartptr[conn2];
			dely=ystartptr[conn3]-ystartptr[conn2];
			delz=zstartptr[conn3]-zstartptr[conn2];
			r=sqrt(delx*delx+dely*dely+delz*delz);
			if (r>=mincutoff) {
				if (r<minvertex) minvertex=r;
				if (r>maxvertex) maxvertex=r;
			}
		}
	}
};

/* Calculate the variable ranges and store them in _varmin and _varmax */
void clTecplot_zone::calculate_varrange(void)
{
	int n_data,num_data,n_var;
	double varmin,varmax;
	const double *dataptr=NULL;

	_varmin.allocate(_num_var,0.0);
	_varmax.allocate(_num_var,0.0);
	if (zonetype_ordered()) {
		int idim=get_idim();
		int jdim=get_jdim();
		int kdim=get_kdim();
		num_data=idim*jdim*kdim;
	}
	for (n_var=0; n_var<_num_var; n_var++) if (get_varsharingzone(n_var)==-1) {
		dataptr=get_datastartptr(n_var);
		varmin=(*dataptr);
		varmax=(*dataptr++);
		if (!zonetype_ordered()) {
			if (get_varlocation(n_var)) num_data=get_num_elem();
			else num_data=get_num_node();
		}
		for (n_data=1; n_data<num_data; n_data++,dataptr++) {
			if ((*dataptr)<varmin) varmin=(*dataptr);
			else if ((*dataptr)>varmax) varmax=(*dataptr);
		}
		_varmin(n_var)=varmin;
		_varmax(n_var)=varmax;
	}
};

/* Copy the header of sourcezone to the current zone */
void clTecplot_zone::copy_zoneheader(const clTecplot_zone* const sourcezone)
{
	int n_var;
	_zonetype=sourcezone->_zonetype;
	_datatype=sourcezone->_datatype;
	_num_var=sourcezone->_num_var;
	_num_varsharing=sourcezone->_num_varsharing;
	_varsharingzone.allocate(_num_var);
	if (!_num_varsharing)
		for (n_var=0; n_var<_num_var; n_var++) _varsharingzone(n_var)=-1;
	else
		for (n_var=0; n_var<_num_var; n_var++) _varsharingzone(n_var)=sourcezone->_varsharingzone(n_var);
	set_title(sourcezone->get_title());
	if (zonetype_ordered()) {
		set_idim(sourcezone->get_idim());
		set_jdim(sourcezone->get_jdim());
		set_kdim(sourcezone->get_kdim());
	} else {
		set_num_node(sourcezone->get_num_node());
		set_num_elem(sourcezone->get_num_elem());
		set_num_conn(sourcezone->get_num_conn());
		set_num_var_cc(sourcezone->get_num_var_cc());
		clVector_namespace::clVector<int> varlocation(_num_var);
		for (n_var=0; n_var<_num_var; n_var++) varlocation(n_var)=sourcezone->get_varlocation(n_var);
		set_varlocation(varlocation);
	}
};

/* Create node to node lookup table for the current zone */
/*   range (>=0.0): 0.0:  Only create the node-to-node table */
/*                  >0.0: Also create the node-to-node distance table, only search for neighbors within the given range */
void clTecplot_zone::create_nton_table(clNton_table **nton_table,const double range) const
{
	const int n_zone=-1;
	create_nton_table(nton_table,n_zone,range);
};
void clTecplot_zone::create_nton_table(clNton_table **nton_table,const int n_zone,const double range) const
{
	const int readwriteflag=0;
	FILE *fileptr=NULL;
	create_nton_table(nton_table,n_zone,fileptr,readwriteflag,range);
};
void clTecplot_zone::create_nton_table(clNton_table **nton_table,FILE *fileptr,
		const int readwriteflag,const double range) const
{
	const int n_zone=-1;
	create_nton_table(nton_table,n_zone,fileptr,readwriteflag,range);
};
void clTecplot_zone::create_nton_table(clNton_table **nton_table,const int n_zone,FILE *fileptr,
	const int readwriteflag,const double range) const
{
	int n_node,num_node,nn_node,nnum_node;
	int *nodesptr=NULL;
	clNton_table *nton_tableptr=NULL;
	clNton_table::nton_entry *nton_entryptr=NULL;
	if (fileptr && readwriteflag>0) {
		fscanf(fileptr,"%*s %*s");
		fscanf(fileptr,"%*s %*s %d",&num_node);
		if (num_node!=get_num_node()) 
			fatal_err("inconsistent num_node in clTecplot_zone::create_nton_table while reading node to node lookup table");
		(*nton_table)=new clNton_table(num_node);
		if (!(*nton_table)) fatal_err(0,"nton_table in clTecplot_zone::create_nton_table");
		nton_tableptr=(*nton_table);
		for (n_node=0; n_node<num_node; n_node++) {
			fscanf(fileptr," %d",&nnum_node);
			nton_entryptr=new clNton_table::nton_entry;
			if (!nton_entryptr) fatal_err(0,"nton_entryptr in clTecplot_zone::create_nton_table");
			nton_tableptr->_nton[n_node]=nton_entryptr;
			nton_tableptr->nton_entry_init(nton_entryptr,nnum_node);
			nodesptr=nton_entryptr->_nodes;
			for (nn_node=0; nn_node<nnum_node; nn_node++) fscanf(fileptr,"%d",nodesptr++);
		}
	} else {
		int gridtype;
		int num_conn=0,num_elem=0;
		if (zonetype_ordered()) {
			const int idim=get_idim();
			const int jdim=get_jdim();
			const int kdim=get_kdim();
			if (jdim!=1 || kdim!=1) gridtype=1;	// Structured grid
			else gridtype=-1;					// Unstructured grid without connectivity info
			num_node=idim*jdim*kdim;
		} else {
			gridtype=0;							// Unstructured grid with connectivity info
			num_node=get_num_node();
			num_elem=get_num_elem();
			num_conn=get_num_conn();
		}
		if (_num_varsharing)
			fatal_err("Variables should be unshared in clTecplot_zone::create_nton_table");
		(*nton_table)=new clNton_table(num_node);
		if (!(*nton_table)) fatal_err(0,"nton_table in clTecplot_zone::create_nton_table");
		nton_tableptr=(*nton_table);
		if (!gridtype) {
			/* Unstructured grid with connectivity info */
			nton_tableptr->construct_nton_table(num_node,num_elem,num_conn,get_connstartptr());
		} else if (gridtype==1) {
			/* Structured grid */
			fatal_err("Structured grid not implemented yet in clTecplot_zone::create_nton_table");
		} else if (gridtype==-1) {
			/* Unstructured grid without connectivity info */
			const int num_var=get_num_var();
			const double **dataptr=new const double* [num_var];
			if (!dataptr) fatal_err(0,"dataptr in clTecplot_zone::create_nton_table");
			for (int n_var=0; n_var<num_var; n_var++) dataptr[n_var]=get_datastartptr(n_var);
			nton_tableptr->construct_nton_table(num_node,get_num_var(),dataptr,range);
		} else fatal_err("Illegal value of gridtype in clTecplot_zone::create_nton_table");
		if (fileptr && readwriteflag<0) {
			fprintf(fileptr,"Zone %d\n",n_zone+1);
			fprintf(fileptr,"num_node = %d\n",num_node);
			for (n_node=0; n_node<num_node; n_node++) {
				nton_entryptr=nton_tableptr->_nton[n_node];
				nnum_node=nton_entryptr->_num_node;
				fprintf(fileptr," %d",nnum_node);
				nodesptr=nton_entryptr->_nodes;
				for (nn_node=0; nn_node<nnum_node; nn_node++) 
					fprintf(fileptr," %d",(*nodesptr++));
				fprintf(fileptr,"\n");
			}
		}
	}
};

/* Find nearest neighbor info for a target Tecplot zone relative to the current zone */
/* Assumes that x,y,z are the first three variables */
/*   eps: The range within which a target node is considered an exact location match with a local node */
void clTecplot_zone::find_nearest_neighbor(const clTecplot_zone* const targetzone,
	const clNton_table* const nton_table,clNton_table **ntondist_table) const
{
	int readwriteflag=0;
	double eps=0.0;
	FILE *fileptr=NULL;
	find_nearest_neighbor(targetzone,nton_table,ntondist_table,eps,fileptr,readwriteflag);
};
void clTecplot_zone::find_nearest_neighbor(const clTecplot_zone* const targetzone,
	const clNton_table* const nton_table,clNton_table **ntondist_table,const double eps) const
{
	int readwriteflag=0;
	FILE *fileptr=NULL;
	find_nearest_neighbor(targetzone,nton_table,ntondist_table,eps,fileptr,readwriteflag);
};
void clTecplot_zone::find_nearest_neighbor(const clTecplot_zone* const targetzone,const clNton_table* const nton_table,
	clNton_table **ntondist_table,FILE *fileptr,const int readwriteflag) const
{
	double eps=0.0;
	find_nearest_neighbor(targetzone,nton_table,ntondist_table,eps,fileptr,readwriteflag);
};
void clTecplot_zone::find_nearest_neighbor(const clTecplot_zone* const targetzone,const clNton_table* const nton_table,
	clNton_table **ntondist_table,const double eps,FILE *fileptr,const int readwriteflag) const
{
	const int numpart=1;
	const double minx=0.0,xmin=0.0,xmax=0.0,deltax=1.e30;

	int num_node;
	int *splitgrid_ltogptr=NULL;
	clVector_namespace::clVector<int> splitgrid_num_node;
	clVector_namespace::clMatrix<int> splitgrid_ltog;
	if (_num_varsharing) fatal_err("Variables should be unshared in clTecplot_zone::find_nearest_neighbor");
	if (zonetype_ordered()) {
		const int idim=get_idim();
		const int jdim=get_jdim();
		const int kdim=get_kdim();
		if (jdim!=1 || kdim!=1) fatal_err("Structured zones not yet implemented in clTecplot_zone::find_nearest_neighbor");	// Structured grid
		num_node=idim*jdim*kdim;
	} else {
		num_node=get_num_node();
	}

	splitgrid_num_node.allocate(numpart);
	splitgrid_num_node(0)=num_node;
	splitgrid_ltog.allocate(numpart,num_node);
	splitgrid_ltogptr=splitgrid_ltog.get_startptr();
	for (int n_node=0; n_node<num_node; n_node++) (*splitgrid_ltogptr++)=n_node;
	find_nearest_neighbor(targetzone,nton_table,ntondist_table,eps,numpart,minx,xmin,xmax,deltax,
		splitgrid_num_node,splitgrid_ltog,fileptr,readwriteflag);
};
void clTecplot_zone::find_nearest_neighbor(const clTecplot_zone* const targetzone,const clNton_table* const nton_table,
	clNton_table **ntondist_table,const double eps,const int numpart,const double minx,
	const double xmin,const double xmax,const double deltax,
	const clVector_namespace::clVector<int> &splitgrid_num_node,const clVector_namespace::clMatrix<int> &splitgrid_ltog,
	FILE *fileptr,const int readwriteflag) const
{
	int npart,n_node_target,num_node_target,nn_node,nnum_node;
	clNton_table *ntondist_tableptr=NULL;
	clNton_table::ntondist_entry *ntondist_entryptr=NULL;
	void *io_ptr=NULL;

	if (fileptr && readwriteflag>0) {
		io_ptr=(void *)&num_node_target;
		if (fread(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(2,"nearest neighbor file");
		if (targetzone->zonetype_ordered()) {
			const int idim=targetzone->get_idim();
			const int jdim=targetzone->get_jdim();
			const int kdim=targetzone->get_kdim();
			if (jdim!=1 || kdim!=1) fatal_err("Structured zones not yet implemented in clTecplot_zone::find_nearest_neighbor");
			if (num_node_target!=idim*jdim*kdim)
				fatal_err("inconsistent num_node_target in clTecplot_zone::find_nearest_neighbor while reading nearest neighbor info");
		} else {
			if (num_node_target!=targetzone->get_num_node())
				fatal_err("inconsistent num_node_target in clTecplot_zone::find_nearest_neighbor while reading nearest neighbor info");
		}
		(*ntondist_table)=new clNton_table(num_node_target,2);
		if (!(*ntondist_table)) fatal_err(0,"ntondist_table in clTecplot_zone::find_nearest_neighbor");
		ntondist_tableptr=(*ntondist_table);
		for (n_node_target=0; n_node_target<num_node_target; n_node_target++) {
			io_ptr=(void *)&nnum_node;
			if (fread(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(2,"nearest neighbor file");
			ntondist_entryptr=new clNton_table::ntondist_entry;
			if (!ntondist_entryptr) fatal_err(0,"ntondist_entryptr in clTecplot_zone::find_nearest_neighbor");
			ntondist_tableptr->_ntondist[n_node_target]=ntondist_entryptr;
			ntondist_tableptr->ntondist_entry_init(ntondist_entryptr,nnum_node);
			io_ptr=(void *)ntondist_entryptr->_nodes;
			if (fread(io_ptr,sizeof(int),nnum_node,fileptr)!=(unsigned)(nnum_node)) fatal_err(2,"nearest neighbor file");
			io_ptr=(void *)ntondist_entryptr->_r2inv;
			if (fread(io_ptr,sizeof(double),nnum_node,fileptr)!=(unsigned)(nnum_node)) fatal_err(2,"nearest neighbor file");
			io_ptr=(void *)&ntondist_entryptr->_sum_r2inv_inv;
			if (fread(io_ptr,sizeof(double),1,fileptr)!=1) fatal_err(2,"nearest neighbor file");
		}
	} else {
		int nn_nodemin;
		int *nodesptr=NULL,*n_nodesptr=NULL;
		const int *ltogptr=NULL;
		double x,y,z,delx,dely,delz,r2,r2min,r2inv,sum_r2inv;
		double *distptr=NULL,epssq=eps*eps;
		const clNton_table *nton_tableptr=NULL;
		const clNton_table::nton_entry *nton_entryptr=NULL;
		const int num_var=get_num_var();
		if (_num_varsharing || targetzone->_num_varsharing) 
			fatal_err("Variables should be unshared in clTecplot_zone::find_nearest_neighbor");
		if (targetzone->zonetype_ordered()) {
			const int idim=targetzone->get_idim();
			const int jdim=targetzone->get_jdim();
			const int kdim=targetzone->get_kdim();
			if (jdim!=1 || kdim!=1) fatal_err("Structured zones not yet implemented in clTecplot_zone::find_nearest_neighbor");	// Structured grid
			num_node_target=idim*jdim*kdim;
		} else {
			num_node_target=targetzone->get_num_node();
			if (get_num_var_cc() || targetzone->get_num_var_cc()) {
				if (get_varlocation(0) || targetzone->get_varlocation(0) ||
					get_varlocation(1) || targetzone->get_varlocation(1) ||
					get_varlocation(2) || targetzone->get_varlocation(2)) 
					fatal_err("Cell-centered x,y or z variables not yet implemented in clTecplot_zone::find_nearest_neighbor");
			}
		}
		nton_tableptr=nton_table;
		(*ntondist_table)=new clNton_table(num_node_target,2);
		if (!(*ntondist_table)) fatal_err(0,"ntondist_table in clTecplot_zone::find_nearest_neighbor");
		ntondist_tableptr=(*ntondist_table);
		const double* const xstartptr=get_datastartptr(0);
		const double* const ystartptr=get_datastartptr(1);
		const double* const zstartptr=get_datastartptr(2);
		const double *xtargetptr=targetzone->get_datastartptr(0);
		const double *ytargetptr=targetzone->get_datastartptr(1);
		const double *ztargetptr=targetzone->get_datastartptr(2);
		if (fileptr && readwriteflag<0) {
			io_ptr=(void *)&num_node_target;
			if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(3,"nearest neighbor file");
		}
		for (n_node_target=0; n_node_target<num_node_target; n_node_target++) {
			x=(*xtargetptr++);
			y=(*ytargetptr++);
			z=(*ztargetptr++);
			if (numpart>1) {
				npart=(int)((x-xmin)/deltax);
				if (npart<0|| npart>=numpart) fatal_err("npart out of bounds in clTecplot_zone::find_nearest_neighbor");
			} else npart=0;
			nnum_node=splitgrid_num_node(npart);
			ltogptr=splitgrid_ltog(npart);
			r2min=deltax*deltax;
			for (nn_node=0,nn_nodemin=-1; nn_node<nnum_node; nn_node++,ltogptr++) {
				delx=(*(xstartptr+(*ltogptr)))-x;
				dely=(*(ystartptr+(*ltogptr)))-y;
				delz=(*(zstartptr+(*ltogptr)))-z;
				r2=delx*delx+dely*dely+delz*delz;
				if (r2<epssq) {
					r2min=0.0;
					nn_nodemin=(*ltogptr);
					break;
				}
				if (r2<r2min) {
					r2min=r2;
					nn_nodemin=(*ltogptr);
					if (!r2min) break;
				}
			}
			if (nn_nodemin<0) fatal_err("unable to find nn_nodemin in clTecplot_zone::find_nearest_neighbor");
			if (r2min>0.0 && (x-xmin)<minx && npart>0) {
				nnum_node=splitgrid_num_node(npart-1);
				ltogptr=splitgrid_ltog(npart-1);
				for (nn_node=0; nn_node<nnum_node; nn_node++,ltogptr++) {
					delx=(*(xstartptr+(*ltogptr)))-x;
					dely=(*(ystartptr+(*ltogptr)))-y;
					delz=(*(zstartptr+(*ltogptr)))-z;
					r2=delx*delx+dely*dely+delz*delz;
					if (r2<r2min) {
						r2min=r2;
						nn_nodemin=(*ltogptr);
					}
				}
			}
			if (r2min>0.0 && (xmax-x)<minx && npart+1<numpart) {
				nnum_node=splitgrid_num_node(npart+1);
				ltogptr=splitgrid_ltog(npart+1);
				for (nn_node=0; nn_node<nnum_node; nn_node++,ltogptr++) {
					delx=(*(xstartptr+(*ltogptr)))-x;
					dely=(*(ystartptr+(*ltogptr)))-y;
					delz=(*(zstartptr+(*ltogptr)))-z;
					r2=delx*delx+dely*dely+delz*delz;
					if (r2<r2min) {
						r2min=r2;
						nn_nodemin=(*ltogptr);
					}
				}
			}
			/* Store nearest neighbor distances in node to node table */
			nton_entryptr=nton_tableptr->_nton[nn_nodemin];
			if (!r2min) nnum_node=0;
			else nnum_node=nton_entryptr->_num_node;
			ntondist_entryptr=new clNton_table::ntondist_entry;
			if (!ntondist_entryptr) fatal_err(0,"ntondist_entryptr in clTecplot_zone::find_nearest_neighbor");
			ntondist_tableptr->_ntondist[n_node_target]=ntondist_entryptr;
			ntondist_tableptr->ntondist_entry_init(ntondist_entryptr,nnum_node+1);
			nodesptr=ntondist_entryptr->_nodes;
			(*nodesptr++)=nn_nodemin;
			n_nodesptr=nton_entryptr->_nodes;
			for (nn_node=0; nn_node<nnum_node; nn_node++,nodesptr++) *nodesptr=(*n_nodesptr++);
			if (r2min) r2inv=1.0/r2min;
			else r2inv=0.0;
			sum_r2inv=r2inv;
			distptr=ntondist_entryptr->_r2inv;
			(*distptr++)=r2inv;
			nodesptr=ntondist_entryptr->_nodes;
			for (nn_node=0; nn_node<nnum_node; nn_node++) {
				delx=(*(xstartptr+(*(++nodesptr))))-x;
				dely=(*(ystartptr+(*(nodesptr))))-y;
				delz=(*(zstartptr+(*(nodesptr))))-z;
				r2inv=1.0/(delx*delx+dely*dely+delz*delz);
				sum_r2inv+=r2inv;
				(*distptr++)=r2inv;
			}
			if (sum_r2inv) sum_r2inv=1.0/sum_r2inv;
			else sum_r2inv=1.0;
			ntondist_entryptr->_sum_r2inv_inv=sum_r2inv;
			if (fileptr && readwriteflag<0) {
				nn_node=nnum_node+1;
				io_ptr=(void *)&nn_node;
				if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(3,"nearest neighbor file");
				io_ptr=(void *)ntondist_entryptr->_nodes;
				if (fwrite(io_ptr,sizeof(int),nn_node,fileptr)!=(unsigned)(nn_node)) fatal_err(3,"nearest neighbor file");
				io_ptr=(void *)ntondist_entryptr->_r2inv;
				if (fwrite(io_ptr,sizeof(double),nn_node,fileptr)!=(unsigned)(nn_node)) fatal_err(3,"nearest neighbor file");
				io_ptr=(void *)&sum_r2inv;
				if (fwrite(io_ptr,sizeof(double),1,fileptr)!=1) fatal_err(3,"nearest neighbor file");
			}
		}
	}
};

/* Find nearest neighbor info for a target Tecplot zone relative to the current zone */
/* Use only for cell-centered data */
/* Assumes that x,y,z are the first three variables */
/*   eps: The range within which a target node is considered an exact location match with a local node */
void clTecplot_zone::find_nearest_neighbor_cc(const clTecplot_zone* const targetzone,const clNton_table* const etoe_table,
	clNton_table **etoedist_table) const
{
	int readwriteflag=0;
	double eps=0.0;
	FILE *fileptr=NULL;
	find_nearest_neighbor_cc(targetzone,etoe_table,etoedist_table,eps,fileptr,readwriteflag);
};
void clTecplot_zone::find_nearest_neighbor_cc(const clTecplot_zone* const targetzone,const clNton_table* const etoe_table,
	clNton_table **etoedist_table,const double eps) const
{
	int readwriteflag=0;
	FILE *fileptr=NULL;
	find_nearest_neighbor_cc(targetzone,etoe_table,etoedist_table,eps,fileptr,readwriteflag);
};
void clTecplot_zone::find_nearest_neighbor_cc(const clTecplot_zone* const targetzone,const clNton_table* const etoe_table,
	clNton_table **etoedist_table,FILE *fileptr,const int readwriteflag) const
{
	double eps=0.0;
	find_nearest_neighbor_cc(targetzone,etoe_table,etoedist_table,eps,fileptr,readwriteflag);
};
void clTecplot_zone::find_nearest_neighbor_cc(const clTecplot_zone* const targetzone,const clNton_table* const etoe_table,
	clNton_table **etoedist_table,const double eps,FILE *fileptr,const int readwriteflag) const
{
	int n_elem_target,num_elem_target,nnum_elem;
	clNton_table *etoedist_tableptr=NULL;
	clNton_table::ntondist_entry *etoedist_entryptr=NULL;
	void *io_ptr=NULL;

	if (fileptr && readwriteflag>0) {
		io_ptr=(void *)&num_elem_target;
		if (fread(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(2,"nearest neighbor file");
		if (num_elem_target!=targetzone->get_num_elem())
			fatal_err("inconsistent num_elem in clTecplot_zone::find_nearest_neighbor_cc while reading nearest neighbor info");
		(*etoedist_table)=new clNton_table(num_elem_target,2);
		if (!(*etoedist_table)) fatal_err(0,"etoedist_table in clTecplot_zone::find_nearest_neighbor_cc");
		etoedist_tableptr=(*etoedist_table);
		for (n_elem_target=0; n_elem_target<num_elem_target; n_elem_target++) {
			io_ptr=(void *)&nnum_elem;
			if (fread(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(2,"nearest neighbor file");
			etoedist_entryptr=new clNton_table::ntondist_entry;
			if (!etoedist_entryptr) fatal_err(0,"etoedist_entryptr in clTecplot_zone::find_nearest_neighbor_cc");
			etoedist_tableptr->_ntondist[n_elem_target]=etoedist_entryptr;
			etoedist_tableptr->ntondist_entry_init(etoedist_entryptr,nnum_elem);
			io_ptr=(void *)etoedist_entryptr->_nodes;
			if (fread(io_ptr,sizeof(int),nnum_elem,fileptr)!=(unsigned)(nnum_elem)) fatal_err(2,"nearest neighbor file");
			io_ptr=(void *)etoedist_entryptr->_r2inv;
			if (fread(io_ptr,sizeof(double),nnum_elem,fileptr)!=(unsigned)(nnum_elem)) fatal_err(2,"nearest neighbor file");
			io_ptr=(void *)&etoedist_entryptr->_sum_r2inv_inv;
			if (fread(io_ptr,sizeof(double),1,fileptr)!=1) fatal_err(2,"nearest neighbor file");
		}
	} else {
		int n_conn,num_conn,n_elem,n_elemmin,nn_elem;
		int *elemsptr=NULL,*eelemsptr=NULL;
		const int *connptr=NULL;
		double x,y,z,delx,dely,delz,r2,r2min,r2inv,sum_r2inv;
		double *distptr=NULL,epssq=eps*eps;
		double *xdataptr=NULL,*ydataptr=NULL,*zdataptr=NULL;
		const double *xxdataptr=NULL,*yydataptr=NULL,*zzdataptr=NULL;
		const double *xdatastartptr=NULL,*ydatastartptr=NULL,*zdatastartptr=NULL;
		const double *xdatatargetptr=NULL,*ydatatargetptr=NULL,*zdatatargetptr=NULL;
		const double* const xstartptr=get_datastartptr(0);
		const double* const ystartptr=get_datastartptr(1);
		const double* const zstartptr=get_datastartptr(2);
		const double* const xtargetstartptr=targetzone->get_datastartptr(0);
		const double* const ytargetstartptr=targetzone->get_datastartptr(1);
		const double* const ztargetstartptr=targetzone->get_datastartptr(2);
		clVector_namespace::clVector<double> xdata,ydata,zdata,xdatatarget,ydatatarget,zdatatarget;
		const clNton_table *etoe_tableptr=NULL;
		const clNton_table::nton_entry *etoe_entryptr=NULL;
		const int num_var=get_num_var();
		const int num_elem=get_num_elem();
		num_elem_target=targetzone->get_num_elem();
		int idummy=get_varlocation(0)+get_varlocation(1)+get_varlocation(2);
		if (!idummy) {
			xdata.allocate(num_elem,0.0);
			ydata.allocate(num_elem,0.0);
			zdata.allocate(num_elem,0.0);
			xdataptr=xdata.get_startptr();
			ydataptr=ydata.get_startptr();
			zdataptr=zdata.get_startptr();
			num_conn=get_num_conn();
			connptr=get_connstartptr();
			for (n_elem=0; n_elem<num_elem; n_elem++) {
				for (n_conn=0; n_conn<num_conn; n_conn++) {
					(*xdataptr)+=(*(xstartptr+(*connptr)));
					(*ydataptr)+=(*(ystartptr+(*connptr)));
					(*zdataptr)+=(*(zstartptr+(*connptr++)));
				}
				(*xdataptr++)/=num_conn;
				(*ydataptr++)/=num_conn;
				(*zdataptr++)/=num_conn;
			}
			xdatastartptr=xdata.get_startptr();
			ydatastartptr=ydata.get_startptr();
			zdatastartptr=zdata.get_startptr();
		} else if (idummy==3) {
			xdatastartptr=xstartptr;
			ydatastartptr=ystartptr;
			zdatastartptr=zstartptr;
		} else fatal_err("Nodal/cell-centered mix of x,y or z variables not allowed in clTecplot_zone::find_nearest_neighbor_cc");
		idummy=targetzone->get_varlocation(0)+targetzone->get_varlocation(1)+targetzone->get_varlocation(2);
		if (!idummy) {
			xdatatarget.allocate(num_elem_target,0.0);
			ydatatarget.allocate(num_elem_target,0.0);
			zdatatarget.allocate(num_elem_target,0.0);
			xdataptr=xdatatarget.get_startptr();
			ydataptr=ydatatarget.get_startptr();
			zdataptr=zdatatarget.get_startptr();
			num_conn=targetzone->get_num_conn();
			connptr=targetzone->get_connstartptr();
			for (n_elem_target=0; n_elem_target<num_elem_target; n_elem_target++) {
				for (n_conn=0; n_conn<num_conn; n_conn++) {
					(*xdataptr)+=(*(xtargetstartptr+(*connptr)));
					(*ydataptr)+=(*(ytargetstartptr+(*connptr)));
					(*zdataptr)+=(*(ztargetstartptr+(*connptr++)));
				}
				(*xdataptr++)/=num_conn;
				(*ydataptr++)/=num_conn;
				(*zdataptr++)/=num_conn;
			}
			xdatatargetptr=xdatatarget.get_startptr();
			ydatatargetptr=ydatatarget.get_startptr();
			zdatatargetptr=zdatatarget.get_startptr();
		} else if (idummy==3) {
			xdatatargetptr=xtargetstartptr;
			ydatatargetptr=ytargetstartptr;
			zdatatargetptr=ztargetstartptr;
		} else fatal_err("Nodal/cell-centered mix of x,y or z variables not allowed in clTecplot_zone::find_nearest_neighbor_cc");
		etoe_tableptr=etoe_table;
		(*etoedist_table)=new clNton_table(num_elem_target,2);
		if (!(*etoedist_table)) fatal_err(0,"etoedist_table in clTecplot_zone::find_nearest_neighbor_cc");
		etoedist_tableptr=(*etoedist_table);
		if (fileptr && readwriteflag<0) {
			io_ptr=(void *)&num_elem_target;
			if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(3,"nearest neighbor file");
		}
		for (n_elem_target=0; n_elem_target<num_elem_target; n_elem_target++) {
			x=(*xdatatargetptr++);
			y=(*ydatatargetptr++);
			z=(*zdatatargetptr++);
			xxdataptr=xdatastartptr;
			yydataptr=ydatastartptr;
			zzdataptr=zdatastartptr;
			for (n_elem=0,n_elemmin=-1; n_elem<num_elem; n_elem++) {
				delx=(*xxdataptr++)-x;
				dely=(*yydataptr++)-y;
				delz=(*zzdataptr++)-z;
				r2=delx*delx+dely*dely+delz*delz;
				if (r2<epssq) {
					r2min=0.0;
					n_elemmin=n_elem;
					break;
				}
				if (n_elemmin<0) {
					r2min=r2;
					n_elemmin=n_elem;
				} else if (r2<r2min) {
					r2min=r2;
					n_elemmin=n_elem;
					if (!r2min) break;
				}
			}
			if (n_elemmin<0) fatal_err("unable to find n_elemmin in clTecplot_zone::find_nearest_neighbor_cc");
			/* Store nearest neighbor distances in elem to elem table */
			if (etoe_tableptr) etoe_entryptr=etoe_tableptr->_nton[n_elemmin];
			if (!r2min || !etoe_entryptr) nnum_elem=0;
			else nnum_elem=etoe_entryptr->_num_node;
			etoedist_entryptr=new clNton_table::ntondist_entry;
			if (!etoedist_entryptr) fatal_err(0,"etoedist_entryptr in clTecplot_zone::find_nearest_neighbor_cc");
			etoedist_tableptr->_ntondist[n_elem_target]=etoedist_entryptr;
			etoedist_tableptr->ntondist_entry_init(etoedist_entryptr,nnum_elem+1);
			elemsptr=etoedist_entryptr->_nodes;
			(*elemsptr++)=n_elemmin;
			if (etoe_entryptr) {
				fatal_err("NOT TESTED YET");
				eelemsptr=etoe_entryptr->_nodes;
				for (nn_elem=0; nn_elem<nnum_elem; nn_elem++,elemsptr++) *elemsptr=(*eelemsptr++);
			}
			if (r2min) r2inv=1.0/r2min;
			else r2inv=0.0;
			sum_r2inv=r2inv;
			distptr=etoedist_entryptr->_r2inv;
			(*distptr++)=r2inv;
			elemsptr=etoedist_entryptr->_nodes+1;
			for (nn_elem=0; nn_elem<nnum_elem; nn_elem++) {
				fatal_err("NOT TESTED YET");
				delx=(*(xdatastartptr+(*elemsptr)))-x;
				dely=(*(ydatastartptr+(*elemsptr)))-y;
				delz=(*(zdatastartptr+(*elemsptr++)))-z;
				r2inv=1.0/(delx*delx+dely*dely+delz*delz);
				sum_r2inv+=r2inv;
				(*distptr++)=r2inv;
			}
			if (sum_r2inv) sum_r2inv=1.0/sum_r2inv;
			else sum_r2inv=1.0;
			etoedist_entryptr->_sum_r2inv_inv=sum_r2inv;
			if (fileptr && readwriteflag<0) {
				nnum_elem++;
				io_ptr=(void *)&nnum_elem;
				if (fwrite(io_ptr,sizeof(int),1,fileptr)!=1) fatal_err(3,"nearest neighbor file");
				io_ptr=(void *)etoedist_entryptr->_nodes;
				if (fwrite(io_ptr,sizeof(int),nnum_elem,fileptr)!=(unsigned)(nnum_elem)) fatal_err(3,"nearest neighbor file");
				io_ptr=(void *)etoedist_entryptr->_r2inv;
				if (fwrite(io_ptr,sizeof(double),nnum_elem,fileptr)!=(unsigned)(nnum_elem)) fatal_err(3,"nearest neighbor file");
				io_ptr=(void *)&sum_r2inv;
				if (fwrite(io_ptr,sizeof(double),1,fileptr)!=1) fatal_err(3,"nearest neighbor file");
			}
		}
	}
};

/* Interpolate the current zone onto a target Tecplot zone using nearest neighbors */
/* Assumes that x,y,z are the first three variables */
/* WARNING: No checking for matching boundaries is performed so serious extrapolation errors can occur */
/*   range (>=0.0): >0.0: only use neighbors within the given range */
/*   eps: The range within which a target node is considered an exact location match with a local node */
/*   usevar: The variable mapping: usevar(n_var_target) is the current variable that maps to the n_var_target's target variable */
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone) const
{
	const int num_var_target=targetzone->get_num_var();
	if (_num_var>num_var_target) fatal_err("Number of variables mismatch in clTecplot_zone::interpolate");
	const double eps=0.0,range=0.0;
	clString_namespace::clString nnfilename;
	clVector_namespace::clVector<int> usevar(num_var_target,-1);
	for (int n_var_target=0; n_var_target<num_var_target; n_var_target++) usevar[n_var_target]=n_var_target;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double range) const
{
	const int num_var_target=targetzone->get_num_var();
	if (_num_var>num_var_target) fatal_err("Number of variables mismatch in clTecplot_zone::interpolate");
	const double eps=0.0;
	clString_namespace::clString nnfilename;
	clVector_namespace::clVector<int> usevar(num_var_target,-1);
	for (int n_var_target=0; n_var_target<num_var_target; n_var_target++) usevar(n_var_target)=n_var_target;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double eps,const char* const nnfilename) const
{
	const int num_var_target=targetzone->get_num_var();
	if (_num_var>num_var_target) fatal_err("Number of variables mismatch in clTecplot_zone::interpolate");
	const double range=0.0;
	clVector_namespace::clVector<int> usevar(num_var_target,-1);
	for (int n_var_target=0; n_var_target<num_var_target; n_var_target++) usevar(n_var_target)=n_var_target;
	clString_namespace::clString nnnfilename(nnfilename);
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double eps,
	const clString_namespace::clString &nnfilename) const
{
	const int num_var_target=targetzone->get_num_var();
	if (_num_var>num_var_target) fatal_err("Number of variables mismatch in clTecplot_zone::interpolate");
	const double range=0.0;
	clVector_namespace::clVector<int> usevar(num_var_target,-1);
	for (int n_var_target=0; n_var_target<num_var_target; n_var_target++) usevar(n_var_target)=n_var_target;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double eps,const double range) const
{
	const int num_var_target=targetzone->get_num_var();
	if (_num_var>num_var_target) fatal_err("Number of variables mismatch in clTecplot_zone::interpolate");
	clString_namespace::clString nnfilename;
	clVector_namespace::clVector<int> usevar(num_var_target,-1);
	for (int n_var_target=0; n_var_target<num_var_target; n_var_target++) usevar(n_var_target)=n_var_target;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double eps,const double range,
	const clString_namespace::clString &nnfilename) const
{
	const int num_var_target=targetzone->get_num_var();
	if (_num_var>num_var_target) fatal_err("Number of variables mismatch in clTecplot_zone::interpolate");
	clVector_namespace::clVector<int> usevar(num_var_target,-1);
	for (int n_var_target=0; n_var_target<num_var_target; n_var_target++) usevar(n_var_target)=n_var_target;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const clVector_namespace::clVector<int> &usevar) const
{
	const double eps=0.0,range=0.0;
	clString_namespace::clString nnfilename;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double range,
	const clVector_namespace::clVector<int> &usevar) const
{
	const double eps=0.0;
	clString_namespace::clString nnfilename;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double eps,const double range,
	const clVector_namespace::clVector<int> &usevar) const
{
	clString_namespace::clString nnfilename;
	interpolate(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double eps,const double range,
	const clVector_namespace::clVector<int> &usevar,const char* const nnfilename) const
{
	clString_namespace::clString nnnfilename(nnfilename);
	interpolate(targetzone,eps,range,usevar,nnnfilename);
};
void clTecplot_zone::interpolate(clTecplot_zone* const targetzone,const double eps,const double range,
	const clVector_namespace::clVector<int> &usevar,const clString_namespace::clString &nnfilename) const
{
	int readflag=0;
	FILE *nnfileptr=NULL;
	if (nnfilename.get_len()) {
		nnfileptr=fopen(nnfilename.get_string(),"rb");
		if (nnfileptr) {
			readflag=1;
			fclose(nnfileptr);
			nnfileptr=NULL;
		}
	}

	/* Create the node to node lookup table for the current zone */
	clNton_table *nton_table=NULL;
	if (!readflag) create_nton_table(&nton_table,range);

	/* Find nearest neighbor current zone info for target node */
	clNton_table *ntondist_table=NULL;
	if (!nnfilename.get_len()) {
		find_nearest_neighbor(targetzone,nton_table,&ntondist_table);
	} else {
		if (readflag) {
			nnfileptr=fopen(nnfilename.get_string(),"rb");
			if (!nnfileptr) fatal_err(1,"nearest neighbor table file in clTecplot_zone::interpolate");
			find_nearest_neighbor(targetzone,nton_table,&ntondist_table,eps,nnfileptr,1);
			fclose(nnfileptr);
		} else {
			nnfileptr=fopen(nnfilename.get_string(),"wb");
			if (!nnfileptr) fatal_err(1,"nearest neighbor table file in clTecplot_zone::interpolate");
			find_nearest_neighbor(targetzone,nton_table,&ntondist_table,eps,nnfileptr,-1);
			fclose(nnfileptr);
		}
	}

	delete nton_table;

	/* Perform interpolation from current to target zone */
	int n_node,n_node_target,num_node_target,nn_node,nnum_node,n_var,n_var_target;
	const int *nodesptr=NULL;
	double data;
	const double *distptr=NULL;
	const clNton_table* const ntondist_tableptr=ntondist_table;
	const clNton_table::ntondist_entry *ntondist_entryptr=NULL;
	if (targetzone->_num_varsharing) 
		fatal_err("Variables in targetzone should be unshared in clTecplot_zone::interpolate");
	if (_num_varsharing) 
		fatal_err("Variables in current zone should be unshared in clTecplot_zone::interpolate");
	const int num_var_target=targetzone->get_num_var();
	for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (usevar(n_var_target)>=0) {
		if (usevar(n_var_target)>=_num_var) fatal_err("Illegal value of usevar in clTecplot_zone::interpolate");
	}
	if (targetzone->zonetype_ordered()) {
		const int idim=targetzone->get_idim();
		const int jdim=targetzone->get_jdim();
		const int kdim=targetzone->get_kdim();
		num_node_target=idim*jdim*kdim;
	} else {
		num_node_target=targetzone->get_num_node();
	}
	const double **dataptr=new const double* [_num_var];
	if (!dataptr) fatal_err(0,"dataptr");
	for (n_var=0; n_var<_num_var; n_var++) dataptr[n_var]=get_datastartptr(n_var);
	double **datatargetptr=new double* [num_var_target];
	if (!datatargetptr) fatal_err(0,"datatargetptr");
	for (n_var_target=0; n_var_target<num_var_target; n_var_target++) 
		datatargetptr[n_var_target]=targetzone->get_datastartptr(n_var_target);
	clVector_namespace::clVector<int> intpolvar(num_var_target,0);
	int num_var_intpol_nodal=0,num_var_intpol_cc=0;
	for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (usevar(n_var_target)>=0) {
		if (targetzone->get_varlocation(n_var_target)) {
			intpolvar(n_var_target)=-1;
			num_var_intpol_cc++;
		} else {
			intpolvar(n_var_target)=1;
			num_var_intpol_nodal++;
		}
	}
	
	/* Interpolate the nodal variables */
	if (num_var_intpol_nodal) for (n_node_target=0; n_node_target<num_node_target; n_node_target++) {
		ntondist_entryptr=ntondist_tableptr->_ntondist[n_node_target];
		nnum_node=ntondist_entryptr->_num_node;
		nodesptr=ntondist_entryptr->_nodes;
		distptr=ntondist_entryptr->_r2inv;
		for (n_var_target=3; n_var_target<num_var_target; n_var_target++) 
			if (intpolvar(n_var_target)>0) *(datatargetptr[n_var_target])=0.0;
		for (nn_node=0; nn_node<nnum_node; nn_node++,nodesptr++,distptr++) {
			for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)>0) {
				n_var=usevar(n_var_target);
				if (n_var>=0) {
					data=(*(dataptr[n_var]+(*nodesptr)));
					if (!(*distptr)) {
						*(datatargetptr[n_var_target])=data;
					} else {
						*(datatargetptr[n_var_target])+=data*(*distptr);
					}
				}
			}
			if (!(*distptr)) break;
		}
		if (nn_node==nnum_node) {
			for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)>0) 
				*(datatargetptr[n_var_target])*=ntondist_entryptr->_sum_r2inv_inv;
		}
		for (n_var_target=0; n_var_target<num_var_target; n_var_target++) 
			if (intpolvar(n_var_target)>0) datatargetptr[n_var_target]++;
	}

	/* Interpolate the cell-centered variables */
	if (num_var_intpol_cc) if (get_num_var_cc() || targetzone->get_num_var_cc()) {

		if (targetzone->zonetype_ordered()) 
			fatal_err("Cell-entered variables for a structured grid not allowed in clTecplot_zone::interpolate");
		if (get_varlocation(0) || targetzone->get_varlocation(0) ||
			get_varlocation(1) || targetzone->get_varlocation(1) ||
			get_varlocation(2) || targetzone->get_varlocation(2)) 
			fatal_err("Cell-centered variables require nodal coordinates in clTecplot_zone::interpolate");
		int n_conn,n_elem;
		const int num_node=get_num_node();
		const int num_conn=get_num_conn();
		if (targetzone->get_num_conn()!=num_conn) 
			fatal_err("targetzone->get_num_conn()!=num_conn in clTecplot_zone::interpolate");
		const int num_elem=get_num_elem();

		int n_data;
		const int num_data=num_conn*num_elem;
		const int *numnodalptr=NULL;
		double *datanodalptr=NULL;
		const double *ddataptr=NULL;
		clVector_namespace::clVector<int> numnodal(num_node,0);
		const int *connptr=get_connstartptr();
		for (n_data=0; n_data<num_data; n_data++) numnodal(*connptr++)++;
		clVector_namespace::clVector<double> *datanodal=new clVector_namespace::clVector<double> [_num_var];
		for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)<0) {
			n_var=usevar(n_var_target);
			if (n_var>=0) {
				datanodal[n_var].allocate(num_node,0.0);
				datanodalptr=datanodal[n_var].get_startptr();
				connptr=get_connstartptr();
				ddataptr=dataptr[n_var];
				for (n_elem=0; n_elem<num_elem; n_elem++) {
					for (n_conn=0; n_conn<num_conn; n_conn++) datanodalptr[*connptr++]+=(*ddataptr);
					ddataptr++;
				}
				numnodalptr=numnodal.get_startptr();
				datanodalptr=datanodal[n_var].get_startptr();
				for (n_node=0; n_node<num_node; n_node++,numnodalptr++,datanodalptr++) 
					if (*numnodalptr) (*datanodalptr)/=(*numnodalptr);
				dataptr[n_var]=datanodal[n_var].get_startptr();
			}
		}
		
		const int num_elem_target=targetzone->get_num_elem();
		const int *conn_targetptr=targetzone->get_connstartptr();
		clVector_namespace::clVector<double> newdata(num_var_target);
		for (n_elem=0; n_elem<num_elem_target; n_elem++) {
			for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)<0) *(datatargetptr[n_var_target])=0.0;
			for (n_conn=0; n_conn<num_conn; n_conn++) {
				ntondist_entryptr=ntondist_tableptr->_ntondist[*conn_targetptr++];
				nnum_node=ntondist_entryptr->_num_node;
				nodesptr=ntondist_entryptr->_nodes;
				distptr=ntondist_entryptr->_r2inv;
				newdata=0.0;
				for (nn_node=0; nn_node<nnum_node; nn_node++,nodesptr++,distptr++) {
					for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)<0) {
						n_var=usevar(n_var_target);
						if (n_var>=0) {
							data=(*(dataptr[n_var]+(*nodesptr)));
							if (!(*distptr)) {
								newdata(n_var_target)=data;
							} else {
								newdata(n_var_target)+=data*(*distptr);
							}
						}
					}
					if (!(*distptr)) break;
				}
				if (nn_node==nnum_node) {
					for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)<0)
						newdata(n_var_target)*=ntondist_entryptr->_sum_r2inv_inv;
				}
				for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)<0) 
					(*(datatargetptr[n_var_target]))+=newdata(n_var_target);
			}
			for (n_var_target=0; n_var_target<num_var_target; n_var_target++) if (intpolvar(n_var_target)<0) {
				*(datatargetptr[n_var_target])/=num_conn;
				datatargetptr[n_var_target]++;
			}
		}
		delete [] datanodal;
	}

	for (n_var=0; n_var<_num_var; n_var++) dataptr[n_var]=NULL;
	delete [] dataptr;
	for (n_var_target=0; n_var_target<num_var_target; n_var_target++) datatargetptr[n_var_target]=NULL;
	delete [] datatargetptr;

	delete ntondist_table;
};

/* Interpolate the current zone onto a different Tecplot zone using nearest neighbor elements */
/* Use only for cell-centered data */
/* Assumes that x,y,z are the first three variables */
/* WARNING: No checking about matching boundaries is performed so serious extrapolation errors can occur */
/*   range (>=0.0): >0.0: only use neighbors within the given range */
/*   eps: The range within which a target node is considered an exact location match with a local node */
/*   usevar: The variable mapping: usevar(n_var_target) is the current variable that maps to the n_var_target's target variable */
void clTecplot_zone::interpolate_cc(clTecplot_zone* const targetzone) const
{
	const int num_var_target=targetzone->get_num_var();
	if (_num_var>num_var_target) fatal_err("Number of variables mismatch in clTecplot_zone::interpolate_cc");
	const double eps=0.0,range=0.0;
	clString_namespace::clString nnfilename;
	clVector_namespace::clVector<int> usevar(num_var_target,-1);
	for (int n_var_target=0; n_var_target<num_var_target; n_var_target++) usevar[n_var_target]=n_var_target;
	interpolate_cc(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate_cc(clTecplot_zone* const targetzone,const clVector_namespace::clVector<int> &usevar) const
{
	const double eps=0.0,range=0.0;
	clString_namespace::clString nnfilename;
	interpolate_cc(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate_cc(clTecplot_zone* const targetzone,const double eps,const double range,
	const clVector_namespace::clVector<int> &usevar) const
{
	clString_namespace::clString nnfilename;
	interpolate_cc(targetzone,eps,range,usevar,nnfilename);
};
void clTecplot_zone::interpolate_cc(clTecplot_zone* const targetzone,const double eps,const double range,
	const clVector_namespace::clVector<int> &usevar,const char* const nnfilename) const
{
	clString_namespace::clString nnnfilename(nnfilename);
	interpolate_cc(targetzone,eps,range,usevar,nnnfilename);
};
void clTecplot_zone::interpolate_cc(clTecplot_zone* const targetzone,const double eps,const double range,
	const clVector_namespace::clVector<int> &usevar,const clString_namespace::clString &nnfilename) const
{
	int readflag=0;
	FILE *nnfileptr=NULL;
	if (nnfilename.get_len()) {
		nnfileptr=fopen(nnfilename.get_string(),"rb");
		if (nnfileptr) {
			readflag=1;
			fclose(nnfileptr);
			nnfileptr=NULL;
		}
	}

	/* Create the element to element lookup table for the current zone */
	clNton_table *etoe_table=NULL;
//	if (!readflag) create_etoe_table(&etoe_table,range);	Currently use only the closest of the nearest neighbors

	/* Find nearest neighbor current zone info for target node */
	clNton_table *etoedist_table=NULL;
	if (!nnfilename.get_len()) {
		find_nearest_neighbor_cc(targetzone,etoe_table,&etoedist_table);
	} else {
		if (readflag) {
			nnfileptr=fopen(nnfilename.get_string(),"rb");
			if (!nnfileptr) fatal_err(1,"nearest neighbor table file in clTecplot_zone::interpolate_cc");
			find_nearest_neighbor_cc(targetzone,etoe_table,&etoedist_table,eps,nnfileptr,1);
			fclose(nnfileptr);
		} else {
			nnfileptr=fopen(nnfilename.get_string(),"wb");
			if (!nnfileptr) fatal_err(1,"nearest neighbor table file in clTecplot_zone::interpolate_cc");
			find_nearest_neighbor_cc(targetzone,etoe_table,&etoedist_table,eps,nnfileptr,-1);
			fclose(nnfileptr);
		}
	}

//	delete etoe_table;

	/* Perform interpolation from current to target zone */
	int n_var,n_var_target,num_var_target;
	const clNton_table* const etoedist_tableptr=etoedist_table;
	const clNton_table::ntondist_entry *etoedist_entryptr=NULL;
	if (_num_varsharing || targetzone->_num_varsharing) 
		fatal_err("Variables in should be unshared in clTecplot_zone::interpolate_cc");
	if (zonetype_ordered() || targetzone->zonetype_ordered()) 
		fatal_err("Zones must be unstructured in clTecplot_zone::interpolate_cc");
	num_var_target=targetzone->get_num_var();
	for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (usevar(n_var_target)>=0) {
		if (usevar(n_var_target)>=_num_var) fatal_err("Illegal value of usevar in clTecplot_zone::interpolate_cc");
	}
	const double **dataptr=new const double* [_num_var];
	if (!dataptr) fatal_err(0,"dataptr");
	for (n_var=0; n_var<_num_var; n_var++) dataptr[n_var]=get_datastartptr(n_var);
	double **datatargetptr=new double* [num_var_target];
	if (!datatargetptr) fatal_err(0,"datatargetptr");
	for (n_var_target=0; n_var_target<num_var_target; n_var_target++) 
		datatargetptr[n_var_target]=targetzone->get_datastartptr(n_var_target);
	for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (usevar(n_var_target)>=0) {
		if (!targetzone->get_varlocation(n_var_target)) 
			fatal_err("Interpolating nodal variables not implemented in clTecplot_zone::interpolate_cc");
	}

	/* Interpolate the cell-centered variables */
	int n_elem_target,nn_elem,nnum_elem;
	const int num_elem=get_num_elem();
	const int num_elem_target=targetzone->get_num_elem();
	const int *elemsptr=NULL;
	double data;
	const double *distptr=NULL;
	clVector_namespace::clVector<double> newdata(num_var_target);
	for (n_elem_target=0; n_elem_target<num_elem_target; n_elem_target++) {
		etoedist_entryptr=etoedist_tableptr->_ntondist[n_elem_target];
		nnum_elem=etoedist_entryptr->_num_node;
		elemsptr=etoedist_entryptr->_nodes;
		distptr=etoedist_entryptr->_r2inv;
		newdata=0.0;
		for (nn_elem=0; nn_elem<nnum_elem; nn_elem++,elemsptr++,distptr++) {
			for (n_var_target=3; n_var_target<num_var_target; n_var_target++) {
				n_var=usevar(n_var_target);
				if (n_var>=0) {
					data=(*(dataptr[n_var]+(*elemsptr)));
					if (!(*distptr)) {
						newdata(n_var_target)=data;
					} else {
						newdata(n_var_target)+=data*(*distptr);
					}
				}
			}
			if (!(*distptr)) break;
		}
		if (nn_elem==nnum_elem) {
			for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (usevar(n_var_target)>=0)
				newdata(n_var_target)*=etoedist_entryptr->_sum_r2inv_inv;
		}
		for (n_var_target=3; n_var_target<num_var_target; n_var_target++) if (usevar(n_var_target)>=0) {
			(*(datatargetptr[n_var_target]))=newdata(n_var_target);
			datatargetptr[n_var_target]++;
		}
	}

	for (n_var=0; n_var<_num_var; n_var++) dataptr[n_var]=NULL;
	delete [] dataptr;
	for (n_var_target=0; n_var_target<num_var_target; n_var_target++) datatargetptr[n_var_target]=NULL;
	delete [] datatargetptr;

	delete etoedist_table;
};

void clTecplot::allocate(const int num_zone,const int num_var)
{
	_version=102;
//	_filetype=0;
	_num_zone=num_zone;
	_num_var=num_var;
	_varname.allocate(_num_var);
	_zones=new clTecplot_zone* [_num_zone];
	if (!_zones) clUtils_namespace::fatal_err(0,"_zones in clTecplot");
	for (int n_zone=0; n_zone<_num_zone; n_zone++) _zones[n_zone]=NULL;
	_firstzone=NULL;
};

/* Remove zones firstzone through and including lastzone */
void clTecplot::remove_zone(const int firstzone,const int lastzone)
{
	for (int n_zone=lastzone; n_zone>=firstzone; n_zone--) remove_zone(n_zone);
};
/* Remove targetzone */
void clTecplot::remove_zone(const int targetzone)
{
	int n_var,n_zone;
	clTecplot_zone *currentzone=NULL;

	/* Check for variable sharing */
	if (!_zonesharing.get_dim()) set_zonesharing();

	if (_num_zonesharing) if (_zonesharing(targetzone))
		fatal_err("removing a shared variable zone not implemented in clTecplot::remove_zone");

	/* Reorganize variable sharing zones */
	if (_num_zonesharing) {
		clVector_namespace::clVector<int> zone_reftable(_num_zone);
		for (n_zone=0; n_zone<targetzone; n_zone++) zone_reftable(n_zone)=n_zone;
		zone_reftable(targetzone)=-1;
		for (n_zone=targetzone+1; n_zone<_num_zone; n_zone++) zone_reftable(n_zone)=n_zone-1;
		for (n_zone=0; n_zone<_num_zone; n_zone++) {
			if (n_zone==targetzone) continue;
			currentzone=_zones[n_zone];
			if (!currentzone->_num_varsharing) continue;
			for (n_var=0; n_var<_num_var; n_var++) 
				currentzone->_varsharingzone(n_var)=zone_reftable(currentzone->_varsharingzone(n_var));
		}
	}

	/* Remove targetzone */
	if (!_num_zone || !_zones) fatal_err("!_num_zone || !_zones in clTecplot::remove_zone");
	if (!_zones[targetzone]) fatal_err("!_zones[targetzone] in clTecplot::remove_zone");
	if (!targetzone) 
		_firstzone=_zones[targetzone+1];
	else if (targetzone==_num_zone-1)
		_zones[targetzone-1]->_next=NULL;
	else if (targetzone>0 && targetzone<_num_zone-1)
		_zones[targetzone-1]->_next=_zones[targetzone+1];
	else 
		fatal_err("targetzone out of range in clTecplot::remove_zone");
	delete _zones[targetzone];
	_num_zone--;

	/* Reorganize zones */
	delete [] _zones;
	_zones=new clTecplot_zone* [_num_zone];
	if (!_zones) fatal_err(0,"!_zones in clTecplot::remove_zone");
	currentzone=_firstzone;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		if (!currentzone) fatal_err("!currentzone in clTecplot::remove_zone");;
		_zones[n_zone]=currentzone;
		currentzone=currentzone->_next;
	}
	set_zonesharing();
};

/* Insert zone at end */
void clTecplot::insert_zone(const int _zonetype)
{
	insert_zone(_num_zone-1,_zonetype);
};
/* Insert zone after targetzone */
void clTecplot::insert_zone(const int targetzone,const int _zonetype)
{
	int n_var,n_zone;
	clTecplot_zone *currentzone=NULL,*newzone=NULL;

	/* Check for variable sharing */
	if (!_zonesharing.get_dim()) set_zonesharing();

	/* Reorganize variable sharing zones */
	if (_num_zonesharing && targetzone<_num_zone-1) {
		clVector_namespace::clVector<int> zone_reftable(_num_zone);
		for (n_zone=0; n_zone<=targetzone; n_zone++) zone_reftable(n_zone)=n_zone;
		for (n_zone=targetzone+1; n_zone<_num_zone; n_zone++) zone_reftable(n_zone)=n_zone+1;
		for (n_zone=0; n_zone<_num_zone; n_zone++) {
			currentzone=_zones[n_zone];
			if (!currentzone->_num_varsharing) continue;
			for (n_var=0; n_var<_num_var; n_var++) 
				currentzone->_varsharingzone(n_var)=zone_reftable(currentzone->_varsharingzone(n_var));
		}
	}

	if (_zonetype==ORDERED) 
		newzone=new clTecplot_zone_structured;
	else
		newzone=new clTecplot_zone_unstructured;
	if (!newzone) fatal_err(0,"newzone in clTecplot::insert_zone");

	/* Insert zone after targetzone */
	if (!_zones[targetzone]) fatal_err("!_zones[targetzone] in clTecplot::insert_zone");
	if (targetzone==_num_zone-1)
		newzone->_next=NULL;
	else if (targetzone>=0 && targetzone<_num_zone-1)
		newzone->_next=_zones[targetzone]->_next;
	else 
		fatal_err("targetzone out of range in clTecplot::insert_zone");
	_zones[targetzone]->_next=newzone;

	_num_zone++;

	/* Reorganize zones */
	delete [] _zones;
	_zones=NULL;
	_zones=new clTecplot_zone* [_num_zone];
	if (!_zones) fatal_err(0,"!_zones in clTecplot::insert_zone");
	currentzone=_firstzone;
	_zones[0]=currentzone;
	for (n_zone=1; n_zone<_num_zone; n_zone++) {
		currentzone=currentzone->_next;
		if (!currentzone) fatal_err("!currentzone in clTecplot::insert_zone");;
		_zones[n_zone]=currentzone;
	}
	set_zonesharing();
};

/* Copy sourcezone to end */
void clTecplot::copy_zone(const clTecplot_zone* const sourcezone)
{
	copy_zone(_num_zone,sourcezone);
};
/* Copy sourcezone after targetzone  */
void clTecplot::copy_zone(const int targetzone,const clTecplot_zone* const sourcezone)
{
	if (_num_var!=sourcezone->get_num_var()) fatal_err("Inconsistent num_var in clTecplot::copy_zone");

	int n_data,num_data,n_var,n_zone;
	int *connptr=NULL;
	const int *sourceconnptr=NULL;
	double *dataptr=NULL;
	const double *sourcedataptr=NULL;
	clTecplot_zone *currentzone=NULL,*newzone=NULL;

	/* Check for variable sharing */
	if (targetzone && !_zonesharing.get_dim()) set_zonesharing();

	/* Reorganize variable sharing zones */
	clVector_namespace::clVector<int> zone_reftable(_num_zone);
	if (_num_zonesharing) {
		for (n_zone=0; n_zone<targetzone; n_zone++) zone_reftable(n_zone)=n_zone;
		for (n_zone=targetzone; n_zone<_num_zone; n_zone++) zone_reftable(n_zone)=n_zone+1;
		for (n_zone=0; n_zone<_num_zone; n_zone++) {
			currentzone=_zones[n_zone];
			if (!currentzone->_num_varsharing) continue;
			for (n_var=0; n_var<_num_var; n_var++) 
				currentzone->_varsharingzone(n_var)=zone_reftable(currentzone->_varsharingzone(n_var));
		}
	}

	/* Create new zone */
	if (sourcezone->zonetype_ordered()) 
		newzone=new clTecplot_zone_structured;
	else
		newzone=new clTecplot_zone_unstructured;
	if (!newzone) fatal_err(0,"newzone in clTecplot::copy_zone");
	newzone->copy_zoneheader(sourcezone);
	if (newzone->_num_varsharing) fatal_err("Adding shared variables not allowed in clTecplot::copy_zone");

	/* Copy zone */
	if (sourcezone->_num_varsharing) fatal_err("Shared zones not implemented yet in clTecplot::copy_zone");
	if (sourcezone->zonetype_ordered()) {
		const int idim=sourcezone->get_idim();
		const int jdim=sourcezone->get_jdim();
		const int kdim=sourcezone->get_kdim();
		num_data=idim*jdim*kdim;
	}
	newzone->allocate_data();
	for (n_var=0; n_var<_num_var; n_var++) {
		sourcedataptr=sourcezone->get_datastartptr(n_var);
		if (!sourcedataptr) fatal_err("Unable to get sourcedataptr in clTecplot::copy_zone");
		dataptr=newzone->get_datastartptr(n_var);
		if (!dataptr) fatal_err("Unable to get dataptr in clTecplot::copy_zone");
		if (!sourcezone->zonetype_ordered()) {
			if (sourcezone->get_varlocation(n_var)) num_data=newzone->get_num_elem();
			else num_data=newzone->get_num_node();
		}
		for (n_data=0; n_data<num_data; n_data++) (*dataptr++)=(*sourcedataptr++);
	}

	if (!sourcezone->zonetype_ordered()) {
		const int num_elem=sourcezone->get_num_elem();
		const int num_conn=sourcezone->get_num_conn();
		clVector_namespace::clMatrix<int> conn(num_elem,num_conn);
		connptr=conn.get_startptr();
		num_data=num_elem*num_conn;
		sourceconnptr=sourcezone->get_connstartptr();
		for (n_data=0; n_data<num_data; n_data++) (*connptr++)=(*sourceconnptr++);
		newzone->set_conn(conn);
	}

	/* Copy zone before targetzone */
	if (!targetzone) {
		_firstzone=newzone;
		if (_num_zone) newzone->_next=_zones[0];
	} else if (targetzone==_num_zone) {
		_zones[_num_zone-1]->_next=newzone;
	} else if (targetzone>0 && targetzone<_num_zone) {
		_zones[targetzone-1]->_next=newzone;
		newzone->_next=_zones[targetzone];
	} else {
		fatal_err("targetzone out of range in clTecplot::copy_zone");
	}

	_num_zone++;

	/* Reorganize zones */
	delete [] _zones;
	_zones=NULL;
	_zones=new clTecplot_zone* [_num_zone];
	if (!_zones) fatal_err(0,"!_zones in clTecplot::copy_zone");
	currentzone=_firstzone;
	_zones[0]=currentzone;
	for (n_zone=1; n_zone<_num_zone; n_zone++) {
		currentzone=currentzone->_next;
		if (!currentzone) fatal_err("!currentzone in clTecplot::copy_zone");;
		_zones[n_zone]=currentzone;
	}
	set_zonesharing();
};

/* Split zone n_zone into numpart parts */
void clTecplot::split_zone(const int n_zone,const int numpart)
{
	clTecplot_zone *currentzone=_zones[n_zone];
	if (currentzone->zonetype_ordered()) 
		fatal_err("clTecplot::split_zone not yet implemented for structured zones");
	else split_zone_unstructured(n_zone,numpart);
};
void clTecplot::split_zone_unstructured(const int n_zone,const int numpart)
{
	int num_conn,n_data,num_data,num_elem,npart,n_var;
	int *connptr=NULL;
	const int *connorgptr=NULL;
	double *dataptr=NULL;
	const double *dataorgptr=NULL;
	const char *currenttitle=NULL;
	char newtitle[NAME_LENGTH];
	clVector_namespace::clMatrix<int> conn;
	clVector_namespace::clVector<double> data;
	clTecplot_zone *currentzone=NULL,*newzone=NULL;

	currentzone=_zones[n_zone];
	if (currentzone->zonetype_ordered()) fatal_err("Zone splitting not implemented for ordered zones");

	num_elem=currentzone->get_num_elem();
	num_conn=currentzone->get_num_conn();
	if (numpart<2 || numpart>num_elem) fatal_err("Illegal value of numpart in clTecplot::split_zone");
	currenttitle=currentzone->get_title();

	int num_elem_part_target=num_elem/numpart;
	if (num_elem%num_elem_part_target) num_elem_part_target++;
	clVector_namespace::clVector<int> num_elem_part(numpart);
	clVector_namespace::clVector<int> lowelem_part(numpart);
	for (npart=0; npart<numpart; npart++) {
		lowelem_part(npart)=npart*num_elem_part_target;
		num_elem_part(npart)=num_elem_part_target;
	}
	num_elem_part(numpart-1)=num_elem-lowelem_part(npart-1);

	/* Create new zones */
	int num_varsharing=currentzone->_num_varsharing;
	for (npart=1; npart<numpart; npart++) {
		insert_zone(n_zone+npart-1,currentzone->_zonetype);
		newzone=_zones[n_zone+npart];
		newzone->copy_zoneheader(currentzone);
		newzone->set_num_elem(num_elem_part(npart));
		newzone->_num_varsharing=_num_var;
		if (!num_varsharing) {
			for (n_var=0; n_var<_num_var; n_var++) if (!newzone->get_varlocation(n_var))
				newzone->set_varsharingzone(n_var,n_zone);
		} else {
			for (n_var=0; n_var<_num_var; n_var++) if (!newzone->get_varlocation(n_var)) {
				if (currentzone->get_varsharingzone(n_var)==-1) 
					newzone->set_varsharingzone(n_var,n_zone);
				else
					newzone->set_varsharingzone(n_var,currentzone->get_varsharingzone(n_var));
			}
		}
		num_data=num_elem_part(npart);
		if (newzone->get_num_var_cc()) {
			newzone->allocate_data();
			for (n_var=0; n_var<_num_var; n_var++) if (newzone->get_varlocation(n_var)) {
				data.allocate(num_data);
				dataptr=data.get_startptr();
				dataorgptr=currentzone->get_datastartptr(n_var)+lowelem_part(npart);
				for (n_data=0; n_data<num_data; n_data++) (*dataptr++)=(*dataorgptr++);
				newzone->set_data(data,n_var);
			}
		}
		conn.allocate(num_data,num_conn);
		connptr=conn.get_startptr();
		connorgptr=currentzone->get_connstartptr()+num_conn*lowelem_part(npart);
		num_data*=num_conn;
		for (n_data=0; n_data<num_data; n_data++) (*connptr++)=(*connorgptr++);
		newzone->set_conn(conn);
		sprintf(newtitle,"%s Part %d",currenttitle,npart+1);
		newzone->set_title(newtitle);
	}

	/* Reduce currentzone */
	num_data=num_elem_part(0);
	if (currentzone->get_num_var_cc()) {
		for (n_var=0; n_var<_num_var; n_var++) if (currentzone->get_varlocation(n_var)) {
			data.allocate(num_data);
			dataptr=data.get_startptr();
			dataorgptr=currentzone->get_datastartptr(n_var);
			for (n_data=0; n_data<num_data; n_data++) (*dataptr++)=(*dataorgptr++);
			currentzone->set_data(data,n_var);
		}
	}
	currentzone->set_num_elem(num_data);
	conn.allocate(num_data,num_conn);
	connptr=conn.get_startptr();
	connorgptr=currentzone->get_connstartptr();
	num_data*=num_conn;
	for (n_data=0; n_data<num_data; n_data++) (*connptr++)=(*connorgptr++);
	currentzone->set_conn(conn);
	sprintf(newtitle,"%s Part %d",currenttitle,1);
	currentzone->set_title(newtitle);

	set_zonesharing();
};

/* Identify zones that share variables */
void clTecplot::set_zonesharing(void)
{
	int n_var,n_zone;
	_num_zonesharing=0;
	clTecplot_zone *currentzone=NULL;
	_zonesharing.deallocate();
	if (!_num_zone) return;
	_zonesharing.allocate(_num_zone,0);
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->_num_varsharing) {
			for (n_var=0; n_var<_num_var; n_var++) 
				if (currentzone->_varsharingzone(n_var)!=-1) _zonesharing(currentzone->_varsharingzone(n_var))++;
		}
	}
	for (n_zone=0; n_zone<_num_zone; n_zone++) if (_zonesharing(n_zone)) _num_zonesharing++;
	if (!_num_zonesharing && _zonesharing.get_dim()) _zonesharing.deallocate();
};

/* Read a binary Tecplot file */
void clTecplot::read(const char* const filename)
{
	int inputfiletype=0;

	/*  Magic number  */
	char readchar;
	FILE *fileptr=fopen(filename,"rb");
	if (!fileptr) fatal_err(1,filename);
	void *io_ptr=(void *)(&readchar);
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (readchar=='#') inputfiletype++;
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (readchar=='!') inputfiletype++;
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (readchar=='T') inputfiletype++;
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (readchar=='D') inputfiletype++;
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (readchar=='V') inputfiletype++;
	fclose(fileptr);

	/* Read Tecplot file */
	if (inputfiletype==5) {
		read_bin(filename);
	} else {
		read_ascii(filename);
	}
};

/* Read an ascii Tecplot file */
void clTecplot::read_ascii(const char* const filename)
{
	/*  Get and initialize the number of variables  */
	_num_var=get_num_var(filename);

	/*  Allocate arrays for coordinate and variable names  */
	if (_num_var) {
		_varname.allocate(_num_var);
	} else {
		fatal_err("Unable to obtain the number of variables in clTecplot::read_ascii");
	}

	/*  Read the header and initialize the title and variable names  */
	read_fileheader(filename);

	/*  Read all zones  */
	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) fatal_err(1,filename);
	int ierror,idim,jdim,kdim,n_var,num_var,n_data,num_data=0;
	int *connptr=NULL;
	double *dataptr=NULL;
	double **ddataptr=NULL;
	clVector_namespace::clVector<int> varlocation(_num_var);
	clVector_namespace::clMatrix<int> *conn=NULL;
	clVector_namespace::clMatrix3D<double> *data=NULL;
	clTecplot_zone *newzone=NULL;
	clTecplot_zone *currentzone=NULL;
	clTecplot_zone_readheader *zoneheader=NULL;
	int n,num=1000,finishflag;
	char c,line[LINE_LENGTH];
	struct dataset {
		int _num;
		clVector_namespace::clMatrix<double> _data;
		struct dataset *_next;
	};
	struct dataset *firstdataset=NULL,*currentdataset=NULL,*newdataset=NULL;
	for (;;) {
		finishflag=0;
		zoneheader=new clTecplot_zone_readheader(_num_var);
		if (!zoneheader) fatal_err(0,"zoneheader in clTecplot::read_ascii");
		ierror=read_zoneheader(fileptr,zoneheader);
		if (!ierror) fatal_err(2,filename);
		else if (ierror==EOF) {
			delete zoneheader;
			break;
		}
		if (zoneheader->zonetype_ordered()) {
			newzone=new clTecplot_zone_structured(_num_var);
			if (!newzone) fatal_err(0,"newzone in clTecplot::read_ascii");
			idim=zoneheader->get_idim();
			jdim=zoneheader->get_jdim();
			kdim=zoneheader->get_kdim();
			if (idim*jdim*kdim==0) {
				int maxdim=idim;
				if (jdim>maxdim) maxdim=jdim;
				if (kdim>maxdim) maxdim=kdim;
				if (maxdim) {
					if (!idim) idim++;
					if (!jdim) jdim++;
					if (!kdim) kdim++;
				}
			}
			newzone->set_idim(idim);
			newzone->set_jdim(jdim);
			newzone->set_kdim(kdim);
		} else {
			newzone=new clTecplot_zone_unstructured(_num_var);
			if (!newzone) fatal_err(0,"newzone in clTecplot::read_ascii");
			newzone->set_num_node(zoneheader->get_num_node());
			newzone->set_num_elem(zoneheader->get_num_elem());
			newzone->set_num_conn(zoneheader->get_num_conn());
			newzone->set_num_var_cc(zoneheader->get_num_var_cc());
			for (n_var=0; n_var<_num_var; n_var++) varlocation(n_var)=zoneheader->get_varlocation(n_var);
			newzone->set_varlocation(varlocation);
		}
		newzone->_zonetype=zoneheader->_zonetype;
		newzone->_num_varsharing=zoneheader->_num_varsharing;
		if (newzone->_num_varsharing) 
			for (n_var=0; n_var<_num_var; n_var++) 
				newzone->_varsharingzone(n_var)=zoneheader->get_varsharingzone(n_var);
		newzone->set_title(zoneheader->get_title());
		newzone->set_datatype(zoneheader->_datatype);
		if (!zoneheader->zonetype_ordered()) if (newzone->get_num_var_cc() && newzone->_num_varsharing) {
			for (n_var=0; n_var<_num_var; n_var++) if (varlocation(n_var) && newzone->get_varsharingzone(n_var)!=-1)
				fatal_err("Illegal cell-centered variable sharing");
		}
		delete zoneheader;
		
		if (++_num_zone==1) _firstzone=newzone;
		else currentzone->_next=newzone;
		currentzone=newzone;

		num_var=_num_var-currentzone->_num_varsharing;
		if (num_var) {
			ierror=0;
			if (currentzone->zonetype_ordered()) {
				idim=currentzone->get_idim();
				jdim=currentzone->get_jdim();
				kdim=currentzone->get_kdim();
				num_data=idim*jdim*kdim;
				if (!num_data) ierror=1;
			} else if (currentzone->datatype_point()) {
				if (currentzone->get_num_var_cc()) {
					fatal_err("Datatype must be BLOCK for cell-centered data in clTecplot::read_ascii");
					if (currentzone->_num_varsharing) 
						fatal_err("Varsharing not implemented yet for cell-centered data in clTecplot::read_ascii");
				}
				num_data=currentzone->get_num_node();
			}
			if (!ierror) {
				currentzone->allocate_data();
				if (currentzone->datatype_point()) {
					ddataptr=new double* [_num_var];
					if (!ddataptr) fatal_err(0,"ddataptr");
					for (n_var=0; n_var<_num_var; n_var++) ddataptr[n_var]=currentzone->get_datastartptr(n_var);
					for (n_data=0; n_data<num_data; n_data++) {
						for (n_var=0; n_var<_num_var; n_var++) if (currentzone->get_varsharingzone(n_var)==-1) {
							if (fscanf(fileptr," %le",ddataptr[n_var])!=1) fatal_err(2,filename);
							ddataptr[n_var]++;
						}
					}
					for (n_var=0; n_var<_num_var; n_var++) ddataptr[n_var]=NULL;
					delete [] ddataptr;
					ddataptr=NULL;
				} else {
					for (n_var=0; n_var<_num_var; n_var++) if (currentzone->get_varsharingzone(n_var)==-1) {
						if (!currentzone->zonetype_ordered()) {
							if (currentzone->get_varlocation(n_var)) num_data=currentzone->get_num_elem();
							else num_data=currentzone->get_num_node();
						}
						dataptr=currentzone->get_datastartptr(n_var);
						for (n_data=0; n_data<num_data; n_data++)
							if (fscanf(fileptr," %le",dataptr++)!=1) fatal_err(2,filename);
					}
				}
			} else {
				if (!currentzone->datatype_point()) 
					fatal_err("currentzone->_datatype!=POINT for unspecified idim,jdim,kdim in clTecplot::read_ascii");
				newdataset=new struct dataset;
				if (!newdataset) fatal_err(0,"newdataset in clTecplot::read_ascii");
				newdataset->_data.allocate(num,num_var);
				newdataset->_num=0;
				newdataset->_next=NULL;
				firstdataset=newdataset;
				currentdataset=newdataset;
				idim=0;
				jdim=1;
				kdim=1;
				for (;;) {
					dataptr=currentdataset->_data.get_startptr();
					for (n=0; n<num; n++) {
						for (;;) {
							c=getc(fileptr);
							if (c==' ' || c=='\n' || c=='\t' || c==10 || c==13) continue;
							if (c=='#') {
								ungetc(c,fileptr);
								if (!fgets(line,LINE_LENGTH,fileptr)) fatal_err(2,filename);
								continue;
							}
							if (c=='+' || c=='-' || c=='.' || (c>='0' && c<='9')) {
								ungetc(c,fileptr);
								break;
							}
							finishflag=1;
							if (c!=EOF) ungetc(c,fileptr);
							break;
						}
						if (finishflag) break;
						for (n_var=0; n_var<_num_var; n_var++) if (currentzone->get_varsharingzone(n_var)==-1) 
							if (fscanf(fileptr," %le",dataptr++)!=1) fatal_err(2,filename);
						(currentdataset->_num)++;
						idim++;
					}
					if (finishflag) break;
					newdataset=new struct dataset;
					if (!newdataset) fatal_err(0,"newdataset in clTecplot::read_ascii");
					newdataset->_data.allocate(num,num_var);
					newdataset->_num=0;
					newdataset->_next=NULL;
					currentdataset->_next=newdataset;
					currentdataset=newdataset;
				}
				data=new clVector_namespace::clMatrix3D<double> [num_var];
				if (!data) fatal_err(0,"data in clTecplot::read_ascii");
				ddataptr=new double* [num_var];
				if (!ddataptr) fatal_err(0,"ddataptr");
				for (n_var=0; n_var<_num_var; n_var++) {
					data[n_var].allocate(kdim,jdim,idim);
					ddataptr[n_var]=data[n_var].get_startptr();
				}
				currentdataset=firstdataset;
				for (;;) {
					dataptr=currentdataset->_data.get_startptr();
					for (n=0; n<currentdataset->_num; n++) {
						for (n_var=0; n_var<_num_var; n_var++) {
							(*ddataptr[n_var])=(*dataptr++);
							ddataptr[n_var]++;
						}
					}
					newdataset=currentdataset->_next;
					delete currentdataset;
					if (!newdataset) break;
					currentdataset=newdataset;
				}
				currentzone->set_data(data);
				currentzone->set_idim(idim);
				currentzone->set_jdim(kdim);
				currentzone->set_kdim(jdim);
				delete [] data;
				data=NULL;
				delete [] ddataptr;
				ddataptr=NULL;
			}
		}
		if (!currentzone->zonetype_ordered()) {
			const int num_elem=currentzone->get_num_elem();
			const int num_conn=currentzone->get_num_conn();
			if (num_conn<=0) fatal_err("Illegal value of num_conn in clTecplot::read_ascii");
			conn=currentzone->get_connptr();
			conn->allocate(num_elem,num_conn);
			connptr=conn->get_startptr();
			num_data=num_elem*num_conn;
			for (n_data=0; n_data<num_data; n_data++) {
				if (fscanf(fileptr," %d",connptr)!=1) fatal_err(2,filename);
				--(*connptr++);
			}
		}
	}
	if (!_num_zone) fatal_err("_num_zone=0 in clTecplot::read_ascii");

	fclose(fileptr);

	/* Consolidate all zones */
	_zones=new clTecplot_zone* [_num_zone];
	if (!_zones) fatal_err(0,"_zones in clTecplot::read_ascii");
	int n_zone;
	currentzone=_firstzone;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		if (!currentzone) fatal_err("currentzone not allocated in clTecplot::read_ascii");
		_zones[n_zone]=currentzone;
		currentzone=currentzone->_next;
	}

	/* Set zone sharing info */
	set_zonesharing();
};

/*  Get the number of variables from an ascii Tecplot file */
int clTecplot::get_num_var(const char* const filename) const
{
	int c,num_var;
	char keyword[LINE_LENGTH],varname[NAME_LENGTH];

	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) fatal_err(1,filename);
	for (num_var=0; !num_var;) {
		if (readletterstring(fileptr,keyword,LINE_LENGTH)<0)
			fatal_err(2,filename);
		uppercase(keyword);
		if (!strcmp(keyword,"VARIABLES")) {
			for (;;) {
				do {c=getc(fileptr);} while (c=='=' || c==',' || c==' ' || c=='\n' || c=='\t' || c==10 || c==13);
				if (c=='"') {
					ungetc(c,fileptr);
					if (readdqstring(fileptr,varname,NAME_LENGTH)<0) fatal_err(2,filename);
					num_var++;
				} else {
					ungetc(c,fileptr);
					if (readletternumberstring(fileptr,varname,NAME_LENGTH)<0) fatal_err(2,filename);
					if (!strcmp(varname,"ZONE")) {
						fclose(fileptr);
						return num_var;
					}
					num_var++;
				}
			}
		} else if (!strcmp(keyword,"ZONE")) {
			break;
		}
	}
	fclose(fileptr);
	return num_var;
};

/*  Read the header and initialize the title and variable names from an ascii Tecplot file */
void clTecplot::read_fileheader(const char* const filename)
{
	int c,n_var=0;
	char keyword[LINE_LENGTH],line[LINE_LENGTH],name[NAME_LENGTH];

	FILE *fileptr=fopen(filename,"r");
	if (!fileptr) fatal_err(1,filename);
	for (;;) {
		c=getc(fileptr);
		if (c==EOF) fatal_err("EOF before completion of reading header");
		if (c=='=' || c==',' || c==' ' || c=='\n' || c=='\t' || c==10 || c==13) continue;
		if (c=='#') {
			ungetc(c,fileptr);
			if (!fgets(line,LINE_LENGTH,fileptr)) fatal_err(2,filename);
		} else {
			ungetc(c,fileptr);
			if (readletterstring(fileptr,keyword,LINE_LENGTH)<0) fatal_err(2,filename);
			uppercase(keyword);
			if (!strcmp(keyword,"TITLE")) {
				if (readdqstring(fileptr,name,NAME_LENGTH)<0) fatal_err(2,filename);
				set_title(name);
			} else if (!strcmp(keyword,"VARIABLES")) {
				for (n_var=0;;) {
					do {c=getc(fileptr);} while (c=='=' || c==',' || c==' ' || c=='\n' || c=='\t' || c==10 || c==13);
					if (c=='"') {
						ungetc(c,fileptr);
						if (readdqstring(fileptr,name,NAME_LENGTH)<0) fatal_err(2,filename);
					} else {
						ungetc(c,fileptr);
						if (readletternumberstring(fileptr,name,NAME_LENGTH)<0) fatal_err(2,filename);
						if (!strcmp(name,"ZONE")) {
							if (n_var!=_num_var) fatal_err("Unable to read all variable names in inputfile");
							fclose(fileptr);
							return;
						}
					}
					if (n_var>=_num_var) fatal_err(2,filename);
					if (!set_varname(n_var,name)) fatal_err(2,filename);
					n_var++;
				}
			} else if (!strcmp(keyword,"ZONE")) {
				break;
			} else {
				continue;
			}
		}
	}
	if (n_var!=_num_var) fatal_err("Unable to read all variable names in inputfile");
	fclose(fileptr);
};

/* Read the zone header from an ascii Tecplot file */
int clTecplot::read_zoneheader(FILE *fileptr,clTecplot_zone_readheader *zoneheader)
{
	int c,idummy;
	double dummy;
	char keyword[LINE_LENGTH],line[LINE_LENGTH],name[NAME_LENGTH];

	for (;;) {
		if (readletterstring(fileptr,keyword,LINE_LENGTH)<0) return EOF;
		uppercase(keyword);
		if (!strcmp(keyword,"ZONE")) {
			for (;;) {
				c=getc(fileptr);
				if (c==EOF) return c;
				if (c=='=' || c==',' || c==' ' || c=='\n' || c=='\t' || c==10 || c==13) continue;
				if (c=='#') {
					ungetc(c,fileptr);
					if (!fgets(line,LINE_LENGTH,fileptr)) return 0;
				} else if ((c>='a' && c<='z') || (c>='A' && c<='Z')) {
					ungetc(c,fileptr);
					if (readletterstring(fileptr,keyword,LINE_LENGTH)<0) return 0;
					uppercase(keyword);
					do {
						if ((c=getc(fileptr))==EOF) return 0;
					} while (c!='=');
					if (!strcmp(keyword,"T")) {
						if (readdqstring(fileptr,name,NAME_LENGTH)<0) return 0;
						zoneheader->set_title(name);
					} else if (!strcmp(keyword,"I")) {
						if (fscanf(fileptr,"%d",&idummy)!=1) return 0;
						zoneheader->set_idim(idummy);
					} else if (!strcmp(keyword,"J")) {
						if (fscanf(fileptr,"%d",&idummy)!=1) return 0;
						zoneheader->set_jdim(idummy);
					} else if (!strcmp(keyword,"K")) {
						if (fscanf(fileptr,"%d",&idummy)!=1) return 0;
						zoneheader->set_kdim(idummy);
					} else if (!strcmp(keyword,"N") || !strcmp(keyword,"NODES")) {
						if (fscanf(fileptr,"%d",&idummy)!=1) return 0;
						zoneheader->set_num_node(idummy);
					} else if (!strcmp(keyword,"E") || !strcmp(keyword,"ELEMENTS")) {
						if (fscanf(fileptr,"%d",&idummy)!=1) return 0;
						zoneheader->set_num_elem(idummy);
					} else if (!strcmp(keyword,"F") || !strcmp(keyword,"DATAPACKING")) {
						if (readletterstring(fileptr,line,LINE_LENGTH)<0) return 0;
						uppercase(line);
						if (!zoneheader->set_datatype(line)) return 0;
					} else if (!strcmp(keyword,"STRANDID")) {
						if (fscanf(fileptr,"%le",&dummy)!=1) return 0;
					} else if (!strcmp(keyword,"SOLUTIONTIME")) {
						if (fscanf(fileptr,"%le",&dummy)!=1) return 0;
					} else if (!strcmp(keyword,"ZONETYPE")) {
						if (readletterstring(fileptr,line,LINE_LENGTH)<0) return 0;
						uppercase(line);
						if (!zoneheader->set_zonetype_num_conn(line)) return 0;
					} else if (!strcmp(keyword,"VARSHARELIST")) {
						int n_var,n_varstart,n_varend,n_zone;
						clVector_namespace::clVector<int> varsharing(_num_var,0);
						do {
							if ((c=getc(fileptr))==EOF) return 0;
						} while (c!='(');
						do {
							if ((c=getc(fileptr))==EOF) return 0;
						} while (c!='[');
						for (;;) {
							if (fscanf(fileptr,"%d",&n_varstart)!=1) return 0;
							if (n_varstart<1 || n_varstart>_num_var) return 0;
							do {
								if ((c=getc(fileptr))==EOF) return 0;
							} while (c!='-' && c!=',' && c!=']');
							if (c=='-') {
								if (fscanf(fileptr,"%d",&n_varend)!=1) return 0;
								if (n_varend<=n_varstart || n_varend>_num_var) return 0;
								n_varstart--;
								for (n_var=n_varstart; n_var<n_varend; n_var++) {
									if (varsharing(n_var)) return 0;
									varsharing(n_var)++;
								}
								zoneheader->_num_varsharing+=n_varend-n_varstart;
								do {
									if ((c=getc(fileptr))==EOF) return 0;
								} while (c!=',' && c!=']');
							} else {
								n_varstart--;
								if (varsharing(n_varstart)) return 0;
								varsharing(n_varstart)++;
								(zoneheader->_num_varsharing)++;
							}
							if (c==']') break;
						}
						do {
							if ((c=getc(fileptr))==EOF) return 0;
						} while (c!='=');
						if (fscanf(fileptr,"%d",&n_zone)!=1) return 0;
						n_zone--;
						if (n_zone<0 || n_zone>=_num_zone) return 0;
						for (n_var=0; n_var<_num_var; n_var++) 
							if (varsharing(n_var)) zoneheader->set_varsharingzone(n_var,n_zone);
						do {
							if ((c=getc(fileptr))==EOF) return 0;
						} while (c!=',' && c!=')');
					} else if (!strcmp(keyword,"VARLOCATION")) {
						int n_var,nn_var,n_varstart,n_varend;
						clVector_namespace::clVector<int> varlocation(_num_var);
						do {
							if ((c=getc(fileptr))==EOF) return 0;
						} while (c!='(');
						for (;;) {
							varlocation=0;
							do {
								if ((c=getc(fileptr))==EOF) return 0;
							} while (c!='[');
							for (;;) {
								if (fscanf(fileptr,"%d",&n_varstart)!=1) return 0;
								if (n_varstart<1 || n_varstart>_num_var) return 0;
								do {
									if ((c=getc(fileptr))==EOF) return 0;
								} while (c!='-' && c!=',' && c!=']');
								if (c=='-') {
									if (fscanf(fileptr,"%d",&n_varend)!=1) return 0;
									if (n_varend<=n_varstart || n_varend>_num_var) return 0;
									n_varstart--;
									for (n_var=n_varstart; n_var<n_varend; n_var++) {
										if (varlocation(n_var)) return 0;
										varlocation(n_var)++;
									}
									do {
										if ((c=getc(fileptr))==EOF) return 0;
									} while (c!=',' && c!=']');
								} else {
									n_varstart--;
									if (varlocation(n_varstart)) return 0;
									varlocation(n_varstart)++;
								}
								if (c==']') break;
							}
							do {
								if ((c=getc(fileptr))==EOF) return 0;
							} while (c!='=');
							if (readletterstring(fileptr,line,LINE_LENGTH)<0) return 0;
							uppercase(line);
							if (!strcmp(line,"CELLCENTERED")) {
								for (n_var=0,nn_var=0; n_var<_num_var; n_var++) if (varlocation(n_var)) {
									zoneheader->set_varlocation(n_var);
									nn_var++;
								}
								zoneheader->incr_num_var_cc(nn_var);
							} else if (strcmp(line,"NODAL")) fatal_err("Illegal VARLOCATION value in input file");
							do {
								if ((c=getc(fileptr))==EOF) return 0;
							} while (c!=',' && c!=')');
							if (c==')') break;
						}
					} else if (!strcmp(keyword,"DT")) {
						if (!fgets(line,LINE_LENGTH,fileptr)) return 0;
					}
				} else {
					ungetc(c,fileptr);
					return 1;
				}
			}
		}
	}
};

/* Write an ascii Tecplot file */
void clTecplot::write_ascii(const char* const filename,const char* const datatype) const
{
//RVRV TODO not working after loading a mixed nodal/face removing the face vars and printing as point
	int n_zone;
	clVector_namespace::clVector<int> swapdatatype(_num_zone,0);
	clTecplot_zone *currentzone=NULL;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (!currentzone) fatal_err("currentzone not allocated in clTecplot::write_ascii");
		if (!currentzone->zonetype_ordered() && currentzone->get_num_var_cc()) continue;
		if (!strcmp(datatype,"POINT")) {
			if (!(currentzone->datatype_point())) swapdatatype(n_zone)++;
			currentzone->set_datatype(POINT);
		} else {
			if (currentzone->datatype_point()) swapdatatype(n_zone)++;
			currentzone->set_datatype(BLOCK);
		}
	}
	write_ascii(filename);
	/* Restore datatype */
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (!currentzone) fatal_err("currentzone not allocated in clTecplot::write_ascii");
		if (!currentzone->zonetype_ordered() && currentzone->get_num_var_cc()) continue;
		if (swapdatatype(n_zone)) {
			if (currentzone->datatype_point()) currentzone->set_datatype(BLOCK);
			else currentzone->set_datatype(POINT);
		}
	}
};
void clTecplot::write_ascii(const char* const filename) const
{
	clVector_namespace::clVector<int> writezone(_num_zone,1);
	clVector_namespace::clVector<int> writevar(_num_var,1);
	write_ascii(writezone,filename,writevar);
};
void clTecplot::write_ascii(const char* const filename,const int n_varlow,const int n_varupp) const
{
	clVector_namespace::clVector<int> writezone(_num_zone,1);
	clVector_namespace::clVector<int> writevar(_num_var,0);
	if (n_varlow<0 || n_varlow>=_num_var) fatal_err("Illegal variable index in clTecplot::write_ascii");
	writevar(n_varlow)=1;
	for (int n_var=n_varlow+1; n_var<=n_varupp; n_var++) writevar(n_var)=1;
	write_ascii(writezone,filename,writevar);
};
void clTecplot::write_ascii(const char* const filename,const clVector_namespace::clVector<int> &writevar) const
{
	clVector_namespace::clVector<int> writezone(_num_zone,1);
	write_ascii(writezone,filename,writevar);
};
void clTecplot::write_ascii(const char* const filename,const clVector_namespace::clVector<int> &writevar,const char* const datatype) const
{
	clTecplot_zone *currentzone=NULL;
	for (int n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (!currentzone) fatal_err("currentzone not allocated in clTecplot::write_ascii");
		if (!currentzone->zonetype_ordered() && currentzone->get_num_var_cc()) continue;
		if (!strcmp(datatype,"BLOCK")) currentzone->set_datatype(BLOCK);
		else currentzone->set_datatype(POINT);
	}
	clVector_namespace::clVector<int> writezone(_num_zone,1);
	write_ascii(writezone,filename,writevar);
};
void clTecplot::write_ascii(const int n_zone_target,const char* const filename) const
{
	write_ascii(n_zone_target,n_zone_target,filename);
};
void clTecplot::write_ascii(const int n_zonelow,const int n_zoneupp,const char* const filename) const
{
	clVector_namespace::clVector<int> writezone(_num_zone,0);
	clVector_namespace::clVector<int> writevar(_num_var,1);
	if (n_zonelow<0 || n_zonelow>=_num_zone) fatal_err("Illegal zone index in clTecplot::write_ascii");
	writezone(n_zonelow)=1;
	for (int n_zone=n_zonelow+1; n_zone<=n_zoneupp; n_zone++) writezone(n_zone)=1;
	write_ascii(writezone,filename,writevar);
};
void clTecplot::write_ascii(const int n_zonelow,const int n_zoneupp,const char* const filename,const int n_varlow,const int n_varupp) const
{
	clVector_namespace::clVector<int> writezone(_num_zone,0);
	clVector_namespace::clVector<int> writevar(_num_var,0);
	if (n_zonelow<0 || n_zonelow>=_num_zone) fatal_err("Illegal zone index in clTecplot::write_ascii");
	if (n_varlow<0 || n_varlow>=_num_var) fatal_err("Illegal variable index in clTecplot::write_ascii");
	writezone(n_zonelow)=1;
	for (int n_zone=n_zonelow+1; n_zone<=n_zoneupp; n_zone++) writezone(n_zone)=1;
	writevar(n_varlow)=1;
	for (int n_var=n_varlow+1; n_var<=n_varupp; n_var++) writevar(n_var)=1;
	write_ascii(writezone,filename,writevar);
};
void clTecplot::write_ascii(const clVector_namespace::clVector<int> &writezone,const char* const filename) const
{
	clVector_namespace::clVector<int> writevar(_num_var,1);
	write_ascii(writezone,filename,writevar);
};
void clTecplot::write_ascii(const clVector_namespace::clVector<int> &writezone,const char* const filename,
	const clVector_namespace::clVector<int> &writevar) const
{
	int n_data,num_data,n_var,num_var,num_varwrite,n_zone;
	const char *varname=NULL;

	clTecplot_zone *currentzone=NULL;

	for (n_var=0,num_varwrite=0; n_var<_num_var; n_var++) if (writevar(n_var)) num_varwrite++;

	FILE *fileptr=fopen(filename,"w");
	if (!fileptr) fatal_err(1,filename);

	if (if_title()) fprintf(fileptr,"TITLE = \"%s\"\n",get_title());
	fprintf(fileptr,"VARIABLES =");
	for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
		varname=get_varname(n_var);
		if (varname) {
			fprintf(fileptr," \"%s\"",varname);
		} else {
			fprintf(fileptr," \"V%d\"",n_var);
		}
	}
	fprintf(fileptr,"\n");

	/* Check for varsharing */
	int shareflag=0;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		if (!writezone(n_zone) && _num_zonesharing) {
			if (_zonesharing(n_zone)) shareflag++;
		}
	}
	if (shareflag) fatal_err("Cannot remove shared zone from writing in clTecplot::write_ascii");

	int idim=0,jdim=0,kdim=0,num_node=0,num_elem=0;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		if (!writezone(n_zone)) continue;
		currentzone=_zones[n_zone];
		num_var=_num_var;
		if (!currentzone->get_num_var_cc()) num_var-=currentzone->_num_varsharing;
		if (currentzone->zonetype_ordered()) {
			idim=currentzone->get_idim();
			jdim=currentzone->get_jdim();
			kdim=currentzone->get_kdim();
			if (currentzone->if_title())
				fprintf(fileptr,"ZONE T=\"%s\"\nI=%d, J=%d, K=%d, ",
					currentzone->get_title(),idim,jdim,kdim);
			else
				fprintf(fileptr,"ZONE\nI=%d, J=%d, K=%d, ",idim,jdim,kdim);
		} else {
			num_node=currentzone->get_num_node();
			num_elem=currentzone->get_num_elem();
			if (currentzone->if_title())
				fprintf(fileptr,"ZONE T=\"%s\" N=%d, E=%d, ",
					currentzone->get_title(),num_node,num_elem);
			else
				fprintf(fileptr,"ZONE N=%d, E=%d, ",num_node,num_elem);
			if (num_varwrite>1) currentzone->print_zonetype(fileptr);
		}
		currentzone->write_datatype(fileptr);
		if (currentzone->_num_varsharing) {
			if (_num_zonesharing>1) fatal_err("Variable sharing from more than one zone not yet implemented in clTecplot::write_ascii");
			int nn_var=0;
			for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
				if (currentzone->_varsharingzone(n_var)!=-1) nn_var++;
			}
			if (nn_var) {
				clVector_namespace::clVector<int> varsharingzone(num_varwrite);
				for (n_var=0,nn_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
					varsharingzone(nn_var++)=currentzone->_varsharingzone(n_var);
				}
				fprintf(fileptr,"VARSHARELIST=([");
				int printflag=0,rangeflag=0,sharingzone=-1;
				for (n_var=0; n_var<num_varwrite; n_var++) {
					if (varsharingzone(n_var)!=-1) {
						if (sharingzone<0) sharingzone=varsharingzone(n_var);
						else if (sharingzone!=varsharingzone(n_var)) fatal_err("Inconsitent variable sharing zone in clTecplot::write_ascii");
						if (!printflag) {
							printflag=1;
							fprintf(fileptr,"%d",n_var+1);
							nn_var=n_var;
						} else {
							if (n_var==nn_var+1) {
								rangeflag=1;
							} else {
								if (rangeflag) fprintf(fileptr,"-%d,%d",nn_var+1,n_var+1);
								else fprintf(fileptr,",%d",n_var+1);
								rangeflag=0;
							}
							nn_var=n_var;
						}
					}
				}
				if (!printflag) fatal_err("Unable to find nodal variable indices in clTecplot::write_ascii");
				if (rangeflag) fprintf(fileptr,"-%d",nn_var+1);
				fprintf(fileptr,"]=%d)\n",sharingzone+1);
			}
		}
		if (currentzone->get_num_var_cc()) {
			int nn_var=0;
			for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
				if (currentzone->get_varlocation(n_var)) nn_var++;
			}
			if (nn_var) {
				clVector_namespace::clVector<int> varlocation(num_varwrite);
				for (n_var=0,nn_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
					varlocation(nn_var++)=currentzone->get_varlocation(n_var);
				}
				fprintf(fileptr,"VARLOCATION=([");
				int printflagnodal=0,rangeflag=0;
				for (n_var=0; n_var<num_varwrite; n_var++) {
					if (!varlocation(n_var)) {
						if (!printflagnodal) {
							printflagnodal=1;
							fprintf(fileptr,"%d",n_var+1);
							nn_var=n_var;
						} else {
							if (n_var==nn_var+1) {
								rangeflag=1;
							} else {
								if (rangeflag) fprintf(fileptr,"-%d,%d",nn_var+1,n_var+1);
								else fprintf(fileptr,",%d",n_var+1);
								rangeflag=0;
							}
							nn_var=n_var;
						}
					}
				}
				if (printflagnodal) {
					if (rangeflag) fprintf(fileptr,"-%d",nn_var+1);
					fprintf(fileptr,"]=NODAL");
				}
				int printflagcc=0;
				rangeflag=0;
				for (n_var=0; n_var<num_varwrite; n_var++) {
					if (varlocation(n_var)) {
						if (!printflagcc) {
							printflagcc=1;
							if (printflagnodal) fprintf(fileptr,", [%d",n_var+1);
							else fprintf(fileptr,"%d",n_var+1);
							nn_var=n_var;
						} else {
							if (n_var==nn_var+1) {
								rangeflag=1;
							} else {
								if (rangeflag) fprintf(fileptr,"-%d,%d",nn_var+1,n_var+1);
								else fprintf(fileptr,",%d",n_var+1);
								rangeflag=0;
							}
							nn_var=n_var;
						}
					}
				}
				if (printflagcc) {
					if (rangeflag) fprintf(fileptr,"-%d",nn_var+1);
					fprintf(fileptr,"]=CELLCENTERED)\n");
				} else if (printflagnodal) {
					fprintf(fileptr,")\n");
				} else {
					fatal_err("Unable to find nodal variable indices in clTecplot::write_ascii");
				}
			}
		}
		if (num_var) {
			if (num_varwrite==1 || currentzone->datatype_point()) {
				if (num_varwrite>1 && currentzone->get_num_var_cc()) 
					fatal_err("Cannot write cell-centered variables in POINT format");
				if (currentzone->zonetype_ordered()) {
					num_data=idim*jdim*kdim;
				} else {
					num_data=num_node;
				}
				const double **dataptr=new const double* [_num_var];
				if (!dataptr) fatal_err(0,"dataptr");
				for (n_var=0; n_var<_num_var; n_var++) dataptr[n_var]=currentzone->get_datastartptr(n_var);
				if (!currentzone->_num_varsharing) {
					for (n_data=0; n_data<num_data; n_data++) {
						for (n_var=0; n_var<_num_var; n_var++) {
							if (writevar(n_var)) {
								fprintf(fileptr," %.9le",*(dataptr[n_var]));
								dataptr[n_var]++;
							}
						}
						fprintf(fileptr,"\n");
					}
				} else {
					int num_varsharewrite=0;
					for (n_var=0; n_var<_num_var; n_var++) 
						if (writevar(n_var) && currentzone->_varsharingzone(n_var)==-1) num_varsharewrite++;
					for (n_data=0; n_data<num_data; n_data++) {
						for (n_var=0; n_var<_num_var; n_var++) if (currentzone->_varsharingzone(n_var)==-1) {
							if (writevar(n_var)) {
								fprintf(fileptr," %.9le",*(dataptr[n_var]));
								dataptr[n_var]++;
							}
						}
						if (num_varsharewrite) fprintf(fileptr,"\n");
					}
				}
				for (n_var=0; n_var<_num_var; n_var++) dataptr[n_var]=NULL;
				delete [] dataptr;
			} else {
				const double *dataptr;
				if (!currentzone->_num_varsharing) {
					for (n_var=0; n_var<_num_var; n_var++) {
						if (writevar(n_var))	{
							if (currentzone->zonetype_ordered()) {
								num_data=idim*jdim*kdim;
							} else {
								if (currentzone->get_varlocation(n_var)) num_data=num_elem;
								else num_data=num_node;
							}
							dataptr=currentzone->get_datastartptr(n_var);
							for (n_data=0; n_data<num_data;) {
								fprintf(fileptr," %.9le",*dataptr++);
								if (++n_data%6==0) fprintf(fileptr,"\n");
							}
							if (n_data%6) fprintf(fileptr,"\n");
						}
					}
				} else {
					for (n_var=0; n_var<_num_var; n_var++) if (currentzone->_varsharingzone(n_var)==-1) {
						if (writevar(n_var))	{
							if (currentzone->zonetype_ordered()) {
								num_data=idim*jdim*kdim;
							} else {
								if (currentzone->get_varlocation(n_var)) num_data=num_elem;
								else num_data=num_node;
							}
							dataptr=currentzone->get_datastartptr(n_var);
							for (n_data=0; n_data<num_data;) {
								fprintf(fileptr," %.9le",*dataptr++);
								if (++n_data%6==0) fprintf(fileptr,"\n");
							}
							if (n_data%6) fprintf(fileptr,"\n");
						}
					}
				}
			}
		}
		if (!currentzone->zonetype_ordered()) {
			int n_elem,num_elem=currentzone->get_num_elem();
			int n_conn,nn_conn,nnn_conn,num_conn=currentzone->get_num_conn();
			if (num_elem>0) {
				const int *connptr=currentzone->get_connstartptr();
				if (num_conn==3) {
					for (n_elem=0; n_elem<num_elem; n_elem++) {
						for (n_conn=0; n_conn<num_conn; n_conn++) fprintf(fileptr," %d",(*connptr++)+1);
						fprintf(fileptr,"\n");
					}
				} else if (num_conn==4) {
					for (n_elem=0; n_elem<num_elem; n_elem++) {
						for (n_conn=0; n_conn<num_conn; n_conn++) {
							nn_conn=(*connptr++);
							if (nn_conn>=0) {
								fprintf(fileptr," %d",nn_conn+1);
							} else {
								if (!n_conn) fprintf(fileptr," %d",(*connptr)+1);
								else fprintf(fileptr," %d",nnn_conn+1);
							}
							nnn_conn=nn_conn;
						}
						fprintf(fileptr,"\n");
					}
				} else fatal_err("Illegal value of num_conn in clTecplot::write_ascii");
			}
		}
	}

	fclose(fileptr);
};

/* Read a binary Tecplot file */
void clTecplot::read_bin(const char* const filename)
{
	using clUtils_namespace::read_binentry;
	using clUtils_namespace::read_binentry_int;

	if (clUtils_namespace::_size_my_int!=4) fatal_err("Incompatible int32 type");
	if (clUtils_namespace::_size_my_float!=4) fatal_err("Incompatible float type");
	if (clUtils_namespace::_size_my_double!=8) fatal_err("Incompatible double type");

	int idummy;
	clUtils_namespace::my_float readfloat;
	char readchar;
	char *readversion=new char [4];

	FILE *fileptr=fopen(filename,"rb");
	if (!fileptr) fatal_err(1,filename);

	/*  Magic number  */
	void *io_ptr=(void *)(&readchar);
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"Magic number=%c",readchar);
	if (readchar!='#') fatal_err("Illegal Magic number in binary inputfile");
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"%c",readchar);
	if (readchar!='!') fatal_err("Illegal Magic number in binary inputfile");
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"%c",readchar);
	if (readchar!='T') fatal_err("Illegal Magic number in binary inputfile");
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"%c",readchar);
	if (readchar!='D') fatal_err("Illegal Magic number in binary inputfile");
	if (fread(io_ptr,sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	if (read_bin_debug) fprintf(stdout,"%c\n",readchar);
	if (readchar!='V') fatal_err("Illegal Magic number in binary inputfile");
	io_ptr=(void *)readversion;
	if (fread(io_ptr,3*sizeof(char),1,fileptr)!=1) fatal_err(2,filename);
	readversion[3]='\0';
	_version=atoi(readversion);
	if (read_bin_debug) fprintf(stdout,"Version=%d\n",_version);
	if (_version!=102 && _version!=112) {
		fprintf(stdout,"readversion=%s\n",readversion);
		fatal_err("Tecplot version of inputfile not implemented");
	}	
	delete [] readversion;
	int ierror=read_binentry_int(fileptr,idummy,1);
	if (!ierror) fatal_err(2,filename);
	else if (ierror==-1) fatal_err("Big/little endian incompatibility");

	if (!read_title(fileptr)) fatal_err(2,filename);

	clTecplot_zone *currentzone=NULL;
	for (;;) {
		if (!read_binentry(fileptr,readfloat)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Marker=%f\n",readfloat);
		if (readfloat==299.0) {
			if (!read_zone(fileptr,&currentzone)) fatal_err(2,filename);
		} else if (readfloat==399.0) read_geometries(fileptr);
		else if (readfloat==499.0) read_text(fileptr);
		else if (readfloat==599.0) read_customlabel(fileptr);
		else if (readfloat==699.0) read_userrec(fileptr);
		else if (readfloat==799.0) {
			if (!read_dataauxdata(fileptr)) fatal_err(2,filename);
		} else if (readfloat==899.0) read_varauxdata(fileptr);
		else if (readfloat==357.0) break;
		else fatal_err("Illegal marker");
	}

	_zones=new clTecplot_zone* [_num_zone];
	if (!_zones) fatal_err(0,"_zones in clTecplot::read_bin");
	int n_zone;
	currentzone=_firstzone;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		if (!currentzone) 
			fatal_err("currentzone not allocated in clTecplot::read_bin");
		if (!read_binentry(fileptr,readfloat)) fatal_err(2,filename);
		if (read_bin_debug) fprintf(stdout,"Marker=%f\n",readfloat);
		if (readfloat!=299.0) fatal_err("Missing zone marker");
		if (!read_datasection(fileptr,currentzone)) fatal_err(2,filename);
		_zones[n_zone]=currentzone;
		currentzone=currentzone->_next;
	}

	/* Set zone sharing info */
	set_zonesharing();

	fclose(fileptr);
};

/* Read the title from a binary Tecplot file */
int clTecplot::read_title(FILE *fileptr)
{
	using clUtils_namespace::read_binentry_int;

	int i,idummy;
	char name[NAME_LENGTH];

	/*  Filetype  */
	if (_version==112) if (!read_binentry_int(fileptr,idummy)) return 0;
	if (idummy) fatal_err("Only full file format (grid+solution) is currently implemented");

	/*  Title  */
	for (i=0; i<NAME_LENGTH; i++) {
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (!i) strcpy(name,(char *)(&idummy));
		else strcat(name,(char *)(&idummy));
		if (!idummy) break;
	}
	if (read_bin_debug) fprintf(stdout,"Title=%s\n",name);
	set_title(name);

	/*  Number of variables  */
	if (!read_binentry_int(fileptr,_num_var)) return 0;
	if (read_bin_debug) fprintf(stdout,"Number of variables=%d\n",_num_var);

	/*  Allocate arrays for coordinate and variable names  */
	if (_num_var) _varname.allocate(_num_var);
	else fatal_err("Unable to obtain the number of variables in clTecplot::read_title");

	/*  Number of variables and variable names  */
	for (int n_var=0; n_var<_num_var; n_var++) {
		for (i=0; i<NAME_LENGTH; i++) {
			if (!read_binentry_int(fileptr,idummy)) return 0;
			if (!i) strcpy(name,(char *)(&idummy));
			else strcat(name,(char *)(&idummy));
			if (!idummy) break;
		}
		if (!set_varname(n_var,name)) return 0;
		if (read_bin_debug) fprintf(stdout,"Name variable %d=%s\n",n_var,name);
	}

	return 1;
};

/* Read a zone from a binary Tecplot file */
int clTecplot::read_zone(FILE *fileptr,clTecplot_zone **currentzone)
{
	using clUtils_namespace::read_binentry;
	using clUtils_namespace::read_binentry_int;

	int i,idummy;
	clUtils_namespace::my_double readdouble;

	/*  Zone title  */
	char title[NAME_LENGTH];
	for (i=0; i<NAME_LENGTH; i++) {
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (!i) strcpy(title,(char *)(&idummy));
		else strcat(title,(char *)(&idummy));
		if (!idummy) break;
	}
	if (read_bin_debug) fprintf(stdout,"Zone title=%s\n",title);
	int zonetype=-1;
	if (_version==102) {
		/*  Not used, set to -1  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy!=-1) fatal_err("Missing entry in Zones");
		/*  ZoneType  */
		if (!read_binentry_int(fileptr,zonetype)) return 0;
		if (read_bin_debug) fprintf(stdout,"Zonetype=%d\n",zonetype);
		/*  Solution time  */
		if (!read_binentry(fileptr,readdouble)) return 0;
		if (read_bin_debug) fprintf(stdout,"Solution time=%f\n",readdouble);
		/*  StrandID  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"StrandID=%d\n",idummy);
	} else if (_version==112) {
		/*  ParentZone  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"ParentZone=%d\n",idummy);
		/*  StrandID  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"StrandID=%d\n",idummy);
		/*  Solution time  */
		if (!read_binentry(fileptr,readdouble)) return 0;
		if (read_bin_debug) fprintf(stdout,"Solution time=%f\n",readdouble);
		/*  Not used, set to -1 */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy!=-1) fatal_err("Illegal entry in Zones");
		/*  ZoneType  */
		if (!read_binentry_int(fileptr,zonetype)) return 0;
		if (read_bin_debug) fprintf(stdout,"ZoneType=%d\n",zonetype);
	}

	/*  Allocate and add new zone  */
	clTecplot_zone *newzone=NULL;
	if (!zonetype) {
		newzone=new clTecplot_zone_structured(_num_var);
		if (!newzone) fatal_err(0,"newzone in clTecplot::read_zone");
	} else {
		newzone=new clTecplot_zone_unstructured(_num_var);
		if (!newzone) fatal_err(0,"newzone in clTecplot::read_zone");
	}
	if (!newzone->set_zonetype_num_conn(zonetype)) return 0;
	set_title(title);

	/*  Var location  */
	if (_version==112) {
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy) {
			const clVector_namespace::clVector<int> varlocation(_num_var);
			for (int n_var=0; n_var<_num_var; n_var++)
				if (!read_binentry_int(fileptr,varlocation(n_var))) return 0;
			newzone->set_varlocation(varlocation);
		}
		/*  Are 1-to-1 face neighbors supplied  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		/*  Number of miscellaneous user defined neighbor connections  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy) {
			if (!read_binentry_int(fileptr,idummy)) return 0;
			if (zonetype) if (!read_binentry_int(fileptr,idummy)) return 0;
		}
	}

	/*  Zone dimensions  */
	if (newzone->zonetype_ordered()) {
		/*  I,J,K dimensions  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"I=%d\n",idummy);
		newzone->set_idim(idummy);
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"J=%d\n",idummy);
		newzone->set_jdim(idummy);
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"K=%d\n",idummy);
		newzone->set_kdim(idummy);
	} else {
		/*  Number of points  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"N=%d\n",idummy);
		newzone->set_num_node(idummy);
		/*  Faces info  */
		if (newzone->_zonetype==FEPOLYGON || newzone->_zonetype==FEPOLYHEDRON) {
			int numfaces;
			if (!read_binentry_int(fileptr,numfaces)) return 0;
			if (read_bin_debug) fprintf(stdout,"Numfaces=%d\n",numfaces);
			if (!read_binentry_int(fileptr,idummy)) return 0;
			if (newzone->_zonetype==FEPOLYGON && idummy!=2*numfaces) 
				fatal_err("Number of face nodes must be twice the number of faces for FEPOLYGON");
			if (!read_binentry_int(fileptr,idummy)) return 0;
			if (!read_binentry_int(fileptr,idummy)) return 0;
		}
		/*  Number of elements  */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (read_bin_debug) fprintf(stdout,"E=%d\n",idummy);
		newzone->set_num_elem(idummy);
		/* For future use */
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy!=0) fatal_err("Illegal ICellDim");
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy!=0) fatal_err("Illegal JCellDim");
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy!=0) fatal_err("Illegal KCellDim");
	}
	/*  Auxiliary data name/values  */
	if (!read_binentry_int(fileptr,idummy)) return 0;
	for (;;) {
		if (!idummy) break;
		for (i=0; i<NAME_LENGTH; i++) {
			if (!read_binentry_int(fileptr,idummy)) return 0;
			if (!idummy) break;
		}
		if (!read_binentry_int(fileptr,idummy)) return 0;
		for (i=0; i<NAME_LENGTH; i++) {
			if (!read_binentry_int(fileptr,idummy)) return 0;
			if (!idummy) break;
		}
		if (!read_binentry_int(fileptr,idummy)) return 0;
	}

	/* Binary files are always in datatype BLOCK */
	newzone->set_datatype(BLOCK);
	
	if (!_num_zone) _firstzone=newzone;
	else (*currentzone)->_next=newzone;
	(*currentzone)=newzone;
	_num_zone++;

	return 1;
};

/* Read the geometries from a binary Tecplot file */
void clTecplot::read_geometries(FILE *fileptr)
{
	fatal_err("read_geometries not yet implemented");
};

/* Read text from a binary Tecplot file */
void clTecplot::read_text(FILE *fileptr)
{
	fatal_err("read_text not yet implemented");
};

/* Read a custom label from a binary Tecplot file */
void clTecplot::read_customlabel(FILE *fileptr)
{
	fatal_err("read_customlabel not yet implemented");
};

/* Read a user record from a binary Tecplot file */
void clTecplot::read_userrec(FILE *fileptr)
{
	fatal_err("read_userrec not yet implemented");
};

/* Read auxiliary data from a binary Tecplot file */
int clTecplot::read_dataauxdata(FILE *fileptr)
{
	using clUtils_namespace::read_binentry_int;

	int i,idummy;
	char auxname[NAME_LENGTH],auxtext[NAME_LENGTH];

	for (i=0; i<NAME_LENGTH; i++) {
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (!i) strcpy(auxname,(char *)(&idummy));
		else strcat(auxname,(char *)(&idummy));
		if (!idummy) break;
	}
	if (read_bin_debug) fprintf(stdout,"Auxiliary variable name = %s\n",auxname);
	if (!read_binentry_int(fileptr,idummy)) return 0;
	if (idummy) fatal_err("Illegal value for auxiliary variable data format");
	for (i=0; i<NAME_LENGTH; i++) {
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (!i) strcpy(auxtext,(char *)(&idummy));
		else strcat(auxtext,(char *)(&idummy));
		if (!idummy) break;
	}
	if (read_bin_debug) fprintf(stdout,"Auxiliary variable text = %s\n",auxtext);
	return 1;
};

/* Read auxiliary variable data from a binary Tecplot file */
void clTecplot::read_varauxdata(FILE *fileptr)
{
	fatal_err("read_varauxdata not yet implemented");
};

/* Read data section from a binary Tecplot file */
int clTecplot::read_datasection(FILE *fileptr,clTecplot_zone *currentzone)
{
	using clUtils_namespace::read_binentry;
	using clUtils_namespace::read_binarray;
	using clUtils_namespace::read_binentry_int;

	int *readintptr=NULL;
	clVector_namespace::clVector<int> readintbuffer(IOBUFFERSIZE);
	clUtils_namespace::my_float readfloat;
	clUtils_namespace::my_float *readfloatptr=NULL;
	clVector_namespace::clVector<clUtils_namespace::my_float> readfloatbuffer(IOBUFFERSIZE);
	clUtils_namespace::my_double readdouble;
	clUtils_namespace::my_double *readdoubleptr=NULL;
	clVector_namespace::clVector<clUtils_namespace::my_double> readdoublebuffer(IOBUFFERSIZE);

	int idummy,n,n_var,num_var,n_data,num_data,nbuffer,numbuffer,buffersize;
	int *connptr=NULL;
	double *dataptr=NULL;
	clVector_namespace::clVector<int> dataformat(_num_var,0);

	/*  Data type  */
	for (n_var=0; n_var<_num_var; n_var++) {
		if (!read_binentry_int(fileptr,dataformat(n_var))) return 0;
		if (dataformat(n_var)!=1 && dataformat(n_var)!=2) fatal_err("Data format not implemented");
		if (read_bin_debug) fprintf(stdout,"Datatype=%d\n",dataformat(n_var));
	}
	/*  Passive variables  */
	if (_version==112) {
		if (!read_binentry_int(fileptr,idummy)) return 0;
		if (idummy) for (n_var=0; n_var<_num_var; n_var++) {
			fatal_err("Passive variables not yet implemented");
			if (!read_binentry_int(fileptr,idummy)) return 0;
		}
	}
	/*  Variable sharing  */
	currentzone->_num_varsharing=0;
	for (n_var=0; n_var<_num_var; n_var++) currentzone->set_varsharingzone(n_var,-1);
	if (!read_binentry_int(fileptr,idummy)) return 0;
	if (idummy) {
		if (!currentzone->_zonetype) fatal_err("Variable sharing not yet implemented for ordered zones");
		for (n_var=0; n_var<_num_var; n_var++) {
			if (!read_binentry_int(fileptr,idummy)) return 0;
			currentzone->set_varsharingzone(n_var,idummy);
		}
	}
	if (read_bin_debug) {
		fprintf(stdout,"Number of variable sharing=%d\n",currentzone->_num_varsharing);
		if (currentzone->_num_varsharing) for (n_var=0; n_var<_num_var; n_var++) 
			fprintf(stdout,"Variable sharing zone (%d)=%d\n",n_var,currentzone->_varsharingzone(n_var));
	}
	if (!read_binentry_int(fileptr,idummy)) return 0;
	if (idummy!=-1) fatal_err("Connectivities sharing not yet implemented");
	num_var=_num_var-currentzone->_num_varsharing;
	/*  Varmin, varmax */
	if (_version==112) {
		currentzone->_varmin.allocate(_num_var);
		currentzone->_varmax.allocate(_num_var);
		for (n_var=0; n_var<_num_var; n_var++) if (currentzone->get_varsharingzone(n_var)==-1) {
			if (!read_binentry(fileptr,readdouble)) return 0;
			currentzone->_varmin(n_var)=(double)(readdouble);
			if (read_bin_debug) fprintf(stdout,"Varmin[%d]=%f\n",n_var,currentzone->_varmin(n_var));
			if (!read_binentry(fileptr,readdouble)) return 0;
			currentzone->_varmax(n_var)=(double)(readdouble);
			if (read_bin_debug) fprintf(stdout,"Varmax[%d]=%f\n",n_var,currentzone->_varmax(n_var));
		}
	}
	/* Data */
	if (num_var) {
		if (currentzone->datatype_point()) fatal_err("datatype should always be BLOCK for binary files");
		if (!currentzone->_zonetype) {
			const int idim=currentzone->get_idim();
			const int jdim=currentzone->get_jdim();
			const int kdim=currentzone->get_kdim();
			num_data=idim*jdim*kdim;
			numbuffer=(int)((num_data-0.5)/IOBUFFERSIZE)+1;
		}
		currentzone->allocate_data();
		for (n_var=0; n_var<_num_var; n_var++) if (currentzone->get_varsharingzone(n_var)==-1) {
			if (!currentzone->zonetype_ordered()) {
				if (currentzone->get_varlocation(n_var)) num_data=currentzone->get_num_elem();
				else num_data=currentzone->get_num_node();
				numbuffer=(int)((num_data-0.5)/IOBUFFERSIZE)+1;
			}
			dataptr=currentzone->get_datastartptr(n_var);
			buffersize=IOBUFFERSIZE;
			for (nbuffer=0; nbuffer<numbuffer; nbuffer++) {
				if (nbuffer==numbuffer-1) buffersize=num_data-nbuffer*IOBUFFERSIZE;
				if (dataformat(n_var)==1) {
					if (!read_binarray(fileptr,readfloatbuffer.get_startptr(),buffersize)) return 0;
					readfloatptr=readfloatbuffer.get_startptr();
					for (n=0; n<buffersize; n++) (*dataptr++)=(double)(*readfloatptr++);
				} else if (dataformat(n_var)==2) {
					if (!read_binarray(fileptr,readdoublebuffer.get_startptr(),buffersize)) return 0;
					readdoubleptr=readdoublebuffer.get_startptr();
					for (n=0; n<buffersize; n++) (*dataptr++)=(double)(*readdoubleptr++);
				}
			}
		}
	}
	if (currentzone->_zonetype) {
		if (currentzone->_zonetype==FEPOLYGON || currentzone->_zonetype==FEPOLYHEDRON) 
			fatal_err("FEPOLYGON and FEPOLYHEDRON not yet implemented in clTecplot::read_datasection");
		int num_elem=currentzone->get_num_elem();
		int num_conn=currentzone->get_num_conn();
		if (num_conn<=0) fatal_err("Illegal value of num_conn in clTecplot::read_datasection");
		clVector_namespace::clMatrix<int> *conn=currentzone->get_connptr();
		conn->allocate(num_elem,num_conn);
		connptr=conn->get_startptr();
		num_data=num_elem*num_conn;
		numbuffer=(int)((num_data-0.5)/IOBUFFERSIZE)+1;
		buffersize=IOBUFFERSIZE;
		for (nbuffer=0; nbuffer<numbuffer; nbuffer++) {
			if (nbuffer==numbuffer-1) buffersize=num_data-nbuffer*IOBUFFERSIZE;
			if (!read_binentry_int(fileptr,readintbuffer.get_startptr(),buffersize)) return 0;
			readintptr=readintbuffer.get_startptr();
			for (n=0; n<buffersize; n++) (*connptr++)=(int)(*readintptr++);
		}
	}

	return 1;
};

/* Write a binary Tecplot file */
void clTecplot::write(const char* const filename)
{
	clVector_namespace::clVector<int> writevar(_num_var,1);
	int writezone=-1;
	write(filename,writezone,writevar);
};
void clTecplot::write(const char* const filename,const int writezone)
{
	clVector_namespace::clVector<int> writevar(_num_var,1);
	write(filename,writezone,writevar);
};
void clTecplot::write(const char* const filename,const int n_varlow,const int n_varupp)
{
	clVector_namespace::clVector<int> writevar(_num_var,0);
	for (int n_var=0; n_var<_num_var; n_var++) 
		if (n_var>=n_varlow && n_var<=n_varupp) writevar(n_var)=1;
	int writezone=-1;
	write(filename,writezone,writevar);
};
void clTecplot::write(const char* const filename,const clVector_namespace::clVector<int> &writevar)
{
	int writezone=-1;
	write(filename,writezone,writevar);
};
void clTecplot::write(const char* const filename,const int writezone,const clVector_namespace::clVector<int> &writevar)
{
	using clUtils_namespace::write;

	clUtils_namespace::my_int writeint;

	FILE *fileptr=fopen(filename,"wb");
	if (!fileptr) fatal_err(1,filename);

	int saveversion=_version;
	_version=112;

	/*  Magic number  */
	unsigned int size;
	const char* const magic="#!TDV";
	void *io_ptr=(void *)magic;
	size=5*sizeof(char);
	if (fwrite(io_ptr,1,size,fileptr)!=size) fatal_err(3,filename);
	char *version=new char [4];
	sprintf(version,"%d",_version);
	version[3]='\0';
	io_ptr=(void *)version;
	size=3*sizeof(char);
	if (fwrite(io_ptr,1,size,fileptr)!=size) fatal_err(3,filename);
	writeint=1;
	if (!write(fileptr,writeint)) fatal_err(3,filename);
	if (!write_title(fileptr,writevar)) fatal_err(3,filename);
	delete [] version;

	int n_zone;
	clUtils_namespace::my_float writefloat;
	clTecplot_zone *currentzone=NULL;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		if (writezone>-1 && n_zone!=writezone) continue;
		currentzone=_zones[n_zone];
		if (!currentzone) 
			fatal_err("currentzone not allocated in clTecplot::write");
		writefloat=299.0;
		if (!write(fileptr,writefloat)) fatal_err(3,filename);
		if (!write_zone(fileptr,currentzone)) fatal_err(3,filename);
//		write_geometries(fileptr);
//		write_text(fileptr);
//		write_customlabel(fileptr);
//		write_userrec(fileptr);
//		write_dataauxdata(fileptr);
//		write_varauxdata(fileptr);
	}
	writefloat=357.0;
	if (!write(fileptr,writefloat)) fatal_err(3,filename);

	currentzone=_firstzone;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
fprintf(stdout,"n_zone=%d writezone=%d\n",n_zone,writezone);
		if (writezone>-1 && n_zone!=writezone) continue;
		currentzone=_zones[n_zone];
		if (!currentzone) 
			fatal_err("currentzone not allocated in clTecplot::write");
		writefloat=299.0;
		if (!write(fileptr,writefloat)) fatal_err(3,filename);
fprintf(stdout,"currentzone=%p\n",currentzone);
writevar.write_ascii(stdout);
		if (!write_datasection(fileptr,currentzone,writevar)) fatal_err(3,filename);
	}

	fclose(fileptr);

	_version=saveversion;

};

/* Write the title to a binary Tecplot file */
int clTecplot::write_title(FILE *fileptr,const clVector_namespace::clVector<int> &writevar) const
{
	using clUtils_namespace::write;

	clUtils_namespace::my_int writeint;

	/*  Filetype  */
	if (_version==112) {
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
	}

	/*  Title  */
	int i,namelength;
	const char *title=get_title();
	namelength=strlen(title);
	if (namelength) {
		for (i=0; i<namelength; i++) {
			writeint=(int)(title[i]);
			if (!write(fileptr,writeint)) return 0;
		}
	}
	writeint=0;
	if (!write(fileptr,writeint)) return 0;

	/*  Number of variables  */
	int n_var,num_var=0;
	for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) num_var++;
	writeint=num_var;
	if (!write(fileptr,writeint)) return 0;

	/*  Number of variables and variable names  */
	const char *varname=NULL;
	for (n_var=0; n_var<_num_var; n_var++) {
		if (!writevar(n_var)) continue;
		varname=get_varname(n_var);
		if (varname) {
			namelength=strlen(varname);
			for (i=0; i<namelength; i++) {
				writeint=(int)(varname[i]);
				if (!write(fileptr,writeint)) return 0;
			}
		}
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
	}

	return 1;
};

/* Write a zone to a binary Tecplot file */
int clTecplot::write_zone(FILE *fileptr,const clTecplot_zone* const currentzone) const
{
	using clUtils_namespace::write;

	clUtils_namespace::my_int writeint;
	clUtils_namespace::my_double writedouble;

	/*  Zone title  */
	int i,namelength;
	const char *title=currentzone->get_title();
	namelength=strlen(title);
	if (namelength) {
		for (i=0; i<namelength; i++) {
			writeint=(int)(title[i]);
			if (!write(fileptr,writeint)) return 0;
		}
	}
	writeint=0;
	if (!write(fileptr,writeint)) return 0;

	if (_version==102) {
		/*  Not used, set to -1  */
		writeint=-1;
		if (!write(fileptr,writeint)) return 0;
		/*  ZoneType  */
		writeint=currentzone->_zonetype;
		if (!write(fileptr,writeint)) return 0;
		/*  Solution time  */
		writedouble=0.0;
		if (!write(fileptr,writedouble)) return 0;
		/*  StrandID  */
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
	} else if (_version==112) {
		/*  ParentZone  */
		writeint=-1;
		if (!write(fileptr,writeint)) return 0;
		/*  StrandID  */
		writeint=-1;
		if (!write(fileptr,writeint)) return 0;
		/*  Solution time  */
		writedouble=0.0;
		if (!write(fileptr,writedouble)) return 0;
		/*  Not used, set to -1 */
		writeint=-1;
		if (!write(fileptr,writeint)) return 0;
		/*  ZoneType  */
		writeint=currentzone->_zonetype;
		if (!write(fileptr,writeint)) return 0;
	}

	/*  Var location  */
	if (_version==112) {
		writeint=currentzone->get_num_var_cc();
		if (writeint) writeint=1;
		if (!write(fileptr,writeint)) return 0;
		if (writeint) {
			if (currentzone->get_varlocationdim()!=_num_var) fatal_err("Inconsistent varlocationdim in clTecplot::write_zone");
			for (int n_var=0; n_var<_num_var; n_var++) {
				writeint=currentzone->get_varlocation(n_var);
				if (!write(fileptr,writeint)) return 0;
			}
		}
		/*  Are 1-to-1 face neighbors supplied  */
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
		/*  Number of miscellaneous user defined neighbor connections  */
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
	}

	/*  Zone dimensions  */
	if (currentzone->zonetype_ordered()) {
		/*  I,J,K dimensions  */
		writeint=currentzone->get_idim();
		if (!write(fileptr,writeint)) return 0;
		writeint=currentzone->get_jdim();
		if (!write(fileptr,writeint)) return 0;
		writeint=currentzone->get_kdim();
		if (!write(fileptr,writeint)) return 0;
	} else {
		/*  Number of points  */
		writeint=currentzone->get_num_node();
		if (!write(fileptr,writeint)) return 0;
		/*  Faces info  */
		if (currentzone->_zonetype==FEPOLYGON || currentzone->_zonetype==FEPOLYHEDRON) 
			fatal_err("FEPOLYGON or FEPOLYHEDRON not yet implemented");
		/*  Number of elements  */
		writeint=currentzone->get_num_elem();
		if (!write(fileptr,writeint)) return 0;
		/* For future use */
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
	}
	/*  Auxiliary data name/values  */
	writeint=0;
	if (!write(fileptr,writeint)) return 0;

	return 1;
};

/* Write geometries to a binary Tecplot file */
void clTecplot::write_geometries(FILE *fileptr) const
{
	fatal_err("write_geometries not yet implemented");
};

/* Write text to a binary Tecplot file */
void clTecplot::write_text(FILE *fileptr) const
{
	fatal_err("write_text not yet implemented");
};

/* Write a custom label to a binary Tecplot file */
void clTecplot::write_customlabel(FILE *fileptr) const
{
	fatal_err("write_customlabel not yet implemented");
};

/* Write a user record to a binary Tecplot file */
void clTecplot::write_userrec(FILE *fileptr) const
{
	fatal_err("write_userrec not yet implemented");
};

/* Write auxiliary data to a binary Tecplot file */
void clTecplot::write_dataauxdata(FILE *fileptr) const
{
	fatal_err("write_dataauxdata not yet implemented");
};

/* Write auxiliary variable data to a binary Tecplot file */
void clTecplot::write_varauxdata(FILE *fileptr) const
{
	fatal_err("read_varauxdata not yet implemented");
};

/* Write the data section to a binary Tecplot file */
int clTecplot::write_datasection(FILE *fileptr,clTecplot_zone* const currentzone,
	const clVector_namespace::clVector<int> &writevar) const
{
	using clUtils_namespace::write;

	clUtils_namespace::my_int writeint;
	clUtils_namespace::my_int *writeintptr=NULL;
	clVector_namespace::clVector<int> writeintbuffer(IOBUFFERSIZE);
	clUtils_namespace::my_double writedouble;
	clUtils_namespace::my_double *writedoubleptr=NULL;
	clVector_namespace::clVector<clUtils_namespace::my_double> writedoublebuffer(IOBUFFERSIZE);

	int n,n_var,num_var,num_varsharing,n_data,num_data,nbuffer,numbuffer,buffersize;
	const int *connptr=NULL;
	const double *dataptr=NULL;

	num_var=0;
	for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) num_var++;
	num_varsharing=0;
	if (currentzone->_num_varsharing) for (n_var=0; n_var<_num_var; n_var++) 
		if (writevar(n_var) && currentzone->_varsharingzone(n_var)!=-1) num_varsharing++;
	num_var-=num_varsharing;

	/*  Data type  */
	for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
		writeint=2;
		if (!write(fileptr,writeint)) return 0;
	}
	/*  Passive variables  */
	if (_version==112) {
		writeint=0;
		if (!write(fileptr,writeint)) return 0;
	}
	/*  Variable sharing  */
	writeint=0;
	if (num_varsharing) writeint=1;
	if (!write(fileptr,writeint)) return 0;
	if (num_varsharing) {
		for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
			writeint=currentzone->_varsharingzone(n_var);
			if (!write(fileptr,writeint)) return 0;
		}
	}
	writeint=-1;
	if (!write(fileptr,writeint)) return 0;
	if (_version==112) {
		if (currentzone->_varmin.get_dim() || !currentzone->_varmax.get_dim())
			currentzone->calculate_varrange();
		for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
			if (!num_varsharing || currentzone->get_varsharingzone(n_var)==-1) {
				writedouble=currentzone->_varmin(n_var);
				if (!write(fileptr,writedouble)) return 0;
				writedouble=currentzone->_varmax[n_var];
				if (!write(fileptr,writedouble)) return 0;
			}
		}
	}
	if (num_var) {
		if (!currentzone->_zonetype) {
			const int idim=currentzone->get_idim();
			const int jdim=currentzone->get_jdim();
			const int kdim=currentzone->get_kdim();
			num_data=idim*jdim*kdim;
			numbuffer=(int)((num_data-0.5)/IOBUFFERSIZE)+1;
		}
		for (n_var=0; n_var<_num_var; n_var++) if (writevar(n_var)) {
			if (!num_varsharing || currentzone->get_varsharingzone(n_var)==-1) {
				if (!currentzone->zonetype_ordered()) {
					if (currentzone->get_varlocation(n_var)) num_data=currentzone->get_num_elem();
					else num_data=currentzone->get_num_node();
					numbuffer=(int)((num_data-0.5)/IOBUFFERSIZE)+1;
				}
				dataptr=currentzone->get_datastartptr(n_var);
				buffersize=IOBUFFERSIZE;
				for (nbuffer=0; nbuffer<numbuffer; nbuffer++) {
					if (nbuffer==numbuffer-1) buffersize=num_data-nbuffer*IOBUFFERSIZE;
					writedoubleptr=writedoublebuffer.get_startptr();
					for (n=0; n<buffersize; n++) (*writedoubleptr++)=(*dataptr++);
					if (!write(fileptr,writedoublebuffer.get_startptr(),buffersize)) return 0;
				}
			}
		}
	}
	if (currentzone->_zonetype) {
		if (currentzone->_zonetype==FEPOLYGON || currentzone->_zonetype==FEPOLYHEDRON) 
			fatal_err("FEPOLYGON and FEPOLYHEDRON not yet implemented in clTecplot::read_datasection");
		int num_elem=currentzone->get_num_elem();
		int num_conn=currentzone->get_num_conn();
		connptr=currentzone->get_connstartptr();
		num_data=num_elem*num_conn;
		numbuffer=(int)((num_data-0.5)/IOBUFFERSIZE)+1;
		buffersize=IOBUFFERSIZE;
		for (nbuffer=0; nbuffer<numbuffer; nbuffer++) {
			if (nbuffer==numbuffer-1) buffersize=num_data-nbuffer*IOBUFFERSIZE;
			writeintptr=writeintbuffer.get_startptr();
			for (n=0; n<buffersize; n++) (*writeintptr++)=(*connptr++);
			if (!write(fileptr,writeintbuffer.get_startptr(),buffersize)) return 0;
		}
	}

	return 1;
};

/* Convert triangular connectivity (set last repeated index to -1) */
void clTecplot::convert_fetriangle_conn(void)
{
	clTecplot_zone *currentzone=NULL;

	/* Check for variable sharing */
	int n_zone;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->_num_varsharing!=0 && currentzone->_num_varsharing!=_num_var) 
			fatal_err("num_varsharing!=0 and num_varsharing!=_num_var not implemented in clTecplot::convert_fetriangle_conn");
	}

	int num_elem,n_elem,num_conn,conn3,conn4;
	int *connptr=NULL;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->_zonetype!=FEQUADRILATERAL) continue;
		num_conn=currentzone->get_num_conn();
		if (num_conn!=4) fatal_err("num_conn!=4 in clTecplot::convert_fetriangle_conn");
		num_elem=currentzone->get_num_elem();
		connptr=currentzone->get_connstartptr();
		for (n_elem=0; n_elem<num_elem; n_elem++) {
			connptr+=2;
			conn3=(*connptr++);
			conn4=(*connptr);
			if (conn3==conn4) (*connptr)=-1;
			connptr++;
		}
	}
};

/* Convert FETRIANGLE zones to FEQUADRILATERAL zones */
void clTecplot::convert_fetriangle_to_fequadrilateral(void)
{
	clTecplot_zone *currentzone=NULL;

	/* Check for variable sharing */
	int n_zone;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->_num_varsharing!=0 && currentzone->_num_varsharing!=_num_var) 
			fatal_err("num_varsharing!=0 and num_varsharing!=_num_var not implemented in clTecplot::convert_fetriangle_to_fequadrilateral");
	}

	int num_conn,num_elem,n_conn,n_elem;
	int *connptr=NULL;
	const int *connorgptr=NULL;
	clVector_namespace::clMatrix<int> conn;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->_zonetype==FEQUADRILATERAL) continue;

		num_conn=currentzone->get_num_conn();
		if (num_conn!=3) fatal_err("num_conn!=3 in clTecplot::convert_fetriangle_to_fequadrilateral");

		num_elem=currentzone->get_num_elem();
		connorgptr=currentzone->get_connstartptr();
		conn.allocate(num_elem,4,-1);
		connptr=conn.get_startptr();
		for (n_elem=0; n_elem<num_elem; n_elem++) {
			for (n_conn=0; n_conn<num_conn; n_conn++) (*connptr++)=(*connorgptr++);
			connptr++;
		}
		currentzone->set_conn(conn);
		currentzone->set_zonetype_num_conn("FEQUADRILATERAL");
		connptr=currentzone->get_connstartptr();
		if (currentzone->get_num_conn()!=4) fatal_err("num_conn not set correctly in convert_fetriangle_to_fequadrilateral");
	}
};

/* Convert FEQUADRILATERAL zones to FETRIANGLE zones */
void clTecplot::convert_fequadrilateral_to_fetriangle(void)
{
	clTecplot_zone *currentzone=NULL;

	/* Check for variable sharing and consistent connectivity types */
	int n_zone,numquad,varsharing;
	for (n_zone=0,varsharing=0,numquad=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->_num_varsharing!=0 && currentzone->_num_varsharing!=_num_var) 
			fatal_err("num_varsharing!=0 and num_varsharing!=_num_var not implemented in clTecplot::convert_fequadrilateral_to_fetriangle");
		if (currentzone->_num_varsharing) varsharing++;
		if (currentzone->_zonetype==FEQUADRILATERAL) numquad++;
	}
	if (varsharing && numquad<_num_zone)
		fatal_err("variable sharing not implemented for numquad<_num_zone in clTecplot::convert_fequadrilateral_to_fetriangle");

	int iconn,jconn,n_conn,num_conn,n_elem,num_elem,found_match;
	int *connptr=NULL;
	const int *connorgptr=NULL;
	clVector_namespace::clMatrix<int> conn;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->_zonetype!=FEQUADRILATERAL) continue;
		num_conn=currentzone->get_num_conn();
		if (num_conn!=4) fatal_err("num_conn!=4 in clTecplot::convert_fequadrilateral_to_fetriangle");
		found_match=0;
		connorgptr=currentzone->get_connstartptr();
		for (iconn=0; iconn<num_conn; iconn++) {
			for (jconn=iconn+1; jconn<num_conn; jconn++) {
				if (connorgptr[iconn]==connorgptr[jconn]) {
					found_match++;
					break;
				}
			}
			if (found_match) break;
		}
		if (!found_match) continue;
		num_elem=currentzone->get_num_elem();
		conn.allocate(num_elem,3);
		connptr=conn.get_startptr();
		for (n_elem=0; n_elem<num_elem; n_elem++) {
			if (connorgptr[iconn]!=connorgptr[jconn]) 
				fatal_err("inconsistent connectivities in clTecplot::convert_fequadrilateral_to_fetriangle");
			for (n_conn=0; n_conn<num_conn; n_conn++) {
				if (n_conn==iconn) connorgptr++;
				else (*connptr++)=(*connorgptr++);
			}
		}
		currentzone->set_conn(conn);
		currentzone->set_zonetype_num_conn("FETRIANGLE");
	}
};

/* Create local to global and global to local lookup tables from connectivities */
int clTecplot::get_ltog_gtol_tables(const int num_node,const int num_elem,const int num_conn,
	const int *connptr,clVector_namespace::clVector<int> &ltog,clVector_namespace::clVector<int> &gtol) const
{
	int n_data,num_data,n_node,n_node_used,num_node_used;
	int *indexptr=NULL,*gtolptr=NULL,*ltogptr=NULL;

	clVector_namespace::clVector<int> index(num_node,0);
	num_data=num_elem*num_conn;
	for (n_data=0; n_data<num_data; n_data++) index(*connptr++)++;
	gtol.allocate(num_node);
	gtolptr=gtol.get_startptr();
	indexptr=index.get_startptr();
	for (n_node=0,num_node_used=0; n_node<num_node; n_node++) {
		if (*indexptr++) num_node_used++;
		(*gtolptr++)=-1;
	}
	ltog.allocate(num_node_used);
	ltogptr=ltog.get_startptr();
	indexptr=index.get_startptr();
	gtolptr=gtol.get_startptr();
	for (n_node=0,n_node_used=0; n_node<num_node; n_node++,gtolptr++) {
		if (*indexptr++) {
			(*ltogptr++)=n_node;
			*gtolptr=n_node_used++;
		}
	}
	if (n_node_used!=num_node_used) fatal_err("n_node_used!=num_node_used in clTecplot::get_ltog_gtol_tables");

	return num_node_used;
}

/* Eliminate varsharing and unreferenced nodes */
void clTecplot::eliminate_varsharing(void)
{
	int n_data,num_data,num_conn,num_elem,n_node,num_node,n_var,n_zone,nn_zone;
	int *gtolptr=NULL,*ltogptr=NULL,*connptr=NULL;
	double *data_datazone=NULL,*dataptr=NULL;
	const double *dataorgptr=NULL;
	clVector_namespace::clVector<int> varsharingzone;
	clVector_namespace::clVector<int> num_node_used;
	clVector_namespace::clVector<int> *gtol=NULL;
	clVector_namespace::clVector<int> *ltog=NULL;
	clVector_namespace::clVector<double> data;
	clTecplot_zone *currentzone=NULL,*datazone=NULL;

	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (currentzone->zonetype_ordered()) fatal_err("clTecplot::eliminate_unused_nodes not implemented for ordered zones");
	}

	gtol=new clVector_namespace::clVector<int> [_num_zone];
	if (!gtol) fatal_err(0,"gtol in clTecplot::eliminate_varsharing");
	ltog=new clVector_namespace::clVector<int> [_num_zone];
	if (!ltog) fatal_err(0,"ltog in clTecplot::eliminate_varsharing");

	/* Find which zones have shared variables */
	if (!_zonesharing.get_dim()) set_zonesharing();
	if (_num_zonesharing) {
		varsharingzone.allocate(_num_zone,-1);
		for (n_zone=0; n_zone<_num_zone; n_zone++) {
			currentzone=_zones[n_zone];
			if (currentzone->_num_varsharing) {
				nn_zone=-1;
				for (n_var=0; n_var<_num_var; n_var++) if (currentzone->_varsharingzone(n_var)!=-1) {
					if (nn_zone<0) {
						nn_zone=currentzone->_varsharingzone(n_var);
						varsharingzone(n_zone)=nn_zone;
					} else {
						if (nn_zone!=currentzone->_varsharingzone(n_var)) 
							fatal_err("Variable sharing with different zones not implemented in clTecplot::eliminate_varsharing");
					}
				}
				if (nn_zone<0) fatal_err("Unable to find sharing zone in clTecplot::eliminate_varsharing");
				if (!_zonesharing(nn_zone)) fatal_err("Inconsistent zone sharing info in clTecplot::eliminate_varsharing");
			}
		}
	}

	/* Get ltog and gtol tables for all zones */
	num_node_used.allocate(_num_zone);
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (_num_zonesharing) if (_zonesharing(n_zone) && currentzone->_num_varsharing)
			fatal_err("Zones cannot share and be shared at the same time in clTecplot::eliminate_varsharing");
		num_node=currentzone->get_num_node();
		num_elem=currentzone->get_num_elem();
		num_conn=currentzone->get_num_conn();
		num_node_used(n_zone)=get_ltog_gtol_tables(num_node,num_elem,num_conn,currentzone->get_connstartptr(),ltog[n_zone],gtol[n_zone]);
	}

	/* Eliminate shared variables */
	if (_num_zonesharing) for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		if (!currentzone->_num_varsharing) continue;
		datazone=_zones[varsharingzone(n_zone)];
		if (!datazone) fatal_err("Unable to get datazone in clTecplot::eliminate_varsharing");
		num_node=num_node_used(n_zone);
		num_elem=currentzone->get_num_elem();
		currentzone->set_num_node(num_node);
		for (n_var=0; n_var<_num_var; n_var++) {
			if (!currentzone->get_varlocation(n_var)) {
				data_datazone=datazone->get_datastartptr(n_var);
				if (!data_datazone) fatal_err("Unable to get data_datazone in clTecplot::eliminate_varsharing");
				data.allocate(num_node);
				dataptr=data.get_startptr();
				ltogptr=ltog[n_zone].get_startptr();
				for (n_node=0; n_node<num_node; n_node++) {
					dataorgptr=data_datazone+(*ltogptr++);
					(*dataptr++)=(*dataorgptr);
				}
				currentzone->set_data(data,n_var);
			}
		}
		num_conn=currentzone->get_num_conn();
		num_data=num_conn*num_elem;
		connptr=currentzone->get_connstartptr();
		gtolptr=gtol[n_zone].get_startptr();;
		for (n_data=0; n_data<num_data; n_data++,connptr++) *connptr=gtolptr[*connptr];
		currentzone->reset_varsharingzone();
	}

	/* Eliminate unused nodes */
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		if (_num_zonesharing) if (!_zonesharing(n_zone)) continue;
		currentzone=_zones[n_zone];
		if (currentzone->_num_varsharing) 
			fatal_err("Shared variables where there should be none in clTecplot::eliminate_varsharing");
		num_node=num_node_used(n_zone);
		num_elem=currentzone->get_num_elem();
		currentzone->set_num_node(num_node);
		for (n_var=0; n_var<_num_var; n_var++) {
			data_datazone=currentzone->get_datastartptr(n_var);
			if (!currentzone->get_varlocation(n_var)) {
				data_datazone=currentzone->get_datastartptr(n_var);
				if (!data_datazone) fatal_err("Unable to get data_datazone in clTecplot::eliminate_varsharing");
				data.allocate(num_node);
				dataptr=data.get_startptr();
				ltogptr=ltog[n_zone].get_startptr();;
				for (n_node=0; n_node<num_node; n_node++) {
					dataorgptr=data_datazone+(*ltogptr++);
					(*dataptr++)=(*dataorgptr);
				}
				currentzone->set_data(data,n_var);
			}
		}
		num_conn=currentzone->get_num_conn();
		num_data=num_conn*num_elem;
		connptr=currentzone->get_connstartptr();
		gtolptr=gtol[n_zone].get_startptr();
		for (n_data=0; n_data<num_data; n_data++,connptr++) *connptr=gtolptr[*connptr];
		currentzone->reset_varsharingzone();
	}

	delete [] gtol;
	delete [] ltog;

	set_zonesharing();
	if (_num_zonesharing || _zonesharing.get_dim()) 
		fatal_err("Failed to eliminate all shared variables in clTecplot::eliminate_varsharing");
};

/* Copy a zone from source to current */
void clTecplot::copy(const clTecplot* const source)
{
	int num_zone=source->get_num_zone();
	const clTecplot_zone *sourcezone=NULL;
	copy_header(source);
	for (int n_zone=0; n_zone<num_zone; n_zone++) {
		sourcezone=source->get_zoneptr(n_zone);
		copy_zone(sourcezone);
	}
};

/* Add data from source to current for each zone */
void clTecplot::add_data(const clTecplot* const source)
{
	const clVector_namespace::clVector<int> varindex(_num_var,1);
	add_data(source,varindex);
};
/* Add data from source to current for each zone for each variable that has varindex set */
void clTecplot::add_data(const clTecplot* const source,const clVector_namespace::clVector<int> &varindex)
{
	int n_zone;
	const clTecplot_zone *sourcezone=NULL;

	if (_num_zone!=source->get_num_zone()) fatal_err("Incompatible num_zone in clTecplot::add_data");
	if (_num_var!=source->get_num_var()) fatal_err("Incompatible num_var in clTecplot::add_data");
	if (_num_zonesharing!=source->get_num_zonesharing()) fatal_err("Incompatible num_zonesharing in clTecplot::add_data");
	if (_num_zonesharing) for (n_zone=0; n_zone<_num_zone; n_zone++) 
		if (_zonesharing(n_zone)!=source->get_zonesharing(n_zone)) fatal_err("Incompatible zonesharing in clTecplot_zone::add_data");
//RV TODO CHECK varlocation here as well
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		sourcezone=source->get_zoneptr(n_zone);
		_zones[n_zone]->add_data(sourcezone,varindex);
	}
};

/* Add delta to the data in each zone */
void clTecplot::translate_data(const double delta)
{
	clVector_namespace::clVector<int> varindex(_num_var,1);
	translate_data(delta,varindex);
};
/* Add delta to variable n_var in each zone */
void clTecplot::translate_data(const double delta,const int n_var)
{
	clVector_namespace::clVector<int> varindex(_num_var,0);
	varindex(n_var)=1;
	translate_data(delta,varindex);
};
/* Add delta to the data in each zone for each variable that has varindex set */
void clTecplot::translate_data(const double delta,const clVector_namespace::clVector<int> &varindex)
{
	for (int n_zone=0; n_zone<_num_zone; n_zone++)
		_zones[n_zone]->translate_data(delta,varindex);
};

/* Multiply the data in each zone by factor */
void clTecplot::multiply_data(const double factor)
{
	clVector_namespace::clVector<int> varindex(_num_var,1);
	multiply_data(factor,varindex);
};
/* Multiply variable n_var in each zone by factor */
void clTecplot::multiply_data(const double factor,const int n_var)
{
	clVector_namespace::clVector<int> varindex(_num_var,0);
	varindex(n_var)=1;
	multiply_data(factor,varindex);
};
/* Multiply the data in each zone by factor for each variable that has varindex set */
void clTecplot::multiply_data(const double factor,const clVector_namespace::clVector<int> &varindex)
{
	for (int n_zone=0; n_zone<_num_zone; n_zone++)
		_zones[n_zone]->multiply_data(factor,varindex);
};

/* Mirror a variable */
void clTecplot::mirror(const int n_var)
{
	int n_data,num_data,num_var,n_zone,mirrorflag;
	double *dataptr=NULL;
	clTecplot_zone *currentzone=NULL;

	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		currentzone=_zones[n_zone];
		num_var=_num_var-currentzone->_num_varsharing;
		if (currentzone->zonetype_ordered()) {
			const int idim=currentzone->get_idim();
			const int jdim=currentzone->get_jdim();
			const int kdim=currentzone->get_kdim();
			num_data=idim*jdim*kdim;
		} else {
			if (currentzone->get_varlocation(n_var)) num_data=currentzone->get_num_elem();
			else num_data=currentzone->get_num_node();
		}
		if (num_var) {
			mirrorflag=0;
			if (!currentzone->_num_varsharing) mirrorflag=1;
			else if (currentzone->_varsharingzone(n_var)==-1) mirrorflag=1;
			if (mirrorflag) {
				dataptr=currentzone->get_datastartptr(n_var);
				for (n_data=0; n_data<num_data; n_data++,dataptr++) (*dataptr)=-(*dataptr);
			}
		}
	}
};

/* Find the smallest and/or largest vertex in each zone */
double clTecplot::find_min_vertex(const double mincutoff) const
{
	double minvertex=1.e30,maxvertex=0.0;
	find_minmax_vertex(minvertex,maxvertex,mincutoff);
	return minvertex;
};
double clTecplot::find_max_vertex(const double mincutoff) const
{
	double minvertex=1.e30,maxvertex=0.0;
	find_minmax_vertex(minvertex,maxvertex,mincutoff);
	return maxvertex;
};
void clTecplot::find_minmax_vertex(double &minvertex,double &maxvertex,const double mincutoff) const
{
	if (minvertex<0.0) minvertex=1.e30;
	if (maxvertex<0.0) maxvertex=0.0;
	for (int n_zone=0; n_zone<_num_zone; n_zone++) {
		if (!_zones[n_zone]) 
			fatal_err("_zones[n_zone] not allocated in clTecplot::find_minmax_vertex");
		if (_zones[n_zone]->_num_varsharing)
			fatal_err("Variables in _zones[n_zone] should be unshared in clTecplot::find_minmax_vertex");
		_zones[n_zone]->find_minmax_vertex(minvertex,maxvertex);
	}
};

void clTecplot_zone_structured::add_variables(const int num_varadd)
{
	if (num_varadd<0) fatal_err("num_varadd<0 in clTecplot_zone_structured::add_variables");
	if (get_num_varsharing()) fatal_err("Shared variables not yet implemented in clTecplot_zone_structured::add_variables");
	if (!_idim || !_jdim|| !_kdim) fatal_err("Illegal zone dimensions in clTecplot::add_variables");

	const int num_var=get_num_var();
	const int num_varnew=num_var+num_varadd;
	clVector_namespace::clMatrix3D<double> *datanew=new clVector_namespace::clMatrix3D<double> [num_varnew];
	int n_var=0;
	if (_data) for (; n_var<num_var; n_var++) datanew[n_var].move(_data[n_var]);
	for (; n_var<num_varnew; n_var++) datanew[n_var].allocate(_kdim,_jdim,_idim,0);
	deallocate_data();
	_data=datanew;
	set_num_var(num_varnew);
};

void clTecplot_zone_unstructured::add_variables(const int num_varadd,const int varlocation)
{
	if (num_varadd<0) fatal_err("num_varadd<0 in clTecplot::add_variables");

	int n_var;
	const int num_var=get_num_var();
	const int num_varnew=num_var+num_varadd;
	clVector_namespace::clVector<int> varlocationnew(num_varnew);
	for (n_var=0; n_var<num_var; n_var++) varlocationnew(n_var)=get_varlocation(n_var);
	clVector_namespace::clVector<double> *datanew=new clVector_namespace::clVector<double> [num_varnew];
	n_var=0;
	if (_data) for (; n_var<num_var; n_var++) datanew[n_var].move(_data[n_var]);
	for (; n_var<num_varnew; n_var++) {
		varlocationnew(n_var)=varlocation;
		if (varlocation) datanew[n_var].allocate(_num_elem,0);
		else datanew[n_var].allocate(_num_node,0);
	}

	clVector_namespace::clVector<int> varsharingzonenew;
	if (get_num_varsharing()) {
		int varsharingzone=-1;
		varsharingzonenew.allocate(num_varnew,-1);
		for (n_var=0; n_var<num_var; n_var++) {
			varsharingzonenew(n_var)=get_varsharingzone(n_var);
			if (varsharingzonenew(n_var)>=0) {
				if (varsharingzone<0) varsharingzone=varsharingzonenew(n_var);
				else if (varsharingzone!=varsharingzonenew(n_var)) fatal_err("varsharingzone!=varsharingzonenew(n_var) in clTecplot::add_variables");
			}
		}
		if (varsharingzone<-1) fatal_err("varsharingzone<-1 in clTecplot::add_variables");
		if (!varlocation) for (; n_var<num_varnew; n_var++) varsharingzonenew(n_var)=varsharingzone;
	}

	set_num_var(num_varnew,varlocationnew,varsharingzonenew);

	deallocate_data();
	_data=datanew;
};

void clTecplot::add_variable(void)
{
	add_variables(1);
};
void clTecplot::add_variables(const int num_varadd,const int varlocation)
{
	if (num_varadd<0) fatal_err("num_varadd<0 in clTecplot::add_variables");

	clTecplot_zone *zoneptr=NULL;
	for (int n_zone=0; n_zone<_num_zone; n_zone++) {
		zoneptr=_zones[n_zone];
		zoneptr->_varmin.deallocate();
		zoneptr->_varmax.deallocate();
		if (zoneptr->zonetype_ordered()) ((clTecplot_zone_structured *)(zoneptr))->add_variables(num_varadd);
		else ((clTecplot_zone_unstructured *)(zoneptr))->add_variables(num_varadd,varlocation);
	}
	int n_var;
	const int num_varnew=_num_var+num_varadd;
	char *varname=NULL,vvarname[NAME_LENGTH];
	clString_namespace::clString_vector varnamenew(num_varnew);
	for (n_var=0; n_var<_num_var; n_var++) {
		varname=get_varname(n_var);
		varnamenew(n_var)=varname;
	}
	for (n_var=_num_var; n_var<num_varnew; n_var++) {
		sprintf(vvarname,"VAR%d",n_var);
		varnamenew(n_var)=vvarname;
	}
	dealloc_varname();
	_num_var=num_varnew;
	for (n_var=0; n_var<_num_var; n_var++) set_varname(n_var,varnamenew(n_var).get_string());
};

void clTecplot::copy_variable(const clTecplot* const source,const int varadd)
{
	copy_variable(source,varadd,varadd);
}
void clTecplot::copy_variable(const clTecplot* const source,const int varadd_source,const int varadd_dest)
{
	int n_data,num_data,n_zone;
	double *dataptr=NULL;
	const double *sourcedataptr=NULL;
	clTecplot_zone *zoneptr=NULL;
	const clTecplot_zone *sourcezoneptr=NULL;

	if (varadd_dest<0 || varadd_dest>=_num_var) fatal_err("varadd_dest<0 || varadd_source>=_num_var in clTecplot::copy_variable");
	if (_num_zonesharing) fatal_err("Shared variables not yet implemented in clTecplot::copy_variable");
	const int num_var_source=source->get_num_var();
	if (varadd_source<0 || varadd_source>=num_var_source) fatal_err("varadd_dest<0 || varadd_dest>=num_var_source in clTecplot::copy_variable");
	if (source->get_num_zonesharing()) fatal_err("Shared variables not yet implemented in clTecplot::copy_variable");
	if (_num_zone!=source->get_num_zone()) fatal_err("inconstent number of zones in clTecplot::copy_variable");

	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		zoneptr=_zones[n_zone];
		if (zoneptr->_num_varsharing) fatal_err("Shared variables not yet implemented in clTecplot::copy_variable");
		sourcezoneptr=source->get_zoneptr(n_zone);
		if (sourcezoneptr->_num_varsharing) fatal_err("Shared variables not yet implemented in clTecplot::copy_variable");
		zoneptr->_varmin.deallocate();
		zoneptr->_varmax.deallocate();
		if (zoneptr->get_varlocation(varadd_dest)!=sourcezoneptr->get_varlocation(varadd_source)) fatal_err("Inconsistent varlocation in clTecplot::copy_variable");
		if (zoneptr->zonetype_ordered()) {
			int idim=zoneptr->get_idim();
			int jdim=zoneptr->get_jdim();
			int kdim=zoneptr->get_kdim();
			if (idim!=sourcezoneptr->get_idim()) fatal_err("Inconsistent idim  in clTecplot::copy_variable");
			if (jdim!=sourcezoneptr->get_jdim()) fatal_err("Inconsistent jdim  in clTecplot::copy_variable");
			if (kdim!=sourcezoneptr->get_kdim()) fatal_err("Inconsistent idim  in clTecplot::copy_variable");
			num_data=idim*jdim*kdim;
		} else {
			if (zoneptr->get_varlocation(varadd_dest)) {
				num_data=zoneptr->get_num_elem();
				if (num_data!=sourcezoneptr->get_num_elem()) fatal_err("Inconsistent number of elements in clTecplot::copy_variable");
			} else {
				num_data=zoneptr->get_num_node();
				if (num_data!=sourcezoneptr->get_num_node()) fatal_err("Inconsistent number of nodes in clTecplot::copy_variable");
			}
		}
		dataptr=zoneptr->get_datastartptr(varadd_dest);
		sourcedataptr=sourcezoneptr->get_datastartptr(varadd_source);
		for (n_data=0; n_data<num_data; n_data++) (*dataptr++)=(*sourcedataptr++);
	}

	if (source->get_varname(varadd_source)) set_varname(varadd_dest,source->get_varname(varadd_source));
};

/* Remove variables varindexstart to and including varindexstart */
void clTecplot_zone_structured::remove_variable(const int varindexstart,const int vvarindexend)
{
	const int varindexend=(vvarindexend<0 ? varindexstart : vvarindexend);

	const int num_var=get_num_var();
	if (get_num_varsharing()) fatal_err("Shared variables not yet implemented in clTecplot::remove_variable");
	if (varindexstart<0 || varindexstart>num_var-1) fatal_err("Illegal varindexstart in clTecplot::remove_variable");
	if (varindexend>num_var-1) fatal_err("Illegal varindexend in clTecplot::remove_variable");
	if (varindexstart>varindexend) fatal_err("varindexstart>varindexend in clTecplot::remove_variable");
	if (!_idim || !_jdim|| !_kdim) fatal_err("Illegal zone dimensions in clTecplot::remove_variable");

	int n_var,n_varnew=0;
	const int num_varnew=num_var-varindexend+varindexstart-1;
	clVector_namespace::clMatrix3D<double> *datanew=new clVector_namespace::clMatrix3D<double> [num_varnew];
	for (n_var=0; n_var<num_var; n_var++) if (n_var<varindexstart || n_var>varindexend) {
		datanew[n_varnew].move(_data[n_var]);
		n_varnew++;
	}
	deallocate_data();
	_data=datanew;
	set_num_var(num_varnew);
};

/* Remove variables varindexstart to and including varindexstart */
void clTecplot_zone_unstructured::remove_variable(const int varindexstart,const int vvarindexend)
{
	const int varindexend=(vvarindexend<0 ? varindexstart : vvarindexend);

	const int num_var=get_num_var();
	if (varindexstart<0 || varindexstart>num_var-1) fatal_err("Illegal varindexstart in clTecplot::remove_variable");
	if (varindexend>num_var-1) fatal_err("Illegal varindexend in clTecplot::remove_variable");
	if (varindexstart>varindexend) fatal_err("varindexstart>varindexend in clTecplot::remove_variable");

	int n_var,n_varnew;
	const int num_varnew=num_var-varindexend+varindexstart-1;
	clVector_namespace::clVector<int> varlocationnew(num_varnew);
	clVector_namespace::clVector<double> *datanew=new clVector_namespace::clVector<double> [num_varnew];
	for (n_var=0,n_varnew=0; n_var<num_var; n_var++) if (n_var<varindexstart || n_var>varindexend) {
		varlocationnew(n_varnew)=get_varlocation(n_var);
		datanew[n_varnew++].move(_data[n_var]);
	}

	clVector_namespace::clVector<int> varsharingzonenew;
	if (get_num_varsharing()) {
		varsharingzonenew.allocate(num_varnew);
		for (n_var=0,n_varnew=0; n_var<num_var; n_var++) if (n_var<varindexstart || n_var>varindexend) {
			varsharingzonenew(n_varnew++)=get_varsharingzone(n_var);
		}
	}

	set_num_var(num_varnew,varlocationnew,varsharingzonenew);

	deallocate_data();
	_data=datanew;
};

/* Remove variables varindexstart to and including varindexstart */
void clTecplot::remove_variable(const int varindexstart,const int vvarindexend)
{
	const int varindexend=(vvarindexend<0 ? varindexstart : vvarindexend);

	if (varindexstart<0 || varindexstart>_num_var-1) fatal_err("Illegal varindexstart in clTecplot::remove_variable");
	if (varindexend>_num_var-1) fatal_err("Illegal varindexend in clTecplot::remove_variable");
	if (varindexstart>varindexend) fatal_err("varindexstart>varindexend in clTecplot::remove_variable");

	clTecplot_zone *zoneptr=NULL;
	for (int n_zone=0; n_zone<_num_zone; n_zone++) {
		zoneptr=_zones[n_zone];
		zoneptr->_varmin.deallocate();
		zoneptr->_varmax.deallocate();
		if (zoneptr->zonetype_ordered()) ((clTecplot_zone_structured *)(zoneptr))->remove_variable(varindexstart,varindexend);
		else ((clTecplot_zone_unstructured *)(zoneptr))->remove_variable(varindexstart,varindexend);
	}

	int n_var,n_varnew=0;
	char *varname=NULL;
	const int num_varnew=_num_var-varindexend+varindexstart-1;
	clString_namespace::clString_vector varnamenew(num_varnew);
	for (n_var=0; n_var<_num_var; n_var++) if (n_var<varindexstart || n_var>varindexend) {
		varname=get_varname(n_var);
		varnamenew(n_varnew)=varname;
		n_varnew++;
	}
	if (n_varnew!=num_varnew) fatal_err("n_varnew!=num_varnew in clTecplot::remove_variable");
	dealloc_varname();
	_num_var=num_varnew;
	for (n_var=0; n_var<_num_var; n_var++) set_varname(n_var,varnamenew(n_var).get_string());
};

/* Replace the current data with the error between the current data and the data in source */
void clTecplot::abserror_variable(const clTecplot* const source,const char* const varname)
{
	int n_var;
	int index=-1,indexsource=-1;
	for (n_var=0; n_var<_num_var; n_var++) if (!strcmp(get_varname(n_var),varname)) index=n_var;
	for (n_var=0; n_var<source->get_num_var(); n_var++) if (!strcmp(source->get_varname(n_var),varname)) indexsource=n_var;
	if (index<0 || indexsource<0) fatal_err("Unable to get variable index in clTecplot::relerror_variable");
	abserror_variable(source,indexsource,index);
};
void clTecplot::abserror_variable(const clTecplot* const source,const int varindexsource,int varindex)
{
	if (_num_zonesharing || source->get_num_zonesharing()) fatal_err("Shared variables not yet implemented in clTecplot::abserror_variable");

	if (varindex<0) varindex=varindexsource;
	if (varindex<0 || varindex>=_num_var) fatal_err("Illegal varindex in clTecplot::abserror_variable");
	if (varindexsource<0 || varindexsource>=source->get_num_var()) fatal_err("Illegal varindexsource in clTecplot::abserror_variable");
	if (source->get_num_zone()!=_num_zone) fatal_err("Inconsistent number of zones in clTecplot::abserror_variable");
	if (strcmp(source->get_varname(varindexsource),get_varname(varindex))) fatal_err("Inconsistent variable name in clTecplot::abserror_variable");

	int n_zone,n_data,num_data,idim,jdim,kdim;
	double *dataptr=NULL;
	const double *sourcedataptr=NULL;
	clTecplot_zone *zoneptr=NULL;
	const clTecplot_zone *sourcezoneptr=NULL;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		zoneptr=_zones[n_zone];
		if (zoneptr->_num_varsharing) fatal_err("Shared variables not yet implemented in clTecplot::abserror_variable");
		sourcezoneptr=source->get_zoneptr(n_zone);
		if (sourcezoneptr->_num_varsharing) fatal_err("Shared variables not yet implemented in clTecplot::abserror_variable");
		zoneptr->_varmin.deallocate();
		zoneptr->_varmax.deallocate();
		if (zoneptr->get_varlocation(varindex)!=sourcezoneptr->get_varlocation(varindexsource)) fatal_err("Inconsistent varlocation in clTecplot::abserror_variable");
		if (zoneptr->zonetype_ordered()) {
			idim=zoneptr->get_idim();
			jdim=zoneptr->get_jdim();
			kdim=zoneptr->get_kdim();
			if (idim!=sourcezoneptr->get_idim()) fatal_err("Inconsistent idim  in clTecplot::abserror_variable");
			if (jdim!=sourcezoneptr->get_jdim()) fatal_err("Inconsistent jdim  in clTecplot::abserror_variable");
			if (kdim!=sourcezoneptr->get_kdim()) fatal_err("Inconsistent idim  in clTecplot::abserror_variable");
			num_data=idim*jdim*kdim;
		} else {
			if (zoneptr->get_varlocation(varindex)) {
				num_data=zoneptr->get_num_elem();
				if (num_data!=sourcezoneptr->get_num_elem()) fatal_err("Inconsistent number of elements in clTecplot::abserror_variable");
			} else {
				num_data=zoneptr->get_num_node();
				if (num_data!=sourcezoneptr->get_num_node()) fatal_err("Inconsistent number of nodes in clTecplot::abserror_variable");
			}
		}
		dataptr=zoneptr->get_datastartptr(varindex);
		sourcedataptr=sourcezoneptr->get_datastartptr(varindexsource);
		for (n_data=0; n_data<num_data; n_data++) (*dataptr++)-=(*sourcedataptr++);
	}
	clString_namespace::clString varname("abserror(");
	varname+=get_varname(varindex);
	varname+=")";
	set_varname(varindex,varname.get_string());
};

/* Replace the current data with the relative error between the current data and the data in source */
void clTecplot::relerror_variable(const clTecplot* const source,const char* const varname)
{
	int n_var;
	int index=-1,indexsource=-1;
	for (n_var=0; n_var<_num_var; n_var++) if (!strcmp(get_varname(n_var),varname)) index=n_var;
	for (n_var=0; n_var<source->get_num_var(); n_var++) if (!strcmp(source->get_varname(n_var),varname)) indexsource=n_var;
	if (index<0 || indexsource<0) fatal_err("Unable to get variable index in clTecplot::relerror_variable");
	relerror_variable(source,indexsource,index);
};
void clTecplot::relerror_variable(const clTecplot* const source,const int varindexsource,int varindex)
{
	if (_num_zonesharing || source->get_num_zonesharing()) fatal_err("Shared variables not yet implemented in clTecplot::relerror_variable");

	if (varindex<0) varindex=varindexsource;
	if (varindex<0 || varindex>=_num_var) fatal_err("Illegal varindex in clTecplot::relerror_variable");
	if (varindexsource<0 || varindexsource>=source->get_num_var()) fatal_err("Illegal varindexsource in clTecplot::relerror_variable");
	if (source->get_num_zone()!=_num_zone) fatal_err("Inconsistent number of zones in clTecplot::relerror_variable");
	if (strcmp(source->get_varname(varindexsource),get_varname(varindex))) fatal_err("Inconsistent variable name in clTecplot::relerror_variable");

	int n_zone,n_data,num_data,idim,jdim,kdim;
	double *dataptr=NULL;
	const double *sourcedataptr=NULL;
	clTecplot_zone *zoneptr=NULL;
	const clTecplot_zone *sourcezoneptr=NULL;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		zoneptr=_zones[n_zone];
		if (zoneptr->_num_varsharing) fatal_err("Shared variables not yet implemented in clTecplot::relerror_variable");
		sourcezoneptr=source->get_zoneptr(n_zone);
		if (sourcezoneptr->_num_varsharing) fatal_err("Shared variables not yet implemented in clTecplot::relerror_variable");
		zoneptr->_varmin.deallocate();
		zoneptr->_varmax.deallocate();
		if (zoneptr->get_varlocation(varindex)!=sourcezoneptr->get_varlocation(varindexsource)) fatal_err("Inconsistent varlocation in clTecplot::relerror_variable");
		if (zoneptr->zonetype_ordered()) {
			idim=zoneptr->get_idim();
			jdim=zoneptr->get_jdim();
			kdim=zoneptr->get_kdim();
			if (idim!=sourcezoneptr->get_idim()) fatal_err("Inconsistent idim  in clTecplot::relerror_variable");
			if (jdim!=sourcezoneptr->get_jdim()) fatal_err("Inconsistent jdim  in clTecplot::relerror_variable");
			if (kdim!=sourcezoneptr->get_kdim()) fatal_err("Inconsistent idim  in clTecplot::relerror_variable");
			num_data=idim*jdim*kdim;
		} else {
			if (zoneptr->get_varlocation(varindex)) {
				num_data=zoneptr->get_num_elem();
				if (num_data!=sourcezoneptr->get_num_elem()) fatal_err("Inconsistent number of elements in clTecplot::relerror_variable");
			} else {
				num_data=zoneptr->get_num_node();
				if (num_data!=sourcezoneptr->get_num_node()) fatal_err("Inconsistent number of nodes in clTecplot::relerror_variable");
			}
		}
		dataptr=zoneptr->get_datastartptr(varindex);
		sourcedataptr=sourcezoneptr->get_datastartptr(varindexsource);
		for (n_data=0; n_data<num_data; n_data++,dataptr++,sourcedataptr++) {
			if (*dataptr) (*dataptr)=(1.0-(*sourcedataptr)/(*dataptr));
			else (*dataptr)-=(*sourcedataptr);
		}
	}
	clString_namespace::clString varname("relerror(");
	varname+=get_varname(varindex);
	varname+=")";
	set_varname(varindex,varname.get_string());
};

/* Replace the current data with the rpderror between the current data and the data in source */
void clTecplot::rpderror_variable(const clTecplot* const source,const char* const varname)
{
	int n_var;
	int index=-1,indexsource=-1;
	for (n_var=0; n_var<_num_var; n_var++) if (!strcmp(get_varname(n_var),varname)) index=n_var;
	for (n_var=0; n_var<source->get_num_var(); n_var++) if (!strcmp(source->get_varname(n_var),varname)) indexsource=n_var;
	if (index<0 || indexsource<0) fatal_err("Unable to get variable index in clTecplot::rpderror_variable");
	rpderror_variable(source,indexsource,index);
};
void clTecplot::rpderror_variable(const clTecplot* const source,const int varindexsource,int varindex)
{
	if (_num_zonesharing || source->get_num_zonesharing()) fatal_err("Shared variables not yet implemented in clTecplot::rpderror_variable");

	if (varindex<0) varindex=varindexsource;
	if (varindex<0 || varindex>=_num_var) fatal_err("Illegal varindex in clTecplot::rpderror_variable");
	if (varindexsource<0 || varindexsource>=source->get_num_var()) fatal_err("Illegal varindexsource in clTecplot::rpderror_variable");
	if (source->get_num_zone()!=_num_zone) fatal_err("Inconsistent number of zones in clTecplot::rpderror_variable");
	if (strcmp(source->get_varname(varindexsource),get_varname(varindex))) fatal_err("Inconsistent variable name in clTecplot::rpderror_variable");

	int n_zone,n_data,num_data,idim,jdim,kdim;
	double d,dd,dummy;
	double *dataptr=NULL;
	const double *sourcedataptr=NULL;
	clTecplot_zone *zoneptr=NULL;
	const clTecplot_zone *sourcezoneptr=NULL;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		zoneptr=_zones[n_zone];
		if (zoneptr->get_num_varsharing()) fatal_err("Shared variables not yet implemented in clTecplot::rpderror_variable");
		sourcezoneptr=source->get_zoneptr(n_zone);
		if (sourcezoneptr->get_num_varsharing()) fatal_err("Shared variables not yet implemented in clTecplot::rpderror_variable");
		zoneptr->_varmin.deallocate();
		zoneptr->_varmax.deallocate();
		if (zoneptr->get_varlocation(varindex)!=sourcezoneptr->get_varlocation(varindexsource)) fatal_err("Inconsistent varlocation in clTecplot::rpderror_variable");
		if (zoneptr->zonetype_ordered()) {
			idim=zoneptr->get_idim();
			jdim=zoneptr->get_jdim();
			kdim=zoneptr->get_kdim();
			if (idim!=sourcezoneptr->get_idim()) fatal_err("Inconsistent idim  in clTecplot::rpderror_variable");
			if (jdim!=sourcezoneptr->get_jdim()) fatal_err("Inconsistent jdim  in clTecplot::rpderror_variable");
			if (kdim!=sourcezoneptr->get_kdim()) fatal_err("Inconsistent idim  in clTecplot::rpderror_variable");
			num_data=idim*jdim*kdim;
		} else {
			if (zoneptr->get_varlocation(varindex)) {
				num_data=zoneptr->get_num_elem();
				if (num_data!=sourcezoneptr->get_num_elem()) fatal_err("Inconsistent number of elements in clTecplot::rpderror_variable");
			} else {
				num_data=zoneptr->get_num_node();
				if (num_data!=sourcezoneptr->get_num_node()) fatal_err("Inconsistent number of nodes in clTecplot::rpderror_variable");
			}
		}
		dataptr=zoneptr->get_datastartptr(varindex);
		sourcedataptr=sourcezoneptr->get_datastartptr(varindexsource);
		for (n_data=0; n_data<num_data; n_data++,dataptr++,sourcedataptr++) {
			d=clUtils_namespace::abs(*dataptr);
			dd=clUtils_namespace::abs(*sourcedataptr);
			dummy=(d>dd ? d : dd);
			if (dummy) (*dataptr)=2.0*((*dataptr)-(*sourcedataptr))/dummy;
			else (*dataptr)=0.0;
		}
	}
	clString_namespace::clString varname("rpderror(");
	varname+=get_varname(varindex);
	varname+=")";
	set_varname(varindex,varname.get_string());
};

/* Replace variable "varnamedest" with the product of "varnamedest" and "varnamemult"*/
void clTecplot::multiply_variable(const char* const varnamemult,const char* const varnamedest)
{
	int n_var;
	int indexmult=-1,indexdest=-1;
	for (n_var=0; n_var<_num_var; n_var++) {
		if (!strcmp(get_varname(n_var),varnamemult)) indexmult=n_var;
		if (!strcmp(get_varname(n_var),varnamedest)) indexdest=n_var;
	}
	if (indexmult<0 || indexdest<0) fatal_err("Unable to get variable index in clTecplot::relerror_variable");
	multiply_variable(indexmult,indexdest);
};
void clTecplot::multiply_variable(const int varindexmult,const int varindexdest)
{
	if (_num_zonesharing) fatal_err("Shared variables not yet implemented in clTecplot::multiply_variable");

	if (varindexmult<0 || varindexmult>=_num_var) fatal_err("Illegal varindexmult in clTecplot::multiply_variable");
	if (varindexdest<0 || varindexdest>=_num_var) fatal_err("Illegal varindexdest in clTecplot::multiply_variable");

	int n_zone,n_data,num_data,idim,jdim,kdim;
	double *destdataptr=NULL;
	const double *multdataptr=NULL;
	clTecplot_zone *zoneptr=NULL;
	for (n_zone=0; n_zone<_num_zone; n_zone++) {
		zoneptr=_zones[n_zone];
		if (zoneptr->_num_varsharing) fatal_err("Shared variables not yet implemented in clTecplot::multiply_variable");
		if (zoneptr->get_varlocation(varindexmult)!=zoneptr->get_varlocation(varindexdest)) 
			fatal_err("Inconsistent varlocation in clTecplot::multiply_variable");
		zoneptr->_varmin.deallocate();
		zoneptr->_varmax.deallocate();
		if (zoneptr->zonetype_ordered()) {
			idim=zoneptr->get_idim();
			jdim=zoneptr->get_jdim();
			kdim=zoneptr->get_kdim();
			num_data=idim*jdim*kdim;
		} else {
			if (zoneptr->get_varlocation(varindexdest)) num_data=zoneptr->get_num_elem();
			else num_data=zoneptr->get_num_node();
		}
		destdataptr=zoneptr->get_datastartptr(varindexdest);
		multdataptr=zoneptr->get_datastartptr(varindexmult);
		for (n_data=0; n_data<num_data; n_data++) (*destdataptr++)*=(*multdataptr++);
	}
	clString_namespace::clString varname(get_varname(varindexdest));
	varname+="*";
	varname+=get_varname(varindexmult);
	set_varname(varindexdest,varname.get_string());
};
