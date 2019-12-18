/*******************************************************************************
  A collection of Tecplot usage functions
  Author: Rolf Verberg
*******************************************************************************/

#ifndef CLTECPLOT_H__INCLUDED_
#define CLTECPLOT_H__INCLUDED_

#include <stdio.h>
#include <string.h>
#include "clMatrix.h"
#include "clMatrix3D.h"
#include "clString.h"
#include "clUtils.h"
#include "clVector.h"

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

class clNton_table
{
	public:

		struct nton_entry {
			int _num_node;
			int *_nodes;
		};

		struct ntondist_entry {
			int _num_node;
			int *_nodes;
			double _sum_r2inv_inv;
//			double *_r;
			double *_r2inv;
		};

		int _num_node;				// Number of nodes in table
		nton_entry **_nton;			// Node to nodes table
		ntondist_entry **_ntondist;	// Node to nodes and distance table

		clNton_table(void) :
			_num_node(0),
			_nton(NULL),
			_ntondist(NULL)
		{ };
		clNton_table(const int);
		clNton_table(const int,const int);

		virtual ~clNton_table(void)	{
			if (_nton) {
				for (int n_node=0; n_node<_num_node; n_node++) if (_nton[n_node]) {
					if (_nton[n_node]->_nodes) delete [] _nton[n_node]->_nodes;
					delete _nton[n_node];
				}
				delete [] _nton;
			}
			if (_ntondist) {
				for (int n_node=0; n_node<_num_node; n_node++) if (_ntondist[n_node]) {
					if (_ntondist[n_node]->_nodes) delete [] _ntondist[n_node]->_nodes;
					if (_ntondist[n_node]->_r2inv) delete [] _ntondist[n_node]->_r2inv;
					delete _ntondist[n_node];
				}
				delete [] _ntondist;
			}
		};
		
		void nton_entry_init(nton_entry*,const int);
		void ntondist_entry_init(ntondist_entry*,const int);
		void construct_nton_table(const int,const int,const int,const int* const);
		void construct_nton_table(const int,const int,const double* const* const,const double range=0.0,const int maxnode=10);
};

class clTecplot_zone
{
	friend class clTecplot;

	private:

		int _zonetype;			// Zone type
		int _datatype;			// Zone data style: POINT or BLOCK
		int _num_var;			// Number of variables
		int _num_varsharing;	// Number of shared variables
		clVector_namespace::clVector<int> _varsharingzone;	// Zone numbers to share variables with
		clVector_namespace::clVector<double> _varmin;		// Minimum values of variables
		clVector_namespace::clVector<double> _varmax;		// Maximum values of variables
		clString_namespace::clString _title;				// Zone title
		clTecplot_zone *_next;	// Pointer to next zone

	public:

		clTecplot_zone(void) :
			_zonetype(0),
			_datatype(0),
			_num_var(0),
			_num_varsharing(0),
			_next(NULL)
		{ };
		clTecplot_zone(const int num_var) :
			_zonetype(0),
			_datatype(0),
			_num_varsharing(0),
			_next(NULL) 
		{
			_num_var=num_var;
			_varsharingzone.allocate(_num_var,-1);
		};

		virtual ~clTecplot_zone(void) { };

		int set_zonetype_num_conn(const char* const);
		int set_zonetype_num_conn(const int);
		int get_zonetype(void) const {return _zonetype;};
		bool zonetype_ordered(void) const;
//		bool zonetype_fequadrilateral(void) const;
//		bool zonetype_fepolygon(void) const;
//		bool zonetype_fepolyhedron(void) const;

		int set_datatype(const char* const);
		void set_datatype(const int datatype) {_datatype=datatype;};
		int get_datatype(void) const {return _datatype;};
		bool datatype_point(void) const;

		void set_num_var(const int num_var,clVector_namespace::clVector<int> &varlocation,
			clVector_namespace::clVector<int> &varsharingzone) {
			if (_num_var && _num_var!=num_var) {
				_varsharingzone.deallocate();
				dealloc_varlocation();
			}
			if (varlocation.get_dim()) {
				if (varlocation.get_dim()!=num_var) 
					clUtils_namespace::fatal_err(0,"inconsistent varlocation in set_num_var::clTecplot_zone");
				set_varlocation(varlocation);
			}
			if (!varsharingzone.get_dim()) {
				_varsharingzone.allocate(num_var,-1);
			} else {
				if (varsharingzone.get_dim()!=num_var) 
					clUtils_namespace::fatal_err(0,"inconsistent varsharingzone in set_num_var::clTecplot_zone");
				_varsharingzone=varsharingzone;
			}
			_num_var=num_var;
		};
		void set_num_var(const int num_var) {
			if (_num_var && _num_var!=num_var) {
				if (_num_varsharing) clUtils_namespace::fatal_err("_num_varsharing in set_num_var::clTecplot_zone");
				if (get_num_var_cc() && get_varlocationdim()!=num_var) 
					clUtils_namespace::fatal_err("illegal _num_var_cc in set_num_var::clTecplot_zone");
				_varsharingzone.deallocate();
				dealloc_varlocation();
			}
			_num_var=num_var;
		};
		int get_num_var(void) const {return _num_var;};

		int set_num_varsharing(void) {
			_num_varsharing=0;
			if (!_num_var) return 0;
			if (!_varsharingzone.get_dim()) return 1;
			for (int n_var=0; n_var<_num_var; n_var++) if (_varsharingzone(n_var)!=-1) _num_varsharing++;
			return 1;
		};
		void set_num_varsharing(const int num_varsharing) {_num_varsharing=num_varsharing;};
		int get_num_varsharing(void) const {return _num_varsharing;};
		int alloc_varsharingzone(void) {
			if (!_num_var) return 0;
			_varsharingzone.allocate(_num_var,-1);
			return 1;
		};
		void dealloc_varsharingzone(void) {_varsharingzone.deallocate();};
		int reset_varsharingzone(void) {
			if (!_num_var || !_varsharingzone.get_dim()) return 0;
			_varsharingzone=-1;
			_num_varsharing=0;
			return 1;
		};
		void set_varsharingzone(const int n_zone) {
			if (!_varsharingzone.get_dim()) alloc_varsharingzone();
			_varsharingzone=n_zone;
			set_num_varsharing();
		};
		int set_varsharingzone(const int n_var,const int zone) {
			if (!_varsharingzone.get_dim()) alloc_varsharingzone();
			_varsharingzone(n_var)=zone;
			set_num_varsharing();
			return 0;
		};
		int* get_varsharingzoneptr(void) {return _varsharingzone.get_startptr();};
		const int* get_varsharingzoneptr(void) const {return _varsharingzone.get_startptr();};
		int get_varsharingzone(const int n_var) const {
			if (!_varsharingzone.get_dim()) return -1;
			else return _varsharingzone(n_var);
		};

		const double* get_varminptr(void) const {return _varmin.get_startptr();};
		double get_varmin(const int n_var) const {return _varmin(n_var);};
		double get_varmin(const int n_var) {
			if (!_varmin.get_dim()) calculate_varrange();
			return _varmin(n_var);
		};
		const double* get_varmaxptr(void) const {return _varmax.get_startptr();};
		double get_varmax(const int n_var) const {return _varmax(n_var);};
		double get_varmax(const int n_var) {
			if (!_varmax.get_dim()) calculate_varrange();
			return _varmax(n_var);
		};

		void set_title(const char* const title) {_title=title;};
		int if_title(void) const {return _title.get_len();};
		const char* get_title(void) const {return _title.get_string();};

		virtual double* get_datastartptr(const int) {return NULL;};
		virtual const double* get_datastartptr(const int) const {return NULL;};

		virtual void set_idim(const int) {};
		virtual int get_idim(void) const {return 0;};
		virtual void set_jdim(const int) {};
		virtual int get_jdim(void) const {return 0;};
		virtual void set_kdim(const int) {};
		virtual int get_kdim(void) const {return 0;};

		virtual void set_num_node(const int) {};
		virtual int get_num_node(void) const {return 0;};
		virtual void set_num_elem(const int) {};
		virtual int get_num_elem(void) const {return 0;};
		virtual void set_num_conn(const int) {};
		virtual int get_num_conn(void) const {return 0;};
		virtual void set_num_var_cc(const int) {};
		virtual int get_num_var_cc(void) const {return 0;};
		virtual void dealloc_varlocation(void) {};
		virtual void set_varlocation(const int) {};
		virtual void set_varlocation(const clVector_namespace::clVector<int>&) {};
		virtual int get_varlocation(const int) const {return 0;};
		virtual int get_varlocationdim(void) const {return 0;};

		virtual void allocate_data(void) {};
		virtual void set_data(clVector_namespace::clMatrix3D<double>*) {};
		virtual void set_data(clVector_namespace::clMatrix3D<double>&,const int) {};
		virtual void set_data(clVector_namespace::clVector<double>*) {};
		virtual void set_data(clVector_namespace::clVector<double>&,const int) {};

		virtual int* get_connstartptr(void) {return NULL;};
		virtual const int* get_connstartptr(void) const {return NULL;};
		virtual clVector_namespace::clMatrix<int> *get_connptr(void) {return NULL;};
		virtual void set_conn(clVector_namespace::clMatrix<int>&) {};

		void add_data(const clTecplot_zone* const);
		void add_data(const clTecplot_zone* const,const clVector_namespace::clVector<int>&);
		void translate_data(const double);
		void translate_data(const double,const int);
		void translate_data(const double,const clVector_namespace::clVector<int>&);
		void multiply_data(const double);
		void multiply_data(const double,const int);
		void multiply_data(const double,const clVector_namespace::clVector<int>&);
		double find_min_vertex(const double mincutoff=0.0) const;
		double find_max_vertex(const double mincutoff=0.0) const;
		void find_minmax_vertex(double&,double&,const double mincutoff=0.0) const;
		void calculate_varrange(void);
		void copy_zoneheader(const clTecplot_zone* const);
		void create_nton_table(clNton_table**,const double range=0.0) const;
		void create_nton_table(clNton_table**,const int,const double range=0.0) const;
		void create_nton_table(clNton_table**,FILE*,const int,const double range=0.0) const;
		void create_nton_table(clNton_table**,const int,FILE*,const int,const double range=0.0) const;
		void find_nearest_neighbor(const clTecplot_zone* const,const clNton_table* const,clNton_table**) const;
		void find_nearest_neighbor(const clTecplot_zone* const,const clNton_table* const,clNton_table**,const double eps) const;
		void find_nearest_neighbor(const clTecplot_zone* const,const clNton_table* const,clNton_table**,
			FILE*,const int) const;
		void find_nearest_neighbor(const clTecplot_zone* const,const clNton_table* const,clNton_table**,
			const double,FILE*,const int) const;
		void find_nearest_neighbor(const clTecplot_zone* const,const clNton_table* const,clNton_table**,
			const double,const int,const double,const double,const double,const double,
			const clVector_namespace::clVector<int>&,const clVector_namespace::clMatrix<int>&,FILE*,const int) const;
		void find_nearest_neighbor_cc(const clTecplot_zone* const,const clNton_table* const,clNton_table**) const;
		void find_nearest_neighbor_cc(const clTecplot_zone* const,const clNton_table* const,clNton_table**,const double) const;
		void find_nearest_neighbor_cc(const clTecplot_zone* const,const clNton_table* const,clNton_table**,FILE*,const int) const;
		void find_nearest_neighbor_cc(const clTecplot_zone* const,const clNton_table* const,clNton_table**,const double,FILE*,const int) const;
		void interpolate(clTecplot_zone* const) const;
		void interpolate(clTecplot_zone* const,const double) const;
		void interpolate(clTecplot_zone* const,const double,const char* const) const;
		void interpolate(clTecplot_zone* const,const double,const clString_namespace::clString&) const;
		void interpolate(clTecplot_zone* const,const double,const double) const;
		void interpolate(clTecplot_zone* const,const double,const double,const clString_namespace::clString&) const;
		void interpolate(clTecplot_zone* const,const clVector_namespace::clVector<int>&) const;
		void interpolate(clTecplot_zone* const,const double,const clVector_namespace::clVector<int>&) const;
		void interpolate(clTecplot_zone* const,const double,const double,const clVector_namespace::clVector<int>&) const;
		void interpolate(clTecplot_zone* const,const double,const double,const clVector_namespace::clVector<int>&,
			const char* const) const;
		void interpolate(clTecplot_zone* const,const double,const double,const clVector_namespace::clVector<int>&,
			const clString_namespace::clString&) const;
		void interpolate_cc(clTecplot_zone* const) const;
		void interpolate_cc(clTecplot_zone* const,const clVector_namespace::clVector<int>&) const;
		void interpolate_cc(clTecplot_zone* const,const double,const double,const clVector_namespace::clVector<int>&) const;
		void interpolate_cc(clTecplot_zone* const,const double,const double,const clVector_namespace::clVector<int>&,
			const char* const) const;
		void interpolate_cc(clTecplot_zone* const,const double,const double,const clVector_namespace::clVector<int>&,
			const clString_namespace::clString&) const;

	private:

		void print_zonetype(FILE*) const;
		void write_datatype(FILE*) const;
};

class clTecplot_zone_structured : public clTecplot_zone
{
	protected:

		int _idim,_jdim,_kdim;	// Zone dimensions
		clVector_namespace::clMatrix3D<double> *_data;	// Variable fields, each variable in a vector

	public:

		clTecplot_zone_structured(void) :
			_idim(0),
			_jdim(0),
			_kdim(0),
			_data(NULL)
		{ };
		clTecplot_zone_structured(const int num_var) : clTecplot_zone(num_var)
		{
			_idim=0;
			_jdim=0;
			_kdim=0;
			_data=NULL;
		};

		virtual ~clTecplot_zone_structured(void) {
			deallocate_data();
		};

		double* get_datastartptr(const int n_var) {
			return ((_data && n_var>=0 && n_var<get_num_var()) ? _data[n_var].get_startptr() : NULL);
		};
		const double* get_datastartptr(const int n_var) const {
			return ((_data && n_var>=0 && n_var<get_num_var()) ? _data[n_var].get_startptr() : NULL);
		};

		void set_idim(const int idim) {_idim=idim;};
		int get_idim(void) const {return _idim;};
		void set_jdim(const int jdim) {_jdim=jdim;};
		int get_jdim(void) const {return _jdim;};
		void set_kdim(const int kdim) {_kdim=kdim;};
		int get_kdim(void) const {return _kdim;};

		void allocate_data(void) {
			deallocate_data();
			const int num_var=get_num_var();
			_data=new clVector_namespace::clMatrix3D<double> [num_var];
			if (!_data) clUtils_namespace::fatal_err(0,"_data in clTecplot_zone_structured::allocate_data");
			for (int n_var=0; n_var<num_var; n_var++) {
				_data[n_var].allocate(_kdim,_jdim,_idim);
			}
		};
		void set_data(clVector_namespace::clMatrix3D<double> *data) {
			const int num_var=get_num_var();
			if (!num_var) clUtils_namespace::fatal_err("_num_var not initialized in clTecplot_zone_structured::set_data");
			if (!_data) {
				_data=new clVector_namespace::clMatrix3D<double> [num_var];
				if (!_data) clUtils_namespace::fatal_err(0,"_data in clTecplot_zone_structured::set_data");
			}
			for (int n_var=0; n_var<num_var; n_var++) {
				_data[n_var].deallocate();
				_data[n_var].move(data[n_var]);
			}
		};
		void set_data(clVector_namespace::clMatrix3D<double> &data,const int n_var) {
			const int num_var=get_num_var();
			if (n_var<0 || n_var>=num_var) clUtils_namespace::fatal_err("Illegal n_var in clTecplot_zone_structured::set_data");
			if (!_data) {
				_data=new clVector_namespace::clMatrix3D<double> [get_num_var()];
				if (!_data) clUtils_namespace::fatal_err(0,"_data in clTecplot_zone_structured::set_data");
			}
			_data[n_var].deallocate();
			_data[n_var].move(data);
		};
		void set_data(clVector_namespace::clVector<double>*) {return;};
		void set_data(clVector_namespace::clVector<double>&,const int) {return;};

		void add_variables(const int);
		void remove_variable(const int,const int vvarindexend=-1);

	private:

		void deallocate_data(void) {
			if (_data) delete [] _data;
			_data=NULL;
		};
};

class clTecplot_zone_unstructured : public clTecplot_zone
{
	protected:

		int _num_node;			// Number of nodes
		int _num_elem;			// Number of elements
		int _num_conn;			// Number of nodes per element
		int _num_var_cc;		// Number of cell-centered variables
		clVector_namespace::clVector<int> _varlocation;	// Variable location 0: nodal (default) or 1: cell-centered
		clVector_namespace::clMatrix<int> _conn;		// Connectivities
		clVector_namespace::clVector<double> *_data;	// Variable fields, each variable in a vector

	public:

		clTecplot_zone_unstructured(void) :
			_num_node(0),
			_num_elem(0),
			_num_conn(0),
			_num_var_cc(0),
			_data(NULL)
		{ };
		clTecplot_zone_unstructured(int num_var) : clTecplot_zone(num_var)
		{
			_num_node=0;
			_num_elem=0;
			_num_conn=0;
			_num_var_cc=0;
//			_varlocation.allocate(num_var,0);
			_data=NULL;
		};

		virtual ~clTecplot_zone_unstructured(void) {
			deallocate_data();
		};

		double* get_datastartptr(const int n_var) {
			return ((_data && n_var>=0 && n_var<get_num_var()) ? _data[n_var].get_startptr() : NULL);
		};
		const double* get_datastartptr(const int n_var) const {
			return ((_data && n_var>=0 && n_var<get_num_var()) ? _data[n_var].get_startptr() : NULL);
		};

		void set_num_node(const int num_node) {_num_node=num_node;};
		int get_num_node(void) const {return _num_node;};
		void set_num_elem(const int num_elem) {_num_elem=num_elem;};
		int get_num_elem(void) const {return _num_elem;};
		void set_num_conn(const int num_conn) {_num_conn=num_conn;};
		int get_num_conn(void) const {return _num_conn;};
		void set_num_var_cc(const int num_var_cc) {_num_var_cc=num_var_cc;};
		int get_num_var_cc(void) const {return _num_var_cc;};
		void dealloc_varlocation(void) {_varlocation.deallocate();};
		void set_varlocation(const int) {
			clUtils_namespace::fatal_err("void set_varlocation(const int) not implemented in clTecplot_zone_unstructured");
		};
		void set_varlocation(const clVector_namespace::clVector<int> &varlocation) {
			_varlocation=varlocation;
			_num_var_cc=0;
			if (!_varlocation.get_dim()) return;
			for (int n_var=0; n_var<_varlocation.get_dim(); n_var++) if (_varlocation(n_var)) _num_var_cc++;
		};
		int get_varlocation(const int n_var) const {return (_varlocation.get_dim() ? _varlocation(n_var) : 0);};
		int get_varlocationdim(void) const {return _varlocation.get_dim();};

		void allocate_data(void) {
			deallocate_data();
			const int num_var=get_num_var();
			if (!_varlocation.get_dim()) _varlocation.allocate(num_var,0);
			else if (_varlocation.get_dim()!=num_var) clUtils_namespace::fatal_err(0,"Inconsistent varlocation size in clTecplot_zone_unstructured::allocate_data");
			_data=new clVector_namespace::clVector<double> [num_var];
			if (!_data) clUtils_namespace::fatal_err(0,"_data in clTecplot_zone_unstructured::allocate_data");
			for (int n_var=0; n_var<num_var; n_var++) {
				if (get_varsharingzone(n_var)==-1) {
					if (_varlocation(n_var)) _data[n_var].allocate(_num_elem,0);
					else _data[n_var].allocate(_num_node,0);
				} else {
					if (_varlocation(n_var)) clUtils_namespace::fatal_err("Shared cell-centered variable in clTecplot_zone_unstructured::allocate_data");
				}
			}
		};
		void set_data(clVector_namespace::clMatrix3D<double>*) {return;};
		void set_data(clVector_namespace::clMatrix3D<double>&,const int) {return;};
		void set_data(clVector_namespace::clVector<double> *data) {
			const int num_var=get_num_var();
			if (!num_var) clUtils_namespace::fatal_err("_num_var not initialized in clTecplot_zone_unstructured::set_data");
			if (!_data) {
				_data=new clVector_namespace::clVector<double> [num_var];
				if (!_data) clUtils_namespace::fatal_err(0,"_data in clTecplot_zone_unstructured::set_data");
			}
			for (int n_var=0; n_var<num_var; n_var++) {
				_data[n_var].deallocate();
				_data[n_var].move(data[n_var]);
			}
		};
		void set_data(clVector_namespace::clVector<double> &data,const int n_var) {
			const int num_var=get_num_var();
			if (n_var<0 || n_var>=num_var) clUtils_namespace::fatal_err("Illegal n_var in clTecplot_zone_unstructured::set_data");
			if (!_data) {
				_data=new clVector_namespace::clVector<double> [get_num_var()];
				if (!_data) clUtils_namespace::fatal_err(0,"_data in clTecplot_zone_unstructured::set_data");
			}
			_data[n_var].deallocate();
			_data[n_var].move(data);
		};

		int* get_connstartptr(void) {return _conn.get_startptr();};
		const int* get_connstartptr(void) const {return _conn.get_startptr();};
		clVector_namespace::clMatrix<int> *get_connptr(void) {return &_conn;};
		void set_conn(clVector_namespace::clMatrix<int> &conn) {
			_conn.deallocate();
			_conn.move(conn);
		};

		void add_variables(const int,const int varlocation=0);
		void remove_variable(const int,const int vvarindexend=-1);

	private:

		void deallocate_data(void) {
			if (_data) delete [] _data;
			_data=NULL;
		};
};

class clTecplot_zone_readheader : public clTecplot_zone
{
	protected:

		int _idim,_jdim,_kdim;	// Zone dimensions
		int _num_node;			// Number of nodes
		int _num_elem;			// Number of elements
		int _num_conn;			// Number of nodes per element
		int _num_var_cc;		// Number of cell-centered variables
		clVector_namespace::clVector<int> _varlocation;	// Variable location 0: nodal (default) or 1: cell-centered

	public:

		clTecplot_zone_readheader(void) :
			_idim(0),
			_jdim(0),
			_kdim(0),
			_num_node(0),
			_num_elem(0),
			_num_conn(0),
			_num_var_cc(0)
		{ };
		clTecplot_zone_readheader(int num_var) : clTecplot_zone(num_var)
		{
			_idim=0;
			_jdim=0;
			_kdim=0;
			_num_node=0;
			_num_elem=0;
			_num_conn=0;
			_num_var_cc=0;
			_varlocation.allocate(num_var,0);
		};

		virtual ~clTecplot_zone_readheader(void) { };

		void set_idim(const int idim) {_idim=idim;};
		int get_idim(void) const {return _idim;};
		void set_jdim(const int jdim) {_jdim=jdim;};
		int get_jdim(void) const {return _jdim;};
		void set_kdim(const int kdim) {_kdim=kdim;};
		int get_kdim(void) const {return _kdim;};
		void set_num_node(const int num_node) {_num_node=num_node;};
		int get_num_node(void) const {return _num_node;};
		void set_num_elem(const int num_elem) {_num_elem=num_elem;};
		int get_num_elem(void) const {return _num_elem;};
		void set_num_conn(const int num_conn) {_num_conn=num_conn;};
		int get_num_conn(void) const {return _num_conn;};
		void incr_num_var_cc(const int numincr=1) {_num_var_cc+=numincr;};
		int get_num_var_cc(void) const {return _num_var_cc;};
		void set_varlocation(const int n_var) {
			const int num_var=get_num_var();
			if (n_var<0 || n_var>=num_var) clUtils_namespace::fatal_err(0,"Illegal n_var in clTecplot_zone_readheader::set_varlocation");
			if (!_varlocation.get_dim()) _varlocation.allocate(num_var,0);
			else if (_varlocation.get_dim()!=num_var) clUtils_namespace::fatal_err(0,"Inconsistent varlocation size in clTecplot_zone_readheader::set_varlocation");
			_varlocation(n_var)=1;
		};
		void set_varlocation(const clVector_namespace::clVector<int>&) {
			clUtils_namespace::fatal_err("void set_varlocation(const clVector_namespace::clVector<int>&) not implemented in clTecplot_zone_readheader");
		};
		int get_varlocation(const int n_var) const {return (_varlocation.get_dim() ? _varlocation(n_var) : 0);};
};

class clTecplot
{
	friend class clTecplot_zone;

	private:

		enum zonetypes {ORDERED,FELINESEG,FETRIANGLE,FEQUADRILATERAL,FETETRAHEDRON,FEBRICK,FEPOLYGON,FEPOLYHEDRON};
		enum datatypes {POINT,BLOCK};

		static const int LINE_LENGTH;
		static const int NAME_LENGTH;
		static const int READBUFFERSIZE;

	private:

		int _version;					// Tecplot version of binary (plt) file
		int _num_zone;					// Number of zones
		int _num_var;					// Number of variables
		int _num_zonesharing;			// Number of zones that share variables
		clString_namespace::clString _title;			// Title
		clVector_namespace::clVector<int> _zonesharing;	// Zones that share variables
		clString_namespace::clString_vector _varname;	// Variable titles
		clTecplot_zone *_firstzone;		// Pointer to first zone
		clTecplot_zone **_zones;		// Array of pointers to all zones

	public:

		clTecplot(void) :
			_version(102),
			_num_zone(0),
			_num_var(0),
			_num_zonesharing(0),
			_firstzone(NULL),
			_zones(NULL)
		{ };
		clTecplot(const int num_zone,const int num_var) :
			_version(102),
			_num_zonesharing(0),
			_firstzone(NULL)
		{
			_num_zone=num_zone;
			_num_var=num_var;
			_varname.allocate(_num_var);
			_zones=new clTecplot_zone* [_num_zone];
			if (!_zones) clUtils_namespace::fatal_err(0,"_zones in clTecplot");
			for (int n_zone=0; n_zone<_num_zone; n_zone++) _zones[n_zone]=NULL;
		};

		virtual ~clTecplot(void) {
			if (_zones) {
				for (int n_zone=0; n_zone<_num_zone; n_zone++) if (_zones[n_zone]) delete _zones[n_zone];
				delete [] _zones;
			}
		};

		/*  Default Copy Constructor  */
		clTecplot(const clTecplot &arg) :
			_version(arg._version),
			_num_zone(0),
			_num_var(0),
			_num_zonesharing(0),
			_firstzone(NULL),
			_zones(NULL)
		{
			copy(&arg);
		};

		/*  Default Assignment Constructor  */
		clTecplot& operator= (const clTecplot &arg)
		{
			if (this==&arg) return *this;

			copy(&arg);

			return *this;
		};

		int get_namelength(void) const {return NAME_LENGTH;};

		void copy_header(const clTecplot* const source)
		{
			_num_var=source->get_num_var();
			set_title(source->get_title());
			alloc_varname();
			for (int n_var=0; n_var<_num_var; n_var++) set_varname(n_var,source->get_varname(n_var));
		};

		void set_num_zone(const int num_zone) {_num_zone=num_zone;};
		int get_num_zone(void) const {return _num_zone;};

		void set_num_var(const int num_var) {
			if (_num_var && _num_var!=num_var) dealloc_varname();
			_num_var=num_var;
		}
		int get_num_var(void) const {return _num_var;};

		int get_num_zonesharing(void) const {return _num_zonesharing;};
		int get_zonesharing(const int n_zone) const {
			if (!_zonesharing.get_dim()) return -1;
			else return _zonesharing(n_zone);
		};
		const int* get_zonesharingptr(void) const {return _zonesharing.get_startptr();};
		void set_zonesharing(void);

		void set_title(const char* const title) {_title=title;};
		int if_title(void) const {return _title.get_len();};
		const char* get_title(void) const {return _title.get_string();};

		int alloc_varname(void) {
			if (!_num_var) return 0;
			_varname.allocate(_num_var);
			return 1;
		};
		int alloc_varname(const int num_var) {
			if (num_var<0) return 0;
			_num_var=num_var;
			return alloc_varname();
		};
		void dealloc_varname(void) {_varname.deallocate();};
		int set_varname(const int n_var,const char* const varname) {
			if (n_var<0 || n_var>=_num_var || varname==NULL) return 0;
			if (!_varname.get_dim()) if (!alloc_varname()) return 0;
			_varname(n_var)=varname;
			return 1;
		};
		char* get_varname(const int n_var) const {
			if (n_var<0 || n_var>=_num_var || !_varname.get_dim()) return NULL;
			else return _varname(n_var).get_string();
		};
		int get_varindex(const char* const varname) const {
			int index=-1;
			int num_var=get_num_var();
			for (int n_var=0; n_var<num_var; n_var++)
				if (!strcmp(get_varname(n_var),varname)) index=n_var;
			return index;
		};

		int alloc_zones(void) {
			if (!_num_zone) return 0;
			_zones=new clTecplot_zone* [_num_zone];
			if (!_zones) return 0;
			for (int n_zone=0; n_zone<_num_zone; n_zone++) _zones[n_zone]=NULL;
			_firstzone=NULL;
			return 1;
		};
		int alloc_zones(const int num_zone) {
			_num_zone=num_zone;
			return alloc_zones();
		};
		void dealloc_zones(void) {
			if (_zones) {
				for (int n_zone=0; n_zone<_num_zone; n_zone++) if (_zones[n_zone]) delete _zones[n_zone];
				if (_zones) delete [] _zones;
			}
			_zones=NULL;
			_firstzone=NULL;
		};
		void set_zone(const int n_zone,clTecplot_zone* const newzone) {
			if (!_num_zone) 
				fprintf(stderr,"FATAL ERROR: _num_zone=0 in clTecplot::set_zone");
			else if (n_zone<0 || n_zone>=_num_zone)
				fprintf(stderr,"FATAL ERROR: Illegal zone index in clTecplot::set_zone");
			if (!_zones) alloc_zones();
			_zones[n_zone]=newzone;
			if (!n_zone) _firstzone=newzone;
			else _zones[n_zone-1]->_next=newzone;
		};
		clTecplot_zone* get_zoneptr(const int n_zone) {return _zones[n_zone];};
		const clTecplot_zone* get_zoneptr(const int n_zone) const {return _zones[n_zone];};
		void remove_zone(const int,const int);
		void remove_zone(const int);
		void insert_zone(const int);
		void insert_zone(const int,const int);
		void copy_zone(const clTecplot_zone* const);
		void copy_zone(const int,const clTecplot_zone* const);
		void split_zone(const int,const int);

		void allocate(const int,const int);
		void read(const char* const);
		void read_ascii(const char* const);
		void write_ascii(const char* const,const char* const) const;
		void write_ascii(const char* const) const;
		void write_ascii(const char* const,const int,const int n_varupp=-1) const;
		void write_ascii(const char* const,const clVector_namespace::clVector<int>&) const;
		void write_ascii(const char* const,const clVector_namespace::clVector<int>&,const char* const) const;
		void write_ascii(const int,const char* const) const;
		void write_ascii(const int,const int,const char* const) const;
		void write_ascii(const int,const int,const char* const,const int,const int) const;
		void write_ascii(const clVector_namespace::clVector<int>&,const char* const) const;
		void write_ascii(const clVector_namespace::clVector<int>&,const char* const,
			const clVector_namespace::clVector<int>&) const;
		void read_bin(const char* const);
		void write(const char* const);
		void write(const char* const,const int);
		void write(const char* const,const int,const int);
		void write(const char* const,const clVector_namespace::clVector<int>&);
		void write(const char* const,const int,const clVector_namespace::clVector<int>&);
		void convert_fetriangle_conn(void);
		void convert_fetriangle_to_fequadrilateral(void);
		void convert_fequadrilateral_to_fetriangle(void);
		void eliminate_varsharing(void);
		void copy(const clTecplot* const);
		void add_data(const clTecplot* const);
		void add_data(const clTecplot* const,const clVector_namespace::clVector<int>&);
		void translate_data(const double);
		void translate_data(const double,const int);
		void translate_data(const double,const clVector_namespace::clVector<int>&);
		void multiply_data(const double);
		void multiply_data(const double,const int);
		void multiply_data(const double,const clVector_namespace::clVector<int>&);
		inline void mirror(const char* const varname)
		{
			int n_var;
			for (n_var=0; n_var<_num_var; n_var++) 
				if (!strcmp(_varname[n_var].get_string(),varname)) break;
			mirror(n_var);
		};
		void mirror(const int);
		double find_min_vertex(const double mincutoff=0.0) const;
		double find_max_vertex(const double mincutoff=0.0) const;
		void find_minmax_vertex(double&,double&,const double mincutoff=0.0) const;
		void create_nton_table(clNton_table**,const double range=0.0) const;
		void create_nton_table(clNton_table**,const char* const,const double range=0.0) const;
		void add_variable(void);
		void add_variables(const int,const int varlocation=0);
		void copy_variable(const clTecplot* const,const int);
		void copy_variable(const clTecplot* const,const int,const int);
		void remove_variable(const int,const int vvarindexend=-1);
		void abserror_variable(const clTecplot* const,const char* const);
		void abserror_variable(const clTecplot* const,const int,int varindexdest=-1);
		void relerror_variable(const clTecplot* const,const char* const);
		void relerror_variable(const clTecplot* const,const int,int varindexdest=-1);
		void rpderror_variable(const clTecplot* const,const char* const);
		void rpderror_variable(const clTecplot* const,const int,int varindexdest=-1);
		void multiply_variable(const char* const,const char* const);
		void multiply_variable(const int,const int);

	private:

		int get_num_var(const char* const) const;
		void read_fileheader(const char* const);
		int read_zoneheader(FILE*,clTecplot_zone_readheader* const);
		int read_title(FILE*);
		int read_zone(FILE*,clTecplot_zone**);
		void split_zone_unstructured(const int,const int);
		void read_geometries(FILE*);
		void read_text(FILE*);
		void read_customlabel(FILE*);
		void read_userrec(FILE*);
		int read_dataauxdata(FILE*);
		void read_varauxdata(FILE*);
		int read_datasection(FILE*,clTecplot_zone* const);
		int write_title(FILE*,const clVector_namespace::clVector<int>&) const;
		int write_zone(FILE*,const clTecplot_zone* const) const;
		void write_geometries(FILE*) const;
		void write_text(FILE*) const;
		void write_customlabel(FILE*) const;
		void write_userrec(FILE*) const;
		void write_dataauxdata(FILE*) const;
		void write_varauxdata(FILE*) const;
		int write_datasection(FILE*,clTecplot_zone* const,const clVector_namespace::clVector<int>&) const;
		int get_ltog_gtol_tables(const int,const int,const int,const int*,
			clVector_namespace::clVector<int>&,clVector_namespace::clVector<int>&) const;
};

#endif // CLTECPLOT_H__INCLUDED_
