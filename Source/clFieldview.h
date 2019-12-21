/*******************************************************************************
  A collection of Fieldview usage functions
  Author: Rolf Verberg
*******************************************************************************/

#ifndef clFieldview_H__INCLUDED_
#define clFieldview_H__INCLUDED_

#include <stdio.h>
#include <string.h>
#include "clUtils.h"
#include "clTecplot.h"

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

class clFieldview_grid
{
	friend class clFieldview;

	private:

		struct standard_face {
			int _num_face;
			int **_conn;
			standard_face *_next;
		};

		int _filetype;							// Filetype: FV_GRIDS_FILE, FV_RESULTS_FILE or FV_COMBINED_FILE
		int _num_facetype;						// Number of boundary face types
		int _num_var;							// Number of function variables
		int _num_bvar;							// Number of boundary variables
		int _num_node;							// Number of nodes
		clUtils_namespace::my_float *_xcoord;	// Nodal x-coordinates
		clUtils_namespace::my_float *_ycoord;	// Nodal y-coordinates
		clUtils_namespace::my_float *_zcoord;	// Nodal z-coordinates
		clUtils_namespace::my_float **_func;	// Function variables
		standard_face **_standard_face;			// Standard faces
		clFieldview_grid *_next;				// Pointer to next grid

	public:

		clFieldview_grid(void) :
			_filetype(0),
			_num_facetype(0),
			_num_var(0),
			_num_bvar(0),
			_num_node(0),
			_xcoord(NULL),
			_ycoord(NULL),
			_zcoord(NULL),
			_func(NULL),
			_standard_face(NULL),
			_next(NULL)
		{ };
		clFieldview_grid(const int filetype,const int num_facetype,const int num_var,const int num_bvar) :
			_num_node(0),
			_xcoord(NULL),
			_ycoord(NULL),
			_zcoord(NULL),
			_func(NULL),
			_next(NULL) 
		{
			_filetype=filetype;
			_num_facetype=num_facetype;
			_num_var=num_var;
			_num_bvar=num_bvar;
			_standard_face=new standard_face* [_num_facetype];
			for (int i=0; i<_num_facetype; i++) _standard_face[i]=NULL;
		};
		virtual ~clFieldview_grid(void) {
			if (_xcoord) delete [] _xcoord;
			if (_ycoord) delete [] _ycoord;
			if (_zcoord) delete [] _zcoord;
			if (_func) {
				delete [] _func[0];
				delete [] _func;
			}
			if (_standard_face) {
				standard_face *standard_face_ptr=NULL;
				standard_face *standard_face_nextptr=NULL;
				for (int i=0; i<_num_facetype; i++) if (_standard_face[i]) {
					standard_face_ptr=_standard_face[i];
					for (;; standard_face_ptr=standard_face_ptr->_next) {
						if (standard_face_ptr->_conn) {
							delete [] standard_face_ptr->_conn[0];
							delete [] standard_face_ptr->_conn;
						}
						standard_face_nextptr=standard_face_ptr->_next;
						delete standard_face_ptr;
						if (!standard_face_nextptr) break;
						standard_face_ptr=standard_face_nextptr;
					}
				}
				delete [] _standard_face;
			}
		};

	private:

		void read_grid(FILE*,const char*);
		void read_standard_face(FILE*,const char*);
		void read_arbitrary_face(FILE*,const char*);
		void read_standard_element(FILE*,const char*);
};

class clFieldview
{
	friend class clFieldview_grid;

	private:

		static const int IO_LINE_LENGTH;

		/* Numeric tags (codes) for FIELDVIEW binary file format. */
		static const int FV_MAGIC;

		/* Content of the file (grid only, results only or combined). */
		static const int FV_GRIDS_FILE;
		static const int FV_RESULTS_FILE;
		static const int FV_COMBINED_FILE;

		static const int FV_NODES;
		static const int FV_FACES;
		static const int FV_ELEMENTS;
		static const int FV_VARIABLES;
		static const int FV_BNDRY_VARS;
		static const int FV_ARB_POLY_FACES;
		static const int FV_ARB_POLY_ELEMENTS;
		static const int FV_ARB_POLY_BNDRY_VARS;

		static const int FV_TET_ELEM_ID;
		static const int FV_HEX_ELEM_ID;
		static const int FV_PRISM_ELEM_ID;
		static const int FV_PYRA_ELEM_ID;
		static const int FV_ARB_POLY_ELEM_ID;

	private:

		int _filetype;					// Filetype: FV_GRIDS_FILE, FV_RESULTS_FILE or FV_COMBINED_FILE
		int _num_grid;					// Number of grids
		int _num_facetype;				// Number of boundary face types
		int _num_var;					// Number of function variables
		int _num_bvar;					// Number of boundary variables
		int *_faceresultflag;			// Boundary face result flag
		int *_facenormalflag;			// Boundary face normal flag
		double _version;				// Fieldview version of data file
		double _time;					// Solution time
		double _mach;					// Solution Mach number
		double _alpha;					// Solution alpha
		double _re;						// Solution Reynolds number
		char **_facetypename;			// Face type names
		char **_varname;				// Function variable names
		char **_bvarname;				// Boundary variable names
		clFieldview_grid *_firstgrid;	// Pointer to first grid
		clFieldview_grid **_grids;		// Array of pointers to all grids

	public:

		clFieldview(void) :
			_filetype(0),
			_num_grid(0),
			_num_facetype(0),
			_num_var(0),
			_num_bvar(0),
			_faceresultflag(NULL),
			_facenormalflag(NULL),
			_version(3.0),
			_time(0.0),
			_mach(0.0),
			_alpha(0.0),
			_re(0.0),
			_facetypename(NULL),
			_varname(NULL),
			_bvarname(NULL),
			_firstgrid(NULL),
			_grids(NULL)
		{ };

		virtual ~clFieldview(void) {
			if (_faceresultflag) delete [] _faceresultflag;
			if (_facenormalflag) delete [] _facenormalflag;
			if (_facetypename) {
				for (int i=0; i<_num_facetype; i++) delete [] _facetypename[i];
				delete [] _facetypename;
			}
			if (_varname) {
				for (int i=0; i<_num_var; i++) delete [] _varname[i];
				delete [] _varname;
			}
			if (_bvarname) {
				for (int i=0; i<_num_bvar; i++) delete [] _bvarname[i];
				delete [] _bvarname;
			}
			if (_grids) {
				for (int i=0; i<_num_grid; i++) if (_grids[i]) delete _grids[i];
				if (_grids) delete [] _grids;
			}
		};

		void read(const char *filename);
//		inline void convert_to_tecplot(clTecplot *tecplotgrid,const int surfvol_flag) {
//			if (!surfvol_flag) convert_to_tecplot_surface(tecplotgrid);
//			else convert_to_tecplot_volume(tecplotgrid);
//		}
		void convert_to_tecplot_surface(clTecplot*);
		void convert_to_tecplot_volume(clTecplot*);
		void add_result_to_tecplot_surface(clTecplot*);

	private:

		int fread_str80(FILE*,char*);
		void read_header(FILE*,const char*);

};

#endif // clFieldview_H__INCLUDED_
