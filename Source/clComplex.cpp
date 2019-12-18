/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/
#include "clComplex.h"

void clComplex_namespace::bcast(clVector_complex &vptr,const int ip)
{
#ifdef DIST_MPI
	int num_proc=clMpi_namespace::get_num_proc();
	int proc_id=clMpi_namespace::get_proc_id();
	if (num_proc<2) return;
	int dim;
	double *vector;
	if (proc_id==ip) dim=vptr.get_dim();
	clMpi_namespace::bcast(&dim,1,ip);
	if (!dim) {
		vptr.deallocate();
		return;
	}
	if (proc_id!=ip) vptr.allocate(dim);
	vector=vptr.get_startptr();
	const int dimvec=dim+dim;
	clMpi_namespace::bcast(vector,dimvec,ip);
#endif
};
