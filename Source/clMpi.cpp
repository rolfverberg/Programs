/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/
#include "clMpi.h"

namespace clMpi_namespace
{
#ifdef DIST_MPI
	MPI_Comm _mpi_comm=MPI_COMM_NULL;
#endif
	int _init_flag=0;
	int _num_proc=1;
	int _proc_id=0;
	int _n_warn=0;
}

void clMpi_namespace::mpi_init(int *argc,char ***argv)
{
	int num_proc=1,proc_id=0;
#ifdef DIST_MPI
	if (_mpi_comm!=MPI_COMM_NULL) fatal_err("MPI namespace cannot support multiple communicators yet: _mpi_comm must be NULL");
	_mpi_comm=MPI_COMM_WORLD;
	MPI_Init(argc,argv);
	MPI_Comm_size(_mpi_comm,&num_proc);
	MPI_Comm_rank(_mpi_comm,&proc_id);
#endif
	_num_proc=num_proc;
	_proc_id=proc_id;
	_init_flag=1;
};

#ifdef DIST_MPI
void clMpi_namespace::mpi_init(const int num_proc,const int proc_id, const MPI_Comm mpi_comm)
{
	if (_mpi_comm!=MPI_COMM_NULL) fatal_err("MPI namespace cannot support multiple communicators yet: _mpi_comm must be NULL");
	if (proc_id >= num_proc) return;

	_mpi_comm=mpi_comm;
	_num_proc=num_proc;
	_proc_id=proc_id;
	_init_flag=1;
	barrier();
};

MPI_Comm clMpi_namespace::get_mpi_comm(void)
{
	if (_mpi_comm==MPI_COMM_NULL) fatal_err("Must call clMpi_namespace::mpi_init before clMpi_namespace::get_mpi_comm");
	return(_mpi_comm);
};

MPI_Comm clMpi_namespace::mpi_create_comm(const int num_procs)
{
	int new_nproc=1,proc_id=0;
	MPI_Group world_group_id, group_ranks_id;
	MPI_Comm group_ranks_comm;
	int *group_ranks;

	MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);
	MPI_Comm_group(MPI_COMM_WORLD, &world_group_id);

	group_ranks=new int[num_procs];
	for (int i=0; i<num_procs; i++) group_ranks[i]=i;

	MPI_Group_incl(world_group_id,num_procs,group_ranks,&group_ranks_id);
	MPI_Comm_create(MPI_COMM_WORLD,group_ranks_id,&group_ranks_comm);

	delete group_ranks;

	if (proc_id<num_procs) {
		MPI_Comm_size(group_ranks_comm,&new_nproc);
		if (new_nproc!=num_procs) fatal_err("Error in create_mpi_comm: new_num_proc!=requested_num_proc");
	}

	return (group_ranks_comm);
};

void clMpi_namespace::mpi_clear(void)
{
	barrier();
	if (_mpi_comm != MPI_COMM_WORLD) mpi_free_comm();
	_mpi_comm=MPI_COMM_NULL;
	_num_proc=1;
	_proc_id=0;
	_init_flag=0;
};

void clMpi_namespace::mpi_free_comm(void)
{
	MPI_Comm_free(&_mpi_comm);
	_mpi_comm=MPI_COMM_NULL;
};
#endif

void clMpi_namespace::mpi_finalize(void)
{
#ifdef DIST_MPI
	MPI_Barrier(_mpi_comm);
	MPI_Finalize();
#endif
};

void clMpi_namespace::mpi_abort(int errorcode)
{
#ifdef DIST_MPI
	MPI_Abort(_mpi_comm,errorcode);
#endif
};

int clMpi_namespace::get_num_proc(void)
{
	if (!_init_flag) fatal_err("Must call clMpi_namespace::mpi_init before clMpi_namespace::get_num_proc");
	return (_num_proc);
};

int clMpi_namespace::get_proc_id(void)
{
	if (!_init_flag) fatal_err("Must call clMpi_namespace::mpi_init before clMpi_namespace::get_proc_id");
	return (_proc_id);
};

void clMpi_namespace::barrier(void)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	MPI_Barrier(_mpi_comm);
#endif
};

void clMpi_namespace::bcast(int* const iptr,const int n,const int ip)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	MPI_Bcast((void *)iptr,n,MPI_INT,ip,_mpi_comm);
#endif
};
void clMpi_namespace::bcast(double* const dptr,const int n,const int ip)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	MPI_Bcast((void *)dptr,n,MPI_DOUBLE,ip,_mpi_comm);
#endif
};
void clMpi_namespace::bcast(char* const cptr,const int n,const int ip)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	MPI_Bcast((void *)cptr,n,MPI_CHAR,ip,_mpi_comm);
#endif
};
/*void clMpi_namespace::bcast(bool* const bptr,const int n,const int ip)
{
#ifdef DIST_MPI
	MPI_Bcast((void *)bptr,n,MPI_C_BOOL,ip,_mpi_comm);
#endif
};*/

void clMpi_namespace::reduce_sum(const double* const dptr_in,double* const dptr_out,const int n,const int ip)
{
#ifdef DIST_MPI
	if (_num_proc>1) 
		MPI_Reduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_SUM,ip,_mpi_comm);
	else
		for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#else
	for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#endif
};
void clMpi_namespace::reduce_sum(double* const dptr_in,const int n,const int ip)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	double *dptr_out=new double [n];
	if (!dptr_out) fatal_err(0,"dptr_out in clMpi_namespace::reduce_sum");
	MPI_Reduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_SUM,ip,_mpi_comm);
	for (int i=0; i<n; i++) dptr_in[i]=dptr_out[i];
	delete [] dptr_out;
#endif
};

void clMpi_namespace::reduce_prod(const double* const dptr_in,double* const dptr_out,const int n,const int ip)
{
#ifdef DIST_MPI
	if (_num_proc>1) 
		MPI_Reduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_PROD,ip,_mpi_comm);
	else
		for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#else
	for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#endif
};
void clMpi_namespace::reduce_prod(double* const dptr_in,const int n,const int ip)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	double *dptr_out=new double [n];
	if (!dptr_out) fatal_err(0,"dptr_out in clMpi_namespace::reduce_prod");
	MPI_Reduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_PROD,ip,_mpi_comm);
	for (int i=0; i<n; i++) dptr_in[i]=dptr_out[i];
	delete [] dptr_out;
#endif
};

void clMpi_namespace::allreduce_sum(const int* const iptr_in,int* const iptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1) 
		MPI_Allreduce((void *)iptr_in,(void *)iptr_out,n,MPI_INT,MPI_SUM,_mpi_comm);
	else
		for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#else
	for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#endif
};
void clMpi_namespace::allreduce_sum(const double* const dptr_in,double* const dptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1) 
		MPI_Allreduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_SUM,_mpi_comm);
	else
		for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#else
	for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#endif
};
void clMpi_namespace::allreduce_sum(int* const iptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	int *iptr_out=new int [n];
	if (!iptr_out) fatal_err(0,"iptr_out in clMpi_namespace::allreduce_sum");
	MPI_Allreduce((void *)iptr_in,(void *)iptr_out,n,MPI_INT,MPI_SUM,_mpi_comm);
	for (int i=0; i<n; i++) iptr_in[i]=iptr_out[i];
	delete [] iptr_out;
#endif
};
void clMpi_namespace::allreduce_sum(double* const dptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	double *dptr_out=new double [n];
	if (!dptr_out) fatal_err(0,"dptr_out in clMpi_namespace::allreduce_sum");
	MPI_Allreduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_SUM,_mpi_comm);
	for (int i=0; i<n; i++) dptr_in[i]=dptr_out[i];
	delete [] dptr_out;
#endif
};

void clMpi_namespace::allreduce_min(const int* const iptr_in,int* const iptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allreduce((void *)iptr_in,(void *)iptr_out,n,MPI_INT,MPI_MIN,_mpi_comm);
	else
		for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#else
	for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#endif
};
void clMpi_namespace::allreduce_min(const double* const dptr_in,double* const dptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allreduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_MIN,_mpi_comm);
	else
		for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#else
	for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#endif
};
void clMpi_namespace::allreduce_min(int* const iptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	int *iptr_out=new int [n];
	if (!iptr_out) fatal_err(0,"iptr_out in clMpi_namespace::allreduce_min");
	MPI_Allreduce((void *)iptr_in,(void *)iptr_out,n,MPI_INT,MPI_MIN,_mpi_comm);
	for (int i=0; i<n; i++) iptr_in[i]=iptr_out[i];
	delete [] iptr_out;
#endif
};
void clMpi_namespace::allreduce_min(double* const dptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	double *dptr_out=new double [n];
	if (!dptr_out) fatal_err(0,"dptr_out in clMpi_namespace::allreduce_min");
	MPI_Allreduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_MIN,_mpi_comm);
	for (int i=0; i<n; i++) dptr_in[i]=dptr_out[i];
	delete [] dptr_out;
#endif
};

void clMpi_namespace::allreduce_max(const int* const iptr_in,int* const iptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allreduce((void *)iptr_in,(void *)iptr_out,n,MPI_INT,MPI_MAX,_mpi_comm);
	else
		for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#else
	for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#endif
};
void clMpi_namespace::allreduce_max(const double* const dptr_in,double* const dptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allreduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_MAX,_mpi_comm);
	else
		for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#else
	for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#endif
};
void clMpi_namespace::allreduce_max(int* const iptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	int *iptr_out=new int [n];
	if (!iptr_out) fatal_err(0,"iptr_out in clMpi_namespace::allreduce_max");
	MPI_Allreduce((void *)iptr_in,(void *)iptr_out,n,MPI_INT,MPI_MAX,_mpi_comm);
	for (int i=0; i<n; i++) iptr_in[i]=iptr_out[i];
	delete [] iptr_out;
#endif
};
void clMpi_namespace::allreduce_max(double* const dptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	double *dptr_out=new double [n];
	if (!dptr_out) fatal_err(0,"dptr_out in clMpi_namespace::allreduce_max");
	MPI_Allreduce((void *)dptr_in,(void *)dptr_out,n,MPI_DOUBLE,MPI_MAX,_mpi_comm);
	for (int i=0; i<n; i++) dptr_in[i]=dptr_out[i];
	delete [] dptr_out;
#endif
};

void clMpi_namespace::allreduce_maxloc(const double_int* const diptr_in,double_int* const diptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allreduce((void *)diptr_in,(void *)diptr_out,n,MPI_DOUBLE_INT,MPI_MAXLOC,_mpi_comm);
	else
		for (int i=0; i<n; i++) diptr_out[i]=diptr_in[i];
#else
	for (int i=0; i<n; i++) diptr_out[i]=diptr_in[i];
#endif
};
void clMpi_namespace::allreduce_maxloc(double_int* const diptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	double_int *diptr_out=new double_int [n];
	if (!diptr_out) fatal_err(0,"diptr_out in clMpi_namespace::allreduce_maxloc");
	MPI_Allreduce((void *)diptr_in,(void *)diptr_out,n,MPI_DOUBLE_INT,MPI_MAXLOC,_mpi_comm);
	for (int i=0; i<n; i++) diptr_in[i]=diptr_out[i];
	delete [] diptr_out;
#endif
};

void clMpi_namespace::allreduce_minloc(const double_int* const diptr_in,double_int* const diptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allreduce((void *)diptr_in,(void *)diptr_out,n,MPI_DOUBLE_INT,MPI_MINLOC,_mpi_comm);
	else
		for (int i=0; i<n; i++) diptr_out[i]=diptr_in[i];
#else
	for (int i=0; i<n; i++) diptr_out[i]=diptr_in[i];
#endif
};
void clMpi_namespace::allreduce_minloc(double_int* const diptr_in,const int n)
{
#ifdef DIST_MPI
	if (_num_proc<2) return;
	double_int *diptr_out=new double_int [n];
	if (!diptr_out) fatal_err(0,"diptr_out in clMpi_namespace::allreduce_minloc");
	MPI_Allreduce((void *)diptr_in,(void *)diptr_out,n,MPI_DOUBLE_INT,MPI_MINLOC,_mpi_comm);
	for (int i=0; i<n; i++) diptr_in[i]=diptr_out[i];
	delete [] diptr_out;
#endif
};

void clMpi_namespace::allgather(const int* const iptr_in,int* const iptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allgather((void *)iptr_in,n,MPI_INT,(void *)iptr_out,n,MPI_INT,_mpi_comm); 
	else
		for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#else
	for (int i=0; i<n; i++) iptr_out[i]=iptr_in[i];
#endif
};
void clMpi_namespace::allgather(const double* const dptr_in,double* const dptr_out,const int n)
{
#ifdef DIST_MPI
	if (_num_proc>1)
		MPI_Allgather((void *)dptr_in,n,MPI_INT,(void *)dptr_out,n,MPI_INT,_mpi_comm); 
	else
		for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#else
	for (int i=0; i<n; i++) dptr_out[i]=dptr_in[i];
#endif
};

void clMpi_namespace::allgatherv(const int* const iptr_in,const int n_in,int* const iptr_out,
	const int* const n_out,const int* const displs)
{
#ifdef DIST_MPI
	if (_num_proc>1) {
		MPI_Allgatherv((void *)iptr_in,n_in,MPI_DOUBLE,(void *)iptr_out,(int *)n_out,(int *)displs,MPI_DOUBLE,_mpi_comm); 
	} else {
		if (n_in!=(*n_out)) fatal_err("n_in!=(*n_out) in clMpi_namespace::allgatherv");
		for (int i=0; i<n_in; i++) iptr_out[i]=iptr_in[i];
	}
#else
	if (n_in!=(*n_out)) fatal_err("n_in!=(*n_out) in clMpi_namespace::allgatherv");
	for (int i=0; i<n_in; i++) iptr_out[i]=iptr_in[i];
#endif
};
void clMpi_namespace::allgatherv(const double* const dptr_in,const int n_in,double* const dptr_out,
	const int* const n_out,const int* const displs)
{
#ifdef DIST_MPI
	if (_num_proc>1) {
		MPI_Allgatherv((void *)dptr_in,n_in,MPI_DOUBLE,(void *)dptr_out,(int *)n_out,(int *)displs,MPI_DOUBLE,_mpi_comm); 
	} else {
		if (n_in!=(*n_out)) fatal_err("n_in!=(*n_out) in clMpi_namespace::allgatherv");
		for (int i=0; i<n_in; i++) dptr_out[i]=dptr_in[i];
	}
#else
	if (n_in!=(*n_out)) fatal_err("n_in!=(*n_out) in clMpi_namespace::allgatherv");
	for (int i=0; i<n_in; i++) dptr_out[i]=dptr_in[i];
#endif
};

void clMpi_namespace::transfer(const int* const iptr_send,int* const iptr_recv,const int n,
	const int ip_send,const int ip_recv,const int itag)
{
	if (_num_proc<2 && ip_send!=ip_recv) {
//		fprintf(stdout,"1: _num_proc=%d _proc_id=%d ip_send=%d ip_recv=%d itag=%d\n",_num_proc,_proc_id,ip_send,ip_recv,itag);
//		warning("1: num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
		fatal_err("num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
	}
	if (_num_proc<2 || ip_send==ip_recv) {
		if (ip_send==_proc_id && iptr_send!=iptr_recv) 
			for (int i=0; i<n; i++) iptr_recv[i]=iptr_send[i];
	} 
#ifdef DIST_MPI
	else if (_proc_id==ip_send || _proc_id==ip_recv) {
		MPI_Status status;
		if (_proc_id==ip_send) {
			MPI_Send((void *)iptr_send,n,MPI_INT,ip_recv,itag,_mpi_comm);
		} else if (_proc_id==ip_recv) {
			MPI_Recv((void *)iptr_recv,n,MPI_INT,ip_send,itag,_mpi_comm,&status);
		}
	}
#endif
};
void clMpi_namespace::transfer(const double* const dptr_send,double* const dptr_recv,const int n,
	const int ip_send,const int ip_recv,const int itag)
{
	if (_num_proc<2 && ip_send!=ip_recv) {
//		fprintf(stdout,"2: _num_proc=%d _proc_id=%d ip_send=%d ip_recv=%d itag=%d\n",_num_proc,_proc_id,ip_send,ip_recv,itag);
//		warning("2: num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
		fatal_err("num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
	}
	if (_num_proc<2 || ip_send==ip_recv) {
		if (ip_send==_proc_id && dptr_send!=dptr_recv) 
			for (int i=0; i<n; i++) dptr_recv[i]=dptr_send[i];
	}
#ifdef DIST_MPI
	else if (_proc_id==ip_send || _proc_id==ip_recv) {
		MPI_Status status;
		if (_proc_id==ip_send) {
			MPI_Send((void *)dptr_send,n,MPI_DOUBLE,ip_recv,itag,_mpi_comm);
		} else if (_proc_id==ip_recv) {
			MPI_Recv((void *)dptr_recv,n,MPI_DOUBLE,ip_send,itag,_mpi_comm,&status);
		}
	}
#endif
};
void clMpi_namespace::transfer(int* const iptr_send,const int n,const int ip_send,const int ip_recv,const int itag)
{
#ifdef DIST_MPI
	if (ip_send!=ip_recv) {
		if (_num_proc<2) {
//			fprintf(stdout,"3: _num_proc=%d _proc_id=%d ip_send=%d ip_recv=%d itag=%d\n",_num_proc,_proc_id,ip_send,ip_recv,itag);
//			warning("3: num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
			fatal_err("num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
		}
		int *iptr_recv=NULL;
		if (_proc_id==ip_recv) iptr_recv=new int [n];
		if (_proc_id==ip_send || _proc_id==ip_recv) {
			MPI_Status status;
			if (_proc_id==ip_send) {
				MPI_Send((void *)iptr_send,n,MPI_INT,ip_recv,itag,_mpi_comm);
			} else if (_proc_id==ip_recv) {
				MPI_Recv((void *)iptr_recv,n,MPI_INT,ip_send,itag,_mpi_comm,&status);
			}
		}
		if (_proc_id==ip_recv) {
			for (int i=0; i<n; i++) iptr_send[i]=iptr_recv[i];
			delete [] iptr_recv;
		}
	}
#endif
};
void clMpi_namespace::transfer(double* const dptr_send,const int n,const int ip_send,const int ip_recv,const int itag)
{
#ifdef DIST_MPI
	if (ip_send!=ip_recv) {
		if (_num_proc<2) {
//			fprintf(stdout,"4: _num_proc=%d _proc_id=%d ip_send=%d ip_recv=%d itag=%d\n",_num_proc,_proc_id,ip_send,ip_recv,itag);
//			warning("4: num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
			fatal_err("num_proc<2 && ip_send!=ip_recv in clMpi_namespace::transfer");
		}
		double *dptr_recv=NULL;
		if (_proc_id==ip_recv) dptr_recv=new double [n];
		if (_proc_id==ip_send || _proc_id==ip_recv) {
			MPI_Status status;
			if (_proc_id==ip_send) {
				MPI_Send((void *)dptr_send,n,MPI_DOUBLE,ip_recv,itag,_mpi_comm);
			} else if (_proc_id==ip_recv) {
				MPI_Recv((void *)dptr_recv,n,MPI_DOUBLE,ip_send,itag,_mpi_comm,&status);
			}
		}
		if (_proc_id==ip_recv) {
			for (int i=0; i<n; i++) dptr_send[i]=dptr_recv[i];
			delete [] dptr_recv;
		}
	}
#endif
};

/*  Print a warning message */
void clMpi_namespace::warning(const char* const warning_msg)
{
	fprintf(stderr,"WARNING: %s.\n",warning_msg);
	_n_warn++;
	if (_n_warn>=MAX_WARN) fatal_err("maximum number of warnings exceeded");
};
void clMpi_namespace::warning(const char* const warning_msg1,const char* const warning_msg2)
{
	fprintf(stderr,"WARNING: %s %s.\n",warning_msg1,warning_msg2);
	_n_warn++;
	if (_n_warn>=MAX_WARN) fatal_err("maximum number of warnings exceeded");
};

/*  Fatal error handler  */
void clMpi_namespace::fatal_err_soft(const char* const error_msg)
{
	fprintf(stderr,"\nFATAL ERROR: %s\n",error_msg);
        fflush(stderr);
        fflush(stdout);
	mpi_finalize();
	exit(1);
};
void clMpi_namespace::fatal_err(const char* const error_msg)
{
	fprintf(stderr,"\nFATAL ERROR: %s\n",error_msg);
        fflush(stderr);
        fflush(stdout);
	mpi_abort(1);
	mpi_finalize();
	exit(1);
};
void clMpi_namespace::fatal_err(const char* const error_msg,const int ierror)
{
	fprintf(stderr,"\nFATAL ERROR: %s ierror=%d\n",error_msg,ierror);
        fflush(stderr);
        fflush(stdout);
	mpi_abort(1);
	mpi_finalize();
	exit(1);
};
void clMpi_namespace::fatal_err(const int flag,const char* const error_msg)
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
	mpi_abort(1);
	mpi_finalize();
	exit(1);
};

double clMpi_namespace::walltime(void)
{
	double t=0.0;
#ifdef DIST_MPI
	t=MPI_Wtime();
#else
	fatal_err("clMpi_namespace::walltime meaningless without MPI");
#endif
	return t;
};

