/*******************************************************************************
  Author: Rolf Verberg
*******************************************************************************/

#if !defined(CLMPI_H__INCLUDED_)
#define CLMPI_H__INCLUDED_

//#ifndef DIST_MPI
//#define DIST_MPI  1
//#endif

#ifdef DIST_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef FALSE
#define FALSE	0
#endif
#ifndef TRUE
#define TRUE	1
#endif
#ifndef NULL
#define NULL	0
#endif

namespace clMpi_namespace
{
#ifdef DIST_MPI
	extern MPI_Comm _mpi_comm;        // MPI communicator (default MPI_COMM_WORLD)
#endif
	extern int _init_flag;		  // Internal flag specifying weather _num_proc and _proc_id are initialized
	extern int _num_proc;		  // Number of processors
	extern int _proc_id;		  // Processor ID or rank
	extern int _n_warn;
//	int _itag;

	const int MAX_WARN=10;
	const int BUFFERSIZE=1000;

	struct double_int {
		double value;
		int rank;
	};

	void mpi_init(int *argc,char ***argv);
#ifdef DIST_MPI
	void mpi_init(const int,const int,const MPI_Comm=MPI_COMM_WORLD);
	MPI_Comm get_mpi_comm(void);
	MPI_Comm mpi_create_comm(const int);
	void mpi_clear(void);
	void mpi_free_comm(void);
#endif
	void mpi_finalize();
	void mpi_abort(int);
	int get_num_proc(void);
	int get_proc_id(void);
	void barrier(void);
	void bcast(int* const,const int n=1,const int ip=0);
	void bcast(double* const,const int n=1,const int ip=0);
	void bcast(char* const,const int n=1,const int ip=0);
//	void bcast(bool* const,const int,const int);
	void reduce_sum(const double* const,double* const,const int n=1,const int ip=0);
	void reduce_sum(double* const,const int n=1,const int ip=0);
	void reduce_prod(const double* const,double* const,const int n=1,const int ip=0);
	void reduce_prod(double* const,const int n=1,const int ip=0);
	void allreduce_sum(const int* const,int* const,const int n=1);
	void allreduce_sum(const double* const,double* const,const int n=1);
	void allreduce_sum(int* const,const int n=1);
	void allreduce_sum(double* const,const int n=1);
	void allreduce_min(const int* const,int* const,const int n=1);
	void allreduce_min(const double* const,double* const,const int n=1);
	void allreduce_min(int* const,const int n=1);
	void allreduce_min(double* const,const int n=1);
	void allreduce_max(const int* const,int* const,const int n=1);
	void allreduce_max(const double* const,double* const,const int n=1);
	void allreduce_max(int* const,const int n=1);
	void allreduce_max(double* const,const int n=1);
	void allreduce_maxloc(const double_int* const,double_int* const,const int n=1);
	void allreduce_maxloc(double_int* const,const int n=1);
	void allreduce_minloc(const double_int* const,double_int* const,const int n=1);
	void allreduce_minloc(double_int* const,const int n=1);
	void allgather(const int* const,int* const,const int);
	void allgather(const double* const,double* const,const int);
	void allgatherv(const int* const,const int,int* const,const int* const,const int* const);
	void allgatherv(const double* const,const int,double* const,const int* const,const int* const);
	void transfer(const int* const,int* const,const int,const int,const int,const int);
	void transfer(const double* const,double* const,const int,const int,const int,const int);
	void transfer(int* const,const int,const int,const int,const int);
	void transfer(double* const,const int,const int,const int,const int);
	void warning(const char* const);
	void warning(const char* const,const char* const);
	void fatal_err_soft(const char* const);
	void fatal_err(const char* const);
	void fatal_err(const char* const,const int);
	void fatal_err(const int,const char* const);
	double walltime(void);

}; // namespace clMpi_namespace

#endif // CLMPI_H__INCLUDED_
