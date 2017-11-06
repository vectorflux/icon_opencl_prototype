//==================================================================================
//
//	This file contains the C/C++ methods called from fortran
//		to access the OpenCL operators
//	The initialization and requests for memory objects should be done when starting the program
//		are static EngineCL and ResourcesCL sufficient?
//		-> another static class is needed to contain the kernels and the memory objects
//	For now, performance outside of the kernel is not important,
//		correctness has priority
//
//==================================================================================

#include "OpenCL_statics.h"
#include <sys/time.h>
#include "math.h"
#include <stdio.h>
#include <assert.h>
#include "Profiler.h"

#include <iostream>
#include <fstream>

static OpenCL_statics ocl;
static Profiler profiler;
static Profiler profilerKernels;
static bool bProfiling = true;

using namespace std;

extern "C" void print_profiling_info_()
{
	profilerKernels.printSummary();
	profiler.printSummary();
}

extern "C" void start_cpp_timer_()
{
	cout << "starting CPP timer...\n";
	profiler.push_start("Fortran");
}

extern "C" void stop_cpp_timer_()
{
	profiler.pop_stop();
	cout << "CPP timer stopped\n";
}

#pragma mark Test Utilities
void checkresultsCPP_(double * resOCL, double * resFortran,
					  int * sx, int * ex, int * nx,
					  int * sy, int * ey, int * ny,
					  int * sz, int * ez, int * nz,
					  char * name)
{
	for (int iz=*sz; iz<*ez; iz++)
		for (int iy=*sy; iy<*ey; iy++)
			for (int ix=*sx; ix<*ex; ix++)
			{
				const int i = ix + iy*(*nx) + iz*(*nx)*(*ny);
				if (fabs(resOCL[i]-resFortran[i])>1.e-8 || isnan(fabs(resOCL[i])) || isnan(fabs(resFortran[i])))
				{
					cout << "ERROR (" << name << "): elements " << ix << "," << iy << "," << iz << "(" << i << ")" << " do not correspond (" << resOCL[i] << " vs " << resFortran[i] << ")" << endl;
					cout << "\tBounds are: [" << *sx << "," << *ex << "[ [" << *sy << "," << *ey << "[ [" << *sz << "," << *ez << "[\n";
					cout << "\tSizes are: " << *nx << " " << *ny << " " << *nz << endl;
					abort();
				}
			}
}

void checkresultsCPP_(int * resOCL, int * resFortran,
					  int * sx, int * ex, int * nx,
					  int * sy, int * ey, int * ny,
					  int * sz, int * ez, int * nz,
					  char * name)
{
	for (int iz=*sz; iz<*ez; iz++)
		for (int iy=*sy; iy<*ey; iy++)
			for (int ix=*sx; ix<*ex; ix++)
			{
				const int i = ix + iy*(*nx) + iz*(*nx)*(*ny);
				if (abs(resOCL[i]-resFortran[i])>0 || isnan(abs(resOCL[i])) || isnan(abs(resFortran[i])))
				{
					cout << "ERROR (" << name << "): elements " << ix << "," << iy << "," << iz << "(" << i << ")" << " do not correspond (" << resOCL[i] << " vs " << resFortran[i] << ")" << endl;
					cout << "\tBounds are: [" << *sx << "," << *ex << "[ [" << *sy << "," << *ey << "[ [" << *sz << "," << *ez << "[\n";
					cout << "\tSizes are: " << *nx << " " << *ny << " " << *nz << endl;
					abort();
				}
			}
}

extern "C" void checkresults_(double * resFortran, double * resOCL,
							  int * sx, int * ex, int * nx,
							  int * sy, int * ey, int * ny,
							  int * sz, int * ez, int * nz,
							  int * idx)
{
	/*
	ofstream fileOCL, fileFort;
	string nameOCL = ocl.lutCheck[*idx]+"_OCL.txt";
	fileOCL.open(nameOCL.c_str());
	fileOCL << resOCL << endl << endl;
	for (int iz=*sz; iz<*ez; iz++)
	{
		for (int iy=*sy; iy<*ey; iy++)
		{
			for (int ix=*sx; ix<*ex; ix++)
			{
				const int i = ix + iy*(*nx) + iz*(*nx)*(*ny);
				fileOCL << i << " " << resOCL[i] << " ";
			}
			fileOCL << endl;
		}
		fileOCL << endl;
	}
	fileOCL.close();
	 
	string nameFort = ocl.lutCheck[*idx]+"_Fort.txt";
	fileFort.open(nameFort.c_str());
	fileFort << resFortran << endl << endl;
	for (int iz=*sz; iz<*ez; iz++)
	{
		for (int iy=*sy; iy<*ey; iy++)
		{
			for (int ix=*sx; ix<*ex; ix++)
			{
				const int i = ix + iy*(*nx) + iz*(*nx)*(*ny);
				fileFort << i << " " << resFortran[i] << " ";
			}
			fileFort << endl;
		}
		fileFort << endl;
	}
	fileFort.close();
	//*/
	
	for (int iz=*sz; iz<*ez; iz++)
		for (int iy=*sy; iy<*ey; iy++)
			for (int ix=*sx; ix<*ex; ix++)
			{
				const int i = ix + iy*(*nx) + iz*(*nx)*(*ny);
				if (fabs(resOCL[i]-resFortran[i])>1.e-6 || isnan(fabs(resOCL[i])) || isnan(fabs(resFortran[i])))
				{
					cout << "ERROR (" << ocl.lutCheck[*idx] << "): elements " << ix << "," << iy << "," << iz << "(" << i << ")" << " do not correspond (" << resFortran[i] << " vs " << resOCL[i] << ")" << endl;
					cout << "\tBounds are: [" << *sx << "," << *ex << "[ [" << *sy << "," << *ey << "[ [" << *sz << "," << *ez << "[\n";
					cout << "\tSizes are: " << *nx << " " << *ny << " " << *nz << endl;
					abort();
				}
			}
}

extern "C" void checkresults4dCPP_(double * resOCL, double * resFortran,
							    int * sx, int * ex, int * nx,
								int * sy, int * ey, int * ny,
							    int * sz, int * ez, int * nz,
							    int * sw, int * ew, int * nw,
							    char * name="")
{
	/*
	 ofstream fileOCL, fileFort;
	 string nameOCL = "4d_OCL.txt";
	 fileOCL.open(nameOCL.c_str());
	 fileOCL << resOCL << endl << endl;
	 for (int iw=*sw; iw<*ew; iw++)
	 {
	 for (int iz=*sz; iz<*ez; iz++)
	 {
	 for (int iy=*sy; iy<*ey; iy++)
	 {
	 for (int ix=*sx; ix<*ex; ix++)
	 {
	 const int i = ix + iy*(*nx) + iz*(*nx)*(*ny) + iw*(*nx)*(*ny)*(*nz);
	 fileOCL << i << " " << resOCL[i] << " ";
	 }
	 fileOCL << endl;
	 }
	 fileOCL << endl;
	 }
	 fileOCL << endl;
	 }
	 fileOCL.close();
	 
	 string nameFort = "4d_Fort.txt";
	 fileFort.open(nameFort.c_str());
	 fileFort << resFortran << endl << endl;
	 for (int iw=*sw; iw<*ew; iw++)
	 {
	 for (int iz=*sz; iz<*ez; iz++)
	 {
	 for (int iy=*sy; iy<*ey; iy++)
	 {
	 for (int ix=*sx; ix<*ex; ix++)
	 {
	 const int i = ix + iy*(*nx) + iz*(*nx)*(*ny) + iw*(*nx)*(*ny)*(*nz);
	 fileFort << i << " " << resFortran[i] << " ";
	 }
	 fileFort << endl;
	 }
	 fileFort << endl;
	 }
	 fileOCL << endl;
	 }
	 fileFort.close();
	 //*/
	
	for (int iw=*sw; iw<*ew; iw++)
		for (int iz=*sz; iz<*ez; iz++)
			for (int iy=*sy; iy<*ey; iy++)
				for (int ix=*sx; ix<*ex; ix++)
				{
					const int i = ix + iy*(*nx) + iz*(*nx)*(*ny) + iw*(*nx)*(*ny)*(*nz);
					//if (iz==0) cout << "ERROR (" << name << "): elements " << ix << "," << iy << "," << iz << "(" << i << ")" << " do not correspond (" << divOCL[i] << " vs " << divFortran[i] << ")" << endl;
					if (fabs(resOCL[i]-resFortran[i])>1.e-6)
					{
						cout << "ERROR (" << name << "): elements " << ix << "," << iy << "," << iz << "," << iw << "(" << i << ")" << " do not correspond (" << resOCL[i] << " vs " << resFortran[i] << ")" << endl;
						cout << "\tBounds are: [" << *sx << "," << *ex << "[ [" << *sy << "," << *ey << "[ [" << *sz << "," << *ez << "[ [" << *sw << "," << *ew << "[\n";
						cout << "\tSizes are: " << *nx << " " << *ny << " " << *nz << " " << *nw << endl;
						abort();
					}
					//cout << fabs(divOCL[i]-divFortran[i]) << endl;
				}
	//abort();
}

extern "C" void checkresults4d_(double * resOCL, double * resFortran,
							    int * sx, int * ex, int * nx,
								int * sy, int * ey, int * ny,
							    int * sz, int * ez, int * nz,
							    int * sw, int * ew, int * nw,
							    int * idx)
{
	for (int iw=*sw; iw<*ew; iw++)
		for (int iz=*sz; iz<*ez; iz++)
			for (int iy=*sy; iy<*ey; iy++)
				for (int ix=*sx; ix<*ex; ix++)
				{
					const int i = ix + iy*(*nx) + iz*(*nx)*(*ny) + iw*(*nx)*(*ny)*(*nz);
					//if (iz==0) cout << "ERROR (" << name << "): elements " << ix << "," << iy << "," << iz << "(" << i << ")" << " do not correspond (" << divOCL[i] << " vs " << divFortran[i] << ")" << endl;
					if (fabs(resOCL[i]-resFortran[i])>1.e-6)
					{
						cout << "ERROR (" << ocl.lutCheck[*idx] << "): elements " << ix << "," << iy << "," << iz << "," << iw << "(" << i << ")" << " do not correspond (" << resOCL[i] << " vs " << resFortran[i] << ")" << endl;
						cout << "\tBounds are: [" << *sx << "," << *ex << "[ [" << *sy << "," << *ey << "[ [" << *sz << "," << *ez << "[ [" << *sw << "," << *ew << "[\n";
						cout << "\tSizes are: " << *nx << " " << *ny << " " << *nz << " " << *nw << endl;
						abort();
					}
					//cout << fabs(divOCL[i]-divFortran[i]) << endl;
				}
	//abort();
}

extern "C" void checkresults4d_int_(int * divOCL, int * divFortran,
							    int * sx, int * ex, int * nx,
								int * sy, int * ey, int * ny,
							    int * sz, int * ez, int * nz,
							    int * sw, int * ew, int * nw,
							    char * name="")
{
	for (int iw=*sw; iw<*ew; iw++)
		for (int iz=*sz; iz<*ez; iz++)
			for (int iy=*sy; iy<*ey; iy++)
				for (int ix=*sx; ix<*ex; ix++)
				{
					const int i = ix + iy*(*nx) + iz*(*nx)*(*ny) + iw*(*nx)*(*ny)*(*nz);
					//if (iz==0) cout << "ERROR (" << name << "): elements " << ix << "," << iy << "," << iz << "(" << i << ")" << " do not correspond (" << divOCL[i] << " vs " << divFortran[i] << ")" << endl;
					if (abs(divOCL[i]-divFortran[i])>0)
					{
						cout << "ERROR (" << name << "): elements " << ix << "," << iy << "," << iz << "," << iw << "(" << i << ")" << " do not correspond (" << divOCL[i] << " vs " << divFortran[i] << ")" << endl;
						cout << "\tBounds are: [" << *sx << "," << *ex << "[ [" << *sy << "," << *ey << "[ [" << *sz << "," << *ez << "[ [" << *sw << "," << *ew << "[\n";
						cout << "\tSizes are: " << *nx << " " << *ny << " " << *nz << " " << *nw << endl;
						abort();
					}
					//cout << fabs(divOCL[i]-divFortran[i]) << endl;
				}
	//abort();
}

#pragma mark Initialization
extern "C" void ocl_init_(int * nblks_c,
						  int * nlev, int * nlevp1,
						  int * nproma,
						  double * wgtfac_c)
{
	// send all constant structures to the GPU with this method
	if (bProfiling) profiler.push_start("ocl_init");
	
	cl_event ev_write_wgtfac_c, ev_kernel_init;
	
	const int size_wgtfac_c     = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_wgtfac_c_out = (*nproma)*(*nlev  )*(*nblks_c);

	const size_t global_init[3] = { *nblks_c, *nlev, *nproma };

	cl_mem wgtfac_c_d      = resourcesCL.getBuffers("wgtfac_c",      1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_wgtfac_c)[0];
	cl_mem wgtfac_c_out1_d = resourcesCL.getBuffers("wgtfac_c_out1", 1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_wgtfac_c_out)[0];
	cl_mem wgtfac_c_out2_d = resourcesCL.getBuffers("wgtfac_c_out2", 1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_wgtfac_c_out)[0];
	
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfac_c_d, CL_FALSE, 0, sizeof(double)*size_wgtfac_c, wgtfac_c, 0, NULL, &ev_write_wgtfac_c) );
	
	CheckCL( clSetKernelArg(ocl.kernel_init, 0, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_init, 1, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_init, 2, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_init, 3, sizeof(cl_mem), (void*)&wgtfac_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_init, 4, sizeof(cl_mem), (void*)&wgtfac_c_out1_d) );
	CheckCL( clSetKernelArg(ocl.kernel_init, 5, sizeof(cl_mem), (void*)&wgtfac_c_out2_d) );
	
	CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_init, 3, NULL, global_init, NULL, 1, &ev_write_wgtfac_c, &ev_kernel_init) );
	
	CheckCL(clFlush(oclhelper.getCommandQueue()));
	CheckCL(clFinish(oclhelper.getCommandQueue()));
	
#ifdef TEST
	cout << "checking init...\n";
	
	double * wgtfac_c_out1 = new double[size_wgtfac_c_out];
	double * wgtfac_c_out2 = new double[size_wgtfac_c_out];
	
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfac_c_out1_d, CL_TRUE, 0, sizeof(cl_double)*size_wgtfac_c_out, wgtfac_c_out1, 0, NULL, NULL) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfac_c_out2_d, CL_TRUE, 0, sizeof(cl_double)*size_wgtfac_c_out, wgtfac_c_out2, 0, NULL, NULL) );

	for (int ib=0; ib<*nblks_c; ib++)
		for (int ik=0; ik<*nlev; ik++)
			for (int ic=0; ic<*nproma; ic++)
			{
				const int idx_in  = ic + ik*(*nproma) + ib*(*nproma)*(*nlevp1);
				const int idx_out = ic + ik*(*nproma) + ib*(*nproma)*(*nlev);
				
				assert(wgtfac_c_out1[idx_out] == wgtfac_c[idx_in]);
				assert(wgtfac_c_out2[idx_out] == wgtfac_c[idx_in+(*nproma)]);
			}
	
	//cout << "test passed!\n";
	//exit(0);
#endif
	
	if (bProfiling) profiler.pop_stop();
}

#pragma mark Solve NH
extern "C" void ocl_exner_(int * i_startblk, int * i_endblk, int * nblks_c,
						   int * nflat, int * nlev, int * nlevp1,
						   int * localStart, int * localEnd, int * nproma,
						   int * istep,
						   double * exner_ref_mc,
						   double * exner_exfac,
						   double * exner_nnow,
						   double * exner_old,
						   double * z_exner_ex_pr,
						   double * z_exner_pr,
						   double * wgtfacq1_c,
						   double * wgtfacq_c,
						   double * wgtfac_c,
						   double * inv_ddqz_z_full,
						   double * z_dexner_dz_c)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_exner");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int i_endblkm1   = *i_endblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_exner_ref_mc, ev_write_exner_exfac, ev_write_exner_nnow, ev_write_exner_old;
	cl_event ev_write_wgtfacq1_c, ev_write_wgtfacq_c;
	cl_event ev_write_wgtfac_c, ev_write_inv_ddqz_z_full;
	cl_event ev_kernel_exner1, ev_kernel_exner2, ev_kernel_exner3a, ev_kernel_exner3b;
	cl_event ev_read_exner_old;
	
	const size_t local_exner_3a[3] = { 1, *nlev, 1 };
	const size_t local_exner_3b[3] = { 1, *nlevp1, 1 };

	const size_t global_exner1[3] = { nBlocks, *nlev, *nproma };
	const size_t global_exner2[2] = { nBlocks, *nproma };
	const size_t global_exner3a[3] = { nBlocks, *nlev, *nproma };
	const size_t global_exner3b[3] = { nBlocks, *nlevp1, *nproma };
	
	const int size_localStart      = nBlocks;
	const int size_localEnd        = nBlocks;
	const int size_3d = (*nproma)*(*nlev)*(*nblks_c);
	const int size_wgtfacq_c       = (*nproma)*3*(*nblks_c);
	const int size_z_exner_ic      = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_inv_ddqz_z_full = (*nproma)*(*nlev)*(*nblks_c);
	const int size_z_dexner_dz_c   = 2*(*nproma)*(*nlev)*(*nblks_c);
	const int size_wgtfac_c        = (*nproma)*(*nlevp1)*(*nblks_c);
	
	cl_mem localStart_d      = resourcesCL.getBuffers("localStart3",     1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d        = resourcesCL.getBuffers("localEnd3",       1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(int)*size_localEnd)[0];
	cl_mem exner_ref_mc_d    = resourcesCL.getBuffers("exner_ref_mc",    1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_3d)[0];
	cl_mem exner_exfac_d     = resourcesCL.getBuffers("exner_exfac",     1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_3d)[0];
	cl_mem exner_nnow_d      = resourcesCL.getBuffers("exner_nnow",      1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_3d)[0];
	cl_mem exner_old_d       = resourcesCL.getBuffers("exner_old",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_exner_ex_pr_d   = resourcesCL.getBuffers("z_exner_ex_pr",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_exner_pr_d      = resourcesCL.getBuffers("z_exner_pr",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem wgtfacq1_c_d      = resourcesCL.getBuffers("wgtfacq1_c",      1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_wgtfacq_c)[0];
	cl_mem wgtfacq_c_d       = resourcesCL.getBuffers("wgtfacq_c",       1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_wgtfacq_c)[0];
	cl_mem wgtfac_c_d        = resourcesCL.getBuffers("wgtfac_c",        1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_wgtfac_c)[0];
	cl_mem z_exner_ic_d      = resourcesCL.getBuffers("z_exner_ic",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_exner_ic)[0];
	cl_mem inv_ddqz_z_full_d = resourcesCL.getBuffers("inv_ddqz_z_full", 1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(double)*size_3d)[0];
	cl_mem z_dexner_dz_c_d   = resourcesCL.getBuffers("z_dexner_dz_c",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_dexner_dz_c)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), exner_ref_mc_d,    CL_FALSE, 0, sizeof(double)*size_3d,    exner_ref_mc,    0, NULL, &ev_write_exner_ref_mc) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), exner_exfac_d,     CL_FALSE, 0, sizeof(double)*size_3d,     exner_exfac,     0, NULL, &ev_write_exner_exfac) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfacq1_c_d,      CL_FALSE, 0, sizeof(double)*size_wgtfacq_c,      wgtfacq1_c,      0, NULL, &ev_write_wgtfacq1_c) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfacq_c_d,       CL_FALSE, 0, sizeof(double)*size_wgtfacq_c,       wgtfacq_c,       0, NULL, &ev_write_wgtfacq_c) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfac_c_d,        CL_FALSE, 0, sizeof(double)*size_wgtfac_c,        wgtfac_c,        0, NULL, &ev_write_wgtfac_c) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_full_d, CL_FALSE, 0, sizeof(double)*size_3d, inv_ddqz_z_full, 0, NULL, &ev_write_inv_ddqz_z_full) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,      CL_FALSE, 0, sizeof(int)*size_localStart,         localStart,      0, NULL, &ev_write_localStart) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,        CL_FALSE, 0, sizeof(int)*size_localEnd,           localEnd,        0, NULL, &ev_write_localEnd) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), exner_nnow_d,      CL_FALSE, 0, sizeof(double)*size_3d, exner_nnow,      0, NULL, &ev_write_exner_nnow) );
	}
	
	if (*istep==1)
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), exner_old_d,       CL_FALSE, 0, sizeof(double)*size_3d, exner_old,       0, NULL, &ev_write_exner_old) );
	
#if TEST
	cl_event ev_write_z_exner_pr, ev_write_z_exner_ex_pr, ev_write_z_dexner_dz_c;
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_exner_ex_pr_d,   CL_FALSE, 0, sizeof(double)*size_3d,   z_exner_ex_pr,   0, NULL, &ev_write_z_exner_ex_pr) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_exner_pr_d,      CL_FALSE, 0, sizeof(double)*size_3d,      z_exner_pr,      0, NULL, &ev_write_z_exner_pr) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_dexner_dz_c_d,   CL_FALSE, 0, sizeof(double)*size_z_dexner_dz_c,   z_dexner_dz_c,   0, NULL, &ev_write_z_dexner_dz_c) );
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	{
		cout << "checking transfers for ocl_exner...\n";
		double * exner_ref_mc_tmp    = new double[size_3d];
		double * exner_exfac_tmp     = new double[size_3d];
		double * wgtfacq1_c_tmp      = new double[size_wgtfacq_c];
		double * wgtfacq_c_tmp       = new double[size_wgtfacq_c];
		double * wgtfac_c_tmp        = new double[size_wgtfac_c];
		double * inv_ddqz_z_full_tmp = new double[size_3d];
		int zero[1] = {0};
		int three[1] = {3};
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_ref_mc_d,    CL_FALSE, 0, sizeof(cl_double)*size_3d,    exner_ref_mc_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_exfac_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d,     exner_exfac_tmp,   0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfacq1_c_d,      CL_FALSE, 0, sizeof(cl_double)*size_wgtfacq_c,      wgtfacq1_c_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfacq_c_d,       CL_FALSE, 0, sizeof(cl_double)*size_wgtfacq_c,       wgtfacq_c_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfac_c_d,        CL_FALSE, 0, sizeof(cl_double)*size_wgtfac_c,        wgtfac_c_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_full_d, CL_FALSE, 0, sizeof(cl_double)*size_3d, inv_ddqz_z_full_tmp,  0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		checkresultsCPP_(exner_ref_mc,exner_ref_mc_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"exner:exner_ref_mc");
		
		checkresultsCPP_(exner_exfac,exner_exfac_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"exner:exner_exfac");
		
		checkresultsCPP_(wgtfacq1_c,wgtfacq1_c_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_c,nblks_c,"exner:wgtfacq1_c");
		
		checkresultsCPP_(wgtfacq_c,wgtfacq_c_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_c,nblks_c,"exner:wgtfacq_c");
		
		checkresultsCPP_(wgtfac_c,wgtfac_c_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"exner:wgtfac_c");
		
		checkresultsCPP_(inv_ddqz_z_full,inv_ddqz_z_full_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"exner:inv_ddqz_z_full");
		
		delete [] exner_ref_mc_tmp;
		delete [] exner_exfac_tmp;
		delete [] wgtfacq1_c_tmp;
		delete [] wgtfacq_c_tmp;
		delete [] wgtfac_c_tmp;
		delete [] inv_ddqz_z_full_tmp;
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("exner1");
	}
#endif
	
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 4, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 5, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 6, sizeof(cl_mem), (void*)&exner_ref_mc_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 7, sizeof(cl_mem), (void*)&exner_exfac_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 8, sizeof(cl_mem), (void*)&exner_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1, 9, sizeof(cl_mem), (void*)&exner_old_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1,10, sizeof(cl_mem), (void*)&z_exner_ex_pr_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner1,11, sizeof(cl_mem), (void*)&z_exner_pr_d) );

	//cl_event ev_wait_list_exner1[8] = { ev_write_localStart, ev_write_localEnd, ev_write_exner_ref_mc, ev_write_exner_exfac, ev_write_exner_nnow, ev_write_exner_old, ev_write_z_exner_ex_pr, ev_write_z_exner_pr };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner1, 3, NULL, global_exner1, NULL, 8, ev_wait_list_exner1, &ev_kernel_exner1) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_exner1[4] = { ev_write_exner_ref_mc, ev_write_exner_exfac, ev_write_exner_nnow, ev_write_exner_old };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner1, 3, NULL, global_exner1, NULL, 4, ev_wait_list_exner1, &ev_kernel_exner1) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_exner1[2] = { ev_write_exner_nnow, ev_write_exner_old };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner1, 3, NULL, global_exner1, NULL, 2, ev_wait_list_exner1, &ev_kernel_exner1) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner1, 1, NULL, global_exner1, NULL, 0, NULL, &ev_kernel_exner1) );
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("exner2");
	}
	
	if (*nflat==1)
	{
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 0, sizeof(int), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 1, sizeof(int), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 2, sizeof(int), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 3, sizeof(int), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 4, sizeof(int), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 5, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 6, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 7, sizeof(cl_mem), (void*)&wgtfacq1_c_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 8, sizeof(cl_mem), (void*)&wgtfacq_c_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a, 9, sizeof(cl_mem), (void*)&z_exner_ex_pr_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2a,10, sizeof(cl_mem), (void*)&z_exner_ic_d) );
		
		//cl_event ev_wait_list_exner2a[7] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq1_c, ev_write_wgtfacq_c, ev_write_z_exner_ex_pr, ev_write_z_exner_pr, ev_kernel_exner1 };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner2a, 2, NULL, global_exner2, NULL, 7, ev_wait_list_exner2a, &ev_kernel_exner2) );
		if (ocl.bFirstIter)
		{
			cl_event ev_wait_list_exner2a[3] = { ev_write_wgtfacq1_c, ev_write_wgtfacq_c, ev_kernel_exner1 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner2a, 2, NULL, global_exner2, NULL, 3, ev_wait_list_exner2a, &ev_kernel_exner2) );
		}
		else
		{
			cl_event ev_wait_list_exner2a[1] = { ev_kernel_exner1 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner2a, 2, NULL, global_exner2, NULL, 1, ev_wait_list_exner2a, &ev_kernel_exner2) );
		}
	}
	else
	{
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 0, sizeof(int), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 1, sizeof(int), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 2, sizeof(int), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 3, sizeof(int), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 4, sizeof(int), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 5, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 6, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 7, sizeof(cl_mem), (void*)&wgtfacq_c_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 8, sizeof(cl_mem), (void*)&z_exner_ex_pr_d) );
		CheckCL( clSetKernelArg(ocl.kernel_exner2b, 9, sizeof(cl_mem), (void*)&z_exner_ic_d) );
		
		//cl_event ev_wait_list_exner2b[6] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq_c, ev_write_z_exner_ex_pr, ev_write_z_exner_pr, ev_kernel_exner1 };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner2b, 2, NULL, global_exner2, NULL, 6, ev_wait_list_exner2b, &ev_kernel_exner2) );
		if (ocl.bFirstIter)
		{
			cl_event ev_wait_list_exner2b[2] = { ev_write_wgtfacq_c, ev_kernel_exner1 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner2b, 2, NULL, global_exner2, NULL, 2, ev_wait_list_exner2b, &ev_kernel_exner2) );
		}
		else
		{
			cl_event ev_wait_list_exner2b[1] = { ev_kernel_exner1 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner2b, 2, NULL, global_exner2, NULL, 1, ev_wait_list_exner2b, &ev_kernel_exner2) );
		}
	}
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("exner3a");
	}
	
	const int slev = max(2,*nflat)-1;
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 0, sizeof(int), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 1, sizeof(int), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 2, sizeof(int), (void*)&slev) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 3, sizeof(int), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 4, sizeof(int), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 5, sizeof(int), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 8, sizeof(cl_mem), (void*)&wgtfac_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a, 9, sizeof(cl_mem), (void*)&z_exner_ex_pr_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3a,10, sizeof(cl_mem), (void*)&z_exner_ic_d) );
	
	//cl_event ev_wait_list_exner3a[5] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfac_c, ev_write_z_exner_ex_pr, ev_kernel_exner2 };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner3a, 3, NULL, global_exner3, NULL, 5, ev_wait_list_exner3a, &ev_kernel_exner3a) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_exner3a[2] = { ev_write_wgtfac_c, ev_kernel_exner2 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner3a, 3, NULL, global_exner3a, local_exner_3a, 2, ev_wait_list_exner3a, &ev_kernel_exner3a) );
	}
	else
	{
		cl_event ev_wait_list_exner3a[1] = { ev_kernel_exner2 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner3a, 3, NULL, global_exner3a, local_exner_3a, 1, ev_wait_list_exner3a, &ev_kernel_exner3a) );
	}
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("exner3b");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 0, sizeof(int), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 1, sizeof(int), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 2, sizeof(int), (void*)&slev) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 3, sizeof(int), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 4, sizeof(int), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 5, sizeof(int), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 8, sizeof(cl_mem), (void*)&inv_ddqz_z_full_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b, 9, sizeof(cl_mem), (void*)&z_exner_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_exner3b,10, sizeof(cl_mem), (void*)&z_dexner_dz_c_d) );
	
	//cl_event ev_wait_list_exner3b[5] = { ev_write_localStart, ev_write_localEnd, ev_write_inv_ddqz_z_full, ev_write_z_dexner_dz_c, ev_kernel_exner3a };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner3b, 3, NULL, global_exner3, NULL, 5, ev_wait_list_exner3b, &ev_kernel_exner3b) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_exner3b[2] = { ev_write_inv_ddqz_z_full, ev_kernel_exner3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner3b, 3, NULL, global_exner3b, local_exner_3b, 2, ev_wait_list_exner3b, &ev_kernel_exner3b) );
	}
	else
	{
		cl_event ev_wait_list_exner3b[1] = { ev_kernel_exner3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_exner3b, 3, NULL, global_exner3b, local_exner_3b, 1, ev_wait_list_exner3b, &ev_kernel_exner3b) );
	}
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());

	if (*istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_old_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d,     exner_old,     0, NULL, &ev_read_exner_old) );
	
#if TEST
	cl_event ev_read_z_exner_ex_pr, ev_read_z_exner_pr, ev_read_z_dexner_dz_c;
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_old_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d,     exner_old,     0, NULL, &ev_read_exner_old) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_exner_ex_pr_d, CL_FALSE, 0, sizeof(cl_double)*size_3d, z_exner_ex_pr, 0, NULL, &ev_read_z_exner_ex_pr) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_exner_pr_d,    CL_FALSE, 0, sizeof(cl_double)*size_3d, z_exner_pr,    0, NULL, &ev_read_z_exner_pr) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_dexner_dz_c_d, CL_FALSE, 0, sizeof(cl_double)*size_z_dexner_dz_c, z_dexner_dz_c, 0, NULL, &ev_read_z_dexner_dz_c) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());

	if (bProfiling)
		profiler.pop_stop();
}

extern "C" void ocl_theta_(int * i_startblk, int * i_endblk, int * nblks_c,
						   int * slev, int * nlev, int * nlevp1,
						   int * localStart, int * localEnd, int * nproma,
						   int * istep,
						   double * theta_v_nvar, int * nameIdx,
						   double * theta_v_nnow,
						   double * theta_ref_mc,
						   double * rho_refcorr_ic,
						   double * wgtfac_c,
						   double * rho_nvar, int * nameIdx2,
						   double * rho_nnow,
						   double * vwind_expl_wgt,
						   double * z_exner_pr,
						   double * ddqz_z_half,
						   double * d_exner_dz_ref_ic,
						   double * rho_ic,
						   double * z_th_ddz_exner_c,
						   double * wgtfacq_c,
						   double * theta_ref_ic,
						   double * theta_v_h,
						   double * d_exner_dz_ref_mc,
						   double * d2_exner_dz2_ref_mc,
						   double * z_dexner_dz_c)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_theta");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int i_endblkm1   = *i_endblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;
	const int slevm1 = *slev-1;

	if (ocl.bFirstIter)
		for (int i=0; i<nBlocks; i++)
		{
			if (localStart[i]!=0 || localEnd[i]!=1)
			{
				cout << "need to use 3d indices\n";
				abort();
			}
		}
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_theta_v_nvar, ev_write_theta_v_nnow, ev_write_theta_ref_mc;
	cl_event ev_write_rho_refcorr_ic/*, ev_write_wgtfac_c*/, ev_write_rho_nvar, ev_write_rho_nnow, ev_write_theta_ref_ic, ev_write_vwind_expl_wgt/*, ev_write_z_exner_pr*/, ev_write_ddqz_z_half, ev_write_d_exner_dz_ref_ic;
	cl_event ev_write_rho_ic, ev_write_theta_v_h/*, ev_write_z_th_ddz_exner_c*/;
	cl_event /*ev_write_wgtfacq_c, */ev_write_d_exner_dz_ref_mc, ev_write_d2_exner_dz2_ref_mc/*, ev_write_z_dexner_dz_c*/;
	cl_event ev_kernel_theta1, ev_kernel_theta2, ev_kernel_theta3a, ev_kernel_theta3b;
	cl_event ev_read_rho_ic, ev_read_theta_v_h/*, ev_read_z_th_ddz_exner_c, ev_read_z_dexner_dz_c*/;
	
	const size_t local_size[3] = { 1, *nlev, 1 };
	
	const size_t global_theta1[1]  = { nBlocks*(*nlev)*(*nproma) };
	const size_t global_theta2[3]  = { nBlocks, *nlev, *nproma };
	const size_t global_theta3a[2] = { nBlocks, *nproma };
	const size_t global_theta3b[3] = { nBlocks, *nlev-1, *nproma };
	
	const int size_localStart     = nBlocks;
	const int size_localEnd       = nBlocks;
	const int size_3d_nlev        = (*nproma)*(*nlev)*(*nblks_c);
	const int size_3d_nlevp1      = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_wgtfacq_c      = (*nproma)*3*(*nblks_c);
	const int size_z_dexner_dz_c  = 2*(*nproma)*(*nlev)*(*nblks_c);
	const int size_vwind_expl_wgt = (*nproma)*(*nblks_c);
	
	cl_mem theta_v_nvar_d        = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx].c_str(),  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem rho_nvar_d            = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx2].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem localStart_d          = resourcesCL.getBuffers("localStartTheta",     1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d            = resourcesCL.getBuffers("localEndTheta",       1, oclhelper.getContext(), CL_MEM_READ_WRITE,  sizeof(int)*size_localEnd)[0];
	cl_mem theta_v_nnow_d        = resourcesCL.getBuffers("theta_v_nnow",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem theta_ref_mc_d        = resourcesCL.getBuffers("theta_ref_mc",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem z_theta_v_pr_mc_d     = resourcesCL.getBuffers("z_theta_v_pr_mc",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem rho_refcorr_ic_d      = resourcesCL.getBuffers("rho_refcorr_ic",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem wgtfac_c_d            = resourcesCL.getBuffers("wgtfac_c",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem rho_nnow_d            = resourcesCL.getBuffers("rho_nnow",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem vwind_expl_wgt_d      = resourcesCL.getBuffers("vwind_expl_wgt",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vwind_expl_wgt)[0];
	cl_mem z_exner_pr_d          = resourcesCL.getBuffers("z_exner_pr",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem ddqz_z_half_d         = resourcesCL.getBuffers("ddqz_z_half",         1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem d_exner_dz_ref_ic_d   = resourcesCL.getBuffers("d_exner_dz_ref_ic",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem rho_ic_d              = resourcesCL.getBuffers("rho_ic",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem z_theta_v_pr_ic_d     = resourcesCL.getBuffers("z_theta_v_pr_ic",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem z_th_ddz_exner_c_d    = resourcesCL.getBuffers("z_th_ddz_exner_c",    1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem wgtfacq_c_d           = resourcesCL.getBuffers("wgtfacq_c",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_wgtfacq_c)[0];
	cl_mem theta_ref_ic_d        = resourcesCL.getBuffers("theta_ref_ic",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem theta_v_h_d           = resourcesCL.getBuffers("theta_v_h",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlevp1)[0];
	cl_mem d_exner_dz_ref_mc_d   = resourcesCL.getBuffers("d_exner_dz_ref_mc",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem d2_exner_dz2_ref_mc_d = resourcesCL.getBuffers("d2_exner_dz2_ref_mc", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_nlev)[0];
	cl_mem z_dexner_dz_c_d       = resourcesCL.getBuffers("z_dexner_dz_c",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_dexner_dz_c)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_ref_mc_d,        CL_FALSE, 0, sizeof(double)*size_3d_nlev,        theta_ref_mc,        0, NULL, &ev_write_theta_ref_mc) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rho_refcorr_ic_d,      CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      rho_refcorr_ic,      0, NULL, &ev_write_rho_refcorr_ic) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfac_c_d,            CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      wgtfac_c,            0, NULL, &ev_write_wgtfac_c) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vwind_expl_wgt_d,      CL_FALSE, 0, sizeof(double)*size_vwind_expl_wgt, vwind_expl_wgt,      0, NULL, &ev_write_vwind_expl_wgt) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddqz_z_half_d,         CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      ddqz_z_half,         0, NULL, &ev_write_ddqz_z_half) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), d_exner_dz_ref_ic_d,   CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      d_exner_dz_ref_ic,   0, NULL, &ev_write_d_exner_dz_ref_ic) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfacq_c_d,           CL_FALSE, 0, sizeof(double)*size_wgtfacq_c,      wgtfacq_c,           0, NULL, &ev_write_wgtfacq_c) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_ref_ic_d,        CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      theta_ref_ic,        0, NULL, &ev_write_theta_ref_ic) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), d_exner_dz_ref_mc_d,   CL_FALSE, 0, sizeof(double)*size_3d_nlev,        d_exner_dz_ref_mc,   0, NULL, &ev_write_d_exner_dz_ref_mc) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), d2_exner_dz2_ref_mc_d, CL_FALSE, 0, sizeof(double)*size_3d_nlev,        d2_exner_dz2_ref_mc, 0, NULL, &ev_write_d2_exner_dz2_ref_mc) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,          CL_FALSE, 0, sizeof(int)*size_localStart,        localStart,          0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,            CL_FALSE, 0, sizeof(int)*size_localEnd,          localEnd,            0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_v_nvar_d,        CL_FALSE, 0, sizeof(double)*size_3d_nlev,        theta_v_nvar,        0, NULL, &ev_write_theta_v_nvar) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_v_nnow_d,        CL_FALSE, 0, sizeof(double)*size_3d_nlev,        theta_v_nnow,        0, NULL, &ev_write_theta_v_nnow) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rho_nvar_d,            CL_FALSE, 0, sizeof(double)*size_3d_nlev,        rho_nvar,            0, NULL, &ev_write_rho_nvar) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rho_nnow_d,            CL_FALSE, 0, sizeof(double)*size_3d_nlev,        rho_nnow,            0, NULL, &ev_write_rho_nnow) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rho_ic_d,              CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      rho_ic,              0, NULL, &ev_write_rho_ic) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_v_h_d,           CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      theta_v_h,           0, NULL, &ev_write_theta_v_h) );
	}
	
#ifdef TEST
	cl_event ev_write_z_th_ddz_exner_c, ev_write_z_dexner_dz_c, ev_write_z_exner_pr;
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_exner_pr_d,          CL_FALSE, 0, sizeof(double)*size_3d_nlev,        z_exner_pr,          0, NULL, &ev_write_z_exner_pr) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_th_ddz_exner_c_d,    CL_FALSE, 0, sizeof(double)*size_3d_nlevp1,      z_th_ddz_exner_c,    0, NULL, &ev_write_z_th_ddz_exner_c) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_dexner_dz_c_d,       CL_FALSE, 0, sizeof(double)*size_z_dexner_dz_c,  z_dexner_dz_c,       0, NULL, &ev_write_z_dexner_dz_c) );

	{
		cout << "checking transfers for ocl_theta...\n";
		double * theta_ref_mc_tmp        = new double[size_3d_nlev];
		double * rho_refcorr_ic_tmp      = new double[size_3d_nlevp1];
		double * vwind_expl_wgt_tmp      = new double[size_vwind_expl_wgt];
		double * ddqz_z_half_tmp         = new double[size_3d_nlevp1];
		double * d_exner_dz_ref_ic_tmp   = new double[size_3d_nlevp1];
		double * theta_ref_ic_tmp        = new double[size_3d_nlevp1];
		double * d_exner_dz_ref_mc_tmp   = new double[size_3d_nlev];
		double * d2_exner_dz2_ref_mc_tmp = new double[size_3d_nlev];
		double * wgtfac_c_tmp            = new double[size_3d_nlevp1];
		double * wgtfacq_c_tmp           = new double[size_wgtfacq_c];
		double * theta_v_nnow_tmp        = new double[size_3d_nlev];
		double * rho_nnow_tmp            = new double[size_3d_nlev];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		int zero[1] = {0};
		int one[1] = {1};
		int three[1] = {3};
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_ref_mc_d,        CL_FALSE, 0, sizeof(cl_double)*size_3d_nlev,        theta_ref_mc_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_refcorr_ic_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,      rho_refcorr_ic_tmp,   0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vwind_expl_wgt_d,      CL_FALSE, 0, sizeof(cl_double)*size_vwind_expl_wgt, vwind_expl_wgt_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddqz_z_half_d,         CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,      ddqz_z_half_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), d_exner_dz_ref_ic_d,   CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,      d_exner_dz_ref_ic_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_ref_ic_d,        CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,      theta_ref_ic_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), d_exner_dz_ref_mc_d,   CL_FALSE, 0, sizeof(cl_double)*size_3d_nlev,        d_exner_dz_ref_mc_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), d2_exner_dz2_ref_mc_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_nlev,        d2_exner_dz2_ref_mc_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfac_c_d,            CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,      wgtfac_c_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfacq_c_d,           CL_FALSE, 0, sizeof(cl_double)*size_wgtfacq_c,      wgtfacq_c_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_nnow_d,        CL_FALSE, 0, sizeof(cl_double)*size_3d_nlev,        theta_v_nnow_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_nnow_d,            CL_FALSE, 0, sizeof(cl_double)*size_3d_nlev,        rho_nnow_tmp,  0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		checkresultsCPP_(theta_ref_mc,theta_ref_mc_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"theta:exner_ref_mc");
		
		checkresultsCPP_(rho_refcorr_ic,rho_refcorr_ic_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"theta:rho_refcorr_ic");
		
		checkresultsCPP_(vwind_expl_wgt,vwind_expl_wgt_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,one,one,"theta:vwind_expl_wgt");
		
		checkresultsCPP_(ddqz_z_half,ddqz_z_half_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"theta:ddqz_z_half");
		
		checkresultsCPP_(d_exner_dz_ref_ic,d_exner_dz_ref_ic_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"theta:d_exner_dz_ref_ic");
		
		checkresultsCPP_(theta_ref_ic,theta_ref_ic_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"theta:theta_ref_ic");
		
		checkresultsCPP_(d_exner_dz_ref_mc,d_exner_dz_ref_mc_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"theta:d_exner_dz_ref_mc");
		
		checkresultsCPP_(d2_exner_dz2_ref_mc,d2_exner_dz2_ref_mc_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"theta:d2_exner_dz2_ref_mc");
		
		checkresultsCPP_(wgtfac_c,wgtfac_c_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"theta:wgtfac_c");
		
		checkresultsCPP_(wgtfacq_c,wgtfacq_c_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_c,nblks_c,"theta:wgtfacq_c");
		
		checkresultsCPP_(theta_v_nnow,theta_v_nnow_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"theta:theta_v_nnow");
		
		checkresultsCPP_(rho_nnow,rho_nnow_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"theta:rho_nnow");
		
		
		delete [] theta_ref_mc_tmp;
		delete [] rho_refcorr_ic_tmp;
		delete [] vwind_expl_wgt_tmp;
		delete [] ddqz_z_half_tmp;
		delete [] d_exner_dz_ref_ic_tmp;
		delete [] theta_ref_ic_tmp;
		delete [] d_exner_dz_ref_mc_tmp;
		delete [] d2_exner_dz2_ref_mc_tmp;
		delete [] wgtfac_c_tmp;
		delete [] wgtfacq_c_tmp;
		delete [] theta_v_nnow_tmp;
		delete [] rho_nnow_tmp;
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}	
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("theta1");
	}
#endif
	
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 4, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 5, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 6, sizeof(cl_mem), (void*)&theta_v_nvar_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 7, sizeof(cl_mem), (void*)&theta_v_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 8, sizeof(cl_mem), (void*)&theta_ref_mc_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta1, 9, sizeof(cl_mem), (void*)&z_theta_v_pr_mc_d) );
	
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_theta1[4] = { ev_write_localStart, ev_write_localEnd, ev_write_theta_v_nvar, ev_write_theta_ref_mc };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta1, 1, NULL, global_theta1, NULL, 4, ev_wait_list_theta1, &ev_kernel_theta1) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_theta1[1] = { ev_write_theta_v_nvar };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta1, 1, NULL, global_theta1, NULL, 1, ev_wait_list_theta1, &ev_kernel_theta1) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta1, 1, NULL, global_theta1, NULL, 0, NULL, &ev_kernel_theta1) );
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("theta2");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 7, sizeof(cl_mem), (void*)&rho_refcorr_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 8, sizeof(cl_mem), (void*)&wgtfac_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2, 9, sizeof(cl_mem), (void*)&rho_nvar_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,10, sizeof(cl_mem), (void*)&rho_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,11, sizeof(cl_mem), (void*)&z_theta_v_pr_mc_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,12, sizeof(cl_mem), (void*)&theta_ref_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,13, sizeof(cl_mem), (void*)&vwind_expl_wgt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,14, sizeof(cl_mem), (void*)&z_exner_pr_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,15, sizeof(cl_mem), (void*)&ddqz_z_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,16, sizeof(cl_mem), (void*)&d_exner_dz_ref_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,17, sizeof(cl_mem), (void*)&rho_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,18, sizeof(cl_mem), (void*)&z_theta_v_pr_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,19, sizeof(cl_mem), (void*)&theta_v_h_d) );
	CheckCL( clSetKernelArg(ocl.kernel_theta2,20, sizeof(cl_mem), (void*)&z_th_ddz_exner_c_d) );

	//cl_event ev_wait_list_theta2[15] = { ev_write_localStart, ev_write_localEnd, ev_write_rho_refcorr_ic, ev_write_wgtfac_c, ev_write_rho_nvar, ev_write_rho_nnow,
	//	                                 ev_write_theta_ref_ic, ev_write_vwind_expl_wgt, ev_write_z_exner_pr, ev_write_ddqz_z_half,
	//	                                 ev_write_d_exner_dz_ref_ic, ev_write_rho_ic, ev_write_theta_v_h, ev_write_z_th_ddz_exner_c, ev_kernel_theta1 };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta2, 3, NULL, global_theta2, NULL, 15, ev_wait_list_theta2, &ev_kernel_theta2) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_theta2[11] = { ev_write_localStart, ev_write_localEnd, ev_write_rho_refcorr_ic, ev_write_rho_nvar,
											 ev_write_theta_ref_ic, ev_write_vwind_expl_wgt, ev_write_ddqz_z_half,
		                                     ev_write_d_exner_dz_ref_ic, ev_write_rho_ic, ev_write_theta_v_h, ev_kernel_theta1 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta2, 3, NULL, global_theta2, local_size, 11, ev_wait_list_theta2, &ev_kernel_theta2) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_theta2[4] = { ev_write_rho_nvar, ev_write_rho_ic, ev_write_theta_v_h, ev_kernel_theta1 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta2, 3, NULL, global_theta2, local_size, 4, ev_wait_list_theta2, &ev_kernel_theta2) );
#endif
	}
	else
	{
		cl_event ev_wait_list_theta2[1] = { ev_kernel_theta1 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta2, 3, NULL, global_theta2, local_size, 1, ev_wait_list_theta2, &ev_kernel_theta2) );
	}
	
	if (*istep==1)
	{
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.pop_stop();
			profilerKernels.push_start("theta3a");
		}
		
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 2, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 3, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 4, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 5, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 6, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 7, sizeof(cl_mem), (void*)&wgtfacq_c_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 8, sizeof(cl_mem), (void*)&z_theta_v_pr_mc_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a, 9, sizeof(cl_mem), (void*)&theta_ref_ic_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a,10, sizeof(cl_mem), (void*)&z_theta_v_pr_ic_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3a,11, sizeof(cl_mem), (void*)&theta_v_h_d) );
		
		//cl_event ev_wait_list_theta3a[6] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq_c, ev_write_theta_ref_ic, ev_write_theta_v_h, ev_kernel_theta2 };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3a, 3, NULL, global_theta3a, NULL, 6, ev_wait_list_theta3a, &ev_kernel_theta3a) );
		if (ocl.bFirstIter)
		{
			cl_event ev_wait_list_theta3a[5] = { ev_write_localStart, ev_write_localEnd, ev_write_theta_ref_ic, ev_write_theta_v_h, ev_kernel_theta2 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3a, 2, NULL, global_theta3a, NULL, 5, ev_wait_list_theta3a, &ev_kernel_theta3a) );
#ifndef TURBO_MODE
		}
		else if (*istep==1)
		{
			cl_event ev_wait_list_theta3a[2] = { ev_write_theta_v_h, ev_kernel_theta2 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3a, 2, NULL, global_theta3a, NULL, 2, ev_wait_list_theta3a, &ev_kernel_theta3a) );
#endif
		}
		else
		{
			cl_event ev_wait_list_theta3a[1] = { ev_kernel_theta2 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3a, 2, NULL, global_theta3a, NULL, 1, ev_wait_list_theta3a, &ev_kernel_theta3a) );
		}
		
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.pop_stop();
			profilerKernels.push_start("theta3b");
		}
		
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 2, sizeof(int32_t), (void*)&slevm1) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 3, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 4, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 5, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 6, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 7, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 8, sizeof(cl_mem), (void*)&theta_v_nnow_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b, 9, sizeof(cl_mem), (void*)&theta_v_h_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b,10, sizeof(cl_mem), (void*)&z_theta_v_pr_ic_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b,11, sizeof(cl_mem), (void*)&d_exner_dz_ref_mc_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b,12, sizeof(cl_mem), (void*)&z_theta_v_pr_mc_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b,13, sizeof(cl_mem), (void*)&d2_exner_dz2_ref_mc_d) );
		CheckCL( clSetKernelArg(ocl.kernel_theta3b,14, sizeof(cl_mem), (void*)&z_dexner_dz_c_d) );
		
		//cl_event ev_wait_list_theta3b[8] = { ev_write_localStart, ev_write_localEnd, ev_write_theta_v_nnow, ev_write_theta_v_h, ev_write_d_exner_dz_ref_mc, ev_write_d2_exner_dz2_ref_mc, ev_write_z_dexner_dz_c, ev_kernel_theta3a };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3b, 3, NULL, global_theta3b, NULL, 8, ev_wait_list_theta3b, &ev_kernel_theta3b) );
		if (ocl.bFirstIter)
		{
			cl_event ev_wait_list_theta3b[6] = { ev_write_localStart, ev_write_localEnd, ev_write_theta_v_h, ev_write_d_exner_dz_ref_mc, ev_write_d2_exner_dz2_ref_mc, ev_kernel_theta3a };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3b, 3, NULL, global_theta3b, NULL, 6, ev_wait_list_theta3b, &ev_kernel_theta3b) );
#ifndef TURBO_MODE
		}
		else if (*istep==1)
		{
			cl_event ev_wait_list_theta3b[2] = { ev_write_theta_v_h, ev_kernel_theta3a };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3b, 3, NULL, global_theta3b, NULL, 2, ev_wait_list_theta3b, &ev_kernel_theta3b) );
#endif
		}
		else
		{
			cl_event ev_wait_list_theta3b[1] = { ev_kernel_theta3a };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_theta3b, 3, NULL, global_theta3b, NULL, 1, ev_wait_list_theta3b, &ev_kernel_theta3b) );
		}
	}
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
#ifndef TURBO_MODE
	if (*istep!=1)
	{
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_ic_d,           CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,     rho_ic,           0, NULL, &ev_read_rho_ic) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_h_d,        CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,     theta_v_h,        0, NULL, &ev_read_theta_v_h) );
	}
#endif
	
#ifdef TEST
	cl_event ev_read_z_dexner_dz_c, ev_read_z_th_ddz_exner_c;
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_ic_d,           CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,     rho_ic,           0, NULL, &ev_read_rho_ic) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_h_d,        CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,     theta_v_h,        0, NULL, &ev_read_theta_v_h) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_th_ddz_exner_c_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_nlevp1,       z_th_ddz_exner_c, 0, NULL, &ev_read_z_th_ddz_exner_c) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_dexner_dz_c_d,    CL_FALSE, 0, sizeof(cl_double)*size_z_dexner_dz_c, z_dexner_dz_c,    0, NULL, &ev_read_z_dexner_dz_c) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}

extern "C" void ocl_grad_(int * i_startblk, int * i_endblk, int * nblks_c, int * nblks_e,
						  int * nflat, int * nflat_gradp, int * nlev, int * nlevp1,
						  int * localStart, int * localEnd, int * nproma,
						  int * igradp_method, int * istep,
						  int * icidx, int * icblk, int * ikidx,
						  double * grav_o_cpd,
						  double * inv_dual_edge_length,
						  double * z_exner_ex_pr,
						  double * z_dexner_dz_c,
						  double * ddxn_z_full,
						  double * c_lin_e,
						  double * zdiff_gradp,
						  double * theta_v_nnow,
						  double * theta_v_h,
						  double * inv_ddqz_z_full,
						  double * z_gradh_exner,
						  double * z_hydro_corr)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_grad");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int i_endblkm1   = *i_endblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;

	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event /*ev_write_icidx, ev_write_icblk, */ev_write_ikidx;
	//cl_event ev_write_inv_dual_edge_length/*, ev_write_z_exner_ex_pr*/;
	cl_event ev_write_ddxn_z_full/*, ev_write_c_lin_e, ev_write_z_dexner_dz_c*/;
	cl_event ev_write_zdiff_gradp;
	cl_event /*ev_write_theta_v_nnow, */ev_write_theta_v_h/*, ev_write_inv_ddqz_z_full*/;
	cl_event ev_kernel_grad1, ev_kernel_grad2, ev_kernel_grad3, ev_kernel_grad4;
	//cl_event ev_read_z_gradh_exner/*, ev_read_z_hydro_corr*/;
	
	const size_t local_grad2[3] = { 1, (*nflat_gradp)-(*nflat)+1, 1 };
	
	const size_t global_grad1[3] = { nBlocks, (*nflat)-1, *nproma };
	const size_t global_grad2[3] = { nBlocks, (*nflat_gradp)-(*nflat)+1, *nproma };
	const size_t global_grad3[3] = { nBlocks, (*nlev)-(*nflat_gradp), *nproma };
	const size_t global_grad4[2] = { nBlocks, *nproma };
	
	const int size_localStart = nBlocks;
	const int size_localEnd   = nBlocks;
	const int size_ic         = (*nproma)*(*nblks_e)*2; // c_lin_e (*nproma)*2*(*nblks_e)
	const int size_2d         = (*nproma)*(*nblks_e);
	const int size_3d_c       = (*nproma)*(*nlev)*(*nblks_c);
	const int size_3d_e       = (*nproma)*(*nlev)*(*nblks_e);
	const int size_3d_c_p1    = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_4d_c       = 2*(*nproma)*(*nlev)*(*nblks_c);
	const int size_4d_e       = 2*(*nproma)*(*nlev)*(*nblks_e);
	
	cl_mem localStart_d           = resourcesCL.getBuffers("localStart4",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d             = resourcesCL.getBuffers("localEnd4",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem icidx_d                = resourcesCL.getBuffers("icidx",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_ic)[0];
	cl_mem icblk_d                = resourcesCL.getBuffers("icblk",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_ic)[0];
	cl_mem ikidx_d                = resourcesCL.getBuffers("ikidx",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_4d_e)[0];
	cl_mem inv_dual_edge_length_d = resourcesCL.getBuffers("inv_dual_edge_length", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_2d)[0];
	cl_mem z_exner_ex_pr_d        = resourcesCL.getBuffers("z_exner_ex_pr",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_c)[0];
	cl_mem z_gradh_exner_d        = resourcesCL.getBuffers("z_gradh_exner",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_e)[0];
	cl_mem c_lin_e_d              = resourcesCL.getBuffers("c_lin_e",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ic)[0];
	cl_mem z_dexner_dz_c_d        = resourcesCL.getBuffers("z_dexner_dz_c",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_4d_c)[0];
	cl_mem zdiff_gradp_d          = resourcesCL.getBuffers("zdiff_gradp",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_4d_e)[0];
	cl_mem theta_v_nnow_d         = resourcesCL.getBuffers("theta_v_nnow",         1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_c)[0];
	cl_mem theta_v_h_d            = resourcesCL.getBuffers("theta_v_h",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_c_p1)[0];
	cl_mem inv_ddqz_z_full_d      = resourcesCL.getBuffers("inv_ddqz_z_full",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_c)[0];
	cl_mem z_hydro_corr_d         = resourcesCL.getBuffers("z_hydro_corr",         1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_2d)[0];
	cl_mem ddxn_z_full_d          = resourcesCL.getBuffers("ddxn_z_full",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_e)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), icidx_d,                CL_FALSE, 0, sizeof(int)*size_ic,         icidx,                0, NULL, &ev_write_icidx) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), icblk_d,                CL_FALSE, 0, sizeof(int)*size_ic,         icblk,                0, NULL, &ev_write_icblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ikidx_d,                CL_FALSE, 0, sizeof(int)*size_4d_e,       ikidx,                0, NULL, &ev_write_ikidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddxn_z_full_d,          CL_FALSE, 0, sizeof(double)*size_3d_e,    ddxn_z_full,          0, NULL, &ev_write_ddxn_z_full) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), zdiff_gradp_d,          CL_FALSE, 0, sizeof(double)*size_4d_e,    zdiff_gradp,          0, NULL, &ev_write_zdiff_gradp) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_full_d,      CL_FALSE, 0, sizeof(double)*size_3d_c,    inv_ddqz_z_full,      0, NULL, &ev_write_inv_ddqz_z_full) );
	}
	
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,           CL_FALSE, 0, sizeof(int)*size_localStart, localStart,           0, NULL, &ev_write_localStart) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,             CL_FALSE, 0, sizeof(int)*size_localEnd,   localEnd,             0, NULL, &ev_write_localEnd) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), inv_dual_edge_length_d, CL_FALSE, 0, sizeof(double)*size_2d,      inv_dual_edge_length, 0, NULL, &ev_write_inv_dual_edge_length) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), c_lin_e_d,              CL_FALSE, 0, sizeof(double)*size_ic,      c_lin_e,              0, NULL, &ev_write_c_lin_e) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_v_nnow_d,         CL_FALSE, 0, sizeof(double)*size_3d_c,    theta_v_nnow,         0, NULL, &ev_write_theta_v_nnow) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_v_h_d,            CL_FALSE, 0, sizeof(double)*size_3d_c_p1, theta_v_h,            0, NULL, &ev_write_theta_v_h) );
	
#ifdef TEST
	cl_event ev_write_z_gradh_exner, ev_write_z_hydro_corr, ev_write_z_dexner_dz_c, ev_write_z_exner_ex_pr;
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_gradh_exner_d,        CL_FALSE, 0, sizeof(double)*size_3d_e,    z_gradh_exner,        0, NULL, &ev_write_z_gradh_exner) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_hydro_corr_d,         CL_FALSE, 0, sizeof(double)*size_2d,      z_hydro_corr,         0, NULL, &ev_write_z_hydro_corr) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_dexner_dz_c_d,        CL_FALSE, 0, sizeof(double)*size_4d_c,    z_dexner_dz_c,        0, NULL, &ev_write_z_dexner_dz_c) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_exner_ex_pr_d,        CL_FALSE, 0, sizeof(double)*size_3d_c,    z_exner_ex_pr,        0, NULL, &ev_write_z_exner_ex_pr) );
	
	
	{
		cout << "checking transfers for ocl_grad...\n";
		int * icidx_tmp             = new int[size_ic];
		int * icblk_tmp             = new int[size_ic];
		int * ikidx_tmp             = new int[size_4d_e];
		double * ddxn_z_full_tmp     = new double[size_3d_e];
		double * zdiff_gradp_tmp     = new double[size_4d_e];
		double * inv_ddqz_z_full_tmp = new double[size_3d_c];
		double * inv_dual_edge_length_tmp = new double[size_2d];
		double * c_lin_e_tmp              = new double[size_ic];
		double * theta_v_nnow_tmp         = new double[size_3d_c];
		double * theta_v_h_tmp            = new double[size_3d_c_p1];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		int zero[1] = {0};
		int one[1] = {1};
		int two[1] = {2};
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), icidx_d,           CL_FALSE, 0, sizeof(cl_int)*size_ic,      icidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), icblk_d,           CL_FALSE, 0, sizeof(cl_int)*size_ic,      icblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ikidx_d,           CL_FALSE, 0, sizeof(cl_int)*size_4d_e,    ikidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddxn_z_full_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d_e, ddxn_z_full_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), zdiff_gradp_d,     CL_FALSE, 0, sizeof(cl_double)*size_4d_e, zdiff_gradp_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_full_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_c, inv_ddqz_z_full_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), inv_dual_edge_length_d, CL_FALSE, 0, sizeof(cl_double)*size_2d,   inv_dual_edge_length_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), c_lin_e_d,              CL_FALSE, 0, sizeof(cl_double)*size_ic,   c_lin_e_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_nnow_d,         CL_FALSE, 0, sizeof(cl_double)*size_3d_c, theta_v_nnow_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_h_d,            CL_FALSE, 0, sizeof(cl_double)*size_3d_c_p1, theta_v_h_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		checkresultsCPP_(icidx,icidx_tmp,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,
						 zero,two,two,"grad:icidx");
		
		checkresultsCPP_(icblk,icblk_tmp,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,
						 zero,two,two,"grad:icblk");
		
		checkresults4d_int_(ikidx,ikidx_tmp,
						zero,two,two,
						zero,nproma,nproma,
						zero,nblks_e,nblks_e,
						zero,two,two,"grad:ikidx");
		
		checkresultsCPP_(ddxn_z_full,ddxn_z_full_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"grad:ddxn_z_full");
		
		checkresults4dCPP_(zdiff_gradp,zdiff_gradp_tmp,
						zero,two,two,
						zero,nproma,nproma,
						zero,nlev,nlev,
						zero,nblks_e,nblks_e,"grad:zdiff_gradp");
		
		checkresultsCPP_(inv_ddqz_z_full,inv_ddqz_z_full_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"grad:inv_ddqz_z_full");
		
		checkresultsCPP_(inv_dual_edge_length,inv_dual_edge_length_tmp,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,
						 zero,one,one,"grad:inv_dual_edge_length");
		
		checkresultsCPP_(c_lin_e,c_lin_e_tmp,
						 zero,nproma,nproma,
						 zero,two,two,
						 zero,nblks_e,nblks_e,"grad:c_lin_e");
		
		checkresultsCPP_(theta_v_nnow,theta_v_nnow_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"grad:theta_v_nnow");
		
		checkresultsCPP_(theta_v_h,theta_v_h_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"grad:theta_v_h");
		
		
		delete [] icidx_tmp;
		delete [] icblk_tmp;
		delete [] ikidx_tmp;
		delete [] ddxn_z_full_tmp;
		delete [] zdiff_gradp_tmp;
		delete [] inv_ddqz_z_full_tmp;
		delete [] inv_dual_edge_length_tmp;
		delete [] c_lin_e_tmp;
		delete [] theta_v_nnow_tmp;
		delete [] theta_v_h_tmp;
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("grad1");
	}
#endif
	
	const int elev1 = *nflat-1;
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 2, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 3, sizeof(int32_t), (void*)&elev1) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 4, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 5, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 8, sizeof(cl_mem), (void*)&icidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1, 9, sizeof(cl_mem), (void*)&icblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1,10, sizeof(cl_mem), (void*)&inv_dual_edge_length_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1,11, sizeof(cl_mem), (void*)&z_exner_ex_pr_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad1,12, sizeof(cl_mem), (void*)&z_gradh_exner_d) );
	
	//cl_event ev_wait_list_grad1[6] = { ev_write_localStart, ev_write_localEnd, ev_write_icidx, ev_write_icblk, ev_write_z_exner_ex_pr, ev_write_z_gradh_exner };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad1, 3, NULL, global_grad1, NULL, 6, ev_wait_list_grad1, &ev_kernel_grad1) );
	//cl_event ev_wait_list_grad1[2] = { ev_write_localStart, ev_write_localEnd };
	CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad1, 3, NULL, global_grad1, NULL, 0, NULL, &ev_kernel_grad1) );
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("grad2");
	}
	
	const int slev2 = *nflat-1;
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 2, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 3, sizeof(int32_t), (void*)&slev2) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 4, sizeof(int32_t), (void*)nflat_gradp) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 5, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 6, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 7, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 8, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2, 9, sizeof(cl_mem), (void*)&icidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2,10, sizeof(cl_mem), (void*)&icblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2,11, sizeof(cl_mem), (void*)&inv_dual_edge_length_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2,12, sizeof(cl_mem), (void*)&z_exner_ex_pr_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2,13, sizeof(cl_mem), (void*)&ddxn_z_full_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2,14, sizeof(cl_mem), (void*)&c_lin_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2,15, sizeof(cl_mem), (void*)&z_dexner_dz_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_grad2,16, sizeof(cl_mem), (void*)&z_gradh_exner_d) );

	//cl_event ev_wait_list_grad2[11] = { ev_write_localStart, ev_write_localEnd, ev_write_icidx, ev_write_icblk, ev_write_inv_dual_edge_length, ev_write_z_exner_ex_pr, ev_write_ddxn_z_full, ev_write_c_lin_e, ev_write_z_dexner_dz_c, ev_write_z_gradh_exner, ev_kernel_grad1 };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad2, 3, NULL, global_grad2, NULL, 11, ev_wait_list_grad2, &ev_kernel_grad2) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_grad2[2] = { ev_write_ddxn_z_full, ev_kernel_grad1 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad2, 3, NULL, global_grad2, local_grad2, 2, ev_wait_list_grad2, &ev_kernel_grad2) );
	}
	else
	{
		cl_event ev_wait_list_grad2[1] = { ev_kernel_grad1 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad2, 3, NULL, global_grad2, local_grad2, 1, ev_wait_list_grad2, &ev_kernel_grad2) );
	}
	
	if (*igradp_method >= 2)
	{	
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.pop_stop();
			profilerKernels.push_start("grad3");
		}
		
		cout << "untested\n";
		const int slev3 = *nflat_gradp;
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 0, sizeof(int32_t), (void*)&i_startblk) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 2, sizeof(int32_t), (void*)nblks_e) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 3, sizeof(int32_t), (void*)&slev3) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 4, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 5, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 6, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 7, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 8, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3, 9, sizeof(cl_mem), (void*)&icidx_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3,10, sizeof(cl_mem), (void*)&icblk_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3,11, sizeof(cl_mem), (void*)&ikidx_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3,12, sizeof(cl_mem), (void*)&inv_dual_edge_length_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3,13, sizeof(cl_mem), (void*)&z_exner_ex_pr_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3,14, sizeof(cl_mem), (void*)&zdiff_gradp_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3,15, sizeof(cl_mem), (void*)&z_dexner_dz_c_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad3,16, sizeof(cl_mem), (void*)&z_gradh_exner_d) );
	
		//cl_event ev_wait_list_grad3[11] = { ev_write_localStart, ev_write_localEnd, ev_write_icidx, ev_write_icblk, ev_write_ikidx, ev_write_inv_dual_edge_length, ev_write_z_exner_ex_pr, ev_write_zdiff_gradp, ev_write_z_dexner_dz_c, ev_write_z_gradh_exner, ev_kernel_grad2 };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad3, 3, NULL, global_grad3, NULL, 11, ev_wait_list_grad3, &ev_kernel_grad3) );
		if (ocl.bFirstIter)
		{
			cl_event ev_wait_list_grad3[3] = { ev_write_ikidx, ev_write_zdiff_gradp, ev_kernel_grad2 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad3, 3, NULL, global_grad3, NULL, 3, ev_wait_list_grad3, &ev_kernel_grad3) );
		}
		else
		{
			cl_event ev_wait_list_grad3[1] = { ev_kernel_grad2 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad3, 3, NULL, global_grad3, NULL, 1, ev_wait_list_grad3, &ev_kernel_grad3) );
		}
	}
	if (*igradp_method == 3)
	{	
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.pop_stop();
			profilerKernels.push_start("grad4");
		}
		cout << "untested\n";
		/*
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 0, sizeof(int32_t), (void*)&i_startblk) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 2, sizeof(int32_t), (void*)nblks_e) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 3, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 4, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 5, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 6, sizeof(int32_t), (void*)grav_o_cpd) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 7, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 8, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4, 9, sizeof(cl_mem), (void*)&icidx_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4,10, sizeof(cl_mem), (void*)&icblk_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4,11, sizeof(cl_mem), (void*)&ikidx_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4,12, sizeof(cl_mem), (void*)&theta_v_nnow_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4,13, sizeof(cl_mem), (void*)&zdiff_gradp_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4,14, sizeof(cl_mem), (void*)&theta_v_h_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4,15, sizeof(cl_mem), (void*)&inv_ddqz_z_full_d) );
		CheckCL( clSetKernelArg(ocl.kernel_grad4,16, sizeof(cl_mem), (void*)&inv_dual_edge_length_d) );
		//CheckCL( clSetKernelArg(ocl.kernel_grad4,17, sizeof(cl_mem), (void*)&z_hydro_corr_d) );
		
		//cl_event ev_wait_list_grad4[12] = { ev_write_localStart, ev_write_localEnd, ev_write_icidx, ev_write_icblk, ev_write_ikidx, ev_write_theta_v, ev_write_zdiff_gradp, ev_write_theta_v_h, ev_write_inv_ddqz_z_full, ev_write_inv_dual_edge_length, ev_write_z_hydro_corr, ev_kernel_grad3 };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad4, 2, NULL, global_grad4, NULL, 12, ev_wait_list_grad4, &ev_kernel_grad4) );
		cl_event ev_wait_list_grad4[11] = { ev_write_localStart, ev_write_localEnd, ev_write_icidx, ev_write_icblk, ev_write_ikidx, ev_write_theta_v, ev_write_zdiff_gradp, ev_write_theta_v_h, ev_write_inv_ddqz_z_full, ev_write_inv_dual_edge_length, ev_kernel_grad3 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_grad4, 2, NULL, global_grad4, NULL, 11, ev_wait_list_grad4, &ev_kernel_grad4) );
		//*/
	}
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
#ifdef TEST
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	cl_event ev_read_z_gradh_exner, ev_read_z_hydro_corr;
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_gradh_exner_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_e, z_gradh_exner, 0, NULL, &ev_read_z_gradh_exner) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_hydro_corr_d,  CL_FALSE, 0, sizeof(cl_double)*size_2d,   z_hydro_corr,  0, NULL, &ev_read_z_hydro_corr) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling) profiler.pop_stop();
}


extern "C" void ocl_vn_(int * i_startblk, int * i_endblk, int * nblks_e,
						int * nlev,
						int * localStart, int * localEnd, int * nproma,
						int * itime_scheme, int * istep,
						int * ntl1, int * ntl2,
						double * cpd, double * dtime,
						double * vn_nnow,
						double * ddt_vn_adv,
						double * ddt_vn_phy,
						double * z_theta_v_e,
						double * z_gradh_exner,
						double * vn_nnew)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_vn");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;
	
	if (ocl.bFirstIter)
		for (int i=0; i<nBlocks; i++)
		{
			if (localStart[i]!=0 || localEnd[i]!=1)
			{
				cout << "need to use 3d indices\n";
				abort();
			}
		}
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_vn_nnow/*, ev_write_ddt_vn_adv*/, ev_write_ddt_vn_phy, ev_write_z_theta_v_e/*, ev_write_z_gradh_exner*/, ev_write_vn_nnew;
	cl_event ev_kernel_vn1;
	cl_event ev_read_vn_nnew;
	
	const size_t global_vn1[1] = { nBlocks*(*nlev)*(*nproma) };
	
	const int size_localStart = nBlocks;
	const int size_localEnd   = nBlocks;
	const int size_2d         = (*nproma)*(*nblks_e);
	const int size_3d         = (*nproma)*(*nlev)*(*nblks_e);
	const int size_4d_3       = (*nproma)*(*nlev)*(*nblks_e)*3;
	
	cl_mem localStart_d           = resourcesCL.getBuffers("localStart4",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d             = resourcesCL.getBuffers("localEnd4",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem vn_nnow_d              = resourcesCL.getBuffers("vn_nnow",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem ddt_vn_adv_d           = resourcesCL.getBuffers("ddt_vn_adv",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_4d_3)[0];
	cl_mem ddt_vn_phy_d           = resourcesCL.getBuffers("ddt_vn_phy",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_theta_v_e_d          = resourcesCL.getBuffers("z_theta_v_e",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_gradh_exner_d        = resourcesCL.getBuffers("z_gradh_exner",         1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem vn_nnew_d              = resourcesCL.getBuffers("vn_nnew",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,           CL_FALSE, 0, sizeof(int)*size_localStart, localStart,         0, NULL, &ev_write_localStart) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,             CL_FALSE, 0, sizeof(int)*size_localEnd,   localEnd,           0, NULL, &ev_write_localEnd) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_nnow_d,              CL_FALSE, 0, sizeof(double)*size_3d,      vn_nnow,            0, NULL, &ev_write_vn_nnow) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_vn_adv_d,           CL_FALSE, 0, sizeof(double)*size_4d_3,    ddt_vn_adv,         0, NULL, &ev_write_ddt_vn_adv) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_vn_phy_d,           CL_FALSE, 0, sizeof(double)*size_3d,      ddt_vn_phy,         0, NULL, &ev_write_ddt_vn_phy) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_theta_v_e_d,          CL_FALSE, 0, sizeof(double)*size_3d,      z_theta_v_e,        0, NULL, &ev_write_z_theta_v_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_gradh_exner_d,        CL_FALSE, 0, sizeof(double)*size_3d,      z_gradh_exner,      0, NULL, &ev_write_z_gradh_exner) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_nnew_d,              CL_FALSE, 0, sizeof(double)*size_3d,      vn_nnew,            0, NULL, &ev_write_vn_nnew) );
	}
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_vn...\n";
		double * ddt_vn_adv_tmp    = new double[size_4d_3];
		double * z_gradh_exner_tmp = new double[size_3d];
		double * z_theta_v_e_tmp   = new double[size_3d];
		int zero[1] = {0};
		int three[1] = {3};
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_vn_adv_d,    CL_FALSE, 0, sizeof(cl_double)*size_4d_3, ddt_vn_adv_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_gradh_exner_d, CL_FALSE, 0, sizeof(cl_double)*size_3d,   z_gradh_exner_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_theta_v_e_d,   CL_FALSE, 0, sizeof(cl_double)*size_3d,   z_theta_v_e_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		checkresults4dCPP_(ddt_vn_adv,ddt_vn_adv_tmp,
						zero,nproma,nproma,
						zero,nlev,nlev,
						zero,nblks_e,nblks_e,
						zero,three,three,"vn:ddt_vn_adv");
		
		checkresultsCPP_(z_gradh_exner,z_gradh_exner_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"vn:z_gradh_exner");
		
		checkresultsCPP_(z_theta_v_e,z_theta_v_e_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"vn:z_theta_v_e");
		
		
		delete [] ddt_vn_adv_tmp;
		delete [] z_gradh_exner_tmp;
		delete [] z_theta_v_e_tmp;
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
	}
#endif
	
	const int ntl1m1 = *ntl1-1;
	const int ntl2m1 = *ntl2-1;
	if (*itime_scheme==6 && *istep==2)
	{
		if (bProfiling)
			profilerKernels.push_start("vn1a");
		
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 2, sizeof(int32_t), (void*)nblks_e) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 3, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 4, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 5, sizeof(int32_t), (void*)&ntl1m1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 6, sizeof(int32_t), (void*)&ntl2m1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 7, sizeof(double), (void*)dtime) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 8, sizeof(double), (void*)cpd) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a, 9, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a,10, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a,11, sizeof(cl_mem), (void*)&vn_nnow_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a,12, sizeof(cl_mem), (void*)&ddt_vn_adv_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a,13, sizeof(cl_mem), (void*)&ddt_vn_phy_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a,14, sizeof(cl_mem), (void*)&z_theta_v_e_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a,15, sizeof(cl_mem), (void*)&z_gradh_exner_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1a,16, sizeof(cl_mem), (void*)&vn_nnew_d) );

		//cl_event ev_wait_list_vn1a[8] = { ev_write_localStart, ev_write_localEnd, ev_write_vn_nnow, ev_write_ddt_vn_adv, ev_write_ddt_vn_phy, ev_write_z_theta_v_e, ev_write_z_gradh_exner, ev_write_vn_nnew };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn1a, 3, NULL, global_vn1, NULL, 8, ev_wait_list_vn1a, &ev_kernel_vn1) );
#ifndef TURBO_MODE
		if (*istep==1)
		{
			cl_event ev_wait_list_vn1a[3] = { ev_write_vn_nnow, ev_write_ddt_vn_phy, ev_write_vn_nnew };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn1a, 1, NULL, global_vn1, NULL, 3, ev_wait_list_vn1a, &ev_kernel_vn1) );
		}
		else
#endif
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn1a, 1, NULL, global_vn1, NULL, 0, NULL, &ev_kernel_vn1) );
	}
	else
	{
		if (bProfiling)
			profilerKernels.push_start("vn1b");
		
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 2, sizeof(int32_t), (void*)nblks_e) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 3, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 4, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 5, sizeof(int32_t), (void*)&ntl1m1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 6, sizeof(double), (void*)dtime) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 7, sizeof(double), (void*)cpd) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 8, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b, 9, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b,10, sizeof(cl_mem), (void*)&vn_nnow_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b,11, sizeof(cl_mem), (void*)&ddt_vn_adv_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b,12, sizeof(cl_mem), (void*)&ddt_vn_phy_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b,13, sizeof(cl_mem), (void*)&z_theta_v_e_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b,14, sizeof(cl_mem), (void*)&z_gradh_exner_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn1b,15, sizeof(cl_mem), (void*)&vn_nnew_d) );
		
		//cl_event ev_wait_list_vn1b[8] = { ev_write_localStart, ev_write_localEnd, ev_write_vn_nnow, ev_write_ddt_vn_adv, ev_write_ddt_vn_phy, ev_write_z_theta_v_e, ev_write_z_gradh_exner, ev_write_vn_nnew };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn1b, 3, NULL, global_vn1, NULL, 8, ev_wait_list_vn1b, &ev_kernel_vn1) );
#ifndef TURBO_MODE
		if (*istep==1)
		{
			cl_event ev_wait_list_vn1b[3] = { ev_write_vn_nnow, ev_write_ddt_vn_phy, ev_write_vn_nnew };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn1b, 1, NULL, global_vn1, NULL, 3, ev_wait_list_vn1b, &ev_kernel_vn1) );
		}
		else
#endif
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn1b, 1, NULL, global_vn1, NULL, 0, NULL, &ev_kernel_vn1) );
	}
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
#ifndef TURBO_MODE
	if (*istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_nnew_d, CL_FALSE, 0, sizeof(cl_double)*size_3d, vn_nnew, 0, NULL, &ev_read_vn_nnew) );
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_nnew_d, CL_FALSE, 0, sizeof(cl_double)*size_3d, vn_nnew, 0, NULL, &ev_read_vn_nnew) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}


extern "C" void ocl_vn_avg_(int * i_startblk, int * i_endblk, int * nblks_e,
							int * nlev,
							int * localStart, int * localEnd, int * nproma,
							int * istep,
							int * iqidx, int * iqblk,
							double * vn_nnew,
							double * e_flx_avg,
							double * z_vn_avg)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_vn_avg");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_iqidx, ev_write_iqblk, ev_write_vn_nnew, ev_write_e_flx_avg;
	cl_event ev_kernel_z_vn_avg;
	
	const size_t local_vn_avg[3] = { 1, *nlev, 1 };
	const size_t global_vn_avg[3] = { nBlocks, *nlev, *nproma };
	
	const int size_localStart = nBlocks;
	const int size_localEnd   = nBlocks;
	const int size_iq         = (*nproma)*(*nblks_e)*4;
	const int size_3d         = (*nproma)*(*nlev)*(*nblks_e);
	const int size_e_flx_avg  = (*nproma)*5*(*nblks_e);
	
	cl_mem localStart_d           = resourcesCL.getBuffers("localStartVnAvg",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d             = resourcesCL.getBuffers("localEndVnAvg",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem iqidx_d                = resourcesCL.getBuffers("iqidx",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iq)[0];
	cl_mem iqblk_d                = resourcesCL.getBuffers("iqblk",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iq)[0];
	cl_mem vn_nnew_d              = resourcesCL.getBuffers("vn_nnew",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem e_flx_avg_d            = resourcesCL.getBuffers("e_flx_avg",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_e_flx_avg)[0];
	cl_mem z_vn_avg_d             = resourcesCL.getBuffers("z_vn_avg",             1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iqidx_d,                CL_FALSE, 0, sizeof(int)*size_iq,           iqidx,              0, NULL, &ev_write_iqidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iqblk_d,                CL_FALSE, 0, sizeof(int)*size_iq,           iqblk,              0, NULL, &ev_write_iqblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,           CL_FALSE, 0, sizeof(int)*size_localStart,   localStart,         0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,             CL_FALSE, 0, sizeof(int)*size_localEnd,     localEnd,           0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_nnew_d,              CL_FALSE, 0, sizeof(double)*size_3d,        vn_nnew,            0, NULL, &ev_write_vn_nnew) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), e_flx_avg_d,            CL_FALSE, 0, sizeof(double)*size_e_flx_avg, e_flx_avg,          0, NULL, &ev_write_e_flx_avg) );
	}
#ifdef TEST
	cl_event ev_write_z_vn_avg;
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_vn_avg_d,             CL_FALSE, 0, sizeof(double)*size_3d,        z_vn_avg,           0, NULL, &ev_write_z_vn_avg) );
	
	{
		cout << "checking transfers for ocl_vn_avg...\n";
		int * iqidx_tmp = new int[size_iq];
		int * iqblk_tmp = new int[size_iq];
		double * vn_nnew_tmp = new double[size_3d];
		
		int zero[1] = {0};
		int four[1] = {4};
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iqidx_d, CL_FALSE, 0, sizeof(cl_int)*size_iq, iqidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iqblk_d, CL_FALSE, 0, sizeof(cl_int)*size_iq, iqblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_nnew_d, CL_FALSE, 0, sizeof(cl_double)*size_3d, vn_nnew_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		checkresultsCPP_(iqidx,iqidx_tmp,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,
						 zero,four,four,"vn_avg:iqidx");
		
		checkresultsCPP_(iqblk,iqblk_tmp,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,
						 zero,four,four,"vn_avg:iqblk");
		
		checkresultsCPP_(vn_nnew,vn_nnew_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"vn_avg:vn_nnew");
		
		delete [] iqidx_tmp;
		delete [] iqblk_tmp;
		delete [] vn_nnew_tmp;
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("vn_avg");
	}
#endif
	
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 2, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 7, sizeof(cl_mem), (void*)&iqidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 8, sizeof(cl_mem), (void*)&iqblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg, 9, sizeof(cl_mem), (void*)&vn_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg,10, sizeof(cl_mem), (void*)&e_flx_avg_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_avg,11, sizeof(cl_mem), (void*)&z_vn_avg_d) );
	
	//cl_event ev_wait_list_vn_avg[7] = { ev_write_localStart, ev_write_localEnd, ev_write_iqidx, ev_write_iqblk, ev_write_vn_nnew, ev_write_e_flx_avg, ev_write_z_vn_avg };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_avg, 3, NULL, global_vn_avg, NULL, 7, ev_wait_list_vn_avg, &ev_kernel_z_vn_avg) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_vn_avg[5] = { ev_write_localStart, ev_write_localEnd, ev_write_iqidx, ev_write_iqblk, ev_write_e_flx_avg };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_avg, 3, NULL, global_vn_avg, local_vn_avg, 5, ev_wait_list_vn_avg, &ev_kernel_z_vn_avg) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_vn_avg[1] = { ev_write_e_flx_avg };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_avg, 3, NULL, global_vn_avg, local_vn_avg, 1, ev_wait_list_vn_avg, &ev_kernel_z_vn_avg) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_avg, 3, NULL, global_vn_avg, local_vn_avg, 0, NULL, &ev_kernel_z_vn_avg) );
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
		//profiler.pop_stop();
	}
#endif
	
#ifdef TEST
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	cl_event ev_read_z_vn_avg;
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_vn_avg_d, CL_FALSE, 0, sizeof(cl_double)*size_3d, z_vn_avg, 0, NULL, &ev_read_z_vn_avg) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling) profiler.pop_stop();
}

extern "C" void ocl_vn_half_(int * i_startblk, int * i_endblk, int * nblks_e,
							 int * nflat, int * nlev, int * nlevp1,
							 int * localStart, int * localEnd, int * nproma,
							 int * istep, bool * l_vert_nested,
							 double * wgtfac_e,
							 double * wgtfacq_e,
							 double * wgtfacq1_e,
							 double * vn_nnew,
							 double * vt,
							 double * vn_half,
							 double * ddxn_z_half,
							 double * ddxt_z_half,
							 double * z_concorr_e)
{
	// read notes in the kernels for details about vn_half5 (computationally different than the original code
	//  but delivers same results with 4-5x performance improvement
	
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("vn_half");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;
	
	if (ocl.bFirstIter)
		for (int i=0; i<nBlocks; i++)
		{
			if (localStart[i]!=0 || localEnd[i]!=1)
			{
				cout << "need to use 3d indices\n";
				abort();
			}
		}
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event /*ev_write_wgtfac_e, ev_write_vn_nnew, */ev_write_vt/*, ev_write_vn_half*/;
	//cl_event ev_write_wgtfacq1_e, ev_write_wgtfacq_e;
	//cl_event ev_write_ddxn_z_half, ev_write_ddxt_z_half;
	cl_event ev_write_z_concorr_e;
	cl_event ev_kernel_vn_half1, ev_kernel_vn_half2, ev_kernel_vn_half3, ev_kernel_vn_half4, ev_kernel_vn_half5;
	cl_event ev_read_vn_half, ev_read_z_concorr_e;
	
	const size_t local_vn_half[3] = { 1, *nlev, 1 };
	
	const size_t global_vn_half1[3] = { nBlocks, (*nlev), *nproma };
	const size_t global_vn_half2[2] = { nBlocks, *nproma };
	const size_t global_vn_half3[3] = { nBlocks, (*nflat)-2, *nproma };
	const size_t global_vn_half4[2] = { nBlocks, *nproma };
	const size_t global_vn_half5[1] = { nBlocks*(*nlevp1)*(*nproma) };
	
	const int size_localStart = nBlocks;
	const int size_localEnd   = nBlocks;
	const int size_3d         = (*nproma)*(*nlev)*(*nblks_e);
	const int size_3d_p1      = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_wgtfacq    = (*nproma)*3*(*nblks_e);
	
	cl_mem localStart_d           = resourcesCL.getBuffers("localStart1",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d             = resourcesCL.getBuffers("localEnd1",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem wgtfac_e_d             = resourcesCL.getBuffers("wgtfac_e",             1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem vn_nnew_d              = resourcesCL.getBuffers("vn_nnew",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem vt_d                   = resourcesCL.getBuffers("vt",                   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem vn_half_d              = resourcesCL.getBuffers("vn_half",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_vt_ie_d              = resourcesCL.getBuffers("z_vt_ie",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem wgtfacq1_e_d           = resourcesCL.getBuffers("wgtfacq1_e",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_wgtfacq)[0];
	cl_mem wgtfacq_e_d            = resourcesCL.getBuffers("wgtfacq_e",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_wgtfacq)[0];
	cl_mem ddxn_z_half_d          = resourcesCL.getBuffers("ddxn_z_half",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem ddxt_z_half_d          = resourcesCL.getBuffers("ddxt_z_half",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_concorr_e_d          = resourcesCL.getBuffers("z_concorr_e",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfac_e_d,             CL_FALSE, 0, sizeof(double)*size_3d_p1,     wgtfac_e,           0, NULL, &ev_write_wgtfac_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfacq1_e_d,           CL_FALSE, 0, sizeof(double)*size_wgtfacq,   wgtfacq1_e,         0, NULL, &ev_write_wgtfacq1_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfacq_e_d,            CL_FALSE, 0, sizeof(double)*size_wgtfacq,   wgtfacq_e,          0, NULL, &ev_write_wgtfacq_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddxn_z_half_d,          CL_FALSE, 0, sizeof(double)*size_3d_p1,     ddxn_z_half,        0, NULL, &ev_write_ddxn_z_half) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddxt_z_half_d,          CL_FALSE, 0, sizeof(double)*size_3d_p1,     ddxt_z_half,        0, NULL, &ev_write_ddxt_z_half) );
	}
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,           CL_FALSE, 0, sizeof(int)*size_localStart,   localStart,         0, NULL, &ev_write_localStart) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,             CL_FALSE, 0, sizeof(int)*size_localEnd,     localEnd,           0, NULL, &ev_write_localEnd) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_nnew_d,              CL_FALSE, 0, sizeof(double)*size_3d,        vn_nnew,            0, NULL, &ev_write_vn_nnew) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vt_d,                   CL_FALSE, 0, sizeof(double)*size_3d,        vt,                 0, NULL, &ev_write_vt) );
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_half_d,              CL_FALSE, 0, sizeof(double)*size_3d_p1,     vn_half,            0, NULL, &ev_write_vn_half) );
	
#ifdef TEST
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_concorr_e_d,          CL_FALSE, 0, sizeof(double)*size_3d_p1,     z_concorr_e,        0, NULL, &ev_write_z_concorr_e) );
	
		cout << "checking transfers for vn_half...\n";
		double * vn_half_tmp     = new double[size_3d_p1];
		double * wgtfac_e_tmp    = new double[size_3d_p1];
		double * wgtfacq_e_tmp   = new double[size_wgtfacq];
		double * wgtfacq1_e_tmp  = new double[size_wgtfacq];
		double * ddxn_z_half_tmp = new double[size_3d_p1];
		double * ddxt_z_half_tmp = new double[size_3d_p1];
		double * vn_nnew_tmp     = new double[size_3d];
		double * z_concorr_e_tmp = new double[size_3d_p1];
		double * vt_tmp          = new double[size_3d];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		int zero[1] = {0};
		int three[1] = {3};
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_half_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d_p1,   vn_half_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfac_e_d,    CL_FALSE, 0, sizeof(cl_double)*size_3d_p1,   wgtfac_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfacq_e_d,   CL_FALSE, 0, sizeof(cl_double)*size_wgtfacq, wgtfacq_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfacq1_e_d,  CL_FALSE, 0, sizeof(cl_double)*size_wgtfacq, wgtfacq1_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddxn_z_half_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_p1,   ddxn_z_half_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddxt_z_half_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_p1,   ddxt_z_half_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_nnew_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d,      vn_nnew_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_concorr_e_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_p1,   z_concorr_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vt_d,          CL_FALSE, 0, sizeof(cl_double)*size_3d,      vt_tmp,  0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		checkresultsCPP_(vn_half,vn_half_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_e,nblks_e,"vn_half:vn_half");
		
		checkresultsCPP_(wgtfac_e,wgtfac_e_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_e,nblks_e,"vn_half:wgtfac_e");
		
		checkresultsCPP_(wgtfacq_e,wgtfacq_e_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_e,nblks_e,"vn_half:wgtfacq_e");
		
		checkresultsCPP_(wgtfacq1_e,wgtfacq1_e_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_e,nblks_e,"vn_half:wgtfacq1_e");
		
		checkresultsCPP_(ddxn_z_half,ddxn_z_half_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_e,nblks_e,"vn_half:ddxn_z_half");
		
		checkresultsCPP_(ddxt_z_half,ddxt_z_half_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_e,nblks_e,"vn_half:ddxt_z_half");
		
		checkresultsCPP_(vn_nnew,vn_nnew_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"vn_half:vn_nnew");
		
		checkresultsCPP_(z_concorr_e,z_concorr_e_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_e,nblks_e,"vn_half:z_concorr_e");
		
		checkresultsCPP_(vt,vt_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"vn_half:vt");
		
		delete [] vn_half_tmp;
		delete [] wgtfac_e_tmp;
		delete [] wgtfacq_e_tmp;
		delete [] wgtfacq1_e_tmp;
		delete [] ddxn_z_half_tmp;
		delete [] ddxt_z_half_tmp;
		delete [] vn_nnew_tmp;
		delete [] z_concorr_e_tmp;
		delete [] vt_tmp;
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("vn_half1");
	}
#endif
	
	const int nflatm1 = *nflat-1;
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 2, sizeof(int32_t), (void*)&nflatm1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 4, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 5, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 8, sizeof(cl_mem), (void*)&wgtfac_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1, 9, sizeof(cl_mem), (void*)&vn_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1,10, sizeof(cl_mem), (void*)&vt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1,11, sizeof(cl_mem), (void*)&vn_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half1,12, sizeof(cl_mem), (void*)&z_vt_ie_d) );
	
	//cl_event ev_wait_list_vn_half1[6] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfac_e, ev_write_vn_nnew, ev_write_vt, ev_write_vn_half };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half1, 3, NULL, global_vn_half1, NULL, 6, ev_wait_list_vn_half1, &ev_kernel_vn_half1) );
	CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half1, 3, NULL, global_vn_half1, local_vn_half, 0, NULL, &ev_kernel_vn_half1) );
	
	if (*istep==1 && !(*l_vert_nested)==1)
	{	
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.pop_stop();
			profilerKernels.push_start("vn_half2");
		}
		
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 2, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 3, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 4, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 5, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 6, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 7, sizeof(cl_mem), (void*)&wgtfacq1_e_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 8, sizeof(cl_mem), (void*)&vn_nnew_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half2, 9, sizeof(cl_mem), (void*)&vn_half_d) );
		
		//cl_event ev_wait_list_vn_half2[5] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq1_e, ev_write_vn_nnew, ev_write_vn_half };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half2, 2, NULL, global_vn_half2, NULL, 5, ev_wait_list_vn_half2, &ev_kernel_vn_half2) );
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half2, 2, NULL, global_vn_half2, NULL, 0, NULL, &ev_kernel_vn_half2) );
	}
	
	if (*istep==1)
	{	
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.pop_stop();
			profilerKernels.push_start("vn_half3");
		}
		
		const int elev = *nflat-1;
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 2, sizeof(int32_t), (void*)&elev) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 3, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 4, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 5, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 6, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 7, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 8, sizeof(cl_mem), (void*)&wgtfac_e_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3, 9, sizeof(cl_mem), (void*)&vn_nnew_d) );
		CheckCL( clSetKernelArg(ocl.kernel_vn_half3,10, sizeof(cl_mem), (void*)&vn_half_d) );
		
		//cl_event ev_wait_list_vn_half3[5] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfac_e, ev_write_vn_nnew, ev_write_vn_half };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half3, 3, NULL, global_vn_half3, NULL, 5, ev_wait_list_vn_half3, &ev_kernel_vn_half3) );
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half3, 3, NULL, global_vn_half3, NULL, 0, NULL, &ev_kernel_vn_half3) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("vn_half4");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 7, sizeof(cl_mem), (void*)&wgtfacq_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 8, sizeof(cl_mem), (void*)&vn_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4, 9, sizeof(cl_mem), (void*)&vt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4,10, sizeof(cl_mem), (void*)&vn_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half4,11, sizeof(cl_mem), (void*)&z_vt_ie_d) );
	
	//cl_event ev_wait_list_vn_half4[6] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq_e, ev_write_vn_nnew, ev_write_vt, ev_write_vn_half };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half4, 2, NULL, global_vn_half4, NULL, 6, ev_wait_list_vn_half4, &ev_kernel_vn_half4) );
	CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half4, 2, NULL, global_vn_half4, NULL, 0, NULL, &ev_kernel_vn_half4) );
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("vn_half5");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 2, sizeof(int32_t), (void*)nflat) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 7, sizeof(cl_mem), (void*)&vn_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 8, sizeof(cl_mem), (void*)&ddxn_z_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5, 9, sizeof(cl_mem), (void*)&z_vt_ie_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5,10, sizeof(cl_mem), (void*)&ddxt_z_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_vn_half5,11, sizeof(cl_mem), (void*)&z_concorr_e_d) );
	
	// some kernels should also be waited
	if (*istep==1 && !(*l_vert_nested)==1)
	{
		cl_event ev_wait_list_vn_half5[4] = { ev_kernel_vn_half1, ev_kernel_vn_half2, ev_kernel_vn_half3, ev_kernel_vn_half4 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half5, 1, NULL, global_vn_half5, NULL, 4, ev_wait_list_vn_half5, &ev_kernel_vn_half5) );
	}
	else if (*istep==1 && (*l_vert_nested)==1)
	{
		cl_event ev_wait_list_vn_half5[3] = { ev_kernel_vn_half1, ev_kernel_vn_half3, ev_kernel_vn_half4 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half5, 1, NULL, global_vn_half5, NULL, 3, ev_wait_list_vn_half5, &ev_kernel_vn_half5) );
	}
	else
	{
		cl_event ev_wait_list_vn_half5[2] = { ev_kernel_vn_half1, ev_kernel_vn_half4 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_vn_half5, 1, NULL, global_vn_half5, NULL, 2, ev_wait_list_vn_half5, &ev_kernel_vn_half5) );
	}
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
#ifndef TURBO_MODE
	if (*istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_half_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, vn_half,     0, NULL, &ev_read_vn_half) );
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_half_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, vn_half,     0, NULL, &ev_read_vn_half) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_concorr_e_d, CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, z_concorr_e, 0, NULL, &ev_read_z_concorr_e) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}

extern "C" void ocl_fluxes_(int * i_startblk, int * i_endblk, int * nblks_e,
							int * nlev,
							int * localStart, int * localEnd, int * nproma,
							int * istep,
							double * z_rho_e,
							double * ptr_vn, int * nameIdx,
							double * ddqz_z_full_e,
							double * z_theta_v_e,
							double * mass_fl_e,
							double * z_theta_v_fl_e)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_fluxes");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;
	
	if (ocl.bFirstIter)
		for (int i=0; i<nBlocks; i++)
		{
			if (localStart[i]!=0 || localEnd[i]!=1)
			{
				cout << "need to use 3d indices\n";
				abort();
			}
		}
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_z_rho_e, ev_write_ptr_vn/*, ev_write_ddqz_z_full_e, ev_write_z_theta_v_e*/, ev_write_mass_fl_e;
	cl_event ev_kernel_fluxes;
	cl_event ev_read_mass_fl_e;
	
	
	const size_t global_fluxes[1] = { nBlocks*(*nlev)*(*nproma) };
	
	const int size_localStart = nBlocks;
	const int size_localEnd   = nBlocks;
	const int size_3d         = (*nproma)*(*nlev)*(*nblks_e);
	
	cl_mem localStart_d           = resourcesCL.getBuffers("localStart",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d             = resourcesCL.getBuffers("localEnd",             1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem z_rho_e_d              = resourcesCL.getBuffers("z_rho_e",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem ptr_vn_d               = resourcesCL.getBuffers("ptr_vn",               1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem ddqz_z_full_e_d        = resourcesCL.getBuffers("ddqz_z_full_e",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_theta_v_e_d          = resourcesCL.getBuffers("z_theta_v_e",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem mass_fl_e_d            = resourcesCL.getBuffers("mass_fl_e",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_theta_v_fl_e_d       = resourcesCL.getBuffers("z_theta_v_fl_e",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,           CL_FALSE, 0, sizeof(int)*size_localStart,   localStart,         0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,             CL_FALSE, 0, sizeof(int)*size_localEnd,     localEnd,           0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_rho_e_d,              CL_FALSE, 0, sizeof(double)*size_3d,        z_rho_e,            0, NULL, &ev_write_z_rho_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddqz_z_full_e_d,        CL_FALSE, 0, sizeof(double)*size_3d,        ddqz_z_full_e,      0, NULL, &ev_write_ddqz_z_full_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_theta_v_e_d,          CL_FALSE, 0, sizeof(double)*size_3d,        z_theta_v_e,        0, NULL, &ev_write_z_theta_v_e) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), mass_fl_e_d,            CL_FALSE, 0, sizeof(double)*size_3d,        mass_fl_e,          0, NULL, &ev_write_mass_fl_e) );
	}
	
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ptr_vn_d,               CL_FALSE, 0, sizeof(double)*size_3d,        ptr_vn,             0, NULL, &ev_write_ptr_vn) );
	
#ifdef TEST
	cl_event ev_write_z_theta_v_fl_e;
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_theta_v_fl_e_d,       CL_FALSE, 0, sizeof(double)*size_3d,        z_theta_v_fl_e,     0, NULL, &ev_write_z_theta_v_fl_e) );

	double * ddqz_z_full_e_tmp = new double[size_3d];
	double * z_theta_v_e_tmp   = new double[size_3d];
	double * z_rho_e_tmp       = new double[size_3d];
	
	int * localStart_tmp = new int[nBlocks];
	int * localEnd_tmp = new int[nBlocks];
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
	for (int i=0; i<nBlocks; i++)
	{
		assert(localStart[i] == localStart_tmp[i]);
		assert(localEnd[i] == localEnd_tmp[i]);
	}
	
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddqz_z_full_e_d, CL_TRUE, 0, sizeof(double)*size_3d, ddqz_z_full_e_tmp, 0, NULL, NULL) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_theta_v_e_d, CL_TRUE, 0, sizeof(double)*size_3d, z_theta_v_e_tmp, 0, NULL, NULL) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_rho_e_d, CL_TRUE, 0, sizeof(double)*size_3d, z_rho_e_tmp, 0, NULL, NULL) );
	
	int zero[1] = {0};
	checkresultsCPP_(ddqz_z_full_e,ddqz_z_full_e_tmp,
					 zero,nproma,nproma,
					 zero,nlev,nlev,
					 zero,nblks_e,nblks_e,"fluxes:ddqz_z_full_e");
	
	
	checkresultsCPP_(z_theta_v_e,z_theta_v_e_tmp,
					 zero,nproma,nproma,
					 zero,nlev,nlev,
					 zero,nblks_e,nblks_e,"fluxes:z_theta_v_e");
	
	
	checkresultsCPP_(z_rho_e,z_rho_e_tmp,
					 zero,nproma,nproma,
					 zero,nlev,nlev,
					 zero,nblks_e,nblks_e,"fluxes:z_rho_e");
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("fluxes");
	}
#endif
	
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 4, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 5, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 6, sizeof(cl_mem), (void*)&z_rho_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 7, sizeof(cl_mem), (void*)&ptr_vn_d) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 8, sizeof(cl_mem), (void*)&ddqz_z_full_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes, 9, sizeof(cl_mem), (void*)&z_theta_v_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes,10, sizeof(cl_mem), (void*)&mass_fl_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_fluxes,11, sizeof(cl_mem), (void*)&z_theta_v_fl_e_d) );
	
	//cl_event ev_wait_list_fluxes[8] = { ev_write_localStart, ev_write_localEnd, ev_write_z_rho_e, ev_write_ptr_vn, ev_write_ddqz_z_full_e, ev_write_z_theta_v_e, ev_write_mass_fl_e, ev_write_z_theta_v_fl_e };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_fluxes, 3, NULL, global_fluxes, NULL, 8, ev_wait_list_fluxes, &ev_kernel_fluxes) );	
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_fluxes[4] = { ev_write_localStart, ev_write_localEnd, ev_write_ptr_vn, ev_write_mass_fl_e };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_fluxes, 1, NULL, global_fluxes, NULL, 4, ev_wait_list_fluxes, &ev_kernel_fluxes) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_fluxes[2] = { ev_write_ptr_vn, ev_write_mass_fl_e };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_fluxes, 1, NULL, global_fluxes, NULL, 2, ev_wait_list_fluxes, &ev_kernel_fluxes) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_fluxes, 1, NULL, global_fluxes, NULL, 0, NULL, &ev_kernel_fluxes) );
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
#ifndef TURBO_MODE
	if (*istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), mass_fl_e_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d, mass_fl_e,      0, NULL, &ev_read_mass_fl_e) );
#endif
	
#ifdef TEST
	cl_event ev_read_z_theta_v_fl_e;
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), mass_fl_e_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d, mass_fl_e,      0, NULL, &ev_read_mass_fl_e) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_theta_v_fl_e_d, CL_FALSE, 0, sizeof(cl_double)*size_3d, z_theta_v_fl_e, 0, NULL, &ev_read_z_theta_v_fl_e) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}


extern "C" void ocl_endphase_(int * i_startblk, int * i_endblk, int * nblks_c,
							  int * nlev, int * nlevp1, int * nrdmax,
							  int * localStart, int * localEnd, int * nproma,
							  int * ntl1, int * ntl2, int * istep,
							  int * itime_scheme,
							  double * dtime, double * cpd, double * cvd, double * rd, double * cvd_o_rd,
							  double * w_nnow,
							  double * w_nnew,
							  double * ddt_w_adv,
							  double * z_th_ddz_exner_c,
							  double * rho_ic,
							  double * w_concorr_c,
							  double * vwind_expl_wgt,
							  double * vwind_impl_wgt,
							  double * exner_nnow,
							  double * exner_nnew,
							  double * rhotheta_v_nnow,
							  double * rhotheta_v_nnew,
							  double * ddqz_z_full,
							  double * theta_v_h,
							  double * ddqz_z_half,
							  double * rho_nnow,
							  double * rho_nnew,
							  double * inv_ddqz_z_full,
							  double * z_mass_fl_div,
							  double * z_exner_pr,
							  double * z_theta_v_fl_div,
							  double * ddt_exner,
							  double * ddt_exner_phy,
							  double * rayleigh_w,
							  double * exner_ref_mc,
							  double * theta_v_nnew)
{
	// to check that the correct number of arguments is given (safety check)
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_endphase");
#endif
	
	const int i_startblkm1 = *i_startblk-1;
	const int nBlocks = *i_endblk-i_startblkm1;
	
	if (ocl.bFirstIter)
		for (int i=0; i<nBlocks; i++)
		{
			if (localStart[i]!=0 || localEnd[i]!=1)
			{
				cout << "need to use 3d indices\n";
				abort();
			}
		}
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_w_nnow, ev_write_ddt_w_adv1, ev_write_ddt_w_adv2, ev_write_ddt_w_adv3, ev_write_w_concorr_c;
	cl_event ev_write_ddqz_z_full, ev_write_vwind_impl_wgt;
	cl_event ev_write_z_mass_fl_div, ev_write_exner_pr, ev_write_z_theta_v_fl_div, ev_write_ddt_exner, ev_write_ddt_exner_phy, ev_write_w_nnew;
	cl_event ev_write_rayleigh_w;
	cl_event ev_write_exner_nnew, ev_write_rho_nnew, ev_write_rhotheta_v_nnow, ev_write_rhotheta_v_nnew, ev_write_theta_v_nnew;
	
	// to check
	cl_event ev_write_vwind_expl_wgt, ev_write_theta_v_h, ev_write_rho_ic, ev_write_z_th_ddz_exner_c, ev_write_exner_nnow, ev_write_ddqz_z_half, ev_write_z_exner_pr, ev_write_exner_ref_mc, ev_write_rho_nnow, ev_write_inv_ddqz_z_full;
	
	
	cl_event ev_kernel_endphase1, ev_kernel_endphase1c, ev_kernel_endphase2a, ev_kernel_endphase2b, ev_kernel_endphase3, ev_kernel_endphase4, ev_kernel_endphase5a, ev_kernel_endphase5b, ev_kernel_endphase5c, ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8;
	cl_event ev_kernel_endphase9a, ev_kernel_endphase9b, ev_kernel_endphase9c, ev_kernel_endphase9d;
	
	cl_event ev_read_w_nnew;
	cl_event ev_read_rho_nnew, ev_read_exner_nnew, ev_read_rhotheta_v_nnew, ev_read_theta_v_nnew;
	
	const size_t local_end_1c[3] = { 1, *nlev, 1 };
	const size_t local_end_2b[3] = { 1, *nlev, 1 };
	const size_t local_end_3[3] = { 1, *nlevp1, 1 };
	const size_t local_end_5[3] = { 1, *nlevp1, 1 };
	const size_t local_end_6[3] = { 1, *nlev,   1 };
	const size_t local_end_7[3] = { 1, 1 };
	const size_t local_end_9[3] = { 1, *nlevp1, 1 };
	
	const size_t global_end1[1] = { nBlocks*(*nlevp1)*(*nproma) }; // lvl 0,nlev is not important as it is not used in end6 (z_w_expl)
	const size_t global_end1c[3] = { nBlocks, *nlev, *nproma };
	const size_t global_end2a[1] = { nBlocks*(*nlev)*(*nproma) };
	const size_t global_end2b[3] = { nBlocks, *nlev, *nproma };
	const size_t global_end3[3] = { nBlocks, *nlevp1, *nproma };
	const size_t global_end4[2] = { nBlocks, *nproma };
	const size_t global_end5a[3] = { nBlocks, 2, *nproma };
	const size_t global_end5b[3] = { nBlocks, *nlevp1, *nproma };
	const size_t global_end5c[1] = { nBlocks*(*nlev)*(*nproma) };
	const size_t global_end6[3] = { nBlocks, *nlev, *nproma };
	const size_t global_end7[2] = { nBlocks, *nproma };
	const size_t global_end8[3] = { nBlocks, *nrdmax-1, *nproma };
	const size_t global_end9[3] = { nBlocks, *nlev, *nproma };
	const size_t global_end9a[3] = { nBlocks, *nlevp1, *nproma };
	const size_t global_end9c[1] = { nBlocks*(*nlev)*(*nproma) };
	
	const int size_localStart = nBlocks;
	const int size_localEnd   = nBlocks;
	const int size_2d_no_jk   = (*nproma)*(*nblks_c);
	const int size_3d         = (*nproma)*(*nlev)*(*nblks_c);
	const int size_3d_p1      = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_4d_p1_3    = (*nproma)*(*nlevp1)*(*nblks_c)*3;
	
	cl_mem localStart_d           = resourcesCL.getBuffers("localStart5",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d             = resourcesCL.getBuffers("localEnd5",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem w_nnow_d               = resourcesCL.getBuffers("w_nnow",               1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem w_nnew_d               = resourcesCL.getBuffers("w_nnew",               1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem ddt_w_adv1_d           = resourcesCL.getBuffers("ddt_w_adv1",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem ddt_w_adv2_d           = resourcesCL.getBuffers("ddt_w_adv2",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem ddt_w_adv3_d           = resourcesCL.getBuffers("ddt_w_adv3",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_th_ddz_exner_c_d     = resourcesCL.getBuffers("z_th_ddz_exner_c",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem rho_ic_d               = resourcesCL.getBuffers("rho_ic",               1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem w_concorr_c_d          = resourcesCL.getBuffers("w_concorr_c",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem vwind_expl_wgt_d       = resourcesCL.getBuffers("vwind_expl_wgt",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_2d_no_jk)[0];
	cl_mem vwind_impl_wgt_d       = resourcesCL.getBuffers("vwind_impl_wgt",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_2d_no_jk)[0];
	cl_mem z_w_expl_d             = resourcesCL.getBuffers("z_w_expl",             1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_contr_w_fl_l_d       = resourcesCL.getBuffers("z_contr_w_fl_l",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem exner_nnow_d           = resourcesCL.getBuffers("exner_nnow",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem exner_nnew_d           = resourcesCL.getBuffers("exner_nnew",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem rhotheta_v_nnow_d      = resourcesCL.getBuffers("rhotheta_v_nnow",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem rhotheta_v_nnew_d      = resourcesCL.getBuffers("rhotheta_v_nnew",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem ddqz_z_full_d          = resourcesCL.getBuffers("ddqz_z_full",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem theta_v_h_d            = resourcesCL.getBuffers("theta_v_h",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_beta_d               = resourcesCL.getBuffers("z_beta",               1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_alpha_d              = resourcesCL.getBuffers("z_alpha",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem ddqz_z_half_d          = resourcesCL.getBuffers("ddqz_z_half",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_gamma_d              = resourcesCL.getBuffers("z_gamma",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_a_d                  = resourcesCL.getBuffers("z_a",                  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_b_d                  = resourcesCL.getBuffers("z_b",                  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_c_d                  = resourcesCL.getBuffers("z_c",                  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_q_d                  = resourcesCL.getBuffers("z_q",                  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem z_g_d                  = resourcesCL.getBuffers("z_g",                  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem rho_nnow_d             = resourcesCL.getBuffers("rho_nnow",             1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem rho_nnew_d             = resourcesCL.getBuffers("rho_nnew",             1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem inv_ddqz_z_full_d      = resourcesCL.getBuffers("inv_ddqz_z_full",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_mass_fl_div_d        = resourcesCL.getBuffers("z_mass_fl_div",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_exner_pr_d           = resourcesCL.getBuffers("z_exner_pr",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_theta_v_fl_div_d     = resourcesCL.getBuffers("z_theta_v_fl_div",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem ddt_exner_d            = resourcesCL.getBuffers("ddt_exner",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem ddt_exner_phy_d        = resourcesCL.getBuffers("ddt_exner_phy",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_rho_expl_d           = resourcesCL.getBuffers("z_rho_expl",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_exner_expl_d         = resourcesCL.getBuffers("z_exner_expl",         1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem rayleigh_w_d           = resourcesCL.getBuffers("rayleigh_w",           1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d_p1)[0];
	cl_mem exner_ref_mc_d         = resourcesCL.getBuffers("exner_ref_mc",         1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem theta_v_nnew_d         = resourcesCL.getBuffers("theta_v_nnew",         1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_contr_w_fl_l_diff_d  = resourcesCL.getBuffers("z_contr_w_fl_l_diff",  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	cl_mem z_th_contr_diff_d      = resourcesCL.getBuffers("z_th_contr_diff",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_3d)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		// futile transfers
		/*{
			CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vwind_expl_wgt_d,       CL_FALSE, 0, sizeof(double)*size_2d_no_jk, vwind_expl_wgt,     0, NULL, &ev_write_vwind_expl_wgt) );
			CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddqz_z_half_d,          CL_FALSE, 0, sizeof(double)*size_3d_p1,    ddqz_z_half,        0, NULL, &ev_write_ddqz_z_half) );
			CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_full_d,      CL_FALSE, 0, sizeof(double)*size_3d,       inv_ddqz_z_full,    0, NULL, &ev_write_inv_ddqz_z_full) );
			CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), exner_ref_mc_d,         CL_FALSE, 0, sizeof(double)*size_3d,       exner_ref_mc,       0, NULL, &ev_write_exner_ref_mc) );
		}*/
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vwind_impl_wgt_d,       CL_FALSE, 0, sizeof(double)*size_2d_no_jk, vwind_impl_wgt,     0, NULL, &ev_write_vwind_impl_wgt) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddqz_z_full_d,          CL_FALSE, 0, sizeof(double)*size_3d,       ddqz_z_full,        0, NULL, &ev_write_ddqz_z_full) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rayleigh_w_d,           CL_FALSE, 0, sizeof(double)*size_3d_p1,    rayleigh_w,         0, NULL, &ev_write_rayleigh_w) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		// in theory, these transfers are futile
		/*{
		 CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_v_h_d,            CL_FALSE, 0, sizeof(double)*size_3d_p1,    theta_v_h,          0, NULL, &ev_write_theta_v_h) );
		 CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rho_nnow_d,             CL_FALSE, 0, sizeof(double)*size_3d,       rho_nnow,           0, NULL, &ev_write_rho_nnow) );
		 CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), exner_nnow_d,           CL_FALSE, 0, sizeof(double)*size_3d,       exner_nnow,         0, NULL, &ev_write_exner_nnow) );
		 CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rho_ic_d,               CL_FALSE, 0, sizeof(double)*size_3d_p1,    rho_ic,             0, NULL, &ev_write_rho_ic) );
		 }*/
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,           CL_FALSE, 0, sizeof(int)*size_localStart,  localStart,         0, NULL, &ev_write_localStart) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,             CL_FALSE, 0, sizeof(int)*size_localEnd,    localEnd,           0, NULL, &ev_write_localEnd) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_nnow_d,               CL_FALSE, 0, sizeof(double)*size_3d_p1,    w_nnow,             0, NULL, &ev_write_w_nnow) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_w_adv1_d,           CL_FALSE, 0, sizeof(double)*size_3d_p1,  ddt_w_adv,                0, NULL, &ev_write_ddt_w_adv1) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_w_adv2_d,           CL_FALSE, 0, sizeof(double)*size_3d_p1,  &ddt_w_adv[size_3d_p1],   0, NULL, &ev_write_ddt_w_adv2) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_w_adv3_d,           CL_FALSE, 0, sizeof(double)*size_3d_p1,  &ddt_w_adv[size_3d_p1*2], 0, NULL, &ev_write_ddt_w_adv3) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_concorr_c_d,          CL_FALSE, 0, sizeof(double)*size_3d_p1,    w_concorr_c,        0, NULL, &ev_write_w_concorr_c) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rhotheta_v_nnow_d,      CL_FALSE, 0, sizeof(double)*size_3d,       rhotheta_v_nnow,    0, NULL, &ev_write_rhotheta_v_nnow) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rhotheta_v_nnew_d,      CL_FALSE, 0, sizeof(double)*size_3d,       rhotheta_v_nnew,    0, NULL, &ev_write_rhotheta_v_nnew) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_mass_fl_div_d,        CL_FALSE, 0, sizeof(double)*size_3d,       z_mass_fl_div,      0, NULL, &ev_write_z_mass_fl_div) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_theta_v_fl_div_d,     CL_FALSE, 0, sizeof(double)*size_3d,       z_theta_v_fl_div,   0, NULL, &ev_write_z_theta_v_fl_div) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_exner_d,            CL_FALSE, 0, sizeof(double)*size_3d,       ddt_exner,          0, NULL, &ev_write_ddt_exner) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_exner_phy_d,        CL_FALSE, 0, sizeof(double)*size_3d,       ddt_exner_phy,      0, NULL, &ev_write_ddt_exner_phy) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_nnew_d,               CL_FALSE, 0, sizeof(double)*size_3d_p1,    w_nnew,             0, NULL, &ev_write_w_nnew) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), exner_nnew_d,           CL_FALSE, 0, sizeof(double)*size_3d,       exner_nnew,         0, NULL, &ev_write_exner_nnew) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rho_nnew_d,             CL_FALSE, 0, sizeof(double)*size_3d,       rho_nnew,           0, NULL, &ev_write_rho_nnew) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), theta_v_nnew_d,         CL_FALSE, 0, sizeof(double)*size_3d,       theta_v_nnew,       0, NULL, &ev_write_theta_v_nnew) );
	}
	
#ifdef TEST
	//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_th_ddz_exner_c_d,     CL_FALSE, 0, sizeof(double)*size_3d,       z_th_ddz_exner_c,   0, NULL, NULL) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_exner_pr_d,           CL_FALSE, 0, sizeof(double)*size_3d,       z_exner_pr,         0, NULL, NULL) );
	{
		cout << "checking transfers for ocl_endphase...\n";
		double * vwind_impl_wgt_tmp  = new double[size_2d_no_jk];
		double * ddqz_z_full_tmp     = new double[size_3d];
		double * inv_ddqz_z_full_tmp = new double[size_3d];
		
		double * exner_ref_mc_tmp   = new double[size_3d];
		double * rho_ic_tmp         = new double[size_3d_p1];
		double * exner_nnow_tmp     = new double[size_3d];
		double * theta_v_h_tmp      = new double[size_3d_p1];
		double * rho_nnow_tmp       = new double[size_3d];
		
		double * z_mass_fl_div_tmp    = new double[size_3d];
		double * z_theta_v_fl_div_tmp = new double[size_3d];
		double * w_concorr_c_tmp      = new double[size_3d_p1];
		
		double * ddt_w_adv_tmp = new double[size_4d_p1_3];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		int zero[1] = {0};
		int one[1] = {0};
		int three[1] = {0};
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vwind_impl_wgt_d,  CL_FALSE, 0, sizeof(cl_double)*size_2d_no_jk, vwind_impl_wgt_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddqz_z_full_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d,       ddqz_z_full_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_full_d, CL_FALSE, 0, sizeof(cl_double)*size_3d,       inv_ddqz_z_full_tmp, 0, NULL, NULL) );

		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_ref_mc_d,   CL_FALSE, 0, sizeof(cl_double)*size_3d,    exner_ref_mc_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_ic_d,         CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, rho_ic_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_nnow_d,     CL_FALSE, 0, sizeof(cl_double)*size_3d,    exner_nnow_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_h_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, theta_v_h_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_nnow_d,       CL_FALSE, 0, sizeof(cl_double)*size_3d,    rho_nnow_tmp, 0, NULL, NULL) );

		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_mass_fl_div_d,    CL_FALSE, 0, sizeof(cl_double)*size_3d,    z_mass_fl_div_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_theta_v_fl_div_d, CL_FALSE, 0, sizeof(cl_double)*size_3d,    z_theta_v_fl_div_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_concorr_c_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, w_concorr_c_tmp, 0, NULL, NULL) );
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv1_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, ddt_w_adv_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv2_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, &ddt_w_adv_tmp[size_3d_p1], 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv3_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, &ddt_w_adv_tmp[size_3d_p1*2], 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		checkresultsCPP_(vwind_impl_wgt,vwind_impl_wgt_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,one,one,"endphase:vwind_impl_wgt");
		
		checkresultsCPP_(ddqz_z_full,ddqz_z_full_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"endphase:ddqz_z_full");
		
		checkresultsCPP_(inv_ddqz_z_full,inv_ddqz_z_full_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"endphase:inv_ddqz_z_full");
		
		checkresultsCPP_(exner_ref_mc,exner_ref_mc_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"endphase:exner_ref_mc");
		
		checkresultsCPP_(rho_ic,rho_ic_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"endphase:rho_ic");
		
		checkresultsCPP_(exner_nnow,exner_nnow_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"endphase:exner_nnow");
		
		checkresultsCPP_(theta_v_h,theta_v_h_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"endphase:theta_v_h");
		
		checkresultsCPP_(rho_nnow,rho_nnow_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"endphase:rho_nnow");
		
		checkresultsCPP_(z_mass_fl_div,z_mass_fl_div_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"endphase:z_mass_fl_div");
		
		checkresultsCPP_(z_theta_v_fl_div,z_theta_v_fl_div_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"endphase:z_theta_v_fl_div");
		
		checkresultsCPP_(w_concorr_c,w_concorr_c_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"endphase:w_concorr_c");
		
		checkresults4dCPP_(ddt_w_adv,ddt_w_adv_tmp,
						   zero,nproma,nproma,
						   zero,nlevp1,nlevp1,
						   zero,nblks_c,nblks_c,
						   zero,three,three,"endphase:ddt_w_adv");
	}
#endif
	
#ifndef REGION_TIMING	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
	}
#endif
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.push_start("endphase1ab");
	}
	
	if (*istep==2 && *itime_scheme==6)
	{
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 2, sizeof(int32_t), (void*)nblks_c) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 3, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 4, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 5, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 6, sizeof(double), (void*)dtime) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 7, sizeof(double), (void*)cpd) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 8, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a, 9, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a,10, sizeof(cl_mem), (void*)&w_nnow_d) );
		if (*ntl1==1)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1a,11, sizeof(cl_mem), (void*)&ddt_w_adv1_d) );
		}
		else if (*ntl1==2)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1a,11, sizeof(cl_mem), (void*)&ddt_w_adv2_d) );
		}
		else if (*ntl1==3)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1a,11, sizeof(cl_mem), (void*)&ddt_w_adv3_d) );
		}
		else abort();
		if (*ntl2==1)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1a,12, sizeof(cl_mem), (void*)&ddt_w_adv1_d) );
		}
		else if (*ntl2==2)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1a,12, sizeof(cl_mem), (void*)&ddt_w_adv2_d) );
		}
		else if (*ntl2==3)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1a,12, sizeof(cl_mem), (void*)&ddt_w_adv3_d) );
		}
		else abort();
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a,13, sizeof(cl_mem), (void*)&z_th_ddz_exner_c_d) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1a,14, sizeof(cl_mem), (void*)&z_w_expl_d) );
		
#ifndef TURBO_MODE
		if (*istep==1)
		{
			cl_event ev_wait_list_end1[4] = { ev_write_w_nnow, ev_write_ddt_w_adv1, ev_write_ddt_w_adv2, ev_write_ddt_w_adv3 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase1a, 1, NULL, global_end1, NULL, 4, ev_wait_list_end1, &ev_kernel_endphase1) );
		}
		else
#endif
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase1a, 1, NULL, global_end1, NULL, 0, NULL, &ev_kernel_endphase1) );
	}
	else
	{
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 2, sizeof(int32_t), (void*)nblks_c) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 3, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 4, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 5, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 6, sizeof(double), (void*)dtime) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 7, sizeof(double), (void*)cpd) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 8, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b, 9, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b,10, sizeof(cl_mem), (void*)&w_nnow_d) );
		if (*ntl1==1)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1b,11, sizeof(cl_mem), (void*)&ddt_w_adv1_d) );
		}
		else if (*ntl1==2)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1b,11, sizeof(cl_mem), (void*)&ddt_w_adv2_d) );
		}
		else if (*ntl1==3)
		{
			CheckCL( clSetKernelArg(ocl.kernel_endphase1b,11, sizeof(cl_mem), (void*)&ddt_w_adv3_d) );
		}
		else abort();
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b,12, sizeof(cl_mem), (void*)&z_th_ddz_exner_c_d) );
		CheckCL( clSetKernelArg(ocl.kernel_endphase1b,13, sizeof(cl_mem), (void*)&z_w_expl_d) );
		
#ifndef TURBO_MODE
		if (*istep==1)
		{
			cl_event ev_wait_list_end1[4] = { ev_write_w_nnow, ev_write_ddt_w_adv1, ev_write_ddt_w_adv2, ev_write_ddt_w_adv3 };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase1b, 1, NULL, global_end1, NULL, 4, ev_wait_list_end1, &ev_kernel_endphase1) );
		}
		else
#endif
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase1b, 1, NULL, global_end1, NULL, 0, NULL, &ev_kernel_endphase1) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase1c");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 2, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 4, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 5, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 8, sizeof(cl_mem), (void*)&w_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c, 9, sizeof(cl_mem), (void*)&rho_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c,10, sizeof(cl_mem), (void*)&w_concorr_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c,11, sizeof(cl_mem), (void*)&vwind_expl_wgt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase1c,12, sizeof(cl_mem), (void*)&z_contr_w_fl_l_d) );
	
#ifndef TURBO_MODE
	if (*istep==1)
	{
		cl_event ev_wait_list_end1c[1] = { ev_write_w_nnow };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase1c, 3, NULL, global_end1c, local_end_1c, 1, ev_wait_list_end1c, &ev_kernel_endphase1c) );
	}
	else
#endif
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase1c, 3, NULL, global_end1c, local_end_1c, 0, NULL, &ev_kernel_endphase1c) );
	
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase2a");
	}
	
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 4, sizeof(double), (void*)dtime) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 5, sizeof(double), (void*)rd) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 6, sizeof(double), (void*)cvd) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 7, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 8, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a, 9, sizeof(cl_mem), (void*)&exner_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a,10, sizeof(cl_mem), (void*)&rhotheta_v_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a,11, sizeof(cl_mem), (void*)&ddqz_z_full_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2a,12, sizeof(cl_mem), (void*)&z_beta_d) );
	
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_end2a[2] = { ev_write_rhotheta_v_nnow, ev_write_ddqz_z_full };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase2a, 1, NULL, global_end2a, NULL, 2, ev_wait_list_end2a, &ev_kernel_endphase2a) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_end2a[1] = { ev_write_rhotheta_v_nnow };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase2a, 1, NULL, global_end2a, NULL, 1, ev_wait_list_end2a, &ev_kernel_endphase2a) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase2a, 1, NULL, global_end2a, NULL, 0, NULL, &ev_kernel_endphase2a) );
	
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase2b");
	}
	
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 7, sizeof(cl_mem), (void*)&vwind_impl_wgt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 8, sizeof(cl_mem), (void*)&theta_v_h_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b, 9, sizeof(cl_mem), (void*)&rho_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase2b,10, sizeof(cl_mem), (void*)&z_alpha_d) );
	
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_end2b[1] = { ev_write_vwind_impl_wgt };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase2b, 3, NULL, global_end2b, local_end_2b, 2, ev_wait_list_end2b, &ev_kernel_endphase2b) );
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase2b, 3, NULL, global_end2b, local_end_2b, 0, NULL, &ev_kernel_endphase2b) );
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase3");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 5, sizeof(double), (void*)dtime) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 6, sizeof(double), (void*)cpd) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 7, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 8, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3, 9, sizeof(cl_mem), (void*)&vwind_impl_wgt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,10, sizeof(cl_mem), (void*)&theta_v_h_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,11, sizeof(cl_mem), (void*)&ddqz_z_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,12, sizeof(cl_mem), (void*)&z_alpha_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,13, sizeof(cl_mem), (void*)&z_beta_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,14, sizeof(cl_mem), (void*)&z_gamma_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,15, sizeof(cl_mem), (void*)&z_a_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,16, sizeof(cl_mem), (void*)&z_b_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,17, sizeof(cl_mem), (void*)&z_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase3,18, sizeof(cl_mem), (void*)&z_q_d) );
	
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_end3[3] = { ev_write_vwind_impl_wgt, ev_kernel_endphase2a, ev_kernel_endphase2b };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase3, 3, NULL, global_end3, local_end_3, 3, ev_wait_list_end3, &ev_kernel_endphase3) );
	}
	else
	{
		cl_event ev_wait_list_end3[2] = { ev_kernel_endphase2a, ev_kernel_endphase2b };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase3, 3, NULL, global_end3, local_end_3, 2, ev_wait_list_end3, &ev_kernel_endphase3) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase4");
	}
	 
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 7, sizeof(cl_mem), (void*)&z_a_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 8, sizeof(cl_mem), (void*)&z_b_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4, 9, sizeof(cl_mem), (void*)&z_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4,10, sizeof(cl_mem), (void*)&z_g_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase4,11, sizeof(cl_mem), (void*)&z_q_d) );
	
	cl_event ev_wait_list_end4[1] = { ev_kernel_endphase3 };
	CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase4, 2, NULL, global_end4, NULL, 1, ev_wait_list_end4, &ev_kernel_endphase4) );
	
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase5a");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 2, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 4, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 5, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 6, sizeof(cl_mem), (void*)&w_concorr_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 7, sizeof(cl_mem), (void*)&z_contr_w_fl_l_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5a, 8, sizeof(cl_mem), (void*)&w_nnew_d) );

#ifndef TURBO_MODE
	if (*istep==1)
	{
		cl_event ev_wait_list_end5a[2] = { ev_write_w_nnew, ev_kernel_endphase1c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase5a, 3, NULL, global_end5a, NULL, 2, ev_wait_list_end5a, &ev_kernel_endphase5a) );
	}
	else
#endif
	{
		cl_event ev_wait_list_end5a[1] = { ev_kernel_endphase1c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase5a, 3, NULL, global_end5a, NULL, 1, ev_wait_list_end5a, &ev_kernel_endphase5a) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase5b");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 5, sizeof(double), (void*)dtime) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 8, sizeof(cl_mem), (void*)&z_contr_w_fl_l_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b, 9, sizeof(cl_mem), (void*)&theta_v_h_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b,10, sizeof(cl_mem), (void*)&z_contr_w_fl_l_diff_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5b,11, sizeof(cl_mem), (void*)&z_th_contr_diff_d) );
	
#ifndef TURBO_MODE
	if (*istep==1)
	{
		cl_event ev_wait_list_end5b[3] = { ev_write_ddt_exner, ev_write_ddt_exner_phy, ev_kernel_endphase1c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase5b, 3, NULL, global_end5b, local_end_5, 3, ev_wait_list_end5b, &ev_kernel_endphase5b) );
	}
	else
#endif
	{
		cl_event ev_wait_list_end5b[1] = { ev_kernel_endphase1c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase5b, 3, NULL, global_end5b, local_end_5, 1, ev_wait_list_end5b, &ev_kernel_endphase5b) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase5c");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 4, sizeof(double), (void*)dtime) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 7, sizeof(cl_mem), (void*)&rho_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 8, sizeof(cl_mem), (void*)&inv_ddqz_z_full_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c, 9, sizeof(cl_mem), (void*)&z_mass_fl_div_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,10, sizeof(cl_mem), (void*)&z_contr_w_fl_l_diff_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,11, sizeof(cl_mem), (void*)&z_exner_pr_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,12, sizeof(cl_mem), (void*)&z_beta_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,13, sizeof(cl_mem), (void*)&z_theta_v_fl_div_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,14, sizeof(cl_mem), (void*)&z_th_contr_diff_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,15, sizeof(cl_mem), (void*)&ddt_exner_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,16, sizeof(cl_mem), (void*)&ddt_exner_phy_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,17, sizeof(cl_mem), (void*)&z_rho_expl_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase5c,18, sizeof(cl_mem), (void*)&z_exner_expl_d) );
	
#ifndef TURBO_MODE
	if (*istep==1)
	{
		cl_event ev_wait_list_end5c[4] = { ev_write_ddt_exner, ev_write_ddt_exner_phy, ev_kernel_endphase2a, ev_kernel_endphase5b };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase5c, 1, NULL, global_end5c, NULL, 4, ev_wait_list_end5c, &ev_kernel_endphase5c) );
	}
	else
#endif
	{
		cl_event ev_wait_list_end5c[2] = { ev_kernel_endphase2a, ev_kernel_endphase5b };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase5c, 1, NULL, global_end5c, NULL, 2, ev_wait_list_end5c, &ev_kernel_endphase5c) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase6");
	}

	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 7, sizeof(cl_mem), (void*)&z_w_expl_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 8, sizeof(cl_mem), (void*)&z_gamma_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6, 9, sizeof(cl_mem), (void*)&z_exner_expl_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6,10, sizeof(cl_mem), (void*)&z_b_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase6,11, sizeof(cl_mem), (void*)&w_nnew_d) );

#ifndef TURBO_MODE
	if (*istep==1)
	{
		cl_event ev_wait_list_end6[5] = { ev_write_w_nnew, ev_kernel_endphase1, ev_kernel_endphase3, ev_kernel_endphase5a, ev_kernel_endphase5c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase6, 3, NULL, global_end6, local_end_6, 5, ev_wait_list_end6, &ev_kernel_endphase6) );
	}
	else
#endif
	{
		cl_event ev_wait_list_end6[4] = { ev_kernel_endphase1, ev_kernel_endphase3, ev_kernel_endphase5a, ev_kernel_endphase5c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase6, 3, NULL, global_end6, local_end_6, 4, ev_wait_list_end6, &ev_kernel_endphase6) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase7");
	}
	 
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 7, sizeof(cl_mem), (void*)&z_a_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 8, sizeof(cl_mem), (void*)&z_g_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7, 9, sizeof(cl_mem), (void*)&z_q_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase7,10, sizeof(cl_mem), (void*)&w_nnew_d) );
	
#ifndef TURBO_MODE
	if (*istep==1)
	{
		cl_event ev_wait_list_end7[5] = { ev_write_w_nnew, ev_kernel_endphase3, ev_kernel_endphase4, ev_kernel_endphase5a, ev_kernel_endphase6 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase7, 2, NULL, global_end7, local_end_7, 5, ev_wait_list_end7, &ev_kernel_endphase7) );
	}
	else
#endif
	{
		cl_event ev_wait_list_end7[4] = { ev_kernel_endphase3, ev_kernel_endphase4, ev_kernel_endphase5a, ev_kernel_endphase6 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase7, 2, NULL, global_end7, local_end_7, 4, ev_wait_list_end7, &ev_kernel_endphase7) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase8");
	}
	 
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 2, sizeof(int32_t), (void*)nrdmax) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 5, sizeof(double), (void*)dtime) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 8, sizeof(cl_mem), (void*)&rayleigh_w_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase8, 9, sizeof(cl_mem), (void*)&w_nnew_d) );

	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_end8[5] = { ev_write_rayleigh_w, ev_write_w_nnew, ev_kernel_endphase5a, ev_kernel_endphase6, ev_kernel_endphase7 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase8, 3, NULL, global_end8, NULL, 5, ev_wait_list_end8, &ev_kernel_endphase8) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_end8[4] = { ev_write_w_nnew, ev_kernel_endphase5a, ev_kernel_endphase6, ev_kernel_endphase7 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase8, 3, NULL, global_end8, NULL, 4, ev_wait_list_end8, &ev_kernel_endphase8) );
#endif
	}
	else
	{
		cl_event ev_wait_list_end8[3] = { ev_kernel_endphase5a, ev_kernel_endphase6, ev_kernel_endphase7 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase8, 3, NULL, global_end8, NULL, 3, ev_wait_list_end8, &ev_kernel_endphase8) );
	}

	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase9a");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 5, sizeof(double), (void*)dtime) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 8, sizeof(cl_mem), (void*)&z_rho_expl_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a, 9, sizeof(cl_mem), (void*)&vwind_impl_wgt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a,10, sizeof(cl_mem), (void*)&inv_ddqz_z_full_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a,11, sizeof(cl_mem), (void*)&rho_ic_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a,12, sizeof(cl_mem), (void*)&w_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9a,13, sizeof(cl_mem), (void*)&rho_nnew_d) );
	
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_end9a[8] = { ev_write_vwind_impl_wgt, ev_write_w_nnew, ev_write_rho_nnew,
			                               ev_kernel_endphase5a, ev_kernel_endphase5c, ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9a, 3, NULL, global_end9a, local_end_9, 8, ev_wait_list_end9a, &ev_kernel_endphase9a) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_end9a[7] = { ev_write_w_nnew, ev_write_rho_nnew,
			                               ev_kernel_endphase5a, ev_kernel_endphase5c, ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9a, 3, NULL, global_end9a, local_end_9, 7, ev_wait_list_end9a, &ev_kernel_endphase9a) );
#endif
	}
	else
	{
		cl_event ev_wait_list_end9a[5] = { ev_kernel_endphase5a, ev_kernel_endphase5c, ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9a, 3, NULL, global_end9a, local_end_9, 5, ev_wait_list_end9a, &ev_kernel_endphase9a) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase9b");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 7, sizeof(cl_mem), (void*)&w_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 8, sizeof(cl_mem), (void*)&z_exner_expl_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b, 9, sizeof(cl_mem), (void*)&exner_ref_mc_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b,10, sizeof(cl_mem), (void*)&z_beta_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b,11, sizeof(cl_mem), (void*)&z_alpha_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9b,12, sizeof(cl_mem), (void*)&exner_nnew_d) );
	
#ifndef TURBO_MODE
	if (ocl.bFirstIter || *istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		cl_event ev_wait_list_end9b[7] = { ev_write_w_nnew, ev_write_exner_nnew,
			                               ev_kernel_endphase5a, ev_kernel_endphase5c, ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9b, 3, NULL, global_end9a, local_end_9, 7, ev_wait_list_end9b, &ev_kernel_endphase9b) );
	}
	else
	{
		cl_event ev_wait_list_end9b[5] = { ev_kernel_endphase5a, ev_kernel_endphase5c, ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9b, 3, NULL, global_end9a, local_end_9, 5, ev_wait_list_end9b, &ev_kernel_endphase9b) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase9c");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 4, sizeof(double), (void*)cvd_o_rd) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 7, sizeof(cl_mem), (void*)&exner_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 8, sizeof(cl_mem), (void*)&exner_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c, 9, sizeof(cl_mem), (void*)&rhotheta_v_nnow_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9c,10, sizeof(cl_mem), (void*)&rhotheta_v_nnew_d) );
	
#ifndef TURBO_MODE
	if (ocl.bFirstIter || *istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		cl_event ev_wait_list_end9c[7] = { ev_write_exner_nnew, ev_write_rhotheta_v_nnow, ev_write_rhotheta_v_nnew,
			                               ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8, ev_kernel_endphase9b };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9c, 1, NULL, global_end9c, NULL, 7, ev_wait_list_end9c, &ev_kernel_endphase9c) );
	}
	else
	{
		cl_event ev_wait_list_end9c[4] = { ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8, ev_kernel_endphase9b };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9c, 1, NULL, global_end9c, NULL, 4, ev_wait_list_end9c, &ev_kernel_endphase9c) );
	}
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("endphase9d");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 3, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 4, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 5, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 6, sizeof(cl_mem), (void*)&rho_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 7, sizeof(cl_mem), (void*)&rhotheta_v_nnew_d) );
	CheckCL( clSetKernelArg(ocl.kernel_endphase9d, 8, sizeof(cl_mem), (void*)&theta_v_nnew_d) );
	
#ifndef TURBO_MODE
	if (ocl.bFirstIter || *istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		cl_event ev_wait_list_end9d[8] = { ev_write_rho_nnew, ev_write_rhotheta_v_nnew, ev_write_theta_v_nnew,
			                                ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8, ev_kernel_endphase9a, ev_kernel_endphase9c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9d, 1, NULL, global_end9c, NULL, 8, ev_wait_list_end9d, &ev_kernel_endphase9d) );
	}
	else
	{
		cl_event ev_wait_list_end9d[5] = { ev_kernel_endphase6, ev_kernel_endphase7, ev_kernel_endphase8, ev_kernel_endphase9a, ev_kernel_endphase9c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_endphase9d, 1, NULL, global_end9c, NULL, 5, ev_wait_list_end9d, &ev_kernel_endphase9d) );
	}
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (*istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_nnew_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d,    exner_nnew,      0, NULL, &ev_read_exner_nnew) );
	
#ifndef TURBO_MODE
	if (*istep!=1)
	{
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_nnew_d,          CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, w_nnew,          0, NULL, &ev_read_w_nnew) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_nnew_d,        CL_FALSE, 0, sizeof(cl_double)*size_3d,    rho_nnew,        0, NULL, &ev_read_rho_nnew) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rhotheta_v_nnew_d, CL_FALSE, 0, sizeof(cl_double)*size_3d,    rhotheta_v_nnew, 0, NULL, &ev_read_rhotheta_v_nnew) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_nnew_d,    CL_FALSE, 0, sizeof(cl_double)*size_3d,    theta_v_nnew,    0, NULL, &ev_read_theta_v_nnew) );
	}
#endif

#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_nnew_d,          CL_FALSE, 0, sizeof(cl_double)*size_3d_p1, w_nnew,          0, NULL, &ev_read_w_nnew) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rho_nnew_d,        CL_FALSE, 0, sizeof(cl_double)*size_3d,    rho_nnew,        0, NULL, &ev_read_rho_nnew) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), exner_nnew_d,      CL_FALSE, 0, sizeof(cl_double)*size_3d,    exner_nnew,      0, NULL, &ev_read_exner_nnew) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rhotheta_v_nnew_d, CL_FALSE, 0, sizeof(cl_double)*size_3d,    rhotheta_v_nnew, 0, NULL, &ev_read_rhotheta_v_nnew) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), theta_v_nnew_d,    CL_FALSE, 0, sizeof(cl_double)*size_3d,    theta_v_nnew,    0, NULL, &ev_read_theta_v_nnew) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	ocl.bFirstIter = false;
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}


#pragma mark Velocity Tendencies
extern "C" void ocl_part1_(int * i_startblk, int * i_endblk,
						   int * nblks_e, int * nblks_c,
						   int * nlev, int * nlevp1,
						   int * localStart, int * localEnd,
						   int * nproma,
						   int * istep,
						   int * nflat,
						   bool * l_vert_nested,
						   double * wgtfac_e, // start of inputs
						   double * wgtfacq1_e,
						   double * wgtfacq_e,
						   double * vn, int * nameIdx,
						   double * vt,
						   double * ddxn_z_half,
						   double * ddxt_z_half,
						   double * vn_half, // start of outputs
						   double * z_kin_hor_e,
						   double * z_concorr_e)
{	
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_part1");
#endif
	
	const int i_startblkm1 = *i_startblk - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_wgtfac_e, ev_write_wgtfacq1_e, ev_write_wgtfacq_e, ev_write_vn, ev_write_vt, ev_write_ddxn_z_half, ev_write_ddxt_z_half;
	cl_event ev_kernel_step1ia, ev_kernel_step1ib, ev_kernel_step1ic, ev_kernel_step1a, ev_kernel_step1b, ev_kernel_stepn;
	cl_event ev_read_vn_half, ev_read_z_kin_hor_e, ev_read_z_concorr_e;
	cl_event ev_write_vn_half/*, ev_write_z_kin_hor_e*/;

	const int size_wgtfac_e    = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_wgtfacq1_e  = (*nproma)*3*(*nblks_e);
	const int size_wgtfacq_e   = (*nproma)*3*(*nblks_e);
	const int size_vn          = (*nproma)*(*nlev)*(*nblks_e);
	const int size_vt          = (*nproma)*(*nlev)*(*nblks_e);
	const int size_ddxn_z_half = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_ddxt_z_half = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_vn_half     = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_z_kin_hor_e = (*nproma)*(*nlev)*(*nblks_e);
	const int size_z_vt_ie     = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_z_concorr_e = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_localStart  = nBlocks;
	const int size_localEnd    = nBlocks;
	
	size_t local_1a[3] = {1, *nlev, 1 };
	
	size_t global_1ia[3] = { nBlocks, *nlev-*nflat, *nproma };
	size_t global_1ib[2] = { nBlocks, *nproma };
	size_t global_1ic[3] = { nBlocks, *nlevp1-*nflat, *nproma };
	size_t global_1a[3]  = { nBlocks, *nlev, *nproma };
	size_t global_1b[2]  = { nBlocks, *nproma };
	size_t global_n[1]   = { nBlocks*(*nlev)*(*nproma) };
	
	cl_mem vn_d          = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vn)[0];
	cl_mem wgtfac_e_d    = resourcesCL.getBuffers("wgtfac_e",    1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_wgtfac_e)[0];
	cl_mem wgtfacq1_e_d  = resourcesCL.getBuffers("wgtfacq1_e",  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_wgtfacq1_e)[0];
	cl_mem wgtfacq_e_d   = resourcesCL.getBuffers("wgtfacq_e",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_wgtfacq_e)[0];
	cl_mem vt_d          = resourcesCL.getBuffers("vt",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vt)[0];
	cl_mem ddxn_z_half_d = resourcesCL.getBuffers("ddxn_z_half", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ddxn_z_half)[0];
	cl_mem ddxt_z_half_d = resourcesCL.getBuffers("ddxt_z_half", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ddxt_z_half)[0];
	cl_mem vn_half_d     = resourcesCL.getBuffers("vn_half",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vn_half)[0];
	cl_mem z_kin_hor_e_d = resourcesCL.getBuffers("z_kin_hor_e", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_kin_hor_e)[0];
	cl_mem z_vt_ie_d     = resourcesCL.getBuffers("z_vt_ie",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_vt_ie)[0];
	cl_mem z_concorr_e_d = resourcesCL.getBuffers("z_concorr_e", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_concorr_e)[0];
	cl_mem localStart_d  = resourcesCL.getBuffers("localStart1", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d    = resourcesCL.getBuffers("localEnd1",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfac_e_d,    CL_FALSE, 0, sizeof(double)*size_wgtfac_e,    wgtfac_e,    0, NULL, &ev_write_wgtfac_e) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfacq1_e_d,  CL_FALSE, 0, sizeof(double)*size_wgtfacq1_e,  wgtfacq1_e,  0, NULL, &ev_write_wgtfacq1_e) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), wgtfacq_e_d,   CL_FALSE, 0, sizeof(double)*size_wgtfacq_e,   wgtfacq_e,   0, NULL, &ev_write_wgtfacq_e) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddxn_z_half_d, CL_FALSE, 0, sizeof(double)*size_ddxn_z_half, ddxn_z_half, 0, NULL, &ev_write_ddxn_z_half) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddxt_z_half_d, CL_FALSE, 0, sizeof(double)*size_ddxt_z_half, ddxt_z_half, 0, NULL, &ev_write_ddxt_z_half) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,  CL_FALSE, 0, sizeof(int)*size_localStart,     localStart,  0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,    CL_FALSE, 0, sizeof(int)*size_localEnd,       localEnd,    0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_d,          CL_FALSE, 0, sizeof(double)*size_vn,          vn,          0, NULL, &ev_write_vn) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vt_d,          CL_FALSE, 0, sizeof(double)*size_vt,          vt,          0, NULL, &ev_write_vt) );
		
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_half_d, CL_FALSE, 0, sizeof(double)*size_vn_half, vn_half, 0, NULL, &ev_write_vn_half) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_vt_ie_d, CL_FALSE, 0, sizeof(double)*size_z_vt_ie, z_vt_ie, 0, NULL, &ev_write_z_vt_ie) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_kin_hor_e_d, CL_FALSE, 0, sizeof(double)*size_z_kin_hor_e, z_kin_hor_e, 0, NULL, &ev_write_z_kin_hor_e) );
	}
	
#ifdef TEST
	cl_event ev_write_z_concorr_e;
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_concorr_e_d, CL_FALSE, 0, sizeof(double)*size_z_concorr_e, z_concorr_e, 0, NULL, &ev_write_z_concorr_e) );
	
	{
		cout << "checking transfers for ocl_part1...\n";
		double * wgtfac_e_tmp    = new double[size_wgtfac_e];
		double * wgtfacq_e_tmp   = new double[size_wgtfacq_e];
		double * wgtfacq1_e_tmp  = new double[size_wgtfacq1_e];
		double * ddxn_z_half_tmp = new double[size_ddxn_z_half];
		double * ddxt_z_half_tmp = new double[size_ddxt_z_half];
		double * vn_tmp          = new double[size_vn];
		double * vt_tmp          = new double[size_vt];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfac_e_d,    CL_FALSE, 0, sizeof(cl_double)*size_wgtfac_e,    wgtfac_e_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfacq_e_d,   CL_FALSE, 0, sizeof(cl_double)*size_wgtfacq_e,   wgtfacq_e_tmp,   0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), wgtfacq1_e_d,  CL_FALSE, 0, sizeof(cl_double)*size_wgtfacq1_e,  wgtfacq1_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddxn_z_half_d, CL_FALSE, 0, sizeof(cl_double)*size_ddxn_z_half, ddxn_z_half_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddxt_z_half_d, CL_FALSE, 0, sizeof(cl_double)*size_ddxt_z_half, ddxt_z_half_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_d,          CL_FALSE, 0, sizeof(cl_double)*size_vn,          vn_tmp,          0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vt_d,          CL_FALSE, 0, sizeof(cl_double)*size_vt,          vt_tmp,          0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		for (int i=0; i<size_wgtfac_e; i++)
			assert(wgtfac_e[i] == wgtfac_e_tmp[i]);
		
		for (int i=0; i<size_wgtfacq_e; i++)
			assert(wgtfacq_e[i] == wgtfacq_e_tmp[i]);
		
		for (int i=0; i<size_wgtfacq1_e; i++)
			assert(wgtfacq1_e[i] == wgtfacq1_e_tmp[i]);
		
		for (int i=0; i<size_ddxn_z_half; i++)
			assert(ddxn_z_half[i] == ddxn_z_half_tmp[i]);
		
		for (int i=0; i<size_ddxt_z_half; i++)
			assert(ddxt_z_half[i] == ddxt_z_half_tmp[i]);
		
		int zero[1] = {0};
		checkresultsCPP_(vn,vn_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"part1:vn");
		
		checkresultsCPP_(vt,vt_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"part1:vt");
		
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
		delete [] wgtfac_e_tmp;
		delete [] wgtfacq_e_tmp;
		delete [] wgtfacq1_e_tmp;
		delete [] ddxn_z_half_tmp;
		delete [] ddxt_z_half_tmp;
		delete [] vn_tmp;
		delete [] vt_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
	}
#endif
	
	if (*istep == 1)
	{
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.push_start("part1:step1a");
		}
		
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 2, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 3, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 4, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 5, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 6, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 7, sizeof(cl_mem), (void*)&wgtfac_e_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 8, sizeof(cl_mem), (void*)&vn_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a, 9, sizeof(cl_mem), (void*)&vt_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a,10, sizeof(cl_mem), (void*)&vn_half_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1a,11, sizeof(cl_mem), (void*)&z_kin_hor_e_d) );
		
		//cl_event ev_wait_list_1a[7] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfac_e, ev_write_vn, ev_write_vt, ev_write_vn_half, ev_write_z_kin_hor_e };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1a, 3, NULL, global_1a, NULL, 7, ev_wait_list_1a, &ev_kernel_step1a) );
		if (ocl.bFirstIter)
		{
			cl_event ev_wait_list_1a[4] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfac_e, ev_write_vn_half };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1a, 3, NULL, global_1a, local_1a, 4, ev_wait_list_1a, &ev_kernel_step1a) );
#ifndef TURBO_MODE
		}
		else if (*istep==1)
		{
			cl_event ev_wait_list_1a[1] = { ev_write_vn_half };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1a, 3, NULL, global_1a, local_1a, 1, ev_wait_list_1a, &ev_kernel_step1a) );
#endif
		}
		else
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1a, 3, NULL, global_1a, local_1a, 0, NULL, &ev_kernel_step1a) );
		
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.pop_stop();
			profilerKernels.push_start("part1:step1b");
		}
		
		int vert_nested = (int)*l_vert_nested;
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 2, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 3, sizeof(int32_t), (void*)nlevp1) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 4, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 5, sizeof(int32_t), (void*)&vert_nested) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 6, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 7, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 8, sizeof(cl_mem), (void*)&wgtfacq1_e_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b, 9, sizeof(cl_mem), (void*)&wgtfacq_e_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b,10, sizeof(cl_mem), (void*)&vn_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b,11, sizeof(cl_mem), (void*)&vt_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b,12, sizeof(cl_mem), (void*)&vn_half_d) );
		CheckCL( clSetKernelArg(ocl.kernel_step1b,13, sizeof(cl_mem), (void*)&z_kin_hor_e_d) );
		
		//cl_event ev_wait_list_1b[8] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq1_e, ev_write_wgtfacq_e, ev_write_vn, ev_write_vt, ev_write_vn_half, ev_write_z_kin_hor_e };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1b, 2, NULL, global_1b, NULL, 8, ev_wait_list_1b, &ev_kernel_step1b) );
		if (ocl.bFirstIter)
		{
			cl_event ev_wait_list_1b[5] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq1_e, ev_write_wgtfacq_e, ev_write_vn_half };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1b, 2, NULL, global_1b, NULL, 5, ev_wait_list_1b, &ev_kernel_step1b) );
#ifndef TURBO_MODE
		}
		else if (*istep==1)
		{
			cl_event ev_wait_list_1b[1] = { ev_write_vn_half };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1b, 2, NULL, global_1b, NULL, 1, ev_wait_list_1b, &ev_kernel_step1b) );
#endif
		}
		else
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1b, 2, NULL, global_1b, NULL, 0, NULL, &ev_kernel_step1b) );
		
		 
		if (ocl.bFirstIter)
		{
			if (bProfiling)
			{
				clFlush(oclhelper.getCommandQueue());
				clFinish(oclhelper.getCommandQueue());
				profilerKernels.pop_stop();
				profilerKernels.push_start("part1:step1ia");
			}
			
			const int i_endblkm1 = *i_endblk-1;
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 0, sizeof(int32_t), (void*)&i_startblkm1) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 1, sizeof(int32_t), (void*)nblks_e) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 2, sizeof(int32_t), (void*)nflat) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 3, sizeof(int32_t), (void*)nlev) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 4, sizeof(int32_t), (void*)nlevp1) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 5, sizeof(int32_t), (void*)nproma) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 6, sizeof(cl_mem), (void*)&localStart_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 7, sizeof(cl_mem), (void*)&localEnd_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 8, sizeof(cl_mem), (void*)&wgtfac_e_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia, 9, sizeof(cl_mem), (void*)&vt_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ia,10, sizeof(cl_mem), (void*)&z_vt_ie_d) );
			
			//cl_event ev_wait_list_1ia[5] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfac_e, ev_write_vt, ev_write_z_vt_ie };
			//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1ia, 3, NULL, global_1ia, NULL, 5, ev_wait_list_1ia, &ev_kernel_step1ia) );
			cl_event ev_wait_list_1ia[3] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfac_e };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1ia, 3, NULL, global_1ia, NULL, 3, ev_wait_list_1ia, &ev_kernel_step1ia) );
			
			if (bProfiling)
			{
				clFlush(oclhelper.getCommandQueue());
				clFinish(oclhelper.getCommandQueue());
				profilerKernels.pop_stop();
				profilerKernels.push_start("part1:step1ib");
			}
			
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 0, sizeof(int32_t), (void*)&i_startblkm1) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 1, sizeof(int32_t), (void*)nblks_e) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 2, sizeof(int32_t), (void*)nflat) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 3, sizeof(int32_t), (void*)nlev) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 4, sizeof(int32_t), (void*)nlevp1) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 5, sizeof(int32_t), (void*)nproma) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 6, sizeof(cl_mem), (void*)&localStart_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 7, sizeof(cl_mem), (void*)&localEnd_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 8, sizeof(cl_mem), (void*)&wgtfacq_e_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib, 9, sizeof(cl_mem), (void*)&vt_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ib,10, sizeof(cl_mem), (void*)&z_vt_ie_d) );
			
			//cl_event ev_wait_list_1ib[5] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq_e, ev_write_vt, ev_write_z_vt_ie };
			//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1ib, 2, NULL, global_1ib, NULL, 5, ev_wait_list_1ib, &ev_kernel_step1ib) );
			cl_event ev_wait_list_1ib[3] = { ev_write_localStart, ev_write_localEnd, ev_write_wgtfacq_e };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1ib, 2, NULL, global_1ib, NULL, 3, ev_wait_list_1ib, &ev_kernel_step1ib) );
			
			if (bProfiling)
			{
				clFlush(oclhelper.getCommandQueue());
				clFinish(oclhelper.getCommandQueue());
				profilerKernels.pop_stop();
				profilerKernels.push_start("part1:step1ic");
			} 
			
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 0, sizeof(int32_t), (void*)&i_startblkm1) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 1, sizeof(int32_t), (void*)i_endblk) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 2, sizeof(int32_t), (void*)nflat) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 3, sizeof(int32_t), (void*)nlev) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 4, sizeof(int32_t), (void*)nlevp1) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 5, sizeof(int32_t), (void*)nproma) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 6, sizeof(cl_mem), (void*)&localStart_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 7, sizeof(cl_mem), (void*)&localEnd_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 8, sizeof(cl_mem), (void*)&vn_half_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic, 9, sizeof(cl_mem), (void*)&z_vt_ie_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic,10, sizeof(cl_mem), (void*)&ddxn_z_half_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic,11, sizeof(cl_mem), (void*)&ddxt_z_half_d) );
			CheckCL( clSetKernelArg(ocl.kernel_step1ic,12, sizeof(cl_mem), (void*)&z_concorr_e_d) );
			
			cl_event ev_wait_list_1ic[9] = { ev_write_localStart, ev_write_localEnd, ev_write_vn_half, ev_write_ddxn_z_half, ev_write_ddxt_z_half, ev_kernel_step1a, ev_kernel_step1b, ev_kernel_step1ia, ev_kernel_step1ib };
			CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_step1ic, 3, NULL, global_1ic, NULL, 9, ev_wait_list_1ic, &ev_kernel_step1ic) );
		}
	}
	else
	{
		if (bProfiling)
		{
			clFlush(oclhelper.getCommandQueue());
			clFinish(oclhelper.getCommandQueue());
			profilerKernels.push_start("part1:stepn");
		}
		
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 0, sizeof(int32_t), (void*)&i_startblkm1) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 1, sizeof(int32_t), (void*)i_endblk) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 2, sizeof(int32_t), (void*)nlev) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 3, sizeof(int32_t), (void*)nproma) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 4, sizeof(cl_mem), (void*)&localStart_d) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 5, sizeof(cl_mem), (void*)&localEnd_d) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 6, sizeof(cl_mem), (void*)&vn_d) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 7, sizeof(cl_mem), (void*)&vt_d) );
		CheckCL( clSetKernelArg(ocl.kernel_stepn, 8, sizeof(cl_mem), (void*)&z_kin_hor_e_d) );
		
		//cl_event ev_wait_list_n[5] = { ev_write_localStart, ev_write_localEnd, ev_write_vn, ev_write_vt, ev_write_z_kin_hor_e };
		//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_stepn, 3, NULL, global_n, NULL, 5, ev_wait_list_n, &ev_kernel_stepn) );
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_stepn, 1, NULL, global_n, NULL, 0, NULL, &ev_kernel_stepn) );
	}
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_half_d,     CL_FALSE, 0, sizeof(cl_double)*size_vn_half,     vn_half,     0, NULL, &ev_read_vn_half) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_kin_hor_e_d, CL_FALSE, 0, sizeof(cl_double)*size_z_kin_hor_e, z_kin_hor_e, 0, NULL, &ev_read_z_kin_hor_e) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_concorr_e_d, CL_FALSE, 0, sizeof(cl_double)*size_z_concorr_e, z_concorr_e, 0, NULL, &ev_read_z_concorr_e) );
#endif
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	int zero[1] = {0};
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
};

extern "C" void ocl_part2_(int * i_startblk,
						   int * i_endblk,
						   int * nblks_e,
						   int * nblks_c,
						   int * slev,
						   int * elev,
						   int * nlev,
						   int * nlevp1,
						   int * nproma,
						   int * localStart,
						   int * localEnd,
						   int * istep,
						   int * icidx,
						   int * icblk,
						   double * vn_half,
						   double * c_lin_e,
						   double * w, int * nameIdx,
						   double * inv_dual_edge_length,
						   double * e_kinh,
						   double * z_vnw,
						   double * z_ddxn_ekin_e)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_part2");
#endif
	
	const int i_startblkm1 = *i_startblk - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	const int slevm1 = *slev-1;
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event ev_write_icidx, ev_write_icblk/*, ev_write_vn_half*/, ev_write_c_lin_e, ev_write_w, ev_write_inv_dual_edge_length, ev_write_e_kinh;
	cl_event ev_kernel_phase2;
	
	
	const int size_icidx                = (*nproma)*(*nblks_e)*2;
	const int size_icblk                = (*nproma)*(*nblks_e)*2;
	const int size_vn_half              = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_c_lin_e              = (*nproma)*2*(*nblks_e);
	const int size_w                    = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_inv_dual_edge_length = (*nproma)*(*nblks_e);
	const int size_e_kinh               = (*nproma)*(*nlev)*(*nblks_c);
	const int size_z_vnw                = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_z_ddxn_ekin_e        = (*nproma)*(*nlev)*(*nblks_e);
	const int size_localStart           = nBlocks;
	const int size_localEnd             = nBlocks;
	
	size_t local_phase2[3] = { 1, *nlev, 1 };
	size_t global_phase2[3] = { nBlocks, *nlev, *nproma };
	
	cl_mem w_d                    = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_w)[0];
	cl_mem icidx_d                = resourcesCL.getBuffers("icidx",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_icidx)[0];
	cl_mem icblk_d                = resourcesCL.getBuffers("icblk",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_icblk)[0];
	cl_mem vn_half_d              = resourcesCL.getBuffers("vn_half",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vn_half)[0];
	cl_mem c_lin_e_d              = resourcesCL.getBuffers("c_lin_e",              1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_c_lin_e)[0];
	cl_mem inv_dual_edge_length_d = resourcesCL.getBuffers("inv_dual_edge_length", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_inv_dual_edge_length)[0];
	cl_mem e_kinh_d               = resourcesCL.getBuffers("e_kinh",               1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_e_kinh)[0];
	cl_mem localStart_d           = resourcesCL.getBuffers("localStart2",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d             = resourcesCL.getBuffers("localEnd2",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem z_vnw_d                = resourcesCL.getBuffers("z_vnw",                1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_vnw)[0];
	cl_mem z_ddxn_ekin_e_d        = resourcesCL.getBuffers("z_ddxn_ekin_e",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_ddxn_ekin_e)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), icidx_d,                CL_FALSE, 0, sizeof(int)*size_icidx,                   icidx,                0, NULL, &ev_write_icidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), icblk_d,                CL_FALSE, 0, sizeof(int)*size_icblk,                   icblk,                0, NULL, &ev_write_icblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,           CL_FALSE, 0, sizeof(int)*size_localStart,              localStart,           0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,             CL_FALSE, 0, sizeof(int)*size_localEnd,                localEnd,             0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_half_d,              CL_FALSE, 0, sizeof(double)*size_vn_half,              vn_half,              0, NULL, &ev_write_vn_half) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), c_lin_e_d,              CL_FALSE, 0, sizeof(double)*size_c_lin_e,              c_lin_e,              0, NULL, &ev_write_c_lin_e) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_d,                    CL_FALSE, 0, sizeof(double)*size_w,                    w,                    0, NULL, &ev_write_w) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), inv_dual_edge_length_d, CL_FALSE, 0, sizeof(double)*size_inv_dual_edge_length, inv_dual_edge_length, 0, NULL, &ev_write_inv_dual_edge_length) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), e_kinh_d,               CL_FALSE, 0, sizeof(double)*size_e_kinh,               e_kinh,               0, NULL, &ev_write_e_kinh) );
	}
	
#ifdef TEST
	{
		cl_event ev_write_z_vnw, ev_write_z_ddxn_ekin_e, ev_read_vn_half;
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_vnw_d,                CL_FALSE, 0, sizeof(double)*size_z_vnw,                z_vnw,                0, NULL, &ev_write_z_vnw) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_ddxn_ekin_e_d,        CL_FALSE, 0, sizeof(double)*size_z_ddxn_ekin_e,        z_ddxn_ekin_e,        0, NULL, &ev_write_z_ddxn_ekin_e) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
	}
	
	{
		cout << "checking transfers for ocl_part2...\n";
		int * icidx_tmp   = new int[size_icidx];
		int * icblk_tmp   = new int[size_icblk];
		double * vn_half_tmp = new double[size_vn_half];
		double * e_kinh_tmp  = new double[size_e_kinh];
		int zero[1] = {0};
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), icidx_d,   CL_FALSE, 0, sizeof(int)*size_icidx,   icidx_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), icblk_d,   CL_FALSE, 0, sizeof(int)*size_icblk,   icblk_tmp,   0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_half_d, CL_FALSE, 0, sizeof(cl_double)*size_vn_half, vn_half_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), e_kinh_d,  CL_FALSE, 0, sizeof(cl_double)*size_e_kinh,  e_kinh_tmp,  0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		for (int i=0; i<size_icidx; i++)
			assert(icidx[i] == icidx_tmp[i]);
		
		for (int i=0; i<size_icblk; i++)
			assert(icblk[i] == icblk_tmp[i]);
		
		checkresultsCPP_(vn_half,vn_half_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_e,nblks_e,"part2:vn_half");
		
		checkresultsCPP_(e_kinh,e_kinh_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"part2:e_kinh");
		
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
		delete [] icidx_tmp;
		delete [] icblk_tmp;
		delete [] vn_half_tmp;
		delete [] e_kinh_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("part2");
	}
#endif
	
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 2, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 3, sizeof(int32_t), (void*)&slevm1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 4, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 5, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 6, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 7, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 8, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2, 9, sizeof(cl_mem), (void*)&icidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,10, sizeof(cl_mem), (void*)&icblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,11, sizeof(cl_mem), (void*)&vn_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,12, sizeof(cl_mem), (void*)&c_lin_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,13, sizeof(cl_mem), (void*)&w_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,14, sizeof(cl_mem), (void*)&inv_dual_edge_length_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,15, sizeof(cl_mem), (void*)&e_kinh_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,16, sizeof(cl_mem), (void*)&z_vnw_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase2,17, sizeof(cl_mem), (void*)&z_ddxn_ekin_e_d) );
	
	//cl_event ev_wait_list[11] = { ev_write_icidx, ev_write_icblk, ev_write_vn_half, ev_write_c_lin_e, ev_write_w, ev_write_inv_dual_edge_length, ev_write_e_kinh, ev_write_localStart, ev_write_localEnd, ev_write_z_vnw, ev_write_z_ddxn_ekin_e };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase2, 3, NULL, global_phase2, NULL, 11, ev_wait_list, &ev_kernel_phase2) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list[7] = { ev_write_icidx, ev_write_icblk, ev_write_c_lin_e, ev_write_w, ev_write_inv_dual_edge_length, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase2, 3, NULL, global_phase2, local_phase2, 7, ev_wait_list, &ev_kernel_phase2) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list[3] = { ev_write_c_lin_e, ev_write_w, ev_write_inv_dual_edge_length };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase2, 3, NULL, global_phase2, local_phase2, 3, ev_wait_list, &ev_kernel_phase2) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase2, 3, NULL, global_phase2, local_phase2, 0, NULL, &ev_kernel_phase2) );
	
	
#ifdef TEST
	cl_event ev_read_z_vnw, ev_read_z_ddxn_ekin_e;
	{
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_vnw_d,         CL_FALSE, 0, sizeof(cl_double)*size_z_vnw,         z_vnw,         1, &ev_kernel_phase2, &ev_read_z_vnw) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_ddxn_ekin_e_d, CL_FALSE, 0, sizeof(cl_double)*size_z_ddxn_ekin_e, z_ddxn_ekin_e, 1, &ev_kernel_phase2, &ev_read_z_ddxn_ekin_e) );
	}
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profilerKernels.pop_stop();
		profiler.pop_stop();
		//profiler.pop_stop();
	}
};

extern "C" void ocl_part3_(int * i_startblk,
						   int * i_endblk,
						   int * nblks_e,
						   int * nblks_c,
						   int * nflat,
						   int * nlev,
						   int * nlevp1,
						   int * nproma,
						   int * localStart,
						   int * localEnd,
						   int * istep,
						   int * ieidx,
						   int * ieblk,
						   double * vn_half,
						   double * geofac_div,
						   double * w, int * nameIdx,
						   double * z_vnw,
						   double * z_hadv_w,
						   double * w_con,
						   double * w_concorr_c,
						   double * z_w_con_c_full)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_part3");
#endif
	
	const int i_startblkm1 = *i_startblk - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	
	cl_event ev_write_localStart, ev_write_localEnd, ev_write_ieidx, ev_write_ieblk/*, ev_write_vn_half*/, ev_write_geofac_div/*, ev_write_w, ev_write_z_vnw*/, ev_write_w_con, ev_write_w_concorr_c;
	cl_event ev_kernel_phase3a, ev_kernel_phase3b, ev_kernel_phase3c;
	cl_event ev_read_z_hadv_w, ev_read_w_con;
	
	size_t local_phase3a[3] = { 1, *nlev, 1 };
	size_t local_phase3c[3] = { 1, *nlevp1, 1 };
	
	size_t global_phase3a[3] = { nBlocks, *nlev, *nproma };
	//size_t global_phase3b[3] = { nBlocks, *nlev-*nflat, *nproma };
	size_t global_phase3c[3] = { nBlocks, *nlevp1, *nproma };
	
	const int size_localStart     = nBlocks;
	const int size_localEnd       = nBlocks;
	const int size_ieidx          = (*nproma)*(*nblks_c)*3;
	const int size_ieblk          = (*nproma)*(*nblks_c)*3;
	const int size_vn_half        = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_geofac_div     = (*nproma)*3*(*nblks_e);
	const int size_w              = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_z_vnw          = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_z_hadv_w       = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_w_con          = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_w_concorr_c    = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_z_w_con_c_full = (*nproma)*(*nlev)*(*nblks_c);
	
	cl_mem w_d              = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_w)[0];
	cl_mem localStart_d     = resourcesCL.getBuffers("localStart3",    1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d       = resourcesCL.getBuffers("localEnd3",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem ieidx_d          = resourcesCL.getBuffers("ieidx",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_ieidx)[0];
	cl_mem ieblk_d          = resourcesCL.getBuffers("ieblk",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ieblk)[0];
	cl_mem vn_half_d        = resourcesCL.getBuffers("vn_half",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vn_half)[0];
	cl_mem geofac_div_d     = resourcesCL.getBuffers("geofac_div",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_geofac_div)[0];
	cl_mem z_vnw_d          = resourcesCL.getBuffers("z_vnw",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_vnw)[0];
	cl_mem w_concorr_c_d    = resourcesCL.getBuffers("w_concorr_c",    1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_w_concorr_c)[0];
	cl_mem z_hadv_w_d       = resourcesCL.getBuffers("z_hadv_w",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_hadv_w)[0];
	cl_mem w_con_d          = resourcesCL.getBuffers("w_con",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_w_con)[0];
	cl_mem z_w_con_c_full_d = resourcesCL.getBuffers("z_w_con_c_full", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_w_con_c_full)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ieidx_d,          CL_FALSE, 0, sizeof(cl_int)*size_ieidx,             ieidx,          0, NULL, &ev_write_ieidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ieblk_d,          CL_FALSE, 0, sizeof(cl_int)*size_ieblk,             ieblk,          0, NULL, &ev_write_ieblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,     CL_FALSE, 0, sizeof(cl_int)*size_localStart,        localStart,     0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,       CL_FALSE, 0, sizeof(cl_int)*size_localEnd,          localEnd,       0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_half_d,        CL_FALSE, 0, sizeof(cl_double)*size_vn_half,        vn_half,        0, NULL, &ev_write_vn_half) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), geofac_div_d,     CL_FALSE, 0, sizeof(cl_double)*size_geofac_div,     geofac_div,     0, NULL, &ev_write_geofac_div) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_d,              CL_FALSE, 0, sizeof(cl_double)*size_w,              w,              0, NULL, &ev_write_w) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_vnw_d,          CL_FALSE, 0, sizeof(cl_double)*size_z_vnw,          z_vnw,          0, NULL, &ev_write_z_vnw) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_concorr_c_d,    CL_FALSE, 0, sizeof(cl_double)*size_w_concorr_c,    w_concorr_c,    0, NULL, &ev_write_w_concorr_c) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_con_d,          CL_FALSE, 0, sizeof(cl_double)*size_w_con,          w_con,          0, NULL, &ev_write_w_con) );
	}
	
#ifdef TEST
	cl_event ev_write_z_hadv_w, ev_write_z_w_con_c_full;
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_hadv_w_d,       CL_FALSE, 0, sizeof(cl_double)*size_z_hadv_w,       z_hadv_w,       0, NULL, &ev_write_z_hadv_w) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_w_con_c_full_d, CL_FALSE, 0, sizeof(cl_double)*size_z_w_con_c_full, z_w_con_c_full, 0, NULL, &ev_write_z_w_con_c_full) );
	
	{
		cout << "checking transfers for ocl_part3...\n";
		int * ieidx_tmp   = new int[size_ieidx];
		int * ieblk_tmp   = new int[size_ieblk];
		double * vn_half_tmp = new double[size_vn_half];
		double * w_tmp       = new double[size_w];
		double * z_vnw_tmp   = new double[size_z_vnw];
		double * w_concorr_c_tmp = new double[size_w_concorr_c];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ieidx_d,   CL_FALSE, 0, sizeof(int)*size_ieidx,   ieidx_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ieblk_d,   CL_FALSE, 0, sizeof(int)*size_ieblk,   ieblk_tmp,   0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_half_d, CL_FALSE, 0, sizeof(cl_double)*size_vn_half, vn_half_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_d,       CL_FALSE, 0, sizeof(cl_double)*size_w,       w_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_vnw_d,   CL_FALSE, 0, sizeof(cl_double)*size_z_vnw,   z_vnw_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_concorr_c_d,   CL_FALSE, 0, sizeof(cl_double)*size_w_concorr_c, w_concorr_c_tmp,  0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		for (int i=0; i<size_ieidx; i++)
			assert(ieidx[i] == ieidx_tmp[i]);
		
		for (int i=0; i<size_ieblk; i++)
			assert(ieblk[i] == ieblk_tmp[i]);
		
		for (int i=0; i<size_vn_half; i++)
			assert(fabs(vn_half[i] - vn_half_tmp[i])<1.e-6 || (isnan(fabs(vn_half[i])) && isnan(fabs(vn_half_tmp[i]))));

		for (int i=0; i<size_z_vnw; i++)
			assert(fabs(z_vnw[i] - z_vnw_tmp[i])<1.e-6 || (isnan(fabs(z_vnw[i])) && isnan(fabs(z_vnw_tmp[i]))));
		
		int zero[1] = {0};
		checkresultsCPP_(w,w_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"part3:w");
		
		checkresultsCPP_(w_concorr_c,w_concorr_c_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"part3:w_concorr_c");
		
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
		delete [] ieidx_tmp;
		delete [] ieblk_tmp;
		delete [] vn_half_tmp;
		delete [] w_tmp;
		delete [] z_vnw_tmp;
		delete [] w_concorr_c_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("part3:phase3a");
	}
#endif
	int argIdx=0;
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(int32_t), (void*)nflat) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&ieidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&ieblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&vn_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&geofac_div_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&w_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&z_vnw_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&z_hadv_w_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&w_concorr_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3a, argIdx++, sizeof(cl_mem), (void*)&w_con_d) );
	
	//cl_event ev_wait_list_phase3a[10] = { ev_write_localStart, ev_write_localEnd, ev_write_ieidx, ev_write_ieblk, ev_write_vn_half, ev_write_geofac_div, ev_write_w, ev_write_z_vnw, ev_write_z_hadv_w, ev_write_w_con };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3a, 3, NULL, global_phase3a, NULL, 10, ev_wait_list_phase3a, &ev_kernel_phase3a) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_phase3a[6] = { ev_write_localStart, ev_write_localEnd, ev_write_ieidx, ev_write_ieblk, ev_write_geofac_div, ev_write_w_con };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3a, 3, NULL, global_phase3a, local_phase3a, 6, ev_wait_list_phase3a, &ev_kernel_phase3a) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_phase3a[2] = { ev_write_geofac_div, ev_write_w_con };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3a, 3, NULL, global_phase3a, local_phase3a, 2, ev_wait_list_phase3a, &ev_kernel_phase3a) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3a, 3, NULL, global_phase3a, local_phase3a, 0, NULL, &ev_kernel_phase3a) );
	/*
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("part3:phase3b");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 2, sizeof(int32_t), (void*)nflat) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 4, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 5, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 6, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 7, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 8, sizeof(cl_mem), (void*)&w_concorr_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3b, 9, sizeof(cl_mem), (void*)&w_con_d) );
	
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_phase3b[4] = { ev_write_localStart, ev_write_localEnd, ev_write_w_con, ev_kernel_phase3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3b, 3, NULL, global_phase3b, NULL, 4, ev_wait_list_phase3b, &ev_kernel_phase3b) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_phase3b[2] = { ev_write_w_con, ev_kernel_phase3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3b, 3, NULL, global_phase3b, NULL, 2, ev_wait_list_phase3b, &ev_kernel_phase3b) );
#endif
	}
	else
	{
		cl_event ev_wait_list_phase3b[1] = { ev_kernel_phase3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3b, 3, NULL, global_phase3b, NULL, 1, ev_wait_list_phase3b, &ev_kernel_phase3b) );
	}
	*/
	
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profilerKernels.push_start("part3:phase3c");
	}
	
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 5, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 6, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 7, sizeof(cl_mem), (void*)&w_con_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase3c, 8, sizeof(cl_mem), (void*)&z_w_con_c_full_d) );
	
	//cl_event ev_wait_list_phase3c[5] = { ev_write_localStart, ev_write_localEnd, ev_write_w_con, ev_write_z_w_con_c_full, ev_kernel_phase3b };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3c, 3, NULL, global_phase3c, NULL, 5, ev_wait_list_phase3c, &ev_kernel_phase3c) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list_phase3c[4] = { ev_write_localStart, ev_write_localEnd, ev_write_w_con, ev_kernel_phase3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3c, 3, NULL, global_phase3c, local_phase3c, 4, ev_wait_list_phase3c, &ev_kernel_phase3c) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list_phase3c[2] = { ev_write_w_con, ev_kernel_phase3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3c, 3, NULL, global_phase3c, local_phase3c, 2, ev_wait_list_phase3c, &ev_kernel_phase3c) );
#endif
	}
	else
	{
		cl_event ev_wait_list_phase3c[1] = { ev_kernel_phase3a };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase3c, 3, NULL, global_phase3c, local_phase3c, 1, ev_wait_list_phase3c, &ev_kernel_phase3c) );
	}
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif

#ifndef TURBO_MODE
	if (*istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_con_d,          CL_FALSE, 0, sizeof(cl_double)*size_w_con,          w_con,          1, &ev_kernel_phase3a, &ev_read_w_con) );
#endif
	
#ifdef TEST
	cl_event ev_read_z_w_con_full;
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_con_d,          CL_FALSE, 0, sizeof(cl_double)*size_w_con,          w_con,          1, &ev_kernel_phase3a, &ev_read_w_con) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_hadv_w_d,       CL_FALSE, 0, sizeof(cl_double)*size_z_hadv_w,       z_hadv_w,       1, &ev_kernel_phase3a, &ev_read_z_hadv_w) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_w_con_c_full_d, CL_FALSE, 0, sizeof(cl_double)*size_z_w_con_c_full, z_w_con_c_full, 1, &ev_kernel_phase3c, &ev_read_z_w_con_full) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
};

extern "C" void ocl_part4_(int * i_startblk,
						   int * i_endblk,
						   int * nblks_e,
						   int * nblks_c,
						   int * nblks_v,
						   int * nlev,
						   int * nlevp1,
						   int * nproma,
						   int * ntnd,
						   int * localStart,
						   int * localEnd,
						   int * istep,
						   int * ividx,
						   int * ivblk,
						   int * icidx,
						   int * icblk,
						   double * z_ddxn_ekin_e,
						   double * vt,
						   double * f_e,
						   double * omega_z,
						   double * c_lin_e,
						   double * z_w_con_c_full,
						   double * vn_half,
						   double * ddqz_z_full_e,
						   double * ddt_vn_adv)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_part4");
#endif
	
	const int i_startblkm1 = *i_startblk - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event /*ev_write_icidx, ev_write_icblk, */ev_write_ividx, ev_write_ivblk;
	cl_event /*ev_write_z_ddxn_ekin_e, ev_write_vt, */ev_write_f_e, ev_write_omega_z/*, ev_write_c_lin_e, ev_write_z_w_con_c_full, ev_write_vn_half*/, ev_write_ddqz_z_full_e, ev_write_ddt_vn_adv;
	cl_event ev_kernel_phase4;
	cl_event ev_read_ddt_vn_adv;
	
	const int size_localStart     = nBlocks;
	const int size_localEnd       = nBlocks;
	const int size_icidx          = (*nproma)*(*nblks_e)*2;
	const int size_icblk          = (*nproma)*(*nblks_e)*2;
	const int size_ividx          = (*nproma)*(*nblks_e)*4;
	const int size_ivblk          = (*nproma)*(*nblks_e)*4;
	const int size_z_ddxn_ekin_e  = (*nproma)*(*nlev)*(*nblks_e);
	const int size_vt             = (*nproma)*(*nlev)*(*nblks_e);
	const int size_f_e            = (*nproma)*(*nblks_e);
	const int size_omega_z        = (*nproma)*(*nlev)*(*nblks_v); // is this correct?
	const int size_c_lin_e        = (*nproma)*2*(*nblks_e);
	const int size_z_w_con_c_full = (*nproma)*(*nlev)*(*nblks_c);
	const int size_vn_half        = (*nproma)*(*nlevp1)*(*nblks_e);
	const int size_ddqz_z_full_e  = (*nproma)*(*nlev)*(*nblks_e);
	const int size_ddt_vn_adv     = (*nproma)*(*nlev)*(*nblks_e)*3;
	
	size_t local_phase4[3] = { 1, *nlevp1, 1 };
	size_t global_phase4[3] = { nBlocks, *nlevp1, *nproma };
	
	cl_mem localStart_d     = resourcesCL.getBuffers("localStart4",    1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d       = resourcesCL.getBuffers("localEnd4",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem icidx_d          = resourcesCL.getBuffers("icidx",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_icidx)[0];
	cl_mem icblk_d          = resourcesCL.getBuffers("icblk",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_icblk)[0];
	cl_mem ividx_d          = resourcesCL.getBuffers("ividx",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_ividx)[0];
	cl_mem ivblk_d          = resourcesCL.getBuffers("ivblk",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_ivblk)[0];
	cl_mem z_ddxn_ekin_e_d  = resourcesCL.getBuffers("z_ddxn_ekin_e",  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_ddxn_ekin_e)[0];
	cl_mem vt_d             = resourcesCL.getBuffers("vt",             1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vt)[0];
	cl_mem f_e_d            = resourcesCL.getBuffers("f_e",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_f_e)[0];
	cl_mem omega_z_d        = resourcesCL.getBuffers("omega_z",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_omega_z)[0];
	cl_mem c_lin_e_d        = resourcesCL.getBuffers("c_lin_e",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_c_lin_e)[0];
	cl_mem z_w_con_c_full_d = resourcesCL.getBuffers("z_w_con_c_full", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_z_w_con_c_full)[0];
	cl_mem vn_half_d        = resourcesCL.getBuffers("vn_half",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vn_half)[0];
	cl_mem ddqz_z_full_e_d  = resourcesCL.getBuffers("ddqz_z_full_e",  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ddqz_z_full_e)[0];
	cl_mem ddt_vn_adv_d     = resourcesCL.getBuffers("ddt_vn_adv",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ddt_vn_adv)[0];

#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), icidx_d,          CL_FALSE, 0, sizeof(int)*size_icidx,             icidx,          0, NULL, &ev_write_icidx) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), icblk_d,          CL_FALSE, 0, sizeof(int)*size_icblk,             icblk,          0, NULL, &ev_write_icblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ividx_d,          CL_FALSE, 0, sizeof(int)*size_ividx,             ividx,          0, NULL, &ev_write_ividx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ivblk_d,          CL_FALSE, 0, sizeof(int)*size_ivblk,             ivblk,          0, NULL, &ev_write_ivblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddqz_z_full_e_d,  CL_FALSE, 0, sizeof(double)*size_ddqz_z_full_e,  ddqz_z_full_e,  0, NULL, &ev_write_ddqz_z_full_e) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,     CL_FALSE, 0, sizeof(int)*size_localStart,        localStart,     0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,       CL_FALSE, 0, sizeof(int)*size_localEnd,          localEnd,       0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_ddxn_ekin_e_d,  CL_FALSE, 0, sizeof(double)*size_z_ddxn_ekin_e,  z_ddxn_ekin_e,  0, NULL, &ev_write_z_ddxn_ekin_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vt_d,             CL_FALSE, 0, sizeof(double)*size_vt,             vt,             0, NULL, &ev_write_vt) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), f_e_d,            CL_FALSE, 0, sizeof(double)*size_f_e,            f_e,            0, NULL, &ev_write_f_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), omega_z_d,        CL_FALSE, 0, sizeof(double)*size_omega_z,        omega_z,        0, NULL, &ev_write_omega_z) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), c_lin_e_d,        CL_FALSE, 0, sizeof(double)*size_c_lin_e,        c_lin_e,        0, NULL, &ev_write_c_lin_e) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), z_w_con_c_full_d, CL_FALSE, 0, sizeof(double)*size_z_w_con_c_full, z_w_con_c_full, 0, NULL, &ev_write_z_w_con_c_full) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vn_half_d,        CL_FALSE, 0, sizeof(double)*size_vn_half,        vn_half,        0, NULL, &ev_write_vn_half) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_vn_adv_d,     CL_FALSE, 0, sizeof(double)*size_ddt_vn_adv,     ddt_vn_adv,     0, NULL, &ev_write_ddt_vn_adv) );
	}
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_part4...\n";
		int * icidx_tmp   = new int[size_icidx];
		int * icblk_tmp   = new int[size_icblk];
		int * ividx_tmp   = new int[size_ividx];
		int * ivblk_tmp   = new int[size_ivblk];
		double * ddqz_z_full_e_tmp = new double[size_ddqz_z_full_e];
		double * z_ddxn_ekin_e_tmp = new double[size_z_ddxn_ekin_e];
		double * vt_tmp      = new double[size_vt];
		double * c_lin_e_tmp = new double[size_c_lin_e];
		double * vn_half_tmp = new double[size_vn_half];
		double * omega_z_tmp = new double[size_omega_z];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), icidx_d,   CL_FALSE, 0, sizeof(int)*size_icidx,   icidx_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), icblk_d,   CL_FALSE, 0, sizeof(int)*size_icblk,   icblk_tmp,   0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ividx_d,   CL_FALSE, 0, sizeof(int)*size_ividx,   ividx_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ivblk_d,   CL_FALSE, 0, sizeof(int)*size_ivblk,   ivblk_tmp,   0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddqz_z_full_e_d, CL_FALSE, 0, sizeof(cl_double)*size_ddqz_z_full_e, ddqz_z_full_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), z_ddxn_ekin_e_d, CL_FALSE, 0, sizeof(cl_double)*size_z_ddxn_ekin_e, z_ddxn_ekin_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vt_d,      CL_FALSE, 0, sizeof(cl_double)*size_vt,      vt_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), c_lin_e_d, CL_FALSE, 0, sizeof(cl_double)*size_c_lin_e, c_lin_e_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vn_half_d, CL_FALSE, 0, sizeof(cl_double)*size_vn_half, vn_half_tmp,  0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), omega_z_d, CL_FALSE, 0, sizeof(cl_double)*size_omega_z, omega_z_tmp,  0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		int zero[1] = {0};
		
		for (int i=0; i<size_icidx; i++)
			assert(icidx[i] == icidx_tmp[i]);
		
		for (int i=0; i<size_icblk; i++)
			assert(icblk[i] == icblk_tmp[i]);
		
		for (int i=0; i<size_ividx; i++)
			assert(ividx[i] == ividx_tmp[i]);
		
		for (int i=0; i<size_ivblk; i++)
			assert(ivblk[i] == ivblk_tmp[i]);
		
		for (int i=0; i<size_ddqz_z_full_e; i++)
			assert(ddqz_z_full_e[i] == ddqz_z_full_e_tmp[i]);
		
		for (int i=0; i<size_vt; i++)
			assert(vt[i] == vt_tmp[i]);
		
		for (int i=0; i<size_c_lin_e; i++)
			assert(c_lin_e[i] == c_lin_e_tmp[i]);
		
		for (int i=0; i<size_vn_half; i++)
			assert(fabs(vn_half[i] - vn_half_tmp[i])<1.e-6 || (isnan(fabs(vn_half[i])) && isnan(fabs(vn_half_tmp[i]))));
		
		checkresultsCPP_(z_ddxn_ekin_e,z_ddxn_ekin_e_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"part4:z_ddxn_ekin_e");
		
		checkresultsCPP_(omega_z,omega_z_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_v,nblks_v,"part4:omega_z");
		
		
		delete [] icidx_tmp;
		delete [] icblk_tmp;
		delete [] ividx_tmp;
		delete [] ivblk_tmp;
		delete [] ddqz_z_full_e_tmp;
		delete [] z_ddxn_ekin_e_tmp;
		delete [] vt_tmp;
		delete [] c_lin_e_tmp;
		delete [] vn_half_tmp;
		delete [] omega_z_tmp;
		
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("part4");
	}
#endif
	
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 5, sizeof(int32_t), (void*)ntnd) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 6, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 7, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 8, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4, 9, sizeof(cl_mem), (void*)&ividx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,10, sizeof(cl_mem), (void*)&ivblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,11, sizeof(cl_mem), (void*)&icidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,12, sizeof(cl_mem), (void*)&icblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,13, sizeof(cl_mem), (void*)&z_ddxn_ekin_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,14, sizeof(cl_mem), (void*)&vt_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,15, sizeof(cl_mem), (void*)&f_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,16, sizeof(cl_mem), (void*)&omega_z_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,17, sizeof(cl_mem), (void*)&c_lin_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,18, sizeof(cl_mem), (void*)&z_w_con_c_full_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,19, sizeof(cl_mem), (void*)&vn_half_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,20, sizeof(cl_mem), (void*)&ddqz_z_full_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase4,21, sizeof(cl_mem), (void*)&ddt_vn_adv_d) );
	
	//cl_event ev_wait_list[15] = { ev_write_localStart, ev_write_localEnd, ev_write_ividx, ev_write_ivblk, ev_write_icidx, ev_write_icblk,
	//	                          ev_write_z_ddxn_ekin_e, ev_write_vt, ev_write_f_e, ev_write_omega_z, ev_write_c_lin_e, ev_write_z_w_con_c_full, ev_write_vn_half, ev_write_ddqz_z_full_e, ev_write_ddt_vn_adv };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase4, 3, NULL, global_phase4, NULL, 15, ev_wait_list, &ev_kernel_phase4) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list[7] = { ev_write_localStart, ev_write_localEnd, ev_write_ividx, ev_write_ivblk,
		                              ev_write_f_e, ev_write_ddqz_z_full_e, ev_write_ddt_vn_adv };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase4, 3, NULL, global_phase4, local_phase4, 7, ev_wait_list, &ev_kernel_phase4) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list[2] = { ev_write_f_e, ev_write_ddt_vn_adv };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase4, 3, NULL, global_phase4, local_phase4, 2, ev_wait_list, &ev_kernel_phase4) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase4, 3, NULL, global_phase4, local_phase4, 0, NULL, &ev_kernel_phase4) );
	
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_vn_adv_d, CL_FALSE, 0, sizeof(cl_double)*size_ddt_vn_adv, ddt_vn_adv, 1, &ev_kernel_phase4, &ev_read_ddt_vn_adv) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
};


extern "C" void ocl_part5_(int * i_startblk,
						   int * i_endblk,
						   int * nblks_c,
						   int * nlev,
						   int * nlevp1,
						   int * nproma,
						   int * ntnd,
						   int * localStart,
						   int * localEnd,
						   int * istep,
						   double * w_con,
						   double * w, int * nameIdx,
						   double * inv_ddqz_z_half2,
						   double * ddt_w_adv)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_part5");
#endif
	
	const int i_startblkm1 = *i_startblk - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	
	cl_event ev_write_localStart, ev_write_localEnd;
	cl_event /*ev_write_w_con, ev_write_w, */ev_write_inv_ddqz_z_half2, ev_write_ddt_w_adv1, ev_write_ddt_w_adv2, ev_write_ddt_w_adv3;
	cl_event ev_kernel_phase5;
	cl_event ev_read_ddt_w_adv1, ev_read_ddt_w_adv2, ev_read_ddt_w_adv3;
	
	const int size_localStart       = nBlocks;
	const int size_localEnd         = nBlocks;
	const int size_w_con            = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_w                = (*nproma)*(*nlevp1)*(*nblks_c);
	const int size_inv_ddqz_z_half2 = (*nproma)*(*nlev)*(*nblks_c);
	const int size_ddt_w_adv        = (*nproma)*(*nlevp1)*(*nblks_c);
	
	size_t local_phase5[3] = { 1, *nlevp1, 1 };
	size_t global_phase5[3] = { nBlocks, *nlevp1, *nproma };
	
	cl_mem w_d                = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_w)[0];
	cl_mem localStart_d       = resourcesCL.getBuffers("localStart5",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localStart)[0];
	cl_mem localEnd_d         = resourcesCL.getBuffers("localEnd5",        1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_localEnd)[0];
	cl_mem w_con_d            = resourcesCL.getBuffers("w_con",            1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_w_con)[0];
	cl_mem inv_ddqz_z_half2_d = resourcesCL.getBuffers("inv_ddqz_z_half2", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_inv_ddqz_z_half2)[0];
	cl_mem ddt_w_adv1_d       = resourcesCL.getBuffers("ddt_w_adv1",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ddt_w_adv)[0];
	cl_mem ddt_w_adv2_d       = resourcesCL.getBuffers("ddt_w_adv2",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ddt_w_adv)[0];
	cl_mem ddt_w_adv3_d       = resourcesCL.getBuffers("ddt_w_adv3",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_ddt_w_adv)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_half2_d, CL_FALSE, 0, sizeof(double)*size_inv_ddqz_z_half2, inv_ddqz_z_half2, 0, NULL, &ev_write_inv_ddqz_z_half2) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d,       CL_FALSE, 0, sizeof(int)*size_localStart,          localStart,       0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,         CL_FALSE, 0, sizeof(int)*size_localEnd,            localEnd,         0, NULL, &ev_write_localEnd) );
	}
	
#ifndef TURBO_MODE
	if (*istep==1)
#else
	if (ocl.bFirstIter)
#endif
	{
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_con_d,            CL_FALSE, 0, sizeof(double)*size_w_con,            w_con,            0, NULL, &ev_write_w_con) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), w_d,                CL_FALSE, 0, sizeof(double)*size_w,                w,                0, NULL, &ev_write_w) );
	}
	
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_w_adv1_d,        CL_FALSE, 0, sizeof(double)*size_ddt_w_adv,        ddt_w_adv,        0, NULL, &ev_write_ddt_w_adv1) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_w_adv2_d,        CL_FALSE, 0, sizeof(double)*size_ddt_w_adv,        &ddt_w_adv[size_ddt_w_adv],        0, NULL, &ev_write_ddt_w_adv2) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ddt_w_adv3_d,        CL_FALSE, 0, sizeof(double)*size_ddt_w_adv,        &ddt_w_adv[2*size_ddt_w_adv],        0, NULL, &ev_write_ddt_w_adv3) );
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_part5...\n";
		cl_event ev_read_w_con, ev_read_w;
		
		double * inv_ddqz_z_half2_tmp   = new double[size_inv_ddqz_z_half2];
		double * w_con_tmp = new double[size_w_con];
		double * w_tmp = new double[size_w];
		
		int * localStart_tmp = new int[nBlocks];
		int * localEnd_tmp = new int[nBlocks];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_TRUE, 0, sizeof(cl_int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_TRUE, 0, sizeof(cl_int)*nBlocks, localEnd_tmp,   0, NULL, NULL) );
		for (int i=0; i<nBlocks; i++)
		{
			assert(localStart[i] == localStart_tmp[i]);
			assert(localEnd[i] == localEnd_tmp[i]);
		}
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), inv_ddqz_z_half2_d, CL_FALSE, 0, sizeof(double)*size_inv_ddqz_z_half2,   inv_ddqz_z_half2_tmp,    0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_con_d, CL_FALSE, 0, sizeof(double)*size_w_con, w_con_tmp,0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), w_d, CL_FALSE, 0, sizeof(double)*size_w, w_tmp, 0, NULL, NULL) );

		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		for (int i=0; i<size_inv_ddqz_z_half2; i++)
			assert(inv_ddqz_z_half2[i] == inv_ddqz_z_half2_tmp[i]);
		
		int zero[1] = {0};
		checkresultsCPP_(w_con,w_con_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"part5:w_con");
		checkresultsCPP_(w,w_tmp,
						 zero,nproma,nproma,
						 zero,nlevp1,nlevp1,
						 zero,nblks_c,nblks_c,"part5:w");
		
		
		delete [] inv_ddqz_z_half2_tmp;
		delete [] w_con_tmp;
		delete [] w_tmp;
		delete [] localStart_tmp;
		delete [] localEnd_tmp;
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("part5");
	}
#endif
	
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 0, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 1, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 2, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 3, sizeof(int32_t), (void*)nlevp1) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 4, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 5, sizeof(int32_t), (void*)ntnd) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 6, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 7, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 8, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5, 9, sizeof(cl_mem), (void*)&w_con_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5,10, sizeof(cl_mem), (void*)&w_d) );
	CheckCL( clSetKernelArg(ocl.kernel_phase5,11, sizeof(cl_mem), (void*)&inv_ddqz_z_half2_d) );
	if (*ntnd==1)
	{
		CheckCL( clSetKernelArg(ocl.kernel_phase5,12, sizeof(cl_mem), (void*)&ddt_w_adv1_d) );
	}
	else if (*ntnd==2)
	{
		CheckCL( clSetKernelArg(ocl.kernel_phase5,12, sizeof(cl_mem), (void*)&ddt_w_adv2_d) );
	}
	else if (*ntnd==3)
	{
		CheckCL( clSetKernelArg(ocl.kernel_phase5,12, sizeof(cl_mem), (void*)&ddt_w_adv3_d) );
	}
	else abort();
	
	//cl_event ev_wait_list[6] = { ev_write_localStart, ev_write_localEnd, ev_write_w_con, ev_write_w, ev_write_inv_ddqz_z_half2, ev_write_ddt_w_adv };
	//CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase5, 3, NULL, global_phase5, NULL, 6, ev_wait_list, &ev_kernel_phase5) );
	if (ocl.bFirstIter)
	{
		cl_event ev_wait_list[6] = { ev_write_localStart, ev_write_localEnd, ev_write_inv_ddqz_z_half2, ev_write_ddt_w_adv1, ev_write_ddt_w_adv2, ev_write_ddt_w_adv3 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase5, 3, NULL, global_phase5, local_phase5, 6, ev_wait_list, &ev_kernel_phase5) );
#ifndef TURBO_MODE
	}
	else if (*istep==1)
	{
		cl_event ev_wait_list[3] = { ev_write_ddt_w_adv1, ev_write_ddt_w_adv2, ev_write_ddt_w_adv3 };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase5, 3, NULL, global_phase5, local_phase5, 3, ev_wait_list, &ev_kernel_phase5) );
#endif
	}
	else
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_phase5, 3, NULL, global_phase5, local_phase5, 0, NULL, &ev_kernel_phase5) );
	
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
#ifndef TURBO_MODE
	if (*istep!=1)
	{
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv1_d, CL_FALSE, 0, sizeof(cl_double)*size_ddt_w_adv, ddt_w_adv, 1, &ev_kernel_phase5, &ev_read_ddt_w_adv1) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv2_d, CL_FALSE, 0, sizeof(cl_double)*size_ddt_w_adv, &ddt_w_adv[size_ddt_w_adv], 1, &ev_kernel_phase5, &ev_read_ddt_w_adv2) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv3_d, CL_FALSE, 0, sizeof(cl_double)*size_ddt_w_adv, &ddt_w_adv[2*size_ddt_w_adv], 1, &ev_kernel_phase5, &ev_read_ddt_w_adv3) );
	}
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv1_d, CL_FALSE, 0, sizeof(cl_double)*size_ddt_w_adv, ddt_w_adv, 1, &ev_kernel_phase5, &ev_read_ddt_w_adv1) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv2_d, CL_FALSE, 0, sizeof(cl_double)*size_ddt_w_adv, &ddt_w_adv[size_ddt_w_adv], 1, &ev_kernel_phase5, &ev_read_ddt_w_adv2) );
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ddt_w_adv3_d, CL_FALSE, 0, sizeof(cl_double)*size_ddt_w_adv, &ddt_w_adv[2*size_ddt_w_adv], 1, &ev_kernel_phase5, &ev_read_ddt_w_adv3) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
};

#pragma mark operators
extern "C" void ocl_div_tri_(double * vec_e, int * nameIdx_in, bool * bOnGPU_in,
							 int * istep,
							 double * geofac_div, /*p_int%geofac_div*/
							 int * iidx, /*p_patch%cells%edge_idx*/
							 int * iblk, /*p_patch%cells%edge_blk*/
							 int * i_startblk,
							 int * i_endblk,
							 int * nproma,
							 int * nlev,
							 int * nblks_c,
							 int * nblks_e,
							 double * div_vec_c, int * nameIdx_out, bool * bOnGPU_out,
							 int * localStart,
							 int * localEnd)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_div_tri");
#endif
	
	cl_event ev_write_div_vec_c, ev_write_vec_e, ev_write_geofac_div, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_kernel, ev_read_div_vec_c;
	
	const int i_startblkm1 = *i_startblk - 1;
	const int nBlocks = *i_endblk-i_startblkm1;
	
	const int i_cell_type = 3;
	const int size_vec_e = (*nproma)*(*nlev)*(*nblks_e);
	const int size_geofac_div = (*nproma)*i_cell_type*(*nblks_c);
	const int size_iidx = (*nproma)*(*nblks_c)*3;
	const int size_iblk = (*nproma)*(*nblks_c)*3;
	const int size_div_vec_c = (*nproma)*(*nlev)*(*nblks_c);
	
	size_t local_size3D[3] = { 1, *nlev, 1 };
	size_t global_size3D[3] = { nBlocks, *nlev, *nproma };
	//cout << global_size3D[0] << " " << global_size3D[1] << " " << global_size3D[2] << endl;
	
	//cout << "Buffers... ";
	cl_mem geofac_div_d = resourcesCL.getBuffers("geofac_div",     1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_geofac_div)[0];
	cl_mem iidx_d =       resourcesCL.getBuffers("ieidx",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iidx)[0];
	cl_mem iblk_d =       resourcesCL.getBuffers("ieblk",          1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iblk)[0];
	cl_mem localStart_d = resourcesCL.getBuffers("localStartCAvg", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*nBlocks)[0];
	cl_mem localEnd_d =   resourcesCL.getBuffers("localEndCAvg",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*nBlocks)[0];
	cl_mem vec_e_d =      resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_in].c_str(),  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vec_e)[0];
	cl_mem div_vec_c_d =  resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_out].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_div_vec_c)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
#ifdef TEST
	double * vec_e_tmp = new double[size_vec_e];
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vec_e_d, CL_TRUE, 0, sizeof(double)*size_vec_e, vec_e_tmp, 0, NULL, NULL) );
	
	int zero[1] = {0};
	checkresultsCPP_(vec_e,vec_e_tmp,
					 zero,nproma,nproma,
					 zero,nlev,nlev,
					 zero,nblks_e,nblks_e,"div:vec_e");
#endif
	
	//cout << "Writing Buffers... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_in && *istep==1)
#else
	if (!*bOnGPU_in && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vec_e_d,      CL_FALSE, 0, sizeof(cl_double)*size_vec_e,      vec_e,      0, NULL, &ev_write_vec_e) );
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), geofac_div_d, CL_FALSE, 0, sizeof(cl_double)*size_geofac_div, geofac_div, 0, NULL, &ev_write_geofac_div) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iidx_d,       CL_FALSE, 0, sizeof(cl_int)*size_iidx,          iidx,       0, NULL, &ev_write_iidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iblk_d,       CL_FALSE, 0, sizeof(cl_int)*size_iblk,          iblk,       0, NULL, &ev_write_iblk) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(cl_int)*nBlocks,            localStart, 0, NULL, &ev_write_localStart) );
		//CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(cl_int)*nBlocks,            localEnd,   0, NULL, &ev_write_localEnd) );
	}
#ifndef TURBO_MODE
	if(!*bOnGPU_out && *istep==1)
#else
	if(!*bOnGPU_out && ocl.bFirstIter)
#endif
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), div_vec_c_d, CL_FALSE, 0, sizeof(cl_double)*size_div_vec_c, div_vec_c, 0, NULL, &ev_write_div_vec_c) );
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_div_tri...\n";
		int * iidx_tmp          = new int[size_iidx];
		int * iblk_tmp          = new int[size_iblk];
		int * localStart_tmp    = new int[nBlocks];
		int * localEnd_tmp      = new int[nBlocks];
		double * geofac_div_tmp = new double[size_geofac_div];
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iidx_d, CL_FALSE, 0, sizeof(int)*size_iidx, iidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iblk_d, CL_FALSE, 0, sizeof(int)*size_iblk, iblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(int)*nBlocks, localEnd_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), geofac_div_d, CL_FALSE, 0, sizeof(double)*size_geofac_div, geofac_div_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		int zero[1] = {0};
		int one[1] = {1};
		int three[1] = {3};
		int nb[1] = {nBlocks};
		
		checkresultsCPP_(iidx,iidx_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,three,three,"div:iidx");
		checkresultsCPP_(iblk,iblk_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,three,three,"div:iblk");
		checkresultsCPP_(localStart,localStart_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"div:localStart");
		checkresultsCPP_(localEnd,localEnd_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"div:localEnd");
		checkresultsCPP_(geofac_div,geofac_div_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_c,nblks_c,"div:geofac_div");
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("div");
	}
#endif
	
	//cout << "Executing Kernel... ";
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 0, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 1, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 2, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 4, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 5, sizeof(cl_mem), (void*)&vec_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 6, sizeof(cl_mem), (void*)&geofac_div_d) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 7, sizeof(cl_mem), (void*)&iidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 8, sizeof(cl_mem), (void*)&iblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri, 9, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri,10, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_div_tri,11, sizeof(cl_mem), (void*)&div_vec_c_d) );
	
	if (ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[5] = { ev_write_vec_e, ev_write_geofac_div, ev_write_iidx, ev_write_iblk, ev_write_div_vec_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 5, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[4] = { ev_write_geofac_div, ev_write_iidx, ev_write_iblk, ev_write_div_vec_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 4, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[4] = { ev_write_vec_e, ev_write_geofac_div, ev_write_iidx, ev_write_iblk };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 4, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[3] = { ev_write_geofac_div, ev_write_iidx, ev_write_iblk };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 3, ev_wait_list, &ev_kernel) );
	}
#ifndef TURBO_MODE
	else if (!ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[2] = { ev_write_vec_e, ev_write_div_vec_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 2, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_div_vec_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_vec_e };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if ((!ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out) || *istep!=1)
#else
	else
#endif
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_div_tri, 3, NULL, global_size3D, local_size3D, 0, NULL, &ev_kernel) );
		
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	//cout << "Downloading Data... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_out && *istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), div_vec_c_d, CL_FALSE, 0, sizeof(cl_double)*size_div_vec_c, div_vec_c, 1, &ev_kernel, &ev_read_div_vec_c) );
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), div_vec_c_d, CL_FALSE, 0, sizeof(cl_double)*size_div_vec_c, div_vec_c, 1, &ev_kernel, &ev_read_div_vec_c) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}


extern "C" void ocl_rot_vertex_atmos_tri_(double * vec_e /*ok*/, int * nameIdx_in, bool * bOnGPU_in,
										  int * istep,
										  double * geofac_rot /*ptr_int%geofac_rot*/,
										  int * iidx /*ptr_patch%verts%edge_idx*/,
										  int * iblk /*ptr_patch%verts%edge_blk*/,
										  int * i_startblk /*ptr_patch%verts%start_blk(rl_start,1)*/,
										  int * nproma /*ok*/,
										  int * slev /*ok*/,
										  int * elev /*ok*/,
										  int * nlev /*ok*/,
										  int * nblks_v /*ptr_patch%nblks_int_v*/,
										  int * nblks_e /*ptr_patch%nblks_int_e*/,
										  double * rot_vec /*ok*/, int * nameIdx_out, bool * bOnGPU_out,
										  int * localStart /*fill with get_indices_v*/,
										  int * localEnd /*fill with get_indices_v*/)
{

#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_rot_vertex_atmos_tri");
#endif
	
	cl_event ev_write_rot_vec, ev_write_vec_e, ev_write_geofac_rot, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_kernel, ev_read_rot_vec;
	
	const int i_size = *nproma;
	const int i_startblkm1 = *i_startblk - 1;
	
	const int i_cell_type = 3;
	const int size_vec_e = (*nproma)*(*nlev)*(*nblks_e);
	const int size_geofac_rot = (*nproma)*6*(*nblks_v);
	const int size_iidx = (*nproma)*(*nblks_v)*6;
	const int size_iblk = (*nproma)*(*nblks_v)*6;
	const int size_rot_vec = (*nproma)*(*nlev)*(*nblks_v);
	
	size_t local_size3D[3] = { 1, *nlev, 1 };
	size_t global_size3D[3] = { *nblks_v, *nlev, i_size };
	//cout << global_size3D[0] << " " << global_size3D[1] << " " << global_size3D[2] << endl;
	
	//cout << "Buffers... ";
	cl_mem geofac_rot_d = resourcesCL.getBuffers("geofac_rot",    1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_geofac_rot)[0];
	cl_mem iidx_d =       resourcesCL.getBuffers("iidxRBF",       1, oclhelper.getContext(), CL_MEM_READ_ONLY, sizeof(int)*size_iidx)[0];
	cl_mem iblk_d =       resourcesCL.getBuffers("iblkRBF",       1, oclhelper.getContext(), CL_MEM_READ_ONLY, sizeof(int)*size_iblk)[0];
	cl_mem localStart_d = resourcesCL.getBuffers("localStartRot", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*(*nblks_v))[0];
	cl_mem localEnd_d =   resourcesCL.getBuffers("localEndRot",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*(*nblks_v))[0];
	cl_mem vec_e_d =      resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_in].c_str(),  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_vec_e)[0];
	cl_mem rot_vec_d =    resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_out].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_rot_vec)[0];
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif

#ifdef TEST
	if (*bOnGPU_in)
	{
		double * vec_e_tmp = new double[size_vec_e];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), vec_e_d, CL_TRUE, 0, sizeof(double)*size_vec_e, vec_e_tmp, 0, NULL, NULL) );
		
		int zero[1] = {0};
		checkresultsCPP_(vec_e,vec_e_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"rot:vec_e");
	}
#endif
	
	//cout << "Writing Buffers... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_in && *istep==1)
#else
	if (!*bOnGPU_in && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), vec_e_d,      CL_FALSE, 0, sizeof(cl_double)*size_vec_e,      vec_e,      0, NULL, &ev_write_vec_e) );
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), geofac_rot_d, CL_FALSE, 0, sizeof(cl_double)*size_geofac_rot, geofac_rot, 0, NULL, &ev_write_geofac_rot) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iidx_d,       CL_FALSE, 0, sizeof(cl_int)*size_iidx,          iidx,       0, NULL, &ev_write_iidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iblk_d,       CL_FALSE, 0, sizeof(cl_int)*size_iblk,          iblk,       0, NULL, &ev_write_iblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(cl_int)*(*nblks_v),         localStart, 0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(cl_int)*(*nblks_v),         localEnd,   0, NULL, &ev_write_localEnd) );
	}
#ifndef TURBO_MODE
	if(!*bOnGPU_out && *istep==1)
#else
	if(!*bOnGPU_out && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), rot_vec_d, CL_FALSE, 0, sizeof(cl_double)*size_rot_vec, rot_vec, 0, NULL, &ev_write_rot_vec) );
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_rot_vertex_atmos_tri...\n";
		int * iidx_tmp          = new int[size_iidx];
		int * iblk_tmp          = new int[size_iblk];
		int * localStart_tmp    = new int[*nblks_v];
		int * localEnd_tmp      = new int[*nblks_v];
		double * geofac_rot_tmp = new double[size_geofac_rot];
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iidx_d, CL_FALSE, 0, sizeof(int)*size_iidx, iidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iblk_d, CL_FALSE, 0, sizeof(int)*size_iblk, iblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(int)*(*nblks_v), localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(int)*(*nblks_v), localEnd_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), geofac_rot_d, CL_FALSE, 0, sizeof(double)*size_geofac_rot, geofac_rot_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		int zero[1] = {0};
		int one[1] = {1};
		int six[1] = {6};
		
		checkresultsCPP_(iidx,iidx_tmp,
						 zero,nproma,nproma,
						 zero,nblks_v,nblks_v,
						 zero,six,six,"rot:iidx");
		checkresultsCPP_(iblk,iblk_tmp,
						 zero,nproma,nproma,
						 zero,nblks_v,nblks_v,
						 zero,six,six,"rot:iblk");
		checkresultsCPP_(localStart,localStart_tmp,
						 zero,nblks_v,nblks_v,
						 zero,one,one,
						 zero,one,one,"rot:localStart");
		checkresultsCPP_(localEnd,localEnd_tmp,
						 zero,nblks_v,nblks_v,
						 zero,one,one,
						 zero,one,one,"rot:localEnd");
		checkresultsCPP_(geofac_rot,geofac_rot_tmp,
						 zero,nproma,nproma,
						 zero,six,six,
						 zero,nblks_v,nblks_v,"rot:geofac_rot");
	}
#endif	

#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("rot");
	}
#endif
	
	//cout << "Executing Kernel... ";
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 0, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 1, sizeof(int32_t), (void*)nblks_v) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 2, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 3, sizeof(int32_t), (void*)slev) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 4, sizeof(int32_t), (void*)elev) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 5, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 6, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 7, sizeof(int32_t), (void*)&i_size) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 8, sizeof(cl_mem), (void*)&vec_e_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri, 9, sizeof(cl_mem), (void*)&geofac_rot_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri,10, sizeof(cl_mem), (void*)&iidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri,11, sizeof(cl_mem), (void*)&iblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri,12, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri,13, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri,14, sizeof(cl_mem), (void*)&rot_vec_d) );

	if (ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[7] = { ev_write_geofac_rot, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 5, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_vec_e, ev_write_geofac_rot, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_geofac_rot, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_rot_vec };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[7] = { ev_write_vec_e, ev_write_geofac_rot, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_rot_vec };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 7, ev_wait_list, &ev_kernel) );
	}
#ifndef TURBO_MODE
	else if (!ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[2] = { ev_write_vec_e, ev_write_rot_vec };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 2, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_rot_vec };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_vec_e };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if ((!ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out) || *istep!=1)
#else
	else
#endif
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rot_vertex_atmos_tri, 3, NULL, global_size3D, local_size3D, 0, NULL, &ev_kernel) );
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	//cout << "Downloading Data... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_out && *istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rot_vec_d, CL_FALSE, 0, sizeof(cl_double)*size_rot_vec, rot_vec, 1, &ev_kernel, &ev_read_rot_vec) );
#endif
	
#ifdef TEST
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	if (*bOnGPU_out || *istep==1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), rot_vec_d, CL_FALSE, 0, sizeof(cl_double)*size_rot_vec, rot_vec, 0, NULL, NULL) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());

	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}

#pragma mark Interpolators
extern "C" void ocl_rbf_vec_interpol_edge_(double * p_vn_in /*ok*/, int * nameIdx_in, bool * bOnGPU_in,
										   int * istep,
										   double * ptr_coeff /*ptr_int%rbf_vec_coeff_e*/,
										   int * iidx /*ptr_int%rbf_vec_idx_e*/,
										   int * iblk /*ptr_int%rbf_vec_blk_e*/,
										   int * i_startblk /*ptr_patch%edges%start_blk(i_rcstartlev,1)*/,
										   int * nproma /*ok*/,
										   int * nlev /*ok*/,
										   int * rbf_vec_dim_e /*from mo_intp_data_strc*/,
										   int * nblks_e /*ptr_patch%nblks_int_e*/,
										   double * p_vt_out /*ok*/, int * nameIdx_out, bool * bOnGPU_out,
										   int * localStart /*fill with get_indices_e*/,
										   int * localEnd /*fill with get_indices_e*/,
										   int * slev /*ok*/,
										   int * elev /*ok*/)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_rbf_vec_interpol_edge");
#endif
	
	cl_event ev_write_p_vt_out, ev_write_p_vn_in, ev_write_ptr_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_kernel, ev_read_p_vt_out;
	
	const int i_startblkm1 = *i_startblk - 1;
	const int nBlocks = (*nblks_e) - i_startblkm1;
	
	const int size_p_vn_in = (*nproma)*(*nlev)*(*nblks_e);
	const int size_ptr_coeff = (*rbf_vec_dim_e)*(*nproma)*(*nblks_e);
	const int size_iidx = (*rbf_vec_dim_e)*(*nproma)*(*nblks_e);
	const int size_iblk = (*rbf_vec_dim_e)*(*nproma)*(*nblks_e);
	const int size_p_vt_out = (*nproma)*(*nlev)*(*nblks_e);
	
	size_t local_size3D[3] = { 1, *nlev, 1 };
	size_t global_size3D[3] = { nBlocks, *nlev, *nproma };
	// cout << global_size3D[0] << " " << global_size3D[1] << " " << global_size3D[2] << endl;
	
	// cout << "Buffers... ";
	cl_mem ptr_coeff_d  = resourcesCL.getBuffers("ptr_coeff",     1, oclhelper.getContext(), CL_MEM_READ_ONLY, sizeof(double)*size_ptr_coeff)[0];
	cl_mem iidx_d       = resourcesCL.getBuffers("iidxRBF",       1, oclhelper.getContext(), CL_MEM_READ_ONLY, sizeof(int)*size_iidx)[0];
	cl_mem iblk_d       = resourcesCL.getBuffers("iblkRBF",       1, oclhelper.getContext(), CL_MEM_READ_ONLY, sizeof(int)*size_iblk)[0];
	cl_mem localStart_d = resourcesCL.getBuffers("localStartRBF", 1, oclhelper.getContext(), CL_MEM_READ_ONLY, sizeof(int)*nBlocks)[0];
	cl_mem localEnd_d   = resourcesCL.getBuffers("localEndRBF",   1, oclhelper.getContext(), CL_MEM_READ_ONLY, sizeof(int)*nBlocks)[0];
	cl_mem p_vn_in_d    = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_in].c_str(),  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_p_vn_in)[0];
	cl_mem p_vt_out_d   = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_out].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_p_vt_out)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
#if TEST
	if (*bOnGPU_in)
	{
		double * p_vn_in_tmp = new double[size_p_vn_in];
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_vn_in_d, CL_TRUE, 0, sizeof(double)*size_p_vn_in, p_vn_in_tmp, 0, NULL, NULL) );
		int zero[1] = {0};
		checkresultsCPP_(p_vn_in,p_vn_in_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_e,nblks_e,"rbf:p_vn_in");
	}
#endif
	
	//cout << "Writing Buffers... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_in && *istep==1)
#else
	if (!*bOnGPU_in && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), p_vn_in_d,    CL_FALSE, 0, sizeof(cl_double)*size_p_vn_in,   p_vn_in,    0, NULL, &ev_write_p_vn_in) );
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), ptr_coeff_d,  CL_FALSE, 0, sizeof(cl_double)*size_ptr_coeff, ptr_coeff,  0, NULL, &ev_write_ptr_coeff) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iidx_d,       CL_FALSE, 0, sizeof(cl_int)*size_iidx,         iidx,       0, NULL, &ev_write_iidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iblk_d,       CL_FALSE, 0, sizeof(cl_int)*size_iblk,         iblk,       0, NULL, &ev_write_iblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(cl_int)*nBlocks,           localStart, 0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(cl_int)*nBlocks,           localEnd,   0, NULL, &ev_write_localEnd) );
	}
#ifndef TURBO_MODE
	if(!*bOnGPU_out && *istep==1)
#else
	if(!*bOnGPU_out && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), p_vt_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_vt_out, p_vt_out, 0, NULL, &ev_write_p_vt_out) );
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_rbf_vec_interpol_edge...\n";
		int * iidx_tmp          = new int[size_iidx];
		int * iblk_tmp          = new int[size_iblk];
		int * localStart_tmp    = new int[nBlocks];
		int * localEnd_tmp      = new int[nBlocks];
		double * ptr_coeff_tmp = new double[size_ptr_coeff];
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iidx_d, CL_FALSE, 0, sizeof(int)*size_iidx, iidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iblk_d, CL_FALSE, 0, sizeof(int)*size_iblk, iblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(int)*nBlocks, localEnd_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), ptr_coeff_d, CL_FALSE, 0, sizeof(double)*size_ptr_coeff, ptr_coeff_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		int zero[1] = {0};
		int one[1] = {1};
		int nb[1] = {nBlocks};
		
		checkresultsCPP_(iidx,iidx_tmp,
						 zero,rbf_vec_dim_e,rbf_vec_dim_e,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,"rbf:iidx");
		checkresultsCPP_(iblk,iblk_tmp,
						 zero,rbf_vec_dim_e,rbf_vec_dim_e,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,"rbf:iblk");
		checkresultsCPP_(localStart,localStart_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"rbf:localStart");
		checkresultsCPP_(localEnd,localEnd_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"rbf:localEnd");
		checkresultsCPP_(ptr_coeff,ptr_coeff_tmp,
						 zero,rbf_vec_dim_e,rbf_vec_dim_e,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,"rbf:ptr_coeff");
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("rbf");
	}
#endif
	
	int argIdx = 0;
	// cout << "Executing Kernel... ";
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(int32_t), (void*)rbf_vec_dim_e) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(int32_t), (void*)slev) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(int32_t), (void*)elev) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(cl_mem), (void*)&p_vn_in_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(cl_mem), (void*)&ptr_coeff_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(cl_mem), (void*)&iidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(cl_mem), (void*)&iblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_rbf_vec_interpol_edge, argIdx++, sizeof(cl_mem), (void*)&p_vt_out_d) );

	if (ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[7] = { ev_write_p_vn_in, ev_write_ptr_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_p_vt_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 7, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_ptr_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_p_vt_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_p_vn_in, ev_write_ptr_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[5] = { ev_write_ptr_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 5, ev_wait_list, &ev_kernel) );
	}
#ifndef TURBO_MODE
	else if (!ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[2] = { ev_write_p_vn_in, ev_write_p_vt_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 2, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_p_vt_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_p_vn_in };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if ((!ocl.bFirstIter &&  bOnGPU_in && *bOnGPU_out) || *istep!=1)
#else
	else
#endif
	  CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_rbf_vec_interpol_edge, 3, NULL, global_size3D, local_size3D, 0, NULL, &ev_kernel) );

#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	//cout << "Downloading Data... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_out && *istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_vt_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_vt_out, p_vt_out, 1, &ev_kernel, &ev_read_p_vt_out) );
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_vt_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_vt_out, p_vt_out, 1, &ev_kernel, &ev_read_p_vt_out) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}

extern "C" void ocl_cells2edges_scalar_(double * p_cell_in, int * nameIdx_in, bool * bOnGPU_in,
										int * istep,
										double * c_int,
										int * iidx, int * iblk,
										int * localStart, int * localEnd, int * nproma,
										int * slev, int * elev, int * nlev,
										int * i_startblk, int * i_endblk, int * nblks_c, int * nblks_e,
										double * p_edge_out, int * nameIdx_out, bool * bOnGPU_out)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_cells2edges_scalar");
#endif
	
	const int i_startblkm1 = *i_startblk - 1;
	const int slevm1 = *slev - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	const int nLevels = *elev - slevm1;
	
	cl_event ev_write_p_edge_out, ev_write_p_cell_in, ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_kernel, ev_read_p_edge_out;
	
	const int size_p_cell_in  = (*nproma)*(*nlev)*(*nblks_c);
	const int size_c_int      = (*nproma)*2*(*nblks_e);
	const int size_iidx       = (*nproma)*(*nblks_e)*2;
	const int size_iblk       = (*nproma)*(*nblks_e)*2;
	const int size_p_edge_out = (*nproma)*(*nlev)*(*nblks_e);
	
	size_t local_size3D[3] = { 1, nLevels, 1 };
	
	size_t global_size3D[3] = { nBlocks, nLevels, *nproma };
	//cout << global_size3D[0] << " " << global_size3D[1] << " " << global_size3D[2] << endl;
	
	//cout << "Buffers... ";
	cl_mem c_int_d      = resourcesCL.getBuffers("c_intC2E",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_c_int)[0];
	cl_mem iidx_d       = resourcesCL.getBuffers("iidxC2E",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iidx)[0];
	cl_mem iblk_d       = resourcesCL.getBuffers("iblkC2E",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iblk)[0];
	cl_mem localStart_d = resourcesCL.getBuffers("localStartC2E", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*(*nblks_e))[0];
	cl_mem localEnd_d   = resourcesCL.getBuffers("localEndC2E",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*(*nblks_e))[0];
	cl_mem p_cell_in_d  = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_in].c_str(),  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_p_cell_in)[0];
	cl_mem p_edge_out_d = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_out].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_p_edge_out)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
	//cout << "Writing Buffers... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_in && *istep==1)
#else
	if (!*bOnGPU_in && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), p_cell_in_d,  CL_FALSE, 0, sizeof(cl_double)*size_p_cell_in, p_cell_in,  0, NULL, &ev_write_p_cell_in) );
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), c_int_d,      CL_FALSE, 0, sizeof(cl_double)*size_c_int,     c_int,      0, NULL, &ev_write_c_int) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iidx_d,       CL_FALSE, 0, sizeof(cl_int)*size_iidx,         iidx,       0, NULL, &ev_write_iidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iblk_d,       CL_FALSE, 0, sizeof(cl_int)*size_iblk,         iblk,       0, NULL, &ev_write_iblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(cl_int)*(*nblks_e),        localStart, 0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(cl_int)*(*nblks_e),        localEnd,   0, NULL, &ev_write_localEnd) );
	}
#ifndef TURBO_MODE
	if(!*bOnGPU_out && *istep==1)
#else
	if(!*bOnGPU_out && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), p_edge_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_edge_out, p_edge_out,  0, NULL, &ev_write_p_edge_out) );
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_cells2edges_scalar...\n";
		int * iidx_tmp          = new int[size_iidx];
		int * iblk_tmp          = new int[size_iblk];
		int * localStart_tmp    = new int[nBlocks];
		int * localEnd_tmp      = new int[nBlocks];
		double * c_int_tmp = new double[size_c_int];
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iidx_d, CL_FALSE, 0, sizeof(int)*size_iidx, iidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iblk_d, CL_FALSE, 0, sizeof(int)*size_iblk, iblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(int)*nBlocks, localEnd_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), c_int_d, CL_FALSE, 0, sizeof(double)*size_c_int, c_int_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		int zero[1] = {0};
		int one[1] = {1};
		int two[1] = {2};
		int three[1] = {3};
		int nb[1] = {nBlocks};
		
		checkresultsCPP_(iidx,iidx_tmp,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,
						 zero,two,two,"c2e:iidx");
		checkresultsCPP_(iblk,iblk_tmp,
						 zero,nproma,nproma,
						 zero,nblks_e,nblks_e,
						 zero,two,two,"c2e:iblk");
		checkresultsCPP_(localStart,localStart_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"c2e:localStart");
		checkresultsCPP_(localEnd,localEnd_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"c2e:localEnd");
		checkresultsCPP_(c_int,c_int_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_c,nblks_c,"c2e:c_int");
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("c2e");
	}
#endif
	
	//cout << "Executing Kernel... ";
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 0, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 1, sizeof(int32_t), (void*)&slevm1) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 2, sizeof(int32_t), (void*)elev) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 4, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 5, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 6, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 7, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 8, sizeof(cl_mem), (void*)&p_cell_in_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar, 9, sizeof(cl_mem), (void*)&c_int_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar,10, sizeof(cl_mem), (void*)&iidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar,11, sizeof(cl_mem), (void*)&iblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar,12, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar,13, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cells2edges_scalar,14, sizeof(cl_mem), (void*)&p_edge_out_d) );
	
	if (ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[7] = { ev_write_p_cell_in, ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_p_edge_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 7, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_p_edge_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_p_cell_in, ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[5] = { ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 5, ev_wait_list, &ev_kernel) );
	}
#ifndef TURBO_MODE
	else if (!ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[2] = { ev_write_p_cell_in, ev_write_p_edge_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 2, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_p_edge_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_p_cell_in };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if ((!ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out) || *istep!=1)
#else
	else
#endif
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cells2edges_scalar, 3, NULL, global_size3D, local_size3D, 0, NULL, &ev_kernel) );
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	//cout << "Downloading Data... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_out && *istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_edge_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_edge_out, p_edge_out, 1, &ev_kernel, &ev_read_p_edge_out) );
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_edge_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_edge_out, p_edge_out, 1, &ev_kernel, &ev_read_p_edge_out) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}

extern "C" void ocl_edges2cells_scalar_tri_(double * p_edge_in /*ok*/, int * nameIdx_in, bool * bOnGPU_in,
											int * istep,
											double * c_int, /*ok*/
											int * iidx, /*p_patch%cells%edge_idx*/
											int * iblk, /*p_patch%cells%edge_idx*/
											int * localStart, /*fill with get_indices_c*/
											int * localEnd, /*fill with get_indices_c*/
											int * nproma, /*ok*/
											int * slev, /*ok*/
											int * elev, /*ok*/
											int * nlev, /*ok*/
											int * i_startblk, /*p_patch%cells%start_blk(rl_start,1)*/
											int * i_endblk, /*p_patch%cells%end_blk(rl_end,MAX(1,p_patch%n_childdom))*/
											int * nblks_c, /*p_patch%nblks_int_c*/
											int * nblks_e, /*p_patch%nblks_int_e*/
											double * p_cell_out /*ok*/, int * nameIdx_out, bool * bOnGPU_out)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_edges2cells_scalar_tri");
#endif
	
	const int i_cell_type = 3;
	const int i_startblkm1 = *i_startblk - 1;
	const int slevm1 = *slev - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	const int nLevels = *elev - slevm1;
	
	cl_event ev_write_p_cell_out, ev_write_p_edge_in, ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_kernel, ev_read_p_cell_out;
	
	const int size_p_edge_in  = (*nproma)*(*nlev)*(*nblks_e);
	const int size_c_int      = (*nproma)*i_cell_type*(*nblks_c);
	const int size_iidx       = (*nproma)*(*nblks_c)*3;
	const int size_iblk       = (*nproma)*(*nblks_c)*3;
	const int size_p_cell_out = (*nproma)*(*nlev)*(*nblks_c);
	
	size_t local_size3D[3] = { 1, nLevels, 1 };
	size_t global_size3D[3] = { nBlocks, nLevels, *nproma };
	
	//cout << "Buffers... ";
	cl_mem c_int_d      = resourcesCL.getBuffers("c_intE2C",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_c_int)[0];
	cl_mem iidx_d       = resourcesCL.getBuffers("iidxE2C",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iidx)[0];
	cl_mem iblk_d       = resourcesCL.getBuffers("iblkE2C",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iblk)[0];
	cl_mem localStart_d = resourcesCL.getBuffers("localStartE2C", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*nBlocks)[0];
	cl_mem localEnd_d   = resourcesCL.getBuffers("localEndE2C",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*nBlocks)[0];
	cl_mem p_edge_in_d  = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_in].c_str(),  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_p_edge_in)[0];
	cl_mem p_cell_out_d = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_out].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_p_cell_out)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
#ifdef TEST
	double * p_edge_in_tmp = new double[size_p_edge_in];
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_edge_in_d, CL_TRUE, 0, sizeof(double)*size_p_edge_in, p_edge_in_tmp, 0, NULL, NULL) );
	
	int zero[1] = {0};
	checkresultsCPP_(p_edge_in,p_edge_in_tmp,
					 zero,nproma,nproma,
					 zero,nlev,nlev,
					 zero,nblks_e,nblks_e,"e2c:p_edge_in");
#endif
	
	//cout << "Writing Buffers... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_in && *istep==1)
#else
	if (!*bOnGPU_in && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), p_edge_in_d,  CL_FALSE, 0, sizeof(cl_double)*size_p_edge_in, p_edge_in,  0, NULL, &ev_write_p_edge_in) );
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), c_int_d,      CL_FALSE, 0, sizeof(cl_double)*size_c_int,     c_int,      0, NULL, &ev_write_c_int) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iidx_d,       CL_FALSE, 0, sizeof(cl_int)*size_iidx,         iidx,       0, NULL, &ev_write_iidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iblk_d,       CL_FALSE, 0, sizeof(cl_int)*size_iblk,         iblk,       0, NULL, &ev_write_iblk) );
	}
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(cl_int)*nBlocks,           localStart, 0, NULL, &ev_write_localStart) );
	CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(cl_int)*nBlocks,           localEnd,   0, NULL, &ev_write_localEnd) );
	
#ifndef TURBO_MODE
	if(!*bOnGPU_out && *istep==1)
#else
	if(!*bOnGPU_out && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), p_cell_out_d,   CL_FALSE, 0, sizeof(cl_double)*size_p_cell_out, p_cell_out,   0, NULL, &ev_write_p_cell_out) );
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_edges2cells_scalar_tri...\n";
		int * iidx_tmp          = new int[size_iidx];
		int * iblk_tmp          = new int[size_iblk];
		int * localStart_tmp    = new int[nBlocks];
		int * localEnd_tmp      = new int[nBlocks];
		double * c_int_tmp = new double[size_c_int];
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iidx_d, CL_FALSE, 0, sizeof(int)*size_iidx, iidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iblk_d, CL_FALSE, 0, sizeof(int)*size_iblk, iblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(int)*nBlocks, localEnd_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), c_int_d, CL_FALSE, 0, sizeof(double)*size_c_int, c_int_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		int zero[1] = {0};
		int one[1] = {1};
		int three[1] = {3};
		int nb[1] = {nBlocks};
		
		checkresultsCPP_(iidx,iidx_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,three,three,"e2c:iidx");
		checkresultsCPP_(iblk,iblk_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,three,three,"e2c:iblk");
		checkresultsCPP_(localStart,localStart_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"e2c:localStart");
		checkresultsCPP_(localEnd,localEnd_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"e2c:localEnd");
		checkresultsCPP_(c_int,c_int_tmp,
						 zero,nproma,nproma,
						 zero,three,three,
						 zero,nblks_c,nblks_c,"e2c:c_int");
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("e2c");
	}
#endif
	
	//cout << "Executing Kernel... ";
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 0, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 1, sizeof(int32_t), (void*)&slevm1) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 2, sizeof(int32_t), (void*)elev) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 4, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 5, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 6, sizeof(int32_t), (void*)nblks_e) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 7, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 8, sizeof(cl_mem), (void*)&p_edge_in_d) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri, 9, sizeof(cl_mem), (void*)&c_int_d) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri,10, sizeof(cl_mem), (void*)&iidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri,11, sizeof(cl_mem), (void*)&iblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri,12, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri,13, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_edges2cells_scalar_tri,14, sizeof(cl_mem), (void*)&p_cell_out_d) );
	
	if (ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[7] = { ev_write_localStart, ev_write_localEnd, ev_write_p_edge_in, ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_p_cell_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 7, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_localStart, ev_write_localEnd, ev_write_c_int, ev_write_iidx, ev_write_iblk, ev_write_p_cell_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_localStart, ev_write_localEnd, ev_write_p_edge_in, ev_write_c_int, ev_write_iidx, ev_write_iblk };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[5] = { ev_write_localStart, ev_write_localEnd, ev_write_c_int, ev_write_iidx, ev_write_iblk };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 5, ev_wait_list, &ev_kernel) );
	}
#ifndef TURBO_MODE
	else if (!ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[2] = { ev_write_p_edge_in, ev_write_p_cell_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 2, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_p_cell_out };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_p_edge_in };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if ((!ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out) || *istep!=1)
#else
	else
#endif
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_edges2cells_scalar_tri, 3, NULL, global_size3D, local_size3D, 0, NULL, &ev_kernel) );
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	//cout << "Downloading Data... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_out && *istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_cell_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_cell_out, p_cell_out, 1, &ev_kernel, &ev_read_p_cell_out) );
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), p_cell_out_d, CL_FALSE, 0, sizeof(cl_double)*size_p_cell_out, p_cell_out, 1, &ev_kernel, &ev_read_p_cell_out) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}


extern "C" void ocl_cell_avg_(double * psi_c /*ok*/, int * nameIdx_in, bool * bOnGPU_in,
							  int * istep,
							  double * avg_coeff, /*ok*/
							  int * iidx, /*p_patch%cells%neighbor_idx*/
							  int * iblk, /*p_patch%cells%neighbor_blk*/
							  int * localStart, /*filled with get_indices_c*/
							  int * localEnd, /*filled with get_indices_c*/
							  int * nproma, /*ok*/
							  int * slev, /*1*/
							  int * elev, /*nlev*/
							  int * nlev, /*ok*/
							  int * i_startblk, /*p_patch%cells%start_blk(4,1)*/
							  int * i_endblk, /*p_patch%cells%end_blk(min_rlcell_int,MAX(1,p_patch%n_childdom))*/
							  int * nblks_c, /*ptr_patch%nblks_int_c*/
							  double * avg_psi_c /*ok*/, int * nameIdx_out, bool * bOnGPU_out)
{
#ifdef REGION_TIMING
	if (bProfiling) profiler.push_start("ocl_cell_avg");
#endif
	
	const int i_cell_type = 3;
	const int i_startblkm1 = *i_startblk - 1;
	const int slevm1 = *slev - 1;
	const int nBlocks = *i_endblk - i_startblkm1;
	const int nLevels = *elev - slevm1;
	
	cl_event ev_write_avg_psi_c, ev_write_psi_c, ev_write_avg_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_kernel, ev_read_avg_psi_c;
	
	const int size_psi_c     = (*nproma)*(*nlev)*(*nblks_c);
	const int size_avg_coeff = (*nproma)*(*nlev)*(*nblks_c);
	const int size_iidx      = (*nproma)*(*nblks_c)*3;
	const int size_iblk      = (*nproma)*(*nblks_c)*3;
	const int size_avg_psi_c = (*nproma)*(*nlev)*(*nblks_c);
	
	size_t local_size3D[3] = { 1, nLevels, 1 };
	size_t global_size3D[3] = { nBlocks, nLevels, *nproma };
	//cout << global_size3D[0] << " " << global_size3D[1] << " " << global_size3D[2] << endl;
	
	//cout << "Buffers... ";
	cl_mem avg_coeff_d  = resourcesCL.getBuffers("avg_coeff",      1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_avg_coeff)[0];
	cl_mem iidx_d       = resourcesCL.getBuffers("iidxCAvg",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iidx)[0];
	cl_mem iblk_d       = resourcesCL.getBuffers("iblkCAvg",       1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*size_iblk)[0];
	cl_mem localStart_d = resourcesCL.getBuffers("localStartCAvg", 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*nBlocks)[0];
	cl_mem localEnd_d   = resourcesCL.getBuffers("localEndCAvg",   1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(int)*nBlocks)[0];
	cl_mem psi_c_d      = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_in].c_str(),  1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_psi_c)[0];
	cl_mem avg_psi_c_d  = resourcesCL.getBuffers(ocl.idx2ObjName[*nameIdx_out].c_str(), 1, oclhelper.getContext(), CL_MEM_READ_WRITE, sizeof(double)*size_avg_psi_c)[0];
	
#ifndef REGION_TIMING
	if (bProfiling) profiler.push_start("H2D");
#endif
	
#ifdef TEST
	double * psi_c_tmp = new double[size_psi_c];
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), psi_c_d, CL_TRUE, 0, sizeof(double)*size_psi_c, psi_c_tmp, 0, NULL, NULL) );
	int zero[1] = {0};
	checkresultsCPP_(psi_c,psi_c_tmp,
					 zero,nproma,nproma,
					 zero,nlev,nlev,
					 zero,nblks_c,nblks_c,"avg:psi_c");
#endif
	
	//cout << "Writing Buffers... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_in && *istep==1)
#else
	if (!*bOnGPU_in && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), psi_c_d,      CL_FALSE, 0, sizeof(cl_double)*size_psi_c,     psi_c,      0, NULL, &ev_write_psi_c) );
	if (ocl.bFirstIter)
	{
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), avg_coeff_d,  CL_FALSE, 0, sizeof(cl_double)*size_avg_coeff, avg_coeff,  0, NULL, &ev_write_avg_coeff) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iidx_d,       CL_FALSE, 0, sizeof(cl_int)*size_iidx,         iidx,       0, NULL, &ev_write_iidx) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), iblk_d,       CL_FALSE, 0, sizeof(cl_int)*size_iblk,         iblk,       0, NULL, &ev_write_iblk) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(cl_int)*nBlocks,           localStart, 0, NULL, &ev_write_localStart) );
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(cl_int)*nBlocks,           localEnd,   0, NULL, &ev_write_localEnd) );
	}
#ifndef TURBO_MODE
	if(!*bOnGPU_out && *istep==1)
#else
	if(!*bOnGPU_out && ocl.bFirstIter)
#endif
		CheckCL( clEnqueueWriteBuffer(oclhelper.getCommandQueue(), avg_psi_c_d, CL_FALSE, 0, sizeof(cl_double)*size_avg_psi_c, avg_psi_c, 0, NULL, &ev_write_avg_psi_c) );
	
#ifdef TEST
	{
		cout << "checking transfers for ocl_cell_avg...\n";
		int * iidx_tmp          = new int[size_iidx];
		int * iblk_tmp          = new int[size_iblk];
		int * localStart_tmp    = new int[nBlocks];
		int * localEnd_tmp      = new int[nBlocks];
		double * avg_coeff_tmp = new double[size_avg_coeff];
		double * psi_c_tmp     = new double[size_psi_c];
		double * avg_psi_c_tmp = new double[size_avg_psi_c];
		
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iidx_d, CL_FALSE, 0, sizeof(int)*size_iidx, iidx_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), iblk_d, CL_FALSE, 0, sizeof(int)*size_iblk, iblk_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localStart_d, CL_FALSE, 0, sizeof(int)*nBlocks, localStart_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), localEnd_d,   CL_FALSE, 0, sizeof(int)*nBlocks, localEnd_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), avg_coeff_d, CL_FALSE, 0, sizeof(double)*size_avg_coeff, avg_coeff_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), psi_c_d,     CL_FALSE, 0, sizeof(double)*size_psi_c,     psi_c_tmp, 0, NULL, NULL) );
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), avg_psi_c_d, CL_FALSE, 0, sizeof(double)*size_avg_psi_c, avg_psi_c_tmp, 0, NULL, NULL) );
		
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		
		int zero[1] = {0};
		int one[1] = {1};
		int three[1] = {3};
		int nb[1] = {nBlocks};
		
		checkresultsCPP_(iidx,iidx_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,three,three,"avg:iidx");
		checkresultsCPP_(iblk,iblk_tmp,
						 zero,nproma,nproma,
						 zero,nblks_c,nblks_c,
						 zero,three,three,"avg:iblk");
		checkresultsCPP_(localStart,localStart_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"avg:localStart");
		checkresultsCPP_(localEnd,localEnd_tmp,
						 zero,nb,nb,
						 zero,one,one,
						 zero,one,one,"avg:localEnd");
		checkresultsCPP_(avg_coeff,avg_coeff_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"avg:avg_coeff");
		checkresultsCPP_(psi_c,psi_c_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"avg:psi_c");
		checkresultsCPP_(avg_psi_c,avg_psi_c_tmp,
						 zero,nproma,nproma,
						 zero,nlev,nlev,
						 zero,nblks_c,nblks_c,"avg:avg_psi_c");
	}
#endif
	
#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profiler.pop_stop();
		profiler.push_start("kernels");
		profilerKernels.push_start("avg");
	}
#endif
	
	//cout << "Executing Kernel... ";
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 0, sizeof(int32_t), (void*)nproma) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 1, sizeof(int32_t), (void*)&slevm1) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 2, sizeof(int32_t), (void*)elev) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 3, sizeof(int32_t), (void*)nlev) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 4, sizeof(int32_t), (void*)&i_startblkm1) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 5, sizeof(int32_t), (void*)i_endblk) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 6, sizeof(int32_t), (void*)nblks_c) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 7, sizeof(cl_mem), (void*)&psi_c_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 8, sizeof(cl_mem), (void*)&avg_coeff_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg, 9, sizeof(cl_mem), (void*)&iidx_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg,10, sizeof(cl_mem), (void*)&iblk_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg,11, sizeof(cl_mem), (void*)&localStart_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg,12, sizeof(cl_mem), (void*)&localEnd_d) );
	CheckCL( clSetKernelArg(ocl.kernel_cell_avg,13, sizeof(cl_mem), (void*)&avg_psi_c_d) );
	
	if (ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[7] = { ev_write_psi_c, ev_write_avg_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_avg_psi_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, local_size3D, 7, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_avg_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd, ev_write_avg_psi_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[6] = { ev_write_psi_c, ev_write_avg_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, local_size3D, 6, ev_wait_list, &ev_kernel) );
	}
	else if (ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[5] = { ev_write_avg_coeff, ev_write_iidx, ev_write_iblk, ev_write_localStart, ev_write_localEnd };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, local_size3D, 5, ev_wait_list, &ev_kernel) );
	}
#ifndef TURBO_MODE
	else if (!ocl.bFirstIter && !*bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[2] = { ev_write_psi_c, ev_write_avg_psi_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, NULL, 2, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && *bOnGPU_in && !*bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_avg_psi_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if (!ocl.bFirstIter && !*bOnGPU_in && *bOnGPU_out && *istep==1)
	{
		cl_event ev_wait_list[1] = { ev_write_psi_c };
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, local_size3D, 1, ev_wait_list, &ev_kernel) );
	}
	else if ((!ocl.bFirstIter && *bOnGPU_in && *bOnGPU_out) || *istep!=1)
#else
	else
#endif
		CheckCL( clEnqueueNDRangeKernel(oclhelper.getCommandQueue(), ocl.kernel_cell_avg, 3, NULL, global_size3D, local_size3D, 0, NULL, &ev_kernel) );

#ifndef REGION_TIMING
	if (bProfiling)
	{
		clFlush(oclhelper.getCommandQueue());
		clFinish(oclhelper.getCommandQueue());
		profilerKernels.pop_stop();
		profiler.pop_stop();
		profiler.push_start("D2H");
	}
#endif
	
	//cout << "Downloading Data... ";
#ifndef TURBO_MODE
	if (!*bOnGPU_out && *istep!=1)
		CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), avg_psi_c_d, CL_FALSE, 0, sizeof(cl_double)*size_avg_psi_c, avg_psi_c, 1, &ev_kernel, &ev_read_avg_psi_c) );
#endif
	
#ifdef TEST
	CheckCL( clEnqueueReadBuffer(oclhelper.getCommandQueue(), avg_psi_c_d, CL_FALSE, 0, sizeof(cl_double)*size_avg_psi_c, avg_psi_c, 1, &ev_kernel, &ev_read_avg_psi_c) );
#endif
	
	clFlush(oclhelper.getCommandQueue());
	clFinish(oclhelper.getCommandQueue());
	
	if (bProfiling)
	{
		profiler.pop_stop();
		//profiler.pop_stop();
	}
}
