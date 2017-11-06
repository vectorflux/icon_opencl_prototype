//===========================================================================
//
//	The purpose of this class is to avoid recalling and recreating
//		OpenCL structures that are needed every call of the methods
//		in the various extern "C" methods called from fortran
//	Currently it is not possible to have buffers as they require size information
//		not available statically
//	A downside is that parameters cannot be passed statically
//
//===========================================================================

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <map>

#include "IncludesCL.h"

//static bool bVerbose = false;

using namespace std;

typedef double Real;

class OpenCL_statics
{
private:
public:
	map<int,string> lutCheck;
	map<int,string> idx2ObjName;
	bool bFirstIter;
	
	cl_command_queue command_queue;
	char * options;
	
	int linesize;
	
	// list of kernel objects
	cl_kernel kernel_div_tri, kernel_div_hex, kernel_div_avg, kernel_rot_vertex_atmos_tri, kernel_rot_vertex_atmos_hex;
	cl_kernel kernel_rbf_vec_interpol_edge, kernel_cells2edges_scalar, kernel_edges2cells_scalar_tri, kernel_edges2cells_scalar_hex, kernel_cell_avg;
	
	cl_kernel kernel_step1ia, kernel_step1ib, kernel_step1ic, kernel_step1a, kernel_step1b, kernel_stepn;
	cl_kernel kernel_phase2;
	cl_kernel kernel_phase3a, kernel_phase3b, kernel_phase3c;
	cl_kernel kernel_phase4;
	cl_kernel kernel_phase5;
	
	cl_kernel kernel_exner1, kernel_exner2a, kernel_exner2b, kernel_exner3a, kernel_exner3b;
	
	cl_kernel kernel_theta1, kernel_theta2, kernel_theta3a, kernel_theta3b;
	
	cl_kernel kernel_grad1, kernel_grad2, kernel_grad3, kernel_grad4;
	
	cl_kernel kernel_vn1a, kernel_vn1b;
	cl_kernel kernel_vn_avg;
	cl_kernel kernel_vn_half1, kernel_vn_half2, kernel_vn_half3, kernel_vn_half4, kernel_vn_half5;

	cl_kernel kernel_fluxes;
	
	cl_kernel kernel_endphase1a, kernel_endphase1b, kernel_endphase1c, kernel_endphase2a, kernel_endphase2b, kernel_endphase3, kernel_endphase4;
	cl_kernel kernel_endphase5a, kernel_endphase5b, kernel_endphase5c;
	cl_kernel kernel_endphase6, kernel_endphase7, kernel_endphase8;
	cl_kernel kernel_endphase9a, kernel_endphase9b, kernel_endphase9c, kernel_endphase9d;
	
	cl_kernel kernel_init;
	
	
	// list of memory objects
	
	// the two maps contain the associations between the data structures
	//	and the addresses of the memory objects on host side
	//	NOT safe for localStart, localEnd as they are allocated and deallocated
	// if there is no alloc/deallc of the arrays in between calls to C++ side,
	//	it should, in principle, be safe
	// asynchronicity of transfers might still cause problems: transfer events need caching as well
	//map<cl_mem, int *>      idxCache;
	//map<cl_mem, double *>   dataCache;
	//map<cl_mem, cl_event *> transferCache;
	
	OpenCL_statics()
	{
		lutCheck[ 0] = "rbf";
		lutCheck[ 1] = "rot";
		lutCheck[ 2] = "p1a";
		lutCheck[ 3] = "p1b";
		lutCheck[ 4] = "p1d";
		lutCheck[ 5] = "e2c";
		lutCheck[ 6] = "p2a";
		lutCheck[ 7] = "p2b";
		lutCheck[ 8] = "p3a";
		lutCheck[ 9] = "p3b";
		lutCheck[10] = "p3c";
		lutCheck[11] = "avg";
		lutCheck[12] = "p4";
		lutCheck[13] = "p5";
		lutCheck[14] = "c2e";
		lutCheck[15] = "z_exner_ex_pr";
		lutCheck[16] = "z_exner_pr";
		lutCheck[17] = "exner_old";
		lutCheck[18] = "z_dexner_dz_c";
		lutCheck[19] = "rho_ic";
		lutCheck[20] = "theta_v_h";
		lutCheck[21] = "z_th_ddz_exner_c";
		lutCheck[22] = "z_gradh_exner";
		lutCheck[23] = "z_hydro_corr";
		lutCheck[24] = "vn";
		lutCheck[25] = "vn_avg";
		lutCheck[26] = "vn_half";
		lutCheck[27] = "z_concorr_e";
		lutCheck[28] = "mass_fl";
		lutCheck[29] = "theta_fl";
		lutCheck[30] = "div";
		lutCheck[31] = "w";
		lutCheck[32] = "rho";
		lutCheck[33] = "exner";
		lutCheck[34] = "rhotheta";
		lutCheck[35] = "theta_v";
		lutCheck[36] = "ddt_w_adv1";
		lutCheck[37] = "ddt_w_adv2";
		lutCheck[38] = "ddt_w_adv3";
		
		idx2ObjName[0] = "vn";
		idx2ObjName[1] = "vt";
		idx2ObjName[2] = "omega_z";
		idx2ObjName[3] = "z_kin_hor_e";
		idx2ObjName[4] = "e_kinh";
		idx2ObjName[5] = "z_concorr_e";
		idx2ObjName[6] = "w_concorr_c";
		idx2ObjName[7] = "z_hadv_w";
		//idx2ObjName[8] = "ddt_w_adv_2";
		idx2ObjName[9] = "rho_nnow";
		idx2ObjName[10] = "z_rho_e";
		idx2ObjName[11] = "theta_v_nnow";
		idx2ObjName[12] = "z_theta_v_e";
		idx2ObjName[13] = "vn_nnew"; // is this the same as 0? depends
		idx2ObjName[14] = "mass_fl_e";
		idx2ObjName[15] = "z_mass_fl_div";
		idx2ObjName[16] = "z_theta_v_fl_e";
		idx2ObjName[17] = "z_theta_v_fl_div";
		idx2ObjName[18] = "vn_nnew";
		idx2ObjName[19] = "vn_nnow";
		idx2ObjName[20] = "w_nnew";
		idx2ObjName[21] = "w_nnow";
		idx2ObjName[22] = "theta_v_nnew";
		idx2ObjName[23] = "theta_v_nnow";
		idx2ObjName[24] = "rho_nnew";
		idx2ObjName[25] = "rho_nnow";
		idx2ObjName[26] = "ddt_w_adv1";
		idx2ObjName[27] = "ddt_w_adv2";
		idx2ObjName[28] = "ddt_w_adv3";
		
		bFirstIter = true;
		
		command_queue = oclhelper.getCommandQueue();
		options = new char[256];
		//sprintf(options, " –cl-fast-relaxed-math ");
		sprintf(options, "");
		oclhelper.readProgram("/users/wsawyer/project/Software/ICON/icon_testbed_base.v0.1_1_opencl/src/OpenCL/kernels_operators.cl", options);
		oclhelper.readProgram("/users/wsawyer/project/Software/ICON/icon_testbed_base.v0.1_1_opencl/src/OpenCL/kernels_interpolation.cl", options);
		oclhelper.readProgram("/users/wsawyer/project/Software/ICON/icon_testbed_base.v0.1_1_opencl/src/OpenCL/kernels_velocity_tendencies.cl", options);
		oclhelper.readProgram("/users/wsawyer/project/Software/ICON/icon_testbed_base.v0.1_1_opencl/src/OpenCL/kernels_solve_nh.cl", options);
		
		kernel_div_tri = oclhelper.getKernel("div_tri");
		kernel_div_avg = oclhelper.getKernel("div_avg");
		kernel_rot_vertex_atmos_tri = oclhelper.getKernel("rot_vertex_atmos_tri");
		kernel_rot_vertex_atmos_hex = oclhelper.getKernel("rot_vertex_atmos_hex");
		
		kernel_rbf_vec_interpol_edge  = oclhelper.getKernel("rbf_vec_interpol_edge");
		kernel_cells2edges_scalar     = oclhelper.getKernel("cells2edges_scalar");
		kernel_edges2cells_scalar_tri = oclhelper.getKernel("edges2cells_scalar_tri");
		kernel_edges2cells_scalar_hex = oclhelper.getKernel("edges2cells_scalar_hex");
		kernel_cell_avg               = oclhelper.getKernel("cell_avg");
		
		kernel_step1ia = oclhelper.getKernel("step1ia");
		kernel_step1ib = oclhelper.getKernel("step1ib");
		kernel_step1ic = oclhelper.getKernel("step1ic");
		kernel_step1a  = oclhelper.getKernel("step1a");
		kernel_step1b  = oclhelper.getKernel("step1b");
		kernel_stepn   = oclhelper.getKernel("stepn");
		kernel_phase2  = oclhelper.getKernel("phase2");
		kernel_phase3a = oclhelper.getKernel("phase3a");
		kernel_phase3b = oclhelper.getKernel("phase3b");
		kernel_phase3c = oclhelper.getKernel("phase3c");
		kernel_phase4  = oclhelper.getKernel("phase4");
		kernel_phase5  = oclhelper.getKernel("phase5");
		
		kernel_exner1  = oclhelper.getKernel("exner1");
		kernel_exner2a = oclhelper.getKernel("exner2a");
		kernel_exner2b = oclhelper.getKernel("exner2b");
		kernel_exner3a = oclhelper.getKernel("exner3a");
		kernel_exner3b = oclhelper.getKernel("exner3b");
		
		kernel_theta1  = oclhelper.getKernel("theta1");
		kernel_theta2  = oclhelper.getKernel("theta2");
		kernel_theta3a = oclhelper.getKernel("theta3a");
		kernel_theta3b = oclhelper.getKernel("theta3b");
		
		kernel_grad1 = oclhelper.getKernel("grad1");
		kernel_grad2 = oclhelper.getKernel("grad2");
		kernel_grad3 = oclhelper.getKernel("grad3");
		kernel_grad4 = oclhelper.getKernel("grad4");
		
		kernel_vn1a = oclhelper.getKernel("vn1a");
		kernel_vn1b = oclhelper.getKernel("vn1b");
		kernel_vn_avg = oclhelper.getKernel("vn_avg");
		kernel_vn_half1 = oclhelper.getKernel("vn_half1");
		kernel_vn_half2 = oclhelper.getKernel("vn_half2");
		kernel_vn_half3 = oclhelper.getKernel("vn_half3");
		kernel_vn_half4 = oclhelper.getKernel("vn_half4");
		kernel_vn_half5 = oclhelper.getKernel("vn_half5");
		
		kernel_fluxes = oclhelper.getKernel("fluxes");
		
		kernel_endphase1a = oclhelper.getKernel("endphase1a");
		kernel_endphase1b = oclhelper.getKernel("endphase1b");
		kernel_endphase1c = oclhelper.getKernel("endphase1c");
		kernel_endphase2a = oclhelper.getKernel("endphase2a");
		kernel_endphase2b = oclhelper.getKernel("endphase2b");
		kernel_endphase3  = oclhelper.getKernel("endphase3");
		kernel_endphase4  = oclhelper.getKernel("endphase4");
		kernel_endphase5a = oclhelper.getKernel("endphase5a");
		kernel_endphase5b = oclhelper.getKernel("endphase5b");
		kernel_endphase5c = oclhelper.getKernel("endphase5c");
		kernel_endphase6  = oclhelper.getKernel("endphase6");
		kernel_endphase7  = oclhelper.getKernel("endphase7");
		kernel_endphase8  = oclhelper.getKernel("endphase8");
		kernel_endphase9a = oclhelper.getKernel("endphase9a");
		kernel_endphase9b = oclhelper.getKernel("endphase9b");
		kernel_endphase9c = oclhelper.getKernel("endphase9c");
		kernel_endphase9d = oclhelper.getKernel("endphase9d");
		
		kernel_init = oclhelper.getKernel("init");
	}
	
	void initializeMemObjects()
	{
	}
	
	~OpenCL_statics()
	{
	}
};

