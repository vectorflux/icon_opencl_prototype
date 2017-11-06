//==================================================
//
//	See comments in kernels_operators.cl
//
//	the choice of iterators might need some rethinking
//
//==================================================

#ifndef __Tahiti__
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#else
#pragma OPENCL EXTENSION cl_amd_fp64: enable
#endif

#ifndef NLEV_P1
#define NLEV_P1 36
#endif

#ifndef NLEV_
#define NLEV_ 35
#endif


__kernel void rbf_vec_interpol_edge(__private int nproma,
					                __private int rbf_vec_dim_e,
					                __private int nblocks_e,
									__private int slev,
									__private int elev,
									__private int nlev,
					                __private int i_startblk,
					                __global double * p_vn_in,
					                __global double * ptr_coeff,
					                __global int * iidx,
					                __global int * iblk,
					                __global int * localStart,
					                __global int * localEnd,
					                __global double * p_vt_out)
{
// can make use of shared memory here
       int line_idx[4];
       double line_ptr_coeff[4];
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[jb] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_local_id(1) < 4 && jb < nblocks_e && je < localEnd[get_global_id(0)])
	{
		const int idx = je*rbf_vec_dim_e + jb*rbf_vec_dim_e*nproma + get_local_id(1);
		
		line_idx[localIdx] = iidx[idx]-1 + (iblk[idx]-1)*nproma*nlev;
		line_ptr_coeff[localIdx] = ptr_coeff[idx];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
		
	if (jk < elev && jb < nblocks_e && je < localEnd[get_global_id(0)])
	{	
		// index computations
		const int idx_out = je + jk*nproma + jb*nproma*nlev;
		const int idx_p_vn_in1 = line_idx[0] + jk*nproma;
		const int idx_p_vn_in2 = line_idx[1] + jk*nproma;
		const int idx_p_vn_in3 = line_idx[2] + jk*nproma;
		const int idx_p_vn_in4 = line_idx[3] + jk*nproma;
		
		p_vt_out[idx_out] = line_ptr_coeff[0] * p_vn_in[idx_p_vn_in1] +
		                    line_ptr_coeff[1] * p_vn_in[idx_p_vn_in2] +
		                    line_ptr_coeff[2] * p_vn_in[idx_p_vn_in3] +
		                    line_ptr_coeff[3] * p_vn_in[idx_p_vn_in4];
	}
}

__kernel void cells2edges_scalar(__private int nproma,
								 __private int slev,
								 __private int elev,
								 __private int nlev,
								 __private int i_startblk,
								 __private int i_endblk,
								 __private int nblocks_e,
								 __private int nblocks_c,
								 __global double * p_cell_in,
								 __global double * c_int,
								 __global int * iidx,
								 __global int * iblk,
								 __global int * localStart,
								 __global int * localEnd,
								 __global double * p_edge_out)
{
	__private int line_idx[2];      // was: __local
	__private double line_c_int[2]; // was: __local

	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[jb] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	for (int localIdx=0; localIdx<2; localIdx++)
	{
		const int idx = je + jb*nproma + localIdx*nproma*nblocks_e;
		const int idx2 = je + jb*nproma*2 + localIdx*nproma;
		line_idx[localIdx] = iidx[idx]-1 + (iblk[idx]-1)*nproma*nlev;
		line_c_int[localIdx] = c_int[idx2];
	}

#if 0
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	if (get_local_id(1) < 2 && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jb*nproma + get_local_id(1)*nproma*nblocks_e;
		const int idx2 = je + jb*nproma*2 + get_local_id(1)*nproma;
		
		line_idx[localIdx] = iidx[idx]-1 + (iblk[idx]-1)*nproma*nlev;
		line_c_int[localIdx] = c_int[idx2];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
		
	if (jk < elev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{	
		// index computations
		const int idx_out = je + jk*nproma + jb*nproma*nlev;
		const int idx_p_cell_in1 = line_idx[0] + jk*nproma;
		const int idx_p_cell_in2 = line_idx[1] + jk*nproma;
		
		p_edge_out[idx_out] = line_c_int[0] * p_cell_in[idx_p_cell_in1] +
							  line_c_int[1] * p_cell_in[idx_p_cell_in2];
	}
}

__kernel void edges2cells_scalar_tri(__private int nproma,
									 __private int slev,
									 __private int elev,
								     __private int nlev,
					                 __private int i_startblk,
									 __private int i_endblk,
					                 __private int nblocks_e,
					                 __private int nblocks_c,
					                 __global double * p_edge_in,
					                 __global double * c_int,
					                 __global int * iidx,
					                 __global int * iblk,
							         __global int * localStart,
					                 __global int * localEnd,
					                 __global double * p_cell_out)
{
	__private int line_idx[3];      // was: __local
	__private double line_c_int[3]; // was: __local

	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[jb] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	for (int localIdx=0; localIdx<3; localIdx++)
	{
		const int idx = jc + jb*nproma + localIdx*nproma*nblocks_c;
		const int idx2 = jc + jb*nproma*3 + localIdx*nproma;
		line_idx[localIdx] = iidx[idx]-1 + (iblk[idx]-1)*nproma*nlev;
		line_c_int[localIdx] = c_int[idx2];
	}

#if 0
// Original code, with expensive barrier
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_local_id(1) < 3 && jb < nblocks_e && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jb*nproma + get_local_id(1)*nproma*nblocks_c;
		const int idx2 = jc + jb*nproma*3 + get_local_id(1)*nproma;
		
		line_idx[localIdx] = iidx[idx]-1 + (iblk[idx]-1)*nproma*nlev;
		line_c_int[localIdx] = c_int[idx2];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
		
	if (jk < elev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{	
		// index computations
		const int idx_out = jc + jk*nproma + jb*nproma*nlev;
		const int idx_p_edge_in1 = line_idx[0] + jk*nproma;
		const int idx_p_edge_in2 = line_idx[1] + jk*nproma;
		const int idx_p_edge_in3 = line_idx[2] + jk*nproma;
		
		p_cell_out[idx_out] = line_c_int[0] * p_edge_in[idx_p_edge_in1] +
							  line_c_int[1] * p_edge_in[idx_p_edge_in2] +
							  line_c_int[2] * p_edge_in[idx_p_edge_in3];
	}
}

__kernel void edges2cells_scalar_hex()
{
	//Interpolate edge values to cells
}

__kernel void cell_avg(__private int nproma,
					   __private int slev,
					   __private int elev,
					   __private int nlev,
					   __private int i_startblk,
					   __private int i_endblk,
					   __private int nblocks_c,
					   __global double * psi_c,
					   __global double * avg_coeff,
					   __global int * iidx,
					   __global int * iblk,
					   __global int * localStart,
					   __global int * localEnd,
					   __global double * avg_psi_c)
{
	__private int line_idx[3];          // was: __local
	__private double line_avg_coeff[4]; // was: __local
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[jb] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	for (int localIdx=0; localIdx<3; localIdx++)
	{
		const int idx = jc + jb*nproma + localIdx*nproma*nblocks_c;
		line_idx[localIdx] = iidx[idx]-1 + (iblk[idx]-1)*nproma*nlev;
	}

	for (int localIdx=0; localIdx<4; localIdx++)
	{
		const int idx2 = jc + jb*nproma*4 + localIdx*nproma;
		line_avg_coeff[localIdx] = avg_coeff[idx2];
	}

#if 0
// Original code with expensive barrier
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_local_id(1) < 4 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx2 = jc + jb*nproma*4 + get_local_id(1)*nproma;
		line_avg_coeff[localIdx] = avg_coeff[idx2];
		
		if (get_local_id(1) < 3)
		{
			const int idx = jc + jb*nproma + get_local_id(1)*nproma*nblocks_c;
			line_idx[localIdx] = iidx[idx]-1 + (iblk[idx]-1)*nproma*nlev;
		}
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
		
	if (jk < elev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{	
		// index computations
		const int idx_out = jc + jk*nproma + jb*nproma*nlev;
		const int idx_psi_c1 = line_idx[0] + jk*nproma;
		const int idx_psi_c2 = line_idx[1] + jk*nproma;
		const int idx_psi_c3 = line_idx[2] + jk*nproma;
		
		avg_psi_c[idx_out] = psi_c[idx_out   ] * line_avg_coeff[0] +
							 psi_c[idx_psi_c1] * line_avg_coeff[1] +
							 psi_c[idx_psi_c2] * line_avg_coeff[2] +
							 psi_c[idx_psi_c3] * line_avg_coeff[3];
	}
}
