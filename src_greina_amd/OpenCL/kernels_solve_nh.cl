//==================================================
//
//	Kernels for the components of solve_nh
//		(mo_solve_nonhydro.f90)
//
//	See comments in kernels_operators.cl
//
//	the choice of iterators might need some rethinking
//
//	the kernels that use a 1d iterator assume that all the elements of the array are accessed!
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


__kernel void init(__private int nblks_c,
                   __private int nlev,
				   __private int nproma,
				   __global double * wgtfac_c,
				   __global double * wgtfac_c_out1,
				   __global double * wgtfac_c_out2)
{
	const int jb = get_global_id(0);
	const int jc = get_global_id(2);
	const int jk = get_global_id(1);
	
	if (jk < nlev && jb < nblks_c && jc < nproma)
	{
		const int idx_in  = jc + jk*nproma + jb*nproma*(nlev+1);
		const int idx_out = jc + jk*nproma + jb*nproma*nlev;
		
		wgtfac_c_out1[idx_out] = wgtfac_c[idx_in       ];
		wgtfac_c_out2[idx_out] = wgtfac_c[idx_in+nproma];
	}
}

__kernel void exner1(__private int i_startblk,
					 __private int i_endblk,
					 __private int nlev,
					 __private int nproma,
					 __global int * localStart, /*(nblks_c)*/
					 __global int * localEnd, /*(nblks_c)*/
					 __global double * exner_ref_mc, /*(nproma,nlev,nblks_c)?*/
					 __global double * exner_exfac, /*(nproma,nlev,nblks_c)*/
					 __global double * exner, /*(nproma,nlev,nblks_c)*/
					 __global double * exner_old, /*(nproma,nlev,nblks_c)*/
					 __global double * z_exner_ex_pr, /*(nproma,nlev,nblks_c)*/
					 __global double * z_exner_pr /*(nproma,nlev,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jk*nproma + jb*nproma*nlev;
		
		const double exner_idx = exner[idx];
		const double exner_ref_mc_idx = exner_ref_mc[idx];
		const double exner_exfac_idx = exner_exfac[idx];
		
		z_exner_ex_pr[idx] = -exner_ref_mc_idx + (1. + exner_exfac_idx)*exner_idx - exner_exfac_idx*exner_old[idx];
		
		z_exner_pr[idx] = exner_idx - exner_ref_mc_idx;
		
		exner_old[idx] = exner_idx;
	}
}

// nflat == 1
__kernel void exner2a(__private int i_startblk,
					  __private int i_endblk,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_c)*/
					  __global int * localEnd, /*(nblks_c)*/
					  __global double * wgtfacq1_c, /*(nproma,3,nblks_c)*/
					  __global double * wgtfacq_c, /*(nproma,3,nblks_c)*/
					  __global double * z_exner_ex_pr, /*(nproma,nlev,nblks_c)*/
					  __global double * z_exner_ic /*(nproma,nlevp1,nblkc_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(1);
	
	if (jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int wgtfacq_c_base = jc + jb*nproma*3;
		const int z_exner_ex_pr_base = jc + jb*nproma*nlev;
		
		z_exner_ic[jc + jb*nproma*nlevp1] = wgtfacq1_c[wgtfacq_c_base           ] * z_exner_ex_pr[z_exner_ex_pr_base           ] +
						                    wgtfacq1_c[wgtfacq_c_base +   nproma] * z_exner_ex_pr[z_exner_ex_pr_base +   nproma] +
						                    wgtfacq1_c[wgtfacq_c_base + 2*nproma] * z_exner_ex_pr[z_exner_ex_pr_base + 2*nproma];
		z_exner_ic[jc + (nlevp1-1)*nproma + jb*nproma*nlevp1] = wgtfacq_c[wgtfacq_c_base           ] * z_exner_ex_pr[z_exner_ex_pr_base + (nlev-1)*nproma] +
						                                        wgtfacq_c[wgtfacq_c_base +   nproma] * z_exner_ex_pr[z_exner_ex_pr_base + (nlev-2)*nproma] +
						                                        wgtfacq_c[wgtfacq_c_base + 2*nproma] * z_exner_ex_pr[z_exner_ex_pr_base + (nlev-3)*nproma];
	}
}

// nflat != 1
__kernel void exner2b(__private int i_startblk,
					  __private int i_endblk,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_c)*/
					  __global int * localEnd, /*(nblks_c)*/
					  __global double * wgtfacq_c, /*(nproma,3,nblks_c)*/
					  __global double * z_exner_ex_pr, /*(nproma,nlev,nblks_c)*/
					  __global double * z_exner_ic /*(nproma,nlevp1,nblkc_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(1);
	
	if (jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int wgtfacq_c_base = jc + jb*nproma*3;
		const int z_exner_ex_pr_base = jc + jb*nproma*nlev;
		
		z_exner_ic[jc + (nlevp1-1)*nproma + jb*nproma*nlevp1] = wgtfacq_c[wgtfacq_c_base           ]*z_exner_ex_pr[z_exner_ex_pr_base + (nlev-1)*nproma] +
						                                        wgtfacq_c[wgtfacq_c_base +   nproma]*z_exner_ex_pr[z_exner_ex_pr_base + (nlev-2)*nproma] +
						                                        wgtfacq_c[wgtfacq_c_base + 2*nproma]*z_exner_ex_pr[z_exner_ex_pr_base + (nlev-3)*nproma];
	}
}

__kernel void exner3a(__private int i_startblk,
					  __private int i_endblk,
					  __private int slev,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_c)*/
					  __global int * localEnd, /*(nblks_c)*/
					  __global double * wgtfac_c, /*(nproma,nlevp1,nblks_c)*/
					  __global double * z_exner_ex_pr, /*(nproma,nlev,nblks_c)*/
					  __global double * z_exner_ic /*(nproma,nlevp1,nblkc_c)*/
)
{
	__private double line_z_exner_ex_pr[NLEV_];  // was: __local

	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);

        for (int ilev=0; ilev < nlev; ilev++)
		line_z_exner_ex_pr[ilev] = z_exner_ex_pr[jc + ilev*nproma + jb*nproma*nlev];

#if 0
// Original code with slow barrier
	if (get_local_id(1) < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
		line_z_exner_ex_pr[get_local_id(1)] = z_exner_ex_pr[jc + get_local_id(1)*nproma + jb*nproma*nlev];
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jk*nproma + jb*nproma*nlevp1;
		
		const double wgtfac_c_tmp = wgtfac_c[idx];
		
		z_exner_ic[idx] = wgtfac_c_tmp*line_z_exner_ex_pr[jk] +
		                 (1.-wgtfac_c_tmp)*line_z_exner_ex_pr[jk-1];
	}
}

__kernel void exner3b(__private int i_startblk,
					  __private int i_endblk,
					  __private int slev,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_c)*/
					  __global int * localEnd, /*(nblks_c)*/
					  __global double * inv_ddqz_z_full, /*(nproma,nlev,nblks_c)*/
					  __global double * z_exner_ic, /*(nproma,nlevp1,nblkc_c)*/
					  __global double * z_dexner_dz_c /*(2,nproma,nlev,nblks_c)*/
)
{
	__private double line_z_exner_ic[NLEV_P1];  // was: __local
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
        for (int ilev=0; ilev < nlevp1; ilev++)
		line_z_exner_ic[ilev] = z_exner_ic[jc + ilev*nproma + jb*nproma*nlevp1];

#if 0
// Original code with slow barrier
	if (get_local_id(1) < nlevp1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
		line_z_exner_ic[get_local_id(1)] = z_exner_ic[jc + get_local_id(1)*nproma + jb*nproma*nlevp1];
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx_out = jc*2 + jk*2*nproma + jb*2*nproma*nlev;
	
		z_dexner_dz_c[idx_out] = (line_z_exner_ic[jk]-line_z_exner_ic[jk+1])*inv_ddqz_z_full[jc+jk*nproma+jb*nproma*nlev];
	}
}

__kernel void theta1(__private int i_startblk,
					 __private int i_endblk,
					 __private int nlev,
					 __private int nproma,
					 __global int * localStart, /*(nblks_c)*/
					 __global int * localEnd, /*(nblks_c)*/
					 __global double * theta_v_nvar, /*(nproma,nlev,nblks_c)*/
					 __global double * theta_v_nnow, /*(nproma,nlev,nblks_c)*/
					 __global double * theta_ref_mc, /*(nproma,nlev,nblks_c)*/
					 __global double * z_theta_v_pr_mc /*(nproma,nlev,nblks_c)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlev*nproma)
		z_theta_v_pr_mc[idx] = .5 * (theta_v_nnow[idx] + theta_v_nvar[idx]) - theta_ref_mc[idx];
}

__kernel void theta2(__private int i_startblk,
					 __private int i_endblk,
					 __private int nlev,
					 __private int nlevp1,
					 __private int nproma,
					 __global int * localStart, /*(nblks_c)*/
					 __global int * localEnd, /*(nblks_c)*/
					 __global double * rho_refcorr_ic, /*(nproma,nlevp1,nblks_c)*/
					 __global double * wgtfac_c, /*(nproma,nlevp1,nblks_c)*/
					 __global double * rho_nvar, /*(nproma,nlev,nblks_c)*/
					 __global double * rho_nnow, /*(nproma,nlev,nblks_c)*/
					 __global double * z_theta_v_pr_mc, /*(nproma,nlev,nblks_c)*/
					 __global double * theta_ref_ic,  /*(nproma,nlevp1,nblks_c)*/
					 __global double * vwind_expl_wgt, /*(nproma,nblks_c)*/
					 __global double * z_exner_pr, /*(nproma,nlev,nblks_c)*/
					 __global double * ddqz_z_half, /*(nproma,nlevp1,nblks_c)*/
					 __global double * d_exner_dz_ref_ic, /*(nproma,nlevp1,nblks_c)*/
					 __global double * rho_ic, /*(nproma,nlevp1,nblks_c)*/
					 __global double * z_theta_v_pr_ic, /*(nproma,nlevp1,nblks_c)*/
					 __global double * theta_v_h, /*(nproma,nlevp1,nblks_c)*/
					 __global double * z_th_ddz_exner_c /*(nproma,nlevp1,nblks_c)*/
)
{
// were: __local
	__private double line_rho_nnow[NLEV_];
	__private double line_rho_nvar[NLEV_];
	__private double line_z_theta_v_pr_mc[NLEV_];
	__private double line_z_exner_pr[NLEV_];

	const int jb = i_startblk + get_local_id(0) + get_group_id(0)*get_local_size(0);//get_global_id(0);
	const int jc = localStart[get_local_id(0) + get_group_id(0)*get_local_size(0)] + get_local_id(2) + get_group_id(2)*get_local_size(2);
	const int jk = 1 + get_local_id(1) + get_group_id(1)*get_local_size(1);
	
	__private double cached_vwind_expl_wgt;  // was: __local

	for (int ilev=0; ilev<nlev; ilev++)
	{
		const int idx1 = jc + ilev*nproma + jb*nproma*nlev;
		line_rho_nnow[ilev] = rho_nnow[idx1];
		line_rho_nvar[ilev] = rho_nvar[idx1];
		line_z_theta_v_pr_mc[ilev] = z_theta_v_pr_mc[idx1];
		line_z_exner_pr[ilev] = z_exner_pr[idx1];
	}
        cached_vwind_expl_wgt = vwind_expl_wgt[jc+jk*nproma];

#if 0
// Original code with slow barrier
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	if (get_local_id(1) + get_group_id(1)*get_local_size(1) < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + (get_local_id(1) + get_group_id(1)*get_local_size(1))*nproma + jb*nproma*nlev;
	
		line_rho_nnow[localIdx] = rho_nnow[idx1];
		line_rho_nvar[localIdx] = rho_nvar[idx1];
		line_z_theta_v_pr_mc[localIdx] = z_theta_v_pr_mc[idx1];
		line_z_exner_pr[localIdx] = z_exner_pr[idx1];
		
		if (get_local_id(1) == 0)
		 cached_vwind_expl_wgt = vwind_expl_wgt[jc+jk*nproma];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_local_id(0) + get_group_id(0)*get_local_size(0)])
	{
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
		
		const double wgtfac_c_idx2 = wgtfac_c[idx2];
							   
		const double z_theta_v_pr_ic_idx2 = wgtfac_c_idx2 * line_z_theta_v_pr_mc[jk] + (1. - wgtfac_c_idx2) * line_z_theta_v_pr_mc[jk-1];
		
		z_theta_v_pr_ic[idx2] = z_theta_v_pr_ic_idx2;
		
		const double theta_v_h_idx2 = theta_ref_ic[idx2] + z_theta_v_pr_ic_idx2;
		
		rho_ic[idx2] = rho_refcorr_ic[idx2] + .5 * (wgtfac_c_idx2 * (line_rho_nnow[jk] + line_rho_nvar[jk]) + 
					  (1. - wgtfac_c_idx2) * (line_rho_nnow[jk-1] + line_rho_nvar[jk-1]));
		
		theta_v_h[idx2] = theta_v_h_idx2;
		
		z_th_ddz_exner_c[idx2] = cached_vwind_expl_wgt * theta_v_h_idx2 *
								(line_z_exner_pr[jk-1] - line_z_exner_pr[jk]) / ddqz_z_half[idx2] +
								 z_theta_v_pr_ic_idx2 * d_exner_dz_ref_ic[idx2];
	}
}

__kernel void theta3a(__private int i_startblk,
					  __private int i_endblk,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_c)*/
					  __global int * localEnd, /*(nblks_c)*/
					  __global double * wgtfacq_c, /*(nproma,3,nblks_c)*/
					  __global double * z_theta_v_pr_mc, /*(nproma,nlev,nblks_c)*/
					  __global double * theta_ref_ic, /*(nproma,nlevp1,nblks_c)*/
					  __global double * z_theta_v_pr_ic, /*(nproma,nlevp1,nblks_c)*/
					  __global double * theta_v_h /*(nproma,nlevp1,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(1);
	
	if (jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		z_theta_v_pr_ic[jc+jb*nproma*nlevp1] = 0.;
		
		const int idx = jc+(nlevp1-1)*nproma+jb*nproma*nlevp1;
		
		z_theta_v_pr_ic[idx] =
			wgtfacq_c[jc         +jb*nproma*3] * z_theta_v_pr_mc[jc+(nlev-1)*nproma+jb*nproma*nlev] +
			wgtfacq_c[jc+  nproma+jb*nproma*3] * z_theta_v_pr_mc[jc+(nlev-2)*nproma+jb*nproma*nlev] +
			wgtfacq_c[jc+2*nproma+jb*nproma*3] * z_theta_v_pr_mc[jc+(nlev-3)*nproma+jb*nproma*nlev];
		
		theta_v_h[idx] = theta_ref_ic[idx] + z_theta_v_pr_ic[idx];
	}
}

__kernel void theta3b(__private int i_startblk,
					  __private int i_endblk,
					  __private int slev,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_c)*/
					  __global int * localEnd, /*(nblks_c)*/
					  __global double * theta_v, /*(nproma,nlev,nblks_c)*/
					  __global double * theta_v_h, /*(nproma,nlevp1,nblks_c)*/
					  __global double * z_theta_v_pr_ic, /*(nproma,nlevp1,nblks_c)*/
					  __global double * d_exner_dz_ref_mc, /*(nproma,nlev,nblks_c)*/
					  __global double * z_theta_v_pr_mc, /*(nproma,nlev,nblks_c)*/
					  __global double * d2_exner_dz2_ref_mc, /*(nproma,nlev,nblks_c)*/
					  __global double * z_dexner_dz_c /*(2,nproma,nlev,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + jk*nproma + jb*nproma*nlev;
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
		
		z_dexner_dz_c[1+2*idx1] = -.5 / theta_v[idx1] *
							   ((z_theta_v_pr_ic[idx2] - z_theta_v_pr_ic[idx2-nproma]) *
							    d_exner_dz_ref_mc[idx1] + z_theta_v_pr_mc[idx1] * d2_exner_dz2_ref_mc[idx1]);
	}
}

__kernel void grad1(__private int i_startblk,
					__private int i_endblk,
					__private int nblks_e,
					__private int elev,
					__private int nlev,
					__private int nproma,
					__global int * localStart, /*(nblks_e)*/
					__global int * localEnd, /*(nblks_e)*/
					__global int * icidx, /*(nproma,nblks_e,2)*/
					__global int * icblk, /*(nproma,nblks_e,2)*/
					__global double * inv_dual_edge_length, /*(nproma,nblks_e)*/
					__global double * z_exner_ex_pr, /*(nproma,nlev,nblks_c)*/
					__global double * z_gradh_exner /*(nproma,nlev,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	if (jk < elev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		// benefits from shared memory for the indexes -> currently the potential gain is too small w.r.t. to whole code (<1%)
		const int idx_out = je + jk*nproma + jb*nproma*nlev;
		const int idx_z_exner_ex_pr1 = icidx[je+jb*nproma               ]-1 + jk*nproma + (icblk[je+jb*nproma               ]-1)*nproma*nlev;
		const int idx_z_exner_ex_pr2 = icidx[je+jb*nproma+nproma*nblks_e]-1 + jk*nproma + (icblk[je+jb*nproma+nproma*nblks_e]-1)*nproma*nlev;
		
		z_gradh_exner[idx_out] = inv_dual_edge_length[je+jb*nproma] * (z_exner_ex_pr[idx_z_exner_ex_pr2] - z_exner_ex_pr[idx_z_exner_ex_pr1]);
	}
}

__kernel void grad2(__private int i_startblk,
					__private int i_endblk,
					__private int nblks_e,
					__private int slev,
					__private int elev,
					__private int nlev,
					__private int nproma,
					__global int * localStart, /*(nblks_e)*/
					__global int * localEnd, /*(nblks_e)*/
					__global int * icidx, /*(nproma,nblks_e,2)*/
					__global int * icblk, /*(nproma,nblks_e,2)*/
					__global double * inv_dual_edge_length, /*(nproma,nblks_e)*/
					__global double * z_exner_ex_pr, /*(nproma,nlev,nblks_c)*/
					__global double * ddxn_z_full, /*(nproma,nlev,nblks_e)*/
					__global double * c_lin_e,  /*(nproma,2,nblks_e)*/
					__global double * z_dexner_dz_c, /*(2,nproma,nlev,nblks_c)*/
					__global double * z_gradh_exner /*(nproma,nlev,nblks_e)*/
)
{
	__private int line_idx[2];        // was: __local
	__private double line_c_lin_e[2]; // was: __local
	__private double cached_inv_dual_edge_length; // was: __local
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);

	for (int ilev=0; ilev<2; ilev++)
	{
		const int idx = je + jb*nproma + ilev*nproma*nblks_e;
		line_idx[ilev] = icidx[idx]-1 + (icblk[idx]-1)*nproma*nlev;
		line_c_lin_e[ilev] = c_lin_e[je + jb*nproma*2 + get_global_id(1)];
	}
	cached_inv_dual_edge_length = inv_dual_edge_length[je+jb*nproma];


#if 0
// Original code with slow barrier
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	__local double cached_inv_dual_edge_length;
	
	if (get_global_id(1) < 2 && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jb*nproma + get_global_id(1)*nproma*nblks_e;
		
		line_idx[localIdx] = icidx[idx]-1 + (icblk[idx]-1)*nproma*nlev;
		line_c_lin_e[localIdx] = c_lin_e[je + jb*nproma*2 + get_global_id(1)];
	
		if (get_global_id(1) == 0)
			cached_inv_dual_edge_length = inv_dual_edge_length[je+jb*nproma];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < elev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jk*nproma + jb*nproma*nlev;
		const int idx_z_exner_ex_pr1 = line_idx[0] + jk*nproma;
		const int idx_z_exner_ex_pr2 = line_idx[1] + jk*nproma;
		
		z_gradh_exner[idx] = cached_inv_dual_edge_length * (z_exner_ex_pr[idx_z_exner_ex_pr2] - z_exner_ex_pr[idx_z_exner_ex_pr1]) -
						     ddxn_z_full[idx] * (line_c_lin_e[0]*z_dexner_dz_c[idx_z_exner_ex_pr1*2] +
												 line_c_lin_e[1]*z_dexner_dz_c[idx_z_exner_ex_pr2*2]);
	}
}

__kernel void grad3(__private int i_startblk,
					__private int i_endblk,
					__private int nblks_e,
					__private int slev,
					__private int elev,
					__private int nlev,
					__private int nproma,
					__global int * localStart, /*(nblks_e)*/
					__global int * localEnd, /*(nblks_e)*/
					__global int * icidx, /*(nproma,nblks_e,2)*/
					__global int * icblk, /*(nproma,nblks_e,2)*/
					__global int * ikidx, /*(2,nproma,nlev,nblks_e)*/
					__global double * inv_dual_edge_length, /*(nproma,nblks_e)*/
					__global double * z_exner_ex_pr, /*(nproma,nlev,nblks_c)*/
					__global double * zdiff_gradp, /*(2,nproma,nlev,nblks_e)*/
					__global double * z_dexner_dz_c, /*(2,nproma,nlev,nblks_c)*/
					__global double * z_gradh_exner /*(nproma,nlev,nblks_e)*/)
{
// unused
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	if (jk < elev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jk*nproma + jb*nproma*nlev;
		const int idx_z_exner_ex_pr1 = icidx[je+jb*nproma               ]-1 + (ikidx[  idx*2]-1)*nproma + (icblk[je+jb*nproma               ]-1)*nproma*nlev;
		const int idx_z_exner_ex_pr2 = icidx[je+jb*nproma+nproma*nblks_e]-1 + (ikidx[1+idx*2]-1)*nproma + (icblk[je+jb*nproma+nproma*nblks_e]-1)*nproma*nlev;
		
		z_gradh_exner[idx] = inv_dual_edge_length[je+jb*nproma] *
		                    (z_exner_ex_pr[idx_z_exner_ex_pr2]     + zdiff_gradp[1+idx*2] *
							(z_dexner_dz_c[  idx_z_exner_ex_pr2*2] + zdiff_gradp[1+idx*2] *
							 z_dexner_dz_c[1+idx_z_exner_ex_pr2*2] -
							(z_exner_ex_pr[idx_z_exner_ex_pr1]     + zdiff_gradp[  idx*2] *
							(z_dexner_dz_c[  idx_z_exner_ex_pr1*2] + zdiff_gradp[  idx*2] *
							 z_dexner_dz_c[1+idx_z_exner_ex_pr1*2]))));
	}
}

__kernel void grad4(__private int i_startblk,
					__private int i_endblk,
					__private int nblks_e,
					__private int nlev,
					__private int nlevp1,
					__private int nproma,
					__private double grav_o_cpd,
					__global int * localStart, /*(nblks_e)*/
					__global int * localEnd, /*(nblks_e)*/
					__global int * icidx, /*(nproma,nblks_e,2)*/
					__global int * icblk, /*(nproma,nblks_e,2)*/
					__global int * ikidx, /*(2,nproma,nlev,nblks_e)*/
					__global double * theta_v, /*(nproma,nlev,nblks_c)*/
					__global double * zdiff_gradp, /*(2,nproma,nlev,nblks_e)*/
					__global double * theta_v_h, /*(nproma,nlevp1,nblks_c)*/
					__global double * inv_ddqz_z_full, /*(nproma,nlev,nblks_c)*/
					__global double * inv_dual_edge_length, /*(nproma,nblks_e)*/
					__global double * z_hydro_corr /*(nproma,nblks_e)*/)
{
// unused
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(1);
	
	if (jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + (nlev-1)*nproma + jb*nproma*nlev;
		const int idx_t1 = icidx[je+jb*nproma               ]-1 + (ikidx[  idx*2]-1)*nproma + (icblk[je+jb*nproma               ]-1)*nproma*nlev;
		const int idx_t2 = icidx[je+jb*nproma+nproma*nblks_e]-1 + (ikidx[1+idx*2]-1)*nproma + (icblk[je+jb*nproma+nproma*nblks_e]-1)*nproma*nlev;
		
		const double z_theta1 = theta_v[idx_t1] + zdiff_gradp[  idx*2] * (theta_v_h[idx_t1] - theta_v_h[idx_t1+nproma]) * inv_ddqz_z_full[idx_t1];
		const double z_theta2 = theta_v[idx_t2] + zdiff_gradp[1+idx*2] * (theta_v_h[idx_t2] - theta_v_h[idx_t2+nproma]) * inv_ddqz_z_full[idx_t2];
		
		z_hydro_corr[je+jb*nproma] = grav_o_cpd * inv_dual_edge_length[je+jb*nproma] * (z_theta2-z_theta1) * 4 / ((z_theta1+z_theta2)*(z_theta1+z_theta2));
	}
}


__kernel void vn1a(__private int i_startblk,
				   __private int i_endblk,
				   __private int nblks_e,
				   __private int nlev,
				   __private int nproma,
				   __private int ntl1, // remember to do -1
				   __private int ntl2, // remember to do -1
				   __private double dtime,
				   __private double cpd,
				   __global int * localStart, /*(nblks_e)*/
				   __global int * localEnd, /*(nblks_e)*/
				   __global double * vn_nnow, /*(nproma,nlev,nblks_e)*/
				   __global double * ddt_vn_adv, /*(nproma,nlev,nblks_e,3)*/
				   __global double * ddt_vn_phy, /*(nproma,nlev,nblks_e)*/
				   __global double * z_theta_v_e, /*(nproma,nlev,nblks_e)*/
				   __global double * z_gradh_exner, /*(nproma,nlev,nblks_e)*/
				   __global double * vn_nnew /*(nproma,nlev,nblks_e)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlev*nproma)
	{
		// benefits from shared memory for ddt_vn_adv
		vn_nnew[idx] = vn_nnow[idx] + dtime * (.5 *
					    (ddt_vn_adv[idx+ntl1*nproma*nlev*nblks_e] + ddt_vn_adv[idx+ntl2*nproma*nlev*nblks_e]) +
					    ddt_vn_phy[idx] - cpd * z_theta_v_e[idx] * z_gradh_exner[idx]);
	}
}

__kernel void vn1b(__private int i_startblk,
				   __private int i_endblk,
				   __private int nblks_e,
				   __private int nlev,
				   __private int nproma,
				   __private int ntl1,
				   __private double dtime,
				   __private double cpd,
				   __global int * localStart, /*(nblks_e)*/
				   __global int * localEnd, /*(nblks_e)*/
				   __global double * vn_nnow, /*(nproma,nlev,nblks_e)*/
				   __global double * ddt_vn_adv, /*(nproma,nlev,nblks_e,3)*/
				   __global double * ddt_vn_phy, /*(nproma,nlev,nblks_e)*/
				   __global double * z_theta_v_e, /*(nproma,nlev,nblks_e)*/
				   __global double * z_gradh_exner, /*(nproma,nlev,nblks_e)*/
				   __global double * vn_nnew /*(nproma,nlev,nblks_e)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlev*nproma)
	{
		vn_nnew[idx] = vn_nnow[idx] + dtime *
					    (ddt_vn_adv[idx+ntl1*nproma*nlev*nblks_e] +
					    ddt_vn_phy[idx] - cpd * z_theta_v_e[idx] * z_gradh_exner[idx]);
	}
}

__kernel void vn_avg(__private int i_startblk,
					 __private int i_endblk,
				     __private int nblks_e,
					 __private int nlev,
					 __private int nproma,
					 __global int * localStart, /*(nblks_e)*/
					 __global int * localEnd, /*(nblks_e)*/
					 __global int * iqidx, /*(nproma,nblks_e,4)*/
					 __global int * iqblk, /*(nproma,nblks_e,4)*/
					 __global double * vn_nnew, /*(nproma,nlev,nblks_e)*/
					 __global double * e_flx_avg, /*(nproma,5,nblks_e)*/
					 __global double * z_vn_avg /*(nproma,nlev,nblks_e)*/
)
{
	 __private double line_e_flx_avg[5];  // was: __local
	 __private int line_idx[4];           // was: __local

	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);

	for (int ilev=0; ilev<5; ilev++)
	{
		const int idx = je + ilev*nproma + jb*nproma*5;
		line_e_flx_avg[ilev] = e_flx_avg[idx];
	}

	for (int ilev=0; ilev<4; ilev++)
	{
		const int idx = je+jb*nproma+ilev*nproma*nblks_e;
		line_idx[ilev] = iqidx[idx]-1 + (iqblk[idx]-1)*nproma*nlev;
	}

#if 0	
// Original code with slow barrier
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_local_id(1) < 5 && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + get_local_id(1)*nproma + jb*nproma*5;
		
		line_e_flx_avg[localIdx] = e_flx_avg[idx];
		
		if (get_local_id(1) < 4)
		{
			const int idx = je+jb*nproma+get_local_id(1)*nproma*nblks_e;
			line_idx[localIdx] = iqidx[idx]-1 + (iqblk[idx]-1)*nproma*nlev;
		}
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jk*nproma + jb*nproma*nlev;
		const int idx_vn1 = line_idx[0] + jk*nproma;
		const int idx_vn2 = line_idx[1] + jk*nproma;
		const int idx_vn3 = line_idx[2] + jk*nproma;
		const int idx_vn4 = line_idx[3] + jk*nproma;
		
		z_vn_avg[idx] = vn_nnew[idx    ] * line_e_flx_avg[0] +
		                vn_nnew[idx_vn1] * line_e_flx_avg[1] +
		                vn_nnew[idx_vn2] * line_e_flx_avg[2] +
		                vn_nnew[idx_vn3] * line_e_flx_avg[3] +
		                vn_nnew[idx_vn4] * line_e_flx_avg[4];
	}
}

__kernel void vn_half1(__private int i_startblk,
					   __private int i_endblk,
					   __private int slev,
					   __private int nlev,
					   __private int nlevp1,
					   __private int nproma,
					   __global int * localStart, /*(nblks_e)*/
					   __global int * localEnd, /*(nblks_e)*/
					   __global double * wgtfac_e, /*(nproma,nlevp1,nblks_e)*/
					   __global double * vn_nnew, /*(nproma,nlev,nblks_e)*/
					   __global double * vt, /*(nproma,nlev,nblks_e)*/
					   __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					   __global double * z_vt_ie /*(nproma,nlevp1,nblks_e)*/
)
{
	__private double line_vn_nnew[NLEV_];  // was: __local
	__private double line_vt[NLEV_];       // was: __local

	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);

	for (int ilev=0; ilev<nlev; ilev++)
	{
		const int idx = je + ilev*nproma + jb*nproma*nlev;
		line_vn_nnew[ilev] = vn_nnew[idx];
		line_vt[ilev] = vt[idx];
	}

#if 0	
// Original code with slow barrier
	if (get_local_id(1) < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + get_local_id(1)*nproma + jb*nproma*nlev;
		
		line_vn_nnew[get_local_id(1)] = vn_nnew[idx];
		line_vt[get_local_id(1)] = vt[idx];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx2 = je + jk*nproma + jb*nproma*nlevp1;
		
		const double wgtfac_e_idx2 = wgtfac_e[idx2];
		
		vn_half[idx2] = wgtfac_e_idx2 * line_vn_nnew[jk] + (1.-wgtfac_e_idx2) * line_vn_nnew[jk-1];
		z_vt_ie[idx2] = wgtfac_e_idx2 *      line_vt[jk] + (1.-wgtfac_e_idx2) *      line_vt[jk-1];
	}
}

__kernel void vn_half2(__private int i_startblk,
					   __private int i_endblk,
					   __private int nlev,
					   __private int nlevp1,
					   __private int nproma,
					   __global int * localStart, /*(nblks_e)*/
					   __global int * localEnd, /*(nblks_e)*/
					   __global double * wgtfacq1_e, /*(nproma,3,nblks_e)*/
					   __global double * vn_nnew, /*(nproma,nlev,nblks_e)*/
					   __global double * vn_half /*(nproma,nlevp1,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(1);
	
	if (jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx1 = je + jb*nproma*nlev;
		const int idx2 = je + jb*nproma*nlevp1;
		
		vn_half[idx2] = wgtfacq1_e[je         +jb*nproma*3] * vn_nnew[idx1         ] +
		                wgtfacq1_e[je+  nproma+jb*nproma*3] * vn_nnew[idx1+  nproma] +
		                wgtfacq1_e[je+2*nproma+jb*nproma*3] * vn_nnew[idx1+2*nproma];
	}
}

__kernel void vn_half3(__private int i_startblk,
					   __private int i_endblk,
					   __private int elev,
					   __private int nlev,
					   __private int nlevp1,
					   __private int nproma,
					   __global int * localStart, /*(nblks_e)*/
					   __global int * localEnd, /*(nblks_e)*/
					   __global double * wgtfac_e, /*(nproma,nlevp1,nblks_e)*/
					   __global double * vn_nnew, /*(nproma,nlev,nblks_e)*/
					   __global double * vn_half /*(nproma,nlevp1,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = 1 + get_global_id(1);
	
	if (jk < elev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		// could benefit from shared memory for vn_nnew
		const int idx1 = je + jk*nproma + jb*nproma*nlev;
		const int idx2 = je + jk*nproma + jb*nproma*nlevp1;
		
		const double wgtfac_e_tmp = wgtfac_e[idx2];
		
		vn_half[idx2] = wgtfac_e_tmp * vn_nnew[idx1] +
		               (1.-wgtfac_e_tmp) * vn_nnew[idx1-nproma];
	}
}

__kernel void vn_half4(__private int i_startblk,
					   __private int i_endblk,
					   __private int nlev,
					   __private int nlevp1,
					   __private int nproma,
					   __global int * localStart, /*(nblks_e)*/
					   __global int * localEnd, /*(nblks_e)*/
					   __global double * wgtfacq_e, /*(nproma,3,nblks_e)*/
					   __global double * vn_nnew, /*(nproma,nlev,nblks_e)*/
					   __global double * vt, /*(nproma,nlev,nblks_e)*/
					   __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					   __global double * z_vt_ie /*(nproma,nlevp1,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(1);
	const int jk = nlevp1-1;
	
	if (jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx2 = je + jk*nproma + jb*nproma*nlevp1;
		
		const double wgt0 = wgtfacq_e[je         +jb*nproma*3];
		const double wgt1 = wgtfacq_e[je+  nproma+jb*nproma*3];
		const double wgt2 = wgtfacq_e[je+2*nproma+jb*nproma*3];
		
		const int idx_vm1 = je+(nlev-1)*nproma+jb*nproma*nlev;
		const int idx_vm2 = je+(nlev-2)*nproma+jb*nproma*nlev;
		const int idx_vm3 = je+(nlev-3)*nproma+jb*nproma*nlev;
		
		vn_half[idx2] = wgt0 * vn_nnew[idx_vm1] +
		                wgt1 * vn_nnew[idx_vm2] +
		                wgt2 * vn_nnew[idx_vm3];
		
		z_vt_ie[idx2] = wgt0 * vt[idx_vm1] +
		                wgt1 * vt[idx_vm2] +
		                wgt2 * vt[idx_vm3];
	}
}

__kernel void vn_half5(__private int i_startblk,
					   __private int i_endblk,
					   __private int slev,
					   __private int nlevp1,
					   __private int nproma,
					   __global int * localStart, /*(nblks_e)*/
					   __global int * localEnd, /*(nblks_e)*/
					   __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					   __global double * ddxn_z_half, /*(nproma,nlevp1,nblks_e)*/
					   __global double * z_vt_ie, /*(nproma,nlevp1,nblks_e)*/
					   __global double * ddxt_z_half, /*(nproma,nlevp1,nblks_e)*/
					   __global double * z_concorr_e /*(nproma,nlevp1,nblks_e)*/)
{
	const int idx = i_startblk*nlevp1*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlevp1*nproma)
		z_concorr_e[idx] = vn_half[idx] * ddxn_z_half[idx] + z_vt_ie[idx] * ddxt_z_half[idx];
		
	/*
	// this code is extremely slow compared to the one above, even though it does a lot less computation
	// note that the additional computation does not interfere with the results as
	//  levels [0-slev[ are filled but not used!
	
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	if (jk < nlevp1 && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jk*nproma + jb*nproma*nlevp1;
		
		z_concorr_e[idx] = vn_half[idx] * ddxn_z_half[idx] + z_vt_ie[idx] * ddxt_z_half[idx];
	}
	*/
}

__kernel void fluxes(__private int i_startblk,
					 __private int i_endblk,
					 __private int nlev,
					 __private int nproma,
					 __global int * localStart, /*(nblks_e)*/
					 __global int * localEnd, /*(nblks_e)*/
					 __global double * z_rho_e, /*(nproma,nlev,nblks_e)*/
					 __global double * ptr_vn, /*(nproma,nlev,nblks_e)*/
					 __global double * ddqz_z_full_e, /*(nproma,nlev,nblks_e)*/
					 __global double * z_theta_v_e, /*(nproma,nlev,nblks_e)*/
					 __global double * mass_fl_e, /*(nproma,nlev,nblks_e)*/
					 __global double * z_theta_v_fl_e /*(nproma,nlev,nblks_e)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlev*nproma)
	{
		mass_fl_e[idx] = z_rho_e[idx] * ptr_vn[idx] * ddqz_z_full_e[idx];
		z_theta_v_fl_e[idx] = mass_fl_e[idx] * z_theta_v_e[idx];
	}
}

__kernel void endphase1a(__private int i_startblk,
						 __private int i_endblk,
						 __private int nblks_c,
						 __private int nlev,
						 __private int nlevp1,
						 __private int nproma,
						 __private double dtime,
						 __private double cpd,
						 __global int * localStart, /*(nblks_c)*/
						 __global int * localEnd, /*(nblks_c)*/
						 __global double * w_nnow, /*(nproma,nlevp1,nblks_c)*/
						 __global double * ddt_w_adv1, /*(nproma,nlevp1,nblks_c)*/
						 __global double * ddt_w_adv2, /*(nproma,nlevp1,nblks_c)*/
						 __global double * z_th_ddz_exner_c, /*(nproma,nlev,nblks_c)*/
						 __global double * z_w_expl /*(nproma,nlevp1,nblks_c)*/)
{
	const int idx = i_startblk*nlevp1*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlevp1*nproma)
	{
		z_w_expl[idx] = w_nnow[idx] + dtime * (.5 * (ddt_w_adv1[idx] + ddt_w_adv2[idx]) - cpd * z_th_ddz_exner_c[idx]);
	}
}

__kernel void endphase1b(__private int i_startblk,
						 __private int i_endblk,
						 __private int nblks_c,
						 __private int nlev,
						 __private int nlevp1,
						 __private int nproma,
						 __private double dtime,
						 __private double cpd,
						 __global int * localStart, /*(nblks_c)*/
						 __global int * localEnd, /*(nblks_c)*/
						 __global double * w_nnow, /*(nproma,nlevp1,nblks_c)*/
						 __global double * ddt_w_adv, /*(nproma,nlevp1,nblks_c)*/
						 __global double * z_th_ddz_exner_c, /*(nproma,nlev,nblks_c)*/
						 __global double * z_w_expl /*(nproma,nlevp1,nblks_c)*/)
{
	const int idx = i_startblk*nlevp1*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlevp1*nproma)
	{
		z_w_expl[idx] = w_nnow[idx] + dtime * (ddt_w_adv[idx] - cpd * z_th_ddz_exner_c[idx]);
	}
}

__kernel void endphase1c(__private int i_startblk,
						 __private int i_endblk,
						 __private int nblks_c,
						 __private int nlev,
						 __private int nlevp1,
						 __private int nproma,
						 __global int * localStart, /*(nblks_c)*/
						 __global int * localEnd, /*(nblks_c)*/
						 __global double * w_nnow, /*(nproma,nlevp1,nblks_c)*/
						 __global double * rho_ic, /*(nproma,nlevp1,nblks_c)*/
						 __global double * w_concorr_c, /*(nproma,nlevp1,nblks_c)*/
						 __global double * vwind_expl_wgt, /*(nproma,nblks_c)*/
						 __global double * z_contr_w_fl_l /*(nproma,nlevp1,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = 1 + get_global_id(1);
	
	__private double cached_vwind_expl_wgt;
	cached_vwind_expl_wgt = vwind_expl_wgt[jc+jb*nproma];

#if 0
	if (jk == 1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
		cached_vwind_expl_wgt = vwind_expl_wgt[jc+jb*nproma];
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jk*nproma + jb*nproma*nlevp1;
		
		z_contr_w_fl_l[idx] = rho_ic[idx] * (-w_concorr_c[idx] + cached_vwind_expl_wgt * w_nnow[idx]);
	}
}

__kernel void endphase2a(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nproma,
						__private double dtime,
						__private double rd,
						__private double cvd,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * exner_nnow, /*(nproma,nlev,nblks_c)*/
						__global double * rhotheta_v_nnow, /*(nproma,nlev,nblks_c)*/
						__global double * ddqz_z_full, /*(nproma,nlev,nblks_c)*/
						__global double * z_beta /*(nproma,nlev,nblks_c)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlev*nproma)
	{	
		z_beta[idx] = dtime * rd * exner_nnow[idx] / (cvd * rhotheta_v_nnow[idx] * ddqz_z_full[idx]);
	}
}

__kernel void endphase2b(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * vwind_impl_wgt, /*(nproma,nblks_c)*/
						__global double * theta_v_h, /*(nproma,nlevp1,nblks_c)*/
						__global double * rho_ic, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_alpha /*(nproma,nlevp1,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jk*nproma + jb*nproma*nlevp1;
		
		z_alpha[idx] = vwind_impl_wgt[jc+jb*nproma] * theta_v_h[idx] * rho_ic[idx];
		
		if (jk==0)
			z_alpha[idx+(nlevp1-1)*nproma] = 0.;
	}
}

__kernel void endphase3(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__private double dtime,
						__private double cpd,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * vwind_impl_wgt, /*(nproma,nblks_c)*/
						__global double * theta_v_h, /*(nproma,nlevp1,nblks_c)*/
						__global double * ddqz_z_half, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_alpha, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_beta, /*(nproma,nlev,nblks_c)*/
						__global double * z_gamma, /*(nproma,nlev,nblks_c)*/
						__global double * z_a, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_b, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_c, /*(nproma,nlev,nblks_c)*/
						__global double * z_q /*(nproma,nlev,nblks_c)*/
)
{
	__private double line_z_alpha[NLEV_P1];
	__private double line_z_beta[NLEV_];
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = 1 + get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);

	for (int ilev=0; ilev<nlevp1; ilev++)
	{
		const int idx2 = jc + ilev*nproma + jb*nproma*nlevp1;
		line_z_alpha[ilev] = z_alpha[idx2];
	}

	for (int ilev=0; ilev<nlev; ilev++)
	{
		const int idx1 = jc + ilev*nproma + jb*nproma*nlev;
		line_z_beta[ilev]  = z_beta[idx1];
	}


#if 0
// Original code with slow barrier
	if (get_global_id(1) < nlevp1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx2 = jc + get_global_id(1)*nproma + jb*nproma*nlevp1;
		line_z_alpha[localIdx] = z_alpha[idx2];
		
		if (get_global_id(1) < nlev)
		{
			const int idx1 = jc + get_global_id(1)*nproma + jb*nproma*nlev;
			line_z_beta[localIdx]  = z_beta[idx1];
		}
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
		
		const double z_gamma_res = dtime * cpd * vwind_impl_wgt[jc+jb*nproma] * theta_v_h[idx2] / ddqz_z_half[idx2];
		const double z_betam1 = line_z_beta[localIdx];
		const double z_beta0  = line_z_beta[localIdx+nproma];
		
		z_gamma[idx2] = z_gamma_res;
		z_a[idx2] = -z_gamma_res * z_betam1 * line_z_alpha[localIdx];
		const double z_c_tmp = z_gamma_res * z_beta0  * line_z_alpha[localIdx+2];
		z_c[idx2] = -z_c_tmp;
		const double z_b_tmp = 1. + z_gamma_res * line_z_alpha[localIdx+1] * (z_betam1 + z_beta0);
		z_b[idx2] = z_b_tmp;

		if (jk==1)
			z_q[idx2] = z_c_tmp/z_b_tmp;
	}
}

__kernel void endphase4(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * z_a, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_b, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_c, /*(nproma,nlev,nblks_c)*/
						__global double * z_g, /*(nproma,nlev,nblks_c)*/
						__global double * z_q /*(nproma,nlev,nblks_c)*/)
{
	// adding nlev more threads in the 3rd dimension to put z_a, z_b, z_c in local memory slows down the kernel
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(1);
	
	__local double z_g_tmp, z_q_tmp;
	
	if (jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jb*nproma*nlevp1;
	
		// there are dependencies here between z_g and z_q!!
		z_q_tmp = z_q[idx+nproma];
		for (int jk=2; jk<nlev; jk++)
		{
			const int idx_loop = idx+jk*nproma;
			z_g_tmp = 1. / (z_b[idx_loop] + z_a[idx_loop] * z_q_tmp);
			z_g[idx_loop] = z_g_tmp;
			z_q_tmp = - z_c[idx_loop] * z_g_tmp;
			z_q[idx_loop] = z_q_tmp;
		}
	}
}

__kernel void endphase5a(__private int i_startblk,
						__private int i_endblk,
						__private int nlevp1,
						__private int nproma,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * w_concorr_c, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_contr_w_fl_l, /*(nproma,nlevp1,nblks_c)*/
						__global double * w_nnew /*(nproma,nlevp1,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1)*(nlevp1-1);
	
	if (get_global_id(1) < 2 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jk*nproma + jb*nproma*nlevp1;
		
		w_nnew[idx] = get_global_id(1)*w_concorr_c[idx];
		z_contr_w_fl_l[idx] = 0.;
	}
}

__kernel void endphase5b(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__private double dtime,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * z_contr_w_fl_l, /*(nproma,nlevp1,nblks_c)*/
						__global double * theta_v_h, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_contr_w_fl_l_diff,
						__global double * z_th_contr_diff
)
{
	__private double line_theta_v_h[NLEV_P1];
	__private double line_z_contr_w_fl_l[NLEV_P1];
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	for (int ilev=0; ilev<nlevp1; ilev++)
	{
		const int idx2 = jc + ilev*nproma + jb*nproma*nlevp1;
		line_theta_v_h[ilev]  = theta_v_h[idx2];
		line_z_contr_w_fl_l[ilev] = z_contr_w_fl_l[idx2];
	}

#if 0
// Original code with slow barrier
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (jk < nlevp1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
	
		line_theta_v_h[localIdx]  = theta_v_h[idx2];
		line_z_contr_w_fl_l[localIdx] = z_contr_w_fl_l[idx2];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif	

	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + jk*nproma + jb*nproma*nlev;
		
		const double z_contr_w_fl_l_idx2  = line_z_contr_w_fl_l[jk];
		const double z_contr_w_fl_l_idx2b = line_z_contr_w_fl_l[jk+1];
		
		z_contr_w_fl_l_diff[idx1] = z_contr_w_fl_l_idx2 - z_contr_w_fl_l_idx2b;
		z_th_contr_diff[idx1] = line_theta_v_h[jk] * z_contr_w_fl_l_idx2 - line_theta_v_h[jk+1] * z_contr_w_fl_l_idx2b;
	}
}

__kernel void endphase5c(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nproma,
						__private double dtime,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * rho_nnow, /*(nproma,nlev,nblks_c)*/
						__global double * inv_ddqz_z_full, /*(nproma,nlev,nblks_c)*/
						__global double * z_mass_fl_div, /*(nproma,nlev,nblks_c)*/
						__global double * z_contr_w_fl_l_diff, /*(nproma,nlev,nblks_c)*/
						__global double * z_exner_pr, /*(nproma,nlev,nblks_c)*/
						__global double * z_beta, /*(nproma,nlev,nblks_c)*/
						__global double * z_theta_v_fl_div, /*(nproma,nlev,nblks_c)*/
						__global double * z_th_contr_diff, /*(nproma,nlev,nblks_c)*/
						__global double * ddt_exner, /*(nproma,nlev,nblks_c)*/
						__global double * ddt_exner_phy, /*(nproma,nlev,nblks_c)*/
						__global double * z_rho_expl, /*(nproma,nlev,nblks_c)*/
						__global double * z_exner_expl /*(nproma,nlev,nblks_c)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);
	
	if (idx < i_endblk*nlev*nproma)
	{
		z_rho_expl[idx] = rho_nnow[idx] - dtime * inv_ddqz_z_full[idx] *
		                  (z_mass_fl_div[idx] + z_contr_w_fl_l_diff[idx]);
		z_exner_expl[idx] = z_exner_pr[idx] - z_beta[idx] *
							(z_theta_v_fl_div[idx] + z_th_contr_diff[idx]) +
						    dtime * (ddt_exner[idx] + ddt_exner_phy[idx]);
	}
}

__kernel void endphase6(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * z_w_expl, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_gamma, /*(nproma,nlev,nblks_c)*/
						__global double * z_exner_expl, /*(nproma,nlev,nblks_c)*/
						__global double * z_b, /*(nproma,nlevp1,nblks_c)*/
						__global double * w_nnew /*(nproma,nlevp1,nblks_c)*/
)
{
	__private double line_z_exner_expl[NLEV_];
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = 1 + get_global_id(1);

	for (int ilev=0; ilev<nlev; ilev++)
	{
		const int idx1 = jc + ilev*nproma + jb*nproma*nlev;
		line_z_exner_expl[ilev] = z_exner_expl[idx1];
	}

#if 0
// Original code with slow barrier
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_global_id(1) < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + get_global_id(1)*nproma + jb*nproma*nlev;
	
		line_z_exner_expl[localIdx] = z_exner_expl[idx1];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
		
		w_nnew[idx2] = (z_w_expl[idx2] - z_gamma[idx2] * (line_z_exner_expl[jk-1] - line_z_exner_expl[jk])) / (jk==1 ? z_b[idx2] : 1.);
	}
}

__kernel void endphase7(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * z_a, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_g, /*(nproma,nlev,nblks_c)*/
						__global double * z_q, /*(nproma,nlev,nblks_c)*/
						__global double * w_nnew /*(nproma,nlevp1,nblks_c)*/
)
{
	__private double line_w_nnew[NLEV_];  // was: __local
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0);

	const int idx = jc + jb*nproma*nlevp1;
	for (int jk=0; jk<nlevp1; jk++)
		line_w_nnew[jk] = w_nnew[idx+jk*nproma];

#if 0
	// TODO:  This code makes no sense; out of bounds problems?? (WS)
	if (jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jb*nproma*nlevp1;
		for (int jk=0; jk<nlevp1; jk++)
			line_w_nnew[localIdx+jk*nproma] = w_nnew[idx+jk*nproma];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jb*nproma*nlevp1;
		
		// there are dependencies here between the various levels!
		for (int jk=2; jk<nlev; jk++)
			line_w_nnew[localIdx+jk*nproma] = (line_w_nnew[localIdx+jk] - z_a[idx+jk*nproma] * line_w_nnew[localIdx+(jk-1)]) * z_g[idx+jk*nproma];
		
		for (int jk=nlev-2; jk>=1; jk--)
			line_w_nnew[localIdx+jk] += line_w_nnew[localIdx+(jk+1)] * z_q[idx+jk*nproma];
		
		for (int jk=0; jk<nlevp1; jk++)
			w_nnew[idx+jk] = line_w_nnew[localIdx+jk];
	}
}

__kernel void endphase8(__private int i_startblk,
						__private int i_endblk,
						__private int elev,
						__private int nlevp1,
						__private int nproma,
						__private double dtime,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * rayleigh_w, /*(nproma,nlevp1,nblks_c)*/
						__global double * w_nnew /*(nproma,nlevp1,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = 1 + get_global_id(1);
	
	if (jk < elev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jk*nproma + jb*nproma*nlevp1;
		
		w_nnew[idx] = (w_nnew[idx] - w_nnew[jc+jb*nproma*nlevp1]) / (1. + dtime * rayleigh_w[idx]);
	}
}

__kernel void endphase9a(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__private double dtime,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * z_rho_expl, /*(nproma,nlev,nblks_c)*/
						__global double * vwind_impl_wgt, /*(nproma,nblks_c)*/
						__global double * inv_ddqz_z_full, /*(nproma,nlev,nblks_c)*/
						__global double * rho_ic, /*(nproma,nlevp1,nblks_c)*/
						__global double * w_nnew, /*(nproma,nlevp1,nblks_c)*/
						__global double * rho_nnew /*(nproma,nlev,nblks_c)*/
)
{
	__private double line_rho_ic[NLEV_P1];
	__private double line_w_nnew[NLEV_P1];
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
        for ( int ilev=0; ilev < nlevp1; ilev++ )
	{
		const int idx2 = jc + ilev*nproma + jb*nproma*nlevp1;
		line_rho_ic[ilev] = rho_ic[idx2];
		line_w_nnew[ilev] = w_nnew[idx2];
	}

#if 0
	if (jk < nlevp1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
	
		line_rho_ic[localIdx] = rho_ic[idx2];
		line_w_nnew[localIdx] = w_nnew[idx2];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
#endif
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + jk*nproma + jb*nproma*nlev;
		
		rho_nnew[idx1] = z_rho_expl[idx1] - vwind_impl_wgt[jc+jb*nproma] * dtime * inv_ddqz_z_full[idx1] *
						(line_rho_ic[localIdx] * line_w_nnew[localIdx] - line_rho_ic[localIdx+nproma] * line_w_nnew[localIdx+nproma]);
	}
}

__kernel void endphase9b(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nlevp1,
						__private int nproma,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * w_nnew, /*(nproma,nlevp1,nblks_c)*/
						__global double * z_exner_expl, /*(nproma,nlev,nblks_c)*/
						__global double * exner_ref_mc, /*(nproma,nlev,nblks_c)*/
						__global double * z_beta, /*(nproma,nlev,nblks_c)*/
						__global double * z_alpha, /*(nproma,nlevp1,nblks_c)*/
						__global double * exner_nnew /*(nproma,nlev,nblks_c)*/
)
{
	__private double line_alpha[NLEV_P1];
	__private double line_w_nnew[NLEV_P1];
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);

	// TODO: check this work
	for ( int ilev=0; ilev < nlevp1; ilev++ )
	{
		const int idx2 = jc + ilev*nproma + jb*nproma*nlevp1;
		line_alpha[ilev]  = z_alpha[idx2];
		line_w_nnew[ilev] = w_nnew[idx2];
	}

#if 0
	if (jk < nlevp1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
	
		line_alpha[localIdx]  = z_alpha[idx2];
		line_w_nnew[localIdx] = w_nnew[idx2];
	}
#endif
	
	barrier(CLK_GLOBAL_MEM_FENCE);
	
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + jk*nproma + jb*nproma*nlev;
			
		exner_nnew[idx1] = z_exner_expl[idx1] + exner_ref_mc[idx1] - z_beta[idx1] *
		                  (line_alpha[localIdx] * line_w_nnew[localIdx] - line_alpha[localIdx+nproma] * line_w_nnew[localIdx+nproma]);
	}
}

__kernel void endphase9c(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nproma,
						__private double cvd_o_rd,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * exner_nnew, /*(nproma,nlev,nblks_c)*/
						__global double * exner_nnow, /*(nproma,nlev,nblks_c)*/
						__global double * rhotheta_v_nnow, /*(nproma,nlev,nblks_c)*/
						__global double * rhotheta_v_nnew /*(nproma,nlev,nblks_c)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);

	if (idx < i_endblk*nlev*nproma)
	{
		rhotheta_v_nnew[idx] = rhotheta_v_nnow[idx] *
		                       ((exner_nnew[idx] / exner_nnow[idx] - 1.) * cvd_o_rd + 1.);
	}
}

__kernel void endphase9d(__private int i_startblk,
						__private int i_endblk,
						__private int nlev,
						__private int nproma,
						__global int * localStart, /*(nblks_c)*/
						__global int * localEnd, /*(nblks_c)*/
						__global double * rho_nnew, /*(nproma,nlev,nblks_c)*/
						__global double * rhotheta_v_nnew, /*(nproma,nlev,nblks_c)*/
						__global double * theta_v_nnew /*(nproma,nlev,nblks_c)*/)
{
	const int idx = i_startblk*nlev*nproma + get_global_id(0);

	if (idx < i_endblk*nlev*nproma)
	{
		theta_v_nnew[idx] = rhotheta_v_nnew[idx] / rho_nnew[idx];
	}
}
