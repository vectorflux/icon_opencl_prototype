//==================================================
//
//	Kernels for the components of velocity_tendencies
//		(mo_solve_nonhydro.f90)
//
//==================================================


#ifndef __Cayman__
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

__kernel void step1ia(__private int i_startblk,
                      __private int nblks_e,
                      __private int nflat,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_e)*/
					  __global int * localEnd, /*(nblks_e)*/
					  __global double * wgtfac_e, /*(nproma,nlevp1,nblks_e)*/
					  __global double * vt, /*(nproma,nlev,nblks_e)*/
					  __global double * z_vt_ie /*(nproma,nlevp1,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = nflat + get_global_id(1);
		
	if (jk < nlev && je < localEnd[get_global_id(0)])
	{
		const int idxp1  = je + jk*nproma + jb*nproma*nlevp1;
		const int idx_vt = je + jk*nproma + jb*nproma*nlev;
		const double wgtfac_e_tmp = wgtfac_e[idxp1];
		z_vt_ie[idxp1] = wgtfac_e_tmp    *vt[idx_vt] +
						(1.-wgtfac_e_tmp)*vt[idx_vt-nproma];
	}
}

__kernel void step1ib(__private int i_startblk,
                      __private int nblks_e,
                      __private int nflat,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_e)*/
					  __global int * localEnd, /*(nblks_e)*/
					  __global double * wgtfacq_e, /*(nproma,3,nblks_e)*/
					  __global double * vt, /*(nproma,nlev,nblks_e)*/
					  __global double * z_vt_ie /*(nproma,nlevp1,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(1);
		
	if (je < localEnd[get_global_id(0)])
	{
		z_vt_ie[je+(nlevp1-1)*nproma + jb*nproma*nlevp1] = wgtfacq_e[je         +jb*nproma*3]*vt[je+(nlev-1)*nproma+jb*nproma*nlev] +
		                                                   wgtfacq_e[je+nproma  +jb*nproma*3]*vt[je+(nlev-2)*nproma+jb*nproma*nlev] +
		                                                   wgtfacq_e[je+nproma*2+jb*nproma*3]*vt[je+(nlev-3)*nproma+jb*nproma*nlev];
	}
}

__kernel void step1ic(__private int i_startblk,
                      __private int i_endblk,
                      __private int nflat,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_e)*/
					  __global int * localEnd, /*(nblks_e)*/
					  __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					  __global double * z_vt_ie, /*(nproma,nlevp1,nblks_e)*/
					  __global double * ddxn_z_half, /*(nproma,nlevp1,nblks_e)*/
					  __global double * ddxt_z_half, /*(nproma,nlevp1,nblks_e)*/
					  __global double * z_concorr_e /*(nproma,nlevp1,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = nflat + get_global_id(1);
		
	if (jb<i_endblk && jk < nlevp1 && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jk*nproma + jb*nproma*nlevp1;
		
		z_concorr_e[idx] = vn_half[idx]*ddxn_z_half[idx] +
		                   z_vt_ie[idx]*ddxt_z_half[idx];
	}
}

__kernel void step1a(__private int i_startblk,
                     __private int i_endblk,
                     __private int nlev,
					 __private int nlevp1,
					 __private int nproma,
					 __global int * localStart, /*(nblks_e)*/
					 __global int * localEnd, /*(nblks_e)*/
					 __global double * wgtfac_e, /*(nproma,nlevp1,nblks_e)*/
					 __global double * vn, /*(nproma,nlev,nblks_e)*/
					 __global double * vt, /*(nproma,nlev,nblks_e)*/
					 __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					 __global double * z_kin_hor_e, /*(nproma,nlev,nblks_e)*/
					 __local double line_vn[NLEV_])
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = 1 + get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_global_id(1) < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx1 = je + (get_global_id(1))*nproma + jb*nproma*nlev;
	
		line_vn[localIdx] = vn[idx1];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
	
	if (jk < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx    = je + jk*nproma + jb*nproma*nlevp1;
		const int idx_vn = je + jk*nproma + jb*nproma*nlev;
		
		const double wgtfac_e_tmp = wgtfac_e[idx];
		const double vn_tmp = line_vn[localIdx+nproma];
		const double vt_tmp = vt[idx_vn];
		
		vn_half[idx] = wgtfac_e_tmp*vn_tmp + (1.-wgtfac_e_tmp)*line_vn[localIdx];
		z_kin_hor_e[idx_vn] = .5*(vn_tmp*vn_tmp + vt_tmp*vt_tmp);
	}
}

__kernel void step1b(__private int i_startblk,
                     __private int i_endblk,
                     __private int nlev,
                     __private int nlevp1,
					 __private int nproma,
					 __private int l_vert_nested,
					 __global int * localStart, /*(nblks_e)*/
					 __global int * localEnd, /*(nblks_e)*/
					 __global double * wgtfacq1_e, /*(nproma,3,nblks_e)*/
					 __global double * wgtfacq_e, /*(nproma,3,nblks_e)*/
					 __global double * vn, /*(nproma,nlev,nblks_e)*/
					 __global double * vt, /*(nproma,nlev,nblks_e)*/
					 __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					 __global double * z_kin_hor_e /*(nproma,nlev,nblks_e)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(1);
	
	if (jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jb*nproma*nlevp1;
		const int idx_wgt = je + jb*nproma*3;
		const int idx_vn = je + jb*nproma*nlev;
		
		const double vn_tmp = vn[idx_vn];
		const double vt_tmp = vt[idx_vn];
		
		// might be better to put the conditional on host side
		//	as l_vert_nested is constant across the threads
		if (l_vert_nested==0)
			vn_half[idx] = wgtfacq1_e[idx_wgt         ]*vn_tmp +
						   wgtfacq1_e[idx_wgt+nproma  ]*vn[idx_vn + nproma  ] +
					       wgtfacq1_e[idx_wgt+nproma*2]*vn[idx_vn + nproma*2];
		
		z_kin_hor_e[idx_vn] = .5*(vn_tmp*vn_tmp+vt_tmp*vt_tmp);
		
		vn_half[idx+(nlevp1-1)*nproma] = wgtfacq_e[idx_wgt         ]*vn[idx_vn + (nlev-1)*nproma] +
		                                 wgtfacq_e[idx_wgt+nproma  ]*vn[idx_vn + (nlev-2)*nproma] +
					                     wgtfacq_e[idx_wgt+nproma*2]*vn[idx_vn + (nlev-3)*nproma];
	}
}

__kernel void stepn(__private int i_startblk,
                    __private int i_endblk,
                    __private int nlev,
					__private int nproma,
					__global int * localStart, /*(nblks_e)*/
					__global int * localEnd, /*(nblks_e)*/
					__global double * vn, /*(nproma,nlev,nblks_e)*/
					__global double * vt, /*(nproma,nlev,nblks_e)*/
					__global double * z_kin_hor_e /*(nproma,nlev,nblks_e)*/)
{
	const int idx = i_startblk*nproma*nlev + get_global_id(0);
		
	if (idx < i_endblk*nproma*nlev)	
	//if (jk < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		//const int idx = je + jk*nproma + jb*nproma*nlev;
		const double vt_tmp = vt[idx];
		const double vn_tmp = vn[idx];
		z_kin_hor_e[idx] = .5*(vn_tmp*vn_tmp+vt_tmp*vt_tmp);
		
	}
}

__kernel void phase2(__private int i_startblk,
                     __private int i_endblk,
					 __private int nblocks_e,
                     __private int slev,
					 __private int nlev,
					 __private int nlevp1,
					 __private int nproma,
					 __global int * localStart, /*(nblks_e)*/
					 __global int * localEnd, /*(nblks_e)*/
					 __global int * icidx, /*(nproma,nblks_e,2)*/
					 __global int * icblk, /*(nproma,nblks_e,2)*/
					 __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					 __global double * c_lin_e, /*(nproma,2,nblks_e)*/
					 __global double * w, /*(nproma,nlevp1,nblks_c)*/
					 __global double * inv_dual_edge_length, /*(nproma,nblks_e)*/
					 __global double * e_kinh, /*(nproma,nlev,nblks_c)*/
					 __global double * z_vnw, /*(nproma,nlevp1,nblks_e)*/
					 __global double * z_ddxn_ekin_e, /*(nproma,nlev,nblks_e)*/
					 __local int line_idx[2],
					 __local int line_idx2[2],
					 __local double line_c_lin_e[2])
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	
	if (get_local_id(1) < 2 && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx = je + jb*nproma + get_local_id(1)*nproma*nblocks_e;
		
		line_idx[localIdx]  = icidx[idx]-1 + (icblk[idx]-1)*nproma*nlev;
		line_idx2[localIdx] = icidx[idx]-1 + (icblk[idx]-1)*nproma*nlevp1;
		line_c_lin_e[localIdx] = c_lin_e[je+jb*nproma*2+get_local_id(1)*nproma];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
		
	if (jk < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx_z_vnw = je + jk*nproma + jb*nproma*nlevp1;
		const int idx_z_ddxn_ekin_e = je + jk*nproma + jb*nproma*nlev;
		
		const int idx_w1 = line_idx2[0] + jk*nproma;
		const int idx_w2 = line_idx2[1] + jk*nproma;
		
		z_vnw[idx_z_vnw] = vn_half[idx_z_vnw]*
		                  (line_c_lin_e[0]*w[idx_w1] +
						   line_c_lin_e[1]*w[idx_w2]);
		
		const int idx_e_kinh1 = line_idx[0] + jk*nproma;
		const int idx_e_kinh2 = line_idx[1] + jk*nproma;
		
		z_ddxn_ekin_e[idx_z_ddxn_ekin_e] = inv_dual_edge_length[je+jb*nproma] * (e_kinh[idx_e_kinh2] - e_kinh[idx_e_kinh1]);
	}
}

__kernel void phase3a(__private int i_startblk,
                      __private int i_endblk,
					  __private int nblks_c,
					  __private int slev,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*(nblks_e)*/
					  __global int * localEnd, /*(nblks_e)*/
					  __global int * ieidx, /*(nproma,nblks_c,3)*/
					  __global int * ieblk, /*(nproma,nblks_c,3)*/
					  __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					  __global double * geofac_div, /*(nproma,3,nblks_e)*/
					  __global double * w, /*(nproma,nlevp1,nblks_c)*/
					  __global double * z_vnw, /*(nproma,nlevp1,nblks_e)*/
					  __global double * z_hadv_w, /*(nproma,nlevp1,nblks_c)*/
					  __global double * w_concorr_c, /*(nproma,nlevp1,nblks_c)*/
					  __global double * w_con, /*(nproma,nlevp1,nblks_c)*/
					  __local int line_idx[3],
					  __local double line_geofac_div[3])
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_local_id(1) < 3 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jb*nproma + get_local_id(1)*nproma*nblks_c;
		
		line_idx[localIdx] = ieidx[idx]-1 + (ieblk[idx]-1)*nproma*nlevp1;
		line_geofac_div[localIdx] = geofac_div[jc+jb*nproma*3+get_local_id(1)*nproma];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
		
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx_out = jc + jk*nproma + jb*nproma*nlevp1;
		const int idx_1 = line_idx[0] + jk*nproma;
		const int idx_2 = line_idx[1] + jk*nproma;
		const int idx_3 = line_idx[2] + jk*nproma;
		
		const double w_tmp = w[idx_out];
		const double geofac_div0 = line_geofac_div[0];
		const double geofac_div1 = line_geofac_div[1];
		const double geofac_div2 = line_geofac_div[2];
		
		z_hadv_w[idx_out] = w_tmp*(vn_half[idx_1]*geofac_div0 +
								   vn_half[idx_2]*geofac_div1 +
								   vn_half[idx_3]*geofac_div2) -
							      (z_vnw[idx_1]*geofac_div0 +
								   z_vnw[idx_2]*geofac_div1 +
								   z_vnw[idx_3]*geofac_div2);
		w_con[idx_out] = w_tmp - (nlev>slev ? w_concorr_c[idx_out] : 0.);
	}
}

__kernel void phase3b(__private int i_startblk,
                      __private int i_endblk,
					  __private int slev,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*nblks_c*/
					  __global int * localEnd, /*nblks_c*/
					  __global double * w_concorr_c, /*(nproma,nlevp1,nblks_c)*/
					  __global double * w_con /*(nproma,nlevp1,nblks_c)*/)
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = slev + get_global_id(1);
		
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + jk*nproma + jb*nproma*nlevp1;
		w_con[idx] -= w_concorr_c[idx];
	}
}

__kernel void phase3c(__private int i_startblk,
                      __private int i_endblk,
					  __private int nlev,
					  __private int nlevp1,
					  __private int nproma,
					  __global int * localStart, /*nblks_c*/
					  __global int * localEnd, /*nblks_c*/
					  __global double * w_con, /*(nproma,nlevp1,nblks_c)*/
					  __global double * z_w_con_c_full, /*(nproma,nlev,nblks_c)*/
					  __local double line_w_con[NLEV_P1])
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_local_id(1) < nlevp1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx = jc + get_local_id(1)*nproma + jb*nproma*nlevp1;
		
		line_w_con[localIdx] = w_con[idx];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
		
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx_out = jc + jk*nproma + jb*nproma*nlev;
		
		z_w_con_c_full[idx_out] = .5 * (line_w_con[localIdx] + line_w_con[localIdx+1]);
	}
}

__kernel void phase4(__private int i_startblk,
					 __private int i_endblk,
					 __private int nlev,
					 __private int nlevp1,
					 __private int nproma,
					 __private int ntnd,
					 __private int nblks_e,
					 __global int * localStart, /*nblks_c*/
					 __global int * localEnd, /*nblks_c*/
					 __global int * ividx, /*(nproma,nblks_e,4)*/
					 __global int * ivblk, /*(nproma,nblks_e,4)*/
					 __global int * icidx, /*(nproma,nblks_e,2)*/
					 __global int * icblk, /*(nproma,nblks_e,2)*/
					 __global double * z_ddxn_ekin_e, /*(nproma,nlev,nblks_e)*/
					 __global double * vt, /*(nproma,nlev,nblks_e)*/
					 __global double * f_e, /*(nproma,nblks_e)*/
					 __global double * omega_z, /*(nproma,nlev,nblks_v or nblks_e)*/
					 __global double * c_lin_e, /*(nproma,2,nblks_e)*/
					 __global double * z_w_con_c_full, /*(nproma,nlev,nblks_c)*/
					 __global double * vn_half, /*(nproma,nlevp1,nblks_e)*/
					 __global double * ddqz_z_full_e, /*(nproma,nlev,nblks_e)*/
					 __global double * ddt_vn_adv, /*(nproma,nlev,nblks_e,3)*/
					 __local int line_ividx[2],
					 __local int line_icidx[2],
					 __local double line_c_lin_e[2],
					 __local double line_vn_half[NLEV_P1])
{
	const int jb = i_startblk + get_global_id(0);
	const int je = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (jk < nlevp1 && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx2 = je + jk*nproma + jb*nproma*nlevp1;
		line_vn_half[localIdx] = vn_half[idx2];
		
		if (jk < 2)
		{
			const int idx = je + jb*nproma+get_local_id(1)*nproma*nblks_e;
			line_ividx[localIdx] = ividx[idx]-1 + (ivblk[idx]-1)*nproma*nlev;
			line_icidx[localIdx] = icidx[idx]-1 + (icblk[idx]-1)*nproma*nlev;
			line_c_lin_e[localIdx] = c_lin_e[je+get_local_id(1)*nproma+jb*nproma*2];
		}
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
		
	if (jk < nlev && jb < i_endblk && je < localEnd[get_global_id(0)])
	{
		const int idx_out = je + jk*nproma + jb*nproma*nlev + (ntnd-1)*nproma*nlev*nblks_e;
		const int idx_in1 = je + jk*nproma + jb*nproma*nlev;
		const int idx_omega_z1 = line_ividx[0] + jk*nproma;
		const int idx_omega_z2 = line_ividx[1] + jk*nproma;
		const int idx_z_w_con_c_full1 = line_icidx[0] + jk*nproma;
		const int idx_z_w_con_c_full2 = line_icidx[1] + jk*nproma;
		
		ddt_vn_adv[idx_out] = -(z_ddxn_ekin_e[idx_in1] +
		                      vt[idx_in1] * (f_e[je + jb*nproma] + .5 *
							  (omega_z[idx_omega_z1] + omega_z[idx_omega_z2])) +
							  (line_c_lin_e[0]*z_w_con_c_full[idx_z_w_con_c_full1] +
							   line_c_lin_e[1]*z_w_con_c_full[idx_z_w_con_c_full2])*
							  (line_vn_half[localIdx] - line_vn_half[localIdx+nproma]) / ddqz_z_full_e[idx_in1]);
	}
}

__kernel void phase5(__private int i_startblk,
					 __private int i_endblk,
					 __private int nlev,
					 __private int nlevp1,
					 __private int nproma,
					 __private int ntnd,
					 __private int nblks_c,
					 __global int * localStart, /*nblks_c*/
					 __global int * localEnd, /*nblks_c*/
					 __global double * w_con, /*(nproma,nlevp1,nblks_c)*/
					 __global double * w, /*(nproma,nlevp1,nblks_c)*/
					 __global double * inv_ddqz_z_half2, /*(nproma,nlev,nblks_c)*/
					 __global double * ddt_w_adv, /*(nproma,nlevp1,nblks_c,3)*/
					 __local double line_w[NLEV_P1])
{
	const int jb = i_startblk + get_global_id(0);
	const int jc = localStart[get_global_id(0)] + get_global_id(2);
	const int jk = 1 + get_global_id(1);
	
	const int localIdx = get_local_id(0) + get_local_id(1)*get_local_size(0) + get_local_id(2)*get_local_size(0)*get_local_size(1);
	
	if (get_global_id(1) < nlevp1 && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + (get_global_id(1))*nproma + jb*nproma*nlevp1;
	
		line_w[localIdx] = w[idx1];
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
		
	if (jk < nlev && jb < i_endblk && jc < localEnd[get_global_id(0)])
	{
		const int idx1 = jc + jk*nproma + jb*nproma*nlev;
		const int idx2 = jc + jk*nproma + jb*nproma*nlevp1;
		
		ddt_w_adv[idx2] -= w_con[idx2]*(line_w[localIdx]-line_w[localIdx+2])*inv_ddqz_z_half2[idx1];
	}
}