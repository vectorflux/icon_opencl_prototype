
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
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri,15, sizeof(int)*6, NULL) );
	CheckCL( clSetKernelArg(ocl.kernel_rot_vertex_atmos_tri,16, sizeof(double)*6, NULL) );

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
