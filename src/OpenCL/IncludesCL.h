#pragma once

#include "CL/cl.h"


#include "OpenCLHelper.h"
static OpenCLHelper oclhelper;

//#include "EngineCL.h"
#include "ResourcesCL.h"
#include "ProfilerCL.h"

//static EngineCL engineCL(1);
static ResourcesCL resourcesCL;
static ProfilerCL profilerCL;