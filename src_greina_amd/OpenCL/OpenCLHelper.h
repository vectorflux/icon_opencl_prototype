/*
 *  OpenCLHelper.h
 *
 *  Created by Christian Conti on 06/04/11.
 *
 *	The purpose of this class is to alleviate the programmer
 *		from the burden of setting up the OpenCL framework
 *		and to help him handle the various resources available
 *
 *	Relations between components:
 *		n platforms
 *		n devices        / 1 platform
 *		n contexts       / 1 platform
 *		n devices        / 1 contexts
 *		n command queues / 1 device
 *
 *	If allowed by the device, the command queues are set with:
 *		CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
 *		CL_QUEUE_PROFILING_ENABLE
 *
 *	By default, the following choices are made:
 *		activePlatform = platforms[0]
 *		activeDevices = devices[activePlatform][0] (or the 1st nDevices if specified
 *		activeContext = contexts[0]
 *
 *	The program generation parses for kernel methods
 *		and add them to its list of kernels automatically
 *
 *	FEATURES TO ADD AND NOTES
 *
 *	The class should be able to handle multiple OpenCL platforms
 *		and eventually switch between them.
 *		
 *	Too bulky?
 *
 *	Should it be kept on a single header file?
 *		-> split it and create a library?
 *		tradeoff: readability vs ease of use
 *
 *	All the information could be dumped automatically on file
 *
 *	Is the possiblity to have multiple kernels of the same code useful?
 *		i.e. what happens if the same kernel is used twice (concurrently)?
 *
 */
#pragma once

#include <assert.h>

#include <vector>
#include <string.h>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "CheckCL.h"



class OpenCLHelper
{
private:
	// structures
	vector<cl_platform_id> platforms;
	map<cl_platform_id,vector<cl_device_id> > devices;
	vector<cl_context> contexts; // the platform is already encoded
	map<cl_device_id,vector<cl_command_queue> > commandQueues;
	map<cl_context,map<string,cl_program> > programs;
	map<cl_context,map<string,cl_kernel> > kernels;
	map<string,cl_mem> memObjects;
	
	cl_platform_id activePlatform;
	vector<cl_device_id> activeDevices;
	cl_context activeContext;
	
	void _launchError(bool condition, string message);
	void _launchWarning(bool condition, string message);
	
	void _initPlatforms();
	void _initDevices(int nDevices);
	void _initContexts();
	void _initCommandQueues(int nCommandQueuesPerDevice);
	
	void _DumpInfoOnFile();
	
	string _loadTextFromFile(string filename);
	char* _loadPrecompiledBinaries(string filename, string source, cl_device_id device, string options);
	
public:
	// Constructors
	OpenCLHelper(int nDevices=1, int nCommandQueuesPerDevice=1);
	OpenCLHelper(const OpenCLHelper& o);
	const OpenCLHelper& operator=(OpenCLHelper& o);
	~OpenCLHelper();
	
	// Compile program and parse for kernels (if necessary)
	void readProgram(string programPath, string options="");
	
	//
	void switchPlatform(cl_platform_id newPlatform);
	
	// Getters
	cl_context getContext(int contextId=0);
	cl_command_queue getCommandQueue(int queueId=0, int device=0);
	cl_kernel getKernel(string kernelName);
	
	// Information printers to help debugging of OpenCL code
	void printPlatformInfo();
	void printDeviceInfo(int device);
	void printContextInfo();
	void printCommandQueueInfo(cl_command_queue commandQueue);
	void printProgramInfo(cl_program program);
	void printKernelInfo(string kernelName);
};

#pragma mark Error/Warning handling methods
void OpenCLHelper::_launchError(bool condition, string message)
{
	if (condition)
	{
		cout << "ERROR: " << message << endl;
		abort();
	}
}

void OpenCLHelper::_launchWarning(bool condition, string message)
{
	if (condition)
		cout << "Warning: " << message << endl;
}

#pragma mark Initialization methods for the constructor
void OpenCLHelper::_initPlatforms()
{
	cl_uint nMaxEntries = 10; // in general, there are 1-3 different platforms on a single machine
	cl_platform_id platformsFound[nMaxEntries];
	cl_uint nEntries = -1;
	
	CheckCL( clGetPlatformIDs(nMaxEntries,platformsFound,&nEntries) );
	
	_launchError(nEntries<1, "No OpenCL platform present on the system!");
	
	platforms.resize(0);
	const int nBuf = 800;
	char buf[nBuf];
	char bufName[nBuf];
	for (int i=0; i<nEntries; i++) {
   	    CheckCL(clGetPlatformInfo(platformsFound[i], CL_PLATFORM_NAME, nBuf, bufName, NULL));
     	    cout << endl << "==== OpenCL Platform found: " << bufName;
	    platforms.push_back(platformsFound[i]);
	}
	assert(platforms.size()==nEntries);
	
	activePlatform = platforms[0];
	
	CheckCL(clGetPlatformInfo(activePlatform, CL_PLATFORM_NAME, nBuf, bufName, NULL));
	cout << endl << "==== OpenCL Platform: " << bufName;
	CheckCL(clGetPlatformInfo(activePlatform, CL_PLATFORM_VERSION, nBuf, buf, NULL));
	cout << " (" << buf << ") ====" << endl;
}

void OpenCLHelper::_initDevices(int nDevices)
{
	const int nMaxEntries=100;
	cl_device_id deviceIds[nMaxEntries];
	cl_uint nEntries = 0;
	
	vector<cl_platform_id>::iterator platformIter;
	for (platformIter = platforms.begin(); platformIter != platforms.end(); platformIter++)
	{
		CheckCL( clGetDeviceIDs(*platformIter, CL_DEVICE_TYPE_GPU, nMaxEntries, deviceIds, &nEntries) );
		printf("*********");
		_launchError(nEntries<1, "No OpenCL device present on the system!");
		
		devices[*platformIter].resize(0);
		for (int i=0; i<nEntries; i++)
			devices[*platformIter].push_back(deviceIds[i]);
		assert(devices[*platformIter].size()==nEntries);
	}
	
	activeDevices.resize(0);
	for (int i=0; i<nDevices && i<devices[activePlatform].size(); i++)
		activeDevices.push_back(devices[activePlatform][i]);
	assert(activeDevices.size()<=nDevices);
	
	_launchWarning(activeDevices.size()<nDevices,"not enough devices found, using all devices found for the current platform");
	
	const int nBuf = 800;
	char bufName[nBuf];
	CheckCL( clGetDeviceInfo(activeDevices[0], CL_DEVICE_NAME, nBuf, bufName, NULL) );
	cout << "==== Running on: " << bufName << " ====" << endl << endl;
}

void OpenCLHelper::_initContexts()
{
	// does it make sense to use clContextFromType? is it safer?
	cl_int status;
	cl_context_properties properties[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)activePlatform, 0 };
	contexts.resize(0);
	contexts.push_back(clCreateContext(properties, activeDevices.size(), &activeDevices.front(), NULL, NULL, &status));
	CheckCL(status);
	assert(contexts.size()==1);
	
	activeContext = contexts[0];
}

void OpenCLHelper::_initCommandQueues(int nCommandQueuesPerDevice)
{
	cl_int status;
	
	cl_command_queue_properties properties;
	for (int i=0; i<activeDevices.size(); i++)
	{
		// Get for supported command queue properties
		clGetDeviceInfo(activeDevices[i], CL_DEVICE_QUEUE_PROPERTIES, sizeof(cl_command_queue_properties), &properties, NULL);
		_launchWarning(!((properties|CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)&CL_QUEUE_PROFILING_ENABLE),"This device does not satisfy the OpenCL 1.1 specifications!");
		
		commandQueues[activeDevices[i]].resize(0);
		for (int n=0; n<nCommandQueuesPerDevice; n++)
		{
			commandQueues[activeDevices[i]].push_back(clCreateCommandQueue(activeContext,activeDevices[i],properties,&status));
			CheckCL(status);
		}
		assert(commandQueues[activeDevices[i]].size()==nCommandQueuesPerDevice);
	}
}

#pragma mark System information dumper
void OpenCLHelper::_DumpInfoOnFile()
{
/*
	int bufSize = 1024;
	char * buf[bufSize];
	
	fstream filestream("/users/cconti/SystemInfo.txt");
	_launchError(filestream==NULL,"Cannot open file to write system informations.");
	
	assert(filestream.is_open());
	
	filestream << "=========================================================================" << endl;
	filestream << endl << "  Platforms information:\n";
	for (int p=0; p<platforms.size(); p++)
	{
		filestream << "\tPlatform " << p << endl;
		
		clGetPlatformInfo(platforms[p], CL_PLATFORM_NAME, bufSize, buf, NULL);
		filestream << "\t\tName:\t\t" << string(buf) << endl;
		
		clGetPlatformInfo(platforms[p], CL_PLATFORM_VENDOR, bufSize, buf, NULL);
		filestream << "\t\tVendor:\t\t" << string(buf) << endl;
		
		clGetPlatformInfo(platforms[p], CL_PLATFORM_VERSION, bufSize, buf, NULL);
		filestream << "\t\tVersion:\t" << string(buf) << endl;
		
		clGetPlatformInfo(platforms[p], CL_PLATFORM_EXTENSIONS, bufSize, buf, NULL);
		filestream << "\t\tExtensions:\t" << string(buf) << endl;
	}
	filestream << "=========================================================================" << endl;

	filestream.close();
*/
}

#pragma mark Program builder helpers
string OpenCLHelper::_loadTextFromFile(string filename)
{
	fstream filestream(filename.c_str());
	_launchError(filestream==NULL,"Cannot open program file " + filename);
	
	string text;
	while (filestream)
	{
		string buf;
		getline(filestream,buf);
		text += buf;
		text += '\n';
	}
	
	return text;
}

char* OpenCLHelper::_loadPrecompiledBinaries(string filename, string source, cl_device_id device, string options)
{
	/*
	fstream filestream(filename.c_str());
	if (filestream==NULL)
		return NULL;
	
	string key = _loadTextFromFile(filename);
	*/
	
	_launchWarning(true,"currently not caching programs on file!");
	return NULL;
}

#pragma mark Constructors
OpenCLHelper::OpenCLHelper(int nDevices, int nCommandQueuesPerDevice)
{
	// setting platforms informations, set first to active
	_initPlatforms();
	
	// setting devices informations
	_initDevices(nDevices);
	
	// create a context for the active devices and the activ e platform
	_initContexts();
	
	// create command queues for the active devices
	_initCommandQueues(nCommandQueuesPerDevice);
	
	// dump system informations
	_DumpInfoOnFile();
	
	//printPlatformInfo();
	//printDeviceInfo(0);
	//printContextInfo();
	//printCommandQueueInfo();
}

OpenCLHelper::OpenCLHelper(const OpenCLHelper& o) :
platforms(o.platforms),devices(o.devices),contexts(o.contexts),
commandQueues(o.commandQueues),programs(o.programs),kernels(o.kernels),
activePlatform(o.activePlatform),activeDevices(o.activeDevices),activeContext(o.activeContext),
memObjects(o.memObjects)
{
}

const OpenCLHelper& OpenCLHelper::operator=(OpenCLHelper& o)
{
	platforms = o.platforms;
	devices = o.devices;
	contexts = o.contexts;
	commandQueues = o.commandQueues;
	programs = o.programs;
	kernels = o.kernels;
	activePlatform = o.activePlatform;
	activeDevices = o.activeDevices;
	activeContext = o.activeContext;
	memObjects = o.memObjects;
	
	return *this;
}

OpenCLHelper::~OpenCLHelper()
{
}

#pragma mark Program
void OpenCLHelper::readProgram(string programPath, string options)
{
	cl_int status;
	cl_program program;
	char sProgramFullName[4096];
	sprintf(sProgramFullName, "%s%s", programPath.c_str(), options.c_str());
	
	if (programs[activeContext].find(sProgramFullName) == programs[activeContext].end())
	{
		string source = _loadTextFromFile(programPath);
		_launchError(source=="\n",programPath + " is empty\n");
		
		char** binaries = new char*[activeDevices.size()];
		for (int i=0; i<activeDevices.size(); i++)
			binaries[i] = _loadPrecompiledBinaries(sProgramFullName, source, activeDevices[i], options);
		
		_launchWarning(true,"no caching of program source is implemented yet!");
		
		const char* src = source.data();
		size_t sourceSize = (size_t)source.size();

		program = clCreateProgramWithSource(activeContext, 1, &src, &sourceSize, &status);
		CheckCL(status);
		
		if(options.size() != 0)
			std::cout << "Build Options are : " << options.c_str() << std::endl;
		
		CheckRelaxCL(status = clBuildProgram(program, activeDevices.size(), &activeDevices.front(), options.c_str(), NULL, NULL));
		if (status != CL_SUCCESS)
		{
			cout << "\n====================================================================\n";
			cout << "\tDumping bugged source code\n\n";
			
			for (int i=0; i<activeDevices.size(); i++)
			{
				cout << "\t\tDevice " << i << endl;
				
				size_t bufferSize = 5*source.size();
				char * buf = new char[bufferSize];
				size_t effectiveSize = 0;
				
				CheckCL(clGetProgramBuildInfo(program, activeDevices[i], CL_PROGRAM_BUILD_LOG, bufferSize, buf, &effectiveSize));
				
				cout << string(buf) << endl << endl;
				
				delete [] buf;
			}
			
			cout << "\n\tDump completed!";
			cout << "\n====================================================================\n";
			abort();
		}
		// store the data in a file!
		
		programs[activeContext][sProgramFullName] = program;
		
		// parse for kernels
		size_t pos = 0, pos2 = 0;
		string parsedKernelName = "";
		const int finalPos = source.rfind("kernel");
		
		while(pos<finalPos)
		{
			pos = source.find("kernel",pos) + 6;
			pos = source.find("void ",pos) + 5;
			pos2 = source.find("(",pos);
			parsedKernelName = source.substr(pos,pos2-pos);
			//_launchWarning(true,"Parsed Kernel Name: " + parsedKernelName);
			
			kernels[activeContext][parsedKernelName] = clCreateKernel(program,parsedKernelName.c_str(),&status);
			CheckCL(status);
		}
	}
}

#pragma mark Switch Platform
void OpenCLHelper::switchPlatform(cl_platform_id newPlatform)
{
	activePlatform = newPlatform;
	_launchError(true,"This method is incomplete - changing platform means changing context, devices and queues as well");
}

#pragma mark Getters
cl_context OpenCLHelper::getContext(int contextId)
{
	_launchError(contextId>=contexts.size(),"This context does not exist!");
	return contexts[contextId];
}

cl_command_queue OpenCLHelper::getCommandQueue(int queueId, int device)
{
	_launchError(activeDevices.size()<=device, "This device does not exist or is not active.");
	_launchError(commandQueues[activeDevices[device]].size()<=queueId, "The requested queue does not exist.");
	
	return commandQueues[activeDevices[device]][queueId];
}

cl_kernel OpenCLHelper::getKernel(string kernelName)
{
	_launchError(kernels[activeContext].find(kernelName)==kernels[activeContext].end(),"Kernel '" + kernelName + "' does not exist in the active context.");
	return kernels[activeContext][kernelName];
}

#pragma mark Information printers
void OpenCLHelper::printPlatformInfo()
{
	const int bufSize = 1024;
	char buf[bufSize];
	
	cout << "\n===================================================================\n\n";
	cout << "  PLATFORM INFORMATION\n\n";
	
	clGetPlatformInfo(activePlatform, CL_PLATFORM_NAME, bufSize, buf, NULL);
	cout << "\tName:\t\t" << string(buf) << endl;
	
	clGetPlatformInfo(activePlatform, CL_PLATFORM_VENDOR, bufSize, buf, NULL);
	cout << "\tVendor:\t\t" << string(buf) << endl;
	
	//clGetPlatformInfo(activePlatform, CL_PLATFORM_PROFILE, bufSize, buf, NULL);
	//cout << "\tProfile:\t" << string(buf) << endl;
	
	clGetPlatformInfo(activePlatform, CL_PLATFORM_VERSION, bufSize, buf, NULL);
	cout << "\tVersion:\t" << string(buf) << endl;
	
	clGetPlatformInfo(activePlatform, CL_PLATFORM_EXTENSIONS, bufSize, buf, NULL);
	cout << "\tExtensions:\t" << string(buf) << endl;

	cout << "\n===================================================================\n\n";
}

void OpenCLHelper::printDeviceInfo(int device)
{
	const int bufSize = 1024;
	char buf[bufSize];
	cl_device_type type;
	cl_uint uintVal;
	
	cout << "\n===================================================================\n\n";
	for (int i=0; i<activeDevices.size(); i++)
	{
		cout << "  DEVICE " << i << " INFORMATION\n\n";
		
		CheckCL( clGetDeviceInfo(activeDevices[i], CL_DEVICE_NAME, bufSize, buf, NULL) );
		cout << "\tDevice Name:\t\t" << string(buf) << endl;
		//CheckCL( clGetDeviceInfo(activeDevices[i], CL_DEVICE_PROFILE, bufSize, buf, NULL) );
		//cout << "\tDevice Profile:\t" << string(buf) << endl;
		CheckCL( clGetDeviceInfo(activeDevices[i], CL_DEVICE_VERSION, bufSize, buf, NULL) );
		cout << "\tDevice Version:\t\t" << string(buf) << endl;
		

		clGetDeviceInfo(activeDevices[i], CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL);
		cout << "\tDevice Type\n";
		
		clGetDeviceInfo(activeDevices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &uintVal, NULL);
		cout << "\tMax Compute Units:\t" << uintVal << endl;
		
		clGetDeviceInfo(activeDevices[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &uintVal, NULL);
		cout << "\tMax Dim for Work-items:\t" << uintVal << endl;
		
		size_t * maxWorkItemSizes = new size_t[3];
		clGetDeviceInfo(activeDevices[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*uintVal, maxWorkItemSizes, NULL);
		cout << "\tMaximum Work-items:\t" << maxWorkItemSizes[0];
		if (uintVal>1) cout << " " << maxWorkItemSizes[1];
		if (uintVal>2) cout << " " << maxWorkItemSizes[2];
		cout << endl;
		
		CheckCL( clGetDeviceInfo(activeDevices[i], CL_DEVICE_EXTENSIONS, bufSize, buf, NULL) );
		cout << "\tDevice Extensions:\t" << buf << endl;
		
		cl_ulong ulongVal;
		clGetDeviceInfo(activeDevices[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &ulongVal, NULL);
		cout << "\tMax Memory Alloc Size (MB):\t" << ulongVal/1024./1024. << endl;
		
		clGetDeviceInfo(activeDevices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &ulongVal, NULL);
		cout << "\tGlobal Memory Size (MB):\t" << ulongVal/1024./1024. << endl;
		
	}
	cout << "\n===================================================================\n\n";
}

void OpenCLHelper::printContextInfo()
{
}

void OpenCLHelper::printCommandQueueInfo(cl_command_queue commandQueue)
{
}

void OpenCLHelper::printProgramInfo(cl_program program)
{
}

void OpenCLHelper::printKernelInfo(string kernelName)
{
	const int bufSize = 1024;
	cl_uint n;
	char buf[bufSize];
	cout << "\n===================================================================\n\n";
	cout << "  KERNEL " << kernelName << " INFORMATION\n\n";
	CheckCL(clGetKernelInfo(kernels[activeContext][kernelName], CL_KERNEL_FUNCTION_NAME, sizeof(char)*bufSize, buf, NULL));
	cout << "\tName:\t\t" << string(buf) << endl;
	CheckCL(clGetKernelInfo(kernels[activeContext][kernelName], CL_KERNEL_NUM_ARGS, sizeof(cl_uint), &n, NULL));
	cout << "\tNumber of Arguments:\t" << n << endl;
	cout << "\n===================================================================\n\n";
}
