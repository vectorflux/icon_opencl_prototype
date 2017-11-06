/*
 *  ResourcesCL.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 11/17/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#include <string.h>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

#include "CheckCL.h"

using namespace std;

class ResourcesCL
{
	int bufferTotalSize;
	map<string, vector<cl_mem> > m_bufObjects;
	map<string, vector<cl_mem> > m_imgObjects;
	map<string, cl_program> m_progObjects;
	map< cl_program,  map<string, cl_kernel> > m_kernelObjects;
	map<string, cl_sampler> m_samplerObjects;

	template<typename T>
	bool _getRes(const map<string, T>& resources, string sKey, T& result)
	{
		typename map<string, T>::const_iterator it = resources.find(sKey);

		if (it == resources.end()) return false;

		//cout<<"Cache hit!"<<endl;

		result = it->second;

		return true;
	}

	inline string _loadTextFileCpp(string sFile)
	{
		string result;

		ifstream in_file( sFile.data() );
		if( !in_file ) {
			cerr << "Couldnot open input file. Path:"<< sFile << endl;
			abort();
		}


		string tmp;
		while(!in_file.eof()) {
			getline(in_file, tmp);
			result += tmp;
			result += "\n";
		}

	//	cout<< "read content:"<< result << endl<< "=============================================="<<endl;
	//	cout<< result.size() << "characters" <<endl;
		return result;
	}
	
	inline string _loadTextFile(string sFile)
	{
		string result;
		
		FILE * f = fopen( sFile.data(), "r" );
		if( f == NULL ) {
			cerr << "Could not open input file. Path:"<< sFile << endl;
			abort();
		}
		
		
		//string tmp;
		char buf[1000];
		while(!feof(f)) {
			memset(buf, 1000, 0);
			if (fgets(buf,1000, f)!=NULL)
				result += string(buf);
			//result += "\n";
		}
		
		fclose(f);
		///cout<< "read content:"<< result << endl<< "=============================================="<<endl;
		//cout<< result.size() << "characters" <<endl;
		return result;
	}

public:
	void printTotalBufferSize()
	{
		const double sizeMB = bufferTotalSize/1024./1024.;
		if (sizeMB > 1000)
			cout << "==== Warning: space allocated for buffers on device is " << sizeMB << " MB ====\n";
	}

	vector<cl_mem> getBuffers(const char* sKey, size_t howmany, cl_context context, cl_mem_flags flags, size_t size)
	{
		vector<cl_mem> result;

		char sKeyDetails[1000];
		sprintf(sKeyDetails, "%s%d", sKey, size);

		if (_getRes(m_bufObjects, sKeyDetails, result))
			if (result.size()>= howmany)
			{
#ifdef TEST
				//cout << "cache hit (" << sKeyDetails << ")\n";
#endif
				return result;
			}
			else 
				howmany -= result.size();

		for(unsigned int i = 0; i < howmany; i++)
		{
			cl_int success;
			cl_mem obj = clCreateBuffer(context, flags, size, NULL, &success);
			if (success != CL_SUCCESS) cout << "checking buffer " << sKeyDetails << endl;
			CheckCL( success );

			result.push_back(obj);
		}
		
		if (m_bufObjects.size()==0)
			bufferTotalSize = 0;
		bufferTotalSize += size;

		m_bufObjects[sKeyDetails] = result;
		printTotalBufferSize();
		
#ifdef TEST
		//cout << "cache miss (" << sKeyDetails << ")\n";
#endif
		
		return result;
	}

	vector<cl_mem> getImages(char* sKey, size_t howmany, cl_context context, cl_mem_flags flags, cl_channel_order order, cl_channel_type channel_type, size_t sizeX, size_t sizeY)
	{
		vector<cl_mem> result;

		char sKeyDetails[1000];
		sprintf(sKeyDetails, "%s%dx%d", sKey, sizeX, sizeY );

		if (_getRes(m_imgObjects, sKeyDetails, result))
			if (result.size()>= howmany)
				return result;
			else 
				howmany -= result.size();

		cl_image_format format;
		format.image_channel_order = order;
		format.image_channel_data_type = channel_type;

		for(unsigned int i=0; i < howmany; i++)
		{
			cl_int success;
			cl_mem mem = clCreateImage2D(context, flags, &format, sizeX, sizeY, 0, NULL, &success);
			CheckCL(success);

			result.push_back(mem);
		}

		m_imgObjects[sKeyDetails] = result;

		return result;
		
	}

	inline cl_sampler getSampler(string sKey, cl_context context, cl_bool bNormalizedCoordinates, cl_addressing_mode address, cl_filter_mode filtering)
	{
		cl_sampler result;

		if (_getRes(m_samplerObjects, sKey, result))
			return result;
		else
		{
			cl_int success;
			result = clCreateSampler(context, bNormalizedCoordinates, address, filtering, &success);
			CheckCL(success);

			m_samplerObjects[sKey] = result;

			return result;
		}
	}

protected:
	
	string _generate_key_prefix(cl_device_id device, string options)
	{
		char bufName[1000];
		CheckCL( clGetDeviceInfo(device, CL_DEVICE_NAME, 1000, bufName, NULL) );
		return string("DEVICE_") +  string(bufName) /*+string("_OPTIONS_") + options */ + string("_");
	}
	
	string _generate_key(string sSource, cl_device_id device, string options)
	{
		return string("_OPTIONS_") + options+ "\n"+ _generate_key_prefix(device, options) + string("\n") + sSource;
	}
	
	string _generate_key_filename(string sFileName, cl_device_id device, string options)
	{
		return sFileName + _generate_key_prefix(device, options) + ".key";
	}
	
	string _generate_binaries_filename(string sFileName, cl_device_id device, string options)
	{
		return sFileName + _generate_key_prefix(device, options) + ".binaries";
	}
	
	void _store_precompiled_binaries(string sFileName, string sSource, cl_program program, string options)
	{
		map<cl_device_id, vector<unsigned char> > device2binaries;
		
		//EXTRACT BINARIES
		{
			cl_uint nof_devices;
			CheckCL(clGetProgramInfo(program, CL_PROGRAM_NUM_DEVICES, sizeof(nof_devices), &nof_devices, NULL));
			
			size_t sizes[nof_devices];
			CheckCL(clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(sizes), &sizes, NULL));
			
			cl_device_id devices[nof_devices];
			CheckCL(clGetProgramInfo(program, CL_PROGRAM_DEVICES, sizeof(devices), &devices, NULL));
			
			vector< unsigned char * > collection_of_binaries(nof_devices);
			for(unsigned int i = 0; i < nof_devices; i++)
				collection_of_binaries[i] = new unsigned char[sizes[i]];

			CheckCL(clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *)*collection_of_binaries.size(), &collection_of_binaries.front(), NULL));
			
			for(unsigned int i = 0; i < nof_devices; i++)
			{
				device2binaries[devices[i]] = vector<unsigned char>(collection_of_binaries[i], collection_of_binaries[i]+sizes[i]);
				
				delete [] collection_of_binaries[i];
			}
		}
		
		//STORE KEYS, BINARIES
		{
			for(map<cl_device_id, vector<unsigned char> >::iterator it=device2binaries.begin(); it!= device2binaries.end(); it++)
			{
				//key
				{
					string sKeyFileName = _generate_key_filename(sFileName, it->first, options);
					ofstream outKey(sKeyFileName.c_str(), ios::out);
					//cout << "Saving key file "<< sKeyFileName << "..." <<endl;
					
					assert(outKey.is_open());
					string sKey = _generate_key(sSource, it->first, options);
					outKey << sKey;
					outKey.close();
				}
				
				//binaries
				{
					string sBinaryFileName = _generate_binaries_filename(sFileName, it->first, options); 
					ofstream outBinaries(sBinaryFileName.c_str(), ios::out | ios::binary);
					//cout << "Saving binary file "<< sBinaryFileName << "..." <<endl;
					assert(outBinaries.is_open());
					outBinaries.write((const char*)&it->second.front(), it->second.size());
					outBinaries.close();
				}
			}
		}
	}
	
	unsigned char* _load_precompiled_binaries(string sFileName, string sSource, cl_device_id device, int &length, string options)
	{
		length = 0;
		string sKeyFileName = _generate_key_filename(sFileName, device, options);
		//cout << "Looking for "<< sKeyFileName << "..."<< endl;
		
		ifstream keyFile (sKeyFileName.c_str());
		
		if (!keyFile.is_open())
		{
			//cout << "well, i did not find it!"<<endl;
			
			return NULL;
		}
		else 
		{
			//cout << "yes, found it!"<<endl;
			keyFile.close();
			
			//cout << "Verifying that the content of "<< sKeyFileName << "is up-to-date..."<< endl;
			string key = _loadTextFile(sKeyFileName);
			string thiskey =  _generate_key(sSource, device, options);
			
			if (key != thiskey)
			{
				//cout << "no, the content was old !"<<endl;
			//	printf("KEY WAS ***%s**\n", key.data());
			//	printf("THIS KEY IS ***%s**\n", thiskey.data());
				return NULL;
			}
			else 
			{
				//cout << "yes, the content was up-to-date !"<<endl;
				
				string sBinariesFileName = _generate_binaries_filename(sFileName, device, options);
				ifstream binariesFile (sBinariesFileName.c_str(), ios::in | ios::binary);
				
				assert(binariesFile.is_open());
				
				//cout << "loading binaries..."<<endl;
				std::vector<char> contents((std::istreambuf_iterator<char>(binariesFile)),
										   std::istreambuf_iterator<char>());
				
				unsigned char * binaries = new unsigned char[contents.size()];
				assert(binaries != NULL);
				memcpy(binaries, &contents.front(), contents.size());
				//cout << "loading binaries is done!"<<endl;
				
				length = contents.size();
				
				return binaries;
			}
		}		
	}
	
	struct AllBinaries
	{
		unsigned char ** binaries;
		size_t * lengths;
		int nof_items;
		
		AllBinaries(size_t * l, unsigned char ** b, int n): lengths(l), binaries(b), nof_items(n) {}
		
		~AllBinaries()
		{
			for(int j=0; j<nof_items; j++)
				delete [] binaries[j];
			
			delete [] binaries;
			delete [] lengths;
		}
	};
	
	AllBinaries* _load_precompiled_binaries(string sFileName, string sSource, vector<cl_device_id > devices, string options)
	{		
		vector<size_t> lengths;
		unsigned char ** binaries = new unsigned char*[devices.size()];
		for(unsigned int i = 0; i < devices.size(); i++)
		{
			int length = 0;
			binaries[i] = _load_precompiled_binaries(sFileName, sSource, devices[i], length, options);
			lengths.push_back(length);
			
			if (binaries[i] == NULL)
			{
				for(unsigned int j=0; j<=i; j++)
					delete [] binaries[j];
				
				delete [] binaries;
				
				return NULL;
			}
		}
		
		size_t * l = new size_t[devices.size()];
		copy(lengths.begin(), lengths.end(), l);
		
		return new AllBinaries(l, binaries, devices.size());
	}
	
public:
	
	cl_program getProgram(char* sFileName, cl_context context, vector<cl_device_id > devices, string sOptions = "")
	{
		cl_program program;

		char sOptions_[2000];
		strcpy(sOptions_,sOptions.c_str());
		char sFileNameDetailed[2000];
		sprintf(sFileNameDetailed, "%s%s", sFileName, sOptions_);

		if (_getRes(m_progObjects, sFileNameDetailed, program))
			return program;
		else
		{
			string source = _loadTextFile(sFileName).data();
			AllBinaries* allbinaries = _load_precompiled_binaries(sFileNameDetailed, source, devices, sOptions);
			
			if (allbinaries != NULL)
			{
				cl_int binary_statuses[devices.size()];
				cl_int errcode_rets[devices.size()];
				program = clCreateProgramWithBinary(context, devices.size(), &devices.front(), allbinaries->lengths, (const unsigned char**)allbinaries->binaries, binary_statuses, errcode_rets);
				
				/*for(int i=0; i<devices.size(); i++)
					CheckCL(errcode_rets[i]);*/
				
				for(unsigned int i = 0; i < devices.size(); i++)
					CheckCL(binary_statuses[i]);
				
				//cout << "Building ..." ;
				CheckCL(clBuildProgram(program, devices.size(), &devices.front(), sOptions.data(), NULL, NULL));
				//cout << "done!" <<endl ;
				
				delete allbinaries;
				
				m_progObjects[sFileNameDetailed] = program;
				
				return program;
			}
			
			const char * source_chr = source.data();
			size_t length = (size_t)source.size();

			cl_int success;
			program = clCreateProgramWithSource(context, 1, &source_chr , &length , &success);
			CheckCL(success);

			CheckRelaxCL(success = clBuildProgram(program, devices.size(), &devices.front(), sOptions.data(), NULL, NULL));

			if(success != CL_SUCCESS)
			{
				cout<<"Press any button."<<endl;
				char c;
				cin>>c;

				for(unsigned int i = 0; i < devices.size(); i++)
				{
					printf("device %d\n", i);

					bool bEnoughSpace = false;
					size_t n=50000, a=0;
					char * buf = NULL;
					do
					{
						n *= 2;
						if (buf != NULL) delete [] buf;
						buf  = new char[n];
						CheckCL( clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_LOG, n, buf, &a) );

						bEnoughSpace = n>a;
					}
					while(!bEnoughSpace);

					cout << "compilation errors: "<< string(buf) << endl;
				//	printf("compilation errors: %s\n(%d chars)\n", buf, nactual);

					delete buf;

					cout<<"Press any button."<<endl;
					cin>>c;
				}
			}
			else 
				_store_precompiled_binaries(sFileNameDetailed, source, program, sOptions);


			m_progObjects[sFileNameDetailed] = program;

			return program;
		}
	}


	cl_kernel getKernel(string sKey, string name, cl_program program)
	{
		cl_kernel kernel;
		
		map<string, cl_kernel>& mymap = m_kernelObjects[program];

		if (_getRes(mymap, name, kernel))
			return kernel;
		else
		{
			cl_int success;
			cout << sKey << endl;
			kernel = clCreateKernel(program, sKey.data(), &success);
			CheckCL(success);

			m_kernelObjects[program][name] = kernel;

			return kernel;
		}
	}

	~ResourcesCL()
	{
		map< string, vector<cl_mem> >::iterator it1;
		if(m_bufObjects.size() != 0)
			for(it1 = m_bufObjects.begin(); it1 != m_bufObjects.end(); it1++)
				for(unsigned int i = 0; i < it1->second.size(); i++)
					clReleaseMemObject(it1->second[i]);

		map< string, vector<cl_mem> >::iterator it2;
		if(m_imgObjects.size() != 0)
			for(it2 = m_imgObjects.begin(); it2 != m_imgObjects.end(); it2++)
				for(unsigned int i = 0; i < it2->second.size(); i++)
					clReleaseMemObject(it2->second[i]);

		map< string, cl_program >::iterator itp;
		if(m_progObjects.size() != 0)
			for(itp = m_progObjects.begin(); itp != m_progObjects.end(); itp++)
				clReleaseProgram(itp->second);

		map< cl_program,  map<string, cl_kernel> >::iterator itPP;
		map< string, cl_kernel >::iterator itk;
		if(m_kernelObjects.size() != 0)
			for(itPP = m_kernelObjects.begin(); itPP != m_kernelObjects.end(); itPP++)
			{
				for(itk = itPP->second.begin(); itk != itPP->second.end(); itk++)
					clReleaseKernel(itk->second);
			}

		map< string, cl_sampler >::iterator its;
		if(m_samplerObjects.size() != 0)
			for(its = m_samplerObjects.begin(); its != m_samplerObjects.end(); its++)
				clReleaseSampler(its->second);
	}
};
