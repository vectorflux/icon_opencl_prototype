/*
 * ProfilerCL.h
 *
 *  Created on: Mar 5, 2010
 *      Author: claurent
 */

#ifndef PROFILERCL_H_
#define PROFILERCL_H_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <math.h>
#include <map>

class ProfileCLAgent
{
protected:
	cl_event* event;
	double time;
	int samples;

public:
	ProfileCLAgent() {};

	ProfileCLAgent(cl_event* event_) : event(event_), time(0.0), samples(0) {};

	void timer()
	{
		cl_ulong start_time;
		cl_ulong stop_time;
		CheckCL( clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start_time, NULL) );
		CheckCL( clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &stop_time, NULL) );
		time += (double)(stop_time - start_time) * pow(10.0,-9);
		samples++;
	}

	double getTime()
	{
		return time;
	}

	void printSummary()
	{
		printf("\nAvgTime:\t%e [s]", time / samples);
	}
};

class ProfilerCL
{
private:
	map<string, ProfileCLAgent> m_mapAgents;

public:
	
	void addTask(string name, cl_event* event)
	{
		ProfileCLAgent agent(event);
		m_mapAgents[name] = agent;
	}

	void getTimes()
	{
		map<string, ProfileCLAgent>::iterator it;
		if(m_mapAgents.size() != 0)
			for(it = m_mapAgents.begin(); it != m_mapAgents.end(); it++)
			{
				//printf("%s\n",it->first.c_str());
				it->second.timer();
			}
	}

	double getTime(string name)
	{
		assert(!m_mapAgents.empty());
		map<string, ProfileCLAgent>::iterator it = m_mapAgents.find(name.c_str());
		assert(it != m_mapAgents.end());
		return it->second.getTime();
	}

	void print()
	{
		//printf("computing total time\n");
		map<string, ProfileCLAgent>::iterator it;
		double totalTime = 0;
		if(m_mapAgents.size() != 0)
			for(it = m_mapAgents.begin(); it != m_mapAgents.end(); it++)
				totalTime += it->second.getTime();
		
		printf("\nTIME ALLOCATION:\n");

		for(it = m_mapAgents.begin(); it != m_mapAgents.end(); it++)
			printf("'%30s':\t%3.3e s\t%2.2f%%\n", it->first.c_str(), it->second.getTime(), it->second.getTime()/totalTime*100);

		printf("END\n");
	}

	void cleanup()
	{
		m_mapAgents.clear();
	}
};

#endif /* PROFILERCL_H_ */
