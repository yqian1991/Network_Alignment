/*
 * helper.cpp
 *
 *  Created on: Oct 1, 2014
 *      Author: yqian33
 */

#include "helper.h"

#define CLOCKS_PER_SEC 1000

void outUsedTime(int flag)
{
	static struct timeval tp_start, tp_end;//struct defined in time.h, used to save the accurate system time
	double time_used;
	static unsigned long start_cpu, end_cpu; // CPU time

	if (flag == 0) {
		gettimeofday(&tp_start, NULL);
		//set start time: CPU time
		start_cpu = clock();
	} else {
		//set over time: CPU time
		end_cpu = clock();
		//print the system time
		gettimeofday(&tp_end, NULL);
		//print the CPU time
		cout<<"Total CPU time used: "<< (double) (end_cpu - start_cpu) / (CLOCKS_PER_SEC) << "seconds." << endl;
		time_used = (1000000 * (tp_end.tv_sec - tp_start.tv_sec)
				+ (double) (tp_end.tv_usec - tp_start.tv_usec)) / 1000000;
		cout<<"Total Used Time By gettimeofday(): "<< time_used <<" Seconds."<<endl;
		gettimeofday(&tp_start, NULL);
		//set start time: CPU time
		start_cpu = clock();
	}
	return;
}
