/*
 * helper.h
 *
 *  Created on: Oct 1, 2014
 *      Author: yqian33
 */

#ifndef HELPER_H_
#define HELPER_H_

#include <sys/time.h> 	//sys/time.h是Linux系统的日期时间头文件，sys/time.h通常会包含include time.h。通常如果代码可以是平台相关的，则只需要include sys/time.h。
#include <time.h>		//编写的代码如果是平台无关的，则需要在代码里include time.h，
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <hash_map>

using namespace std;

void outUsedTime(int flag);

#endif /* HELPER_H_ */
