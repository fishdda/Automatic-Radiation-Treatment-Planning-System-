/*
 * Util.cpp
 *
 *  Created on: Feb 2, 2012
 *      Author: user
 */

#include "Util.h"

/** Return time expressed as a free-running microsecond clock
 *
 *  Uses the gettimeofday system call, but converts result to
 *  simple microsecond clock for greater convenience.
 */

uint32_t Util::getTime() {
        // note use of static variables
        static uint32_t now;
        static struct timeval prevTimeval = { 0, 0 };

        if (prevTimeval.tv_sec == 0 && prevTimeval.tv_usec == 0) {
                // first call to getTime(); initialize and return 0
                if (gettimeofday(&prevTimeval, NULL) < 0)
                        printf("Util::getTime: gettimeofday failure");
                now = 0;
                return 0;
        }
        // normal case
        struct timeval nowTimeval;
        if (gettimeofday(&nowTimeval, NULL) < 0)
                printf("Util::getTime: gettimeofday failure");
	uint32_t dsec = nowTimeval.tv_sec; dsec -= prevTimeval.tv_sec;
	uint32_t dusec = nowTimeval.tv_usec - prevTimeval.tv_usec;
	if (nowTimeval.tv_usec < prevTimeval.tv_usec) {
		dusec = nowTimeval.tv_usec + (1000000 - prevTimeval.tv_usec);
		dsec--;
	}
	now += 1000000*dsec + dusec;
        prevTimeval = nowTimeval;

        return now;
}
