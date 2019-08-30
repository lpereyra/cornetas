#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
#include <time.h>

/*
 use the macro OMPP_TIMESTAMP with a *double* 
 parameter argument. the returned value is time
 in seconds passed since some point of time 
 in the past
*/
#define TIMER( time_ )  	  			       \
{                                        \
  struct timeval tv;                     \
  gettimeofday( &tv, NULL );             \
  time_=tv.tv_sec+tv.tv_usec*1.0e-6;     \
}

#endif /* TIMER_H */
