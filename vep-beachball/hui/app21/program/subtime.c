/**************************************************************************
* FILE: subtime.c
* Author: Huai Zhang
* DESCRIPTION:
*  Called by fortran version of timer subroutines to obtain
*  seconds and microseconds from gettimeofday, and system and user
*  seconds from getrusage, and parallel user time from MCLOCK().
*  
***************************************************************************/
/**************************************************************************/
/* void function gettimeofday to get wall time of present time */
/* times[0]=elapsed time(integer part), times[1]=elapsed time(extra part) */
/**************************************************************************/
#include <sys/resource.h> 
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>

void fwtiming_(int times[2]) 
{
  struct timeval tv;

  gettimeofday(&tv, (struct timezone*)0);	  
  times[0] = tv.tv_sec;
  times[1] = tv.tv_usec;
  return;
}
/**********************************************************************/


/**********************************************************************/
/* void function getrusage to get user time and system consuming      */
/* tarray[0]=User time, tarray[1]=Sys time,  all expressed in seconds */
/**********************************************************************/
void futiming_(double *tarray) 
{
  struct rusage RU;

  getrusage(RUSAGE_SELF, &RU);

  tarray[0]=RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6;
  tarray[1]=RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
  return;
} 


/********************************************************************/
/* This function returns various types of elapsed time              */
/* tarray[0]=User time, tarray[1]=Sys time                          */
/* tarray[2]=Greenwich time                                         */
/* all expressed in seconds                                         */
/**********************************************************************/

#include <sys/resource.h>
#include <sys/time.h>

void gettime_(double *tarray)
{
  struct rusage RU;
  struct timeval tv;

  getrusage(RUSAGE_SELF, &RU);
  gettimeofday(&tv, (struct timezone *)0);

  tarray[0]=RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6;
  tarray[1]=RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
  tarray[2]=tv.tv_sec+(double)tv.tv_usec*1e-6;
  return;
}



/**********************************************************************/
/* this function return the absolute path of current runtime program  */
/**********************************************************************/
#include <unistd.h>
  void fgetcwd_(int *lenth, char *abspath) 
{
        size_t size;
        size = 80;

        getcwd(abspath, size);
        *lenth = strlen(abspath);
/*        printf("current path lenth, max lenth is %d %d \n", *lenth, size);*/
/*        printf("current path is %-1.80s \n", abspath);                    */
        return;
} 
/****************************************************************************/



/****************************************************************************/
/* Functions list belowing are used in IBM RS6000 system                    */
/****************************************************************************/
/* void function gettimeofday to get wall time of present time              */
/* times[0]=elapsed time(integer part), times[1]=elapsed time(extra part)   */
/****************************************************************************/
#include <sys/resource.h> 
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>

void fwtiming(int times[2]) 
{
  struct timeval tv;

  gettimeofday(&tv, (struct timezone*)0);	  
  times[0] = tv.tv_sec;
  times[1] = tv.tv_usec;
  return;
}
/****************************************************************************/


/**********************************************************************/
/* void function getrusage to get user time and system consuming      */
/* tarray[0]=User time, tarray[1]=Sys time,  all expressed in seconds */
/**********************************************************************/
void futiming(double *tarray) 
{
  struct rusage RU;

  getrusage(RUSAGE_SELF, &RU);

  tarray[0]=RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6;
  tarray[1]=RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
  return;
} 
/**********************************************************************/


/********************************************************************/
/* This function returns various types of elapsed time              */
/* tarray[0]=User time, tarray[1]=Sys time                          */
/* tarray[2]=Greenwich time                                         */
/* all expressed in seconds                                         */
/********************************************************************/

#include <sys/resource.h>
#include <sys/time.h>

void gettime(double *tarray)
{
  struct rusage RU;
  struct timeval tv;

  getrusage(RUSAGE_SELF, &RU);
  gettimeofday(&tv, (struct timezone *)0);

  tarray[0]=RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6;
  tarray[1]=RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
  tarray[2]=tv.tv_sec+(double)tv.tv_usec*1e-6;
  return;
}
/********************************************************************/



/************************************************************************/
/* this function return the absolute path of current runtime program    */
/************************************************************************/
#include <unistd.h>
  void fgetcwd(int *lenth, char *abspath) 
{
        size_t size;
        size = 80;

        getcwd(abspath, size);
        *lenth = strlen(abspath);
/*    printf("current path lenth, max lenth is %d %d \n", *lenth, size);*/
    printf("Your current runtime path is %-1.80s \n", abspath);
        return;
}
/************************************************************************/


/*************************************************************************/
/* another method to get system and user time, current is not in use     */
/*************************************************************************/
void fctiming_(char times) {
        time_t  ltime;

        ltime =  time(NULL);
/*        times = *ctime(&ltime);                                        */
        times = ctime(&ltime);
        return;
}


/*************************************************************************/
/* This function returns the seconds, currently is not in use            */
/*************************************************************************/
#include <time.h>
void getseconds_ (double *seconds)
{   clock_t clock(void);
	*seconds = clock()/CLOCKS_PER_SEC;
	printf("current time is, %d",*seconds);
	return;	
}
/*************************************************************************/
