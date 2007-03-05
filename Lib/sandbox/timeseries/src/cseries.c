#include <Python.h>
#include <datetime.h>
#include <structmember.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "arrayobject.h"

static char cseries_doc[] = "Speed sensitive time series operations";

#define FR_ANN  1000  /* Annual */
#define FR_QTR  2000  /* Quarterly */
#define FR_MTH  3000  /* Monthly */

#define FR_WK   4000  /* Weekly */
#define FR_WKSUN FR_WK /* Weekly - Sunday end of week */
#define FR_WKMON 4001 /* Weekly - Monday end of week */
#define FR_WKTUE 4002 /* Weekly - Tuesday end of week */
#define FR_WKWED 4003 /* Weekly - Wednesday end of week */
#define FR_WKTHU 4004 /* Weekly - Thursday end of week */
#define FR_WKFRI 4005 /* Weekly - Friday end of week */
#define FR_WKSAT 4006 /* Weekly - Saturday end of week */

#define FR_BUS  5000  /* Business days */
#define FR_DAY  6000  /* Daily */
#define FR_HR   7000  /* Hourly */
#define FR_MIN  8000  /* Minutely */
#define FR_SEC  9000  /* Secondly */
#define FR_UND  -10000 /* Undefined */

#define ADD_INT_TO_DICT(dict, key, val) \
    {PyObject *pyval = PyInt_FromLong(val); \
     PyDict_SetItemString(dict, key, pyval); \
     Py_DECREF(pyval); }

#define DINFO_ERR -99

#define CHECK_ASFREQ(result) if ((result) == DINFO_ERR) return NULL

struct asfreq_info{
    int to_week_end; //day the week ends on in the "to" frequency
    int from_week_end; //day the week ends on in the "from" frequency
};

static struct asfreq_info NULL_AF_INFO;

//DERIVED FROM mx.DateTime
/*
=====================================================
== Functions in the following section are borrowed ==
== from mx.DateTime, and in many cases slightly    ==
== modified                                        ==
=====================================================
*/

#define Py_AssertWithArg(x,errortype,errorstr,a1) {if (!(x)) {PyErr_Format(errortype,errorstr,a1);goto onError;}}
#define Py_Error(errortype,errorstr) {PyErr_SetString(errortype,errorstr);goto onError;}

static PyObject *DateCalc_Error; /* Error Exception object */
static PyObject *DateCalc_RangeError; /* Error Exception object */



#define GREGORIAN_CALENDAR 0
#define JULIAN_CALENDAR 1

#define SECONDS_PER_DAY ((double) 86400.0)

/* Table with day offsets for each month (0-based, without and with leap) */
static int month_offset[2][13] = {
    { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 },
    { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 }
};

/* Table of number of days in a month (0-based, without and with leap) */
static int days_in_month[2][12] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
};

struct date_info {
    long absdate;
    double abstime;

    double second;
    int minute;
    int hour;
    int day;
    int month;
    int quarter;
    int year;
    int day_of_week;
    int day_of_year;
    int calendar;
};


/* Return 1/0 iff year points to a leap year in calendar. */
static
int dInfoCalc_Leapyear(register long year,
            int calendar)
{
    if (calendar == GREGORIAN_CALENDAR) {
        return (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0));
    } else {
        return (year % 4 == 0);
    }
}

static
int dInfoCalc_ISOWeek(struct date_info dinfo)
{
    int week;

    /* Estimate */
    week = (dinfo.day_of_year-1) - dinfo.day_of_week + 3;
    if (week >= 0) week = week / 7 + 1;

    /* Verify */
    if (week < 0) {
        /* The day lies in last week of the previous year */
        if ((week > -2) ||
            (week == -2 && dInfoCalc_Leapyear(dinfo.year-1, dinfo.calendar)))
            week = 53;
        else
            week = 52;
    } else if (week == 53) {
    /* Check if the week belongs to year or year+1 */
        if (31-dinfo.day + dinfo.day_of_week < 3) {
            week = 1;
        }
    }

    return week;
}


/* Return the day of the week for the given absolute date. */
static
int dInfoCalc_DayOfWeek(register long absdate)
{
    int day_of_week;

    if (absdate >= 1) {
        day_of_week = (absdate - 1) % 7;
    } else {
        day_of_week = 6 - ((-absdate) % 7);
    }
    return day_of_week;
}

/* Return the year offset, that is the absolute date of the day
   31.12.(year-1) in the given calendar.

   Note:
   For the Julian calendar we shift the absdate (which is measured
   using the Gregorian Epoch) value by two days because the Epoch
   (0001-01-01) in the Julian calendar lies 2 days before the Epoch in
   the Gregorian calendar. */
static
int dInfoCalc_YearOffset(register long year,
              int calendar)
{
    year--;
    if (calendar == GREGORIAN_CALENDAR) {
    if (year >= 0 || -1/4 == -1)
        return year*365 + year/4 - year/100 + year/400;
    else
        return year*365 + (year-3)/4 - (year-99)/100 + (year-399)/400;
    }
    else if (calendar == JULIAN_CALENDAR) {
    if (year >= 0 || -1/4 == -1)
        return year*365 + year/4 - 2;
    else
        return year*365 + (year-3)/4 - 2;
    }
    Py_Error(DateCalc_Error, "unknown calendar");
 onError:
    return -1;
}


/* Set the instance's value using the given date and time. calendar
   may be set to the flags: GREGORIAN_CALENDAR,
   JULIAN_CALENDAR to indicate the calendar to be used. */

static
int dInfoCalc_SetFromDateAndTime(struct date_info *dinfo,
                  int year,
                  int month,
                  int day,
                  int hour,
                  int minute,
                  double second,
                  int calendar)
{

    /* Calculate the absolute date */
    {
        int leap;
        long yearoffset,absdate;

        /* Range check */
        Py_AssertWithArg(year > -(INT_MAX / 366) && year < (INT_MAX / 366),
                 DateCalc_RangeError,
                 "year out of range: %i",
                 year);

        /* Is it a leap year ? */
        leap = dInfoCalc_Leapyear(year,calendar);

        /* Negative month values indicate months relative to the years end */
        if (month < 0) month += 13;
        Py_AssertWithArg(month >= 1 && month <= 12,
                 DateCalc_RangeError,
                 "month out of range (1-12): %i",
                 month);

        /* Negative values indicate days relative to the months end */
        if (day < 0) day += days_in_month[leap][month - 1] + 1;
        Py_AssertWithArg(day >= 1 && day <= days_in_month[leap][month - 1],
                 DateCalc_RangeError,
                 "day out of range: %i",
                 day);

        yearoffset = dInfoCalc_YearOffset(year,calendar);
        if (yearoffset == -1 && PyErr_Occurred()) goto onError;

        absdate = day + month_offset[leap][month - 1] + yearoffset;

        dinfo->absdate = absdate;

        dinfo->year = year;
        dinfo->month = month;
        dinfo->quarter = ((month-1)/3)+1;
        dinfo->day = day;

        dinfo->day_of_week = dInfoCalc_DayOfWeek(absdate);
        dinfo->day_of_year = (short)(absdate - yearoffset);

        dinfo->calendar = calendar;
    }

    /* Calculate the absolute time */
    {
    Py_AssertWithArg(hour >= 0 && hour <= 23,
             DateCalc_RangeError,
             "hour out of range (0-23): %i",
             hour);
    Py_AssertWithArg(minute >= 0 && minute <= 59,
             DateCalc_RangeError,
             "minute out of range (0-59): %i",
             minute);
    Py_AssertWithArg(second >= (double)0.0 &&
             (second < (double)60.0 ||
              (hour == 23 && minute == 59 &&
               second < (double)61.0)),
             DateCalc_RangeError,
             "second out of range (0.0 - <60.0; <61.0 for 23:59): %f",
             second);

    dinfo->abstime = (double)(hour*3600 + minute*60) + second;

    dinfo->hour = hour;
    dinfo->minute = minute;
    dinfo->second = second;
    }
    return 0;
 onError:
    return -1;
}


/* Sets the date part of the date_info struct using the indicated
   calendar.

   XXX This could also be done using some integer arithmetics rather
       than with this iterative approach... */
static
int dInfoCalc_SetFromAbsDate(register struct date_info *dinfo,
                  long absdate,
                  int calendar)
{
    register long year;
    long yearoffset;
    int leap,dayoffset;
    int *monthoffset;

    /* Approximate year */
    if (calendar == GREGORIAN_CALENDAR) {
        year = (long)(((double)absdate) / 365.2425);
    } else if (calendar == JULIAN_CALENDAR) {
        year = (long)(((double)absdate) / 365.25);
    } else {
        Py_Error(DateCalc_Error, "unknown calendar");
    }
    if (absdate > 0) year++;

    /* Apply corrections to reach the correct year */
    while (1) {
        /* Calculate the year offset */
        yearoffset = dInfoCalc_YearOffset(year,calendar);
        if (yearoffset == -1 && PyErr_Occurred())
            goto onError;

        /* Backward correction: absdate must be greater than the
           yearoffset */
        if (yearoffset >= absdate) {
            year--;
            continue;
        }

        dayoffset = absdate - yearoffset;
        leap = dInfoCalc_Leapyear(year,calendar);

        /* Forward correction: non leap years only have 365 days */
        if (dayoffset > 365 && !leap) {
            year++;
            continue;
        }
        break;
    }

    dinfo->year = year;
    dinfo->calendar = calendar;

    /* Now iterate to find the month */
    monthoffset = month_offset[leap];
    {
        register int month;

        for (month = 1; month < 13; month++) {
            if (monthoffset[month] >= dayoffset)
            break;
        }

        dinfo->month = month;
        dinfo->quarter = ((month-1)/3)+1;
        dinfo->day = dayoffset - month_offset[leap][month-1];
    }


    dinfo->day_of_week = dInfoCalc_DayOfWeek(absdate);
    dinfo->day_of_year = dayoffset;
    dinfo->absdate = absdate;

    return 0;

 onError:
    return -1;
}

/* Sets the time part of the DateTime object. */
static
int dInfoCalc_SetFromAbsTime(struct date_info *dinfo,
                  double abstime)
{
    int inttime;
    int hour,minute;
    double second;

    inttime = (int)abstime;
    hour = inttime / 3600;
    minute = (inttime % 3600) / 60;
    second = abstime - (double)(hour*3600 + minute*60);

    dinfo->hour = hour;
    dinfo->minute = minute;
    dinfo->second = second;

    dinfo->abstime = abstime;

    return 0;
}

/* Set the instance's value using the given date and time. calendar
   may be set to the flags: GREGORIAN_CALENDAR, JULIAN_CALENDAR to
   indicate the calendar to be used. */
static
int dInfoCalc_SetFromAbsDateTime(struct date_info *dinfo,
                  long absdate,
                  double abstime,
                  int calendar)
{

    /* Bounds check */
    Py_AssertWithArg(abstime >= 0.0 && abstime <= SECONDS_PER_DAY,
             DateCalc_Error,
             "abstime out of range (0.0 - 86400.0): %f",
             abstime);

    /* Calculate the date */
    if (dInfoCalc_SetFromAbsDate(dinfo,
                  absdate,
                  calendar))
    goto onError;

    /* Calculate the time */
    if (dInfoCalc_SetFromAbsTime(dinfo,
                  abstime))
    goto onError;

    return 0;
 onError:
    return -1;
}

/*
====================================================
== End of section borrowed from mx.DateTime       ==
====================================================
*/


//////////////////////////////////////////////////////////

static long minval_D_toHighFreq = 719163;

///////////////////////////////////////////////////////////////////////

static long absdatetime_hour(long absdate, long time) {

}

///////////////////////////////////////////////////////////////////////

// helpers for frequency conversion routines //

static long DtoB_weekday(long fromDate) { return (((fromDate) / 7) * 5) + (fromDate)%7; }

static long DtoB_WeekendToMonday(struct date_info dinfo) {

    long absdate = dinfo.absdate;
    if (dinfo.day_of_week > 4) {
        //change to Monday after weekend
        absdate += (7 - dinfo.day_of_week);
    }
    return DtoB_weekday(absdate);
}

static long DtoB_WeekendToFriday(struct date_info dinfo) {

    long absdate = dinfo.absdate;
    if (dinfo.day_of_week > 4) {
        //change to friday before weekend
        absdate -= (dinfo.day_of_week - 4);
    }
    return DtoB_weekday(absdate);
}

static long absdate_from_ymd(int y, int m, int d) {
    struct date_info tempDate;
    if (dInfoCalc_SetFromDateAndTime(&tempDate, y, m, d, 0, 0, 0, GREGORIAN_CALENDAR)) return DINFO_ERR;
    return tempDate.absdate;
}


///////////////////////////////////////////////

// frequency specifc conversion routines
// each function must take an integer fromDate and a char relation ('B' or 'A' for 'BEFORE' or 'AFTER')

//************ FROM DAILY ***************

static long asfreq_DtoA(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return DINFO_ERR;
    return (long)(dinfo.year);
}

static long asfreq_DtoQ(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return DINFO_ERR;
    return (long)((dinfo.year - 1) * 4 + dinfo.quarter);
}

static long asfreq_DtoM(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return DINFO_ERR;
    return (long)((dinfo.year - 1) * 12 + dinfo.month);
}

static long asfreq_DtoW(long fromDate, char relation, struct asfreq_info af_info) {
    return (fromDate - (1 + af_info.to_week_end))/7 + 1;
}

static long asfreq_DtoB(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return DINFO_ERR;

    if (relation == 'B') {
        return DtoB_WeekendToFriday(dinfo);
    } else {
        return DtoB_WeekendToMonday(dinfo);
    }
}

static long asfreq_DtoB_forConvert(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return DINFO_ERR;

    if (dinfo.day_of_week > 4) {
        return -1;
    } else {
        return DtoB_weekday(fromDate);
    }
}

// needed for getDateInfo function
static long asfreq_DtoD(long fromDate, char relation, struct asfreq_info af_info) { return fromDate; }

static long asfreq_DtoHIGHFREQ(long fromDate, char relation, long periodsPerDay) {
    if (fromDate >= minval_D_toHighFreq) {
        if (relation == 'B') { return (fromDate - minval_D_toHighFreq)*(periodsPerDay) + 1; }
        else                 { return (fromDate - minval_D_toHighFreq + 1)*(periodsPerDay); }
    } else { return -1; }
}

static long asfreq_DtoH(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24); }
static long asfreq_DtoT(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24*60); }
static long asfreq_DtoS(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24*60*60); }

//************ FROM SECONDLY ***************

static long asfreq_StoD(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1)/(60*60*24) + minval_D_toHighFreq; }

static long asfreq_StoA(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoA(asfreq_StoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_StoQ(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoQ(asfreq_StoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_StoM(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoM(asfreq_StoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_StoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_StoD(fromDate, relation, NULL_AF_INFO), relation, af_info); }
static long asfreq_StoB(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoB(asfreq_StoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_StoB_forConvert(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoB_forConvert(asfreq_StoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_StoT(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1)/60 + 1; }
static long asfreq_StoH(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1)/(60*60) + 1; }

//************ FROM MINUTELY ***************

static long asfreq_TtoD(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1)/(60*24) + minval_D_toHighFreq; }

static long asfreq_TtoA(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoA(asfreq_TtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_TtoQ(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoQ(asfreq_TtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_TtoM(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoM(asfreq_TtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_TtoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_TtoD(fromDate, relation, NULL_AF_INFO), relation, af_info); }
static long asfreq_TtoB(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoB(asfreq_TtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_TtoB_forConvert(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoB_forConvert(asfreq_TtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_TtoH(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1)/60 + 1; }
static long asfreq_TtoS(long fromDate, char relation, struct asfreq_info af_info) {
    if (relation == 'B') {  return fromDate*60 - 59; }
    else                 {  return fromDate*60;      }}

//************ FROM HOURLY ***************

static long asfreq_HtoD(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1)/24 + minval_D_toHighFreq; }
static long asfreq_HtoA(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoA(asfreq_HtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_HtoQ(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoQ(asfreq_HtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_HtoM(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoM(asfreq_HtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_HtoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_HtoD(fromDate, relation, NULL_AF_INFO), relation, af_info); }
static long asfreq_HtoB(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoB(asfreq_HtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_HtoB_forConvert(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoB_forConvert(asfreq_HtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

// calculation works out the same as TtoS, so we just call that function for HtoT
static long asfreq_HtoT(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_TtoS(fromDate, relation, NULL_AF_INFO); }
static long asfreq_HtoS(long fromDate, char relation, struct asfreq_info af_info) {
    if (relation == 'B') {  return fromDate*60*60 - 60*60 + 1; }
    else                 {  return fromDate*60*60;             }}

//************ FROM BUSINESS ***************

static long asfreq_BtoD(long fromDate, char relation, struct asfreq_info af_info)
    { return ((fromDate-1)/5)*7 + (fromDate-1)%5 + 1; }

static long asfreq_BtoA(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoA(asfreq_BtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_BtoQ(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoQ(asfreq_BtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_BtoM(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoM(asfreq_BtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_BtoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_BtoD(fromDate, relation, NULL_AF_INFO), relation, af_info); }

static long asfreq_BtoH(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoH(asfreq_BtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_BtoT(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoT(asfreq_BtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

static long asfreq_BtoS(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoS(asfreq_BtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

//************ FROM WEEKLY ***************

static long asfreq_WtoD(long fromDate, char relation, struct asfreq_info af_info) {
    if (relation == 'B') { return fromDate * 7 - 6 + af_info.from_week_end;}
    else                 { return fromDate * 7 + af_info.from_week_end; }}

static long asfreq_WtoA(long fromDate, char relation, struct asfreq_info af_info) {
    return asfreq_DtoA(asfreq_WtoD(fromDate, 'A', af_info), relation, NULL_AF_INFO); }
static long asfreq_WtoQ(long fromDate, char relation, struct asfreq_info af_info) {
    return asfreq_DtoQ(asfreq_WtoD(fromDate, 'A', af_info), relation, NULL_AF_INFO); }
static long asfreq_WtoM(long fromDate, char relation, struct asfreq_info af_info) {
    return asfreq_DtoM(asfreq_WtoD(fromDate, 'A', af_info), relation, NULL_AF_INFO); }

static long asfreq_WtoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_WtoD(fromDate, relation, af_info), relation, NULL_AF_INFO); }

static long asfreq_WtoB(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_WtoD(fromDate, relation, af_info),
                    GREGORIAN_CALENDAR)) return DINFO_ERR;

    if (relation == 'B') { return DtoB_WeekendToMonday(dinfo); }
    else                 { return DtoB_WeekendToFriday(dinfo); }
}

static long asfreq_WtoH(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoH(asfreq_WtoD(fromDate, relation, af_info), relation, NULL_AF_INFO); }
static long asfreq_WtoT(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoT(asfreq_WtoD(fromDate, relation, af_info), relation, NULL_AF_INFO); }
static long asfreq_WtoS(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoS(asfreq_WtoD(fromDate, relation, af_info), relation, NULL_AF_INFO); }

//************ FROM MONTHLY ***************

static void MtoD_ym(long fromDate, long *y, long *m) {
    *y = (fromDate - 1) / 12 + 1;
    *m = fromDate - 12 * (*y) - 1;
}

static long asfreq_MtoD(long fromDate, char relation, struct asfreq_info af_info) {

    long y, m, absdate;

    if (relation == 'B') {
        MtoD_ym(fromDate, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == DINFO_ERR) return DINFO_ERR;
        return absdate;
    } else {
        MtoD_ym(fromDate+1, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == DINFO_ERR) return DINFO_ERR;
        return absdate-1;
    }
}

static long asfreq_MtoA(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1) / 12 + 1; }
static long asfreq_MtoQ(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1) / 3 + 1; }

static long asfreq_MtoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_MtoD(fromDate, relation, NULL_AF_INFO), relation, af_info); }

static long asfreq_MtoB(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_MtoD(fromDate, relation, NULL_AF_INFO),
                    GREGORIAN_CALENDAR)) return DINFO_ERR;

    if (relation == 'B') { return DtoB_WeekendToMonday(dinfo); }
    else                 { return DtoB_WeekendToFriday(dinfo); }
}

static long asfreq_MtoH(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoH(asfreq_MtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_MtoT(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoT(asfreq_MtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_MtoS(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoS(asfreq_MtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }

//************ FROM QUARTERLY ***************

static void QtoD_ym(long fromDate, long *y, long *m) {
    *y = (fromDate - 1) / 4 + 1;
    *m = (fromDate + 4) * 3 - 12 * (*y) - 2;
}

static long asfreq_QtoD(long fromDate, char relation, struct asfreq_info af_info) {

    long y, m, absdate;

    if (relation == 'B') {
        QtoD_ym(fromDate, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == DINFO_ERR) return DINFO_ERR;
        return absdate;
    } else {
        QtoD_ym(fromDate+1, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == DINFO_ERR) return DINFO_ERR;
        return absdate - 1;
    }
}

static long asfreq_QtoA(long fromDate, char relation, struct asfreq_info af_info)
    { return (fromDate - 1)/ 4 + 1; }

static long asfreq_QtoM(long fromDate, char relation, struct asfreq_info af_info) {
    if (relation == 'B') { return fromDate * 3 - 2; }
    else {                 return fromDate  * 3; }
}

static long asfreq_QtoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_QtoD(fromDate, relation, NULL_AF_INFO), relation, af_info); }

static long asfreq_QtoB(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_QtoD(fromDate, relation, NULL_AF_INFO),
                    GREGORIAN_CALENDAR)) return DINFO_ERR;

    if (relation == 'B') { return DtoB_WeekendToMonday(dinfo); }
    else                 { return DtoB_WeekendToFriday(dinfo); }
}


static long asfreq_QtoH(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoH(asfreq_QtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_QtoT(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoT(asfreq_QtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_QtoS(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoS(asfreq_QtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }


//************ FROM ANNUAL ***************

static long asfreq_AtoD(long fromDate, char relation, struct asfreq_info af_info) {
    long absdate;
    if (relation == 'B') {
        if ((absdate = absdate_from_ymd(fromDate,1,1)) == DINFO_ERR) return DINFO_ERR;
        return absdate;
    } else {
        if ((absdate = absdate_from_ymd(fromDate+1,1,1)) == DINFO_ERR) return DINFO_ERR;
        return absdate - 1;
    }
}

static long asfreq_AtoQ(long fromDate, char relation, struct asfreq_info af_info) {
    if (relation == 'B') { return fromDate * 4 - 3; }
    else {                 return fromDate * 4; }
}

static long asfreq_AtoM(long fromDate, char relation, struct asfreq_info af_info) {
    if (relation == 'B') { return fromDate * 12 - 11; }
    else {                 return fromDate * 12; }
}

static long asfreq_AtoW(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoW(asfreq_AtoD(fromDate, relation, NULL_AF_INFO), relation, af_info); }

static long asfreq_AtoB(long fromDate, char relation, struct asfreq_info af_info) {

    struct date_info dailyDate;

    if (relation == 'B') {
        if (dInfoCalc_SetFromDateAndTime(&dailyDate,
                        fromDate,1,1, 0, 0, 0,
                        GREGORIAN_CALENDAR)) return DINFO_ERR;
        return DtoB_WeekendToMonday(dailyDate);
    } else {
        long absdate;

        if (dInfoCalc_SetFromDateAndTime(&dailyDate,
                       fromDate+1,1,1, 0, 0, 0,
                       GREGORIAN_CALENDAR)) return DINFO_ERR;

        absdate = dailyDate.absdate - 1;

        if(dInfoCalc_SetFromAbsDate(&dailyDate,
                      absdate,
                      GREGORIAN_CALENDAR)) return DINFO_ERR;

        return DtoB_WeekendToFriday(dailyDate);
    }
}

static long asfreq_AtoH(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoH(asfreq_AtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_AtoT(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoT(asfreq_AtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }
static long asfreq_AtoS(long fromDate, char relation, struct asfreq_info af_info)
    { return asfreq_DtoS(asfreq_AtoD(fromDate, relation, NULL_AF_INFO), relation, NULL_AF_INFO); }


static long nofunc(long fromDate, char relation, struct asfreq_info af_info) { return -1; }

// end of frequency specific conversion routines

// return a pointer to appropriate conversion function
static long (*get_asfreq_func(int fromFreq, int toFreq, int forConvert))(long, char, struct asfreq_info) {

    int fromGroup = (fromFreq/1000)*1000;
    int toGroup = (toFreq/1000)*1000;

    switch(fromGroup)
    {
        case FR_ANN:
            switch(toGroup)
            {
                case FR_QTR: return &asfreq_AtoQ;
                case FR_MTH: return &asfreq_AtoM;
                case FR_WK: return &asfreq_AtoW;
                case FR_BUS: return &asfreq_AtoB;
                case FR_DAY: return &asfreq_AtoD;
                case FR_HR: return &asfreq_AtoH;
                case FR_MIN: return &asfreq_AtoT;
                case FR_SEC: return &asfreq_AtoS;
                default: return &nofunc;
            }

        case FR_QTR:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_QtoA;
                case FR_MTH: return &asfreq_QtoM;
                case FR_WK: return &asfreq_QtoW;
                case FR_BUS: return &asfreq_QtoB;
                case FR_DAY: return &asfreq_QtoD;
                case FR_HR: return &asfreq_QtoH;
                case FR_MIN: return &asfreq_QtoT;
                case FR_SEC: return &asfreq_QtoS;
                default: return &nofunc;
            }

        case FR_MTH:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_MtoA;
                case FR_QTR: return &asfreq_MtoQ;
                case FR_WK: return &asfreq_MtoW;
                case FR_BUS: return &asfreq_MtoB;
                case FR_DAY: return &asfreq_MtoD;
                case FR_HR: return &asfreq_MtoH;
                case FR_MIN: return &asfreq_MtoT;
                case FR_SEC: return &asfreq_MtoS;
                default: return &nofunc;
            }

        case FR_WK:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_WtoA;
                case FR_QTR: return &asfreq_WtoQ;
                case FR_MTH: return &asfreq_WtoM;
                case FR_WK: return &asfreq_WtoW;
                case FR_BUS: return &asfreq_WtoB;
                case FR_DAY: return &asfreq_WtoD;
                case FR_HR: return &asfreq_WtoH;
                case FR_MIN: return &asfreq_WtoT;
                case FR_SEC: return &asfreq_WtoS;
                default: return &nofunc;
            }

        case FR_BUS:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_BtoA;
                case FR_QTR: return &asfreq_BtoQ;
                case FR_MTH: return &asfreq_BtoM;
                case FR_WK: return &asfreq_BtoW;
                case FR_DAY: return &asfreq_BtoD;
                case FR_HR: return &asfreq_BtoH;
                case FR_MIN: return &asfreq_BtoT;
                case FR_SEC: return &asfreq_BtoS;
                default: return &nofunc;
            }

        case FR_DAY:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_DtoA;
                case FR_QTR: return &asfreq_DtoQ;
                case FR_MTH: return &asfreq_DtoM;
                case FR_WK: return &asfreq_DtoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_DtoB_forConvert; }
                    else            { return &asfreq_DtoB; }
                case FR_DAY: return &asfreq_DtoD;
                case FR_HR: return &asfreq_DtoH;
                case FR_MIN: return &asfreq_DtoT;
                case FR_SEC: return &asfreq_DtoS;
                default: return &nofunc;
            }

        case FR_HR:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_HtoA;
                case FR_QTR: return &asfreq_HtoQ;
                case FR_MTH: return &asfreq_HtoM;
                case FR_WK: return &asfreq_HtoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_HtoB_forConvert; }
                    else            { return &asfreq_HtoB; }
                case FR_DAY: return &asfreq_HtoD;
                case FR_MIN: return &asfreq_HtoT;
                case FR_SEC: return &asfreq_HtoS;
                default: return &nofunc;
            }

        case FR_MIN:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_TtoA;
                case FR_QTR: return &asfreq_TtoQ;
                case FR_MTH: return &asfreq_TtoM;
                case FR_WK: return &asfreq_TtoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_TtoB_forConvert; }
                    else            { return &asfreq_TtoB; }
                case FR_DAY: return &asfreq_TtoD;
                case FR_HR: return &asfreq_TtoH;
                case FR_SEC: return &asfreq_TtoS;
                default: return &nofunc;
            }

        case FR_SEC:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_StoA;
                case FR_QTR: return &asfreq_StoQ;
                case FR_MTH: return &asfreq_StoM;
                case FR_WK: return &asfreq_StoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_StoB_forConvert; }
                    else            { return &asfreq_StoB; }
                case FR_DAY: return &asfreq_StoD;
                case FR_HR: return &asfreq_StoH;
                case FR_MIN: return &asfreq_StoT;
                default: return &nofunc;
            }
        default: return &nofunc;
    }
}

/*
Helper function for cseries_convert:
    determine the size of the second dimension for the resulting
    converted array
*/
static long get_height(int fromFreq, int toFreq) {

    int maxBusDaysPerYear, maxBusDaysPerQuarter, maxBusDaysPerMonth;
    int maxDaysPerYear, maxDaysPerQuarter, maxDaysPerMonth;

    int fromGroup = (fromFreq/1000)*1000;
    int toGroup = (toFreq/1000)*1000;

    maxBusDaysPerYear = 262;
    maxBusDaysPerQuarter = 66;
    maxBusDaysPerMonth = 23;

    maxDaysPerYear = 366;
    maxDaysPerQuarter = 92;
    maxDaysPerMonth = 31;



    switch(fromGroup)
    {
        case FR_ANN: return 1;
        case FR_QTR:
            switch(toGroup)
            {
                case FR_ANN: return 4;
                default: return 1;
            }
        case FR_MTH: //monthly
            switch(toGroup)
            {
                case FR_ANN: return 12;
                case FR_QTR: return 3;
                default: return 1;
            }
        case FR_WK: //weekly
            switch(toGroup)
            {
                case FR_ANN: return 53;
                case FR_QTR: return 13;
                case FR_MTH: return 4;
                default: return 1;
            }
        case FR_BUS: //business
            switch(toGroup)
            {
                case FR_ANN: return maxBusDaysPerYear;;
                case FR_QTR: return maxBusDaysPerQuarter;
                case FR_MTH: return maxBusDaysPerMonth;
                case FR_WK: return 5;
                default: return 1;
            }
        case FR_DAY: //daily
            switch(toGroup)
            {
                case FR_ANN: return maxDaysPerYear;;
                case FR_QTR: return maxDaysPerQuarter;
                case FR_MTH: return maxDaysPerMonth;
                case FR_WK: return 7;
                default: return 1;
            }
        case FR_HR: //hourly
            switch(toGroup)
            {
                case FR_ANN: return 24 * maxDaysPerYear;;
                case FR_QTR: return 24 * maxDaysPerQuarter;
                case FR_MTH: return 24 * maxDaysPerMonth;
                case FR_WK: return 24 * 7;
                case FR_DAY: return 24;
                case FR_BUS: return 24;
                default: return 1;
            }
        case FR_MIN: //minutely
            switch(toGroup)
            {
                case FR_ANN: return 24 * 60 * maxDaysPerYear;;
                case FR_QTR: return 24 * 60 * maxDaysPerQuarter;
                case FR_MTH: return 24 * 60 * maxDaysPerMonth;
                case FR_WK: return 24 * 60 * 7;
                case FR_DAY: return 24 * 60;
                case FR_BUS: return 24 * 60;
                case FR_HR: return 60;
                default: return 1;
            }
        case FR_SEC: //minutely
            switch(toGroup)
            {
                case FR_ANN: return 24 * 60 * 60 * maxDaysPerYear;;
                case FR_QTR: return 24 * 60 * 60 * maxDaysPerQuarter;
                case FR_MTH: return 24 * 60 * 60 * maxDaysPerMonth;
                case FR_WK: return 24 * 60 * 60 * 7;
                case FR_DAY: return 24 * 60 * 60;
                case FR_BUS: return 24 * 60 * 60;
                case FR_HR: return 60 * 60;
                case FR_MIN: return 60;
                default: return 1;
            }
        default: return 1;
    }
}

static void *get_asfreq_info(int fromFreq, int toFreq, struct asfreq_info *af_info) {

    int maxBusDaysPerYear, maxBusDaysPerQuarter, maxBusDaysPerMonth;
    int maxDaysPerYear, maxDaysPerQuarter, maxDaysPerMonth;

    int fromGroup = (fromFreq/1000)*1000;
    int toGroup = (toFreq/1000)*1000;

    switch(fromGroup)
    {
        case FR_WK: {
            af_info->from_week_end = fromFreq - fromGroup;
        } break;
    }

    switch(toGroup)
    {
        case FR_WK: {
            af_info->to_week_end = toFreq - toGroup;
        } break;
    }

}

static char cseries_convert_doc[] = "";
static PyObject *
cseries_convert(PyObject *self, PyObject *args)
{
    PyObject *arrayTest;
    PyArrayObject *array, *newArray;
    PyArrayObject *mask, *newMask;

    PyObject *returnVal = NULL;
    PyObject *start_index_retval;

    long startIndex;
    long newStart, newStartTemp;
    long newEnd, newEndTemp;
    long newLen, newHeight;
    long i, currIndex, prevIndex;
    long nd;
    npy_intp *dim, *newIdx;
    long currPerLen;
    char *position;
    int fromFreq, toFreq;
    char relation;
    struct asfreq_info af_info;

    PyObject *val, *valMask;

    long (*asfreq_main)(long, char, struct asfreq_info) = NULL;
    long (*asfreq_endpoints)(long, char, struct asfreq_info) = NULL;
    long (*asfreq_reverse)(long, char, struct asfreq_info) = NULL;

    returnVal = PyDict_New();

    if (!PyArg_ParseTuple(args, "OiislO:convert(array, fromfreq, tofreq, position, startIndex, mask)", &array, &fromFreq, &toFreq, &position, &startIndex, &mask)) return NULL;

    if (toFreq == fromFreq)
    {
        PyDict_SetItemString(returnVal, "values", (PyObject*)array);
        PyDict_SetItemString(returnVal, "mask", (PyObject*)mask);

        Py_DECREF(array);
        Py_DECREF(mask);

        return returnVal;
    }

    switch(position[0])
    {
        case 'S':
            // start -> before
            relation = 'B';
            break;
        case 'E':
            // end -> after
            relation = 'A';
            break;
        default:
            return NULL;
            break;
    }

    get_asfreq_info(fromFreq, toFreq, &af_info);

    asfreq_main = get_asfreq_func(fromFreq, toFreq, 1);
    asfreq_endpoints = get_asfreq_func(fromFreq, toFreq, 0);

    //convert start index to new frequency
    CHECK_ASFREQ(newStartTemp = asfreq_main(startIndex, 'B', af_info));
    if (newStartTemp < 1) {
        CHECK_ASFREQ(newStart = asfreq_endpoints(startIndex, 'A', af_info));
    }
    else { newStart = newStartTemp; }

    //convert end index to new frequency
    CHECK_ASFREQ(newEndTemp = asfreq_main(startIndex+array->dimensions[0]-1, 'A', af_info));
    if (newEndTemp < 1) {
        CHECK_ASFREQ(newEnd = asfreq_endpoints(startIndex+array->dimensions[0]-1, 'B', af_info));
    }
    else { newEnd = newEndTemp; }

    if (newStart < 1) {
        PyErr_SetString(PyExc_ValueError, "start_date outside allowable range for destination frequency");
        return NULL;
    }

    newLen = newEnd - newStart + 1;
    newHeight = get_height(fromFreq, toFreq);

    if (newHeight > 1) {
        long tempval;
        asfreq_reverse = get_asfreq_func(toFreq, fromFreq, 0);
        CHECK_ASFREQ(tempval = asfreq_reverse(newStart, 'B', af_info));
        currPerLen = startIndex - tempval;

        nd = 2;
        dim = PyDimMem_NEW(nd);
        dim[0] = (npy_intp)newLen;
        dim[1] = (npy_intp)newHeight;
    } else {
        nd = 1;
        dim = PyDimMem_NEW(nd);
        dim[0] = (npy_intp)newLen;
    }

    newIdx = PyDimMem_NEW(nd);
    arrayTest = PyArray_SimpleNew(nd, dim, array->descr->type_num);
    if (arrayTest == NULL) { return NULL; }
    newArray = (PyArrayObject*)arrayTest;
    newMask  = (PyArrayObject*)PyArray_SimpleNew(nd, dim, mask->descr->type_num);

    PyDimMem_FREE(dim);

    PyArray_FILLWBYTE(newArray,0);
    PyArray_FILLWBYTE(newMask,1);

    prevIndex = newStart;

    //set values in the new array
    for (i = 0; i < array->dimensions[0]; i++) {

        val = PyArray_GETITEM(array, PyArray_GetPtr(array, &i));
        valMask = PyArray_GETITEM(mask, PyArray_GetPtr(mask, &i));

        CHECK_ASFREQ(currIndex = asfreq_main(startIndex + i, relation, af_info));
        newIdx[0] = currIndex-newStart;

        if (newHeight > 1) {

                if (currIndex != prevIndex)
                {
                    //reset period length
                    currPerLen = 0;
                    prevIndex = currIndex;
                }

                newIdx[1] = currPerLen;
                currPerLen++;
        }

        if (newIdx[0] > -1) {
            PyArray_SETITEM(newArray, PyArray_GetPtr(newArray, newIdx), val);
            PyArray_SETITEM(newMask, PyArray_GetPtr(newMask, newIdx), valMask);
        }

        Py_DECREF(val);
        Py_DECREF(valMask);
    }

    PyDimMem_FREE(newIdx);

    start_index_retval = (PyObject*)PyInt_FromLong(newStart);

    PyDict_SetItemString(returnVal, "values", (PyObject*)newArray);
    PyDict_SetItemString(returnVal, "mask", (PyObject*)newMask);
    PyDict_SetItemString(returnVal, "startindex", start_index_retval);

    Py_DECREF(newArray);
    Py_DECREF(newMask);
    Py_DECREF(start_index_retval);

    return returnVal;
}

static char cseries_asfreq_doc[] = "";
static PyObject *
cseries_asfreq(PyObject *self, PyObject *args)
{
    PyArrayObject *fromDates, *toDates;
    PyArrayIterObject *iterFrom, *iterTo;
    PyObject *fromDateObj, *toDateObj;
    char *relation;
    int fromFreq, toFreq;
    long fromDate, toDate;
    long (*asfreq_main)(long, char, struct asfreq_info) = NULL;
    struct asfreq_info af_info;

    if (!PyArg_ParseTuple(args, "Oiis:asfreq(fromDates, fromfreq, tofreq, relation)", &fromDates, &fromFreq, &toFreq, &relation)) return NULL;

    get_asfreq_info(fromFreq, toFreq, &af_info);

    asfreq_main = get_asfreq_func(fromFreq, toFreq, 0);

    toDates = (PyArrayObject *)PyArray_Copy(fromDates);

    iterFrom = (PyArrayIterObject *)PyArray_IterNew((PyObject *)fromDates);
    if (iterFrom == NULL) return NULL;

    iterTo = (PyArrayIterObject *)PyArray_IterNew((PyObject *)toDates);
    if (iterTo == NULL) return NULL;

    while (iterFrom->index < iterFrom->size) {

        fromDateObj = PyArray_GETITEM(fromDates, iterFrom->dataptr);
        fromDate = PyInt_AsLong(fromDateObj);
        CHECK_ASFREQ(toDate = asfreq_main(fromDate, relation[0], af_info));
        toDateObj = PyInt_FromLong(toDate);

        PyArray_SETITEM(toDates, iterTo->dataptr, toDateObj);

        Py_DECREF(fromDateObj);
        Py_DECREF(toDateObj);

        PyArray_ITER_NEXT(iterFrom);
        PyArray_ITER_NEXT(iterTo);
    }

    Py_DECREF(iterFrom);
    Py_DECREF(iterTo);

    return (PyObject *)toDates;

}

static int dInfo_year(struct date_info dateObj)    { return dateObj.year; }
static int dInfo_quarter(struct date_info dateObj) { return dateObj.quarter; }
static int dInfo_month(struct date_info dateObj)   { return dateObj.month; }
static int dInfo_day(struct date_info dateObj)     { return dateObj.day; }
static int dInfo_day_of_year(struct date_info dateObj) { return dateObj.day_of_year; }
static int dInfo_day_of_week(struct date_info dateObj) { return dateObj.day_of_week; }
static int dInfo_week(struct date_info dateObj)     { return dInfoCalc_ISOWeek(dateObj); }
static int dInfo_hour(struct date_info dateObj)     { return dateObj.hour; }
static int dInfo_minute(struct date_info dateObj)     { return dateObj.minute; }
static int dInfo_second(struct date_info dateObj)     { return (int)dateObj.second; }

static double getAbsTime(int freq, long dailyDate, long originalDate) {

    long startOfDay, periodsPerDay;

    switch(freq)
    {
        case FR_HR:
            periodsPerDay = 24;
            break;
        case FR_MIN:
            periodsPerDay = 24*60;
            break;
        case FR_SEC:
            periodsPerDay = 24*60*60;
            break;
        default:
            return 0;
    }

    startOfDay = asfreq_DtoHIGHFREQ(dailyDate, 'B', periodsPerDay);
    return (24*60*60)*((double)(originalDate - startOfDay))/((double)periodsPerDay);
}



static char cseries_getDateInfo_doc[] = "";
static PyObject *
cseries_getDateInfo(PyObject *self, PyObject *args)
{
    int freq;
    char *info;

    PyArrayObject *array;
    PyArrayObject *newArray;
    PyArrayIterObject *iterSource, *iterResult;
    struct date_info convDate;

    PyObject *val;
    long dateNum, dInfo;
    long absdate;
    double abstime;

    long (*toDaily)(long, char, struct asfreq_info) = NULL;
    int (*getDateInfo)(struct date_info) = NULL;
    struct asfreq_info af_info;

    if (!PyArg_ParseTuple(args, "Ois:getDateInfo(array, freq, info)", &array, &freq, &info)) return NULL;
    newArray = (PyArrayObject *)PyArray_Copy(array);

    iterSource = (PyArrayIterObject *)PyArray_IterNew((PyObject *)array);
    iterResult = (PyArrayIterObject *)PyArray_IterNew((PyObject *)newArray);

    toDaily = get_asfreq_func(freq, FR_DAY, 0);
    get_asfreq_info(freq, FR_DAY, &af_info);

    switch(*info)
    {
        case 'Y': //year
            getDateInfo = &dInfo_year;
            break;
        case 'Q': //quarter
            getDateInfo = &dInfo_quarter;
            break;
        case 'M': //month
            getDateInfo = &dInfo_month;
            break;
        case 'D': //day
            getDateInfo = &dInfo_day;
            break;
        case 'R': //day of year
            getDateInfo = &dInfo_day_of_year;
            break;
        case 'W': //day of week
            getDateInfo = &dInfo_day_of_week;
            break;
        case 'I': //week of year
            getDateInfo = &dInfo_week;
            break;
        case 'H': //hour
            getDateInfo = &dInfo_hour;
            break;
        case 'T': //minute
            getDateInfo = &dInfo_minute;
            break;
        case 'S': //second
            getDateInfo = &dInfo_second;
            break;
        default:
            return NULL;
    }

    while (iterSource->index < iterSource->size) {

        val = PyArray_GETITEM(array, iterSource->dataptr);
        dateNum = PyInt_AsLong(val);
        Py_DECREF(val);
        CHECK_ASFREQ(absdate = toDaily(dateNum, 'B', af_info));
        abstime = getAbsTime(freq, absdate, dateNum);

        if(dInfoCalc_SetFromAbsDateTime(&convDate, absdate, abstime, GREGORIAN_CALENDAR)) return NULL;
        dInfo = getDateInfo(convDate);

        PyArray_SETITEM(newArray, iterResult->dataptr, PyInt_FromLong(dInfo));

        PyArray_ITER_NEXT(iterSource);
        PyArray_ITER_NEXT(iterResult);
    }

    Py_DECREF(iterSource);
    Py_DECREF(iterResult);

    return (PyObject *) newArray;
}

static char *str_replace(const char *s, const char *old, const char *new)
{
    char *ret;
    int i, count = 0;
    size_t newlen = strlen(new);
    size_t oldlen = strlen(old);

    for (i = 0; s[i] != '\0'; i++) {
        if (strstr(&s[i], old) == &s[i]) {
           count++;
           i += oldlen - 1;
        }
    }

    ret = malloc(i + 1 + count * (newlen - oldlen));
    if (ret == NULL) return NULL;

    i = 0;
    while (*s) {
        if (strstr(s, old) == s) {
            strcpy(&ret[i], new);
            i += newlen;
            s += oldlen;
        } else {
            ret[i++] = *s++;
        }
    }
    ret[i] = '\0';

    return ret;
}

static char cseries_strfmt_doc[] = "";
static PyObject *
cseries_strfmt(PyObject *self, PyObject *args)
{

    char *orig_fmt_str, *fmt_str, *q_loc;
    char *result;
    char place_holder[] = "^`";
    struct tm c_date;
    struct date_info tempDate;
    int result_len;
    PyObject *date, *py_result;

    if (!PyArg_ParseTuple(args, "Os:strfmt(datetime, fmt_str)", &date, &orig_fmt_str)) return NULL;

    if (dInfoCalc_SetFromDateAndTime(&tempDate,
                                    PyDateTime_GET_YEAR(date),
                                    PyDateTime_GET_MONTH(date),
                                    PyDateTime_GET_DAY(date),
                                    PyDateTime_DATE_GET_HOUR(date),
                                    PyDateTime_DATE_GET_MINUTE(date),
                                    PyDateTime_DATE_GET_SECOND(date),
                                    GREGORIAN_CALENDAR)) return NULL;

    /* We need to modify the fmt_str passed in to handle our special syntax for quarters.
       We can't modify the string passed in directly, so we must make a copy. */
    fmt_str = malloc((strlen(orig_fmt_str) + 1)*sizeof(char));
    strcpy(fmt_str, orig_fmt_str);

    if ((q_loc = strstr(fmt_str,"%q")) != NULL) {
        q_loc = strstr(fmt_str,"%q");
        strncpy (q_loc,place_holder,2);
    }

    c_date.tm_sec = (int)tempDate.second;
    c_date.tm_min = tempDate.minute;
    c_date.tm_hour = tempDate.hour;
    c_date.tm_mday = tempDate.day;
    c_date.tm_mon = tempDate.month - 1;
    c_date.tm_year = tempDate.year - 1900;
    c_date.tm_wday = tempDate.day_of_week;
    c_date.tm_yday = tempDate.day_of_year;
    c_date.tm_isdst = -1;

    result_len = strlen(orig_fmt_str) + 50;

    result = malloc(result_len * sizeof(char));

    strftime(result, result_len, fmt_str, &c_date);

    if (q_loc != NULL) {
        char *alt_result;
        char q_str[2];

        sprintf(q_str, "%i", tempDate.quarter);
        alt_result = str_replace(result, place_holder, q_str);
        py_result = PyString_FromString(alt_result);
        free(result);
        free(alt_result);
    } else {
        py_result = PyString_FromString(result);
        free(result);
    }

    return py_result;

}


///////////////////////////////////////////////////////////////////////

static PyMethodDef cseries_methods[] = {
    {"strfmt", cseries_strfmt, METH_VARARGS, cseries_strfmt_doc},
    {"convert", cseries_convert, METH_VARARGS, cseries_convert_doc},
    {"asfreq", cseries_asfreq, METH_VARARGS, cseries_asfreq_doc},
    {"getDateInfo", cseries_getDateInfo, METH_VARARGS, cseries_getDateInfo_doc},
    {NULL, NULL}
};

PyMODINIT_FUNC
initcseries(void)
{
    PyObject *m, *TSER_CONSTANTS;
    m = Py_InitModule3("cseries", cseries_methods, cseries_doc);
    import_array();
    PyDateTime_IMPORT;

    TSER_CONSTANTS = PyDict_New();

    // Add all the frequency constants to a python dictionary
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_ANN", FR_ANN);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_QTR", FR_QTR);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_MTH", FR_MTH);

    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WK", FR_WK);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WKSUN", FR_WKSUN);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WKSAT", FR_WKSAT);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WKFRI", FR_WKFRI);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WKTHU", FR_WKTHU);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WKWED", FR_WKWED);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WKTUE", FR_WKTUE);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_WKMON", FR_WKMON);

    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_BUS", FR_BUS);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_DAY", FR_DAY);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_HR", FR_HR);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_MIN", FR_MIN);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_SEC", FR_SEC);
    ADD_INT_TO_DICT(TSER_CONSTANTS, "FR_UND", FR_UND);

    PyModule_AddObject(m, "TSER_CONSTANTS", TSER_CONSTANTS);
}