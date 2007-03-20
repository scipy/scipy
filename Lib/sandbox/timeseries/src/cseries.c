#include <Python.h>
#include <datetime.h>
#include <structmember.h>
#include <time.h>
#include "arrayobject.h"

static char cseries_doc[] = "Speed sensitive time series operations";

#define FR_ANN  1000  /* Annual */
#define FR_ANNDEC  FR_ANN  /* Annual - December year end*/
#define FR_ANNJAN  1001  /* Annual - January year end*/
#define FR_ANNFEB  1002  /* Annual - February year end*/
#define FR_ANNMAR  1003  /* Annual - March year end*/
#define FR_ANNAPR  1004  /* Annual - April year end*/
#define FR_ANNMAY  1005  /* Annual - May year end*/
#define FR_ANNJUN  1006  /* Annual - June year end*/
#define FR_ANNJUL  1007  /* Annual - July year end*/
#define FR_ANNAUG  1008  /* Annual - August year end*/
#define FR_ANNSEP  1009  /* Annual - September year end*/
#define FR_ANNOCT  1010  /* Annual - October year end*/
#define FR_ANNNOV  1011  /* Annual - November year end*/

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

#define INT_ERR_CODE -999

#define HIGHFREQ_ORIG 719163

#define CHECK_ASFREQ(result) if ((result) == INT_ERR_CODE) return NULL

static int get_freq_group(int freq) {
    return (freq/1000)*1000;
}

struct asfreq_info{
    int from_week_end; //day the week ends on in the "from" frequency
    int to_week_end; //day the week ends on in the "to" frequency

    int from_a_year_end; //month the year ends on in the "from" frequency
    int to_a_year_end; //month the year ends on in the "to" frequency

};

static struct asfreq_info NULL_AF_INFO;

/*********************************************************
** Python callbacks. These functions must be called by  **
** the module __init__ script                           **
*********************************************************/
static PyObject *
set_callback(PyObject *args, PyObject **callback)
{
    PyObject *result = NULL;
    PyObject *temp;

    if (PyArg_ParseTuple(args, "O:set_callback", &temp)) {

        if (!PyCallable_Check(temp)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return NULL;
        }

        Py_XINCREF(temp);        // Add a reference to new callback
        Py_XDECREF(*callback);  // Dispose of previous callback
        *callback = temp;       // Remember new callback
        // Boilerplate to return "None"
        Py_INCREF(Py_None);
        result = Py_None;
    }
    return result;
}

static PyObject *DateFromString = NULL;
static char set_callback_DateFromString_doc[] =
"set DateFromString function python callback";
static PyObject *
set_callback_DateFromString(PyObject *dummy, PyObject *args) {
    return set_callback(args, &DateFromString);
}

static PyObject *DateTimeFromString = NULL;
static char set_callback_DateTimeFromString_doc[] =
"set DateTimeFromString function python callback";
static PyObject *
set_callback_DateTimeFromString(PyObject *dummy, PyObject *args) {
    return set_callback(args, &DateTimeFromString);
}

/*********************************************************/

static char *
str_uppercase(char *str) {
    if (str) {
        int i, len=strlen(str);
        char *result;
        if((result = malloc((len + 1)*sizeof(char))) == NULL) {
            return (char *)PyErr_NoMemory();
        }
        strcpy(result, str);

        for (i=0;i<len;i++) {
            switch(result[i])
            {
                case 'a': { result[i]='A'; break; }
                case 'b': { result[i]='B'; break; }
                case 'c': { result[i]='C'; break; }
                case 'd': { result[i]='D'; break; }
                case 'e': { result[i]='E'; break; }
                case 'f': { result[i]='F'; break; }
                case 'g': { result[i]='G'; break; }
                case 'h': { result[i]='H'; break; }
                case 'i': { result[i]='I'; break; }
                case 'j': { result[i]='J'; break; }
                case 'k': { result[i]='K'; break; }
                case 'l': { result[i]='L'; break; }
                case 'm': { result[i]='M'; break; }
                case 'n': { result[i]='N'; break; }
                case 'o': { result[i]='O'; break; }
                case 'p': { result[i]='P'; break; }
                case 'q': { result[i]='Q'; break; }
                case 'r': { result[i]='R'; break; }
                case 's': { result[i]='S'; break; }
                case 't': { result[i]='T'; break; }
                case 'u': { result[i]='U'; break; }
                case 'v': { result[i]='V'; break; }
                case 'w': { result[i]='W'; break; }
                case 'x': { result[i]='X'; break; }
                case 'y': { result[i]='Y'; break; }
                case 'z': { result[i]='Z'; break; }
            }
        }

        return result;
    } else { return NULL; }
}

static char *
str_replace(const char *s, const char *old, const char *new) {
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
    if (ret == NULL) {return (char *)PyErr_NoMemory();}

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
int dInfoCalc_ISOWeek(struct date_info *dinfo)
{
    int week;

    /* Estimate */
    week = (dinfo->day_of_year-1) - dinfo->day_of_week + 3;
    if (week >= 0) week = week / 7 + 1;

    /* Verify */
    if (week < 0) {
        /* The day lies in last week of the previous year */
        if ((week > -2) ||
            (week == -2 && dInfoCalc_Leapyear(dinfo->year-1, dinfo->calendar)))
            week = 53;
        else
            week = 52;
    } else if (week == 53) {
    /* Check if the week belongs to year or year+1 */
        if (31-dinfo->day + dinfo->day_of_week < 3) {
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





///////////////////////////////////////////////////////////////////////

// helpers for frequency conversion routines //

static long DtoB_weekday(long fromDate) { return (((fromDate) / 7) * 5) + (fromDate)%7; }

static long DtoB_WeekendToMonday(long absdate, int day_of_week) {

    if (day_of_week > 4) {
        //change to Monday after weekend
        absdate += (7 - day_of_week);
    }
    return DtoB_weekday(absdate);
}

static long DtoB_WeekendToFriday(long absdate, int day_of_week) {

    if (day_of_week > 4) {
        //change to friday before weekend
        absdate -= (day_of_week - 4);
    }
    return DtoB_weekday(absdate);
}

static long absdate_from_ymd(int y, int m, int d) {
    struct date_info tempDate;
    if (dInfoCalc_SetFromDateAndTime(&tempDate, y, m, d, 0, 0, 0, GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    return tempDate.absdate;
}


///////////////////////////////////////////////

// frequency specifc conversion routines
// each function must take an integer fromDate and a char relation ('B' or 'A' for 'BEFORE' or 'AFTER')

//************ FROM DAILY ***************

static long asfreq_DtoA(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    if (dinfo.month > af_info->to_a_year_end) { return (long)(dinfo.year + 1); }
    else { return (long)(dinfo.year); }
}


static long asfreq_DtoQ(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    return (long)((dinfo.year - 1) * 4 + dinfo.quarter);
}

static long asfreq_DtoM(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    return (long)((dinfo.year - 1) * 12 + dinfo.month);
}

static long asfreq_DtoW(long fromDate, char relation, struct asfreq_info *af_info) {
    return (fromDate - (1 + af_info->to_week_end))/7 + 1;
}

static long asfreq_DtoB(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'B') {
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week);
    } else {
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week);
    }
}

static long asfreq_DtoB_forConvert(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (dinfo.day_of_week > 4) {
        return -1;
    } else {
        return DtoB_weekday(fromDate);
    }
}

// needed for getDateInfo function
static long asfreq_DtoD(long fromDate, char relation, struct asfreq_info *af_info) { return fromDate; }

static long asfreq_DtoHIGHFREQ(long fromDate, char relation, long periodsPerDay) {
    if (fromDate >= HIGHFREQ_ORIG) {
        if (relation == 'B') { return (fromDate - HIGHFREQ_ORIG)*(periodsPerDay) + 1; }
        else                 { return (fromDate - HIGHFREQ_ORIG + 1)*(periodsPerDay); }
    } else { return -1; }
}

static long asfreq_DtoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24); }
static long asfreq_DtoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24*60); }
static long asfreq_DtoS(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24*60*60); }

//************ FROM SECONDLY ***************

static long asfreq_StoD(long fromDate, char relation, struct asfreq_info *af_info)
    { return (fromDate - 1)/(60*60*24) + HIGHFREQ_ORIG; }

static long asfreq_StoA(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_StoQ(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_StoM(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_StoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_StoB(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_StoB_forConvert(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoB_forConvert(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_StoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return (fromDate - 1)/60 + 1; }
static long asfreq_StoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return (fromDate - 1)/(60*60) + 1; }

//************ FROM MINUTELY ***************

static long asfreq_TtoD(long fromDate, char relation, struct asfreq_info *af_info)
    { return (fromDate - 1)/(60*24) + HIGHFREQ_ORIG; }

static long asfreq_TtoA(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_TtoQ(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_TtoM(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_TtoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_TtoB(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_TtoB_forConvert(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoB_forConvert(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_TtoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return (fromDate - 1)/60 + 1; }
static long asfreq_TtoS(long fromDate, char relation, struct asfreq_info *af_info) {
    if (relation == 'B') {  return fromDate*60 - 59; }
    else                 {  return fromDate*60;      }}

//************ FROM HOURLY ***************

static long asfreq_HtoD(long fromDate, char relation, struct asfreq_info *af_info)
    { return (fromDate - 1)/24 + HIGHFREQ_ORIG; }
static long asfreq_HtoA(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_HtoQ(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_HtoM(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_HtoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_HtoB(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_HtoB_forConvert(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoB_forConvert(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

// calculation works out the same as TtoS, so we just call that function for HtoT
static long asfreq_HtoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_TtoS(fromDate, relation, &NULL_AF_INFO); }
static long asfreq_HtoS(long fromDate, char relation, struct asfreq_info *af_info) {
    if (relation == 'B') {  return fromDate*60*60 - 60*60 + 1; }
    else                 {  return fromDate*60*60;             }}

//************ FROM BUSINESS ***************

static long asfreq_BtoD(long fromDate, char relation, struct asfreq_info *af_info)
    { return ((fromDate-1)/5)*7 + (fromDate-1)%5 + 1; }

static long asfreq_BtoA(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_BtoQ(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_BtoM(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_BtoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_BtoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_BtoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_BtoS(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM WEEKLY ***************

static long asfreq_WtoD(long fromDate, char relation, struct asfreq_info *af_info) {
    if (relation == 'B') { return fromDate * 7 - 6 + af_info->from_week_end;}
    else                 { return fromDate * 7 + af_info->from_week_end; }
}

static long asfreq_WtoA(long fromDate, char relation, struct asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_WtoD(fromDate, 'A', af_info), relation, af_info); }
static long asfreq_WtoQ(long fromDate, char relation, struct asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_WtoD(fromDate, 'A', af_info), relation, &NULL_AF_INFO); }
static long asfreq_WtoM(long fromDate, char relation, struct asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_WtoD(fromDate, 'A', af_info), relation, &NULL_AF_INFO); }

static long asfreq_WtoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_WtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }

static long asfreq_WtoB(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_WtoD(fromDate, relation, af_info),
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'B') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static long asfreq_WtoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_WtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_WtoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_WtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_WtoS(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_WtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }

//************ FROM MONTHLY ***************

static void MtoD_ym(long fromDate, long *y, long *m) {
    *y = (fromDate - 1) / 12 + 1;
    *m = fromDate - 12 * (*y) - 1;
}

static long asfreq_MtoD(long fromDate, char relation, struct asfreq_info *af_info) {

    long y, m, absdate;

    if (relation == 'B') {
        MtoD_ym(fromDate, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate;
    } else {
        MtoD_ym(fromDate+1, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate-1;
    }
}

static long asfreq_MtoA(long fromDate, char relation, struct asfreq_info *af_info) {
    long year = (fromDate - 1) / 12 + 1;
    long month = fromDate - (year - 1)*12;
    if (month > af_info->to_a_year_end) { return year + 1; }
    else { return year; }
}

static long asfreq_MtoQ(long fromDate, char relation, struct asfreq_info *af_info)
    { return (fromDate - 1) / 3 + 1; }

static long asfreq_MtoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_MtoB(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_MtoD(fromDate, relation, &NULL_AF_INFO),
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'B') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static long asfreq_MtoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_MtoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_MtoS(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM QUARTERLY ***************

static void QtoD_ym(long fromDate, long *y, long *m) {
    *y = (fromDate - 1) / 4 + 1;
    *m = (fromDate + 4) * 3 - 12 * (*y) - 2;
}

static long asfreq_QtoD(long fromDate, char relation, struct asfreq_info *af_info) {

    long y, m, absdate;

    if (relation == 'B') {
        QtoD_ym(fromDate, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate;
    } else {
        QtoD_ym(fromDate+1, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate - 1;
    }
}

static long asfreq_QtoA(long fromDate, char relation, struct asfreq_info *af_info) {
    long year = (fromDate - 1) / 4 + 1;
    long quarter = fromDate - (year - 1)*4;
    if (quarter > af_info->to_a_year_end) { return year + 1; }
    else { return year; }
}

static long asfreq_QtoM(long fromDate, char relation, struct asfreq_info *af_info) {
    if (relation == 'B') { return fromDate * 3 - 2; }
    else {                 return fromDate  * 3; }
}

static long asfreq_QtoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_QtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_QtoB(long fromDate, char relation, struct asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_QtoD(fromDate, relation, &NULL_AF_INFO),
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'B') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}


static long asfreq_QtoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_QtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_QtoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_QtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_QtoS(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_QtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }


//************ FROM ANNUAL ***************

static long asfreq_AtoD(long fromDate, char relation, struct asfreq_info *af_info) {
    long absdate, year, final_adj;
    int month = (af_info->from_a_year_end + 1) % 12;

    if (relation == 'B') {
        if (af_info->from_a_year_end == 12) {year = fromDate;}
        else {year = fromDate - 1;}
        final_adj = 0;
    } else {
        if (af_info->from_a_year_end == 12) {year = fromDate+1;}
        else {year = fromDate;}
        final_adj = -1;
    }
    absdate = absdate_from_ymd(year, month, 1);
    if (absdate  == INT_ERR_CODE) return INT_ERR_CODE;
    return absdate + final_adj;
}

static long asfreq_AtoQ(long fromDate, char relation, struct asfreq_info *af_info) {

    long quarter = asfreq_MtoQ(af_info->from_a_year_end, 'B', &NULL_AF_INFO);

    if (relation == 'B') {
        return fromDate * 4 - 3 - (4 - quarter);
    } else {
        return fromDate * 4 - (4 - quarter);
    }
}

static long asfreq_AtoM(long fromDate, char relation, struct asfreq_info *af_info) {
    if (relation == 'B') {
        return fromDate * 12 - 11 - (12 - af_info->from_a_year_end);
    } else {
        return fromDate * 12 - (12 - af_info->from_a_year_end);
    }
}

static long asfreq_AtoW(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_AtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_AtoB(long fromDate, char relation, struct asfreq_info *af_info) {

    long year;
    int month = (af_info->from_a_year_end + 1) % 12;

    struct date_info dailyDate;

    if (relation == 'B') {

        if (af_info->from_a_year_end == 12) {year = fromDate;}
        else {year = fromDate - 1;}

        if (dInfoCalc_SetFromDateAndTime(&dailyDate,
                        year,month,1, 0, 0, 0,
                        GREGORIAN_CALENDAR)) return INT_ERR_CODE;
        return DtoB_WeekendToMonday(dailyDate.absdate, dailyDate.day_of_week);
    } else {
        long absdate;

        if (af_info->from_a_year_end == 12) {year = fromDate+1;}
        else {year = fromDate;}

        if (dInfoCalc_SetFromDateAndTime(&dailyDate,
                       year,month,1, 0, 0, 0,
                       GREGORIAN_CALENDAR)) return INT_ERR_CODE;

        absdate = dailyDate.absdate - 1;

        if(dInfoCalc_SetFromAbsDate(&dailyDate,
                      absdate,
                      GREGORIAN_CALENDAR)) return INT_ERR_CODE;

        return DtoB_WeekendToFriday(dailyDate.absdate, dailyDate.day_of_week);
    }
}

static long asfreq_AtoH(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_AtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_AtoT(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_AtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_AtoS(long fromDate, char relation, struct asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_AtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }

static long nofunc(long fromDate, char relation, struct asfreq_info *af_info) { return -1; }

// end of frequency specific conversion routines

// return a pointer to appropriate conversion function
static long (*get_asfreq_func(int fromFreq, int toFreq, int forConvert))(long, char, struct asfreq_info*) {

    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    if (fromGroup == FR_UND) { fromGroup = FR_DAY; }

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

    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    if (fromGroup == FR_UND) { fromGroup = FR_DAY; }

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

static int calc_a_year_end(int freq, int group) {
    int result = freq - group;
    if (result == 0) {return 12;}
    else {return result;}
}

static int calc_week_end(int freq, int group) {
    return freq - group;
}

static void get_asfreq_info(int fromFreq, int toFreq, struct asfreq_info *af_info) {

    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    switch(fromGroup)
    {
        case FR_WK: {
            af_info->from_week_end = calc_week_end(fromFreq, fromGroup);
        } break;
        case FR_ANN: {
            af_info->from_a_year_end = calc_a_year_end(fromFreq, fromGroup);
        } break;
    }

    switch(toGroup)
    {
        case FR_WK: {
            af_info->to_week_end = calc_week_end(toFreq, toGroup);
        } break;
        case FR_ANN: {
            af_info->to_a_year_end = calc_a_year_end(toFreq, toGroup);
        } break;
    }

}

static int dInfo_year(struct date_info *dateObj)    { return dateObj->year; }
static int dInfo_quarter(struct date_info *dateObj) { return dateObj->quarter; }
static int dInfo_month(struct date_info *dateObj)   { return dateObj->month; }
static int dInfo_day(struct date_info *dateObj)     { return dateObj->day; }
static int dInfo_day_of_year(struct date_info *dateObj) { return dateObj->day_of_year; }
static int dInfo_day_of_week(struct date_info *dateObj) { return dateObj->day_of_week; }
static int dInfo_week(struct date_info *dateObj)     { return dInfoCalc_ISOWeek(dateObj); }
static int dInfo_hour(struct date_info *dateObj)     { return dateObj->hour; }
static int dInfo_minute(struct date_info *dateObj)     { return dateObj->minute; }
static int dInfo_second(struct date_info *dateObj)     { return (int)dateObj->second; }

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

/************************************************************
** Date type definition
************************************************************/

typedef struct {
    PyObject_HEAD
    int freq; /* frequency of date */
    int value; /* integer representation of date */
    PyObject* cached_vals;
} DateObject;

/* Forward declarations */
static PyTypeObject DateType;
#define DateObject_Check(op) PyObject_TypeCheck(op, &DateType)

static void
DateObject_dealloc(DateObject* self) {
    Py_XDECREF(self->cached_vals);
    self->ob_type->tp_free((PyObject*)self);
}


static PyObject *freq_dict, *freq_dict_rev, *freq_constants;

#define DICT_SETINT_STRKEY(dict, key, val) \
    {PyObject *pyval = PyInt_FromLong(val); \
     PyDict_SetItemString(dict, key, pyval); \
     Py_DECREF(pyval); }

#define ADD_FREQ_CONSTANT(const_name, val) \
    DICT_SETINT_STRKEY(freq_constants, const_name, val)

#define INIT_FREQ(const_name, key, aliases) \
    {PyObject *pykey = PyInt_FromLong(key); \
     PyDict_SetItem(freq_dict, pykey, aliases); \
     PyDict_SetItemString(freq_constants, const_name, pykey); \
     Py_DECREF(pykey); \
     Py_DECREF(aliases); }


static int init_freq_group(int num_items, int num_roots, int group_const,
                            char item_abbrevs[][2][10], char group_prefixes[][10],
                            char item_const_names[][10]) {

    int i;

    for (i = 0; i < num_items; i++) {

        PyObject *aliases;
        int j, size, k;

        if (i == 0) { k = 3; } else { k = 2; }

        size = num_roots * k;

        aliases = PyTuple_New(size);

        for (j = 0; j < num_roots; j++) {
            PyObject *alias_v1, *alias_v2;
            char *root, *alt;

            if ((root = malloc((26) * sizeof(char))) == NULL) return INT_ERR_CODE;
            if ((alt = malloc((26) * sizeof(char))) == NULL) return INT_ERR_CODE;

            strcpy(root, group_prefixes[j]);
            strcpy(alt, group_prefixes[j]);

            if (i == 0) {
                PyObject *alias = PyString_FromString(root);
                PyTuple_SET_ITEM(aliases, j*k + 2, alias);
            }

            strcat(root, "-");
            strcat(root, item_abbrevs[i][0]);
            strcat(alt, "-");
            strcat(alt, item_abbrevs[i][1]);

            alias_v1 = PyString_FromString(root);
            alias_v2 = PyString_FromString(alt);

            free(root);
            free(alt);

            PyTuple_SET_ITEM(aliases, j*k, alias_v1);
            PyTuple_SET_ITEM(aliases, j*k + 1, alias_v2);
        }

        INIT_FREQ(item_const_names[i], group_const+i, aliases);
    }

    return 0;
}

/* take a dictionary with integer keys and tuples of strings for values,
   and populate a dictionary with all the strings as keys and integers
   for values */
static int reverse_dict(PyObject *source, PyObject *dest) {

    PyObject *key, *value;

    int pos = 0;

    while (PyDict_Next(source, &pos, &key, &value)) {
        PyObject *tuple_iter;
        PyObject *item;

        if((tuple_iter = PyObject_GetIter(value)) == NULL) return INT_ERR_CODE;

        while (item = PyIter_Next(tuple_iter)) {
            PyDict_SetItem(dest, item, key);
            Py_DECREF(item);
        }
        Py_DECREF(tuple_iter);
    }
    return 0;
}

static int build_freq_dict() {

    char ANN_prefixes[8][10] = { "A", "Y", "ANN", "ANNUAL", "ANNUALLY",
                                 "YR", "YEAR", "YEARLY" };
    char WK_prefixes[4][10] =  { "W", "WK", "WEEK", "WEEKLY" };

    /* Note: order of this array must match up with how the Annual
       frequency constants are lined up */
    char month_names[12][2][10] = {
        { "DEC", "DECEMBER" },
        { "JAN", "JANUARY" },
        { "FEB", "FEBRUARY" },
        { "MAR", "MARCH" },
        { "APR", "APRIL" },
        { "MAY", "MAY" },
        { "JUN", "JUNE" },
        { "JUL", "JULY" },
        { "AUG", "AUGUST" },
        { "SEP", "SEPTEMBER" },
        { "OCT", "OCTOBER" },
        { "NOV", "NOVEMBER" }};

    char ANN_const_names[12][10] = {
        "FR_ANNDEC",
        "FR_ANNJAN",
        "FR_ANNFEB",
        "FR_ANNMAR",
        "FR_ANNAPR",
        "FR_ANNMAY",
        "FR_ANNJUN",
        "FR_ANNJUL",
        "FR_ANNAUG",
        "FR_ANNSEP",
        "FR_ANNOCT",
        "FR_ANNNOV"};

    char day_names[7][2][10] = {
        { "SUN", "SUNDAY" },
        { "MON", "MONDAY" },
        { "TUE", "TUESDAY" },
        { "WED", "WEDNESDAY" },
        { "THU", "THURSDAY" },
        { "FRI", "FRIDAY" },
        { "SAT", "SATURDAY" }};

    char WK_const_names[7][10] = {
        "FR_WKSUN",
        "FR_WKMON",
        "FR_WKTUE",
        "FR_WKWED",
        "FR_WKTHU",
        "FR_WKFRI",
        "FR_WKSAT"};

    PyObject *aliases;

    freq_dict = PyDict_New();
    freq_dict_rev = PyDict_New();
    freq_constants = PyDict_New();

    aliases = Py_BuildValue("(ssss)", "Q", "QTR", "QUARTER", "QUARTERLY");
    INIT_FREQ("FR_QTR", FR_QTR, aliases);

    aliases = Py_BuildValue("(ssss)", "M", "MTH", "MONTH", "MONTHLY");
    INIT_FREQ("FR_MTH", FR_MTH, aliases);

    aliases = Py_BuildValue("(ssss)", "B", "BUS", "BUSINESS", "BUSINESSLY");
    INIT_FREQ("FR_BUS", FR_BUS, aliases);

    aliases = Py_BuildValue("(ssss)", "D", "DAY", "DLY", "DAILY");
    INIT_FREQ("FR_DAY", FR_DAY, aliases);

    aliases = Py_BuildValue("(sssss)", "H", "HR", "HOUR", "HRLY", "HOURLY");
    INIT_FREQ("FR_HR", FR_HR, aliases);

    aliases = Py_BuildValue("(ssss)", "T", "MIN", "MINUTE", "MINUTELY");
    INIT_FREQ("FR_MIN", FR_MIN, aliases);

    aliases = Py_BuildValue("(ssss)", "S", "SEC", "SECOND", "SECONDLY");
    INIT_FREQ("FR_SEC", FR_SEC, aliases);

    aliases = Py_BuildValue("(ssss)", "U", "UND", "UNDEF", "UNDEFINED");
    INIT_FREQ("FR_UND", FR_UND, aliases);

    ADD_FREQ_CONSTANT("FR_ANN", FR_ANN);

    if(init_freq_group(12, 8, FR_ANN,
        month_names, ANN_prefixes, ANN_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    ADD_FREQ_CONSTANT("FR_WK", FR_WK);

    if(init_freq_group(7, 4, FR_WK,
                    day_names, WK_prefixes, WK_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    if(reverse_dict(freq_dict, freq_dict_rev) == INT_ERR_CODE) {
        return INT_ERR_CODE;
    }

    return 0;
}


/* take user specified frequency and convert to int representation
   of the frequency */
static int check_freq(PyObject *freq_spec) {

    if (PyInt_Check(freq_spec)) {
        return (int)PyInt_AsLong(freq_spec);
    } else if (PyString_Check(freq_spec)) {
        char *freq_str, *freq_str_uc;
        PyObject *freq_val;

        freq_str = PyString_AsString(freq_spec);
        if((freq_str_uc = str_uppercase(freq_str)) == NULL) {return INT_ERR_CODE;}

        freq_val = PyDict_GetItemString(freq_dict_rev, freq_str_uc);

        free(freq_str_uc);

        if (freq_val == NULL) {
            PyErr_SetString(PyExc_ValueError, "invalid frequency specification");
            return INT_ERR_CODE;
        } else {
            int ret_val = (int)PyInt_AsLong(freq_val);
            return ret_val;
        }
    } else if (freq_spec == Py_None) {
        return FR_UND;
    } else {
        int retval = (int)PyInt_AsLong(freq_spec);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError, "invalid frequency specification");
            return INT_ERR_CODE;
        } else { return retval; }
    }

}

static PyObject *
DateObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {

    DateObject *self;

    self = (DateObject*)type->tp_alloc(type, 0);
    if (self != NULL) {
        // initialize attributes that need initializing in here
        self->freq = FR_UND;
        self->value = -1;
    }

    return (PyObject *)self;
}

/* for use in C code */
static DateObject *
DateObject_New() {
    PyObject *dummy;
    return (DateObject*)DateObject_new(&DateType, dummy, dummy);
}

#define INIT_ERR(errortype, errmsg) PyErr_SetString(errortype,errmsg);return -1

static int
DateObject_init(DateObject *self, PyObject *args, PyObject *kwds) {

    PyObject *freq=NULL, *value=NULL, *datetime=NULL, *string=NULL;
    char *INSUFFICIENT_MSG = "insufficient parameters to initialize Date";

    int def_info=INT_ERR_CODE;

    int year=def_info, month=def_info, day=def_info, quarter=def_info,
        hour=def_info, minute=def_info, second=def_info;

    int free_dt=0;

    static char *kwlist[] = {"freq", "value", "string",
                             "year", "month", "day", "quarter",
                             "hour", "minute", "second",
                             "datetime", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|OOiiiiiiiO", kwlist,
                                      &freq, &value, &string,
                                      &year, &month, &day, &quarter,
                                      &hour, &minute, &second,
                                      &datetime)) return -1;

    if (PyObject_HasAttrString(freq, "freq")) {
        PyObject *freq_attr = PyObject_GetAttrString(freq, "freq");
        self->freq = PyInt_AS_LONG(freq_attr);
        Py_DECREF(freq_attr);
    } else {
        if((self->freq = check_freq(freq)) == INT_ERR_CODE) return -1;
    }

    if ((value && PyString_Check(value)) || string) {

        PyObject *string_arg = PyTuple_New(1);
        int freq_group = get_freq_group(self->freq);

        free_dt = 1;

        if (!string) {
            string = value;
        }

        PyTuple_SET_ITEM(string_arg, 0, string);
        Py_INCREF(string);

        if (freq_group == FR_HR ||
            freq_group == FR_MIN ||
            freq_group == FR_SEC)
             { datetime = PyEval_CallObject(DateTimeFromString, string_arg); }
        else { datetime = PyEval_CallObject(DateFromString, string_arg); }

        Py_DECREF(string_arg);

        value = NULL;
    }

    if (value) {
        self->value = PyInt_AsLong(value);
    } else {

        int freq_group = get_freq_group(self->freq);

        if (datetime) {
            year=PyDateTime_GET_YEAR(datetime);
            month=PyDateTime_GET_MONTH(datetime);
            day=PyDateTime_GET_DAY(datetime);
            hour=PyDateTime_DATE_GET_HOUR(datetime);
            minute=PyDateTime_DATE_GET_MINUTE(datetime);
            second=PyDateTime_DATE_GET_SECOND(datetime);
        }

        if (!datetime) {

            // First, some basic checks.....
            if (year == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
            if (self->freq == FR_BUS ||
               self->freq == FR_DAY ||
               self->freq == FR_WK ||
               self->freq == FR_UND) {
                if (month == def_info || day == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }

                // if FR_BUS, check for week day

            } else if (self->freq == FR_MTH) {
                if (month == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
            } else if (self->freq == FR_QTR) {
                if (quarter == def_info && month == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                } else if (month != def_info) {
                    quarter = (int)asfreq_MtoQ(month, '-', &NULL_AF_INFO);
                }
            } else if (self->freq == FR_SEC) {
                if (month == def_info ||
                    day == def_info ||
                    second == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
                if (hour == def_info) {
                    hour = second/3600;
                    minute = (second % 3600)/60;
                    second = second % 60;
                } else if (minute == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
            } else if (self->freq == FR_MIN) {
                if (month == def_info ||
                    day == def_info ||
                    minute == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
                if (hour == def_info) {
                    hour = minute/60;
                    minute = minute % 60;
                }
            } else if (self->freq == FR_HR) {
                if (month == def_info ||
                    day == def_info ||
                    hour == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
            }

        }

        if (self->freq == FR_SEC) {
            long absdays, delta;
            absdays = absdate_from_ymd(year, month, day);
            delta = (absdays - HIGHFREQ_ORIG);
            self->value = (int)(delta*86400 + hour*3600 + minute*60 + second + 1);
        } else if (self->freq == FR_MIN) {
            long absdays, delta;
            absdays = absdate_from_ymd(year, month, day);
            delta = (absdays - HIGHFREQ_ORIG);
            self->value = (int)(delta*1440 + hour*60 + minute + 1);
        } else if (self->freq == FR_HR) {
            long absdays, delta;
            absdays = absdate_from_ymd(year, month, day);
            delta = (absdays - HIGHFREQ_ORIG);
            self->value = (int)(delta*24 + hour + 1);
        } else if (self->freq == FR_DAY) {
            self->value = (int)absdate_from_ymd(year, month, day);
        } else if (self->freq == FR_UND) {
            self->value = (int)absdate_from_ymd(year, month, day);
        } else if (self->freq == FR_BUS) {
            long weeks, days;
            days = absdate_from_ymd(year, month, day);
            weeks = days/7;
            self->value = (int)(days - weeks*2);
        } else if (freq_group == FR_WK) {
            int adj_ordinal;
            int ordinal = (int)absdate_from_ymd(year, month, day);
            int day_adj = (7 - (self->freq - FR_WK)) % 7;

            adj_ordinal = ordinal + ((7 - day_adj) - ordinal % 7) % 7;
            self->value = adj_ordinal/7;
        } else if (self->freq == FR_MTH) {
            self->value = (year-1)*12 + month;
        } else if (self->freq == FR_QTR) {
            self->value = (year-1)*4 + quarter;
        } else if (freq_group == FR_ANN) {
            self->value = year;
        }

    }

    if (free_dt) { Py_DECREF(datetime); }

    return 0;
}

static PyMemberDef DateObject_members[] = {
    {"freq", T_INT, offsetof(DateObject, freq), 0,
     "frequency"},
    {"value", T_INT, offsetof(DateObject, value), 0,
     "integer representation of the Date"},
    {NULL}  /* Sentinel */
};

static char DateObject_toordinal_doc[] =
"Return the proleptic Gregorian ordinal of the date, where January 1 of\n"
"year 1 has ordinal 1";
static PyObject *
DateObject_toordinal(DateObject* self)
{
    if (self->freq == FR_DAY) {
        return PyInt_FromLong(self->value);
    } else {
        long (*toDaily)(long, char, struct asfreq_info*) = NULL;
        struct asfreq_info af_info;

        toDaily = get_asfreq_func(self->freq, FR_DAY, 0);
        get_asfreq_info(self->freq, FR_DAY, &af_info);

        return PyInt_FromLong(toDaily(self->value, 'A', &af_info));
    }
}

static char DateObject_asfreq_doc[] =
"Returns a date converted to a specified frequency.\n\n"
":Parameters:\n"
"   - freq : string/int\n"
"       Frequency to convert the Date to. Accepts any valid frequency\n"
"       specification (string or integer)\n"
"   - relation :string *['After']*\n"
"       Applies only when converting a lower frequency Date to a higher\n"
"       frequency Date, or when converting a weekend Date to a business\n"
"       frequency Date. Valid values are 'before', 'after', 'b', and 'a'.";
static PyObject *
DateObject_asfreq(DateObject *self, PyObject *args, PyObject *kwds)
{

    PyObject *freq=NULL;
    char *relation_raw=NULL;
    char *relation_uc;
    char relation;
    int invalid_relation=0;
    int toFreq;
    int result_val;
    DateObject *result = DateObject_New();

    static char *kwlist[] = {"freq", "relation", NULL};

    long (*asfreq_func)(long, char, struct asfreq_info*) = NULL;
    struct asfreq_info af_info;

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|s", kwlist,
                                      &freq, &relation_raw)) return NULL;

    if(relation_raw) {
        if (strlen(relation_raw) > 0) {
            if((relation_uc = str_uppercase(relation_raw)) == NULL)
            {return PyErr_NoMemory();}

            if (strcmp(relation_uc, "BEFORE") == 0 ||
                strcmp(relation_uc, "B") == 0 ||
                strcmp(relation_uc, "AFTER") == 0 ||
                strcmp(relation_uc, "A") == 0) {
                 relation = relation_uc[0];
            } else { invalid_relation=1; }
        } else {
            invalid_relation=1;
        }

        if (invalid_relation) {
            PyErr_SetString(PyExc_ValueError,"Invalid relation specification");
            return NULL;
        }
    } else {
        relation = 'A';
    }

    if ((toFreq = check_freq(freq)) == INT_ERR_CODE) return NULL;

    get_asfreq_info(self->freq, toFreq, &af_info);
    asfreq_func = get_asfreq_func(self->freq, toFreq, 0);

    result_val = asfreq_func(self->value, relation, &af_info);

    result->freq = toFreq;
    result->value = result_val;

    return (PyObject*)result;

}

static char DateObject_strfmt_doc[] =
"Returns string representation of Date object according to format specified.\n\n"
":Parameters:\n"
"   - fmt : string\n"
"       Formatting string. Uses the same directives as in the time.strftime\n"
"       function in the standard Python time module. In addition, a %q directive\n"
"       directive is recognized which represents the 'quarter' of the date";
static PyObject *
DateObject_strfmt(DateObject *self, PyObject *args)
{

    char *orig_fmt_str, *fmt_str;
    char *result;

    char extra_fmts[1][2][10] = {{"%q", "^`XZ`^"}};
    int extra_fmts_found[1] = {0};
    int extra_fmts_found_one = 0;
    struct tm c_date;
    struct date_info tempDate;
    long absdate;
    double abstime;
    int i, result_len, special_found=0;
    PyObject *py_result;

    long (*toDaily)(long, char, struct asfreq_info*) = NULL;
    struct asfreq_info af_info;

    if (!PyArg_ParseTuple(args, "s:strfmt(fmt)", &orig_fmt_str)) return NULL;

    toDaily = get_asfreq_func(self->freq, FR_DAY, 0);
    get_asfreq_info(self->freq, FR_DAY, &af_info);

    absdate = toDaily(self->value, 'A', &af_info);
    abstime = getAbsTime(self->freq, absdate, self->value);

    if(dInfoCalc_SetFromAbsDateTime(&tempDate, absdate, abstime,
                                    GREGORIAN_CALENDAR)) return NULL;

    // populate standard C date struct with info from our date_info struct
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
    if ((result = malloc(result_len * sizeof(char))) == NULL) {return PyErr_NoMemory();}

    fmt_str = orig_fmt_str;

    // replace any special format characters with their place holder
    for(i=0; i < 1; i++) {
        char *special_loc;
        if ((special_loc = strstr(fmt_str,extra_fmts[i][0])) != NULL) {
            char *tmp_str = fmt_str;
            fmt_str = str_replace(fmt_str, extra_fmts[i][0],
                                           extra_fmts[i][1]);
            /* only free the previous loop value if this is not the first
               special format string found */
            if (extra_fmts_found_one) { free(tmp_str); }

            if (fmt_str == NULL) {return NULL;}

            extra_fmts_found[i] = 1;
            extra_fmts_found_one = 1;
        }
    }

    strftime(result, result_len, fmt_str, &c_date);
    if (extra_fmts_found_one) { free(fmt_str); }

    // replace any place holders with the appropriate value
    for(i=0; i < 1; i++) {
        if (extra_fmts_found[i]) {
            char *tmp_str = result;
            char *extra_str;

            if(strcmp(extra_fmts[i][0], "%q") == 0) {
                if ((extra_str = malloc(2 * sizeof(char))) == NULL) {
                    free(tmp_str);
                    return PyErr_NoMemory();
                }

                sprintf(extra_str, "%i", tempDate.quarter);
            } else {
                PyErr_SetString(PyExc_RuntimeError,"Unrecogized fmt string");
                return NULL;
            }

            result = str_replace(result, extra_fmts[i][1], extra_str);
            free(tmp_str);
            if (result == NULL) { return NULL; }
        }
    }

    py_result = PyString_FromString(result);
    free(result);

    return py_result;
}

static PyObject *
DateObject___str__(DateObject* self)
{

    int freq_group = get_freq_group(self->freq);
    PyObject *string_arg, *retval;

    if (freq_group == FR_ANN) { string_arg = Py_BuildValue("(s)", "%Y"); }
    else if (freq_group == FR_QTR) { string_arg = Py_BuildValue("(s)", "%YQ%q"); }
    else if (freq_group == FR_MTH) { string_arg = Py_BuildValue("(s)", "%b-%Y"); }
    else if (freq_group == FR_DAY ||
             freq_group == FR_BUS ||
             freq_group == FR_WK ||
             freq_group == FR_UND) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y"); }
    else if (freq_group == FR_HR) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:00"); }
    else if (freq_group == FR_MIN) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M"); }
    else if (freq_group == FR_SEC) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M:%S"); }

    if (string_arg == NULL) { return NULL; }

    retval = DateObject_strfmt(self, string_arg);
    Py_DECREF(string_arg);

    return retval;
}

static PyObject *
DateObject_freqstr(DateObject *self, void *closure) {
    PyObject *key = PyInt_FromLong(self->freq);
    PyObject *freq_aliases = PyDict_GetItem(freq_dict, key);
    PyObject *main_alias = PyTuple_GET_ITEM(freq_aliases, 0);
    Py_DECREF(key);
    Py_INCREF(main_alias);
    return main_alias;
}


static PyObject *
DateObject___repr__(DateObject* self)
{
    PyObject *py_str_rep, *py_freqstr, *py_repr;
    char *str_rep, *freqstr, *repr;
    int repr_len;

    py_str_rep = DateObject___str__(self);
    if (py_str_rep == NULL) { return NULL; }

    py_freqstr = DateObject_freqstr(self, NULL);

    str_rep = PyString_AsString(py_str_rep);
    freqstr = PyString_AsString(py_freqstr);

    repr_len = strlen(str_rep) + strlen(freqstr) + 6;

    if((repr = malloc((repr_len + 1) * sizeof(char))) == NULL)
    { return PyErr_NoMemory(); }

    strcpy(repr, "<");
    strcat(repr, freqstr);
    strcat(repr, " : ");
    strcat(repr, str_rep);
    strcat(repr, ">");

    py_repr = PyString_FromString(repr);

    Py_DECREF(py_str_rep);
    Py_DECREF(py_freqstr);

    free(repr);

    return py_repr;
}

/******************************
   These methods seem rather useless. May or may not implement them.
fromordinal(self, ordinal):
    return Date(self.freq, datetime=dt.datetime.fromordinal(ordinal))
tostring(self):
    return str(self)
toobject(self):
    return self
isvalid(self):
    return True
*******************************/


static DateObject *
DateObject_FromFreqAndValue(int freq, int value) {

    DateObject *result = DateObject_New();

    PyObject *args = PyTuple_New(0);
    PyObject *kw = PyDict_New();
    PyObject *py_freq = PyInt_FromLong(freq);
    PyObject *py_value = PyInt_FromLong(value);

    PyDict_SetItemString(kw, "freq", py_freq);
    PyDict_SetItemString(kw, "value", py_value);

    Py_DECREF(py_freq);
    Py_DECREF(py_value);

    DateObject_init(result, args, kw);

    Py_DECREF(args);
    Py_DECREF(kw);

    return result;
}

static PyObject *
DateObject_date_plus_int(PyObject *date, PyObject *pyint) {
    DateObject *dateobj = (DateObject*)date;
    if (DateObject_Check(pyint)) {
        PyErr_SetString(PyExc_TypeError, "Cannot add two Date objects");
        return NULL;
    }

    return (PyObject*)DateObject_FromFreqAndValue(dateobj->freq, PyInt_AsLong(pyint) + dateobj->value);
}

static PyObject *
DateObject___add__(PyObject *left, PyObject *right)
{
    if (DateObject_Check(left)) {
        return DateObject_date_plus_int(left, right);
    } else {
        return DateObject_date_plus_int(right, left);
    }
}

static PyObject *
DateObject___subtract__(PyObject *left, PyObject *right)
{
    int result;
    DateObject *dleft;
    if (!DateObject_Check(left)) {
        PyErr_SetString(PyExc_ValueError, "Cannot subtract Date from non-Date value");
        return NULL;
    }

    dleft = (DateObject*)left;

    if (DateObject_Check(right)) {
        DateObject *dright = (DateObject*)right;
        if (dleft->freq != dright->freq) {
            PyErr_SetString(PyExc_ValueError, "Cannot subtract Dates with different frequency");
            return NULL;
        }
        result = dleft->value - dright->value;
        return PyInt_FromLong(result);
    } else {
        result = dleft->value - PyInt_AsLong(right);
        return (PyObject*)DateObject_FromFreqAndValue(dleft->freq, result);
    }
}

static int
DateObject___compare__(DateObject * obj1, DateObject * obj2)
{
    if (obj1->freq != obj2->freq) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot compare dates with different frequency");
        return -1;
    }

    if (obj1->value < obj2->value) return -1;
    if (obj1->value > obj2->value) return 1;
    if (obj1->value == obj2->value) return 0;

}

static long
DateObject___hash__(DateObject *self)
{
    register int freq_group = get_freq_group(self->freq);

    /* within a given frequency, hash values are guaranteed to be unique
       for different dates. For different frequencies, we make a reasonable
       effort to ensure hash values will be unique, but it is not guaranteed */
    if (freq_group == FR_BUS) {
        return self->value + 10000000;
    } else if (freq_group == FR_WK) {
        return self->value + 100000000;
    } else { return self->value; }
}

static PyObject *
DateObject___int__(DateObject *self)
{
    return PyInt_FromLong(self->value);
}

static PyObject *
DateObject___float__(DateObject *self)
{
    return PyFloat_FromDouble((double)(self->value));
}

/***************************************************
           ====== Date Properties ======
****************************************************/

// helper function for date property funcs
static int
DateObject_set_date_info(DateObject *self, struct date_info *dinfo) {
    PyObject *daily_obj = DateObject_toordinal(self);
    long absdate = PyInt_AsLong(daily_obj);

    Py_DECREF(daily_obj);

    if(dInfoCalc_SetFromAbsDate(dinfo, absdate,
                                GREGORIAN_CALENDAR)) return -1;

    return 0;
}

// helper function for date property funcs
static int
DateObject_set_date_info_wtime(DateObject *self, struct date_info *dinfo) {
    PyObject *daily_obj = DateObject_toordinal(self);
    long absdate = PyInt_AsLong(daily_obj);
    double abstime;

    Py_DECREF(daily_obj);

    abstime = getAbsTime(self->freq, absdate, self->value);

    if(dInfoCalc_SetFromAbsDateTime(dinfo, absdate, abstime,
                                    GREGORIAN_CALENDAR)) return -1;

    return 0;
}

static PyObject *
DateObject_year(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.year);
}

static PyObject *
DateObject_quarter(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.quarter);
}

static PyObject *
DateObject_month(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.month);
}

static PyObject *
DateObject_day(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day);
}

static PyObject *
DateObject_day_of_week(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day_of_week);
}

static PyObject *
DateObject_day_of_year(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day_of_year);
}

static PyObject *
DateObject_week(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dInfoCalc_ISOWeek(&dinfo));
}

static PyObject *
DateObject_hour(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.hour);
}

static PyObject *
DateObject_minute(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.minute);
}

static PyObject *
DateObject_second(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong((int)dinfo.second);
}

static PyObject *
DateObject_datetime(DateObject *self, void *closure) {
    PyObject *datetime;
    struct date_info dinfo;
    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    datetime = PyDateTime_FromDateAndTime(dinfo.year, dinfo.month,
                                          dinfo.day, dinfo.hour,
                                          dinfo.minute, (int)dinfo.second, 0);
    return datetime;
}

static int
DateObject_ReadOnlyErr(DateObject *self, PyObject *value, void *closure) {
   PyErr_SetString(PyExc_AttributeError, "Cannot set read-only property");
   return -1;
}

static PyGetSetDef DateObject_getseters[] = {
    {"year", (getter)DateObject_year, (setter)DateObject_ReadOnlyErr,
            "Returns the year.", NULL},
    {"quarter", (getter)DateObject_quarter, (setter)DateObject_ReadOnlyErr,
            "Returns the quarter.", NULL},
    {"month", (getter)DateObject_month, (setter)DateObject_ReadOnlyErr,
            "Returns the month.", NULL},
    {"week", (getter)DateObject_week, (setter)DateObject_ReadOnlyErr,
            "Returns the week.", NULL},
    {"day", (getter)DateObject_day, (setter)DateObject_ReadOnlyErr,
            "Returns the day of month.", NULL},
    {"day_of_week", (getter)DateObject_day_of_week, (setter)DateObject_ReadOnlyErr,
            "Returns the day of week.", NULL},
    {"day_of_year", (getter)DateObject_day_of_year, (setter)DateObject_ReadOnlyErr,
            "Returns the day of year.", NULL},
    {"second", (getter)DateObject_second, (setter)DateObject_ReadOnlyErr,
            "Returns the second.", NULL},
    {"minute", (getter)DateObject_minute, (setter)DateObject_ReadOnlyErr,
            "Returns the minute.", NULL},
    {"hour", (getter)DateObject_hour, (setter)DateObject_ReadOnlyErr,
            "Returns the hour.", NULL},

    {"freqstr", (getter)DateObject_freqstr, (setter)DateObject_ReadOnlyErr,
            "Returns the string representation of frequency.", NULL},
    {"datetime", (getter)DateObject_datetime, (setter)DateObject_ReadOnlyErr,
            "Returns the Date object converted to standard python datetime object",
            NULL},

    {NULL}  /* Sentinel */
};


static PyNumberMethods DateObject_as_number = {
    (binaryfunc)DateObject___add__,      /* nb_add */
    (binaryfunc)DateObject___subtract__, /* nb_subtract */
    0,                                   /* nb_multiply */
    0,                                   /* nb_divide */
    0,                                   /* nb_remainder */
    0,                                   /* nb_divmod */
    0,                                   /* nb_power */
    0,                                   /* nb_negative */
    0,                                   /* nb_positive */
    0,                                   /* nb_absolute */
    0,                                   /* nb_nonzero */
    0,				                     /* nb_invert */
    0,	                                 /* nb_lshift */
    0,	                                 /* nb_rshift */
    0,	                                 /* nb_and */
    0,	                                 /* nb_xor */
    0,	                                 /* nb_or */
    0,		                             /* nb_coerce */
    (unaryfunc)DateObject___int__,		 /* nb_int */
    (unaryfunc)0,				         /* nb_long */
    (unaryfunc)DateObject___float__,	 /* nb_float */
    (unaryfunc)0,				         /* nb_oct */
    (unaryfunc)0,				         /* nb_hex */
};

static PyMethodDef DateObject_methods[] = {
    {"toordinal", (PyCFunction)DateObject_toordinal, METH_NOARGS,
     DateObject_toordinal_doc},
    {"strfmt", (PyCFunction)DateObject_strfmt, METH_VARARGS,
     DateObject_strfmt_doc},
    {"asfreq", (PyCFunction)DateObject_asfreq, METH_VARARGS | METH_KEYWORDS,
     DateObject_asfreq_doc},
    {NULL}  /* Sentinel */
};


static PyTypeObject DateType = {
    PyObject_HEAD_INIT(NULL)
    0,                               /* ob_size */
    "timeseries.Date",               /* tp_name */
    sizeof(DateObject),              /* tp_basicsize */
    0,                               /* tp_itemsize */
    (destructor)DateObject_dealloc,  /* tp_dealloc */
    0,                               /* tp_print */
    0,                               /* tp_getattr */
    0,                               /* tp_setattr */
    (cmpfunc)DateObject___compare__, /* tp_compare */
    (reprfunc)DateObject___repr__,   /* tp_repr */
    &DateObject_as_number,           /* tp_as_number */
    0,                               /* tp_as_sequence */
    0,                               /* tp_as_mapping */
    (hashfunc)DateObject___hash__,   /* tp_hash */
    0,                               /* tp_call*/
    (reprfunc)DateObject___str__,    /* tp_str */
    0,                               /* tp_getattro */
    0,                               /* tp_setattro */
    0,                               /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |             /* tp_flags */
    Py_TPFLAGS_CHECKTYPES |
    Py_TPFLAGS_BASETYPE,
    "Defines a Date object, as the combination of a date and a frequency.\n"
    "Several options are available to construct a Date object explicitly:\n\n"
    "- Give appropriate values to the `year`, `month`, `day`, `quarter`, `hours`,\n"
    "  `minutes`, `seconds` arguments.\n\n"
    "  >>> td.Date(freq='Q',year=2004,quarter=3)\n"
    "  >>> td.Date(freq='D',year=2001,month=1,day=1)\n\n"
    "- Use the `string` keyword. This method uses a modified version of the\n"
    "  mx.DateTime parser submodule. More information is available in its\n"
    "  documentation.\n\n"
    "  >>> ts.Date('D', '2007-01-01')\n\n"
    "- Use the `datetime` keyword with an existing datetime.datetime object.\n\n"
    "  >>> td.Date('D', datetime=datetime.datetime.now())",  /* tp_doc */
    0,                               /* tp_traverse */
    0,                               /* tp_clear */
    0,                               /* tp_richcompare */
    0,                               /* tp_weaklistoffset */
    0,                               /* tp_iter */
    0,                               /* tp_iternext */
    DateObject_methods,              /* tp_methods */
    DateObject_members,              /* tp_members */
    DateObject_getseters,            /* tp_getset */
    0,                               /* tp_base */
    0,                               /* tp_dict */
    0,                               /* tp_descr_get */
    0,                               /* tp_descr_set */
    0,                               /* tp_dictoffset */
    (initproc)DateObject_init,       /* tp_init */
    0,                               /* tp_alloc */
    DateObject_new,                  /* tp_new */
};


///////////////////////////////////////////////////////////////////////

static char cseries_check_freq_doc[] =
"translate user specified frequency into frequency constant";
static PyObject *
cseries_check_freq(PyObject *self, PyObject *args) {

    PyObject *freq;
    int freq_val;

    if (!PyArg_ParseTuple(args, "O:check_freq(freq)", &freq)) return NULL;
    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) return NULL;

    return PyInt_FromLong(freq_val);
}

static char cseries_check_freq_str_doc[] =
"translate user specified frequency into standard string representation";
static PyObject *
cseries_check_freq_str(PyObject *self, PyObject *args) {

    PyObject *alias_tuple, *result, *freq_key;
    int freq_val;

    if ((freq_key = cseries_check_freq(self, args)) == NULL) return NULL;

    alias_tuple = PyDict_GetItem(freq_dict, freq_key);
    result = PyTuple_GET_ITEM(alias_tuple, 0);

    Py_INCREF(result);

    Py_DECREF(freq_key);

    return result;
}

static char cseries_get_freq_group_doc[] =
"translate user specified frequency into frequency group constant";
static PyObject *
cseries_get_freq_group(PyObject *self, PyObject *args) {

    PyObject *freq;
    int freq_val;

    if (!PyArg_ParseTuple(args, "O:get_freq_group(freq)", &freq)) return NULL;
    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) return NULL;

    return PyInt_FromLong(get_freq_group(freq_val));
}

static char cseries_thisday_doc[] =
"Returns today's date, at the given frequency\n\n"
":Parameters:\n"
"   - freq : string/int\n"
"       Frequency to convert the Date to. Accepts any valid frequency\n"
"       specification (string or integer)\n";
static PyObject *
cseries_thisday(PyObject *self, PyObject *args) {

    PyObject *freq, *init_args, *init_kwargs;
    time_t rawtime;
    struct tm *timeinfo;
    int freq_val;

    DateObject *secondly_date;

    if (!PyArg_ParseTuple(args, "O:thisday(freq)", &freq)) return NULL;

    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) return NULL;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    init_args = PyTuple_New(0);
    init_kwargs = PyDict_New();

    DICT_SETINT_STRKEY(init_kwargs, "freq", FR_SEC);
    DICT_SETINT_STRKEY(init_kwargs, "year", timeinfo->tm_year+1900);
    DICT_SETINT_STRKEY(init_kwargs, "month", timeinfo->tm_mon+1);
    DICT_SETINT_STRKEY(init_kwargs, "day", timeinfo->tm_mday);
    DICT_SETINT_STRKEY(init_kwargs, "hour", timeinfo->tm_hour);
    DICT_SETINT_STRKEY(init_kwargs, "minute", timeinfo->tm_min);
    DICT_SETINT_STRKEY(init_kwargs, "second", timeinfo->tm_sec);

    secondly_date = DateObject_New();
    DateObject_init(secondly_date, init_args, init_kwargs);

    Py_DECREF(init_args);
    Py_DECREF(init_kwargs);

    if (freq_val != FR_SEC) {
        DateObject *result = DateObject_New();

        long (*asfreq_func)(long, char, struct asfreq_info*) = NULL;
        struct asfreq_info af_info;

        int date_val;

        get_asfreq_info(FR_SEC, freq_val, &af_info);
        asfreq_func = get_asfreq_func(FR_SEC, freq_val, 0);

        date_val = asfreq_func(secondly_date->value, 'A', &af_info);

        Py_DECREF(secondly_date);

        result->freq = freq_val;
        result->value = date_val;

        return (PyObject*)result;

    } else { return (PyObject*)secondly_date; }
}


static char TimeSeries_convert_doc[] = "";
static PyObject *
TimeSeries_convert(PyObject *self, PyObject *args)
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
    PyObject *fromFreq_arg, *toFreq_arg;
    int fromFreq, toFreq;
    char relation;
    struct asfreq_info af_info;

    PyObject *val, *valMask;

    long (*asfreq_main)(long, char, struct asfreq_info*) = NULL;
    long (*asfreq_endpoints)(long, char, struct asfreq_info*) = NULL;
    long (*asfreq_reverse)(long, char, struct asfreq_info*) = NULL;

    returnVal = PyDict_New();

    if (!PyArg_ParseTuple(args,
        "OOOslO:convert(array, fromfreq, tofreq, position, startIndex, mask)",
        &array, &fromFreq_arg, &toFreq_arg,
        &position, &startIndex, &mask)) return NULL;

    if((fromFreq = check_freq(fromFreq_arg)) == INT_ERR_CODE) return NULL;
    if((toFreq = check_freq(toFreq_arg)) == INT_ERR_CODE) return NULL;

    if (toFreq == fromFreq)
    {
        PyObject *sidx;
        newArray = (PyArrayObject *)PyArray_Copy(array);
        newMask = (PyArrayObject *)PyArray_Copy(mask);
        sidx = PyInt_FromLong(startIndex);

        PyDict_SetItemString(returnVal, "values", (PyObject*)newArray);
        PyDict_SetItemString(returnVal, "mask", (PyObject*)newMask);
        PyDict_SetItemString(returnVal, "startindex", sidx);

        Py_DECREF(newArray);
        Py_DECREF(newMask);
        Py_DECREF(sidx);

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
    CHECK_ASFREQ(newStartTemp = asfreq_main(startIndex, 'B', &af_info));
    if (newStartTemp < 1) {
        CHECK_ASFREQ(newStart = asfreq_endpoints(startIndex, 'A', &af_info));
    }
    else { newStart = newStartTemp; }

    //convert end index to new frequency
    CHECK_ASFREQ(newEndTemp = asfreq_main(startIndex+array->dimensions[0]-1, 'A', &af_info));
    if (newEndTemp < 1) {
        CHECK_ASFREQ(newEnd = asfreq_endpoints(startIndex+array->dimensions[0]-1, 'B', &af_info));
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
        struct asfreq_info af_info_rev;

        get_asfreq_info(toFreq, fromFreq, &af_info_rev);
        asfreq_reverse = get_asfreq_func(toFreq, fromFreq, 0);

        CHECK_ASFREQ(tempval = asfreq_reverse(newStart, 'B', &af_info_rev));
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

        CHECK_ASFREQ(currIndex = asfreq_main(startIndex + i, relation, &af_info));

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

static char DateArray_asfreq_doc[] = "";
static PyObject *
DateArray_asfreq(PyObject *self, PyObject *args)
{
    PyArrayObject *fromDates, *toDates;
    PyArrayIterObject *iterFrom, *iterTo;
    PyObject *fromDateObj, *toDateObj;
    char *relation;
    int fromFreq, toFreq;
    long fromDate, toDate;
    long (*asfreq_main)(long, char, struct asfreq_info*) = NULL;
    struct asfreq_info af_info;

    if (!PyArg_ParseTuple(args,
                "Oiis:asfreq(fromDates, fromfreq, tofreq, relation)",
                &fromDates, &fromFreq, &toFreq, &relation)) return NULL;

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
        CHECK_ASFREQ(toDate = asfreq_main(fromDate, relation[0], &af_info));
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

static char DateArray_getDateInfo_doc[] = "";
static PyObject *
DateArray_getDateInfo(PyObject *self, PyObject *args)
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

    long (*toDaily)(long, char, struct asfreq_info*) = NULL;
    int (*getDateInfo)(struct date_info*) = NULL;
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
        CHECK_ASFREQ(absdate = toDaily(dateNum, 'B', &af_info));
        abstime = getAbsTime(freq, absdate, dateNum);

        if(dInfoCalc_SetFromAbsDateTime(&convDate, absdate, abstime, GREGORIAN_CALENDAR)) return NULL;
        dInfo = getDateInfo(&convDate);

        PyArray_SETITEM(newArray, iterResult->dataptr, PyInt_FromLong(dInfo));

        PyArray_ITER_NEXT(iterSource);
        PyArray_ITER_NEXT(iterResult);
    }

    Py_DECREF(iterSource);
    Py_DECREF(iterResult);

    return (PyObject *) newArray;
}

static PyMethodDef cseries_methods[] = {

    {"TS_convert", TimeSeries_convert, METH_VARARGS, TimeSeries_convert_doc},

    {"DA_asfreq", DateArray_asfreq, METH_VARARGS, DateArray_asfreq_doc},
    {"DA_getDateInfo", DateArray_getDateInfo, METH_VARARGS, DateArray_getDateInfo_doc},

    {"thisday", cseries_thisday, METH_VARARGS, cseries_thisday_doc},
    {"check_freq", cseries_check_freq, METH_VARARGS, cseries_check_freq_doc},
    {"check_freq_str", cseries_check_freq_str, METH_VARARGS, cseries_check_freq_str_doc},
    {"get_freq_group", cseries_get_freq_group, METH_VARARGS, cseries_get_freq_group_doc},

    {"set_callback_DateFromString", set_callback_DateFromString, METH_VARARGS,
     set_callback_DateFromString_doc},
    {"set_callback_DateTimeFromString", set_callback_DateTimeFromString, METH_VARARGS,
     set_callback_DateTimeFromString_doc},

    {NULL, NULL}
};

PyMODINIT_FUNC
initcseries(void)
{
    PyObject *m;

    if (PyType_Ready(&DateType) < 0) return;

    m = Py_InitModule3("cseries", cseries_methods, cseries_doc);
    if (m == NULL)
      return;

    import_array();
    PyDateTime_IMPORT;

    Py_INCREF(&DateType);
    PyModule_AddObject(m, "Date", (PyObject *)&DateType);

    if(build_freq_dict(m) == INT_ERR_CODE) {
    	PyErr_SetString(					\
    		PyExc_ImportError,				\
    		"initialization of module timeseries.cseries failed");
        return;
    };

    PyModule_AddObject(m, "freq_dict", freq_dict);
    PyModule_AddObject(m, "freq_dict_rev", freq_dict_rev);
    PyModule_AddObject(m, "freq_constants", freq_constants);
}