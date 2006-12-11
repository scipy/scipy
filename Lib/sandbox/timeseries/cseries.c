#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <string.h>
#include "mxDateTime.h"
#include "arrayobject.h"

static char cseries_doc[] = "Speed sensitive time series operations";

///////////////////////////////////////////////////////////////////////


static //PyArrayObject *
setArrayItem(PyArrayObject **theArray, long index, PyObject *newVal)
{
	char *setptr;

	if (index >= 0)
	{
		//set value in array
		setptr = (*theArray)->data + (index) * (*theArray)->strides[0];
		PyArray_SETITEM(*theArray,setptr,newVal);
	}

	//return theArray;
}

static int
freqVal(char freq)
{
	switch(freq)
	{
		case 'A':
			//annual
			return 1;
		case 'Q':
			//quarterly
			return 2;
		case 'M':
			//monthly
			return 3;
		case 'B':
			//business
			return 4;
		case 'D':
			//daily
			return 5;
		default:
			return 0;
	}
}

static long
toDaily(long fromDate, char fromFreq)
{
    long absdate, origin, secondorigin;
    int y,m,d;

	mxDateTimeObject *theDate;

    origin = 675333;
    secondorigin = 722814;

	//convert fromDate to days since (0 AD - 1 day)
    switch(fromFreq)
    {
        case 'D':
            absdate = fromDate;
            break;
        case 'B':
            absdate = (fromDate/5)*7 + fromDate%5;
            break;
        case 'M':
			y = fromDate/12 + 1;
			m = fromDate%12;

			if (m == 0)
			{
				m = 12;
				y--;
			}
			d=1;
			break;
        case 'Q':
        	y = fromDate/4 + 1;
        	m = (fromDate%4) * 3 - 2;

			if (m < 1)
			{
				m += 12;
				y--;
			}
			else if (m == 12)
			{
				m = 1;
				y++;
			}
			d=1;
			break;
        case 'A':
        	y = fromDate-1;
        	m = 1;
        	d = 1;
        	break;
        default:
            return -1;
    }

	if (freqVal(fromFreq) < 4)
	{
		//switch to years from 0 for mxDateTime
		y+= 1849;

		theDate = (mxDateTimeObject *)mxDateTime.DateTime_FromDateAndTime(y,m,d,0,0,0);
		absdate = (long)(theDate->absdate);
	}
	else
	{
		//days from 0 for mxDateTime
		absdate += origin;
	}

	return absdate;

}


static long
getDateInfo_sub(long dateNum, char freq, char info) {

	long monthNum;
	mxDateTimeObject *convDate;
	convDate = (mxDateTimeObject *)mxDateTime.DateTime_FromAbsDateAndTime(toDaily(dateNum,freq),0);

    switch(info)
    {
        case 'Y': //year

			return (long)(convDate->year);

        case 'Q': //quarter
			monthNum = (long)(convDate->month);
			return ((monthNum-1)/3)+1;

        case 'M': //month
			return (long)(convDate->month);

        case 'D': //day
			return (long)(convDate->day);

        case 'W': //day of week
			return (long)(convDate->day_of_week);
        default:
            return -1;
    }
}


static char cseries_getDateInfo_doc[] = "";
static PyObject *
cseries_getDateInfo(PyObject *self, PyObject *args)
{
	char *freq;
	char *info;

    PyArrayObject *array;
    PyArrayObject *tempArray;
    PyArrayObject *newArray;

    char *getptr;
    PyObject *val;
    long i, lngVal, dInfo, dim;

	if (!PyArg_ParseTuple(args, "Oss:getDateInfo(array, freq, info)", &tempArray, &freq, &info)) return NULL;

    array = PyArray_GETCONTIGUOUS(tempArray);

	dim = array->dimensions[0];

	//initialize new array
    newArray = (PyArrayObject*)PyArray_SimpleNew(array->nd, &dim, array->descr->type_num);

    for (i = 0; i < array->dimensions[0]; i++)
    {
		getptr = array->data + i*array->strides[0];
		val = PyArray_GETITEM(array, getptr);
		lngVal = PyInt_AsLong(val);
		dInfo = getDateInfo_sub(lngVal, *freq, *info);

		setArrayItem(&newArray, i, PyInt_FromLong(dInfo));
	}

	return (PyObject *) newArray;

}


//fromDate is periods since Dec 31, 1849
static long
convert(long fromDate, char fromFreq, char toFreq, int notStartInd)
{
    long absdate, origin, secondorigin, secsInDay;
    long converted;
    int rem;
    int y,m,d;

	mxDateTimeObject *convDate;

    origin = 675333;
    secondorigin = 722814;
    secsInDay = 86400;

	absdate = toDaily(fromDate, fromFreq);

	convDate = (mxDateTimeObject *)mxDateTime.DateTime_FromAbsDateAndTime(absdate,0);

	//switch back to days and years since 1849 for pyTSA Date
	absdate -= origin;
	y = convDate->year - 1849;
	m = convDate->month;

	//convert convDate to appropriate # of periods according to toFreq
    switch(toFreq)
    {
        case 'D':
            converted = absdate;
            break;
        case 'B':
        	rem = absdate%7;
            if (rem > 4) //is weekend day
            {
				if (notStartInd == 1 && freqVal(fromFreq) > 4)
				{
					return -1;
				}
				else
				{
					d = convDate->day;
					d -= rem - 4;	//change to friday before weekend
					if (d < 1) d += 3; //if friday was prev. month, change to monday instead
					absdate = absdate - convDate->day + d;
					converted = (long)((absdate / 7 * 5.0) + absdate%7);
				}
			}
			else
			{
				converted = (long)((absdate / 7 * 5.0) + rem);
			}
            break;
        case 'M':
        	converted = (long)((y-1)*12 + m);
        	break;
       	case 'Q':
       		converted = (long)((y-1)*4 + ((m-1)/3) + 1);
       		break;
       	case 'A':
       		converted = (long)(y+1);
       		break;
        default:
            return -1;
    }

    return converted;
}


static long
expand(long oldSize, char fromFr, char toFr)
{
	long newSize;
	int fromFreq, toFreq;

	if (fromFr == toFr) return oldSize;

	fromFreq = freqVal(fromFr);
	toFreq = freqVal(toFr);
	if (fromFreq*toFreq == 0) return oldSize; //invalid frequency

	newSize = oldSize;

	while (toFreq > fromFreq)
	{
		if (fromFreq == 1)	//Annual
		{
			newSize *= 4; //quarters in year
			fromFreq++;
		}
		else if (fromFreq == 2) //Quarterly
		{
			newSize *= 3; //months in quarter
			fromFreq++;
		}
		else if (fromFreq == 3)	//Monthly
		{
			newSize *= 31; //max days in month
			fromFreq++;
		}
		else if (fromFreq == 4)	//Business
		{
			newSize *= 2; //max d days for each b days
			fromFreq++;
		}
	}


	return newSize;
}


///////////////////////////////////////////////////////////////////////
/*
OBSERVED

from lower freq to higher freq
----------------------

summed --	all values in period set as lower freq's value / # of values

rest --		all values in period set as lower freq's value

from higher freq to lower freq
----------------------
begin - 	lower freq's value set as first value in period
end - 		lower freq's value set as end value in period
summed -	lower freq's value set as sum of all values in period
averaged -	lower freq's value set as average of all values in period
high - 		lower freq's value set as largest value in period
low - 		lower freq's value set as smallest value in period

*/
///////////////////////////////////////////////////////////////////////

static void
adjValForObsSet(PyArrayObject *theArray, char obs, PyObject **newVal, PyObject **newValMask, PyObject *val, PyObject *valMask,  long curPerLen)
{
	double dblVal;
	long lngValMask, lngAllMasked;

	lngValMask = PyInt_AsLong(valMask);
	lngAllMasked = PyInt_AsLong(*newValMask);

	if (!lngValMask) {

		// if any value is not masked, then we shall not mask the aggregated result
		*newValMask = valMask;

		if (obs == 'B')
		{
			if (lngAllMasked) {
				*newVal = val;
			}
		}
		else if ( PyArray_ISFLOAT(theArray) && (obs=='S' || obs=='A') )
		{

			if (obs == 'S')
			{
				//observed is summed

				dblVal = PyFloat_AsDouble(*newVal);
				dblVal += PyFloat_AsDouble(val);
				*newVal = PyFloat_FromDouble(dblVal);
			}
			else
			{
				//observed is averaged

				dblVal = PyFloat_AsDouble(*newVal);
				dblVal *= (curPerLen-1);
				dblVal += PyFloat_AsDouble(val);
				dblVal /= curPerLen;
				*newVal = PyFloat_FromDouble(dblVal);
			}

		}
		else if ( PyArray_ISNUMBER(theArray) && (obs=='H' || obs=='L') )
		{

			if (obs == 'H')
			{
				//observed is high

				if (PyFloat_AsDouble(val) > PyFloat_AsDouble(*newVal)) *newVal = val;
			}
			else if (obs == 'L')
			{
				//observed is low

				if (PyFloat_AsDouble(val) < PyFloat_AsDouble(*newVal)) *newVal = val;
			}

		}
		else
		{
			//observed is not beginning and
			//val is string or (val is date and observed is summed/averaged)
			//or observed is end or not supported

			*newVal = val;
		}
	}

}


static char cseries_reindex_doc[] = "";
static PyObject *
cseries_reindex(PyObject *self, PyObject *args)
{
    PyArrayObject *array;
    PyArrayObject *tempArray;
    PyArrayObject *newArray;

    PyArrayObject *mask;
    PyArrayObject *tempMask;
    PyArrayObject *newMask;

    PyObject *returnVal = NULL;

    int notStartInd;
    long startIndex, newStart;
    long i, curPerInd, nextPerInd, prevIndex, curIndex;
    long dim;
    long curPerLen;
    long lngValMask;
    char *fromFreq, *toFreq, *observed;

    char *getptr;
    PyObject *val, *newVal;

    char *getptrMask;
    PyObject *valMask, *newValMask;

    int toFrVal, fromFrVal;

	returnVal = PyDict_New();

	if (!PyArg_ParseTuple(args, "OssslO:reindex(array, fromfreq, tofreq, observed, startIndex,mask)", &tempArray, &fromFreq, &toFreq, &observed, &startIndex, &tempMask)) return NULL;

    if (toFreq[0] == fromFreq[0])
    {

        PyDict_SetItemString(returnVal, "values", (PyObject*)tempArray);
        PyDict_SetItemString(returnVal, "mask", (PyObject*)tempMask);

        return returnVal;
    }

    array = PyArray_GETCONTIGUOUS(tempArray);
    mask = PyArray_GETCONTIGUOUS(tempMask);

	//expand size to fit new values if needed
	dim = expand(array->dimensions[0], fromFreq[0], toFreq[0]);

	//initialize new array
    newArray = (PyArrayObject*)PyArray_SimpleNew(array->nd, &dim, array->descr->type_num);
    newMask  = (PyArrayObject*)PyArray_SimpleNew(mask->nd, &dim, mask->descr->type_num);

    for (i = 0; i < dim; i++)
    {
		setArrayItem(&newArray, i, PyInt_FromLong(1));
		setArrayItem(&newMask, i, PyInt_FromLong(1));
	}

	//convert start index to new frequency
	notStartInd = 0;
    newStart = convert(startIndex, fromFreq[0], toFreq[0], notStartInd);

	//initialize prevIndex
	prevIndex = newStart - 1;

	notStartInd = 1;

	//set values in the new array
    for (i = 0; i < array->dimensions[0]; i++)
    {
		//find index for start of current period in new frequency
        curPerInd = convert(startIndex + i, fromFreq[0], toFreq[0], notStartInd);

		//get frequency numeric mapping
		fromFrVal = freqVal(fromFreq[0]);
		toFrVal = freqVal(toFreq[0]);

		//get value from old array
		getptr = array->data + i*array->strides[0];
		val = PyArray_GETITEM(array,getptr);

		//get the mask corresponding to the old value
		getptrMask = mask->data + i*mask->strides[0];
		valMask = PyArray_GETITEM(mask,getptrMask);

        if (fromFrVal < toFrVal)
        {
			//from lower freq to higher freq

			newVal = val;
			newValMask = valMask;

			//find index for start of next period in new frequency
			nextPerInd = convert(startIndex + i + 1, fromFreq[0], toFreq[0], notStartInd);

			//adjust for observed setting
			if (observed[0] == 'S' && PyArray_ISFLOAT(array) && !( (fromFrVal == 4 && toFrVal == 5) || (fromFrVal == 5 && toFrVal == 4) ) )
			{
				//summed

				//all values in period set as old array's value / # of values
				newVal = PyFloat_FromDouble( PyFloat_AsDouble(val) / (nextPerInd - curPerInd) );
			}

			//set each value in period
			for (curIndex = curPerInd; curIndex < nextPerInd; curIndex++)
			{
				setArrayItem(&newArray, curIndex-newStart, newVal);
				setArrayItem(&newMask, curIndex-newStart, newValMask);
			}
		}
		else
		{

			lngValMask = PyInt_AsLong(valMask);

			//from higher freq to lower freq

			if (curPerInd != prevIndex)
			{
				//starting new period in old array


				//set value in the new array
				setArrayItem(&newArray, prevIndex-newStart, newVal);
				setArrayItem(&newMask, prevIndex-newStart, newValMask);

				//reset period length
				curPerLen = 0;



				if (!lngValMask) {
					curPerLen++;
				}



				//store current index and value
				prevIndex = curPerInd;
				newVal = val;
				newValMask = valMask;

			}
			else
			{
				//still in same period



				if (!lngValMask) {
					curPerLen++;
				}

				//adjust new value according to observed setting
				adjValForObsSet(array, observed[0], &newVal, &newValMask, val, valMask, curPerLen);
			}

		}

    }

	//set value of last item in the new array
	setArrayItem(&newArray, curPerInd-newStart, newVal);
	setArrayItem(&newMask, curPerInd-newStart, newValMask);

	PyDict_SetItemString(returnVal, "values", (PyObject*)newArray);
	PyDict_SetItemString(returnVal, "mask", (PyObject*)newMask);

	return returnVal;

}


static char cseries_convert_doc[] = "";
static PyObject *
cseries_convert(PyObject *self, PyObject *args)
{
    long fromDate;
    char* fromFreq;
    char* toFreq;
    int notStartInd;

    if (!PyArg_ParseTuple(args, "lss:convert(fromDate, fromfreq, tofreq)", &fromDate, &fromFreq, &toFreq)) return NULL;

	//always want start of period (only matters when converting from lower freq to higher freq ie. m -> d)
	notStartInd = 0;

    return PyInt_FromLong(convert(fromDate, fromFreq[0], toFreq[0], notStartInd));
}


///////////////////////////////////////////////////////////////////////

static PyMethodDef cseries_methods[] = {
    {"reindex", cseries_reindex, METH_VARARGS, cseries_reindex_doc},
    {"convert", cseries_convert, METH_VARARGS, cseries_convert_doc},
	{"getDateInfo", cseries_getDateInfo, METH_VARARGS, cseries_getDateInfo_doc},
    {NULL, NULL}
};

PyMODINIT_FUNC
initcseries(void)
{
    Py_InitModule3("cseries", cseries_methods, cseries_doc);
    mxDateTime_ImportModuleAndAPI();
    import_array();
}