
#ifdef _MULTIARRAYMODULE

typedef struct {
        PyObject_HEAD
        npy_bool obval;
} PyBoolScalarObject;


static unsigned int PyArray_GetNDArrayCVersion (void);
static PyTypeObject PyBigArray_Type;
static PyTypeObject PyArray_Type;
static PyTypeObject PyArrayDescr_Type;
static PyTypeObject PyArrayFlags_Type;
static PyTypeObject PyArrayIter_Type;
static PyTypeObject PyArrayMapIter_Type;
static PyTypeObject PyArrayMultiIter_Type;
static int NPY_NUMUSERTYPES=0;
static PyTypeObject PyBoolArrType_Type;
static PyBoolScalarObject _PyArrayScalar_BoolValues[2];

static PyTypeObject PyGenericArrType_Type;
static PyTypeObject PyNumberArrType_Type;
static PyTypeObject PyIntegerArrType_Type;
static PyTypeObject PySignedIntegerArrType_Type;
static PyTypeObject PyUnsignedIntegerArrType_Type;
static PyTypeObject PyInexactArrType_Type;
static PyTypeObject PyFloatingArrType_Type;
static PyTypeObject PyComplexFloatingArrType_Type;
static PyTypeObject PyFlexibleArrType_Type;
static PyTypeObject PyCharacterArrType_Type;
static PyTypeObject PyByteArrType_Type;
static PyTypeObject PyShortArrType_Type;
static PyTypeObject PyIntArrType_Type;
static PyTypeObject PyLongArrType_Type;
static PyTypeObject PyLongLongArrType_Type;
static PyTypeObject PyUByteArrType_Type;
static PyTypeObject PyUShortArrType_Type;
static PyTypeObject PyUIntArrType_Type;
static PyTypeObject PyULongArrType_Type;
static PyTypeObject PyULongLongArrType_Type;
static PyTypeObject PyFloatArrType_Type;
static PyTypeObject PyDoubleArrType_Type;
static PyTypeObject PyLongDoubleArrType_Type;
static PyTypeObject PyCFloatArrType_Type;
static PyTypeObject PyCDoubleArrType_Type;
static PyTypeObject PyCLongDoubleArrType_Type;
static PyTypeObject PyObjectArrType_Type;
static PyTypeObject PyStringArrType_Type;
static PyTypeObject PyUnicodeArrType_Type;
static PyTypeObject PyVoidArrType_Type;
static int PyArray_SetNumericOps \
       (PyObject *);
static PyObject * PyArray_GetNumericOps \
       (void);
static int PyArray_INCREF \
       (PyArrayObject *);
static int PyArray_XDECREF \
       (PyArrayObject *);
static void PyArray_SetStringFunction \
       (PyObject *, int);
static PyArray_Descr * PyArray_DescrFromType \
       (int);
static PyObject * PyArray_TypeObjectFromType \
       (int);
static char * PyArray_Zero \
       (PyArrayObject *);
static char * PyArray_One \
       (PyArrayObject *);
static PyObject * PyArray_CastToType \
       (PyArrayObject *, PyArray_Descr *, int);
static int PyArray_CastTo \
       (PyArrayObject *, PyArrayObject *);
static int PyArray_CastAnyTo \
       (PyArrayObject *, PyArrayObject *);
static int PyArray_CanCastSafely \
       (int, int);
static npy_bool PyArray_CanCastTo \
       (PyArray_Descr *, PyArray_Descr *);
static int PyArray_ObjectType \
       (PyObject *, int);
static PyArray_Descr * PyArray_DescrFromObject \
       (PyObject *, PyArray_Descr *);
static PyArrayObject ** PyArray_ConvertToCommonType \
       (PyObject *, int *);
static PyArray_Descr * PyArray_DescrFromScalar \
       (PyObject *);
static PyArray_Descr * PyArray_DescrFromTypeObject \
       (PyObject *);
static npy_intp PyArray_Size \
       (PyObject *);
static PyObject * PyArray_Scalar \
       (void *, PyArray_Descr *, PyObject *);
static PyObject * PyArray_FromScalar \
       (PyObject *, PyArray_Descr *);
static void PyArray_ScalarAsCtype \
       (PyObject *, void *);
static int PyArray_CastScalarToCtype \
       (PyObject *, void *, PyArray_Descr *);
static int PyArray_CastScalarDirect \
       (PyObject *, PyArray_Descr *, void *, int);
static PyObject * PyArray_ScalarFromObject \
       (PyObject *);
static PyArray_VectorUnaryFunc * PyArray_GetCastFunc \
       (PyArray_Descr *, int);
static PyObject * PyArray_FromDims \
       (int, int *, int);
static PyObject * PyArray_FromDimsAndDataAndDescr \
       (int, int *, PyArray_Descr *, char *);
static PyObject * PyArray_FromAny \
       (PyObject *, PyArray_Descr *, int, int, int, PyObject *);
static PyObject * PyArray_EnsureArray \
       (PyObject *);
static PyObject * PyArray_EnsureAnyArray \
       (PyObject *);
static PyObject * PyArray_FromFile \
       (FILE *, PyArray_Descr *, npy_intp, char *);
static PyObject * PyArray_FromString \
       (char *, npy_intp, PyArray_Descr *, npy_intp, char *);
static PyObject * PyArray_FromBuffer \
       (PyObject *, PyArray_Descr *, npy_intp, npy_intp);
static PyObject * PyArray_FromIter \
       (PyObject *, PyArray_Descr *, npy_intp);
static PyObject * PyArray_Return \
       (PyArrayObject *);
static PyObject * PyArray_GetField \
       (PyArrayObject *, PyArray_Descr *, int);
static int PyArray_SetField \
       (PyArrayObject *, PyArray_Descr *, int, PyObject *);
static PyObject * PyArray_Byteswap \
       (PyArrayObject *, npy_bool);
static PyObject * PyArray_Resize \
       (PyArrayObject *, PyArray_Dims *, int, NPY_ORDER);
static int PyArray_MoveInto \
       (PyArrayObject *, PyArrayObject *);
static int PyArray_CopyInto \
       (PyArrayObject *, PyArrayObject *);
static int PyArray_CopyAnyInto \
       (PyArrayObject *, PyArrayObject *);
static int PyArray_CopyObject \
       (PyArrayObject *, PyObject *);
static PyObject * PyArray_NewCopy \
       (PyArrayObject *, NPY_ORDER);
static PyObject * PyArray_ToList \
       (PyArrayObject *);
static PyObject * PyArray_ToString \
       (PyArrayObject *, NPY_ORDER);
static int PyArray_ToFile \
       (PyArrayObject *, FILE *, char *, char *);
static int PyArray_Dump \
       (PyObject *, PyObject *, int);
static PyObject * PyArray_Dumps \
       (PyObject *, int);
static int PyArray_ValidType \
       (int);
static void PyArray_UpdateFlags \
       (PyArrayObject *, int);
static PyObject * PyArray_New \
       (PyTypeObject *, int, npy_intp *, int, npy_intp *, void *, int, int, PyObject *);
static PyObject * PyArray_NewFromDescr \
       (PyTypeObject *, PyArray_Descr *, int, npy_intp *, npy_intp *, void *, int, PyObject *);
static PyArray_Descr * PyArray_DescrNew \
       (PyArray_Descr *);
static PyArray_Descr * PyArray_DescrNewFromType \
       (int);
static double PyArray_GetPriority \
       (PyObject *, double);
static PyObject * PyArray_IterNew \
       (PyObject *);
static PyObject * PyArray_MultiIterNew \
       (int, ...);
static int PyArray_PyIntAsInt \
       (PyObject *);
static npy_intp PyArray_PyIntAsIntp \
       (PyObject *);
static int PyArray_Broadcast \
       (PyArrayMultiIterObject *);
static void PyArray_FillObjectArray \
       (PyArrayObject *, PyObject *);
static int PyArray_FillWithScalar \
       (PyArrayObject *, PyObject *);
static npy_bool PyArray_CheckStrides \
       (int, int, npy_intp, npy_intp, npy_intp *, npy_intp *);
static PyArray_Descr * PyArray_DescrNewByteorder \
       (PyArray_Descr *, char);
static PyObject * PyArray_IterAllButAxis \
       (PyObject *, int *);
static PyObject * PyArray_CheckFromAny \
       (PyObject *, PyArray_Descr *, int, int, int, PyObject *);
static PyObject * PyArray_FromArray \
       (PyArrayObject *, PyArray_Descr *, int);
static PyObject * PyArray_FromInterface \
       (PyObject *);
static PyObject * PyArray_FromStructInterface \
       (PyObject *);
static PyObject * PyArray_FromArrayAttr \
       (PyObject *, PyArray_Descr *, PyObject *);
static NPY_SCALARKIND PyArray_ScalarKind \
       (int, PyArrayObject **);
static int PyArray_CanCoerceScalar \
       (int, int, NPY_SCALARKIND);
static PyObject * PyArray_NewFlagsObject \
       (PyObject *);
static npy_bool PyArray_CanCastScalar \
       (PyTypeObject *, PyTypeObject *);
static int PyArray_CompareUCS4 \
       (npy_ucs4 *, npy_ucs4 *, register size_t);
static int PyArray_RemoveSmallest \
       (PyArrayMultiIterObject *);
static int PyArray_ElementStrides \
       (PyObject *);
static void PyArray_Item_INCREF \
       (char *, PyArray_Descr *);
static void PyArray_Item_XDECREF \
       (char *, PyArray_Descr *);
static PyObject * PyArray_FieldNames \
       (PyObject *);
static PyObject * PyArray_Transpose \
       (PyArrayObject *, PyArray_Dims *);
static PyObject * PyArray_TakeFrom \
       (PyArrayObject *, PyObject *, int, PyArrayObject *, NPY_CLIPMODE);
static PyObject * PyArray_PutTo \
       (PyArrayObject *, PyObject*, PyObject *, NPY_CLIPMODE);
static PyObject * PyArray_PutMask \
       (PyArrayObject *, PyObject*, PyObject*);
static PyObject * PyArray_Repeat \
       (PyArrayObject *, PyObject *, int);
static PyObject * PyArray_Choose \
       (PyArrayObject *, PyObject *, PyArrayObject *, NPY_CLIPMODE);
static int PyArray_Sort \
       (PyArrayObject *, int, NPY_SORTKIND);
static PyObject * PyArray_ArgSort \
       (PyArrayObject *, int, NPY_SORTKIND);
static PyObject * PyArray_SearchSorted \
       (PyArrayObject *, PyObject *, NPY_SEARCHSIDE);
static PyObject * PyArray_ArgMax \
       (PyArrayObject *, int, PyArrayObject *);
static PyObject * PyArray_ArgMin \
       (PyArrayObject *, int, PyArrayObject *);
static PyObject * PyArray_Reshape \
       (PyArrayObject *, PyObject *);
static PyObject * PyArray_Newshape \
       (PyArrayObject *, PyArray_Dims *, NPY_ORDER);
static PyObject * PyArray_Squeeze \
       (PyArrayObject *);
static PyObject * PyArray_View \
       (PyArrayObject *, PyArray_Descr *, PyTypeObject *);
static PyObject * PyArray_SwapAxes \
       (PyArrayObject *, int, int);
static PyObject * PyArray_Max \
       (PyArrayObject *, int, PyArrayObject *);
static PyObject * PyArray_Min \
       (PyArrayObject *, int, PyArrayObject *);
static PyObject * PyArray_Ptp \
       (PyArrayObject *, int, PyArrayObject *);
static PyObject * PyArray_Mean \
       (PyArrayObject *, int, int, PyArrayObject *);
static PyObject * PyArray_Trace \
       (PyArrayObject *, int, int, int, int, PyArrayObject *);
static PyObject * PyArray_Diagonal \
       (PyArrayObject *, int, int, int);
static PyObject * PyArray_Clip \
       (PyArrayObject *, PyObject *, PyObject *, PyArrayObject *);
static PyObject * PyArray_Conjugate \
       (PyArrayObject *, PyArrayObject *);
static PyObject * PyArray_Nonzero \
       (PyArrayObject *);
static PyObject * PyArray_Std \
       (PyArrayObject *, int, int, PyArrayObject *, int);
static PyObject * PyArray_Sum \
       (PyArrayObject *, int, int, PyArrayObject *);
static PyObject * PyArray_CumSum \
       (PyArrayObject *, int, int, PyArrayObject *);
static PyObject * PyArray_Prod \
       (PyArrayObject *, int, int, PyArrayObject *);
static PyObject * PyArray_CumProd \
       (PyArrayObject *, int, int, PyArrayObject *);
static PyObject * PyArray_All \
       (PyArrayObject *, int, PyArrayObject *);
static PyObject * PyArray_Any \
       (PyArrayObject *, int, PyArrayObject *);
static PyObject * PyArray_Compress \
       (PyArrayObject *, PyObject *, int, PyArrayObject *);
static PyObject * PyArray_Flatten \
       (PyArrayObject *, NPY_ORDER);
static PyObject * PyArray_Ravel \
       (PyArrayObject *, NPY_ORDER);
static npy_intp PyArray_MultiplyList \
       (register npy_intp *, register int);
static int PyArray_MultiplyIntList \
       (register int *, register int);
static void * PyArray_GetPtr \
       (PyArrayObject *, register npy_intp*);
static int PyArray_CompareLists \
       (npy_intp *, npy_intp *, int);
static int PyArray_AsCArray \
       (PyObject **, void *, npy_intp *, int, PyArray_Descr*);
static int PyArray_As1D \
       (PyObject **, char **, int *, int);
static int PyArray_As2D \
       (PyObject **, char ***, int *, int *, int);
static int PyArray_Free \
       (PyObject *, void *);
static int PyArray_Converter \
       (PyObject *, PyObject **);
static int PyArray_IntpFromSequence \
       (PyObject *, npy_intp *, int);
static PyObject * PyArray_Concatenate \
       (PyObject *, int);
static PyObject * PyArray_InnerProduct \
       (PyObject *, PyObject *);
static PyObject * PyArray_MatrixProduct \
       (PyObject *, PyObject *);
static PyObject * PyArray_CopyAndTranspose \
       (PyObject *);
static PyObject * PyArray_Correlate \
       (PyObject *, PyObject *, int);
static int PyArray_TypestrConvert \
       (int, int);
static int PyArray_DescrConverter \
       (PyObject *, PyArray_Descr **);
static int PyArray_DescrConverter2 \
       (PyObject *, PyArray_Descr **);
static int PyArray_IntpConverter \
       (PyObject *, PyArray_Dims *);
static int PyArray_BufferConverter \
       (PyObject *, PyArray_Chunk *);
static int PyArray_AxisConverter \
       (PyObject *, int *);
static int PyArray_BoolConverter \
       (PyObject *, npy_bool *);
static int PyArray_ByteorderConverter \
       (PyObject *, char *);
static int PyArray_OrderConverter \
       (PyObject *, NPY_ORDER *);
static unsigned char PyArray_EquivTypes \
       (PyArray_Descr *, PyArray_Descr *);
static PyObject * PyArray_Zeros \
       (int, npy_intp *, PyArray_Descr *, int);
static PyObject * PyArray_Empty \
       (int, npy_intp *, PyArray_Descr *, int);
static PyObject * PyArray_Where \
       (PyObject *, PyObject *, PyObject *);
static PyObject * PyArray_Arange \
       (double, double, double, int);
static PyObject * PyArray_ArangeObj \
       (PyObject *, PyObject *, PyObject *, PyArray_Descr *);
static int PyArray_SortkindConverter \
       (PyObject *, NPY_SORTKIND *);
static PyObject * PyArray_LexSort \
       (PyObject *, int);
static PyObject * PyArray_Round \
       (PyArrayObject *, int, PyArrayObject *);
static unsigned char PyArray_EquivTypenums \
       (int, int);
static int PyArray_RegisterDataType \
       (PyArray_Descr *);
static int PyArray_RegisterCastFunc \
       (PyArray_Descr *, int, PyArray_VectorUnaryFunc *);
static int PyArray_RegisterCanCast \
       (PyArray_Descr *, int, NPY_SCALARKIND);
static void PyArray_InitArrFuncs \
       (PyArray_ArrFuncs *);
static PyObject * PyArray_IntTupleFromIntp \
       (int, npy_intp *);
static int PyArray_TypeNumFromName \
       (char *);
static int PyArray_ClipmodeConverter \
       (PyObject *, NPY_CLIPMODE *);
static int PyArray_OutputConverter \
       (PyObject *, PyArrayObject **);
static PyObject * PyArray_BroadcastToShape \
       (PyObject *, npy_intp *, int);
static void _PyArray_SigintHandler \
       (int);
static void* _PyArray_GetSigintBuf \
       (void);
static int PyArray_DescrAlignConverter \
       (PyObject *, PyArray_Descr **);
static int PyArray_DescrAlignConverter2 \
       (PyObject *, PyArray_Descr **);
static int PyArray_SearchsideConverter \
       (PyObject *, void *);

#else

#if defined(PY_ARRAY_UNIQUE_SYMBOL)
#define PyArray_API PY_ARRAY_UNIQUE_SYMBOL
#endif

#if defined(NO_IMPORT) || defined(NO_IMPORT_ARRAY)
extern void **PyArray_API;
#else
#if defined(PY_ARRAY_UNIQUE_SYMBOL)
void **PyArray_API;
#else
static void **PyArray_API=NULL;
#endif
#endif

#define PyArray_GetNDArrayCVersion (*(unsigned int (*)(void)) PyArray_API[0])
#define PyBigArray_Type (*(PyTypeObject *)PyArray_API[1])
#define PyArray_Type (*(PyTypeObject *)PyArray_API[2])
#define PyArrayDescr_Type (*(PyTypeObject *)PyArray_API[3])
#define PyArrayFlags_Type (*(PyTypeObject *)PyArray_API[4])
#define PyArrayIter_Type (*(PyTypeObject *)PyArray_API[5])
#define PyArrayMultiIter_Type (*(PyTypeObject *)PyArray_API[6])
#define NPY_NUMUSERTYPES (*(int *)PyArray_API[7])
#define PyBoolArrType_Type (*(PyTypeObject *)PyArray_API[8])
#define _PyArrayScalar_BoolValues ((PyBoolScalarObject *)PyArray_API[9])

#define PyGenericArrType_Type (*(PyTypeObject *)PyArray_API[10])
#define PyNumberArrType_Type (*(PyTypeObject *)PyArray_API[11])
#define PyIntegerArrType_Type (*(PyTypeObject *)PyArray_API[12])
#define PySignedIntegerArrType_Type (*(PyTypeObject *)PyArray_API[13])
#define PyUnsignedIntegerArrType_Type (*(PyTypeObject *)PyArray_API[14])
#define PyInexactArrType_Type (*(PyTypeObject *)PyArray_API[15])
#define PyFloatingArrType_Type (*(PyTypeObject *)PyArray_API[16])
#define PyComplexFloatingArrType_Type (*(PyTypeObject *)PyArray_API[17])
#define PyFlexibleArrType_Type (*(PyTypeObject *)PyArray_API[18])
#define PyCharacterArrType_Type (*(PyTypeObject *)PyArray_API[19])
#define PyByteArrType_Type (*(PyTypeObject *)PyArray_API[20])
#define PyShortArrType_Type (*(PyTypeObject *)PyArray_API[21])
#define PyIntArrType_Type (*(PyTypeObject *)PyArray_API[22])
#define PyLongArrType_Type (*(PyTypeObject *)PyArray_API[23])
#define PyLongLongArrType_Type (*(PyTypeObject *)PyArray_API[24])
#define PyUByteArrType_Type (*(PyTypeObject *)PyArray_API[25])
#define PyUShortArrType_Type (*(PyTypeObject *)PyArray_API[26])
#define PyUIntArrType_Type (*(PyTypeObject *)PyArray_API[27])
#define PyULongArrType_Type (*(PyTypeObject *)PyArray_API[28])
#define PyULongLongArrType_Type (*(PyTypeObject *)PyArray_API[29])
#define PyFloatArrType_Type (*(PyTypeObject *)PyArray_API[30])
#define PyDoubleArrType_Type (*(PyTypeObject *)PyArray_API[31])
#define PyLongDoubleArrType_Type (*(PyTypeObject *)PyArray_API[32])
#define PyCFloatArrType_Type (*(PyTypeObject *)PyArray_API[33])
#define PyCDoubleArrType_Type (*(PyTypeObject *)PyArray_API[34])
#define PyCLongDoubleArrType_Type (*(PyTypeObject *)PyArray_API[35])
#define PyObjectArrType_Type (*(PyTypeObject *)PyArray_API[36])
#define PyStringArrType_Type (*(PyTypeObject *)PyArray_API[37])
#define PyUnicodeArrType_Type (*(PyTypeObject *)PyArray_API[38])
#define PyVoidArrType_Type (*(PyTypeObject *)PyArray_API[39])
#define PyArray_SetNumericOps \
        (*(int (*)(PyObject *)) \
         PyArray_API[40])
#define PyArray_GetNumericOps \
        (*(PyObject * (*)(void)) \
         PyArray_API[41])
#define PyArray_INCREF \
        (*(int (*)(PyArrayObject *)) \
         PyArray_API[42])
#define PyArray_XDECREF \
        (*(int (*)(PyArrayObject *)) \
         PyArray_API[43])
#define PyArray_SetStringFunction \
        (*(void (*)(PyObject *, int)) \
         PyArray_API[44])
#define PyArray_DescrFromType \
        (*(PyArray_Descr * (*)(int)) \
         PyArray_API[45])
#define PyArray_TypeObjectFromType \
        (*(PyObject * (*)(int)) \
         PyArray_API[46])
#define PyArray_Zero \
        (*(char * (*)(PyArrayObject *)) \
         PyArray_API[47])
#define PyArray_One \
        (*(char * (*)(PyArrayObject *)) \
         PyArray_API[48])
#define PyArray_CastToType \
        (*(PyObject * (*)(PyArrayObject *, PyArray_Descr *, int)) \
         PyArray_API[49])
#define PyArray_CastTo \
        (*(int (*)(PyArrayObject *, PyArrayObject *)) \
         PyArray_API[50])
#define PyArray_CastAnyTo \
        (*(int (*)(PyArrayObject *, PyArrayObject *)) \
         PyArray_API[51])
#define PyArray_CanCastSafely \
        (*(int (*)(int, int)) \
         PyArray_API[52])
#define PyArray_CanCastTo \
        (*(npy_bool (*)(PyArray_Descr *, PyArray_Descr *)) \
         PyArray_API[53])
#define PyArray_ObjectType \
        (*(int (*)(PyObject *, int)) \
         PyArray_API[54])
#define PyArray_DescrFromObject \
        (*(PyArray_Descr * (*)(PyObject *, PyArray_Descr *)) \
         PyArray_API[55])
#define PyArray_ConvertToCommonType \
        (*(PyArrayObject ** (*)(PyObject *, int *)) \
         PyArray_API[56])
#define PyArray_DescrFromScalar \
        (*(PyArray_Descr * (*)(PyObject *)) \
         PyArray_API[57])
#define PyArray_DescrFromTypeObject \
        (*(PyArray_Descr * (*)(PyObject *)) \
         PyArray_API[58])
#define PyArray_Size \
        (*(npy_intp (*)(PyObject *)) \
         PyArray_API[59])
#define PyArray_Scalar \
        (*(PyObject * (*)(void *, PyArray_Descr *, PyObject *)) \
         PyArray_API[60])
#define PyArray_FromScalar \
        (*(PyObject * (*)(PyObject *, PyArray_Descr *)) \
         PyArray_API[61])
#define PyArray_ScalarAsCtype \
        (*(void (*)(PyObject *, void *)) \
         PyArray_API[62])
#define PyArray_CastScalarToCtype \
        (*(int (*)(PyObject *, void *, PyArray_Descr *)) \
         PyArray_API[63])
#define PyArray_CastScalarDirect \
        (*(int (*)(PyObject *, PyArray_Descr *, void *, int)) \
         PyArray_API[64])
#define PyArray_ScalarFromObject \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[65])
#define PyArray_GetCastFunc \
        (*(PyArray_VectorUnaryFunc * (*)(PyArray_Descr *, int)) \
         PyArray_API[66])
#define PyArray_FromDims \
        (*(PyObject * (*)(int, int *, int)) \
         PyArray_API[67])
#define PyArray_FromDimsAndDataAndDescr \
        (*(PyObject * (*)(int, int *, PyArray_Descr *, char *)) \
         PyArray_API[68])
#define PyArray_FromAny \
        (*(PyObject * (*)(PyObject *, PyArray_Descr *, int, int, int, PyObject *)) \
         PyArray_API[69])
#define PyArray_EnsureArray \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[70])
#define PyArray_EnsureAnyArray \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[71])
#define PyArray_FromFile \
        (*(PyObject * (*)(FILE *, PyArray_Descr *, npy_intp, char *)) \
         PyArray_API[72])
#define PyArray_FromString \
        (*(PyObject * (*)(char *, npy_intp, PyArray_Descr *, npy_intp, char *)) \
         PyArray_API[73])
#define PyArray_FromBuffer \
        (*(PyObject * (*)(PyObject *, PyArray_Descr *, npy_intp, npy_intp)) \
         PyArray_API[74])
#define PyArray_FromIter \
        (*(PyObject * (*)(PyObject *, PyArray_Descr *, npy_intp)) \
         PyArray_API[75])
#define PyArray_Return \
        (*(PyObject * (*)(PyArrayObject *)) \
         PyArray_API[76])
#define PyArray_GetField \
        (*(PyObject * (*)(PyArrayObject *, PyArray_Descr *, int)) \
         PyArray_API[77])
#define PyArray_SetField \
        (*(int (*)(PyArrayObject *, PyArray_Descr *, int, PyObject *)) \
         PyArray_API[78])
#define PyArray_Byteswap \
        (*(PyObject * (*)(PyArrayObject *, npy_bool)) \
         PyArray_API[79])
#define PyArray_Resize \
        (*(PyObject * (*)(PyArrayObject *, PyArray_Dims *, int, NPY_ORDER)) \
         PyArray_API[80])
#define PyArray_MoveInto \
        (*(int (*)(PyArrayObject *, PyArrayObject *)) \
         PyArray_API[81])
#define PyArray_CopyInto \
        (*(int (*)(PyArrayObject *, PyArrayObject *)) \
         PyArray_API[82])
#define PyArray_CopyAnyInto \
        (*(int (*)(PyArrayObject *, PyArrayObject *)) \
         PyArray_API[83])
#define PyArray_CopyObject \
        (*(int (*)(PyArrayObject *, PyObject *)) \
         PyArray_API[84])
#define PyArray_NewCopy \
        (*(PyObject * (*)(PyArrayObject *, NPY_ORDER)) \
         PyArray_API[85])
#define PyArray_ToList \
        (*(PyObject * (*)(PyArrayObject *)) \
         PyArray_API[86])
#define PyArray_ToString \
        (*(PyObject * (*)(PyArrayObject *, NPY_ORDER)) \
         PyArray_API[87])
#define PyArray_ToFile \
        (*(int (*)(PyArrayObject *, FILE *, char *, char *)) \
         PyArray_API[88])
#define PyArray_Dump \
        (*(int (*)(PyObject *, PyObject *, int)) \
         PyArray_API[89])
#define PyArray_Dumps \
        (*(PyObject * (*)(PyObject *, int)) \
         PyArray_API[90])
#define PyArray_ValidType \
        (*(int (*)(int)) \
         PyArray_API[91])
#define PyArray_UpdateFlags \
        (*(void (*)(PyArrayObject *, int)) \
         PyArray_API[92])
#define PyArray_New \
        (*(PyObject * (*)(PyTypeObject *, int, npy_intp *, int, npy_intp *, void *, int, int, PyObject *)) \
         PyArray_API[93])
#define PyArray_NewFromDescr \
        (*(PyObject * (*)(PyTypeObject *, PyArray_Descr *, int, npy_intp *, npy_intp *, void *, int, PyObject *)) \
         PyArray_API[94])
#define PyArray_DescrNew \
        (*(PyArray_Descr * (*)(PyArray_Descr *)) \
         PyArray_API[95])
#define PyArray_DescrNewFromType \
        (*(PyArray_Descr * (*)(int)) \
         PyArray_API[96])
#define PyArray_GetPriority \
        (*(double (*)(PyObject *, double)) \
         PyArray_API[97])
#define PyArray_IterNew \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[98])
#define PyArray_MultiIterNew \
        (*(PyObject * (*)(int, ...)) \
         PyArray_API[99])
#define PyArray_PyIntAsInt \
        (*(int (*)(PyObject *)) \
         PyArray_API[100])
#define PyArray_PyIntAsIntp \
        (*(npy_intp (*)(PyObject *)) \
         PyArray_API[101])
#define PyArray_Broadcast \
        (*(int (*)(PyArrayMultiIterObject *)) \
         PyArray_API[102])
#define PyArray_FillObjectArray \
        (*(void (*)(PyArrayObject *, PyObject *)) \
         PyArray_API[103])
#define PyArray_FillWithScalar \
        (*(int (*)(PyArrayObject *, PyObject *)) \
         PyArray_API[104])
#define PyArray_CheckStrides \
        (*(npy_bool (*)(int, int, npy_intp, npy_intp, npy_intp *, npy_intp *)) \
         PyArray_API[105])
#define PyArray_DescrNewByteorder \
        (*(PyArray_Descr * (*)(PyArray_Descr *, char)) \
         PyArray_API[106])
#define PyArray_IterAllButAxis \
        (*(PyObject * (*)(PyObject *, int *)) \
         PyArray_API[107])
#define PyArray_CheckFromAny \
        (*(PyObject * (*)(PyObject *, PyArray_Descr *, int, int, int, PyObject *)) \
         PyArray_API[108])
#define PyArray_FromArray \
        (*(PyObject * (*)(PyArrayObject *, PyArray_Descr *, int)) \
         PyArray_API[109])
#define PyArray_FromInterface \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[110])
#define PyArray_FromStructInterface \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[111])
#define PyArray_FromArrayAttr \
        (*(PyObject * (*)(PyObject *, PyArray_Descr *, PyObject *)) \
         PyArray_API[112])
#define PyArray_ScalarKind \
        (*(NPY_SCALARKIND (*)(int, PyArrayObject **)) \
         PyArray_API[113])
#define PyArray_CanCoerceScalar \
        (*(int (*)(int, int, NPY_SCALARKIND)) \
         PyArray_API[114])
#define PyArray_NewFlagsObject \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[115])
#define PyArray_CanCastScalar \
        (*(npy_bool (*)(PyTypeObject *, PyTypeObject *)) \
         PyArray_API[116])
#define PyArray_CompareUCS4 \
        (*(int (*)(npy_ucs4 *, npy_ucs4 *, register size_t)) \
         PyArray_API[117])
#define PyArray_RemoveSmallest \
        (*(int (*)(PyArrayMultiIterObject *)) \
         PyArray_API[118])
#define PyArray_ElementStrides \
        (*(int (*)(PyObject *)) \
         PyArray_API[119])
#define PyArray_Item_INCREF \
        (*(void (*)(char *, PyArray_Descr *)) \
         PyArray_API[120])
#define PyArray_Item_XDECREF \
        (*(void (*)(char *, PyArray_Descr *)) \
         PyArray_API[121])
#define PyArray_FieldNames \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[122])
#define PyArray_Transpose \
        (*(PyObject * (*)(PyArrayObject *, PyArray_Dims *)) \
         PyArray_API[123])
#define PyArray_TakeFrom \
        (*(PyObject * (*)(PyArrayObject *, PyObject *, int, PyArrayObject *, NPY_CLIPMODE)) \
         PyArray_API[124])
#define PyArray_PutTo \
        (*(PyObject * (*)(PyArrayObject *, PyObject*, PyObject *, NPY_CLIPMODE)) \
         PyArray_API[125])
#define PyArray_PutMask \
        (*(PyObject * (*)(PyArrayObject *, PyObject*, PyObject*)) \
         PyArray_API[126])
#define PyArray_Repeat \
        (*(PyObject * (*)(PyArrayObject *, PyObject *, int)) \
         PyArray_API[127])
#define PyArray_Choose \
        (*(PyObject * (*)(PyArrayObject *, PyObject *, PyArrayObject *, NPY_CLIPMODE)) \
         PyArray_API[128])
#define PyArray_Sort \
        (*(int (*)(PyArrayObject *, int, NPY_SORTKIND)) \
         PyArray_API[129])
#define PyArray_ArgSort \
        (*(PyObject * (*)(PyArrayObject *, int, NPY_SORTKIND)) \
         PyArray_API[130])
#define PyArray_SearchSorted \
        (*(PyObject * (*)(PyArrayObject *, PyObject *, NPY_SEARCHSIDE)) \
         PyArray_API[131])
#define PyArray_ArgMax \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[132])
#define PyArray_ArgMin \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[133])
#define PyArray_Reshape \
        (*(PyObject * (*)(PyArrayObject *, PyObject *)) \
         PyArray_API[134])
#define PyArray_Newshape \
        (*(PyObject * (*)(PyArrayObject *, PyArray_Dims *, NPY_ORDER)) \
         PyArray_API[135])
#define PyArray_Squeeze \
        (*(PyObject * (*)(PyArrayObject *)) \
         PyArray_API[136])
#define PyArray_View \
        (*(PyObject * (*)(PyArrayObject *, PyArray_Descr *, PyTypeObject *)) \
         PyArray_API[137])
#define PyArray_SwapAxes \
        (*(PyObject * (*)(PyArrayObject *, int, int)) \
         PyArray_API[138])
#define PyArray_Max \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[139])
#define PyArray_Min \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[140])
#define PyArray_Ptp \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[141])
#define PyArray_Mean \
        (*(PyObject * (*)(PyArrayObject *, int, int, PyArrayObject *)) \
         PyArray_API[142])
#define PyArray_Trace \
        (*(PyObject * (*)(PyArrayObject *, int, int, int, int, PyArrayObject *)) \
         PyArray_API[143])
#define PyArray_Diagonal \
        (*(PyObject * (*)(PyArrayObject *, int, int, int)) \
         PyArray_API[144])
#define PyArray_Clip \
        (*(PyObject * (*)(PyArrayObject *, PyObject *, PyObject *, PyArrayObject *)) \
         PyArray_API[145])
#define PyArray_Conjugate \
        (*(PyObject * (*)(PyArrayObject *, PyArrayObject *)) \
         PyArray_API[146])
#define PyArray_Nonzero \
        (*(PyObject * (*)(PyArrayObject *)) \
         PyArray_API[147])
#define PyArray_Std \
        (*(PyObject * (*)(PyArrayObject *, int, int, PyArrayObject *, int)) \
         PyArray_API[148])
#define PyArray_Sum \
        (*(PyObject * (*)(PyArrayObject *, int, int, PyArrayObject *)) \
         PyArray_API[149])
#define PyArray_CumSum \
        (*(PyObject * (*)(PyArrayObject *, int, int, PyArrayObject *)) \
         PyArray_API[150])
#define PyArray_Prod \
        (*(PyObject * (*)(PyArrayObject *, int, int, PyArrayObject *)) \
         PyArray_API[151])
#define PyArray_CumProd \
        (*(PyObject * (*)(PyArrayObject *, int, int, PyArrayObject *)) \
         PyArray_API[152])
#define PyArray_All \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[153])
#define PyArray_Any \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[154])
#define PyArray_Compress \
        (*(PyObject * (*)(PyArrayObject *, PyObject *, int, PyArrayObject *)) \
         PyArray_API[155])
#define PyArray_Flatten \
        (*(PyObject * (*)(PyArrayObject *, NPY_ORDER)) \
         PyArray_API[156])
#define PyArray_Ravel \
        (*(PyObject * (*)(PyArrayObject *, NPY_ORDER)) \
         PyArray_API[157])
#define PyArray_MultiplyList \
        (*(npy_intp (*)(register npy_intp *, register int)) \
         PyArray_API[158])
#define PyArray_MultiplyIntList \
        (*(int (*)(register int *, register int)) \
         PyArray_API[159])
#define PyArray_GetPtr \
        (*(void * (*)(PyArrayObject *, register npy_intp*)) \
         PyArray_API[160])
#define PyArray_CompareLists \
        (*(int (*)(npy_intp *, npy_intp *, int)) \
         PyArray_API[161])
#define PyArray_AsCArray \
        (*(int (*)(PyObject **, void *, npy_intp *, int, PyArray_Descr*)) \
         PyArray_API[162])
#define PyArray_As1D \
        (*(int (*)(PyObject **, char **, int *, int)) \
         PyArray_API[163])
#define PyArray_As2D \
        (*(int (*)(PyObject **, char ***, int *, int *, int)) \
         PyArray_API[164])
#define PyArray_Free \
        (*(int (*)(PyObject *, void *)) \
         PyArray_API[165])
#define PyArray_Converter \
        (*(int (*)(PyObject *, PyObject **)) \
         PyArray_API[166])
#define PyArray_IntpFromSequence \
        (*(int (*)(PyObject *, npy_intp *, int)) \
         PyArray_API[167])
#define PyArray_Concatenate \
        (*(PyObject * (*)(PyObject *, int)) \
         PyArray_API[168])
#define PyArray_InnerProduct \
        (*(PyObject * (*)(PyObject *, PyObject *)) \
         PyArray_API[169])
#define PyArray_MatrixProduct \
        (*(PyObject * (*)(PyObject *, PyObject *)) \
         PyArray_API[170])
#define PyArray_CopyAndTranspose \
        (*(PyObject * (*)(PyObject *)) \
         PyArray_API[171])
#define PyArray_Correlate \
        (*(PyObject * (*)(PyObject *, PyObject *, int)) \
         PyArray_API[172])
#define PyArray_TypestrConvert \
        (*(int (*)(int, int)) \
         PyArray_API[173])
#define PyArray_DescrConverter \
        (*(int (*)(PyObject *, PyArray_Descr **)) \
         PyArray_API[174])
#define PyArray_DescrConverter2 \
        (*(int (*)(PyObject *, PyArray_Descr **)) \
         PyArray_API[175])
#define PyArray_IntpConverter \
        (*(int (*)(PyObject *, PyArray_Dims *)) \
         PyArray_API[176])
#define PyArray_BufferConverter \
        (*(int (*)(PyObject *, PyArray_Chunk *)) \
         PyArray_API[177])
#define PyArray_AxisConverter \
        (*(int (*)(PyObject *, int *)) \
         PyArray_API[178])
#define PyArray_BoolConverter \
        (*(int (*)(PyObject *, npy_bool *)) \
         PyArray_API[179])
#define PyArray_ByteorderConverter \
        (*(int (*)(PyObject *, char *)) \
         PyArray_API[180])
#define PyArray_OrderConverter \
        (*(int (*)(PyObject *, NPY_ORDER *)) \
         PyArray_API[181])
#define PyArray_EquivTypes \
        (*(unsigned char (*)(PyArray_Descr *, PyArray_Descr *)) \
         PyArray_API[182])
#define PyArray_Zeros \
        (*(PyObject * (*)(int, npy_intp *, PyArray_Descr *, int)) \
         PyArray_API[183])
#define PyArray_Empty \
        (*(PyObject * (*)(int, npy_intp *, PyArray_Descr *, int)) \
         PyArray_API[184])
#define PyArray_Where \
        (*(PyObject * (*)(PyObject *, PyObject *, PyObject *)) \
         PyArray_API[185])
#define PyArray_Arange \
        (*(PyObject * (*)(double, double, double, int)) \
         PyArray_API[186])
#define PyArray_ArangeObj \
        (*(PyObject * (*)(PyObject *, PyObject *, PyObject *, PyArray_Descr *)) \
         PyArray_API[187])
#define PyArray_SortkindConverter \
        (*(int (*)(PyObject *, NPY_SORTKIND *)) \
         PyArray_API[188])
#define PyArray_LexSort \
        (*(PyObject * (*)(PyObject *, int)) \
         PyArray_API[189])
#define PyArray_Round \
        (*(PyObject * (*)(PyArrayObject *, int, PyArrayObject *)) \
         PyArray_API[190])
#define PyArray_EquivTypenums \
        (*(unsigned char (*)(int, int)) \
         PyArray_API[191])
#define PyArray_RegisterDataType \
        (*(int (*)(PyArray_Descr *)) \
         PyArray_API[192])
#define PyArray_RegisterCastFunc \
        (*(int (*)(PyArray_Descr *, int, PyArray_VectorUnaryFunc *)) \
         PyArray_API[193])
#define PyArray_RegisterCanCast \
        (*(int (*)(PyArray_Descr *, int, NPY_SCALARKIND)) \
         PyArray_API[194])
#define PyArray_InitArrFuncs \
        (*(void (*)(PyArray_ArrFuncs *)) \
         PyArray_API[195])
#define PyArray_IntTupleFromIntp \
        (*(PyObject * (*)(int, npy_intp *)) \
         PyArray_API[196])
#define PyArray_TypeNumFromName \
        (*(int (*)(char *)) \
         PyArray_API[197])
#define PyArray_ClipmodeConverter \
        (*(int (*)(PyObject *, NPY_CLIPMODE *)) \
         PyArray_API[198])
#define PyArray_OutputConverter \
        (*(int (*)(PyObject *, PyArrayObject **)) \
         PyArray_API[199])
#define PyArray_BroadcastToShape \
        (*(PyObject * (*)(PyObject *, npy_intp *, int)) \
         PyArray_API[200])
#define _PyArray_SigintHandler \
        (*(void (*)(int)) \
         PyArray_API[201])
#define _PyArray_GetSigintBuf \
        (*(void* (*)(void)) \
         PyArray_API[202])
#define PyArray_DescrAlignConverter \
        (*(int (*)(PyObject *, PyArray_Descr **)) \
         PyArray_API[203])
#define PyArray_DescrAlignConverter2 \
        (*(int (*)(PyObject *, PyArray_Descr **)) \
         PyArray_API[204])
#define PyArray_SearchsideConverter \
        (*(int (*)(PyObject *, void *)) \
         PyArray_API[205])

#if !defined(NO_IMPORT_ARRAY) && !defined(NO_IMPORT)
static int
_import_array(void)
{
  PyObject *numpy = PyImport_ImportModule("numpy.core.multiarray");
  PyObject *c_api = NULL;
  if (numpy == NULL) return -1;
  c_api = PyObject_GetAttrString(numpy, "_ARRAY_API");
  if (c_api == NULL) {Py_DECREF(numpy); return -1;}
  if (PyCObject_Check(c_api)) {
      PyArray_API = (void **)PyCObject_AsVoidPtr(c_api);
  }
  Py_DECREF(c_api);
  Py_DECREF(numpy);
  if (PyArray_API == NULL) return -1;
  /* Perform runtime check of C API version */
  if (NPY_VERSION != PyArray_GetNDArrayCVersion()) {
    PyErr_Format(PyExc_RuntimeError, "module compiled against "\
        "version %x of C-API but this version of numpy is %x", \
        (int) NPY_VERSION, (int) PyArray_GetNDArrayCVersion());
    return -1;
  }
  return 0;
}

#define import_array() {if (_import_array() < 0) {PyErr_Print(); PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import"); return; } }

#define import_array1(ret) {if (_import_array() < 0) {PyErr_Print(); PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import"); return ret; } }

#define import_array2(msg, ret) {if (_import_array() < 0) {PyErr_Print(); PyErr_SetString(PyExc_ImportError, msg); return ret; } }

#endif

#endif
