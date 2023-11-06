
#include "hhgrid.h"

#include <Python.h>
#include <stdlib.h> // setenv, strdup, strtok, NULL
#include <assert.h> // assert
#include <unistd.h> // access

// BOF - Helper Functions (mostly for Fortran)
void python_initialize()
{
    Py_Initialize();
    if (!Py_IsInitialized())
    {
        printf("ERROR: Failed to initialize Python interpreter\n");
        exit(1);
    }
};

void python_decref(PyObject* grid)
{
    Py_XDECREF(grid);
};

void python_finalize()
{
    Py_Finalize();
};

void python_printinfo()
{
    const wchar_t* programFullPath = Py_GetProgramFullPath();
    const char* getVersion = Py_GetVersion();
    const wchar_t* getPythonHome = Py_GetPythonHome();
    const wchar_t* getPath = Py_GetPath();

    printf("== Python Parameters ==\n");
    printf("Py_GetProgramFullPath: %ls\n", programFullPath);
    printf("Py_GetVersion: %s\n", getVersion);
    printf("Py_GetPythonHome: %ls\n", getPythonHome);
    printf("Py_GetPath: %ls\n", getPath);
    printf("\n");
};
// EOF - Helper Functions (mostly for Fortran)

void combine_grids(const char* grid_temp, double cHHH, double ct, double ctt, double cg, double cgg)
{
    char* values[23] = {"+5.000000E-01_+4.782609E-01_+1.000000E+00_+6.875000E-01_+8.888889E-01",
  "-1.000000E+00_-2.500000E+00_+2.857143E-01_+4.666667E-01_+1.818182E-01",
  "+4.545455E-02_-5.263158E-02_-6.875000E-01_+2.941176E-01_+6.923077E-01",
  "+1.250000E-01_+1.111111E-01_-3.333333E-01_+2.173913E-01_-8.000000E-01",
  "-1.142857E+00_+4.444444E-01_+9.090909E-02_+7.272727E-01_+6.428571E-01",
  "+8.333333E-01_+4.545455E-01_-4.285714E-01_+9.523810E-02_+5.263158E-01",
  "+1.000000E+00_+1.428571E-01_-3.333333E-01_+5.217391E-01_-1.500000E+00",
  "-5.263158E-02_-7.142857E-02_+2.083333E-01_+6.923077E-01_+9.000000E-01",
  "-3.750000E-01_+4.705882E-01_-3.750000E-01_-6.666667E-01_+4.166667E-01",
  "-6.470588E-01_-8.333333E-01_+5.000000E-01_-1.333333E+00_+5.000000E-01",
  "-6.666667E-01_+1.666667E-01_+1.100000E+00_+5.000000E-01_+1.000000E+00",
  "+2.000000E-01_+6.000000E-01_+8.695652E-02_-3.478261E-01_-1.333333E+00",
  "-4.210526E-01_+4.444444E-01_-2.222222E-01_+3.809524E-01_-5.714286E-01",
  "-1.176471E-01_+2.200000E+00_-7.692308E-02_-1.875000E-01_+5.555556E-01",
  "+2.000000E-01_-8.571429E-01_-1.000000E+00_+3.125000E-01_-1.166667E+00",
  "+3.000000E-01_+1.111111E-01_+1.285714E+00_+1.285714E+00_-4.615385E-01",
  "-4.347826E-01_-8.000000E-01_+1.111111E-01_-6.315789E-01_+4.347826E-02",
  "-1.142857E+00_-3.333333E-01_-5.000000E-01_-5.000000E-01_+4.117647E-01",
  "+2.250000E+00_-6.666667E-01_+2.727273E-01_+3.571429E-01_-1.000000E+00",
  "+6.111111E-01_+2.777778E-01_+1.111111E-01_-8.000000E-01_+2.272727E-01",
  "+2.173913E-01_+3.000000E+00_-5.263158E-01_+4.761905E-02_-3.809524E-01",
  "+4.545455E-01_+4.000000E-01_-1.500000E+00_+5.454545E-01_+6.428571E-01",
  "+2.500000E-01_+1.111111E+00_-4.166667E-01_-4.444444E-01_+5.000000E-02"};

    int search_paths = 1;
    char* delims = ":";
    char* path_sep = "/";
    char* grid_file_path;
    size_t len_path_sep = strlen(path_sep);
    char* pythonpath = strdup(getenv("PYTHONPATH"));
    char* result = strtok( pythonpath, delims );

    char grid_tmp[10];
    strncpy(grid_tmp, grid_temp, 9); // Only take the first characters to look for the three basic cHHH grids
    grid_tmp[9]='\0';

    for (int i=0; i<(sizeof(values) / sizeof(values[0])); ++i) {
      size_t len_grid_name = strlen(grid_tmp) + strlen("_") + strlen(values[i]) + strlen(".grid") + 1;  // +14 for _***.grid
      char* grid_name = (char*) malloc(len_grid_name);
      memcpy(grid_name, grid_tmp, strlen(grid_tmp));
      memcpy(grid_name + strlen(grid_tmp), "_", strlen("_"));
      memcpy(grid_name + strlen(grid_tmp) + strlen("_"), values[i], strlen(values[i]));
      memcpy(grid_name + strlen(grid_tmp) + strlen("_") + strlen(values[i]), ".grid", strlen(".grid"));
      grid_name[len_grid_name-1] = '\0';

    // Check if grid_name is accessible as-is    
    if( access(grid_name, F_OK) != -1 && access(grid_name, R_OK) != -1 )
    {
        printf("Looking for %s in current directory. Found\n", grid_name);
        grid_file_path = (char*) malloc(len_grid_name + 1); // +1 for null terminator
        memcpy(grid_file_path, grid_name, len_grid_name + 1);
        search_paths = 0;
        setenv("PYTHONPATH", ".", 1); // Set PYTHONPATH to look here
    }

    // Else search PYTHONPATH for grid_name
    while( search_paths == 1 && result != NULL )
    {
        size_t len_result = strlen(result);
        grid_file_path = (char*) malloc (len_result + len_path_sep + len_grid_name + 1); // +1 for null terminator
        memcpy(grid_file_path, result, len_result);
        memcpy(grid_file_path + len_result, path_sep, len_path_sep);
        memcpy(grid_file_path + len_result + len_path_sep, grid_name, len_grid_name + 1); // +1 for null terminator
        printf("Searching for %s in: %s ", grid_name, grid_file_path);
        if( access(grid_file_path, F_OK) != -1 && access(grid_file_path, R_OK) != -1 )
        {
            printf("found\n");
            search_paths = 0;
        } else
        {
            printf("not found\n");
            result = strtok( NULL, delims );
            if(result == NULL)
            {
                printf("ERROR: Failed to find grid\n");
                exit(1);
            }
        }
    }
    }

    PyObject* pModule = PyImport_ImportModule("creategrid");
    if(pModule == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to load creategrid.py: please check that you have numpy and scipy installed\n");
    }
    assert(pModule != NULL);

    PyObject* pFct = PyObject_GetAttrString(pModule, "combinegrids");
    if(pFct == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to locate CreateGrid class: please check that you have the latest version of creategrid.py and Python 3.x\n");
    }
    assert(pFct != NULL);

    PyObject* pGridName = PyUnicode_FromString(grid_temp);
    if(pGridName == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python string from grid_temp: please check that grid_temp is a valid string\n");
    }
    assert(pGridName != NULL);
    
    PyObject* pcHHHValue = PyFloat_FromDouble(cHHH);
    if(pcHHHValue == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python double from cHHH: please check that cHHH is a valid double\n");
    }
    assert(pcHHHValue != NULL);

    PyObject* pctValue = PyFloat_FromDouble(ct);
    if(pctValue == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python double from ct: please check that ct is a valid double\n");
    }
    assert(pctValue != NULL);

    PyObject* pcttValue = PyFloat_FromDouble(ctt);
    if(pcttValue == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python double from ctt: please check that ctt is a valid double\n");
    }
    assert(pcttValue != NULL);

    PyObject* pcgValue = PyFloat_FromDouble(cg);
    if(pcgValue == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python double from cg: please check that cg is a valid double\n");
    }
    assert(pcgValue != NULL);

    PyObject* pcggValue = PyFloat_FromDouble(cgg);
    if(pcggValue == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python double from cgg: please check that cgg is a valid double\n");
    }
    assert(pcggValue != NULL);

    PyObject* pGridNameTuple = PyTuple_Pack(6,pGridName,pcHHHValue,pctValue,pcttValue,pcgValue,pcggValue);
    if(pGridNameTuple == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python tuple: please check that your Python version is 3.x\n");
    }
    assert(pGridNameTuple != NULL);

    PyObject* pFunct = PyObject_CallObject(pFct, pGridNameTuple);
    if(pFunct == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to call combinegrids: please check that you have the latest version of creategrid.py and Python 3.x\n");
    }
    assert(pFunct != NULL);

    // Cleanup
    free(pythonpath);
    free(grid_file_path);
    Py_XDECREF(pModule);
    Py_XDECREF(pFct);
    Py_XDECREF(pGridName);
    Py_XDECREF(pcHHHValue);
    Py_XDECREF(pGridNameTuple);
    Py_XDECREF(pFunct);
};

PyObject* grid_initialize(const char* grid_name)
{
    int search_paths = 1;
    char* delims = ":";
    char* path_sep = "/";
    char* grid_file_path;
    size_t len_path_sep = strlen(path_sep);
    size_t len_grid_name = strlen(grid_name);
    char* pythonpath = strdup(getenv("PYTHONPATH"));
    char* result = strtok( pythonpath, delims );

    // Check if grid_name is accessible as-is
    if( access(grid_name, F_OK) != -1 && access(grid_name, R_OK) != -1 )
    {
	printf("Looking for %s in current directory. Found\n", grid_name);
        grid_file_path = (char*) malloc(len_grid_name + 1); // +1 for null terminator
        memcpy(grid_file_path, grid_name, len_grid_name + 1);
        search_paths = 0;
        setenv("PYTHONPATH", ".", 1); // Set PYTHONPATH to look here
    }

    // Else search PYTHONPATH for grid_name
    while( search_paths == 1 && result != NULL )
    {
        size_t len_result = strlen(result);
        grid_file_path = (char*) malloc (len_result + len_path_sep + len_grid_name + 1); // +1 for null terminator
        memcpy(grid_file_path, result, len_result);
        memcpy(grid_file_path + len_result, path_sep, len_path_sep);
        memcpy(grid_file_path + len_result + len_path_sep, grid_name, len_grid_name + 1); // +1 for null terminator
        printf("Searching for %s in: %s ", grid_name, grid_file_path);
        if( access(grid_file_path, F_OK) != -1 && access(grid_file_path, R_OK) != -1 )
        {
            printf("found\n");
            search_paths = 0;
        } else
        {
            printf("not found\n");
            result = strtok( NULL, delims );
            if(result == NULL)
            {
                printf("ERROR: Failed to find grid\n");
                exit(1);
            }
        }
    }

    printf("PYTHONPATH: %s\n", getenv("PYTHONPATH"));
    printf("Grid Path: %s\n", grid_file_path);

    PyObject* pModule = PyImport_ImportModule("creategrid");
    if(pModule == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to load creategrid.py: please check that you have numpy and scipy installed\n");
    }
    assert(pModule != NULL);

    PyObject* pClass = PyObject_GetAttrString(pModule, "CreateGrid");
    if(pClass == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to locate CreateGrid class: please check that you have the latest version of creategrid.py and Python 3.x\n");
    }
    assert(pClass != NULL);

    if(!PyCallable_Check(pClass))
    {
        PyErr_Print();
        printf("ERROR: CreateGrid is not callable: please check that you have the latest version of creategrid.py and Python 3.x\n");
    }

    PyObject* pGridName = PyUnicode_FromString(grid_file_path);
    if(pGridName == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python string from grid_file_path: please check that grid_file_path is a valid string\n");
    }
    assert(pGridName != NULL);

    PyObject* pGridNameTuple = PyTuple_Pack(1,pGridName);
    if(pGridNameTuple == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create Python tuple: please check that your Python version is 3.x\n");
    }
    assert(pGridNameTuple != NULL);

    PyObject* pInstance = PyObject_CallObject(pClass, pGridNameTuple);
    if(pInstance == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to create instance of CreateGrid: please check that you have the latest version of creategrid.py and Python 3.x\n");
    }
    assert(pInstance != NULL);

    // Cleanup
    free(pythonpath);
    free(grid_file_path);
    Py_XDECREF(pModule);
    Py_XDECREF(pClass);
    Py_XDECREF(pGridName);
    Py_XDECREF(pGridNameTuple);

    return pInstance;
};

double grid_virt(PyObject* grid, double s, double t)
{
    PyObject* pResult = PyObject_CallMethod(grid, "GetAmplitude", "(ff)", s, t);
    if(pResult == NULL)
    {
        PyErr_Print();
        printf("ERROR: Failed to compute result\n");
    }
    assert(pResult != NULL);

    double result = PyFloat_AsDouble(pResult);

    // Cleanup
    Py_XDECREF(pResult);

    return result;
};
