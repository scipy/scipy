/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HToyIO.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifdef __cplusplus
// Includes for C++ compiler
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
using std::string;
#endif

// Includes for C compiler [C++ compiler ignores them]
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include "io/HToyIO.h"
#include "lp_data/HConst.h"

int readToy_LP_cpp(const char* filename, int* m_p, int* n_p, int* maxmin,
                   double* offset, double** A, double** b, double** c,
                   double** lb, double** ub) {
  int* integerColumn;
  int RtCd = readToy_MIP_cpp(filename, m_p, n_p, maxmin, offset, A, b, c, lb,
                             ub, &integerColumn);
  return RtCd;
}
extern "C" int readToy_LP_c(const char* filename, int* m_p, int* n_p,
                            int* maxmin, double* offset, double** A, double** b,
                            double** c, double** lb, double** ub) {
  int RtCd =
      readToy_LP_cpp(filename, m_p, n_p, maxmin, offset, A, b, c, lb, ub);
  return RtCd;
}

int readToy_MIP_cpp(const char* filename, int* m_p, int* n_p, int* maxmin,
                    double* offset, double** A, double** b, double** c,
                    double** lb, double** ub, int** integerColumn) {
  const int buff_ln = 255;
  char buff[buff_ln];
  char* fgets_rt;
  int n, m;
  int start, end;
  int rp_rd = 0;
  double rl_v;
  int i_v;
  char* term;
  int i, j, k;

#ifdef HiGHSDEV
  printf("readMPS: Trying to open file %s\n", filename);
#endif
  FILE* file_p = fopen(filename, "r");
  if (file_p == 0) {
#ifdef HiGHSDEV
    printf("readInput: Not opened file OK\n");
#endif
    return 1;
  }
#ifdef HiGHSDEV
  printf("readInput: Opened file  OK\n");
#endif
  fgets_rt = fgets(buff, buff_ln, file_p);

  sscanf(buff, "%d %d %d", m_p, n_p, maxmin);
  if (rp_rd > 0)
    printf("Read %d %d %d for m_p, n_p, maxmin\n", *(m_p), *(n_p), *(maxmin));
  n = *(n_p);
  m = *(m_p);

  // Allocate arrays
  *c = (double*)malloc(sizeof(double) * n);
  *A = (double*)malloc(sizeof(double) * (m * n));
  *b = (double*)malloc(sizeof(double) * m);
  *lb = (double*)malloc(sizeof(double) * n);
  *ub = (double*)malloc(sizeof(double) * n);
  *integerColumn = (int*)malloc(sizeof(int) * n);

  fgets_rt = fgets(buff, buff_ln, file_p);
  start = 0;
  while (buff[start] == ' ' && buff[start] != '\0') {
    //      printf("buff[%3d] = %s\n", start, buff[start]);
    start++;
  }

  // Read cost and offset
  j = 0;
  start = 0;
  while (j <= n) {
    if (rp_rd > 0) printf("\n Col %d\n", j);
    while (buff[start] == ' ') start++;

    end = start;
    while (buff[end] != ' ' && buff[end] != '\0') end++;

    int len = end - start;
    if (rp_rd > 0)
      printf("buff: start =%2d, len =%2d, end =%2d\n", start, len, end);
    if (len > 0) {
      term = (char*)malloc(sizeof(char) * len);
      for (k = 0; k < len; k++) term[k] = buff[start + k];
    } else {
      // JAJH: QY - LEN<=0 means that the end of the line has been reached
      printf("Error reading file: not enough numbers on row of objective.\n");
      return 1;
    }

    rl_v = atof(term);
    if (j < n) {
      (*c)[j] = rl_v;
      if (rp_rd > 0) printf("Read %s giving %g for (*c)[%d]\n", term, rl_v, j);
      j++;
      if (buff[end] == '\0') {
        printf("Error reading file: not enough numbers on row of objective.\n");
        return 1;
      }
      start = end + 1;
    } else {
      *offset = rl_v;
      if (rp_rd > 0) printf("Read %s giving %g for offset\n", term, rl_v);
      j++;
    }
  }

  // Read matrix and RHS
  i = 0;
  while (i < m) {
    fgets_rt = fgets(buff, buff_ln, file_p);

    start = 0;
    while (buff[start] == ' ' && buff[start] != '\0') start++;

    // JAJH: QY - This doesn't seem to do what was intended
    if (buff[start] == '\0')
      printf("Error reading file: not enough numbers on row %d of A.\n", i);

    j = 0;
    while (j < n) {
      while (buff[start] == ' ') start++;

      end = start;
      while (buff[end] != ' ' && buff[end] != '\0') end++;

      int len = end - start;

      if (len > 0) {
        term = (char*)malloc(sizeof(char) * len);
        for (k = 0; k < len; k++) term[k] = buff[start + k];
      } else {
        // JAJH: QY - LEN<=0 means that the end of the line has been reached
        printf("Error reading file: not enough numbers on row for RHS.\n");
        return 1;
      }

      rl_v = atof(term);

      // Row-wise
      //	int ind = (i * n) + j;
      // Column-wise
      int ind = (j * m) + i;
      (*A)[ind] = rl_v;
      if (rp_rd > 0)
        printf("Read %s giving %g for (*A)[%d]\n", term, rl_v, ind);
      j++;
      if (buff[end] == '\0') {
        printf("Error reading file: not enough numbers on constraint row.\n");
        return 1;
      }
      start = end + 1;
    }

    // RHS
    while (buff[start] == ' ') start++;

    end = start;
    while (buff[end] != ' ' && buff[end] != '\0') end++;

    int len = end - start;
    if (len > 0) {
      term = (char*)malloc(sizeof(char) * len);
      for (k = 0; k < len; k++) term[k] = buff[start + k];
    } else {
      // JAJH: QY - LEN<=0 means that the end of the line has been reached
      printf("Error reading file: not enough numbers on constraint row.\n");
      return 1;
    }

    rl_v = atof(term);
    (*b)[i] = rl_v;
    if (rp_rd > 0) printf("Read %s giving %g for (*b)[%d]\n", term, rl_v, i);
    i++;
  }
  // Set up default values of lower bounds, upper bounds and integer information
  for (int c_n = 0; c_n < n; c_n++) {
    (*lb)[c_n] = 0;
    (*ub)[c_n] = HIGHS_CONST_INF;
    (*integerColumn)[c_n] = 0;
  }

  // Read lower bounds on columns
  fgets_rt = fgets(buff, buff_ln, file_p);
  if (rp_rd > 0) printf("LB: fgets return value = %s\n", fgets_rt);
  if (fgets_rt == NULL) {
    fclose(file_p);
    return 0;
  }
  start = 0;
  while (buff[start] == ' ' && buff[start] != '\0') {
    //      printf("buff[%3d] = %s\n", start, buff[start]);
    start++;
  }

  j = 0;
  start = 0;
  while (j < n) {
    if (rp_rd > 0) printf("\n Col %d\n", j);
    while (buff[start] == ' ') start++;

    end = start;
    while (buff[end] != ' ' && buff[end] != '\0') end++;

    int len = end - start;
    if (rp_rd > 0)
      printf("buff: start =%2d, len =%2d, end =%2d\n", start, len, end);
    if (len > 0) {
      term = (char*)malloc(sizeof(char) * len);
      for (k = 0; k < len; k++) term[k] = buff[start + k];
    } else {
      // JAJH: QY - LEN<=0 means that the end of the line has been reached
      printf("Comment reading file: no numbers on LB row.\n");
      return 0;
    }

    rl_v = atof(term);
    (*lb)[j] = rl_v;
    if (rp_rd > 0) printf("Read %s giving %g for (*lb)[%d]\n", term, rl_v, j);
    j++;
    if (j < n && buff[end] == '\0') {
      printf("Error reading file: not enough numbers on LB row.\n");
      return 1;
    }
    start = end + 1;
  }

  // Read upper bounds on columns
  fgets_rt = fgets(buff, buff_ln, file_p);
  if (rp_rd > 0) printf("UB: fgets return value = %s\n", fgets_rt);
  if (fgets_rt == NULL) {
    fclose(file_p);
    return 0;
  }
  start = 0;
  while (buff[start] == ' ' && buff[start] != '\0') {
    //      printf("buff[%3d] = %s\n", start, buff[start]);
    start++;
  }

  j = 0;
  start = 0;
  while (j < n) {
    if (rp_rd > 0) printf("\n Col %d\n", j);
    while (buff[start] == ' ') start++;

    end = start;
    while (buff[end] != ' ' && buff[end] != '\0') end++;

    int len = end - start;
    if (rp_rd > 0)
      printf("buff: start =%2d, len =%2d, end =%2d\n", start, len, end);
    if (len > 0) {
      term = (char*)malloc(sizeof(char) * len);
      for (k = 0; k < len; k++) term[k] = buff[start + k];
    } else {
      // JAJH: QY - LEN<=0 means that the end of the line has been reached
      printf("Comment reading file: no numbers on UB row.\n");
      return 0;
    }

    rl_v = atof(term);
    (*ub)[j] = rl_v;
    if (rp_rd > 0) printf("Read %s giving %g for (*ub)[%d]\n", term, rl_v, j);
    j++;
    if (j < n && buff[end] == '\0') {
      printf("Error reading file: not enough numbers on UB row.\n");
      return 1;
    }
    start = end + 1;
  }

  // Read integer information for columns
  fgets_rt = fgets(buff, buff_ln, file_p);
  if (rp_rd > 0) printf("IC: fgets return value = %s\n", fgets_rt);
  if (fgets_rt == NULL) {
    fclose(file_p);
    return 0;
  }
  start = 0;
  while (buff[start] == ' ' && buff[start] != '\0') {
    //      printf("buff[%3d] = %s\n", start, buff[start]);
    start++;
  }

  j = 0;
  start = 0;
  while (j < n) {
    if (rp_rd > 0) printf("\n Col %d\n", j);
    while (buff[start] == ' ') start++;

    end = start;
    while (buff[end] != ' ' && buff[end] != '\0') end++;

    int len = end - start;
    if (rp_rd > 0)
      printf("buff: start =%2d, len =%2d, end =%2d\n", start, len, end);
    if (len > 0) {
      term = (char*)malloc(sizeof(char) * len);
      for (k = 0; k < len; k++) term[k] = buff[start + k];
    } else {
      // JAJH: QY - LEN<=0 means that the end of the line has been reached
      printf("Comment reading file: no numbers on integerColumn row.\n");
      return 0;
    }

    i_v = atoi(term);
    (*integerColumn)[j] = i_v;
    if (rp_rd > 0)
      printf("Read %s giving %d for (*integerColumn)[%d]\n", term, i_v, j);
    j++;
    if (j < n && buff[end] == '\0') {
      printf("Error reading file: not enough numbers on integerColumn row.\n");
      return 1;
    }
    start = end + 1;
  }
  fclose(file_p);
  return 0;
}
