/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#ifndef __SUPERLU_DIST_PZBRIDGE
#define __SUPERLU_DIST_PZBRIDGE

#include "superlu_zdefs.h"

typedef struct {
    superlu_dist_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    zScalePermstruct_t ScalePermstruct;
    zLUstruct_t LUstruct;
    zSOLVEstruct_t SOLVEstruct;
    gridinfo_t grid;
} slu_handle_z;

#ifdef __cplusplus
extern "C" {
#endif

extern void pzbridge_init(int algo3d,
                          int_t m,
                          int_t n,
                          int_t nnz,
                          int_t *rowind,
                          int_t *colptr,
                          doublecomplex *nzval,
                          void **pyobj,
                          int argc,
                          char *argv[]);
extern void pzbridge_factor(void **pyobj);
extern void pzbridge_solve(void **pyobj, int nrhs, doublecomplex *b_global);
extern void pzbridge_free(void **pyobj);

extern void pzbridge_init2d(int_t m,
                            int_t n,
                            int_t nnz,
                            int_t *rowind,
                            int_t *colptr,
                            doublecomplex *nzval,
                            void **pyobj,
                            int argc,
                            char *argv[]);
extern void pzbridge_factor2d(void **pyobj);
extern void pzbridge_solve2d(void **pyobj, int nrhs, doublecomplex *b_global);
extern void pzbridge_free2d(void **pyobj);

int zcreate_matrix_from_csc(SuperMatrix *A,
                            int_t m,
                            int_t n,
                            int_t nnz,
                            int_t *rowind0,
                            int_t *colptr0,
                            doublecomplex *nzval0,
                            gridinfo_t *grid);

#ifdef __cplusplus
}
#endif

#endif
