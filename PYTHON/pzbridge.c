/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "pzbridge.h"

#include <stdio.h>
#include <stdlib.h>

static void parse_options_2d(superlu_dist_options_t *options,
                             int *nprow,
                             int *npcol,
                             int argc,
                             char *argv[]) {
    int lookahead = -1;
    int rowperm = -1;
    int colperm = -1;
    int ir = -1;
    int symbfact = -1;
    int sympattern = 0;
    int printstat = 0;
    int tinyp = 0;

    char **cpp;
    char c;

    for (cpp = argv + 1; *cpp; ++cpp) {
        if (**cpp == '-') {
            c = *(*cpp + 1);
            ++cpp;
            switch (c) {
                case 'r':
                    *nprow = atoi(*cpp);
                    break;
                case 'c':
                    *npcol = atoi(*cpp);
                    break;
                case 'l':
                    lookahead = atoi(*cpp);
                    break;
                case 'p':
                    rowperm = atoi(*cpp);
                    break;
                case 'q':
                    colperm = atoi(*cpp);
                    break;
                case 's':
                    symbfact = atoi(*cpp);
                    break;
                case 'i':
                    ir = atoi(*cpp);
                    break;
                case 'm':
                    sympattern = atoi(*cpp);
                    break;
                case 't':
                    printstat = atoi(*cpp);
                    break;
                case 'n':
                    tinyp = atoi(*cpp);
                    break;
                default:
                    break;
            }
        }
    }

    if (rowperm != -1) options->RowPerm = rowperm;
    if (colperm != -1) options->ColPerm = colperm;
    if (lookahead != -1) options->num_lookaheads = lookahead;
    if (ir != -1) options->IterRefine = ir;
    if (symbfact != -1) options->ParSymbFact = symbfact;
    if (sympattern == 1) options->SymPattern = YES;
    if (tinyp == 1) options->ReplaceTinyPivot = YES;
    if (printstat == 1) options->PrintStat = YES;
}

static void compute_local_rows_2d(const gridinfo_t *grid,
                                  int_t m,
                                  int_t *m_loc,
                                  int_t *fst_row) {
    int_t m_loc_fst;

    *m_loc = m / (grid->nprow * grid->npcol);
    m_loc_fst = *m_loc;
    *fst_row = grid->iam * m_loc_fst;

    if ((*m_loc * grid->nprow * grid->npcol) != m) {
        if (grid->iam == (grid->nprow * grid->npcol - 1)) {
            *m_loc = m - (*m_loc) * (grid->nprow * grid->npcol - 1);
        }
    }
}

void pzbridge_init(int algo3d,
                   int_t m,
                   int_t n,
                   int_t nnz,
                   int_t *rowind,
                   int_t *colptr,
                   doublecomplex *nzval,
                   void **pyobj,
                   int argc,
                   char *argv[]) {
    if (algo3d != 0) {
        ABORT("pzbridge currently supports 2D only (algo3d must be 0)");
    }
    pzbridge_init2d(m, n, nnz, rowind, colptr, nzval, pyobj, argc, argv);
}

void pzbridge_factor(void **pyobj) {
    pzbridge_factor2d(pyobj);
}

void pzbridge_solve(void **pyobj, int nrhs, doublecomplex *b_global) {
    pzbridge_solve2d(pyobj, nrhs, b_global);
}

void pzbridge_free(void **pyobj) {
    pzbridge_free2d(pyobj);
}

void pzbridge_init2d(int_t m,
                     int_t n,
                     int_t nnz,
                     int_t *rowind,
                     int_t *colptr,
                     doublecomplex *nzval,
                     void **pyobj,
                     int argc,
                     char *argv[]) {
    slu_handle_z *slu_obj = (slu_handle_z *) malloc(sizeof(slu_handle_z));
    int nprow = 1;
    int npcol = 1;
    int iam;

    set_default_options_dist(&(slu_obj->options));
    (slu_obj->options).PrintStat = NO;
    (slu_obj->options).Algo3d = NO;

    parse_options_2d(&(slu_obj->options), &nprow, &npcol, argc, argv);

    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &(slu_obj->grid));

    iam = (slu_obj->grid).iam;
    if ((iam >= nprow * npcol) || (iam == -1)) {
        *pyobj = (void *) slu_obj;
        return;
    }

    zcreate_matrix_from_csc(&(slu_obj->A), m, n, nnz, rowind, colptr, nzval, &(slu_obj->grid));

    zScalePermstructInit(m, n, &(slu_obj->ScalePermstruct));
    zLUstructInit(n, &(slu_obj->LUstruct));
    PStatInit(&(slu_obj->stat));

    *pyobj = (void *) slu_obj;
}

void pzbridge_factor2d(void **pyobj) {
    slu_handle_z *slu_obj = (slu_handle_z *) (*pyobj);
    int iam = (slu_obj->grid).iam;
    int nprow = (slu_obj->grid).nprow;
    int npcol = (slu_obj->grid).npcol;
    int info;
    int_t m = (slu_obj->A).nrow;
    int_t m_loc, fst_row;
    int nrhs = 1;
    double *berr;
    doublecomplex *b;
    int_t i;

    if ((iam >= nprow * npcol) || (iam == -1)) {
        *pyobj = (void *) slu_obj;
        return;
    }

    compute_local_rows_2d(&(slu_obj->grid), m, &m_loc, &fst_row);

    b = doublecomplexMalloc_dist((size_t) m_loc * (size_t) nrhs);
    if (!b) ABORT("Malloc fails for rhs[].");
    for (i = 0; i < (int_t) m_loc * (int_t) nrhs; ++i) {
        b[i].r = 0.0;
        b[i].i = 0.0;
    }

    berr = doubleMalloc_dist((size_t) nrhs);
    if (!berr) ABORT("Malloc fails for berr[].");

    pzgssvx(&(slu_obj->options),
            &(slu_obj->A),
            &(slu_obj->ScalePermstruct),
            b,
            m_loc,
            nrhs,
            &(slu_obj->grid),
            &(slu_obj->LUstruct),
            &(slu_obj->SOLVEstruct),
            berr,
            &(slu_obj->stat),
            &info);

    if ((slu_obj->options).PrintStat == YES) {
        PStatPrint(&(slu_obj->options), &(slu_obj->stat), &(slu_obj->grid));
    }

    (slu_obj->options).Fact = FACTORED;

    SUPERLU_FREE(b);
    SUPERLU_FREE(berr);
    *pyobj = (void *) slu_obj;
}

void pzbridge_solve2d(void **pyobj, int nrhs, doublecomplex *b_global) {
    slu_handle_z *slu_obj = (slu_handle_z *) (*pyobj);
    int iam = (slu_obj->grid).iam;
    int nprow = (slu_obj->grid).nprow;
    int npcol = (slu_obj->grid).npcol;
    int procs = nprow * npcol;
    int info;
    int ldb;
    int_t m = (slu_obj->A).nrow;
    int_t m_loc, fst_row;
    double *berr;
    doublecomplex *b;
    doublecomplex *btmp = NULL;
    doublecomplex *b_global_tmp = NULL;

    MPI_Bcast(&nrhs, 1, MPI_INT, 0, (slu_obj->grid).comm);

    if ((iam >= nprow * npcol) || (iam == -1)) {
        *pyobj = (void *) slu_obj;
        return;
    }

    compute_local_rows_2d(&(slu_obj->grid), m, &m_loc, &fst_row);

    b = doublecomplexMalloc_dist((size_t) m_loc * (size_t) nrhs);
    if (!b) ABORT("Malloc fails for rhs[].");

    if (procs > 1) {
        int *counts = NULL;
        int *displs = NULL;
        int m_loc1 = (int) (m_loc * nrhs);
        int fst_row1 = (int) (fst_row * nrhs);
        int i, j;

        btmp = doublecomplexMalloc_dist((size_t) m_loc * (size_t) nrhs);
        if (!btmp) ABORT("Malloc fails for btmp[].");

        if (iam == 0) {
            counts = (int *) intMalloc_dist(procs);
            displs = (int *) intMalloc_dist(procs);
            b_global_tmp = doublecomplexMalloc_dist((size_t) m * (size_t) nrhs);
            if (!b_global_tmp) ABORT("Malloc fails for b_global_tmp[].");

            for (j = 0; j < nrhs; ++j) {
                for (i = 0; i < m; ++i) {
                    b_global_tmp[j + i * nrhs] = b_global[(size_t) j * (size_t) m + (size_t) i];
                }
            }
        }

        MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, (slu_obj->grid).comm);
        MPI_Gather(&fst_row1, 1, MPI_INT, displs, 1, MPI_INT, 0, (slu_obj->grid).comm);

        MPI_Scatterv(b_global_tmp,
                     counts,
                     displs,
                     SuperLU_MPI_DOUBLE_COMPLEX,
                     btmp,
                     m_loc1,
                     SuperLU_MPI_DOUBLE_COMPLEX,
                     0,
                     (slu_obj->grid).comm);

        for (j = 0; j < nrhs; ++j) {
            int_t irow;
            for (irow = 0; irow < m_loc; ++irow) {
                b[(size_t) j * (size_t) m_loc + (size_t) irow] =
                    btmp[(size_t) j + (size_t) irow * (size_t) nrhs];
            }
        }

        if (iam == 0) {
            SUPERLU_FREE(counts);
            SUPERLU_FREE(displs);
            SUPERLU_FREE(b_global_tmp);
        }
        SUPERLU_FREE(btmp);
    } else {
        int j;
        for (j = 0; j < nrhs; ++j) {
            int_t irow;
            for (irow = 0; irow < m_loc; ++irow) {
                int_t row = fst_row + irow;
                b[(size_t) j * (size_t) m_loc + (size_t) irow] =
                    b_global[(size_t) j * (size_t) m + (size_t) row];
            }
        }
    }

    ldb = (int) m_loc;

    berr = doubleMalloc_dist((size_t) nrhs);
    if (!berr) ABORT("Malloc fails for berr[].");

    pzgssvx(&(slu_obj->options),
            &(slu_obj->A),
            &(slu_obj->ScalePermstruct),
            b,
            ldb,
            nrhs,
            &(slu_obj->grid),
            &(slu_obj->LUstruct),
            &(slu_obj->SOLVEstruct),
            berr,
            &(slu_obj->stat),
            &info);

    if (procs > 1) {
        int *counts = NULL;
        int *displs = NULL;
        int m_loc1 = (int) (m_loc * nrhs);
        int fst_row1 = (int) (fst_row * nrhs);
        int i, j;

        btmp = doublecomplexMalloc_dist((size_t) m_loc * (size_t) nrhs);
        if (!btmp) ABORT("Malloc fails for btmp[].");

        for (j = 0; j < nrhs; ++j) {
            int_t irow;
            for (irow = 0; irow < m_loc; ++irow) {
                btmp[(size_t) j + (size_t) irow * (size_t) nrhs] =
                    b[(size_t) j * (size_t) m_loc + (size_t) irow];
            }
        }

        if (iam == 0) {
            counts = (int *) intMalloc_dist(procs);
            displs = (int *) intMalloc_dist(procs);
            b_global_tmp = doublecomplexMalloc_dist((size_t) m * (size_t) nrhs);
            if (!b_global_tmp) ABORT("Malloc fails for b_global_tmp[].");
        }

        MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, (slu_obj->grid).comm);
        MPI_Gather(&fst_row1, 1, MPI_INT, displs, 1, MPI_INT, 0, (slu_obj->grid).comm);

        MPI_Gatherv(btmp,
                    m_loc1,
                    SuperLU_MPI_DOUBLE_COMPLEX,
                    b_global_tmp,
                    counts,
                    displs,
                    SuperLU_MPI_DOUBLE_COMPLEX,
                    0,
                    (slu_obj->grid).comm);

        if (iam == 0) {
            for (j = 0; j < nrhs; ++j) {
                for (i = 0; i < m; ++i) {
                    b_global[(size_t) j * (size_t) m + (size_t) i] =
                        b_global_tmp[(size_t) j + (size_t) i * (size_t) nrhs];
                }
            }
            SUPERLU_FREE(counts);
            SUPERLU_FREE(displs);
            SUPERLU_FREE(b_global_tmp);
        }

        SUPERLU_FREE(btmp);
    } else {
        int j;
        for (j = 0; j < nrhs; ++j) {
            int_t irow;
            for (irow = 0; irow < m_loc; ++irow) {
                int_t row = fst_row + irow;
                b_global[(size_t) j * (size_t) m + (size_t) row] =
                    b[(size_t) j * (size_t) m_loc + (size_t) irow];
            }
        }
    }

    if ((slu_obj->options).PrintStat == YES) {
        PStatPrint(&(slu_obj->options), &(slu_obj->stat), &(slu_obj->grid));
    }

    SUPERLU_FREE(b);
    SUPERLU_FREE(berr);
    *pyobj = (void *) slu_obj;
}

void pzbridge_free2d(void **pyobj) {
    slu_handle_z *slu_obj = (slu_handle_z *) (*pyobj);
    int iam = (slu_obj->grid).iam;
    int m = (slu_obj->A).nrow;

    Destroy_CompRowLoc_Matrix_dist(&(slu_obj->A));
    zScalePermstructFree(&(slu_obj->ScalePermstruct));
    zDestroy_LU(m, &(slu_obj->grid), &(slu_obj->LUstruct));
    zLUstructFree(&(slu_obj->LUstruct));
    zSolveFinalize(&(slu_obj->options), &(slu_obj->SOLVEstruct));

    superlu_gridexit(&(slu_obj->grid));
    if (iam != -1) {
        PStatFree(&(slu_obj->stat));
    }

    *pyobj = (void *) slu_obj;
}

int zcreate_matrix_from_csc(SuperMatrix *A,
                            int_t m,
                            int_t n,
                            int_t nnz,
                            int_t *rowind0,
                            int_t *colptr0,
                            doublecomplex *nzval0,
                            gridinfo_t *grid) {
    SuperMatrix GA;
    int_t *rowind;
    int_t *colptr;
    int_t *colind;
    int_t *rowptr;
    int_t m_loc;
    int_t m_loc_fst;
    int_t fst_row;
    int_t nnz_loc;
    int_t *marker;
    doublecomplex *nzval;
    doublecomplex *nzval_loc;

    int_t chunk = 2000000000;
    int_t Nchunk;
    int_t remainder;
    int count;

    int iam = grid->iam;
    int_t i, j, row, relpos;

    if (!iam) {
        zallocateA_dist(n, nnz, &nzval, &rowind, &colptr);
        for (i = 0; i < nnz; ++i) {
            nzval[i] = nzval0[i];
            rowind[i] = rowind0[i];
        }
        for (i = 0; i < n + 1; ++i) {
            colptr[i] = colptr0[i];
        }

        MPI_Bcast(&m, 1, mpi_int_t, 0, grid->comm);
        MPI_Bcast(&n, 1, mpi_int_t, 0, grid->comm);
        MPI_Bcast(&nnz, 1, mpi_int_t, 0, grid->comm);

        Nchunk = CEILING(nnz, chunk);
        remainder = nnz % chunk;
        MPI_Bcast(&Nchunk, 1, mpi_int_t, 0, grid->comm);
        MPI_Bcast(&remainder, 1, mpi_int_t, 0, grid->comm);

        for (i = 0; i < Nchunk; ++i) {
            int_t idx = i * chunk;
            count = (i == Nchunk - 1 && remainder != 0) ? (int) remainder : (int) chunk;
            MPI_Bcast(&nzval[idx], count, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid->comm);
            MPI_Bcast(&rowind[idx], count, mpi_int_t, 0, grid->comm);
        }

        MPI_Bcast(colptr, n + 1, mpi_int_t, 0, grid->comm);
    } else {
        MPI_Bcast(&m, 1, mpi_int_t, 0, grid->comm);
        MPI_Bcast(&n, 1, mpi_int_t, 0, grid->comm);
        MPI_Bcast(&nnz, 1, mpi_int_t, 0, grid->comm);

        MPI_Bcast(&Nchunk, 1, mpi_int_t, 0, grid->comm);
        MPI_Bcast(&remainder, 1, mpi_int_t, 0, grid->comm);

        zallocateA_dist(n, nnz, &nzval, &rowind, &colptr);

        for (i = 0; i < Nchunk; ++i) {
            int_t idx = i * chunk;
            count = (i == Nchunk - 1 && remainder != 0) ? (int) remainder : (int) chunk;
            MPI_Bcast(&nzval[idx], count, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid->comm);
            MPI_Bcast(&rowind[idx], count, mpi_int_t, 0, grid->comm);
        }

        MPI_Bcast(colptr, n + 1, mpi_int_t, 0, grid->comm);
    }

    m_loc = m / (grid->nprow * grid->npcol);
    m_loc_fst = m_loc;
    if ((m_loc * grid->nprow * grid->npcol) != m) {
        if (iam == (grid->nprow * grid->npcol - 1)) {
            m_loc = m - m_loc * (grid->nprow * grid->npcol - 1);
        }
    }

    zCreate_CompCol_Matrix_dist(&GA, m, n, nnz, nzval, rowind, colptr, SLU_NC, SLU_Z, SLU_GE);

    rowptr = (int_t *) intMalloc_dist(m_loc + 1);
    marker = (int_t *) intCalloc_dist(n);

    for (i = 0; i < n; ++i) {
        for (j = colptr[i]; j < colptr[i + 1]; ++j) {
            ++marker[rowind[j]];
        }
    }

    rowptr[0] = 0;
    fst_row = iam * m_loc_fst;
    for (j = 0; j < m_loc; ++j) {
        row = fst_row + j;
        rowptr[j + 1] = rowptr[j] + marker[row];
        marker[j] = rowptr[j];
    }
    nnz_loc = rowptr[m_loc];

    nzval_loc = doublecomplexMalloc_dist(nnz_loc);
    colind = (int_t *) intMalloc_dist(nnz_loc);

    for (i = 0; i < n; ++i) {
        for (j = colptr[i]; j < colptr[i + 1]; ++j) {
            row = rowind[j];
            if ((row >= fst_row) && (row < fst_row + m_loc)) {
                row = row - fst_row;
                relpos = marker[row];
                colind[relpos] = i;
                nzval_loc[relpos] = nzval[j];
                ++marker[row];
            }
        }
    }

    Destroy_CompCol_Matrix_dist(&GA);

    zCreate_CompRowLoc_Matrix_dist(A,
                                   m,
                                   n,
                                   nnz_loc,
                                   m_loc,
                                   fst_row,
                                   nzval_loc,
                                   colind,
                                   rowptr,
                                   SLU_NR_loc,
                                   SLU_Z,
                                   SLU_GE);

    SUPERLU_FREE(marker);
    return 0;
}
