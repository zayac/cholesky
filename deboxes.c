#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "C4SNet.h"
#include "memfun.h"

/* A tile is an array of B*B doubles.
 * It is indexed column major.  */
typedef double *tile_t;

static bool verbose;
static bool compute;
static struct timespec begin;

static double delta_time(struct timespec t0, struct timespec t1)
{
    return t1.tv_sec - t0.tv_sec + t1.tv_nsec * 1e-9 - t0.tv_nsec * 1e-9;
}

/* Print an error message and terminate. */
static void pexit(const char *msg)
{
    perror(msg);
    exit(1);
}

/* Open a file or terminate. */
static FILE* xfopen(const char *fname, const char *mode)
{
    FILE *fp;
    if (!strcmp(fname, "-")) {
        fp = (*mode == 'r') ? stdin : stdout;
    }
    else if ((fp = fopen(fname, mode)) == NULL) {
        pexit(fname);
    }
    return fp;
}

static void xfclose(FILE *fp)
{
    if (fp != stdin && fp != stdout) fclose(fp);
}

/* Examine the optional environment variable DECOMPOSE for options. */
static void get_options(void)
{
    char        *env, *str, *save, *tok;

    compute = true;
    if ((env = getenv("DECOMPOSE")) != NULL) {
        str = strdup(env);
        save = NULL;
        for (tok = strtok_r(str, ",", &save);
             tok != NULL;
             tok = strtok_r(NULL, ",", &save)) {
            if (!strcmp(tok, "v")) {
                verbose = true;
            }
            else if (!strcmp(tok, "c")) {
                compute = false;
            }
        }
        free(str);
    }
}

/* Convert an array of chars to a nul-terminated string. */
static char* chars_to_string(c4snet_data_t *c4data)
{
  size_t size = C4SNetArraySize(c4data);
  char* str = SNetMemAlloc(size + 1);

  memcpy(str, C4SNetGetData(c4data), size);
  str[size] = '\0';
  return str;
}

/* Read the lower triangle of a matrix from file and copy to upper triangle. */
static void read_matrix(int n, double *A, const char *infile)
{
    FILE *fp = xfopen(infile, "r");
    int r, c;

    for (r = 0; r < n; ++r) {
        for (c = 0; c <= r; ++c) {
            if (fscanf(fp, " %lf", &A[r * n + c]) != 1) {
                fprintf(stderr, "Failed to read array from %s\n", infile);
                exit(1);
            }
            A[c * n + r] = A[r * n + c];
        }
    }
    xfclose(fp);
}

/* Setup the S-Net computation: read matrix data from file and distribute. */
void *start(void *hnd, c4snet_data_t *InFile, c4snet_data_t *OutFile, int N, int B)
{
    int         P, X;
    int         j, i;
    void       *result_data;
    char       *infile;
    c4snet_data_t *result;
    double     *array_data;

    get_options();
    /* Matrix size N must be a multiple of the block size B. */
    if (B < 1 || N < 1 || N % B) {
        fprintf(stderr, "%s: matrix size %d is not a multiple of the block size %d\n",
                __func__, N, B);
        exit(1);
    }
    /* P gives the number of tiles in either dimension. */
    P = N / B;
    /* X gives the number of tiles needed to construct the output. */
    X = (P * P + P) / 2;
    /* Get the input filename as a string. */
    infile = chars_to_string(InFile);
    C4SNetFree(InFile);
    if (verbose) {
        printf("%s: N=%d, B=%d, P=%d, X=%d, I=%s\n", __func__, N, B, P, X, infile);
    }

    result = C4SNetAlloc(CTYPE_long, P * P, &result_data);
    memset(result_data, 0, sizeof(long) * P * P);
    array_data = SNetMemAlloc(N * N * sizeof(double));
    if (compute) {
        read_matrix(N, array_data, infile);
    }
    free(infile);

    if (clock_gettime(CLOCK_REALTIME, &begin)) {
        pexit("clock_gettime");
    }

    /* Loop over the rows as j. */
    for (j = 0; j < P; ++j) {
        /* Loop over the columns up to the diagonal as i. */
        for (i = 0; i <= j; ++i) {
            void *tile_ptr;
            c4snet_data_t *tile = C4SNetAlloc(CTYPE_double, B * B, &tile_ptr);
            if (compute) {
                tile_t tile_data = (tile_t) tile_ptr;
                int A_y, T_y, A_x, T_x;
                for (A_y = j * B, T_y = 0; T_y < B; ++A_y, ++T_y) {
                    for (A_x = i * B, T_x = 0; T_x < B; ++A_x, ++T_x) {
                        tile_data[T_y + T_x * B] = array_data[A_y * N + A_x];
                    }
                }
            } else {
                memset(tile_ptr, 0, B * B * sizeof(double));
            }
            if (j == 0 && i == 0) {
                if (verbose) {
                    printf("%s: -> Fac, N=%d, B=%d, (1, , )\n", __func__, N, B);
                }
                C4SNetOut(hnd, 2, tile, N, B, 1);
            }
            else if (i == 0) {
                if (verbose) {
                    printf("%s: -> Tri2, (1, %d, )\n", __func__, j);
                }
                C4SNetOut(hnd, 3, tile, 1, j);
            }
            else {
                if (verbose) {
                    printf("%s: -> Sym3, (1, %d, %d)\n", __func__, j, i);
                }
                C4SNetOut(hnd, 4, tile, 1, j, i);
            }
        }
    }

    C4SNetOut(hnd, 1, OutFile, result, X);
    free(array_data);

    return hnd;
}

/* Do the InitialFactorization from input tile Fac to output tile Tri. */
void *facto(void *hnd, c4snet_data_t *Fac, int N, int B, int K)
{
    int         j;
    const int   P = N / B;
    c4snet_data_t *Tri;

    if (verbose) {
        printf("%s: %p, N=%d, B=%d, (%d, , )\n",
                __func__, Fac, N, B, K);
    }
    if (compute) {
        tile_t  A = (tile_t) C4SNetGetData(Fac);
        void    *tile_ptr;
        c4snet_data_t *tile = C4SNetAlloc(CTYPE_double, B * B, &tile_ptr);
        tile_t  L = (tile_t) tile_ptr;
        int     k_b, j_b, j_bb, i_b;

        memset(L, 0, B * B * sizeof(double));
        for (k_b = 0; k_b < B; ++k_b) {
            if (A[k_b + k_b * B] <= 0) {
                fprintf(stderr, "Not a symmetric positive definite (SPD) matrix\n");
                fprintf(stderr, "K = %d, k_b = %d, B = %d, A[k_b + k_b * B] = %f\n",
                        K, k_b, B, A[k_b + k_b * B]);
                exit(1);
            }
            L[k_b + k_b * B] = sqrt(A[k_b + k_b * B]);
            for (j_b = k_b + 1; j_b < B; ++j_b) {
                L[j_b + k_b * B] = A[j_b + k_b * B] / L[k_b + k_b * B];
            }
            for (j_bb = k_b + 1; j_bb < B; ++j_bb) {
                for (i_b = j_bb; i_b < B; ++i_b) {
                    A[i_b + j_bb * B] -= L[i_b + k_b * B] * L[j_bb + k_b * B];
                }
            }
        }
        Tri = tile;
        C4SNetFree(Fac);
    } else {
        Tri = Fac;
    }
    for (j = K; j < P; ++j) {
        if (verbose) {
            printf("%s: -> Tri, N=%d, B=%d, (%d, %d, )\n", __func__, N, B, K, j);
        }
        C4SNetOut(hnd, 1, C4SNetShallowCopy(Tri), N, B, K, j);
    }
    /* printf("%s: -> Out, N=%d, B=%d, (%d, %d, %d)\n",
           __func__, N, B, K, K - 1, K - 1); */
    C4SNetOut(hnd, 2, Tri, N, B, K, K - 1, K - 1);
    return hnd;
}

/* Do the TriangularSolve from input tiles Tri and Tri2. */
void *trisol(void *hnd, c4snet_data_t *Tri, c4snet_data_t *Tri2,
             int N, int B, int K, int J)
{
    c4snet_data_t *Out;
    c4snet_data_t *Sym;
    c4snet_data_t *Sym2;
    const int   P = N / B;
    int         i, j;

    if (verbose) {
        printf("%s: %p, %p, N=%d, B=%d, (%d, %d, )\n",
                __func__, Tri, Tri2, N, B, K, J);
    }
    if (compute) {
        tile_t  A = C4SNetGetData(Tri2);
        tile_t  Li = C4SNetGetData(Tri);
        void    *tile_ptr;
        c4snet_data_t *tile = C4SNetAlloc(CTYPE_double, B * B, &tile_ptr);
        tile_t  Lo = tile_ptr;
        int     k_b, i_b, j_b;

        for (k_b = 0; k_b < B; ++k_b) {
            for (i_b = 0; i_b < B; ++i_b) {
                Lo[i_b + k_b * B] = A[i_b + k_b * B] / Li[k_b + k_b * B];
            }
            for (j_b = k_b + 1; j_b < B; ++j_b) {
                for (i_b = 0; i_b < B; ++i_b) {
                    A[i_b + j_b * B] -= Li[j_b + k_b * B] * Lo[i_b + k_b * B];
                }
            }
        }
        Out = Sym = Sym2 = tile;
        C4SNetFree(Tri);
        C4SNetFree(Tri2);
    } else {
        Out = Sym = Sym2 = Tri;
        C4SNetFree(Tri2);
    }

    for (i = K + 1; i <= J; ++i) {
        if (verbose) {
            printf("%s: -> Sym(%d,%d)a, N=%d, B=%d, (%d, %d, %d)\n",
                   __func__, K, J, N, B, K, J, i);
        }
        C4SNetOut(hnd, 1, C4SNetShallowCopy(Sym), N, B, K, J, i);
        if (i == J && i != K) {
            if (verbose) {
                printf("%s: -> Sym2(%d,%d)a, (%d, %d, %d)\n",
                       __func__, K, J, K, J, i);
            }
            C4SNetOut(hnd, 2, C4SNetShallowCopy(Sym2), K, J, i);
        }
    }
    for (j = J; j < P; ++j) {
        if (j == J && j != K) {
            if (verbose) {
                printf("%s: -> Sym(%d,%d)b, N=%d, B=%d, (%d, %d, %d)\n",
                       __func__, K, J, N, B, K, j, K);
            }
            C4SNetOut(hnd, 1, C4SNetShallowCopy(Sym2), N, B, K, j, K);
        } else {
            if (verbose) {
                printf("%s: -> Sym2(%d,%d)c, (%d, %d, %d)\n",
                       __func__, K, J, K, j, J);
            }
            C4SNetOut(hnd, 2, C4SNetShallowCopy(Sym2), K, j, J);
        }
        if (j == K) {
            if (verbose) {
                printf("%s: -> Sym(%d,%d)d, N=%d, B=%d, (%d, %d, %d)\n",
                       __func__, K, J, N, B, K, j, K);
            }
            C4SNetOut(hnd, 1, C4SNetShallowCopy(Sym), N, B, K, j, K);
        }
    }
    C4SNetOut(hnd, 3, Out, N, B, K, J, K - 1);

    return hnd;
}

void *trinot(void *hnd, int K, int J)
{
    printf("%s: Invalid record (%d, %d, )\n", __func__, K, J);
    exit(1);
    return hnd;
}

/* Do the SymmetricRankUpdate step from input tiles Sym, Sym2, Sym3. */
void *rankup(void *hnd, c4snet_data_t *Sym, c4snet_data_t *Sym2,
             c4snet_data_t *Sym3, int N, int B, int K, int J, int I)
{
    c4snet_data_t *Fac = Sym3;
    c4snet_data_t *Tri2 = Sym3;
    
    if (verbose) {
        printf("%s: %p, %p, %p, N=%d, B=%d, (%d, %d, %d)\n",
                __func__, Sym, Sym2, Sym3, N, B, K, J, I);
    }
    if (compute) {
        tile_t      L1, L2, A = C4SNetGetData(Sym3);
        int         i_b, j_b, k_b;

        assert( K <= I && I <= J );

        if (I == J) {
            /* L2 = get_tile(ctx, K + 1, J, K); */
            L2 = C4SNetGetData(Sym);
            L1 = NULL;
        } else {
            /* L2 = get_tile(ctx, K + 1, I, K); */
            L2 = C4SNetGetData(Sym2);
            /* L1 = get_tile(ctx, K + 1, J, K); */
            L1 = C4SNetGetData(Sym);
        }

        for (j_b = 0; j_b < B; ++j_b) {
            for (k_b = 0; k_b < B; ++k_b) {
                double temp = -L2[j_b + k_b * B];
                if (I != J) {
                    for (i_b = 0; i_b < B; ++i_b) {
                        A[i_b + j_b * B] += temp * L1[i_b + k_b * B];
                    }
                } else {
                    for (i_b = j_b; i_b < B; ++i_b) {
                        A[i_b + j_b * B] += temp * L2[i_b + k_b * B];
                    }
                }
            }
        }
        C4SNetFree(Sym);
        C4SNetFree(Sym2);
    }
    if (K == J && J == I) {
        if (verbose) {
            printf("%s: -> Fac, N=%d, B=%d, (%d, , )\n", __func__, N, B, K + 1);
        }
        C4SNetOut(hnd, 1, Fac, N, B, K + 1);
    }
    else if (K == I) {
        if (verbose) {
            printf("%s: -> Tri2, (%d, %d, )\n", __func__, K + 1, J);
        }
        C4SNetOut(hnd, 2, Tri2, K + 1, J);
    }
    else {
        if (verbose) {
            printf("%s: -> Sym3, (%d, %d, %d)\n", __func__, K + 1, J, I);
        }
        C4SNetOut(hnd, 3, Sym3, K + 1, J, I);
    }
    return hnd;
}

void *trilog(void *hnd, int K, int J)
{
    printf("%s: (%d, %d, )\n", __func__, K, J);
    C4SNetOut(hnd, 1, K, J);
    return hnd;
}

void *symlog(void *hnd, int K, int J, int I)
{
    printf("%s: (%d, %d, %d)\n", __func__, K, J, I);
    C4SNetOut(hnd, 1, K, J, I);
    return hnd;
}

void *symnot(void *hnd, int K, int J, int I)
{
    printf("%s: Invalid record (%d, %d, %d)\n", __func__, K, J, I);
    exit(1);
    return hnd;
}

static void fnum(double a, char *buf)
{
    char        *s;

    sprintf(buf, " %f", a);
    if (strchr(buf, '.')) {
        for (s = buf + strlen(buf) - 1; s >= buf && *s == '0'; --s) {
            *s = '\0';
        }
        if (*s == '.') {
            *s = '\0';
        }
    }
}

static void write_matrix(c4snet_data_t *Result, int N, int B, c4snet_data_t *OutFile)
{
    const int   P = N / B;
    int         r, c, r_b, c_b;
    c4snet_data_t **result = C4SNetGetData(Result);
    char       *outf = C4SNetGetString(OutFile);
    FILE       *fp = xfopen(outf, "w");

    for (r = 0; r < P; ++r) {
        for (r_b = 0; r_b < B; ++r_b) {
            for (c = 0; c <= r; ++c) {
                c4snet_data_t *T = result[r + c * P];
                tile_t tmp = C4SNetGetData(T);
                for (c_b = 0; (r != c) ? (c_b < B) : (c_b <= r_b); ++c_b) {
                    double a = tmp[r_b + c_b * B];
                    char buf[40];
                    fnum(a, buf);
                    fputs(buf, fp);
                }
            }
            fprintf(fp, "\n");
        }
    }
    xfclose(fp);
    free(outf);
    for (r = 0; r < P; ++r) {
        for (c = 0; c <= r; ++c) {
            c4snet_data_t *Out = result[r + c * P];
            C4SNetFree(Out);
        }
    }
}

void *merge(void *hnd, c4snet_data_t *OutFile, c4snet_data_t *Result,
            int X, c4snet_data_t *Out, int N, int B, int K, int J, int I)
{
    const int   P = N / B;

    if (verbose) {
        printf("%s: %p, %p, X=%d, %p, N=%d, B=%d, (%d, %d, %d)\n",
                __func__, OutFile, Result, X, Out, N, B, K, J, I);
    }
    if (compute) {
        c4snet_data_t **result = C4SNetGetData(Result);
        assert(result[J + I * P] == NULL);
        result[J + I * P] = Out;
    }
    if (X > 1) {
        C4SNetOut(hnd, 1, OutFile, Result, X - 1);
    }
    else if (compute) {
        struct timespec end;
        if (clock_gettime(CLOCK_REALTIME, &end)) {
            pexit("clock_gettime");
        }
        printf("Time for size %d x %d : %lf sec\n", N, N, delta_time(begin, end));
        write_matrix(Result, N, B, OutFile);
        C4SNetOut(hnd, 2, C4SNetCreate(CTYPE_char, 5, "Done."));
        C4SNetFree(OutFile);
        C4SNetFree(Result);
    }
    return hnd;
}

