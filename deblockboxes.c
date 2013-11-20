#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "C4SNet.h"
#include "memfun.h"

typedef double *tile;

static struct timespec begin;
static char* outfile;

static double delta_time(struct timespec t0, struct timespec t1)
{
    return t1.tv_sec - t0.tv_sec + t1.tv_nsec * 1e-9 - t0.tv_nsec * 1e-9;
}

// compute_s1

void * compute_s1 (void *hnd, c4snet_data_t * atiles, c4snet_data_t * ltiles,
    int bs, int p, int k)
{
    int k_b, j_b, j, i_b, j_bb;
    tile * a_tiles = *((tile **) C4SNetGetData (atiles));
    tile * l_tiles = *((tile **) C4SNetGetData (ltiles));
    for (k_b = 0; k_b < bs; k_b++) {
            l_tiles[k * p + k][k_b * bs + k_b] =
                sqrt (a_tiles[k * p + k][k_b * bs + k_b]);

        for (j_b = k_b + 1; j_b < bs; j_b++)
            l_tiles[k * p + k][k_b * bs + j_b] =
                a_tiles[k * p + k][k_b * bs + j_b] 
              / l_tiles[k * p + k][k_b * bs + k_b];

        for (j_bb = k_b + 1; j_bb < bs; j_bb++)
            for (i_b = j_bb; i_b < bs; i_b++)
                a_tiles[k * p + k][i_b * bs + j_bb] =
                    a_tiles[k * p + k][i_b * bs + j_bb]
                  - l_tiles[k * p + k][k_b * bs + i_b]
                  * l_tiles[k * p + k][k_b * bs + j_bb];
    }

    if (k + 1 < p)
        for (j = k + 1; j < p; j++)
            C4SNetOut (hnd, 1, atiles, ltiles, bs, p, k, j);
    else {
        int i, j;
        C4SNetOut (hnd, 2, ltiles, bs, p);
        /* deallocate a_tiles array and remove atiles field.  */
        for (i = 0; i < p; i++)
            for (j = 0; j <= i; j++)
                free (a_tiles[j * p + i]);
        free (a_tiles);
        C4SNetFree (atiles);
    }
    return hnd;
}

// compute_s2

void * solve_s2 (void *hnd, c4snet_data_t * atiles, c4snet_data_t * ltiles,
    int bs, int p, int k, int j)
{
    int k_b, j_b, i_b;
    tile * a_tiles = *((tile **) C4SNetGetData (atiles));
    tile * l_tiles = *((tile **) C4SNetGetData (ltiles));
    memset (l_tiles[k * p + j], 0, sizeof (double) * bs * bs);
    assert (j != k);
    for (k_b = 0; k_b < bs; k_b++) {
        for (i_b = 0; i_b < bs; i_b++)
            l_tiles[k * p + j][k_b * bs + i_b] =
                        a_tiles[k * p + j][k_b * bs + i_b]
                      / l_tiles[k * p + k][k_b * bs + k_b];
        for (j_b = k_b + 1; j_b < bs; j_b++)
            for (i_b = 0; i_b < bs; i_b++)
                a_tiles[k * p + j][j_b * bs + i_b] =
                        a_tiles[k * p + j][j_b * bs + i_b]
                      - l_tiles[k * p + k][k_b * bs + j_b] 
                      * l_tiles[k * p + j][k_b * bs + i_b];
    }
    C4SNetOut (hnd, 1, atiles, ltiles, bs, p, k);
    return hnd;
}

// compute_s3

void * distribute (void *hnd, c4snet_data_t * fatiles, c4snet_data_t * fltiles,
    int bs, int p, int k)
{
    int i, j;
    for (j = k + 1; j < p; j++)
        for (i = k + 1; i <= j; i++)
            C4SNetOut (hnd, 1, fatiles, fltiles, bs, p, k, i * p + j);
    return hnd;
}

void * update (void *hnd, c4snet_data_t * fatiles, c4snet_data_t * fltiles,
    int bs, int p, int k, int index)
{
    tile * a_tiles = *((tile **) C4SNetGetData (fatiles));
    tile * l_tiles = *((tile **) C4SNetGetData (fltiles));
    int j_b, k_b, i_b;
    int i = index % p, j = index / p;
    tile l1_block = NULL, l2_block = NULL;
    // Diagonal tile
    if (i == j)
        l2_block = l_tiles[k * p + j];
    else {
        l1_block = l_tiles[k * p + i];
        l2_block = l_tiles[k * p + j];
    }
    for (j_b = 0; j_b < bs; j_b++)
        for (k_b = 0; k_b < bs; k_b++) {
            double temp = -1 * l2_block[k_b * bs + j_b];
            if (i != j)
                for (i_b = 0; i_b < bs; i_b++)
                    a_tiles[j * p + i][j_b * bs + i_b] =
                        a_tiles[j * p + i][j_b * bs + i_b]
                      + temp * l1_block[k_b * bs + i_b];
            else
                for (i_b = j_b; i_b < bs; i_b++)
                    a_tiles[j * p + i][j_b * bs + i_b] =
                        a_tiles[j * p + i][j_b * bs + i_b]
                      + temp * l2_block[k_b * bs + i_b];
        }
    C4SNetOut (hnd, 1, fatiles, fltiles, bs, p, k);
    return hnd;
}

// decompose

/* Convert an array of chars to a nul-terminated string. */
static char* chars_to_string(c4snet_data_t *c4data)
{
  size_t size = C4SNetArraySize(c4data);
  char* str = SNetMemAlloc(size + 1);

  memcpy(str, C4SNetGetData(c4data), size);
  str[size] = '\0';
  return str;
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

void * decompose (void *hnd, c4snet_data_t *InFile, c4snet_data_t *OutFile,
    int a_size, int bs)
{
    int p, i, j;
    double *array;
    tile *atiles, *ltiles;
    char *infile;

    if (bs <= 0) {
        fprintf (stderr, "A block size must be greater than 0\n");
        exit (1);
    }
    if (a_size <= 0) {
        fprintf (stderr, "The dimension of matrix must be greater than 0\n");
        exit (1);
    }
    if (a_size % bs) {
        fprintf(stderr, 
            "matrix size %d is not a multiple of the block size %d\n",
            a_size, bs);
        exit(1);
    }

    p = a_size / bs;

    /* Get the input filename as a string. */
    infile = chars_to_string(InFile);
    outfile = chars_to_string(OutFile);
    C4SNetFree(InFile);
    C4SNetFree(OutFile);

    array = SNetMemAlloc(a_size * a_size * sizeof(double));
    read_matrix(a_size, array, infile);
    free(infile);

    if (clock_gettime(CLOCK_REALTIME, &begin)) {
        pexit("clock_gettime");
    }

    atiles = (tile *) malloc (sizeof (tile) * p * p);
    ltiles = (tile *) malloc (sizeof (tile) * p * p);
    memset (atiles, 0, sizeof (tile) * p * p);
    memset (ltiles, 0, sizeof (tile) * p * p);
    for (i = 0; i < p; i++)
        for (j = 0; j <= i; j++) {
            atiles[j * p + i] = (double *) malloc (sizeof (double) * bs * bs);
            ltiles[j * p + i] = (double *) malloc (sizeof (double) * bs * bs);
            int ai, aj, ti, tj;
            for (ai = i * bs, ti = 0; ti < bs; ai++, ti++)
                for (aj = j * bs, tj = 0; tj < bs; aj++, tj++)
                    atiles[j * p + i][tj * bs + ti] = array[aj * a_size + ai];
        }
    C4SNetOut (hnd, 1, C4SNetCreate (CTYPE_char, sizeof (void *), &atiles),
               C4SNetCreate (CTYPE_char, sizeof (void *), &ltiles), bs, p, 0);
    free(array);
    return hnd;
}

// merger

void * gen_counter (void *hnd)
{
    C4SNetOut (hnd, 1, 0);
    return hnd;
}

void * merge (void *hnd, c4snet_data_t * atiles, c4snet_data_t * ltiles,
    int counter, int bs, int p, int k)
{
    if (++counter < p - k - 1)
        C4SNetOut (hnd, 1, counter);
    else
        C4SNetOut (hnd, 2, atiles, ltiles, bs, p, k);
    return hnd;
}

static void write_matrix(tile *tiles, int p, int bs)
{
    FILE *fp = xfopen(outfile, "w");
    int size = p * bs;
    int i, j;
    double *array = (double *) malloc (sizeof (double) * size * size);
    memset (array, 0.0, sizeof (double) * sizeof (size * size));
    for (i = 0; i < p; i++) {
        for (j = 0; j <= i; j++) {
            int k;
            for (k = 0; k < bs * bs; k++) {
                int px = j * bs + (k / bs);
                int py = i * bs + (k % bs);
                array[py * size + px] = tiles[j * p + i][k];
            }
            free (tiles[j * p + i]);
        }
    }
    for (i = 0; i < size; i++) {
        for(j = 0; j <= i; j++)
            fprintf(fp, "%lf ", array[i * size + j]);
        fprintf(fp, "\n");
    }

    xfclose(fp);
    free (outfile);
    free (array);
    free (tiles);
}

void *
finalize (void *hnd, c4snet_data_t * fltiles, int bs, int p)
{
    struct timespec end;
    tile * tiles = *((tile **) C4SNetGetData (fltiles));

    if (clock_gettime(CLOCK_REALTIME, &end)) {
        pexit("clock_gettime");
    }
    printf("Time for size %d x %d : %lf sec\n", bs * p, bs * p, delta_time(begin, end));
    
    write_matrix(tiles, p, bs);

    C4SNetOut (hnd, 1, C4SNetCreate (CTYPE_char, 5, "Done."));
    C4SNetFree (fltiles);
    return hnd;
}

void * sync (void *hnd, c4snet_data_t * fatiles, c4snet_data_t * fltiles,
    int counter, int bs, int p, int k)
{
    if (++counter < (p - k) * (p - k - 1) / 2)
        C4SNetOut (hnd, 1, counter);
    else
        C4SNetOut (hnd, 2, fatiles, fltiles, bs, p, k);
    return hnd;
}
