#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>

static void pexit(const char *msg)
{
    perror(msg);
    exit(1);
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

static void write_matrix(double *A, const int n, const char *outfile)
{
    int         j, k;
    char        buf[40];
    FILE       *fp = !strcmp(outfile, "-") ? stdout : fopen(outfile, "w");

    if (!fp) {
        pexit(outfile);
    }
    for (j = 0; j < n; ++j) {
        for (k = 0; k <= j; ++k) {
            fnum(A[k * n + j], buf);
            fputs(buf, fp);
        }
        fprintf(fp, "\n");
    }
    if (fp != stdout) {
        fclose(fp);
    }
}

static void symposdef(double *A, const int n, const char *trifile)
{
    double      *L = calloc(sizeof(double), n * n);
    double      *LT = calloc(sizeof(double), n * n);
    int          j, k, i;

    /* Construct lower triangle matrix L. */
    for (j = 0; j < n; ++j) {
        for (k = 0; k <= j; ++k) {
            if (k < j) {
                L[k * n + j] = (((j * k) / (j + 1.0)
                               / (k + 2.0) * 2.0) - 1.0) / n;
            } else if (k == j) {
                L[k * n + j] = 1.0;
            }
        }
    }
    /* Write lower triangle. */
    if (trifile) {
        write_matrix(L, n, trifile);
    }
    /* Construct transpose of L. */
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            LT[j * n + i] = L[i * n + j];
        }
    }
    /* Multiply L with its transpose and assign to A. */
    for (i = 0; i < n; ++i) {
        double *a = &A[i * n];
        double *l = &L[i * n];
        for (k = 0; k < n; ++k) {
            double *lt = &LT[k * n];
            for (j = 0; j < n; ++j) {
                a[j] += l[k] * lt[j];
            }
        }
    }

    free(L);
    free(LT);
}

static void usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [N] [-d] [-v] [-o outfile] [-t trifile]\n", prog);
    exit(1);
}

int main(int argc, char **argv)
{
    int         n = 16;
    int         argi = 0;
    const char *outfile = NULL;
    const char *trifile = NULL;
    double     *A;
    bool        debug = false;
    bool        verbose = false;

    while (++argi < argc) {
        if (argv[argi][0] == '-') {
            if (!strcmp(argv[argi], "-d")) {
                debug = true;
            }
            else if (!strcmp(argv[argi], "-v")) {
                verbose = true;
            }
            else if (!strcmp(argv[argi], "-o")) {
                outfile = argv[++argi];
            }
            else if (!strcmp(argv[argi], "-t")) {
                trifile = argv[++argi];
            }
            else {
                usage(*argv);
            }
        }
        else if (sscanf(argv[argi], "%d", &n) != 1 || n <= 0) {
            fprintf(stderr, "%s: Invalid array dimension %s\n",
                    *argv, argv[argi]);
            exit(1);
        }
    }

    A = calloc(sizeof(double), n * n);
    symposdef(A, n, trifile);
    if (outfile) {
        write_matrix(A, n, outfile);
    }
    else if (debug || verbose) {
        write_matrix(A, n, "-");
    }
    free(A);

    return 0;
}

