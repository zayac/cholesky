#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

#define SQR(x)          ((x) * (x))
#define TILI(k,j,i)     ((i) + ((j) * p) + ((k) * p * p))
#define TILE(k,j,i)     ctx->T[TILI(k,j,i)]
#define N               ctx->n
#define B               ctx->b
#define P               ctx->p

typedef double *tile_t;

typedef struct ctx {
    tile_t     *T;
    int         n;
    int         b;
    int         p;
} ctx_t;

static bool             opt_debug;
static bool             opt_verbose;
static FILE*            opt_graph;
static const char*      color;

static void pexit(const char *msg)
{
    perror(msg);
    exit(1);
}

static struct timespec begin;

static double delta_time(struct timespec t0, struct timespec t1)
{
    return t1.tv_sec - t0.tv_sec + t1.tv_nsec * 1e-9 - t0.tv_nsec * 1e-9;
}

static FILE* xfopen(const char *fname, const char *mode)
{
    FILE *fp = fopen(fname, mode);
    if (!fp) pexit(fname);
    return fp;
}

static void edge(int p, int src, int end)
{
    if (src < p * p) {
        fprintf(opt_graph, "    struct1:%d -> %d;\n", src, end);
    } else {
        fprintf(opt_graph, "    %d -> %d;\n", src, end);
    }
}

static void store_tile(ctx_t *ctx, int k, int j, int i, tile_t T)
{
    const int   p = P;

    if (T && opt_debug) {
        printf("put %2d, %2d, %2d\n", k, j, i);
    }
    assert((unsigned) k <= p);
    assert((unsigned) j <  p);
    assert((unsigned) i <  p);
    assert(!TILE(k,j,i) || !T);
    TILE(k,j,i) = T;

    if (opt_graph && k && T) {
        char *group = ((k == j + 1 && j == i) ||
                       (k == j && j == i + 1) ||
                       (k == j && j == i)) ? ",group=left" : "";
        fprintf(opt_graph, "    %d [label=\"%d,%d,%d\"%s%s];\n",
                TILI(k,j,i), k, j, i, color, group);
    }
}

static tile_t get_tile(ctx_t *ctx, int k, int j, int i)
{
    const int   p = P;
    tile_t      T;

    if (opt_debug) {
        printf("    %2d, %2d, %2d\n", k, j, i);
    }
    assert((unsigned) k <= p);
    assert((unsigned) j <  p);
    assert((unsigned) i <  p);
    T = TILE(k,j,i);
    assert(T);

    return T;
}

static tile_t new_tile(ctx_t *ctx)
{
    tile_t      T = (tile_t) malloc(B * B * sizeof(double));
    return T;
}

static ctx_t* make_context(const int n, const int b, const int p,
                           const double *A)
{
    ctx_t       *ctx = malloc(sizeof(ctx_t));
    int          j, i;

    N = n;
    B = b;
    P = p;
    ctx->T = (tile_t *) calloc(1, (p + 1) * SQR(p) * sizeof(tile_t));

    if (opt_graph) {
        fprintf(opt_graph, "subgraph top {\n");
        color = ",group=\"top\",fillcolor=\"white\"";
        fprintf(opt_graph,
                "struct1 [shape=record,group=left,"
                "height=0.3,margin=0.15,fillcolor=\"white\",label=\"");
    }

    for (i = 0; i < p; ++i) {
        for (j = 0; j <= i; ++j) {
            tile_t T = new_tile(ctx);
            int A_i, T_i, A_j, T_j;
            for (A_i = i * b, T_i = 0; T_i < b; ++A_i, ++T_i) {
                for (A_j = j * b, T_j = 0; T_j < b; ++A_j, ++T_j) {
                    T[T_i + T_j * b] = A[A_i * N + A_j];
                }
            }
            store_tile(ctx, 0, i, j, T);
            if (opt_graph) {
                int loc = TILI(0,i,j);
                char sep = loc ? '|' : ' ';
                fprintf(opt_graph, "%c <%d> 0,%d,%d ", sep, loc, i, j);
            }
        }
    }

    if (opt_graph) {
        fprintf(opt_graph, "\"];\n");
        color = "";
        fprintf(opt_graph, "}\n");
    }

    return ctx;
}

static void free_context(ctx_t *ctx)
{
    const int   p = P;
    int         i, j, k;
    int         a = 0, b = 0;

    for (k = 0; k <= p; ++k) {
        for (j = 0; j < p; ++j) {
            for (i = 0; i < p; ++i) {
                b = i + (j * p) + (k * p * p);
                if (b) {
                    assert(!a || b > a);
                    free(ctx->T[b]);
                    ctx->T[b] = NULL;
                    a = b;
                }
            }
        }
    }

    free(ctx->T);
    free(ctx);
}

/* Compute tile k+1,j,i from k,j,i and k+1,j,k and k+1,i,k, where k<i<=j<p. */
static void S3_compute(ctx_t *ctx, const int k, const int j, const int i)
{
    tile_t      A = get_tile(ctx, k, j, i);
    tile_t      L2;
    tile_t      L1;
    const int   b = B;
    int         i_b, j_b, k_b;

    assert( k < i && i <= j );

    if (i == j) {
        L2 = get_tile(ctx, k + 1, j, k);
        L1 = NULL;
    } else {
        L2 = get_tile(ctx, k + 1, i, k);
        L1 = get_tile(ctx, k + 1, j, k);
    }
 
   for (j_b = 0; j_b < b; ++j_b) {
        for (k_b = 0; k_b < b; ++k_b) {
            double temp = -L2[j_b + k_b * b];
            if (i != j) {
                for (i_b = 0; i_b < b; ++i_b) {
                    A[i_b + j_b * b] += temp * L1[i_b + k_b * b];
                }
            } else {
                for (i_b = j_b; i_b < b; ++i_b) {
                    A[i_b + j_b * b] += temp * L2[i_b + k_b * b];
                }
            }
        }
    }

    store_tile(ctx, k, j, i, NULL);
    store_tile(ctx, k + 1, j, i, A);

    if (opt_graph) {
        const int p = P;
        edge(p, TILI(k,j,i), TILI(k+1,j,i));
        edge(p, TILI(k+1,j,k), TILI(k+1,j,i));
        if (i != j) {
            edge(p, TILI(k+1,i,k), TILI(k+1,j,i));
        }
    }
}

static void kji_compute(ctx_t *ctx, const int k, const int j)
{
    int         i;
	#pragma omp parallel for
    for (i = k + 1; i <= j; ++i) {
        S3_compute(ctx, k, j, i);
    }
}

/* Compute tile k+1,j,k from tiles k,j,k and k+1,k,k, where k<j<p. */
static void S2_compute(ctx_t *ctx, const int k, const int j)
{
    const int   b = B;
    tile_t      A = get_tile(ctx, k, j, k);
    tile_t      Li = get_tile(ctx, k + 1, k, k);
    tile_t      Lo = new_tile(ctx);
    int         k_b, i_b, j_b;

    for (k_b = 0; k_b < b; ++k_b) {
        for (i_b = 0; i_b < b; ++i_b) {
            Lo[i_b + k_b * b] = A[i_b + k_b * b] / Li[k_b + k_b * b];
        }
        for (j_b = k_b + 1; j_b < b; ++j_b) {
            for (i_b = 0; i_b < b; ++i_b) {
                A[i_b + j_b * b] -= Li[j_b + k_b * b] * Lo[i_b + k_b * b];
            }
        }
    }

    color = ",fillcolor=\"lightblue\"";
    store_tile(ctx, k + 1, j, k, Lo);
    color = "";

    if (opt_graph) {
        const int p = P;
        edge(p, TILI(k,j,k), TILI(k+1,j,k));
        edge(p, TILI(k+1,k,k), TILI(k+1,j,k));
    }

    kji_compute(ctx, k, j);
}

static void kj_compute(ctx_t *ctx, const int k)
{
    const int   p = P;
    int         j;
    for (j = k + 1; j < p; ++j) {
        S2_compute(ctx, k, j);
    }
}

/* Compute tile k+1,k,k from tile k,k,k */
static void S1_compute(ctx_t *ctx, const int k)
{
    tile_t      A = get_tile(ctx, k, k, k);
    tile_t      L = new_tile(ctx);
    const int   b = B;
    int         k_b, j_b, j_bb, i_b;

    memset(L, 0, b * b * sizeof(double));
    for (k_b = 0; k_b < b; ++k_b) {
        if (A[k_b + k_b * b] <= 0) {
            fprintf(stderr, "Not a symmetric positive definite (SPD) matrix\n");
            fprintf(stderr, "k = %d, k_b = %d, b = %d, A[k_b + k_b * b] = %f\n",
                    k, k_b, b, A[k_b + k_b * b]);
            exit(1);
        }
        L[k_b + k_b * b] = sqrt(A[k_b + k_b * b]);
        for (j_b = k_b + 1; j_b < b; ++j_b) {
            L[j_b + k_b * b] = A[j_b + k_b * b] / L[k_b + k_b * b];
        }
        for (j_bb = k_b + 1; j_bb < b; ++j_bb) {
            for (i_b = j_bb; i_b < b; ++i_b) {
                A[i_b + j_bb * b] -= L[i_b + k_b * b] * L[j_bb + k_b * b];
            }
        }
    }

    color = ",fillcolor=\"orange\"";
    store_tile(ctx, k + 1, k, k, L);
    color = "";

    if (opt_graph) {
        const int p = P;
        edge(p, TILI(k,k,k), TILI(k+1,k,k));
    }

    if (k < P) {
        kj_compute(ctx, k);
    }
}

static void k_compute(ctx_t *ctx)
{
    const int   p = P;
    int         k;
//	#pragma omp parallel for
    for (k = 0; k < p; ++k) {
        S1_compute(ctx, k);
    }
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

static void write_matrix(ctx_t *ctx, const char *outfile)
{
    FILE       *fp = strcmp(outfile, "-") ? xfopen(outfile, "w") : stdout;
    int         r, c, r_b, c_b, k;
    const int   p = P;
    const int   b = B;
    char        buf[40];

    if (opt_graph) {
        fprintf(opt_graph, "out [shape=Msquare,group=left,label=\"output\\nresult\"];\n");
    }

    for (r = 0; r < p; ++r) {
        for (r_b = 0; r_b < b; ++r_b) {
            k = 1;
            for (c = 0; c <= r; ++c) {
                tile_t tmp = get_tile(ctx, k, r, c);
                for (c_b = 0; (r != c) ? (c_b < b) : (c_b <= r_b); ++c_b) {
                    double a = tmp[r_b + c_b * b];
                    fnum(a, buf);
                    fputs(buf, fp);
                }
                if (opt_graph && !r_b) {
                    int src = TILI(k, r, c);
                    fprintf(opt_graph, "    %d -> out;\n", src);
                }
                ++k;
            }
            fprintf(fp, "\n");
        }
    }

    if (fp != stdout) {
        fclose(fp);
    }
}

static void decompose(const int n, const int b, double *A,
                      const char *outfile)
{
    const int   p = n / b;
    ctx_t      *ctx = make_context(n, b, p, A);
    struct timespec end;
    k_compute(ctx);

    if (clock_gettime(CLOCK_REALTIME, &end)) {
        pexit("clock_gettime");
	}
	printf("Time for size %d x %d : %lf sec\n", n, n, delta_time(begin, end));

    if (outfile) {
        write_matrix(ctx, outfile);
    }

    free_context(ctx);
}

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
            if (r != c) {
                A[c * n + r] = A[r * n + c];
            }
        }
    }
    fclose(fp);
}

static void dump_array(double *A, int n, FILE *fp)
{
    int         r, c;
    char        buf[40];

    for (r = 0; r < n; ++r) {
        for (c = 0; c < n; ++c) {
            fnum(A[r * n + c], buf);
            fprintf(fp, "%3s", buf);
        }
        fprintf(fp, "\n");
    }       
}

static void gen_matrix(int n, double *A)
{
    int r, c;
    for (r = 0; r < n; ++r) {
        for (c = 0; c < n; ++c) {
            A[r * n + c] = (r == c);
        }
    }
}

static void doit(int n, int b, const char *infile, const char *outfile,
                 const char *testfile)
{
    double *A = calloc(sizeof(double), n * n);

    if (infile) {
        read_matrix(n, A, infile);
        if (testfile) {
            FILE *fp = strcmp(testfile, "-") ? xfopen(testfile, "w") : stdout;
            dump_array(A, n, fp);
            if (fp != stdout) fclose(fp);
        }
    } else {
        gen_matrix(n, A);
    }

	if (clock_gettime(CLOCK_REALTIME, &begin)) {
		pexit("clock_gettime");
	}
    decompose(n, b, A, outfile);
    free(A);
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s N B [-i input] [-o output] [-w test] [-g graph] [-d] [-v]\n",
            prog);
    exit(1);
}

int main(int argc, char **argv)
{
    int         n = 0;
    int         b = 0;
    int         argi = 0;
    const char *infile = NULL;
    const char *outfile = NULL;
    const char *testfile = NULL;
    const char *graphfile = NULL;

    while (++argi < argc) {
        const char *args = argv[argi];
        if (args[0] == '-') {
            if (!strcmp(args, "-d")) opt_debug = true;
            else if (!strcmp(args, "-v")) opt_verbose = true;
            else if (!strcmp(args, "-i")) infile = argv[++argi];
            else if (!strcmp(args, "-o")) outfile = argv[++argi];
            else if (!strcmp(args, "-w")) testfile = argv[++argi];
            else if (!strcmp(args, "-g")) graphfile = argv[++argi];
            else usage(*argv);
        }
        else if (n == 0) {
            if ((n = atoi(args)) <= 0) {
                fprintf(stderr, "%s: Invalid matrix dimension \"%s\"\n",
                        *argv, args);
                usage(*argv);
            }
        }
        else if (b == 0) {
            if ((b = atoi(args)) <= 0) {
                fprintf(stderr, "%s: Invalid block size \"%s\"\n",
                        *argv, args);
                usage(*argv);
            }
        }
        else {
            usage(*argv);
        }
    }
    if (n == 0) {
        n = 16;
    }
    if (b == 0) {
        b = (int) (sqrt(n) + 0.1);
    }
    if (b < 1 || n % b) {
        fprintf(stderr, "Matrix size %d is not divisible by block size %d.\n",
                        n, b);
        exit(1);
    }

    if (graphfile) {
        opt_graph = xfopen(graphfile, "w");
        fprintf(opt_graph, "digraph serial_graph_dump {\n");
    }

    doit(n, b, infile, outfile, testfile);

    if (graphfile) {
        fprintf(opt_graph, "}\n");
        fclose(opt_graph);
    }

    return 0;
}

