/* HeightMapAlgorithm.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>

#include "HeightMapAlgorithm.h"

/* 
 * A note on memory usage:
 *
 * This implementation of the HeightMapAlgorithm is rather straight-forward but
 * consequently not designed for optimal memory efficiency. For example, it grows 
 * the output array A one row at a time, potentially requiring two copies of it 
 * to be present in memory at the same time. Instead, the output array A could be 
 * preallocated to its theoretical maximum size of 1/4 n^2 rows. 
 *   In addition, the output array A needs to be transformed from its (integer) 
 * canonical form to double-valued coordinates, again requiring two copies of the
 * array to be present in memory at the same time. An in-place transformation
 * is possible, but would require some trickery to interpret values of A either
 * as integers or as doubles.
 *   These memory issues are only a concern when the algorithm is run on datasets 
 * where the size of the output array may reach the computer's memory capacity. In
 * such cases, the above suggestions may be used to make the algorithm more memory
 * efficient.
 */

/* internal representation of an endpoint */
typedef struct { double x; int c, k, i; } EndPoint;

/* type definition of a rectangle with double coordinates */
typedef struct { double x1, x2, y1, y2; } DoubleRect;

/* type definition of a rectangle with integer coordinates */
typedef struct { int x1, x2, y1, y2; } IntRect;

/* type definition of rectangle boundaries (open/closed) */
typedef struct { int cx1, cx2, cy1, cy2; } RectBounds;

/* CompareEndpoints(A,B)
 
 pseudocode:
   1. c_A := (c_x_(k,i) = 1)            { boolean indicating A is a closed endpoint }           
   2. c_B := (c_x_(l,j) = 1)            { boolean indicating B is a closed endpoint }           
   3. r_A := (k = 2)                    { boolean indicating A is a right endpoint }            
   4. r_B := (l = 2)                    { boolean indicating B is a right endpoint }           
   5. if (x_(k,i) \ne x_(l,j)) then     { if the endpoints have different coordinates }         
   6.   return ( x_(k,i) < x_(l,j) )    { ...then let their coordinates determine their order } 
   ... 
*/
static int CompareEndpoints(const EndPoint *A, const EndPoint *B)
{
    int ca = (A->c == 1);
    int cb = (B->c == 1);
    int ra = (A->k == 2);
    int rb = (B->k == 2);

    if (A->x != B->x) 
        return (A->x < B->x);

    if (ra==rb && ca==cb)
        return (A->i < B->i);

    if (ra!=rb && ca!=cb)
        return ra;

    return (ra!=ca);
}

static int SortEndpoints(const void *elem1, const void *elem2)
{
    if (CompareEndpoints((EndPoint*) elem1, (EndPoint*) elem2))
        return -1;
    else
        return 1;
}

static void VerifyInput(SEXP R, SEXP B)
{
    int i,n;
    
    if (!isMatrix(R) || ncols(R)!=4)
        error("invalid first argument\n");

    if ((isMatrix(B) && !((nrows(B)==1 && (ncols(B)==2 || ncols(B)==4)) || (nrows(B)==nrows(R) && ncols(B)==4)))
     || (isVector(B) && length(B)!=2 && length(B)!=4))
        error("invalid second argument\n");

    n = nrows(R);
    for (i=0; i<n; i++)
    {
        if (REAL(R)[i]>REAL(R)[i+n])
            error("x1 not less than or equal to x2 in rectangle %d\n", i+1);
        if (REAL(R)[i+2*n]>REAL(R)[i+3*n])
            error("y1 not less than or equal to y2 in rectangle %d\n", i+1);
    }

    for (i=0; i<length(B); i++)
        if (INTEGER(B)[i]!=0 && INTEGER(B)[i]!=1)
            error("second argument may only contain 0's and 1's\n");
}

static IntRect* ObsRectToIntRect(SEXP RR, SEXP BB)
{
    /* retrieve n from R-object RR */
    int n = nrows(RR);
    int i;

    EndPoint *XEndPoints = Calloc(2*n, EndPoint);
    EndPoint *YEndPoints = Calloc(2*n, EndPoint);
	IntRect  *R          = Calloc(n+1, IntRect);

    double   *RRdata = REAL(RR);
    int      *BBdata = INTEGER(BB);
    int       BBexplicit;
    int       BBvalue[4];

    BBexplicit = (isMatrix(BB) && nrows(BB)==nrows(RR));
    if (!BBexplicit)
    {
        BBvalue[0] = BBdata[0];
        BBvalue[1] = BBdata[1];
        BBvalue[2] = BBdata[length(BB)-2+0];
        BBvalue[3] = BBdata[length(BB)-2+1];
    }

    /* separate the observation rectangles R into X and Y endpoints */
    for (i=0; i<n; i++)
    {
        XEndPoints[i*2].x = RRdata[i+0*n];
        XEndPoints[i*2].c = BBexplicit ? BBdata[i+0*n] : BBvalue[0];
        XEndPoints[i*2].k = 1;
        XEndPoints[i*2].i = i+1;

        YEndPoints[i*2].x = RRdata[i+2*n];
        YEndPoints[i*2].c = BBexplicit ? BBdata[i+2*n] : BBvalue[2];
        YEndPoints[i*2].k = 1;
        YEndPoints[i*2].i = i+1;

        XEndPoints[i*2+1].x = RRdata[i+1*n];
        XEndPoints[i*2+1].c = BBexplicit ? BBdata[i+1*n] : BBvalue[1];
        XEndPoints[i*2+1].k = 2;
        XEndPoints[i*2+1].i = i+1;

        YEndPoints[i*2+1].x = RRdata[i+3*n];
        YEndPoints[i*2+1].c = BBexplicit ? BBdata[i+3*n] : BBvalue[3];
        YEndPoints[i*2+1].k = 2;
        YEndPoints[i*2+1].i = i+1;
    }

    /* sort endpoint arrays using quicksort with CompareEndpoints (indirectly) */
    qsort(XEndPoints, n*2, sizeof(EndPoint), SortEndpoints);
    qsort(YEndPoints, n*2, sizeof(EndPoint), SortEndpoints);

    /* transform the sorted arrays into rectangles in canonical form */
    for (i=0; i<2*n; i++)
    {
        if (XEndPoints[i].k == 1)
            R[XEndPoints[i].i].x1 = i+1;
        else
            R[XEndPoints[i].i].x2 = i+1;

        if (YEndPoints[i].k == 1)
            R[YEndPoints[i].i].y1 = i+1;
        else
            R[YEndPoints[i].i].y2 = i+1;
    }

    Free(XEndPoints);
    Free(YEndPoints);

    return R;
}

static int SortDoubles(const void *elem1, const void *elem2)
{
    double d1 = *((double*) elem1);
    double d2 = *((double*) elem2);

    if (d1 == d2)
        return 0;
    if (d1 < d2)
        return -1;
    return 1;
}

static SEXP IntRectToDoubleRect(int m, IntRect *A, SEXP R)
{
    /* retrieve n from R-object R */
    int n = nrows(R);

    double *xcoords = Calloc(2*n, double);
    double *ycoords = Calloc(2*n, double);
    double *Rdata   = REAL(R);
    SEXP    Aout;
    double *Aoutdata;

    int i;

    /* allocate m x 4 matrix of reals in R workspace */
    PROTECT(Aout = allocMatrix(REALSXP, m, 4));
    Aoutdata = REAL(Aout);

    for (i=0; i<n; i++)
    {
        xcoords[2*i]   = Rdata[i+0*n];
        xcoords[2*i+1] = Rdata[i+1*n];

        ycoords[2*i]   = Rdata[i+2*n];
        ycoords[2*i+1] = Rdata[i+3*n];
    }

    qsort(xcoords, n*2, sizeof(double), SortDoubles);
    qsort(ycoords, n*2, sizeof(double), SortDoubles);

    for (i=1; i<=m; i++)
    {
        Aoutdata[i-1+0*m] = xcoords[A[i].x1-1];
        Aoutdata[i-1+1*m] = xcoords[A[i].x2-1];
        Aoutdata[i-1+2*m] = ycoords[A[i].y1-1];
        Aoutdata[i-1+3*m] = ycoords[A[i].y2-1];
    }

    Free(xcoords);
    Free(ycoords);

    return Aout;
}

static IntRect *AddToA(IntRect *A, int m, int x, int j, int b, int k)
{
    /* increment size of A array by 1 */
    IntRect *Anew = Realloc(A, m+1, IntRect);

    /* add new element to array */
    Anew[m].x1 = x;
    Anew[m].x2 = j;
    Anew[m].y1 = b;
    Anew[m].y2 = k;

    return Anew;
}

#define X1(i) (R[i].x1)
#define X2(i) (R[i].x2)
#define Y1(i) (R[i].y1)
#define Y2(i) (R[i].y2)

SEXP HeightMapAlgorithm(SEXP RR, SEXP BB)
{
    int b, i, j, k, m, n;
    int *h, *e, *r, *lb;

    IntRect *A = NULL;
	IntRect *R = NULL;
    SEXP AA;

    /* 0. initialization */
    VerifyInput(RR,BB);
    n  = nrows(RR);
    h  = Calloc( (2*n+1), int);
    e  = Calloc( (2*n+1), int);
    r  = Calloc( (2*n+1), int);
    lb = Calloc( (2*n+1), int);

    /* 1. transform to canonical rectangles */
    R = ObsRectToIntRect(RR, BB);

    /* 2. sort x_(1,i) and x_(2,i) (i=1...n) in ascending order, 
          and store their indices i in r_j (j=1...2n) */
    for (i=1; i<=n; i++)
    {
        r[X1(i)]  = i;
        lb[X1(i)] = 1;

        r[X2(i)]  = i;
        lb[X2(i)] = 0;
    }

    /* 3. m := 0 */
    m = 0;

    /* 4. and 5. set h_1..h_2n and e_1..e_2n to zero */
    for (i=1; i<=2*n; i++)
    {
        h[i] = 0; e[i] = 0;
    }

    /* 6. for j = 1 to 2n do */
    for (j=1; j<=2*n; j++)
    {
        /* 7. if r_j is a left boundary */
        if (lb[j])
        {
            /* 8. for k = y_(1,r_j) + 1 to y_(2,r_y) do */
            for (k=Y1(r[j]) + 1; k<=Y2(r[j]); k++)
            {
                /* 9. h_k := h_k + 1; e_k := r_j */
                h[k] += 1;
                e[k] = r[j];
            }
        }
        /* 10. else */
        else 
        {
            /* 11. b := y_(1,r_j) */
            b = Y1(r[j]);

            /* 12. for k = y_(1,r_j) + 1 to y_(2,r_j) - 1 do */
            for (k=Y1(r[j]) + 1; k<=Y2(r[j]) - 1; k++)
            {
                /* 13. if (h_(k+1) < h_k and b > 0) then */
                if ((h[k+1]<h[k]) && (b>0))
                {
                    /* 14. if (e_k > 0) then */
                    if (e[k]>0)
                    {
                        /* 15. m:=m+1; A_m := (x_(1,e_k), j, b, k) */
                        m += 1;
                        A = AddToA(A, m, X1(e[k]), j, b, k);
                    }
                    /* 16. b := 0 */
                    b = 0;
                }
                /* 17. if (h_(k+1) > h_k) then */
                if (h[k+1]>h[k])
                {
                    /* 18. b := k */
                    b = k;
                }
            }

            /* 19. if (b > 0 and e_k > 0) then */
            if ((b>0) && (e[k]>0))
            {
                /* 20. m := m + 1; A_m := (x_(1,e_k), j, b, k) */
                m += 1;
                A = AddToA(A, m, X1(e[k]), j, b, k);
            }

            /* 21. k := y_(1,r_j) */
            k = Y1(r[j]);

            /* 22. while (k >= 1 and h_k <= h_(k+1) ) do  */
            while ((k>=1) && (h[k]<=h[k+1]))
            {
                /* 23. e_k = 0; k := k - 1 */
                e[k] = 0;
                k -= 1;
            }

            /* 24. k := y_(2,r_j) + 1 */
            k = Y2(r[j]) + 1;

            /* 25. while (k <= 2n and h_k <= h_(k-1)) do */
            while ((k<=2*n) && (h[k]<=h[k-1]))
            {
                /* 26. e_k := 0; k := k + 1 */
                e[k] = 0;
                k += 1;
            }

            /* 27. for k = y_(1,r_j) + 1 to y_(2,r_j) do */
            for (k=Y1(r[j]) + 1; k<=Y2(r[j]); k++)
            {
                /* 28. h_k := h_k - 1; e_k := 0 */
                h[k] -= 1;
                e[k] = 0;
            }
        }
    }

    /* 29. transform the canonical maximal intersections A_1..A_m 
           back to the original coordinates */
    AA = IntRectToDoubleRect(m, A, RR);

    /* release memory of all intermediate data structures */
    Free(h);
    Free(e);
    Free(r);
    Free(lb);
    Free(A);
    Free(R);

    /* 30. return A_1..A_m */
    UNPROTECT(1);
    return AA;
}

/*
 * vim: ts=4 sw=4 et fileformat=unix
 */
