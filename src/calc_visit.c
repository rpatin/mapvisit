#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>



/* ***********************************************************************
 *                                                                       *
 *                          Declaration of functions                     *
 *                                                                       *
 * ********************************************************************* */


/* from adehabitatHR. Author : C. Calenge
 xi and yi are the center of the circle of radius r
 x1 and y1 are the coordinates of the location inside the circle
 x2 and y2 are the coordinates of the location outside
 */
double interpLoc(double xi, double yi, double x1, double y1,
                 double x2, double y2, double r)
{
  double a, u, v, p, d, g;

  d = hypot(x2-x1, y2-y1);
  a = atan2(y2-y1, x2-x1);
  u = ((xi-x1)*cos(a))+((yi-y1)*sin(a));
  v = ((yi-y1)*cos(a))-((xi-x1)*sin(a));
  g = sqrt(R_pow(r, 2.0) - R_pow(v, 2.0));
  p = (u + g)/d;
  return(p);
}


/* Functions modified from adehabitatHR nvisits and HRresidtime */

SEXP calc_visit_c(SEXP xyt, SEXP xygrid, SEXP distr, SEXP maxt)
{
  /* declaring the variables */
  int n,i, j, *deds, sortie, *nvisi, ngrid;
  double *xr, *yr, *tr, *xgri, *ygri, dist, maxtr;
  double refti, refti2, p, *residtimer;
  SEXP x, y, t, dedsr, nvisit, residtime, returnlist, list_names, xgrid, ygrid;

  /* coercing the arguments */
  PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
  PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
  PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));
  PROTECT(xgrid = coerceVector(VECTOR_ELT(xygrid,0), REALSXP));
  PROTECT(ygrid = coerceVector(VECTOR_ELT(xygrid,1), REALSXP));
  n = length(x); /* number of relocations */
  ngrid = length(xgrid); /* number of points on the grid */
  PROTECT(dedsr = allocVector(INTSXP, n)); /* will be used to check whether
   the relocation j
   is within the distance distr
   of the grid point i */
  PROTECT(nvisit = allocVector(INTSXP, ngrid)); /* the output vector */
  PROTECT(residtime = allocVector(REALSXP, ngrid)); /* the output vector */


  char *names[2] = {"nvisit", "residtime"};

  PROTECT(list_names = allocVector(STRSXP,2));
  for(i = 0; i < 2; i++)
    SET_STRING_ELT(list_names,i,mkChar(names[i]));

  PROTECT(returnlist = allocVector(VECSXP, 2));

  /* for the ease of manipulation: gets the pointers */
  xr = REAL(x);
  yr = REAL(y);
  tr = REAL(t);
  xgri = REAL(xgrid);
  ygri = REAL(ygrid);
  deds = INTEGER(dedsr);
  nvisi = INTEGER(nvisit);
  residtimer = REAL(residtime);

    /* get the three constants passed as arguments */
  maxtr = REAL(maxt)[0];
  dist = REAL(distr)[0];
  /* Now, calculate the residence time for each relocation */
  for (i = 0; i < ngrid; i++) {
    /* No visit at the beginning */
    nvisi[i] = 0;

    /* checks which relocations are within the distance dist from grid point i  */
    for (j = 0; j < n; j++) {
      if (hypot(xgri[i]-xr[j], ygri[i]-yr[j])<=dist) {
        deds[j] = 1;
      } else {
        deds[j] = 0;
      }
    }

    /* forward time */
    sortie = 1; /* = 0 when the animal is still
     inside the circle; =1 otherwise */
    refti2 = tr[n]+maxtr;

    residtimer[i] = 0;

    /* if this is not the last relocation (otherwise, n forward visit = 0.0) */
    /* for all next relocations */
    for (j = 0; j<n; j++) {

      /* if the relocation is outside the circle */
      if (deds[j]==0) {

        /* if this is the first relocation outside */
        if (sortie==0) {
          sortie = 1;
          /* interpolating the time when the animal come out
           of the circle. Same procedure as above
           */
          p = interpLoc(xgri[i], ygri[i], xr[j-1], yr[j-1],
                        xr[j], yr[j], dist);
          refti2 = tr[j-1] + p*(tr[j]-tr[j-1]);
          residtimer[i] += refti2 - refti;
        }
      } else { /* if the relocation is inside the circle */
        if (sortie>0) {
          sortie = 0;
          p = interpLoc(xgri[i], ygri[i], xr[j], yr[j],
                        xr[j-1], yr[j-1], dist);
          refti = tr[j] + p*(tr[j-1]-tr[j]);

          if (fabs(refti2 - refti) > maxtr) {
            nvisi[i]++;
          }
        }
      }
    }
  }


  // attaching myint vector to list:
  SET_VECTOR_ELT(returnlist, 0, nvisit);
  // attaching mydouble vector to list:
  SET_VECTOR_ELT(returnlist, 1, residtime);
  // and attaching the vector names:
  setAttrib(returnlist, R_NamesSymbol, list_names);
  /* output */
  UNPROTECT(10);

  return(returnlist);

}
