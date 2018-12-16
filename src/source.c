/* History Sep 02 2015  Round cutpoints down
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


#define DOUBLE_MISS 999999999.9
#define DOUBLE_MISS_TEST 999999999.0
#define CHECK_MEM(obj) if (obj == NULL) {Rprintf("ERROR: allocating memory \n"); error("1");}
#define MIN_LOG_ARG 1e-300
#define LOG0 -690.7755
#define SMALL 1e-12
#define LOG10 2.3025850929940459011
#define MYINF 1.0e100

void ref_risk1(double *ref_risks, double *betavec, double *zmat, double *refmat, int *n_beta,\
   int *nr_z, int *nc_z, int *nr_ref, int *nc_ref, double *probs, int *n_probs, int *pdebug,\
   double *weights, double *ref_full_LP, double *retvec, int *retflag, double *retlps);

void ref_risk2(double *betavec, double *zmat, double *refmat, int *n_beta, int *nr_z, int *nc_z,\
   int *nr_ref, int *nc_ref, double *probs, int *n_probs, int *pdebug, int *age_new, int *age_int,\
   int *n_lambda, double *popSubFromLP, double *lambda, double *compRates, double *weights,\
   double *ref_full_LP, double *retvec, int *retflag, double *retlps);

static const R_CMethodDef callMethods[] = {
  {"ref_risk1", (DL_FUNC)&ref_risk1, 17},
  {"ref_risk2", (DL_FUNC)&ref_risk2, 22},
  {NULL, NULL, 0}
};


/****************** Memory allocation ************************/
/* Function to allocate memory for a double vector */
double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  /*ret = (double *) malloc(n*sizeof(double));*/
  ret = (double *) R_Calloc(n, double);
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }
  
  return(ret);

} /* END: dVec_alloc */

/* Function to allocate memory for a integer vector */
int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  /*ret = (int *) malloc(n*sizeof(int));*/
  ret = (int *) R_Calloc(n, int);
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }
  
  return(ret);

} /* END: iVec_alloc */

/* Function to allocate an integer matrix */
double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  /*mat = (double **) malloc(nrow*sizeof(double *));*/
  mat = (double **) R_Calloc(nrow, double *);
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to free a matrix */
void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) R_Free(x[i]);
  }
  R_Free(x);

} /* END: matrix_free */

/***************************************************************************/


/* Function to get non missing indices in a double vector */
static int get_non_miss(vec, n, ovec)
double *vec;
int n, *ovec;
{
  double *pd;
  int i, *pint, ret;

  ret  = 0;
  pint = ovec;
  for (i=0, pd=vec; i<n; i++, pd++) {
    if (!(*pd > DOUBLE_MISS_TEST)) {
      *pint++ = i;
      ret++;
    }
  }

  return(ret);

} /* END: get_non_miss */

/* Function to compute Xy for matrix X and vector y based on index vector */
static void X_y(X, nr, nc, y, index, ni, ret)
double **X, *y, *ret;
int nr, nc, *index, ni;
{
  int i, j, k, *pint;
  double sum, *pret, *vec, **row;

  for (j=0, pret=ret, row=X; j<nr; j++, pret++, row++) {
    sum = 0.0;
    vec = *row;
    for (i=0, pint=index; i<ni; i++, pint++) {
      k = *pint;
      sum += vec[k]*y[k];
    }
    *pret = sum; 
  }

} /* END: X_y */

/* Function to compute X[temp]y-z[temp] for matrix X and vector y, z based on index vector */
static void X_y2(X, nc, y, z, index, ni, ret)
double **X, *y, *z, *ret;
int nc, *index, ni;
{
  int i, j, *pint, row, ok;
  double sum, *p, *pret, *vec, *py;

  for (i=0, pint=index, pret=ret; i<ni; i++, pint++, pret++) {
    row  = *pint;
    vec  = X[row];
    sum  = 0.0;
    ok   = 1;
    for (j=0, p=vec, py=y; j<nc; j++, p++, py++) {
      if (!(*p > DOUBLE_MISS_TEST)) {
        sum += *p * *py;
      } else {
        ok = 0;
        break;
      }
    }
    if (ok) { 
      *pret = sum - z[row];
    } else {
      *pret = DOUBLE_MISS;
    } 
  }

} /* END: X_y2 */

/* Function to compute x*y for vectors x and y based on index vector */
static double vec_mult(x, y, index, ni)
double *x, *y;
int *index, ni;
{
  int i, k, *pint;
  double sum;

  sum = 0.0;
  for (i=0, pint=index; i<ni; i++, pint++) {
    k = *pint;
    sum += x[k]*y[k];
  }
 
  return(sum);

} /* END: vec_mult */

/* Function to get the missing percentile */
static int vec_lteq_c(vec, n, c)
double *vec, c; /* vec must be sorted in ascending order */
int n;
{
  int j, ret;
  double *p;

  ret = 0;
  for (j=0, p=vec; j<n; j++, p++) {
    if (*p <= c) {
      ret++;
    } else {
      break;
    }
  }
  
  return(ret);

} /* vec_lteq_c */

/* Function to compute mean of vector based on vec between 2 values */
static double mean_2val(vec, subvec, nsub, weights)
double *vec, *weights;
int *subvec, nsub;
{
  int i, *pii;
  double sum, *pw;

  sum = 0.0;
  for (i=0, pw=weights, pii=subvec; i<nsub; i++, pw++, pii++) {
    sum += vec[*pii]* *pw;
  }

  return(sum);

} /* mean_2val */

/* Function to create matrix from a vector */
static void mat_from_vec(vec, nr, nc, mat)
double *vec, **mat; /* vec must be the columns stacked */
int nr, nc;
{
  int i, j;
  double *p;

  p = vec;
  for (i=0; i<nc; i++) {
    for (j=0; j<nr; j++) mat[j][i] = *p++;
  }

} /* END: mat_from_vec */

/* Function to get quantile vectors */
static int quant_vecs(probs, nprobs, n, ivec, lo, hi, h)
double *probs, *h;
int n, *ivec, *lo, *hi, nprobs;
{
  int i, *pi1, *plo, *phi, low, n_i;
  double *p, index, *ph, nm1;

  n_i = 0;
  nm1 = (double)n - 1.0;
  pi1 = ivec;
  ph  = h;
  for (i=0, p=probs, phi=hi, plo=lo; i<nprobs; i++, p++, phi++, plo++) {
    index = nm1* *p; /* Subtract 1 since arrays start at 0 in C */
    low   = (int) floor(index);
    *plo  = low;
    *phi  = (int) ceil(index);
    if (index > low) {
      n_i++;
      *pi1++ = i; 
      *ph++  = index - low;
    }
  }

  return(n_i);

} /* END: quant_vecs */

/* Function to compute quantiles */
static void quantile(x, n, probs, nprobs, lo, hi, h, ivec, ni, qs)
double *x, *probs, *h, *qs;
int n, nprobs, *ivec, ni, *lo, *hi;
{
  int last, flag, j, *q_ivec, *pi1, i, *q_lo, *q_hi;
  double *pqs, *q_h, *p, temp;

  last = n-1;
  flag = 0;

  R_qsort(x, 1, n);

  if (x[last] > DOUBLE_MISS_TEST) {

    /* Get the number of non-missing */
    last = 0;
    for (i=0, p=x; i<n; i++, p++) {
      if (!(*p > DOUBLE_MISS_TEST)) last++;
    }
    flag   = 1;
    q_lo   = iVec_alloc(nprobs, 0, 0);
    q_hi   = iVec_alloc(nprobs, 0, 0);
    q_h    = dVec_alloc(nprobs, 0, 0);
    q_ivec = iVec_alloc(nprobs, 0, 0);
    ni     = quant_vecs(probs, nprobs, last, q_ivec, q_lo, q_hi, q_h);
  } else {
    q_lo   = lo;
    q_hi   = hi;
    q_h    = h;
    q_ivec = ivec;
  }

  /* qs <- x[lo] */
  for (j=0, pqs=qs, pi1=q_lo; j<nprobs; j++, pqs++, pi1++) {
    *pqs = x[*pi1];
  }
  
  /* qs[i] <- (1 - h) * qs[i] + h * x[hi[i]] */
  for (j=0, pi1=q_ivec, p=q_h; j<ni; j++, pi1++, p++) {
   i     = *pi1;
   temp  = *p;
   qs[i] = (1.0 - temp)*qs[i] + temp*x[q_hi[i]]; 
  }

  if (flag) {
    R_Free(q_lo);
    R_Free(q_hi);
    R_Free(q_h);
    R_Free(q_ivec);
  }

} /* END: quantile */

/* Function to compute the mean and variance of a double vector */
static double var_dvec(vec, n)
double *vec;
int n;
{
  double sum=0.0, sum2=0.0, ret, *p, temp, variance;
  int i;

  for (i=0, p=vec; i<n; i++, p++) {
    temp = *p;
    sum += temp;
    sum2 += temp*temp;
  }
  ret = sum/n;
  variance = (sum2 - n*ret*ret)/(n - 1.0);

  return(variance);

} /* END: var_dvec */

/* Function to get the unique cutpoints rounded */
static int get_cutpoints(quantiles, nq, ref_LP, nref)
double *quantiles, *ref_LP;
int nq, nref;
{
  int i, nuniq;
  double se, prec, temp, digits, *p, temp2;

  /* Get the se of ref_LP */
  se   = sqrt(var_dvec(ref_LP, nref));
  if (se < SMALL) se = SMALL;
  prec = 1.0/(0.001*se);

  /* Get the number of digits of precision */
  digits = 0.0;
  temp   = 10.0;
  while (1) {
    if (prec >= temp) {
      digits += 1.0;
      temp *= 10.0;
    } else {
      break;
    }
  }

  /* Round quantiles. 02Sep2015 round down  */
  temp  = pow(10.0, (double) digits);
  temp2 = 1.0/temp;
  for (i=0, p=quantiles; i<nq; i++, p++) *p = temp2*floor(*p * temp);
  
  /* Get unique quantiles for the cutpoints */
  p     = quantiles;
  nuniq = 1;
  for (i=1; i<nq; i++) {
    if (!(fabs(p[i] - p[i-1]) < SMALL)) {
      p[nuniq] = p[i];
      nuniq++; 
    } 
  }

  /* Add infinity */
  p[nuniq] = MYINF;
  nuniq++; 

  return(nuniq);

} /* END: get_cutpoints */

/* Function to get the missing percentile */
static int get_miss_perc(cutpoints, ncuts, miss_LP)
double *cutpoints, miss_LP;
int ncuts;
{
  int miss_perc, ncm1;  

  ncm1 = ncuts - 1;

  /* Get the missing percentile */
  miss_perc = vec_lteq_c(cutpoints, ncuts, miss_LP) - 1;
  if (miss_perc < 0) {
    miss_perc = 0;
  } else if (miss_perc == ncm1) {
    miss_perc--;
  }
  
  return(miss_perc);

} /* END: get_miss_perc */

/* Function to determine if any vector elements between 2 values */
static int vec_ab(x, n, a, b)
double *x, a, b;
int n;
{
  int i;
  double *p, temp;

  for (i=0, p=x; i<n; i++, p++) {
    temp = *p;
    if ((temp >= a) && (temp < b)) return(1);
  }
  
  return(0);

} /* END: vec_ab */

/* Function to get new endpoints */
static void get_newEndpoints(cutpoints, ncuts, miss_LP, miss_perc, ref_LP, nref, rcut1, rcut2)
double *cutpoints, miss_LP, *ref_LP, *rcut1, *rcut2;
int miss_perc, ncuts, nref;
{
  int i, ncm1, left, right, m;
  double cut1, cut2;

  ncm1  = ncuts - 1;
  cut1  = cutpoints[miss_perc];
  cut2  = cutpoints[miss_perc + 1];
  left  = miss_perc - 1;
  right = miss_perc + 2;
  for (i=0; i<ncuts; i++) {
    if (left < 0) left = 0;
    if (right > ncm1) right = ncm1;
    cut1 = cutpoints[left];
    cut2 = cutpoints[right];
    m    = vec_ab(ref_LP, nref, cut1, cut2);
    if (m) break;
  }

  *rcut1 = cut1;
  *rcut2 = cut2;

} /* END: get_newEndpoints */

/* Function to get subjects between 2 cutpoints */
static int get_subs_2cut(x, n, a, b, ret) 
double *x, a, b;
int *ret, n;
{
  int i, *pret, j;
  double *p, val;

  j    = 0;
  pret = ret;
  for (i=0, p=x; i<n; i++, p++) {
    val = *p;
    if ((val >= a) && (val < b)) {
      *pret++ = i;
      j++;
    }
  }

  return(j);

} /* END: get_subs_2cut */

/* Function to get vector of weights that sum to 1 */
static void get_wgts(weights, subvec, nsub, retvec)
double *weights, *retvec;
int *subvec, nsub;
{
  double sumWgt, *p1, tempd;
  int j, *pii;

  sumWgt = 0.0;
  for (j=0, p1=retvec, pii=subvec; j<nsub; j++, p1++, pii++) {
    tempd  = weights[*pii];
    *p1    = tempd;
    sumWgt += tempd;
  }
  for (j=0, p1=retvec; j<nsub; j++, p1++) *p1 = *p1/sumWgt;

} /* END: get_wgts */

/* Only pass in missing rows for Z_new */
void ref_risk1(ref_risks, betavec, zmat, refmat, n_beta, nr_z, nc_z, nr_ref, nc_ref,
   probs, n_probs, pdebug, weights, ref_full_LP, retvec, retflag, retlps)
double *ref_risks, *betavec, *zmat, *refmat, *retvec, *probs, *weights, *ref_full_LP, *retlps;
int *retflag, *n_beta, *nr_z, *nc_z, *nr_ref, *nc_ref, *n_probs, *pdebug;
{
  double **REF, **ZMAT, *Zvec, miss_LP, *cutpoints, *ref_LP, *ref_LP2, *q_h, *p1, *p2;
  int ii, zmat_nrow, zmat_ncol, ref_nrow, ref_ncol, nprobs, ncuts, DEBUG;
  int *present, nonmiss_n, j, *q_ivec, n_ivec, *q_hi, *q_lo, ncutpoints, miss_perc, *subvec, nsub;
  double cut1, cut2, *wgts;

  *retflag = 1;

  DEBUG     = *pdebug;
  zmat_nrow = *nr_z;
  zmat_ncol = *nc_z;
  ref_nrow  = *nr_ref;
  ref_ncol  = *nc_ref;
  nprobs    = *n_probs;
  ncuts     = nprobs;
  for (j=0; j<zmat_nrow; j++) retvec[j] = DOUBLE_MISS;
  
  /* Allocate memory */
  ZMAT      = dMat_alloc(zmat_nrow, zmat_ncol, 0, 0.0);
  REF       = dMat_alloc(ref_nrow, ref_ncol, 0, 0.0);
  present   = iVec_alloc(zmat_ncol, 0, 0);
  ref_LP    = dVec_alloc(ref_nrow, 0, 0);
  ref_LP2   = dVec_alloc(ref_nrow, 0, 0);
  cutpoints = dVec_alloc(ncuts+1, 0, 0); /* The cutpoint Inf is added later */
  q_lo      = iVec_alloc(nprobs, 0, 0);
  q_hi      = iVec_alloc(nprobs, 0, 0);
  q_h       = dVec_alloc(nprobs, 0, 0);
  q_ivec    = iVec_alloc(nprobs, 0, 0);
  subvec    = iVec_alloc(ref_nrow, 0, 0);
  wgts      = dVec_alloc(ref_nrow, 0, 0);

  /* Fill the matrices */
  mat_from_vec(zmat, zmat_nrow, zmat_ncol, ZMAT);
  mat_from_vec(refmat, ref_nrow, ref_ncol, REF);

  /* Quantile vectors when there are no missing values */
  n_ivec = quant_vecs(probs, nprobs, ref_nrow, q_ivec, q_lo, q_hi, q_h);

  /* Loop over each missing */
  for (ii=0; ii<zmat_nrow; ii++) {
    if (DEBUG) Rprintf("iter = %d\n", ii);
    Zvec = ZMAT[ii]; 

    /* Get non-missing variables in zmat */
    nonmiss_n = get_non_miss(Zvec, zmat_ncol, present);

    /* Compute LP for the reference set */
    X_y(REF, ref_nrow, ref_ncol, betavec, present, nonmiss_n, ref_LP);
    
    /* Copy ref_LP into ref_LP2, because it is sorted in quantile function */
    for (j=0, p1=ref_LP, p2=ref_LP2; j<ref_nrow; j++, p1++, p2++) *p2 = *p1;

    /* Compute missing LP */
    miss_LP = vec_mult(Zvec, betavec, present, nonmiss_n);

    /* Compute the cutpoints */
    quantile(ref_LP2, ref_nrow, probs, nprobs, q_lo, q_hi, q_h, q_ivec, n_ivec, cutpoints);

    /* Get the rounded, unique cutpoints */
    ncutpoints = get_cutpoints(cutpoints, ncuts, ref_LP, ref_nrow);

    /* Get the missing percentile and cutpoint interval */
    miss_perc = get_miss_perc(cutpoints, ncutpoints, miss_LP);
    cut1      = cutpoints[miss_perc];
    cut2      = cutpoints[miss_perc + 1];

    /* Get the subjects between the 2 cutpoints */
    nsub      = get_subs_2cut(ref_LP, ref_nrow, cut1, cut2, subvec); 
    if (!nsub) {
      /* Get new endpoints */
      get_newEndpoints(cutpoints, ncutpoints, miss_LP, miss_perc, ref_LP, ref_nrow, &cut1, &cut2);
      nsub = get_subs_2cut(ref_LP, ref_nrow, cut1, cut2, subvec); 
    }

    if (DEBUG) Rprintf("ncutpoints=%d, cut1=%20.15f, cut2=%20.15f\n", ncutpoints, cut1, cut2);
    /*if ((cut1 > DOUBLE_MISS_TEST) || (cut2 > DOUBLE_MISS_TEST)) continue;*/

    /* Get the weights. wgts will have length nsub */
    get_wgts(weights, subvec, nsub, wgts);

    retvec[ii] = mean_2val(ref_risks, subvec, nsub, wgts);
    retlps[ii] = mean_2val(ref_full_LP, subvec, nsub, wgts);
  }

  matrix_free((void **) ZMAT, zmat_nrow);
  matrix_free((void **) REF, ref_nrow);
  R_Free(present);
  R_Free(ref_LP);
  R_Free(ref_LP2);
  R_Free(cutpoints);
  R_Free(q_lo);
  R_Free(q_hi);
  R_Free(q_h);
  R_Free(q_ivec);
  R_Free(subvec);
  R_Free(wgts);

  *retflag = 0;

  return;

} /* END: ref_risk1 */


static void get_int(nsub, age_new, t, aexp_rr, lambda, comp, zbeta, holder)
int nsub, age_new, t;
double aexp_rr, *lambda, *comp, *zbeta, *holder;
{
  int u, j;
  double hold1, *ph, *p, c1, lambdau, temp, compu;

  /* Initialize */
  c1 = 1.0/aexp_rr;
  for (j=0, ph=holder; j<nsub; j++, ph++) *ph = 0.0;

  for (u=age_new; u<t; u++) {
    lambdau = lambda[u];
    temp    = lambdau*c1;
    compu   = comp[u];
    for (j=0, ph=holder, p=zbeta; j<nsub; j++, ph++, p++) {
      hold1 = temp* *p;
      *ph   = *ph + hold1 + compu;  

    }


  }

} /* END: get_int */

static double finalVal(nsub, zbeta, expzbeta, age_new, age_int, lambda, loglambda, compRates, holder, tempv, weights)
double *zbeta, *expzbeta, *lambda, *loglambda, *compRates, *holder, *tempv, *weights; 
int nsub, age_new, age_int;
{
  int i, age2, t, j;
  double ret, this_year, log_lam_t, *p, *ph, *pzb, val, *pw;
  
  age2 = age_new + age_int;

  for (i=0, p=tempv; i<nsub; i++, p++) {
    *p = 0.0;
  }

  for (t=age_new; t<age2; t++) {
    get_int(nsub, age_new, t, 1.0, lambda, compRates, expzbeta, holder);
    log_lam_t = loglambda[t];
    for (i=0, p=tempv, ph=holder, pzb=zbeta; i<nsub; i++, p++, ph++, pzb++) {
      this_year = exp(log_lam_t + *pzb - *ph);
      *p = *p + this_year;
    }
  }

  /* Compute the weighted mean */
  ret = 0.0;
  j   = 0;
  for (i=0, p=tempv, pw=weights; i<nsub; i++, p++, pw++) {
    val  = *p;
    if (!(val > DOUBLE_MISS_TEST)) {
      ret += val* *pw;
      j++;
    } 
  }
  if (!j) ret = DOUBLE_MISS;

  return(ret);

} /* END: finalVal */

/* Only pass in missing rows for Z_new */
void ref_risk2(betavec, zmat, refmat, n_beta, nr_z, nc_z, nr_ref, nc_ref,
   probs, n_probs, pdebug, age_new, age_int, n_lambda, popSubFromLP, lambda, compRates,
   weights, ref_full_LP, retvec, retflag, retlps)
double *betavec, *zmat, *refmat, *retvec, *probs, *ref_full_LP, *retlps;
double *popSubFromLP, *lambda, *compRates, *weights;
int *retflag, *n_beta, *nr_z, *nc_z, *nr_ref, *nc_ref, *n_probs, *pdebug;
int *age_new, *age_int, *n_lambda;
{
  double **REF, **ZMAT, *Zvec, miss_LP, *cutpoints, *ref_LP, *ref_LP2, *q_h, *p1, *p2;
  int ii, zmat_nrow, zmat_ncol, ref_nrow, ref_ncol, nprobs, ncuts, DEBUG;
  int *present, nonmiss_n, j, *q_ivec, n_ivec, *q_hi, *q_lo, nlambda;
  double cut1, cut2, *zbeta, tempd, *holder, *loglambda, *expzbeta, *wgts;
  int *subvec, nsub, ncutpoints, miss_perc;
  
  *retflag = 1;

  DEBUG     = *pdebug;
  zmat_nrow = *nr_z;
  zmat_ncol = *nc_z;
  ref_nrow  = *nr_ref;
  ref_ncol  = *nc_ref;
  nprobs    = *n_probs;
  ncuts     = nprobs;
  nlambda   = *n_lambda;
  for (j=0; j<zmat_nrow; j++) retvec[j] = DOUBLE_MISS;
  
  /* Allocate memory */
  ZMAT      = dMat_alloc(zmat_nrow, zmat_ncol, 0, 0.0);
  REF       = dMat_alloc(ref_nrow, ref_ncol, 0, 0.0);
  present   = iVec_alloc(zmat_ncol, 0, 0);
  ref_LP    = dVec_alloc(ref_nrow, 0, 0);
  ref_LP2   = dVec_alloc(ref_nrow, 0, 0);
  cutpoints = dVec_alloc(ncuts+1, 0, 0); /* The cutpoint Inf is added later */
  q_lo      = iVec_alloc(nprobs, 0, 0);
  q_hi      = iVec_alloc(nprobs, 0, 0);
  q_h       = dVec_alloc(nprobs, 0, 0);
  q_ivec    = iVec_alloc(nprobs, 0, 0);
  subvec    = iVec_alloc(ref_nrow, 0, 0);
  zbeta     = dVec_alloc(ref_nrow, 0, 0);
  expzbeta  = dVec_alloc(ref_nrow, 0, 0);
  loglambda = dVec_alloc(nlambda, 0, 0);
  holder    = dVec_alloc(ref_nrow, 0, 0);
  wgts      = dVec_alloc(ref_nrow, 0, 0);

  /* Fill the matrices */
  mat_from_vec(zmat, zmat_nrow, zmat_ncol, ZMAT);
  mat_from_vec(refmat, ref_nrow, ref_ncol, REF);

  /* Quantile vectors when there are no missing values */
  n_ivec = quant_vecs(probs, nprobs, ref_nrow, q_ivec, q_lo, q_hi, q_h);

  /* Compute log(lambda) */
  for (j=0, p1=lambda, p2=loglambda; j<nlambda; j++, p1++, p2++) {
    tempd = *p1;
    if (tempd > MIN_LOG_ARG) {
      *p2 = log(tempd);
    } else {
      *p2 = LOG0;
    }
  }

  /* Loop over each missing */
  for (ii=0; ii<zmat_nrow; ii++) {
    if (DEBUG) Rprintf("iter = %d\n", ii);

    Zvec = ZMAT[ii]; 

    /* Get non-missing variables in zmat */
    nonmiss_n = get_non_miss(Zvec, zmat_ncol, present);

    /* Compute LP for the reference set */
    X_y(REF, ref_nrow, ref_ncol, betavec, present, nonmiss_n, ref_LP);
    
    /* Copy ref_LP into ref_LP2, because it is sorted in quantile function */
    for (j=0, p1=ref_LP, p2=ref_LP2; j<ref_nrow; j++, p1++, p2++) *p2 = *p1;

    /* Compute missing LP */
    miss_LP = vec_mult(Zvec, betavec, present, nonmiss_n);

    /* Compute the quantiles */
    quantile(ref_LP2, ref_nrow, probs, nprobs, q_lo, q_hi, q_h, q_ivec, n_ivec, cutpoints);

    /* Get the rounded, unique cutpoints */
    ncutpoints = get_cutpoints(cutpoints, ncuts, ref_LP, ref_nrow);

    /* Get the missing percentile and cutpoint interval */
    miss_perc = get_miss_perc(cutpoints, ncutpoints, miss_LP);
    cut1      = cutpoints[miss_perc];
    cut2      = cutpoints[miss_perc + 1];

    /* Get the subjects between the 2 cutpoints */
    nsub      = get_subs_2cut(ref_LP, ref_nrow, cut1, cut2, subvec); 
    if (!nsub) {
      /* Get new endpoints */
      get_newEndpoints(cutpoints, ncutpoints, miss_LP, miss_perc, ref_LP, ref_nrow, &cut1, &cut2);
      nsub = get_subs_2cut(ref_LP, ref_nrow, cut1, cut2, subvec); 
    }

    if (DEBUG) Rprintf("nsub = %d, ncutpoints=%d, cut1=%20.15f, cut2=%20.15f\n", nsub, ncutpoints, cut1, cut2);
    /*if ((cut1 > DOUBLE_MISS_TEST) || (cut2 > DOUBLE_MISS_TEST)) continue;*/

    /* Compute Zbeta-sub */
    X_y2(REF, ref_ncol, betavec, popSubFromLP, subvec, nsub, zbeta);

    /* Compute exp(zbeta) */
    for (j=0, p1=zbeta, p2=expzbeta; j<nsub; j++, p1++, p2++) *p2 = exp(*p1);

    /* Get the weights. wgts will have length nsub */
    get_wgts(weights, subvec, nsub, wgts);

    /* Compute final value. ref_LP2 is for scratch space */
    retvec[ii] = finalVal(nsub, zbeta, expzbeta, age_new[ii], age_int[ii], lambda, loglambda, compRates, holder, ref_LP2, wgts);
    retlps[ii] = mean_2val(ref_full_LP, subvec, nsub, wgts);
  }

  matrix_free((void **) ZMAT, zmat_nrow);
  matrix_free((void **) REF, ref_nrow);
  R_Free(present);
  R_Free(ref_LP);
  R_Free(ref_LP2);
  R_Free(cutpoints);
  R_Free(q_lo);
  R_Free(q_hi);
  R_Free(q_h);
  R_Free(q_ivec);
  R_Free(subvec);
  R_Free(zbeta);
  R_Free(expzbeta);
  R_Free(loglambda);
  R_Free(holder);
  R_Free(wgts);

  *retflag = 0;

  return;

} /* END: ref_risk2 */

void R_init_iCARE(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
}




