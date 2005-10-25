/**************************************************************************
*                                                                         *
*               Swiss Federal Institute of Technology (ETH)               *
*                       CH-8092 Zuerich, Switzerland                      *
*                                                                         *
*                       (C) 1999 All Rights Reserved                      *
*                                                                         *
*                                NOTICE                                   *
*                                                                         *
*  Permission to use, copy, modify, and distribute this software and      *
*  its documentation for any purpose and without fee is hereby granted    *
*  provided that the above copyright notice appear in all copies and      *
*  that both the copyright notice and this permission notice appear in    *
*  supporting documentation.                                              *
*                                                                         *
*  Neither the Swiss Federal Institute of Technology nor the author make  *
*  any representations about the suitability of this software for any     *
*  purpose.  This software is provided ``as is'' without express or       *
*  implied warranty.                                                      *
*                                                                         *
***************************************************************************
*
*  Main routine of eigensolver
*
*  $Id: jdsym.c,v 1.2 2005/01/04 14:28:01 hamsel Exp $
*
*
**************************************************************************/


/****************************************************************************
 *                                                                          *
 * Prototypes of static functions                                           *
 *                                                                          *
 ****************************************************************************/
static void print_status(int clvl, int it, int k, int j, int kmax, 
			 int blksize, int actblksize,
			 double *s, double *resnrm, int *actcorrits);
static void quicksort(int n, double arr[], int idx[]);
static void sorteig(int j, double S[], double U[], int ldu, double tau,
		    double dtemp[], int idx1[], int idx2[], int strategy);


/****************************************************************************
 *                                                                          *
 * Main eigensolver routine                                                 *
 *                                                                          *
 ****************************************************************************/

static int
jdsym (int n, double tau, double jdtol, 
       int kmax, int jmax, int jmin, int itmax,
       int blksize, int blkwise, 
       PyArrayObject *V0, 
       PyObject *linsolver,
       int optype, 
       int linitmax, double eps_tr, double toldecay,
       int strategy, int clvl,
       int *k_conv, double *Q, double *lambda, int *it, int *it_inner,
       PyObject *amat,
       PyObject *mmat,
       PyObject *prec,
       PyObject *proj) {
  
  /****************************************************************************
   *                                                                          *
   * Local variables                                                          *
   *                                                                          *
   ****************************************************************************/

  /* Correction equation object */
  CorrEqObject *correq = NULL;

  /* allocatables: 
   * initialize with NULL, so we can free even unallocated ptrs */
  double *V = NULL, *Vtmp = NULL, *U = NULL, *Qm = NULL, *Y = NULL, 
    *s = NULL, *Res = NULL, *resnrm = NULL, 
    *M = NULL, *H = NULL, *Hlu = NULL, 
    *eigwork = NULL, *temp1 = NULL, *temp2 = NULL;

  int *Hpiv = NULL, *idx1 = NULL, *idx2 = NULL, 
    *convind = NULL, *keepind = NULL, *solvestep = NULL, 
    *actcorrits = NULL;

  /* non-allocated ptrs */
  double *q, *v, *u, *r, *y, *qm;
  double *matdummy, *vecdummy;

  /* scalar vars */
  double theta, alpha, it_tol;

  int eigworklen, found, conv, keep, act, cnt, idummy, info, ret;
  double linres;		/* return values from linear solver */
  int linit;
  int k;			/* number of already converged eigenpairs */
  int j;			/* dimension of serach subspace (matrix V) */
  int nof_ce_its = 0;		/* number of inner iterations spent in CE */
  int nof_ce = 0;		/* number of CE solved */
  int actblksize;		/* current block size (may be smaller than blksize) */
  int i;			/* look variable */

  /* variables for random number generator */
  int IDIST = 1;
  int ISEED[4] = {2, 3, 5, 7};

  /****************************************************************************
   *                                                                          *
   * Execution starts here...                                                 *
   *                                                                          *
   ****************************************************************************/

  /* print info header */
  if (clvl >= 1) {
    printf("JDSYM     Solving  %s  %s preconditioning.\n\n",
	   !mmat ? "A*x = lambda*x" : "A*x = lambda*M*x",
	   !prec ? "without" : "with");
    printf("  N=      %10d  ITMAX=%4d\n", n, itmax);
    printf("  KMAX=%3d  JMIN=%3d  JMAX=%3d  V0DIM=%3d\n", 
	   kmax, jmin, jmax, 
	   V0 == NULL ? 0 : (V0->nd == 1 ? 1 : V0->dimensions[1]));
    printf("  BLKSIZE=        %2d  BLKWISE=      %5s\n", 
	   blksize, blkwise ? "TRUE" : "FALSE");
    printf("  TAU=   %11.4e  JDTOL=  %11.4e  STRATEGY= %8d\n", 
	   tau, jdtol, strategy);
    printf("  OPTYPE=%12s\n", 
	   optype == JDSYM_OP_UNSYM ? "UNSYM" : "SYM");
    printf("  LINITMAX=    %5d  EPS_TR=  %10.3e  TOLDECAY=%9.2e\n", 
	   linitmax, eps_tr, toldecay);
    printf("\n");
  }

  /* validate input parameters */
  assert(0 < jdtol);
  assert(0 < kmax && kmax <= n);
  assert(0 < jmax && jmax <= n);
  assert(0 < jmin && jmin < jmax);
  assert(0 <= itmax);
  assert(blksize <= jmin);
  assert(blksize <= jmax - jmin);
  assert(0 < blksize && blksize <= kmax);
  assert(blkwise == 0 || blkwise == 1);
  assert(optype == JDSYM_OP_UNSYM || optype == JDSYM_OP_SYM);
  assert(0 <= linitmax);
  assert(0.0 <= eps_tr);
  assert(1.0 < toldecay);
  
  /* Get hardware-dependent values:
   * Opt size of workspace for DSYEV is (NB+2)*j, where NB is the opt.
   * block size... */
  eigworklen = (2 + F77(ilaenv)(&ONE, "dsytrd", "vu", &jmax, &MONE, &MONE, &MONE, 6, 2)) * jmax;

  /*
   * Allocate memory
   */

  V = (double *)malloc(n * jmax * sizeof(double));
  U = (double *)malloc(jmax * jmax * sizeof(double));
  s = (double *)malloc(jmax * sizeof(double));
  Res = (double *)malloc(n * blksize * sizeof(double));
  resnrm = (double *)malloc(blksize * sizeof(double));
  M = (double *)malloc(jmax * jmax * sizeof(double));
  Vtmp = (double *)malloc(jmax * jmax * sizeof(double));
  idx1 = (int *)malloc(jmax * sizeof(int));
  idx2 = (int *)malloc(jmax * sizeof(int));
  convind = (int *)malloc(blksize * sizeof(int));
  keepind = (int *)malloc(blksize * sizeof(int));
  solvestep = (int *)malloc(blksize * sizeof(int));
  actcorrits = (int *)malloc(blksize * sizeof(int));
  eigwork = (double *)malloc(eigworklen * sizeof(double));
  temp1 = (double *)malloc(n * sizeof(double));
  temp2 = (double *)malloc(n * sizeof(double));

  if (!(V && U && s && Res && resnrm && M && Vtmp && idx1 && idx2 && convind && 
	keepind && solvestep && actcorrits && eigwork && temp1 && temp2)) {
    PyErr_NoMemory();
    goto fail;
  }
    
  /* Allocate matrices H, Y and G only if necessary */
  if (prec) {
    H = (double *)malloc(kmax * kmax * sizeof(double));
    Hlu = (double *)malloc(kmax * kmax * sizeof(double));
    Hpiv = (int *)malloc(kmax * sizeof(int));
    Y = (double *)malloc(n * kmax * sizeof(double));
    if (!(H && Hlu && Hpiv && Y)) {
      PyErr_NoMemory();
      goto fail;
    }
  }
  
  /* Analogous for Qm only if necessary */
  if (mmat) { 
    Qm = (double *)malloc(n * kmax * sizeof(double));
    if (!Qm) {
      PyErr_NoMemory();
      goto fail;
    }
  }

  /* Create CorrEqObject */
  correq = (CorrEqObject *)newCorrEqObject(optype, n, 
					   amat, mmat, prec,
					   Q, Qm, Y, Hpiv, Hlu, kmax);
  if (!correq)
    goto fail;

  /**************************************************************************
   *                                                                        *
   * Generate initial search subspace V. Vectors are taken from V0 and if   *
   * necessary randomly generated.                                          *
   *                                                                        *
   **************************************************************************/

  /* copy V0 to V */
  if (V0 != NULL) {
    if (V0->nd == 1) {
      j = 1;
      for (i = 0; i < n; i ++)
	V[i] = *(double *)(V0->data + i*V0->strides[0]);
    }
    else {
      j = V0->dimensions[1];
      if (j > jmax)
	j = jmax;
      for (k = 0; k < j; k ++)
	for (i = 0; i < n; i ++)
	  V[n*k + i] = *(double *)(V0->data + i*V0->strides[0] + k*V0->strides[1]);
    }
  } else
    j = 0;
  
  /* if j < blksize: generate additional random vectors */
  if (j < blksize) {
    idummy = (blksize - j)*n; /* nof random numbers */
    F77(dlarnv)(&IDIST, ISEED, &idummy, V + j*n);
    j = blksize;
  }
  /* Project into user subspace */
  if (proj)
    for (cnt = 0; cnt < j; cnt ++)
      Jdsym_Proj(proj, n, V + cnt*n);
  /* (M-)orthogonalize columns of V */
  if (!mmat) {
    for (cnt = 0; cnt < j; cnt ++) {
      mgs(V + cnt*n, n, cnt, V);
      alpha = 1.0 / F77(dnrm2)(&n, V + cnt*n, &ONE);
      F77(dscal)(&n, &alpha, V + cnt*n, &ONE);
    }
  }
  else {
    for (cnt = 0; cnt < j; cnt ++) {
      icgsm(V + cnt*n, &alpha, n, cnt, V, mmat, temp1, temp2);
      alpha = 1.0/alpha;
      assert(alpha > 0.0);
      F77(dscal)(&n, &alpha, V + cnt*n, &ONE);
    }
  }
  
  /* Generate interaction matrix M = V'*A*V. Only the upper triangle
     is computed. */
  for (cnt = 0; cnt < j; cnt++){
    ret = SpMatrix_Matvec(amat, n, V+cnt*n, n, temp1);
    assert(ret == 0);
    idummy = cnt+1;
    F77(dgemv)("t", &n, &idummy, &DONE, V, &n, temp1, &ONE, 
	       &DZER, M+cnt*jmax, &ONE, 1);
  }

  /* Other initializations */
  k = 0; *it = 0;
  actblksize = blksize; 
  for(act = 0; act < blksize; act ++)
    solvestep[act] = 1;


  /****************************************************************************
   *                                                                          *
   * Main JD-iteration loop                                                   *
   *                                                                          *
   ****************************************************************************/

  while(*it < itmax) {

    /****************************************************************************
     *                                                                          *
     * Solving the projected eigenproblem                                       *
     *                                                                          *
     * M*u = V'*A*V*u = s*u                                                     *
     * M is symmetric, only the upper triangle is stored                        *
     *                                                                          *
     ****************************************************************************/

    F77(dlacpy)("u", &j, &j, M, &jmax, U, &jmax, 1);
    F77(dsyev)("v", "u", &j, U, &jmax, s, eigwork, &eigworklen, &info, 1, 1);
    if (info != 0) {
      printf("jdsym: error solving the projected eigenproblem.");
      printf(" dsyev: info = %d\n", info);
    }
    assert(info == 0);
  
    /* sort eigenpairs, such that |S(i)-tau| <= |S(i+1)-tau| for i=1..j-1 */
    sorteig(j, s, U, jmax, tau, temp1, idx1, idx2, strategy);

    /****************************************************************************
     *                                                                          *
     * Convergence/Restart Check                                                *
     *                                                                          *
     * In case of convergence, strip off a whole block or just the converged    *
     * ones and put 'em into Q.  Update the matrices Q, V, U, s and if          *
     * necessary Qm, Y, and H.                                                  *
     *                                                                          *
     * In case of a restart update the V, U and M matrices and recompute the    *
     * Eigenvectors                                                             *
     *                                                                          *
     ****************************************************************************/

    found = 1;
    while(found) {

      /* conv/keep = Number of converged/non-converged Approximations */
      conv = 0; keep = 0;

      for(act=0; act < actblksize; act++) {

	/* Setting pointers for single vectors */
	q = Q + (act+k)*n; 
	u = U + act*jmax; 
	r = Res + act*n; 
	qm = Qm + (act+k)*n; 
	y = Y + (act+k)*n;
	
	/* Compute Ritz-Vector Q[:,k+cnt1]=V*U[:,cnt1] */
	theta = s[act];
	F77(dgemv)("n", &n, &j, &DONE, V, &n, u, &ONE, &DZER, q, &ONE, 1);

	/* Compute the residual and update the matrix Qm if necessary. */
	if (!mmat){ 
	  /* M is Identity */
	  SpMatrix_Matvec(amat, n, q, n, r);
	  theta = -theta;
	  F77(daxpy)(&n, &theta, q, &ONE, r, &ONE);
	}
	else{ 
	  /* M is NOT Identity */
	  SpMatrix_Matvec(mmat, n, q, n, qm);
	  SpMatrix_Matvec(amat, n, q, n, r);
	  theta = -theta;
	  F77(daxpy)(&n, &theta, qm, &ONE, r, &ONE);
	}

	/* Finally update matrices H, Y  if necessary */
	if (prec){ 
	  
	  if (mmat) {
	    matdummy = Qm; vecdummy = qm;    /* If M exists, then also Qm does */
	  }
	  else {
	    matdummy = Q; vecdummy = q;      /* Without M, no Qm */
	  }

	  /* Calculate y = inv(K)*qm */
	  SpMatrix_Precon(prec, n, vecdummy, y);

	  /* update H, starting with the column ... */
	  idummy=k+act+1;
	  F77(dgemv)("t", &n, &idummy, &DONE, matdummy, &n, y, &ONE, &DZER, 
		     H+(k+act)*kmax, &ONE, 1);

	  /* ... and then the row */
	  F77(dgemv)("t", &n, &idummy, &DONE, Y, &n, vecdummy, &ONE, &DZER, 
		     H+(k+act), &kmax, 1);
	}

	/* Compute norm of the residual and update arrays convind/keepind*/
	resnrm[act] = F77(dnrm2)(&n, r, &ONE);
	if (resnrm[act] < jdtol)
	  { convind[conv] = act; conv = conv + 1; }
	else
	  { keepind[keep] = act; keep = keep + 1; }
	
      }  /* for(act = 0; act < actblksize; act ++) */

      /* Check whether the blkwise-mode is chosen and ALL the
	 approximations converged, or whether the strip-off mode is
	 active and SOME of the approximations converged */

      found = ((blkwise==1 && conv==actblksize) || (blkwise==0 && conv!=0)) 
	&& (j > actblksize || k == kmax - actblksize);
      
      /***************************************************************************
	*                                                                        *
	* Convergence Case                                                       *
	*                                                                        *
	* In case of convergence, strip off a whole block or just the converged  *
	* ones and put 'em into Q.  Update the matrices Q, V, U, s and if        *
	* necessary Qm, Y, and H.                                                *
	*                                                                        *
	**************************************************************************/

      if (found) {

	/* Store Eigenvalues */
	for(act = 0; act < conv; act++)
	  lambda[k+act] = s[convind[act]];
	 
	/* Re-use non approximated Ritz-Values */
	for(act = 0; act < keep; act++)
	  s[act] = s[keepind[act]];

	/* Shift the others in the right position */
	for(act = 0; act < (j-actblksize); act ++)
	  s[act+keep] = s[act+actblksize];

	/* Update V. Re-use the V-Vectors not looked at yet. */
	idummy = j - actblksize;
	for (act = 0; act < n; act = act + jmax) {
	  cnt = act + jmax > n ? n-act : jmax;
	  F77(dlacpy)("a", &cnt, &j, V+act, &n, Vtmp, &jmax, 1);
	  F77(dgemm)("n", "n", &cnt, &idummy, &j, &DONE, Vtmp, 
		     &jmax, U+actblksize*jmax, &jmax, &DZER, V+act+keep*n, &n, 1, 1);
	}

	/* Insert the not converged approximations as first columns in V */
	for(act = 0; act < keep; act++){
	  F77(dlacpy)("a",&n,&ONE,Q+(k+keepind[act])*n,&n,V+act*n,&n,1);
	}

	/* Store Eigenvectors */
	for(act = 0; act < conv; act++)
	  F77(dlacpy)("a",&n,&ONE,Q+(k+convind[act])*n,&n,Q+(k+act)*n,&n,1);

	/* Update Qm if necessary */
	if (mmat){
	  for(act = 0; act < conv; act++)
	    F77(dlacpy)("a",&n,&ONE,Qm+(k+convind[act])*n,&n,Qm+(k+act)*n,&n,1);
	}

	/* Update H and Y if necessary */
	if (prec){
	  for(act = 0; act < conv; act++)   /* Y */
	    F77(dlacpy)("a",&n,&ONE,Y+(k+convind[act])*n,&n,Y+(k+act)*n,&n,1);

	  for(act=0; act < conv; act++){   /* H */
	    idummy = k + act;
	    /* Copy column ... */
	    F77(dlacpy)("a",&idummy,&ONE,H+(k+convind[act])*kmax,&kmax,H+idummy*kmax,&kmax,1);
	    /* ... diagonalelement ... */
	    H[idummy*(kmax+1)]=H[(k+convind[act])*(kmax+1)];
	    /* ... and row */
	    F77(dlacpy)("a",&ONE,&idummy,H+(k+convind[act]),&kmax,H+idummy,&kmax,1);
	  }
	}

	/* Update SearchSpaceSize j */
	j = j - conv;

	/* Let M become a diagonal matrix with the Ritzvalues as entries ... */ 
	F77(dlaset)("u", &j, &j, &DZER, &DZER, M, &jmax, 1);
	for (act = 0; act < j; act++)
	  M[act*jmax + act] = s[act];
	
	/* ... and U the Identity(jnew,jnew) */
	F77(dlaset)("a", &j, &j, &DZER, &DONE, U, &jmax, 1);

	/* Avoid computation of zero eigenvalues:

	   If STRATEGY == 1: set tau to the largest of the now
	   converged eigenvalues.

	   Warning: This may not work well if BLKSIZE > 1. */
	if (strategy == 1)
	  for(act = 0; act < conv; act ++)
	    if (lambda[k+act] > tau)
	      tau = lambda[k+act];

	/* Update Converged-Eigenpair-counter and Pro_k */
	k = k + conv;

	/* Update the new blocksize */
	actblksize = blksize < kmax-k ? blksize : kmax-k;

	/* Exit main iteration loop when kmax eigenpairs have been
           approximated */
	if (k == kmax)
	  goto end;

	/* Counter for the linear-solver-accuracy */
	for(act = 0; act < keep; act++)
	  solvestep[act] = solvestep[keepind[act]];

	for(act = keep; act < blksize; act++)
	  solvestep[act] = 1;

      } /* if(found) */
      
      /**************************************************************************
       *                                                                        *
       * Restart                                                                *
       *                                                                        *
       * The Eigenvector-Aproximations corresponding to the first jmin          *
       * Petrov-Vectors are kept.  if (j+actblksize > jmax) {                   *
       *                                                                        *
       **************************************************************************/
      if (j+actblksize > jmax) {

	idummy = j; j = jmin;

	for (act = 0; act < n; act = act + jmax) { /* V = V * U(:,1:j) */
	  cnt = act+jmax > n ? n-act : jmax;
	  F77(dlacpy)("a", &cnt, &idummy, V+act, &n, Vtmp, &jmax, 1);
	  F77(dgemm)("n", "n", &cnt, &j, &idummy, &DONE, Vtmp, 
		     &jmax, U, &jmax, &DZER, V+act, &n, 1, 1);
	}
	  
	F77(dlaset)("a", &j, &j, &DZER, &DONE, U, &jmax, 1);
	F77(dlaset)("u", &j, &j, &DZER, &DZER, M, &jmax, 1);
	for (act = 0; act < j; act++)
	  M[act*jmax + act] = s[act];
      }

    } /* while(found) */    


    /****************************************************************************
     *                                                                          *
     * Solving the correction equations                                         *
     *                                                                          *
     * Depending on the input-arguments we choose an appropriate lin.solver.    *
     *                                                                          *
     ****************************************************************************/


    /* calculate Hlu (LU-factorization), if necessary */
    if (prec){
      idummy = k + actblksize;
      F77(dlacpy)("a", &idummy, &idummy, H, &kmax, Hlu, &kmax, 1); 
      F77(dgetrf)(&idummy, &idummy, Hlu, &kmax, Hpiv, &info);
      if (info != 0)
	printf("jdsym: factorization of H failed: info=%d\n", info);
      assert(info == 0);
    }

    /* Solve actblksize times the correction equation ... */
    for (act = 0; act < actblksize; act ++) {      

      /* Setting start-value for vector v as zeros(n,1). Guarantees
         (M-)orthogonality */
      v = V + j*n;
      for (cnt = 0; cnt < n; cnt ++) 
	v[cnt] = 0.0;

      /* Adaptive accuracy and shift for the lin.solver. In case the
	 residual is big, we don't need a too precise solution for the
	 correction equation, since even in exact arithmetic the
	 solution wouldn't be too useful for the Eigenproblem. */
      r = Res + act*n;

      if (resnrm[act] < eps_tr)
	correq->update(correq, k + actblksize, s[act]);
      else
	correq->update(correq, k + actblksize, tau);
      
      it_tol = pow(toldecay, (double)(-solvestep[act]));
      solvestep[act] = solvestep[act] + 1;

      /* Form the right hand side of the correction equation */
      correq->right(correq, r);
      
      /* Solve the correction equation ... */
      ret = ItSolvers_Solve(linsolver, (PyObject *)correq, n, 
			    r, v, it_tol, linitmax, (PyObject *)correq, 
			    &info, &linit, &linres);
      if (ret == -1)
	goto fail;
      
      /* Actualizing profiling data */
/*#warning error check after iteative solver?*/
      nof_ce ++;
      nof_ce_its += linit;
      actcorrits[act] = linit;

      /* (M-)orthonormalize v to Q, project into user subspace and
         finally (M-)orthonormalize to V

	 (M-)orthonormalize v to Q is necessary, because the implicit
	 orthogonalization in the solvers may be too inaccurate.
	 IteratedCGS is used to prevent numerical breakdown */
      if (mmat) {
	mgsm(v, n, k+actblksize, Q, Qm);
	if (proj) 
	  Jdsym_Proj(proj, n, v);
	icgsm(v, &alpha, n, j, V, mmat, temp1, temp2);
      }
      else {
	mgs(v, n, k+actblksize, Q);
	if (proj) 
	  Jdsym_Proj(proj, n, v);
	icgs(v, &alpha, n, j, V, temp1);
      }

      alpha = 1.0 / alpha;
      F77(dscal)(&n, &alpha, v, &ONE);
      
      /* update interaction matrix M */
      SpMatrix_Matvec(amat, n, v, n, temp1);
      idummy = j+1;
      F77(dgemv)("t", &n, &idummy, &DONE, V, &n, temp1, &ONE, 
		 &DZER, M+j*jmax, &ONE, 1);
      
      /* Increasing SearchSpaceSize j */
      j ++;
    }   /* for (act = 0;act < actblksize; act ++) */    

    /* Print information line */
    print_status(clvl, *it, k, j - blksize, kmax, blksize, actblksize, 
		 s, resnrm, actcorrits);    

    /* Increase iteration-counter for outer loop  */
    (*it) ++;

  } /* Main iteration loop */
  
 end:

  /******************************************************************
   *                                                                *
   * Eigensolutions converged or iteration limit reached            *
   *                                                                *
   * Print statistics. Free memory. Return.                         *
   *                                                                *
   ******************************************************************/

  *k_conv = k;
  *it_inner = nof_ce_its;
  if (clvl >= 1) {
    printf("\nJDSYM execution statistics\n\n");
    printf("IT_OUTER=%d   IT_INNER_TOT=%d   IT_INNER_AVG=%8.2f   IT_INNER_PER_OUTER=%8.2f\n",
	   (*it), nof_ce_its, (double)nof_ce_its/(*it), (double)nof_ce_its/nof_ce);
    printf("\nConverged eigensolutions in order of convergence:\n");
    printf("\n  I              LAMBDA(I)      RES(I)\n");
    printf("---------------------------------------\n");
    
    for (act = 0; act < *k_conv; act ++) {
      /* Compute the residual for solution act */
      q = Q + act*n; qm = Qm + act*n;
      theta = -lambda[act];
      ret = SpMatrix_Matvec(amat, n, q, n, r);
      assert(ret == 0);
      if (!mmat)
	F77(daxpy)(&n, &theta, q, &ONE, r, &ONE);
      else { 
	ret = SpMatrix_Matvec(mmat, n, q, n, qm);
	assert(ret == 0);
	F77(daxpy)(&n, &theta, qm, &ONE, r, &ONE);
      }
      printf("%3d %22.15e %12.5e\n", act+1, lambda[act],
	     F77(dnrm2)(&n, r, &ONE));
    }
    printf("\n");
  }

  free(V); free(Vtmp); free(U); free(Qm); free(Y); 
  free(s); free(Res); free(resnrm); 
  free(M); free(H); free(Hlu); 
  free(eigwork); free(temp1); free(temp2);
  free(Hpiv); free(idx1); free(idx2); 
  free(convind); free(keepind); free(solvestep); free(actcorrits);
  Py_DECREF(correq);
  return 0;

 fail:
  free(V); free(Vtmp); free(U); free(Qm); free(Y); 
  free(s); free(Res); free(resnrm); 
  free(M); free(H); free(Hlu); 
  free(eigwork); free(temp1); free(temp2);
  free(Hpiv); free(idx1); free(idx2); 
  free(convind); free(keepind); free(solvestep); free(actcorrits);
  Py_XDECREF(correq);
  return -1;
} /* jdsym(.....) */


/****************************************************************************
 *                                                                          *
 * Supporting functions                                                     *
 *                                                                          *
 ****************************************************************************/

/* PRINT_STATUS - print status line (called for each outer iteration)
 */
static void print_status(int clvl, int it, int k, int j, int kmax, 
			 int blksize, int actblksize,
			 double *s, double *resnrm, int *actcorrits) {
  const int max_vals = 5;

  int i, idummy;

  if (clvl >= 1) {
    if (blksize == 1) {
      if (it == 0) {
	printf("  IT   K   J       RES CGIT RITZVALS(1:5)\n");
	idummy = 27 + ( 13 > max_vals*10 ? 13 : max_vals*10);
	for (i = 0; i < idummy; i ++)
	  putchar('-');
	printf("\n");
      }
      printf("%4d %3d %3d %9.2e %4d", it + 1, k, j, resnrm[0], actcorrits[0]);
      for (i = 0; i < (j < max_vals ? j : max_vals); i ++)
	printf(" %9.2e", s[i]);
      printf("\n");
    }
    else {			/* blksize > 1 */
      if (it == 0) {
	printf("  IT   K   J RITZVALS ");
	for (i = 1; i < actblksize; i ++)
	  printf("          ");
	printf("   RES      ");
	for (i = 1; i < actblksize; i ++)
	  printf("          ");
	printf("   CGIT\n");
	idummy = 12 + 4 + blksize*(10 + 10 + 5);
	for (i = 0; i < idummy; i ++)
	  putchar('-');
	printf("\n");
      }
      printf("%4d %3d %3d", it + 1, k, j);
      for (i = 0; i < blksize; i ++)
	if (i < actblksize)
	  printf(" %9.2e", s[i]);
	else
	  printf("          ");
      printf("  ");
      for (i = 0; i < blksize; i ++)
	if (i < actblksize)
	  printf(" %9.2e", resnrm[i]);
	else
	  printf("          ");
      printf("  ");
      for (i = 0; i < blksize; i ++)
	if (i < actblksize)
	  printf(" %4d", actcorrits[i]);
	else
	  printf("     ");
      printf("\n");
    }
  }
}

/*
 * SORTEIG
 *
 * Default behaviour (strategy == 0):
 *
 *   Sort eigenpairs (S(i),U(:,i)), such that 
 *
 *       |S(i) - tau| <= |S(i+1) -tau| for i=1..j-1.
 *
 *     j  : dimension of S
 *     ldu: leading dimension of U
 *   dtemp: double array of length j
 *     idx: int array of length j
 *
 * Alternate behaviour (strategy == 1):
 *
 *   Same as above but put all S(i) < tau to the end. This is used to
 *   avoid computation of zero eigenvalues.
 */

static void sorteig(int j, double S[], double U[], int ldu, double tau,
		    double dtemp[], int idx1[], int idx2[], int strategy)
{
  int i;

  /* setup vector to be sorted and index vector */
  switch (strategy) {
  case 0:
    for (i = 0; i < j; i ++)
      dtemp[i] = fabs(S[i] - tau);
    break;
  case 1:
    for (i = 0; i < j; i ++)
      if (S[i] < tau)
	dtemp[i] = DBL_MAX;
      else
	dtemp[i] = fabs(S[i] - tau);
    break;
  default:
    assert(0);
  }
  for (i = 0; i < j; i ++)
    idx1[i] = i;

  /* sort dtemp in ascending order carrying itemp along */
  quicksort(j, dtemp, idx1);

  /* compute 'inverse' index vector */
  for (i = 0; i < j; i ++)
    idx2[idx1[i]] = i;

  /* sort eigenvalues */
  memcpy(dtemp, S, j * sizeof(double));
  for (i = 0; i < j; i ++)
    S[i] = dtemp[idx1[i]];

  /* sort eigenvectors (in place) */
  for (i = 0; i < j; i ++) {
    if (i != idx1[i]) {
      memcpy(dtemp, U+i*ldu, j*sizeof(double));
      memcpy(U+i*ldu, U+idx1[i]*ldu, j*sizeof(double));
      memcpy(U+idx1[i]*ldu, dtemp, j*sizeof(double));
      idx1[idx2[i]] = idx1[i];
      idx2[idx1[i]] = idx2[i];
    }
  }
}


/* 
 * QUICKSORT 
 *
 * Sorts a double array using a non-recursive quicksort algorithm in
 * ascending order carrying along an int array.
 *  
 */

static void quicksort(int n, double arr[], int idx[])
{
  double v, td;
  int i, j, l, r, ti, tos, stack[32];
  
  l = 0; r = n-1; tos = -1;
  for (;;)
    {
      while (r > l)
	{ 
	  v = arr[r]; i = l; j = r-1;
	  for (;;)
	    { 
	      while (arr[i] < v) i ++;
	      /* j > l prevents underflow */
	      while (arr[j] >= v && j > l) j --;
	      if (i >= j) break;
	      td = arr[i]; arr[i] = arr[j]; arr[j] = td;
	      ti = idx[i]; idx[i] = idx[j]; idx[j] = ti;
	    }
	  td = arr[i]; arr[i] = arr[r]; arr[r] = td;
	  ti = idx[i]; idx[i] = idx[r]; idx[r] = ti;
	  if (i-l > r-i)
	    { stack[++tos] = l; stack[++tos] = i-1; l = i+1; }
	  else
	    { stack[++tos] = i+1; stack[++tos] = r; r = i-1; }
	  assert(tos < 32);
	} 
      if (tos == -1) break;
      r = stack[tos--]; l = stack[tos--]; 
    }
}
