#include "supermatrix.h"
#include "util.h"

void
sp_preorder(char *refact,  SuperMatrix *A, int *perm_c, 
	    int *etree, SuperMatrix *AC)
{
/*
 * Purpose
 * =======
 *
 * sp_preorder() permutes the columns of the original matrix. It performs
 * the following steps:
 *
 *    1. Apply column permutation perm_c[] to A's column pointers to form AC;
 *
 *    2. If refact = 'N', then
 *       (1) Compute column elimination tree etree[] of AC'AC;
 *       (2) Post order etree[] to get a postordered elimination tree etree[],
 *           and a postorder permutation post[];
 *       (3) Apply post[] permutation to columns of AC;
 *       (4) Overwrite perm_c[] with the product perm_c * post.
 *
 * Arguments
 * =========
 *
 * refact (input) char *
 *         Specifies whether or not the elimination tree will be re-used.
 *         = 'N': first time factor A, etree is computed and output.
 *         = 'Y': re-factor A, etree is input, unchanged on exit.
 *
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of the linear equations is A->nrow. Currently, the type of A can be:
 *         Stype = NC or NCP; Dtype = _D; Mtype = GE. In the future,
 *         more general A can be handled.
 *
 * perm_c  (input/output) int*
 *	   Column permutation vector of size A->ncol, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 *
 * etree   (input/output) int*
 *         Elimination tree of Pc'*A'*A*Pc, dimension A->ncol.
 *         If fact is not 'F' and refact = 'Y', etree is an input argument,
 *         otherwise it is an output argument.
 *         Note: etree is a vector of parent pointers for a forest whose
 *         vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
 *
 * AC      (output) SuperMatrix*
 *         The resulting matrix after applied the column permutation
 *         perm_c[] to matrix A. The type of AC can be:
 *         Stype = NCP; Dtype = D; Mtype = GE.
 *
 */

    NCformat  *Astore;
    NCPformat *ACstore;
    int       *iwork, *post;
    register  int n, i;

    n = A->ncol;
    iwork = (int*) SUPERLU_MALLOC((n+1)*sizeof(int)); 
    if ( !iwork ) ABORT("SUPERLU_MALLOC fails for iwork[]");
    
    /* Apply column permutation perm_c to A's column pointers so to
       obtain NCP format in AC = A*Pc.  */
    AC->Stype       = NCP;
    AC->Dtype       = A->Dtype;
    AC->Mtype       = A->Mtype;
    AC->nrow        = A->nrow;
    AC->ncol        = A->ncol;
    Astore          = A->Store;
    ACstore = AC->Store = (void *) SUPERLU_MALLOC( sizeof(NCPformat) );
    if ( !ACstore ) ABORT("SUPERLU_MALLOC fails for ACstore");
    ACstore->nnz    = Astore->nnz;
    ACstore->nzval  = Astore->nzval;
    ACstore->rowind = Astore->rowind;
    ACstore->colbeg = (int*) SUPERLU_MALLOC(n*sizeof(int));
    if ( !(ACstore->colbeg) ) ABORT("SUPERLU_MALLOC fails for ACstore->colbeg");
    ACstore->colend = (int*) SUPERLU_MALLOC(n*sizeof(int));
    if ( !(ACstore->colend) ) ABORT("SUPERLU_MALLOC fails for ACstore->colend");

#ifdef DEBUG
    print_int_vec("pre_order:", n, perm_c);
    check_perm("Initial perm_c", n, perm_c);
#endif      

    for (i = 0; i < n; i++) {
	ACstore->colbeg[perm_c[i]] = Astore->colptr[i]; 
	ACstore->colend[perm_c[i]] = Astore->colptr[i+1];
    }
	
    if ( lsame_(refact, "N") ) {
	
	/* Compute the column elimination tree */
	sp_coletree(ACstore->colbeg, ACstore->colend, ACstore->rowind,
		    A->nrow, A->ncol, etree);
#ifdef DEBUG	
	print_int_vec("etree:", n, etree);
#endif	
	
	/* Post order etree */
	post = (int *) TreePostorder(n, etree);
	/* for (i = 0; i < n+1; ++i) inv_post[post[i]] = i;
	   iwork = post; */

#ifdef DEBUG
	print_int_vec("post:", n+1, post);
	check_perm("post", n, post);	
#endif	

	/* Renumber etree in postorder */
	for (i = 0; i < n; ++i) iwork[post[i]] = post[etree[i]];
	for (i = 0; i < n; ++i) etree[i] = iwork[i];

#ifdef DEBUG	
	print_int_vec("postorder etree:", n, etree);
#endif
	
	/* Postmultiply A*Pc by post[] */
	for (i = 0; i < n; ++i) iwork[post[i]] = ACstore->colbeg[i];
	for (i = 0; i < n; ++i) ACstore->colbeg[i] = iwork[i];
	for (i = 0; i < n; ++i) iwork[post[i]] = ACstore->colend[i];
	for (i = 0; i < n; ++i) ACstore->colend[i] = iwork[i];

	for (i = 0; i < n; ++i)
	    iwork[i] = post[perm_c[i]];  /* product of perm_c and post */
	for (i = 0; i < n; ++i) perm_c[i] = iwork[i];

#ifdef DEBUG
	print_int_vec("Pc*post:", n, perm_c);
	check_perm("final perm_c", n, perm_c);	
#endif

	SUPERLU_FREE (post);

    } /* if refact = 'N' */

    SUPERLU_FREE (iwork);

}

check_perm(char *what, int n, int *perm)
{
    register int i;
    int          *marker;
    marker = (int *) calloc(n, sizeof(int));

    for (i = 0; i < n; ++i) {
	if ( marker[perm[i]] == 1 || perm[i] >= n ) {
	    printf("%s: Not a valid PERM[%d] = %d\n", what, i, perm[i]);
	    ABORT("check_perm");
	} else {
	    marker[perm[i]] = 1;
	}
    }

    SUPERLU_FREE(marker);
    return 0;
}
