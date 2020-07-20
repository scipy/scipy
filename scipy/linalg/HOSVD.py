def HOSVD(A):
    #A is n(1) x ... x n(d) tensor
    #U is d cell array with U(k) being the left singular vector of a's mode-k unfolding
    #S is n(1) x ... x n(d) tensor : A x1 U(1) x2 U(2) ... xd U(d)
    S = A

    for k in range(A.size()):
        C = A[k,:].numpy() ##may need to make sure columns are arranged in increasing order
        [U[k],Sig,V] = svd(C)
        t_U(k) = numpyu.transpose(U(k))
        S = numpy.tensordot(S,t_U(k),k)
    return S, U

def SHOSVD(A):
    #A is n(1) x ... x n(d) tensor
    #U is d cell array with U(k) being the left singular vector of a's mode-k unfolding
    #S is n(1) x ... x n(d) tensor : A x1 U(1) x2 U(2) ... xd U(d)
    S = A

    for k in range(A.size()):
        C = A[k,:].numpy() ##may need to make sure columns are arranged in increasing order
        [U[k],Sig,V] = svds(C)
        t_U(k) = numpyu.transpose(U(k))
        S = numpy.tensordot(S,t_U(k),k)
    return S, U