# anova() and its supporting functions have been removed from scipy.stats to
# this module. No care has been taken to ensure that it works in its current
# state. At minimum, you will have to figure out what it needs to import to work.

from numpy import sum

def anova(data, effects=['A','B','C','D','E','F','G','H','I','J','K']):
    """
Prints the results of single-variable between- and within-subject ANOVA
designs.  The function can only handle univariate ANOVAs with a single
random factor.  The random factor is coded in column 0 of the input
list/array (see below) and the measured variable is coded in the last
column of the input list/array. The following were used as references
when writing the code:

Maxwell, SE, Delaney HD (1990)  Designing Experiments and Analyzing
    Data, Wadsworth: Belmont, CA.
Lindman, HR (1992) Analysis of Variance in Experimental Design,
    Springer-Verlag: New York.

TO DO:  Increase Current Max Of 10 Levels Per W/I-Subject Factor
        Consolidate Between-Subj Analyses For Between And Within/Between
        Front-end for different input data-array shapes/organization
        Axe mess of 'global' statements (particularly for d_restrict fcns)

Usage:   anova(data,                         data = |Stat format
               effects=['A','B','C','D','E','F','G','H','I','J','K'])

Note: |Stat format is as follows ... one datum per row, first element of
row is the subject identifier, followed by all within/between subject
variable designators, and the measured data point as the last element in the
row.  Thus, [1, 'short', 'drugY', 2, 14.7] represents subject 1 when measured
in the short / drugY / 2 condition, and subject 1 gave a measured value of
14.7 in this combination of conditions.  Thus, all input lists are '2D'
lists-of-lists.
"""
    global alluniqueslist, Nlevels, Nfactors, Nsubjects, Nblevels, Nallsources
    global Bscols, Bbetweens, SSlist, SSsources, Bwonly_sources, D
    global alleffects, alleffsources
    outputlist = []
    #SSbtw = []
    #SSbtwsources = []
    #SSwb = []
    #SSwbsources = []
    alleffects = []
    alleffsources = []
    SSlist = []
    SSsources = []

    print
    variables = 1       # this function only handles one measured variable
    data = asarray(data)
    if not isscalar(data):
        data = data.tolist()

## Create a list of all unique values in each column, and a list of these Ns
    alluniqueslist = [0]*(len(data[0])-variables) # all cols but data cols
    Nlevels = [0]*(len(data[0])-variables)        # (as above)
    for column in range(len(Nlevels)):
        alluniqueslist[column] = _support.unique(_support.colex(data,column))
        Nlevels[column] = len(alluniqueslist[column])

    #Ncells = multiply.reduce(Nlevels[1:]) # total num cells (w/i AND btw)
    Nfactors = len(Nlevels[1:])             # total num factors
    Nallsources = 2**(Nfactors+1)  # total no. possible sources (factor-combos)
    Nsubjects = len(alluniqueslist[0])  # total # subj in study (# of diff. subj numbers in column 0)

## Within-subj factors defined as those where there are fewer subj than
## scores in the first level of a factor (quick and dirty; findwithin() below)
    Bwithins = findwithin(data)         # binary w/i subj factors (excl. col 0)
    Bbetweens = ~Bwithins & (Nallsources-1) - 1

    Wcolumns = makelist(Bwithins,Nfactors+1)  # get list of cols of w/i factors
    Wscols = [0] + Wcolumns                   # w/i subj columns INCL col 0
    Bscols = makelist(Bbetweens+1,Nfactors+1) #list of btw-subj cols,INCL col 0
    Nwifactors = len(Wscols) - 1 # WAS len(Wcolumns)
    #Nwlevels = take(array(Nlevels),Wscols,axis=0) # no.lvls for each w/i subj fact
    #Nbtwfactors = len(Bscols) - 1 # WASNfactors - Nwifactors + 1
    Nblevels = take(array(Nlevels),Bscols,0)

    Nwsources = 2**Nwifactors - 1 # num within-subject factor-combos
    #Nbsources = Nallsources - Nwsources

    #
    # CALC M-VARIABLE (LIST) and Marray/Narray VARIABLES (ARRAY OF CELL MNS/NS)
    #
    # Eliminate replications for the same subject in same condition as well as
    # within-subject repetitions, keep as list
    M = _support.collapse(data,Bscols,-1,0,0)
    # Create an arrays of Nblevels shape (excl. subj dim)
    Marray = zeros(Nblevels[1:],'f')
    Narray = zeros(Nblevels[1:],'f')
    # Fill arrays by looping through all scores in the (collapsed) M
    for row in M:
        idx = []
        for i in range(len(row[:-1])):
            idx.append(alluniqueslist[Bscols[i]].index(row[i]))
        idx = idx[1:]
        Marray[idx] = Marray[idx] + row[-1]
        Narray[idx] = Narray[idx] + 1
    Marray = Marray / Narray

    #
    # CREATE DATA ARRAY, DA, FROM ORIGINAL INPUT DATA
    # (this is an unbelievably bad, wasteful data structure, but it makes lots
    # of tasks much easier; should nevertheless be fixed someday)

    # This limits the within-subject level count to 10!
    coefflist =[[[1]],
                [[-1,1]],
                [[-1,0,1],[1,-2,1]],
                [[-3,-1,1,3],[1,-1,-1,1],[-1,3,-3,1]],
                [[-2,-1,0,1,2],[2,-1,-2,-1,2],[-1,2,0,-2,1],[1,-4,6,-4,1]],
                [[-5,-3,-1,1,3,5],[5,-1,-4,-4,-1,5],[-5,7,4,-4,-7,5],
                 [1,-3,2,2,-3,1],[-1,5,-10,10,-5,1]],
                [[-3,-2,-1,0,1,2,3],[5,0,-3,-4,-3,0,5],[-1,1,1,0,-1,-1,1],
                 [3,-7,1,6,1,-7,3],[-1,4,-5,0,5,-4,1],[1,-6,15,-20,15,-6,1]],
                [[-7,-5,-3,-1,1,3,5,7],[7,1,-3,-5,-5,-3,1,7],
                 [-7,5,7,3,-3,-7,-5,7],[7,-13,-3,9,9,-3,-13,7],
                 [-7,23,-17,-15,15,17,-23,7],[1,-5,9,-5,-5,9,-5,1],
                 [-1,7,-21,35,-35,21,-7,1]],
                [[-4,-3,-2,-1,0,1,2,3,4],[28,7,-8,-17,-20,-17,-8,7,28],
                 [-14,7,13,9,0,-9,-13,-7,14],[14,-21,-11,9,18,9,-11,-21,14],
                 [-4,11,-4,-9,0,9,4,-11,4],[4,-17,22,1,-20,1,22,-17,4],
                 [-1,6,-14,14,0,-14,14,-6,1],[1,-8,28,-56,70,-56,28,-8,1]],
                [[-9,-7,-5,-3,-1,1,3,5,7,9],[6,2,-1,-3,-4,-4,-3,-1,2,6],
                 [-42,14,35,31,12,-12,-31,-35,-14,42],
                 [18,-22,-17,3,18,18,3,-17,-22,18],
                 [-6,14,-1,-11,-6,6,11,1,-14,6],[3,-11,10,6,-8,-8,6,10,-11,3],
                 [9,-47,86,-42,-56,56,42,-86,47,-9],
                 [1,-7,20,-28,14,14,-28,20,-7,1],
                 [-1,9,-36,84,-126,126,-84,36,-9,1]]]

    dindex = 0
    # Prepare a list to be filled with arrays of D-variables, array per within-
    # subject combo (i.e., for 2 w/i subj factors E and F ... E, F, ExF)
    NDs = [0]* Nwsources
    for source in range(Nwsources):
        if subset(source,Bwithins):
            NDs[dindex] = numlevels(source,Nlevels)
            dindex = dindex + 1

    # Collapse multiple repetitions on the same subject and same condition
    cdata = _support.collapse(data,range(Nfactors+1),-1,0,0)

    # Find a value that's not a data score with which to fill the array DA
    dummyval = -1
    datavals = _support.colex(data,-1)
    while dummyval in datavals:  # find a value that's not a data score
        dummyval = dummyval - 1
    DA = ones(Nlevels,'f')*dummyval # create plenty of data-slots to fill

    if len(Bscols) == 1: # ie., if no btw-subj factors
        # 1 (below) needed because we need 2D array even w/ only 1 group of subjects
        subjslots = ones((Nsubjects,1))
    else: # create array to hold 1s (subj present) and 0s (subj absent)
        subjslots = zeros(Nblevels)
    for i in range(len(data)): # for every datapoint given as input
        idx = []
        for j in range(Nfactors+1): # get n-D bin idx for this datapoint
            new = alluniqueslist[j].index(data[i][j])
            idx.append(new)
        DA[idx] = data[i][-1] # put this data point in proper place in DA
        btwidx = take(idx,array(Bscols),0)
        subjslots[btwidx] = 1
    # DONE CREATING DATA ARRAY, DA ... #dims = numfactors+1, dim 0=subjects
    # dim -1=measured values, dummyval = values used to fill empty slots in DA

    # PREPARE FOR MAIN LOOP
    dcount = -1     # prepare for pre-increment of D-variable pointer
    Bwsources = []  # binary #s, each=source containing w/i subj factors
    Bwonly_sources = [] # binary #s, each=source of w/i-subj-ONLY factors
    D = zeros(Nwsources,PyObject) # one slot for each Dx,2**Nwifactors
    DM = [0] *Nwsources # Holds arrays of cell-means
    DN = [0] *Nwsources # Holds arrays of cell-ns

    # BEGIN MAIN LOOP!!!!!
    # BEGIN MAIN LOOP!!!!!
    # BEGIN MAIN LOOP!!!!!
    for source in range(3,Nallsources,2): # all sources that incl. subjects
        if ((source-1) & Bwithins) != 0: # 1 or more w/i subj sources?
            Bwsources.append(source-1)   # add it to a list
        #
        # WITHIN-SUBJECT-ONLY TERM?  IF SO ... NEED TO CALCULATE NEW D-VARIABLE
        # (per Maxwell & Delaney pp.622-4)
        if subset((source-1),Bwithins):
            # Keep track of which D-var set we're working with (De, Df, Def, etc.)
            dcount = dcount + 1
            Bwonly_sources.append(source-1) #add source, minus subj,to list
            dwsc = 1.0 * DA       # get COPY of w/i-subj data array
            # Find all non-source columns, note ~source alone (below) -> negative number
            Bnonsource = (Nallsources-1) & ~source
            Bwscols = makebin(Wscols) # make a binary version of Wscols
            # Figure out which cols from the ORIGINAL (input) data matrix are both non-
            # source and also within-subj vars (excluding subjects col)
            Bwithinnonsource = Bnonsource & Bwscols

            # Next, make a list of the above.  The list is a list of axes in DA
            # because DA has the same number of axes as there are factors
            # (including subjects), but with extra dummyval='-1' values the original
            # data array (assuming between-subj vars exist)
            Lwithinnonsource = makelist(Bwithinnonsource,Nfactors+1)

            # Collapse all non-source, w/i subj dims, FROM THE END (otherwise the
            # dim-numbers change as you collapse).  THIS WORKS BECAUSE WE'RE
            # COLLAPSING ACROSS W/I SUBJECT AXES, WHICH WILL ALL HAVE THE
            # SAME SUBJ IN THE SAME ARRAY LOCATIONS (i.e., dummyvals will still exist
            # but should remain the same value through the mean(,axis=0) function
            for i in range(len(Lwithinnonsource)-1,-1,-1):
                dwsc = mean(dwsc,Lwithinnonsource[i])
            mns = dwsc

            # NOW, ACTUALLY COMPUTE THE D-VARIABLE ENTRIES FROM DA
            # CREATE LIST OF COEFF-COMBINATIONS TO DO (len=e-1, f-1, (e-1)*(f-1), etc...)
            #
            # Figure out which cols are both source and within-subjects, including col 0
            Bwithinsource = source & Bwscols
            # Make a list of within-subj cols, incl subjects col (0)
            Lwithinsourcecol = makelist(Bwithinsource, Nfactors+1)
            # Make a list of cols that are source within-subj OR btw-subj
            Lsourceandbtws = makelist(source | Bbetweens, Nfactors+1)
            if Lwithinnonsource != []:
                Lwithinsourcecol = map(Lsourceandbtws.index,Lwithinsourcecol)
                # Now indxlist should hold a list of indices into the list of possible
                # coefficients, one row per combo of coefficient. Next line PRESERVES dummyval
            dvarshape = array(take(mns.shape,Lwithinsourcecol[1:],0)) -1
            idxarray = indices(dvarshape)
            newshape = array([idxarray.shape[0],
                                multiply.reduce(idxarray.shape[1:])])
            indxlist = swapaxes(reshape(idxarray,newshape),0,1)

            # The following is what makes the D-vars 2D.  It takes an n-dim array
            # and retains the first (num of factors) dim while making the 2nd dim
            # equal to the total number of source within-subject cells.

            #
            # CREATE ALL D-VARIABLES FOR THIS COMBINATION OF FACTORS
            #
            for i in range(len(indxlist)):
                #
                # FILL UP COEFFMATRIX (OF SHAPE = MNS) WITH CORRECT COEFFS FOR 1 D-VAR
                #
                coeffmatrix = ones(mns.shape,Float) # fewer dims than DA (!!)
                # Make a list of dim #s that are both in source AND w/i subj fact, incl subj
                #Wsourcecol = makelist(Bwscols&source,Nfactors+1)
                # Fill coeffmatrix with a complete set of coeffs (1 per w/i-source factor)
                for wfactor in range(len(Lwithinsourcecol[1:])):
                    #put correct coeff. axis as first axis, or "swap it up"
                    coeffmatrix = swapaxes(coeffmatrix,0,
                                             Lwithinsourcecol[wfactor+1])
                    # Find appropriate ROW of (static) coefflist we need
                    nlevels = coeffmatrix.shape[0]
                    # Get the next coeff in that row
                    try:
                        nextcoeff = coefflist[nlevels-1][indxlist[i,wfactor]]
                    except IndexError:
                        raise IndexError, "anova() can only handle up to 10 levels on a within-subject factors"
                    for j in range(nlevels):
                        coeffmatrix[j] = coeffmatrix[j] * nextcoeff[j]
                    # Swap it back to where it came from
                    coeffmatrix = swapaxes(coeffmatrix,0,
                                             Lwithinsourcecol[wfactor+1])

                # CALCULATE D VARIABLE
                scratch = coeffmatrix * mns
                # Collapse all axes EXCEPT subjects dim (dim 0)
                for j in range(len(coeffmatrix.shape[1:])):
                    scratch = add.reduce(scratch,1)
                if len(scratch.shape) == 1:
                    scratch.shape = list(scratch.shape)+[1]
                try:
                    # Tack this column onto existing ones
                    #tmp = D[dcount].shape
                    D[dcount] = _support.abut(D[dcount],scratch)
                except AttributeError: # i.e., D[dcount]=integer/float
                    # If this is the first, plug it in
                    D[dcount] = scratch


            # Big long thing to create DMarray (list of DM variables) for this source
            variables = D[dcount].shape[1] # Num variables for this source
            tidx = range(1,len(subjslots.shape)) + [0] # [0] = Ss dim
            tsubjslots = transpose(subjslots,tidx) # put Ss in last dim
            DMarray = zeros(list(tsubjslots.shape[0:-1]) +
                              [variables],'f') # btw-subj dims, then vars
            DNarray = zeros(list(tsubjslots.shape[0:-1]) +
                              [variables],'f') # btw-subj dims, then vars
            idx = [0] *len(tsubjslots.shape[0:-1])
            idx[0] = -1
            loopcap = array(tsubjslots.shape[0:-1]) -1
            while incr(idx,loopcap) != -1:
                DNarray[idx] = float(sum(tsubjslots[idx],axis=0))
                thismean =  (add.reduce(tsubjslots[idx] * # 1=subj dim
                                          transpose(D[dcount]),1) /
                             DNarray[idx])
                thismean = array(thismean,PyObject)
                DMarray[idx] = thismean
            DM[dcount] = DMarray
            DN[dcount] = DNarray

        #
        # DONE CREATING M AND D VARIABLES ... TIME FOR SOME SS WORK
        # DONE CREATING M AND D VARIABLES ... TIME FOR SOME SS WORK
        #
        if Bscols[1:] != []:
            BNs = _support.colex([Nlevels],Bscols[1:])
        else:
            BNs = [1]
            #
            # FIGURE OUT WHICH VARS TO RESTRICT, see p.680 (Maxwell&Delaney)
            #
            # BETWEEN-SUBJECTS VARIABLES ONLY, use M variable for analysis
            #
        if ((source-1) & Bwithins) == 0:  # btw-subjects vars only?
            #sourcecols = makelist(source-1,Nfactors+1)

            # Determine cols (from input list) required for n-way interaction
            Lsource = makelist((Nallsources-1)&Bbetweens,Nfactors+1)
            # NOW convert this list of between-subject column numbers to a list of
            # AXES in M, since M has fewer dims than the original data array
            # (assuming within-subj vars exist); Bscols has list of between-subj cols
            # from input list, the indices of which correspond to that var's loc'n in M
            btwcols = map(Bscols.index,Lsource)
            # Obviously-needed loop to get cell means is embedded in the collapse fcn, -1
            # represents last (measured-variable) column, None=std, 1=retain Ns

            #hn = ahmean(Narray,-1) # -1=unravel first

            # CALCULATE SSw ... SUBTRACT APPROPRIATE CELL MEAN FROM EACH SUBJ SCORE
            SSw = 0.0
            #idxlist = _support.unique(_support.colex(M,btwcols))
            for row in M:
                idx = []
                for i in range(len(row[:-1])):
                    idx.append(alluniqueslist[Bscols[i]].index(row[i]))
                idx = idx[1:]   # Strop off Ss col/dim
                newval = row[-1] - Marray[idx]
                SSw = SSw + (newval)**2

            # Determine which cols from input are required for this source
            Lsource = makelist(source-1,Nfactors+1)
            # NOW convert this list of between-subject column numbers to a list of
            # AXES in M, since M has fewer dims than the original data array
            # (assuming within-subj vars exist); Bscols has list of between-subj cols
            # from input list, the indices of which correspond to that var's loc'n in M
            #btwsourcecols = (array(map(Bscols.index,Lsource))-1).tolist()

            # Average Marray and get harmonic means of Narray OVER NON-SOURCE DIMS
            Bbtwnonsourcedims = ~source & Bbetweens
            Lbtwnonsourcedims = makelist(Bbtwnonsourcedims,Nfactors+1)
            btwnonsourcedims = (array(map(Bscols.index,Lbtwnonsourcedims))-1).tolist()

    ## Average Marray over non-source axes (1=keep squashed dims)
            sourceMarray = apply_over_axes(mean, Marray,btwnonsourcedims)

    ## Calculate harmonic means for each level in source
            sourceNarray = apply_over_axes(hmean, Narray,btwnonsourcedims)

    ## Calc grand average (ga,axis=0), used for ALL effects
            ga = sum((sourceMarray*sourceNarray)/ \
                     sum(sourceNarray,axis=0),axis=0)
            ga = reshape(ga,ones(len(Marray.shape)))

    ## If GRAND interaction, use harmonic mean of ALL cell Ns
            if source == Nallsources-1:
                sourceNarray = hmean(Narray, None)

    ## Calc all SUBSOURCES to be subtracted from sourceMarray (M&D p.320)
            sub_effects = 1.0 * ga # start with grand mean
            for subsource in range(3,source,2):
        ## Make a list of the non-subsource axes
                if subset(subsource-1,source-1):
                    sub_effects = (sub_effects +
                                   alleffects[alleffsources.index(subsource)])
        ## Calc this effect (a(j)'s, b(k)'s, ab(j,k)'s, whatever)
            effect = sourceMarray - sub_effects

        ## Save it so you don't have to calculate it again next time
            alleffects.append(effect)
            alleffsources.append(source)

    ## Calc and save sums of squares for this source
            SS = sum((effect**2 *sourceNarray) *
                      multiply.reduce(take(Marray.shape,btwnonsourcedims,0)),axis=0)
        ## Save it so you don't have to calculate it again next time
            SSlist.append(SS)
            SSsources.append(source)

            collapsed = _support.collapse(M,btwcols,-1,0,1)
            # Obviously needed for-loop to get source cell-means embedded in collapse fcns
            #contrastmns = _support.collapse(collapsed,btwsourcecols,-2,1,1)
            # Collapse again, this time SUMMING instead of averaging (to get cell Ns)
            #contrastns = _support.collapse(collapsed,btwsourcecols,-1,0,0,
            #                            sum)

            # Collapse again, this time calculating hmeans (for hns)
            #contrasthns = _support.collapse(collapsed,btwsourcecols,-1,0,0,
            #                             hmean)
            # CALCULATE *BTW-SUBJ* dfnum, dfden
            sourceNs = _support.colex([Nlevels],makelist(source-1,Nfactors+1))
            dfnum = multiply.reduce(ravel(array(sourceNs)-1))
            dfden = Nsubjects - multiply.reduce(ravel(BNs))

            # CALCULATE MS, MSw, F AND PROB FOR ALL-BETWEEN-SUBJ SOURCES ONLY
            MS = SS / dfnum
            MSw = SSw / dfden
            if MSw != 0:
                f = MS / MSw
            else:
                f = 0  # i.e., absolutely NO error in the full model

            if f >= 0:
                prob = fprob(dfnum, dfden, f)
            else:
                prob = 1.0
    # Now this falls thru to output stage

    #
    # SOME WITHIN-SUBJECTS FACTORS TO DEAL WITH ... use appropriate D variable
    #
        else:  # Source has some w/i subj factors
            # FIGURE OUT WHICH D-VAR TO USE BASED ON WHICH W/I-SUBJ FACTORS ARE IN SOURCE
            # Determine which w/i-subj factors are in this source
            sourcewithins = (source-1) & Bwithins
            # Use D-var that was created for that w/i subj combo (the position of that
            # source within Bwsources determines the index of that D-var in D)
            workd = asarray(D[Bwonly_sources.index(sourcewithins)])

            # CALCULATE Er, Ef
    ## Set up workd and subjslots for upcoming calcs
            if len(workd.shape)==1:
                workd = workd[:,newaxis]
            if len(subjslots.shape)==1:
                subjslots = subjslots[:,newaxis]

    ## Calculate full-model sums of squares
            ef = d_full_model(workd,subjslots) # Uses cell-means model

            #
            # **ONLY** WITHIN-SUBJECT VARIABLES TO CONSIDER
            #
            if subset((source-1),Bwithins):
                # restrict grand mean, as per M&D p.680
                er = d_restrict_mean(workd,subjslots)
        #
        # **BOTH** WITHIN- AND BETWEEN-SUBJECTS VARIABLES TO CONSIDER
        #
            else:
                er = d_restrict_source(workd,subjslots,source) + ef
            SSw = linalg.det(ef)
            SS = linalg.det(er) - SSw

        # CALCULATE *W/I-SUBJ* dfnum, dfden
            sourceNs = _support.colex([Nlevels],makelist(source,Nfactors+1))
            # Calculation of dfnum is straightforward regardless
            dfnum = multiply.reduce(ravel(array(sourceNs)-1)[1:])
            # If only within-subject factors are involved, dfden is straightforward
            if subset(source-1,Bwithins):
                dfden = Nsubjects -multiply.reduce(ravel(BNs))-dfnum +1
                MS = SS / dfnum
                MSw = SSw / dfden
                if MSw != 0:
                    f = MS / MSw
                else:
                    f = 0  # i.e., absolutely NO error in full model

                if f >= 0:
                    prob = fprob(dfnum, dfden, f)
                else:
                    prob = 1.0

            # If combined within-between source, must use Rao's approximation for dfden
            # from Tatsuoka, MM (1988) Multivariate Analysis (2nd Ed), MacMillan: NY p93
            else: # it's a within-between combo source
                try:
                    p = workd.shape[1]
                except IndexError:
                    p = 1
                k = multiply.reduce(ravel(BNs))
                m = Nsubjects -1 -(p+k)/2.0
                d_en = float(p**2 + (k-1)**2 - 5)
                if d_en == 0.0:
                    s = 1.0
                else:
                    s = math.sqrt(((p*(k-1))**2-4) / d_en)
                dfden = m*s - dfnum/2.0 + 1

                # Given a within-between combined source, Wilk's Lambda is appropriate
                if linalg.det(er) != 0:
                    lmbda = linalg.det(ef) / linalg.det(er)
                    W = math.pow(lmbda,(1.0/s))
                    f = ((1.0-W)/W) * (dfden/dfnum)
                else:
                    f = 0  # i.e., absolutely NO error in restricted model

                if f >= 0:
                    prob = fprob(dfnum,dfden,f)
                else:
                    prob = 1.0

        #
        # CREATE STRING-LIST FOR RESULTS FROM THIS PARTICULAR SOURCE
        #
        suffix = ''                       # for *s after the p-value
        if  prob < 0.001:  suffix = '***'
        elif prob < 0.01:  suffix = '**'
        elif prob < 0.05:  suffix = '*'
        adjsourcecols = array(makelist(source-1,Nfactors+1)) -1
        thiseffect = ''
        for col in adjsourcecols:
            if len(adjsourcecols) > 1:
                thiseffect = thiseffect + effects[col][0]
            else:
                thiseffect = thiseffect + (effects[col])
        outputlist = (outputlist
        # These terms are for the numerator of the current effect/source
                      + [[thiseffect, round4(SS),dfnum,
                          round4(SS/float(dfnum)),round4(f),
                          round4(prob),suffix]]
        # These terms are for the denominator for the current effect/source
                      + [[thiseffect+'/w', round4(SSw),dfden,
                          round4(SSw/float(dfden)),'','','']]
                      + [['\n']])

        #
        # PRINT OUT ALL MEANS AND Ns FOR THIS SOURCE (i.e., this combo of factors)
        #
        Lsource = makelist(source-1,Nfactors+1)
        collapsed = _support.collapse(cdata,Lsource,-1,1,1)

        # First, get the list of level-combos for source cells
        prefixcols = range(len(collapsed[0][:-3]))
        outlist = _support.colex(collapsed,prefixcols)
        # Start w/ factor names (A,B,C, or ones input to anova())
        eff = []
        for col in Lsource:
            eff.append(effects[col-1])
        # Add in the mean and N labels for printout
        for item in ['MEAN','STDERR','N']:
            eff.append(item)
        # To the list of level-combos, abut the corresp. means and Ns
        outlist = _support.abut(outlist,
                             map(round4,_support.colex(collapsed,-3)),
                             map(round4,_support.colex(collapsed,-2)),
                             map(round4,_support.colex(collapsed,-1)))
        outlist = [eff] + outlist # add titles to the top of the list
        _support.printcc(outlist)    # print it in customized columns
        print


###
### OUTPUT FINAL RESULTS (ALL SOURCES TOGETHER)
### Note: All 3 types of source-calcs fall through to here
###
    print
    title = [['FACTORS: ','RANDOM'] + effects[:Nfactors]]
    title = title + [['LEVELS:  ']+Nlevels]
    facttypes = ['BETWEEN']*Nfactors
    for i in range(len(Wscols[1:])):
        facttypes[Wscols[i+1]-1] = 'WITHIN'
    title = title + [['TYPE:    ','RANDOM']+facttypes]
    _support.printcc(title)
    print

    title = [['Effect','SS','DF','MS','F','p','sig']] + ['dashes']
    outputlist = title + outputlist
    _support.printcc(outputlist)
    return


def d_full_model(workd, subjslots):
    """
    RESTRICTS NOTHING (i.e., FULL MODEL CALCULATION).  Subtracts D-variable
    cell-mean for each between-subj group and then calculates the SS array.
    """
    workd = subtr_cellmeans(workd,subjslots)
    sserr = multivar_SScalc(workd)
    return sserr


def d_restrict_mean(workd, subjslots):
    """
    RESTRICTS GRAND MEA  Subtracts D-variable cell-mean for each between-
    subj group, and then adds back each D-variable's grand mean.
    """
    # subtract D-variable cell-mean for each (btw-subj) group
    errors = subtr_cellmeans(workd,subjslots)

    # add back in appropriate grand mean from individual scores
    grandDmeans = expand_dims(mean(workd,0),0)
    errors = errors + transpose(grandDmeans) # errors has reversed dims!!
    # SS for mean-restricted model is calculated below.  Note: already put
    # subj as last dim because later code expects this code here to leave
    # workd that way
    sserr = multivar_SScalc(errors)
    return sserr


def d_restrict_source(workd, subjslots, source):
    """
Calculates error for a given model on array workd.  Subjslots is an
array of 1s and 0s corresponding to whether or not the subject is a
member of that between-subjects variable combo.  source is the code
for the type of model to calculate.  source=-1 means no restriction;
source=0 means to restrict workd's grand mean; source>0 means to
restrict the columns of the main data array, DA, specified (in binary)
by the source-value.

Returns: SS array for multivariate F calculation
"""
###
### RESTRICT COLUMNS/AXES SPECIFIED IN source (BINARY)
### (i.e., is the value of source not equal to 0 or -1?)
###
    global D
    if source > 0:
        sourcewithins = (source-1) & Bwithins
        sourcebetweens = (source-1) & Bbetweens
        dindex = Bwonly_sources.index(sourcewithins)
        all_cellmeans = transpose(DM[dindex],[-1]+range(0,len(DM[dindex].shape)-1))
        all_cellns = transpose(DN[dindex],[-1]+range(0,len(DN[dindex].shape)-1))
        hn = hmean(all_cellns, None)

        levels = D[dindex].shape[1]  # GENERAL, 'cause each workd is always 2D
        SSm = zeros((levels,levels),'f') #called RCm=SCm in Lindman,p.317-8
        tworkd = transpose(D[dindex])

        ## Calculate SSw, within-subj variance (Lindman approach)
        RSw = zeros((levels,levels),'f')
        RSinter = zeros((levels,levels),PyObject)
        for i in range(levels):
            for j in range(i,levels):
                RSw[i,j] = RSw[j,i] = sum(tworkd[i]*tworkd[j],axis=0)
                cross = all_cellmeans[i] * all_cellmeans[j]
                multfirst = sum(cross*all_cellns[i],axis=0)
                RSinter[i,j] = RSinter[j,i] = asarray(multfirst)
                SSm[i,j] = SSm[j,i] = (mean(all_cellmeans[i],None) *
                                        mean(all_cellmeans[j],None) *
                                        len(all_cellmeans[i]) *hn)
        #SSw = RSw - RSinter

### HERE BEGINS THE MAXWELL & DELANEY APPROACH TO CALCULATING SS
        Lsource = makelist(sourcebetweens,Nfactors+1)
        #btwsourcecols = (array(map(Bscols.index,Lsource))-1).tolist()
        Bbtwnonsourcedims = ~source & Bbetweens
        Lbtwnonsourcedims = makelist(Bbtwnonsourcedims,Nfactors+1)
        btwnonsourcedims = (array(map(Bscols.index,Lbtwnonsourcedims))-1).tolist()

        # Average Marray over non-source axes
        sourceDMarray = DM[dindex] *1.0
        for dim in btwnonsourcedims: # collapse all non-source dims
            if dim == len(DM[dindex].shape)-1:
                raise ValueError, "Crashing ... shouldn't ever collapse ACROSS variables"
            sourceDMarray = expand_dims(mean(sourceDMarray,dim),dim)

        # Calculate harmonic means for each level in source
        sourceDNarray = apply_over_axes(hmean, DN[dindex],btwnonsourcedims)

        # Calc grand average (ga,axis=0), used for ALL effects
        variableNs = apply_over_axes(sum, sourceDNarray,
                                     range(len(sourceDMarray.shape)-1))
        ga = apply_over_axes(sum, (sourceDMarray*sourceDNarray) / \
                             variableNs,
                             range(len(sourceDMarray.shape)-1))

        # If GRAND interaction, use harmonic mean of ALL cell Ns
        if source == Nallsources-1:
            sourceDNarray = hmean(DN[dindex],
                                  range(len(sourceDMarray.shape)-1))

        # Calc all SUBSOURCES to be subtracted from sourceMarray (M&D p.320)
        sub_effects = ga *1.0   # start with grand mean
        for subsource in range(3,source-2,2):
            # Make a list of the non-subsource axes
            #subsourcebtw = (subsource-1) & Bbetweens
            if (propersubset(subsource-1,source-1) and
                (subsource-1)&Bwithins == (source-1)&Bwithins and
                (subsource-1) != (source-1)&Bwithins):
                sub_effects = (sub_effects +
                                alleffects[alleffsources.index(subsource)])

        # Calc this effect (a(j)'s, b(k)'s, ab(j,k)'s, whatever)
        effect = sourceDMarray - sub_effects

        # Save it so you don't have to calculate it again next time
        alleffects.append(effect)
        alleffsources.append(source)

        # Calc and save sums of squares for this source
        SS = zeros((levels,levels),'f')
        SS = sum((effect**2 *sourceDNarray) *
            multiply.reduce(take(DM[dindex].shape,btwnonsourcedims,0)),
            range(len(sourceDMarray.shape)-1))
        # Save it so you don't have to calculate it again next time
        SSlist.append(SS)
        SSsources.append(source)

        return SS


def multivar_sscalc(workd):
###
### DO SS CALCS ON THE OUTPUT FROM THE SOURCE=0 AND SOURCE=-1 CASES
###
    # this section expects workd to have subj. in LAST axis!!!!!!
    if len(workd.shape) == 1:
        levels = 1
    else:
        levels = workd.shape[0] # works because workd is always 2D

    sserr = zeros((levels,levels),'f')
    for i in range(levels):
        for j in range(i,levels):
            ssval = add.reduce(workd[i]*workd[j])
            sserr[i,j] = ssval
            sserr[j,i] = ssval
    return sserr


def subtr_cellmeans(workd, subjslots):
    """
    Subtract all cell means when within-subjects factors are present ...
    i.e., calculate full-model using a D-variable.
    """
    # Get a list of all dims that are source and between-subj
    sourcedims = makelist(Bbetweens,Nfactors+1)

    # Now, fix this list by mapping the dims from the original source
    # to dims for a between-subjects variable (namely, subjslots)
    transidx = range(len(subjslots.shape))[1:] + [0] # put subj dim at end
    tsubjslots = transpose(subjslots,transidx) # get all Ss for this idx
    tworkd = transpose(workd) # swap subj. and variable dims
    errors = 1.0 * tworkd

    if len(sourcedims) == 0:
        idx = [-1]
        loopcap = [0]
    if len(sourcedims) != 0:
        btwsourcedims = map(Bscols.index,sourcedims)
        idx = [0] * len(btwsourcedims)
        idx[0] = -1 # compensate for pre-increment of 1st slot in incr()

        # Get a list of the maximum values each factor can handle
        loopcap = take(array(Nlevels),sourcedims,0)-1

### WHILE STILL MORE GROUPS, CALCULATE GROUP MEAN FOR EACH D-VAR
    while incr(idx,loopcap) != -1:  # loop through source btw level-combos
        mask = tsubjslots[idx]
        thisgroup = tworkd*mask[newaxis,:]
        groupmns = mean(compress(mask,thisgroup,axis=-1),1)

### THEN SUBTRACT THEM FROM APPROPRIATE SUBJECTS
        errors = errors - multiply.outer(groupmns,mask)
    return errors

def member(factor, source):
    return (1 << factor) & source != 0

def setsize(source):
    size = 0
    for bit in source:
        if bit == 1:
            size = size + 1
    return size

def subset(a, b):
    return (a&b)==a

def propersubset(a, b):
    sub = ((a&b)==a)
    if a==b:
        sub = 0
    return sub

def numlevels(source, Nlevels):
    for i in range(30): # find the biggest i such that 2**i >= source
        if 1<<i >= source:
            break
    levelcount = 1
    for j in range(i): # loop up through each bit
        if subset(1<<j,source):
            levelcount = levelcount * Nlevels[j] - 1
    return levelcount

def numbitson(a):
    numon = 0
    while a>0:
        numon = numon + a%2
        a = a>>1
    return numon

def makebin(sourcelist):
    outbin = 0
    for item in sourcelist:
        outbin = outbin + 2**item
    return outbin

def makelist(source, ncols):
    levellist = []
    for j in range(ncols):
        if subset(1<<j,source):
            levellist.append(j)
    return levellist

def round4(num):
    try:
        return around(num,4)
    except:
        return 'N/A'


def findwithin(data):
    """
Returns a binary vector, 1=within-subject factor, 0=between.  Input
equals the entire data array (i.e., column 1=random factor, last
column = measured values.

"""
    numfact = len(data[0])-2
    withinvec = [0]*numfact
    for col in range(1,numfact+1):
        # get 1 level of this factor
        rows = _support.linexand(data,col,_support.unique(_support.colex(data,1))[0])  
        # if fewer subjects than scores on this factor
        if len(_support.unique(_support.colex(rows,0))) < len(rows):   
            withinvec[col-1] = 1
    return withinvec
