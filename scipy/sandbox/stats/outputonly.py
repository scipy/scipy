def outputfstats(Enum, Eden, dfnum, dfden, f, prob):
    Enum = around(Enum,3)
    Eden = around(Eden,3)
    dfnum = around(Enum,3)
    dfden = around(dfden,3)
    f = around(f,3)
    prob = around(prob,3)
    suffix = ''                       # for *s after the p-value
    if  prob < 0.001:  suffix = '  ***'
    elif prob < 0.01:  suffix = '  **'
    elif prob < 0.05:  suffix = '  *'
    title = [['EF/ER','DF','Mean Square','F-value','prob','']]
    lofl = title+[[Enum, dfnum, around(Enum/float(dfnum),3), f, prob, suffix],
                  [Eden, dfden, around(Eden/float(dfden),3),'','','']]
    _support.printcc(lofl)
    return

def paired(x, y):
    """
Interactively determines the type of data in x and y, and then runs the
appropriated statistic for paired group data.

Returns: appropriate statistic name, value, and probability
"""
    x,y = map(asarray, (x,y))
    samples = ''
    while samples not in ['i','r','I','R','c','C']:
        print '\nIndependent or related samples, or correlation (i,r,c): ',
        samples = raw_input()

    if samples in ['i','I','r','R']:
        print '\nComparing variances ...',
# USE O'BRIEN'S TEST FOR HOMOGENEITY OF VARIANCE, Maxwell & delaney, p.112
        r = obrientransform(x,y)
        f,p = f_oneway(_support.colex(r,0),_support.colex(r,1))
        if p<0.05:
            vartype='unequal, p='+str(around(p,4))
        else:
            vartype='equal'
        print vartype
        if samples in ['i','I']:
            if vartype[0]=='e':
                t,p = ttest_ind(x,y,None,0)
                print '\nIndependent samples t-test:  ', around(t,4),around(p,4)
            else:
                if len(x)>20 or len(y)>20:
                    z,p = ranksums(x,y)
                    print '\nRank Sums test (NONparametric, n>20):  ', around(z,4),around(p,4)
                else:
                    u,p = mannwhitneyu(x,y)
                    print '\nMann-Whitney U-test (NONparametric, ns<20):  ', around(u,4),around(p,4)

        else:  # RELATED SAMPLES
            if vartype[0]=='e':
                t,p = ttest_rel(x,y,0)
                print '\nRelated samples t-test:  ', around(t,4),around(p,4)
            else:
                t,p = ranksums(x,y)
                print '\nWilcoxon T-test (NONparametric):  ', around(t,4),around(p,4)
    else:  # CORRELATION ANALYSIS
        corrtype = ''
        while corrtype not in ['c','C','r','R','d','D']:
            print '\nIs the data Continuous, Ranked, or Dichotomous (c,r,d): ',
            corrtype = raw_input()
        if corrtype in ['c','C']:
            m,b,r,p,see = linregress(x,y)
            print '\nLinear regression for continuous variables ...'
            lol = [['Slope','Intercept','r','Prob','SEestimate'],
                   [around(m,4),around(b,4),around(r,4),around(p,4),around(see,4)]]
            _support.printcc(lol)
        elif corrtype in ['r','R']:
            r,p = spearmanr(x,y)
            print '\nCorrelation for ranked variables ...'
            print "Spearman's r: ",around(r,4),around(p,4)
        else: # DICHOTOMOUS
            r,p = pointbiserialr(x,y)
            print '\nAssuming x contains a dichotomous variable ...'
            print 'Point Biserial r: ',around(r,4),around(p,4)
    print '\n\n'
    return None


def writecc (listoflists, file, writetype='w', extra=2):
    """
Writes a list of lists to a file in columns, customized by the max
size of items within the columns (max size of items in col, +2 characters)
to specified file.  File-overwrite is the default.

Usage:   writecc (listoflists,file,writetype='w',extra=2)
Returns: None
"""
    if not isinstance(listoflists[0], (list, tuple)):
        listoflists = [listoflists]
    outfile = open(file,writetype)
    rowstokill = []
    list2print = copy.deepcopy(listoflists)
    for i in range(len(listoflists)):
        if listoflists[i] == ['\n'] or listoflists[i]=='\n' or listoflists[i]=='dashes':
            rowstokill = rowstokill + [i]
    rowstokill.reverse()
    for row in rowstokill:
        del list2print[row]
    maxsize = [0]*len(list2print[0])
    for col in range(len(list2print[0])):
        items = _support.colex(list2print,col)
        items = map(_support.makestr,items)
        maxsize[col] = max(map(len,items)) + extra
    for row in listoflists:
        if row == ['\n'] or row == '\n':
            outfile.write('\n')
        elif row == ['dashes'] or row == 'dashes':
            dashes = [0]*len(maxsize)
            for j in range(len(maxsize)):
                dashes[j] = '-'*(maxsize[j]-2)
            outfile.write(_support.lineincustcols(dashes,maxsize))
        else:
            outfile.write(_support.lineincustcols(row,maxsize))
        outfile.write('\n')
    outfile.close()
    return None


def outputpairedstats(fname,writemode,name1,n1,m1,se1,min1,max1,name2,n2,m2,se2,min2,max2,statname,stat,prob):
    """
Prints or write to a file stats for two groups, using the name, n,
mean, sterr, min and max for each group, as well as the statistic name,
its value, and the associated p-value.

Usage:   outputpairedstats(fname,writemode,
                           name1,n1,mean1,stderr1,min1,max1,
                           name2,n2,mean2,stderr2,min2,max2,
                           statname,stat,prob)
Returns: None
"""
    suffix = ''                       # for *s after the p-value
    try:
        x = prob.shape
        prob = prob[0]
    except:
        pass
    if  prob < 0.001:  suffix = '  ***'
    elif prob < 0.01:  suffix = '  **'
    elif prob < 0.05:  suffix = '  *'
    title = [['Name','N','Mean','SD','Min','Max']]
    lofl = title+[[name1,n1,round(m1,3),round(math.sqrt(se1),3),min1,max1],
                  [name2,n2,round(m2,3),round(math.sqrt(se2),3),min2,max2]]
    if not isinstance(fname, basestring) or len(fname) == 0:
        print
        print statname
        print
        _support.printcc(lofl)
        print
        try:
            if stat.shape == ():
                stat = stat[0]
            if prob.shape == ():
                prob = prob[0]
        except:
            pass
        print 'Test statistic = ',round(stat,3),'   p = ',round(prob,3),suffix
        print
    else:
        file = open(fname,writemode)
        file.write('\n'+statname+'\n\n')
        file.close()
        writecc(lofl,fname,'a')
        file = open(fname,'a')
        try:
            if stat.shape == ():
                stat = stat[0]
            if prob.shape == ():
                prob = prob[0]
        except:
            pass
        file.write(_support.list2string(['\nTest statistic = ',round(stat,4),'   p = ',round(prob,4),suffix,'\n\n']))
        file.close()
    return None
