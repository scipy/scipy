import numpy.testing as npt
import numpy as np
import nose

from scipy import stats

"""
Test all continuous distributions.

Parameters were chosen for those distributions that pass the
Kolmogorov-Smirnov test.  This provides safe parameters for each
distributions so that we can perform further testing of class methods.

These tests currently check only/mostly for serious errors and exceptions,
not for numerically exact results.


TODO: 
* make functioning test for skew and kurtosis
  still known failures - skip for now


"""

#currently not used
DECIMAL = 5 # specify the precision of the tests  # increased from 0 to 5
DECIMAL_kurt = 0

distcont = [
    ['alpha', (3.5704770516650459,)],
    ['anglit', ()],
    ['arcsine', ()],
    ['beta', (2.3098496451481823, 0.62687954300963677)],
    ['betaprime', (5, 6)],   # avoid unbound error in entropy with (100, 86)],
    ['bradford', (0.29891359763170633,)],
    ['burr', (10.5, 4.3)],    #incorrect mean and var for(0.94839838075366045, 4.3820284068855795)],
    ['cauchy', ()],
    ['chi', (78,)],
    ['chi2', (55,)],
    ['cosine', ()],
    ['dgamma', (1.1023326088288166,)],
    ['dweibull', (2.0685080649914673,)],
    ['erlang', (20,)],    #correction numargs = 1
    ['expon', ()],
    ['exponpow', (2.697119160358469,)],
    ['exponweib', (2.8923945291034436, 1.9505288745913174)],
    ['f', (29, 18)],
    ['fatiguelife', (29,)],   #correction numargs = 1
    ['fisk', (3.0857548622253179,)],
    ['foldcauchy', (4.7164673455831894,)],
    ['foldnorm', (1.9521253373555869,)],
    ['frechet_l', (3.6279911255583239,)],
    ['frechet_r', (1.8928171603534227,)],
    ['gamma', (1.9932305483800778,)],
    ['gausshyper', (13.763771604130699, 3.1189636648681431,
                    2.5145980350183019, 5.1811649903971615)],  #veryslow
    ['genexpon', (9.1325976465418908, 16.231956600590632, 3.2819552690843983)],
    ['genextreme', (-0.1,)],  # sample mean test fails for (3.3184017469423535,)],
    ['gengamma', (4.4162385429431925, 3.1193091679242761)],
    ['genhalflogistic', (0.77274727809929322,)],
    ['genlogistic', (0.41192440799679475,)],
    ['genpareto', (0.1,)],   # use case with finite moments
    ['gilbrat', ()],
    ['gompertz', (0.94743713075105251,)],
    ['gumbel_l', ()],
    ['gumbel_r', ()],
    ['halfcauchy', ()],
    ['halflogistic', ()],
    ['halfnorm', ()],
    ['hypsecant', ()],
    ['invgamma', (2.0668996136993067,)],
    ['invnorm', (0.14546264555347513,)],
    ['invweibull', (10.58,)], # sample mean test fails at(0.58847112119264788,)]
    ['johnsonsb', (4.3172675099141058, 3.1837781130785063)],
    ['johnsonsu', (2.554395574161155, 2.2482281679651965)],
    ['ksone', (1000,)],  #replace 22 by 100 to avoid failing range, ticket 956
    ['kstwobign', ()],
    ['laplace', ()],
    ['levy', ()],
    ['levy_l', ()],
#    ['levy_stable', (0.35667405469844993,
#                     -0.67450531578494011)], #NotImplementedError
    #           rvs not tested
    ['loggamma', (0.41411931826052117,)],
    ['logistic', ()],
    ['loglaplace', (3.2505926592051435,)],
    ['lognorm', (0.95368226960575331,)],
    ['lomax', (1.8771398388773268,)],
    ['maxwell', ()],
    ['mielke', (10.4, 3.6)], # sample mean test fails for (4.6420495492121487, 0.59707419545516938)],
                             # mielke: good results if 2nd parameter >2, weird mean or var below
    ['nakagami', (4.9673794866666237,)],
    ['ncf', (27, 27, 0.41578441799226107)],
    ['nct', (14, 0.24045031331198066)],
    ['ncx2', (21, 1.0560465975116415)],
    ['norm', ()],
    ['pareto', (2.621716532144454,)],
    ['powerlaw', (1.6591133289905851,)],
    ['powerlognorm', (2.1413923530064087, 0.44639540782048337)],
    ['powernorm', (4.4453652254590779,)],
    ['rayleigh', ()],
    ['rdist', (0.9,)],   # feels also slow
#    ['rdist', (3.8266985793976525,)],  #veryslow, especially rvs
    #['rdist', (541.0,)],   # from ticket #758    #veryslow
    ['recipinvgauss', (0.63004267809369119,)],
    ['reciprocal', (0.0062309367010521255, 1.0062309367010522)],
    ['rice', (0.7749725210111873,)],
    ['semicircular', ()],
    ['t', (2.7433514990818093,)],
    ['triang', (0.15785029824528218,)],
    ['truncexpon', (4.6907725456810478,)],
    ['truncnorm', (-1.0978730080013919, 2.7306754109031979)],
    ['tukeylambda', (3.1321477856738267,)],
    ['uniform', ()],
    ['vonmises', (3.9939042581071398,)],
    ['wald', ()],
    ['weibull_max', (2.8687961709100187,)],
    ['weibull_min', (1.7866166930421596,)],
    ['wrapcauchy', (0.031071279018614728,)]]

# for testing only specific functions
##distcont = [
##    ['erlang', (20,)],    #correction numargs = 1
##    ['fatiguelife', (29,)],   #correction numargs = 1
##    ['loggamma', (0.41411931826052117,)]]

# for testing ticket:767
##distcont = [
##    ['genextreme', (3.3184017469423535,)],
##    ['genextreme', (0.01,)],
##    ['genextreme', (0.00001,)],
##    ['genextreme', (0.0,)],
##    ['genextreme', (-0.01,)]
##    ]

##distcont = [['gumbel_l', ()],
##            ['gumbel_r', ()],
##            ['norm', ()]
##            ]

##distcont = [['norm', ()]]

distmissing = ['wald', 'gausshyper', 'genexpon', 'rv_continuous',
    'loglaplace', 'rdist', 'semicircular', 'invweibull', 'ksone',
    'cosine', 'kstwobign', 'truncnorm', 'mielke', 'recipinvgauss', 'levy',
    'johnsonsu', 'levy_l', 'powernorm', 'wrapcauchy',
    'johnsonsb', 'truncexpon', 'rice', 'invnorm', 'invgamma',
    'powerlognorm']

distmiss = [[dist,args] for dist,args in distcont if dist in distmissing]
distslow = ['rdist', 'gausshyper', 'recipinvgauss', 'ksone', 'genexpon',
            'vonmises', 'rice', 'mielke', 'semicircular', 'cosine', 'invweibull',
            'powerlognorm', 'johnsonsu', 'kstwobign']
#distslow are sorted by speed (very slow to slow)

def test_cont_basic():
    # this test skips slow distributions
    for distname, arg in distcont[:]:
        if distname in distslow: continue
        distfn = getattr(stats, distname)
        np.random.seed(765456)
        sn = 1000
        rvs = distfn.rvs(size=sn,*arg)
        sm = rvs.mean()
        sv = rvs.var()
        skurt = stats.kurtosis(rvs)
        sskew = stats.skew(rvs)
        m,v = distfn.stats(*arg)
        
        yield check_sample_meanvar_, distfn, arg, m, v, sm, sv, sn, distname + \
              'sample mean test'
        # the sample skew kurtosis test has known failures, not very good distance measure
        #yield check_sample_skew_kurt, distfn, arg, sskew, skurt, distname
        yield check_moment, distfn, arg, m, v, distname
        yield check_cdf_ppf, distfn, arg, distname
        yield check_sf_isf, distfn, arg, distname
        yield check_pdf, distfn, arg, distname
        if distname in distmissing:
            alpha = 0.01
            yield check_distribution_rvs, dist, args, alpha, rvs


@npt.dec.slow
def test_cont_basic_slow():
    # same as above for slow distributions
    for distname, arg in distcont[:]:
        if distname not in distslow: continue
        distfn = getattr(stats, distname)
        np.random.seed(765456)
        sn = 1000
        rvs = distfn.rvs(size=sn,*arg)
        sm = rvs.mean()
        sv = rvs.var()
        skurt = stats.kurtosis(rvs)
        sskew = stats.skew(rvs)
        m,v = distfn.stats(*arg)
        yield check_sample_meanvar_, distfn, arg, m, v, sm, sv, sn, distname + \
              'sample mean test'
        # the sample skew kurtosis test has known failures, not very good distance measure
        #yield check_sample_skew_kurt, distfn, arg, sskew, skurt, distname
        yield check_moment, distfn, arg, m, v, distname
        yield check_cdf_ppf, distfn, arg, distname
        yield check_sf_isf, distfn, arg, distname
        yield check_pdf, distfn, arg, distname
        #yield check_oth, distfn, arg # is still missing
        if distname in distmissing:
            alpha = 0.01
            yield check_distribution_rvs, dist, args, alpha, rvs




def check_moment(distfn, arg, m, v, msg):
    m1  = distfn.moment(1,*arg)
    m2  = distfn.moment(2,*arg)
    if not np.isinf(m):
        npt.assert_almost_equal(m1, m, decimal=10, err_msg= msg + \
                            ' - 1st moment')
    else:                     # or np.isnan(m1), 
        assert np.isinf(m1), \
               msg + ' - 1st moment -infinite, m1=%s' % str(m1)
        #np.isnan(m1) temporary special treatment for loggamma
    if not np.isinf(v):
        npt.assert_almost_equal(m2-m1*m1, v, decimal=10, err_msg= msg + \
                            ' - 2ndt moment')
    else:                     #or np.isnan(m2), 
        assert np.isinf(m2), \
               msg + ' - 2nd moment -infinite, m2=%s' % str(m2)
        #np.isnan(m2) temporary special treatment for loggamma

def check_sample_meanvar_(distfn, arg, m, v, sm, sv, sn, msg):
    #this did not work, skipped silently by nose
    #check_sample_meanvar, sm, m, msg + 'sample mean test'
    #check_sample_meanvar, sv, v, msg + 'sample var test'
    if not np.isinf(m):
        check_sample_mean(sm, sv, sn, m)
    if not np.isinf(v):
        check_sample_var(sv, sn, v)
##    check_sample_meanvar( sm, m, msg + 'sample mean test')
##    check_sample_meanvar( sv, v, msg + 'sample var test')

def check_sample_mean(sm,v,n, popmean):
    """
from stats.stats.ttest_1samp(a, popmean):
Calculates the t-obtained for the independent samples T-test on ONE group
of scores a, given a population mean.

Returns: t-value, two-tailed prob
"""
##    a = asarray(a)
##    x = np.mean(a)
##    v = np.var(a, ddof=1)
##    n = len(a)
    df = n-1
    svar = ((n-1)*v) / float(df)    #looks redundant
    t = (sm-popmean)/np.sqrt(svar*(1.0/n))
    prob = stats.betai(0.5*df,0.5,df/(df+t*t))

    #return t,prob
    assert prob>0.01, 'mean fail, t,prob = %f, %f, m,sm=%f,%f' % (t,prob,popmean,sm)

def check_sample_var(sv,n, popvar):
    '''
two-sided chisquare test for sample variance equal to hypothesized variance
    '''
    df = n-1
    chi2 = (n-1)*popvar/float(popvar)
    pval = stats.chisqprob(chi2,df)*2
    assert pval>0.01, 'var fail, t,pval = %f, %f, v,sv=%f,%f' % (chi2,pval,popvar,sv)
    

    
def check_sample_skew_kurt(distfn, arg, ss, sk, msg):
    skew,kurt = distfn.stats(moments='sk',*arg)
##    skew = distfn.stats(moment='s',*arg)[()]
##    kurt = distfn.stats(moment='k',*arg)[()]
    check_sample_meanvar( sk, kurt, msg + 'sample kurtosis test')
    check_sample_meanvar( ss, skew, msg + 'sample skew test')

def check_sample_meanvar(sm,m,msg):
    if not np.isinf(m) and not np.isnan(m):
        npt.assert_almost_equal(sm, m, decimal=DECIMAL, err_msg= msg + \
                                ' - finite moment')
##    else:
##        assert abs(sm) > 10000, 'infinite moment, sm = ' + str(sm)

def check_cdf_ppf(distfn,arg,msg):
    npt.assert_almost_equal(distfn.cdf(distfn.ppf([0.001,0.5,0.999], *arg), *arg),
                            [0.001,0.5,0.999], decimal=DECIMAL, err_msg= msg + \
                            ' - cdf-ppf roundtrip')

def check_sf_isf(distfn,arg,msg):
    npt.assert_almost_equal(distfn.sf(distfn.isf([0.1,0.5,0.9], *arg), *arg),
                            [0.1,0.5,0.9], decimal=DECIMAL, err_msg= msg + \
                            ' - sf-isf roundtrip')
    npt.assert_almost_equal(distfn.cdf([0.1,0.9], *arg),
                            1.0-distfn.sf([0.1,0.9], *arg),
                            decimal=DECIMAL, err_msg= msg + \
                            ' - cdf-sf relationship')

def check_pdf(distfn, arg, msg):
    # compares pdf at median with numerical derivative of cdf
    median = distfn.ppf(0.5, *arg)
    eps = 1e-6
    pdfv = distfn.pdf(median, *arg)
    if (pdfv < 1e-4) or (pdfv > 1e4):
        # avoid checking a case where pdf is close to zero or huge (singularity)
        median = median + 0.1
        pdfv = distfn.pdf(median, *arg)
    cdfdiff = (distfn.cdf(median + eps, *arg) -
               distfn.cdf(median - eps, *arg))/eps/2.0
    #replace with better diff and better test (more points),
    #actually, this works pretty well
    npt.assert_almost_equal(pdfv, cdfdiff,
                decimal=DECIMAL, err_msg= msg + ' - cdf-pdf relationship')


def check_distribution_rvs(dist, args, alpha, rvs):
    #test from scipy.stats.tests
    #this version reuses existing random variables
    D,pval = stats.kstest(rvs, dist, args=args, N=1000)
    if (pval < alpha):
        D,pval = stats.kstest(dist,'',args=args, N=1000)
        assert (pval > alpha), "D = " + str(D) + "; pval = " + str(pval) + \
               "; alpha = " + str(alpha) + "\nargs = " + str(args)

if __name__ == "__main__":
    #nose.run(argv=['', __file__])
    nose.runmodule(argv=[__file__,'-s'], exit=False)

