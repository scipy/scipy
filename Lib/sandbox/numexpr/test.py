from numpy import *
from numexpr import evaluate

tests = [
('MISC', ['b*c+d*e', 
          '2*a+3*b', 
          '2*a + (cos(3)+5)*sinh(cos(b))',
          '2*a + arctan2(a, b)',
          'where(0.1*a > arctan2(a, b), 2*a, arctan2(a,b))',
          'where(a, 2, b)',
          'where(a-10, a, 2)',
          'cos(1+1)',
          '1+1',
          '1',
          'cos(a2)'])]
cmptests = []
for op in ['<', '<=', '==', '>=', '>', '!=']:
    cmptests.append("a/2+5 %s b" % op)
tests.append(('COMPARISONS', cmptests))
func1tests = []
for func in ['sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh']:
    func1tests.append("a + %s(b+c)" % func)
tests.append(('1-ARG FUNCS', func1tests))
func2tests = []
for func in ['arctan2', 'fmod']:
    func2tests.append("a + %s(b+c, d+1)" % func)
    func2tests.append("a + %s(b+c, 1)" % func)
    func2tests.append("a + %s(1, d+1)" % func)
tests.append(('2-ARG FUNCS', func2tests))


array_size = 1e2
a = arange(array_size)
a2 = zeros([array_size, array_size])
b = arange(array_size)
c = arange(array_size)
d = arange(array_size)
e = arange(array_size)

def check(expr):
    try:
        npval = eval(expr)
        neval = evaluate(expr)
        if not (shape(npval) == shape(neval) and alltrue(ravel(npval) == ravel(neval))):
            print ">>> %s: FAILED" % expr
            print "  GOT      :", neval[:3], "..."
            print "  EXPECTED :", npval[:3], "..."
            return False
        else:
            print "    %s" % expr#, evaluate(expr)[:3], "..."
            return True
    except StandardError, err:
        print ">>> %s: ERROR(%s)" % (expr, err)
        return False
        
def test():
    all_pass = True
    for category, expressions in tests:
        print "%s:" % category
        for x in expressions:
            all_pass &= check(x)
    print "*******************"
    if not all_pass:
        print "THERE WERE FAILURES"
    else:
        print "ALL TESTS PASSED"
    print "*******************"







if __name__ == '__main__':
    test()
