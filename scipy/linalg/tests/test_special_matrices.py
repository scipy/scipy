"""Tests for functions in special_matrices.py."""

from numpy import arange, add, array, eye, copy
from numpy.testing import TestCase, run_module_suite, assert_raises, \
    assert_equal, assert_array_equal, assert_array_almost_equal, \
    assert_allclose

from scipy.misc import comb
from scipy.linalg import toeplitz, hankel, circulant, hadamard, leslie, \
                            companion, tri, triu, tril, kron, block_diag, \
                            hilbert, invhilbert, pascal
from numpy.linalg import cond


def get_mat(n):
    data = arange(n)
    data = add.outer(data,data)
    return data


class TestTri(TestCase):
    def test_basic(self):
        assert_equal(tri(4),array([[1,0,0,0],
                                   [1,1,0,0],
                                   [1,1,1,0],
                                   [1,1,1,1]]))
        assert_equal(tri(4,dtype='f'),array([[1,0,0,0],
                                                [1,1,0,0],
                                                [1,1,1,0],
                                                [1,1,1,1]],'f'))
    def test_diag(self):
        assert_equal(tri(4,k=1),array([[1,1,0,0],
                                       [1,1,1,0],
                                       [1,1,1,1],
                                       [1,1,1,1]]))
        assert_equal(tri(4,k=-1),array([[0,0,0,0],
                                        [1,0,0,0],
                                        [1,1,0,0],
                                        [1,1,1,0]]))
    def test_2d(self):
        assert_equal(tri(4,3),array([[1,0,0],
                                     [1,1,0],
                                     [1,1,1],
                                     [1,1,1]]))
        assert_equal(tri(3,4),array([[1,0,0,0],
                                     [1,1,0,0],
                                     [1,1,1,0]]))
    def test_diag2d(self):
        assert_equal(tri(3,4,k=2),array([[1,1,1,0],
                                         [1,1,1,1],
                                         [1,1,1,1]]))
        assert_equal(tri(4,3,k=-2),array([[0,0,0],
                                          [0,0,0],
                                          [1,0,0],
                                          [1,1,0]]))

class TestTril(TestCase):
    def test_basic(self):
        a = (100*get_mat(5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[k,l] = 0
        assert_equal(tril(a),b)

    def test_diag(self):
        a = (100*get_mat(5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[k,l] = 0
        assert_equal(tril(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[k,l] = 0
        assert_equal(tril(a,k=-2),b)


class TestTriu(TestCase):
    def test_basic(self):
        a = (100*get_mat(5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[l,k] = 0
        assert_equal(triu(a),b)

    def test_diag(self):
        a = (100*get_mat(5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[l,k] = 0
        assert_equal(triu(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[l,k] = 0
        assert_equal(triu(a,k=-2),b)


class TestToeplitz(TestCase):

    def test_basic(self):
        y = toeplitz([1,2,3])
        assert_array_equal(y,[[1,2,3],[2,1,2],[3,2,1]])
        y = toeplitz([1,2,3],[1,4,5])
        assert_array_equal(y,[[1,4,5],[2,1,4],[3,2,1]])

    def test_complex_01(self):
        data = (1.0 + arange(3.0)) * (1.0 + 1.0j)
        x = copy(data)
        t = toeplitz(x)
        # Calling toeplitz should not change x.
        assert_array_equal(x, data)
        # According to the docstring, x should be the first column of t.
        col0 = t[:,0]
        assert_array_equal(col0, data)
        assert_array_equal(t[0,1:], data[1:].conj())

    def test_scalar_00(self):
        """Scalar arguments still produce a 2D array."""
        t = toeplitz(10)
        assert_array_equal(t, [[10]])
        t = toeplitz(10, 20)
        assert_array_equal(t, [[10]])

    def test_scalar_01(self):
        c = array([1,2,3])
        t = toeplitz(c, 1)
        assert_array_equal(t, [[1],[2],[3]])

    def test_scalar_02(self):
        c = array([1,2,3])
        t = toeplitz(c, array(1))
        assert_array_equal(t, [[1],[2],[3]])

    def test_scalar_03(self):
        c = array([1,2,3])
        t = toeplitz(c, array([1]))
        assert_array_equal(t, [[1],[2],[3]])

    def test_scalar_04(self):
        r = array([10,2,3])
        t = toeplitz(1, r)
        assert_array_equal(t, [[1,2,3]])


class TestHankel(TestCase):
    def test_basic(self):
        y = hankel([1,2,3])
        assert_array_equal(y, [[1,2,3], [2,3,0], [3,0,0]])
        y = hankel([1,2,3], [3,4,5])
        assert_array_equal(y, [[1,2,3], [2,3,4], [3,4,5]])


class TestCirculant(TestCase):
    def test_basic(self):
        y = circulant([1,2,3])
        assert_array_equal(y, [[1,3,2], [2,1,3], [3,2,1]])


class TestHadamard(TestCase):

    def test_basic(self):

        y = hadamard(1)
        assert_array_equal(y, [[1]])

        y = hadamard(2, dtype=float)
        assert_array_equal(y, [[1.0, 1.0], [1.0, -1.0]])

        y = hadamard(4)
        assert_array_equal(y, [[1,1,1,1], [1,-1,1,-1], [1,1,-1,-1], [1,-1,-1,1]])

        assert_raises(ValueError, hadamard, 0)
        assert_raises(ValueError, hadamard, 5)


class TestLeslie(TestCase):

    def test_bad_shapes(self):
        assert_raises(ValueError, leslie, [[1,1],[2,2]], [3,4,5])
        assert_raises(ValueError, leslie, [3,4,5], [[1,1],[2,2]])
        assert_raises(ValueError, leslie, [1,2], [1,2])
        assert_raises(ValueError, leslie, [1], [])

    def test_basic(self):
        a = leslie([1, 2, 3], [0.25, 0.5])
        expected = array([
            [1.0,  2.0, 3.0],
            [0.25, 0.0, 0.0],
            [0.0,  0.5, 0.0]])
        assert_array_equal(a, expected)


class TestCompanion(TestCase):

    def test_bad_shapes(self):
        assert_raises(ValueError, companion, [[1,1],[2,2]])
        assert_raises(ValueError, companion, [0,4,5])
        assert_raises(ValueError, companion, [1])
        assert_raises(ValueError, companion, [])

    def test_basic(self):
        c = companion([1, 2, 3])
        expected = array([
            [-2.0, -3.0],
            [ 1.0,  0.0]])
        assert_array_equal(c, expected)

        c = companion([2.0, 5.0, -10.0])
        expected = array([
            [-2.5, 5.0],
            [ 1.0, 0.0]])
        assert_array_equal(c, expected)


class TestBlockDiag:
    def test_basic(self):
        x = block_diag(eye(2), [[1,2], [3,4], [5,6]], [[1, 2, 3]])
        assert_array_equal(x, [[1, 0, 0, 0, 0, 0, 0],
                               [0, 1, 0, 0, 0, 0, 0],
                               [0, 0, 1, 2, 0, 0, 0],
                               [0, 0, 3, 4, 0, 0, 0],
                               [0, 0, 5, 6, 0, 0, 0],
                               [0, 0, 0, 0, 1, 2, 3]])

    def test_dtype(self):
        x = block_diag([[1.5]])
        assert_equal(x.dtype, float)

        x = block_diag([[True]])
        assert_equal(x.dtype, bool)

    def test_scalar_and_1d_args(self):
        a = block_diag(1)
        assert_equal(a.shape, (1,1))
        assert_array_equal(a, [[1]])

        a = block_diag([2,3], 4)
        assert_array_equal(a, [[2, 3, 0], [0, 0, 4]])

    def test_bad_arg(self):
        assert_raises(ValueError, block_diag, [[[1]]])

    def test_no_args(self):
        a = block_diag()
        assert_equal(a.ndim, 2)
        assert_equal(a.nbytes, 0)


class TestKron:

    def test_basic(self):

        a = kron(array([[1, 2], [3, 4]]), array([[1, 1, 1]]))
        assert_array_equal(a, array([[1, 1, 1, 2, 2, 2],
                                     [3, 3, 3, 4, 4, 4]]))

        m1 = array([[1, 2], [3, 4]])
        m2 = array([[10], [11]])
        a = kron(m1, m2)
        expected = array([[ 10, 20 ],
                          [ 11, 22 ],
                          [ 30, 40 ],
                          [ 33, 44 ]])
        assert_array_equal(a, expected)


class TestHilbert(TestCase):

    def test_basic(self):
        h3 = array([[1.0,  1/2., 1/3.],
                    [1/2., 1/3., 1/4.],
                    [1/3., 1/4., 1/5.]])
        assert_array_almost_equal(hilbert(3), h3)

        assert_array_equal(hilbert(1), [[1.0]])

        h0 = hilbert(0)
        assert_equal(h0.shape, (0,0))


class TestInvHilbert(TestCase):

    def test_basic(self):
        invh1 = array([[1]])
        assert_array_equal(invhilbert(1, exact=True), invh1)
        assert_array_equal(invhilbert(1), invh1)

        invh2 = array([[ 4, -6],
                       [-6, 12]])
        assert_array_equal(invhilbert(2, exact=True), invh2)
        assert_array_almost_equal(invhilbert(2), invh2)

        invh3 = array([[  9,  -36,  30],
                       [-36,  192, -180],
                        [30, -180,  180]])
        assert_array_equal(invhilbert(3, exact=True), invh3)
        assert_array_almost_equal(invhilbert(3), invh3)

        invh4 = array([[  16,  -120,   240,  -140],
                       [-120,  1200, -2700,  1680],
                       [ 240, -2700,  6480, -4200],
                       [-140,  1680, -4200,  2800]])
        assert_array_equal(invhilbert(4, exact=True), invh4)
        assert_array_almost_equal(invhilbert(4), invh4)

        invh5 = array([[   25,   -300,    1050,  -1400,     630],
                       [ -300,   4800,  -18900,   26880, -12600],
                       [ 1050, -18900,   79380, -117600,  56700],
                       [-1400,  26880, -117600,  179200, -88200],
                       [  630, -12600,   56700,  -88200,  44100]])
        assert_array_equal(invhilbert(5, exact=True), invh5)
        assert_array_almost_equal(invhilbert(5), invh5)

        invh17 = array([
            [289, -41616, 1976760, -46124400, 629598060, -5540462928,
             33374693352, -143034400080, 446982500250, -1033026222800,
             1774926873720, -2258997839280, 2099709530100, -1384423866000,
             613101997800, -163493866080, 19835652870],
            [-41616, 7990272, -426980160, 10627061760, -151103534400, 1367702848512,
             -8410422724704, 36616806420480, -115857864064800, 270465047424000,
             -468580694662080, 600545887119360, -561522320049600, 372133135180800,
             -165537539406000, 44316454993920, -5395297580640],
            [1976760, -426980160, 24337869120, -630981792000, 9228108708000,
             -85267724461920, 532660105897920, -2348052711713280, 7504429831470000,
             -17664748409880000, 30818191841236800, -39732544853164800,
             37341234283298400, -24857330514030000, 11100752642520000,
             -2982128117299200, 364182586693200],
            [-46124400, 10627061760, -630981792000, 16826181120000,
             -251209625940000, 2358021022156800, -14914482965141760,
             66409571644416000, -214015221119700000, 507295338950400000,
             -890303319857952000, 1153715376477081600, -1089119333262870000,
             727848632044800000, -326170262829600000, 87894302404608000,
             -10763618673376800],
            [629598060, -151103534400, 9228108708000,
             -251209625940000, 3810012660090000, -36210360321495360,
             231343968720664800, -1038687206500944000, 3370739732635275000,
             -8037460526495400000, 14178080368737885600, -18454939322943942000,
             17489975175339030000, -11728977435138600000, 5272370630081100000,
             -1424711708039692800, 174908803442373000],
            [-5540462928, 1367702848512, -85267724461920, 2358021022156800,
             -36210360321495360, 347619459086355456, -2239409617216035264,
             10124803292907663360, -33052510749726468000, 79217210949138662400,
             -140362995650505067440, 183420385176741672960, -174433352415381259200,
             117339159519533952000, -52892422160973595200, 14328529177999196160,
             -1763080738699119840],
            [33374693352, -8410422724704, 532660105897920,
             -14914482965141760, 231343968720664800, -2239409617216035264,
             14527452132196331328, -66072377044391477760, 216799987176909536400,
             -521925895055522958000, 928414062734059661760, -1217424500995626443520,
             1161358898976091015200, -783401860847777371200, 354015418167362952000,
             -96120549902411274240, 11851820521255194480],
            [-143034400080, 36616806420480, -2348052711713280, 66409571644416000,
             -1038687206500944000, 10124803292907663360, -66072377044391477760,
             302045152202932469760, -995510145200094810000, 2405996923185123840000,
             -4294704507885446054400, 5649058909023744614400,
             -5403874060541811254400, 3654352703663101440000,
             -1655137020003255360000, 450325202737117593600, -55630994283442749600],
            [446982500250, -115857864064800, 7504429831470000, -214015221119700000,
             3370739732635275000, -33052510749726468000, 216799987176909536400,
             -995510145200094810000, 3293967392206196062500,
             -7988661659013106500000, 14303908928401362270000,
             -18866974090684772052000, 18093328327706957325000,
             -12263364009096700500000, 5565847995255512250000,
             -1517208935002984080000, 187754605706619279900],
            [-1033026222800, 270465047424000, -17664748409880000,
             507295338950400000, -8037460526495400000, 79217210949138662400,
             -521925895055522958000, 2405996923185123840000,
             -7988661659013106500000, 19434404971634224000000,
             -34894474126569249192000, 46141453390504792320000,
             -44349976506971935800000, 30121928988527376000000,
             -13697025107665828500000, 3740200989399948902400,
             -463591619028689580000],
            [1774926873720, -468580694662080,
             30818191841236800, -890303319857952000, 14178080368737885600,
             -140362995650505067440, 928414062734059661760, -4294704507885446054400,
             14303908928401362270000, -34894474126569249192000,
             62810053427824648545600, -83243376594051600326400,
             80177044485212743068000, -54558343880470209780000,
             24851882355348879230400, -6797096028813368678400, 843736746632215035600],
            [-2258997839280, 600545887119360, -39732544853164800,
             1153715376477081600, -18454939322943942000, 183420385176741672960,
             -1217424500995626443520, 5649058909023744614400,
             -18866974090684772052000, 46141453390504792320000,
             -83243376594051600326400, 110552468520163390156800,
             -106681852579497947388000, 72720410752415168870400,
             -33177973900974346080000, 9087761081682520473600,
             -1129631016152221783200],
            [2099709530100, -561522320049600, 37341234283298400,
             -1089119333262870000, 17489975175339030000, -174433352415381259200,
             1161358898976091015200, -5403874060541811254400,
             18093328327706957325000, -44349976506971935800000,
             80177044485212743068000, -106681852579497947388000,
             103125790826848015808400, -70409051543137015800000,
             32171029219823375700000, -8824053728865840192000,
             1098252376814660067000],
            [-1384423866000, 372133135180800,
             -24857330514030000, 727848632044800000, -11728977435138600000,
             117339159519533952000, -783401860847777371200, 3654352703663101440000,
             -12263364009096700500000, 30121928988527376000000,
             -54558343880470209780000, 72720410752415168870400,
             -70409051543137015800000, 48142941226076592000000,
             -22027500987368499000000, 6049545098753157120000,
             -753830033789944188000],
            [613101997800, -165537539406000,
             11100752642520000, -326170262829600000, 5272370630081100000,
             -52892422160973595200, 354015418167362952000, -1655137020003255360000,
             5565847995255512250000, -13697025107665828500000,
             24851882355348879230400, -33177973900974346080000,
             32171029219823375700000, -22027500987368499000000,
             10091416708498869000000, -2774765838662800128000, 346146444087219270000],
            [-163493866080, 44316454993920, -2982128117299200, 87894302404608000,
             -1424711708039692800, 14328529177999196160, -96120549902411274240,
             450325202737117593600, -1517208935002984080000, 3740200989399948902400,
             -6797096028813368678400, 9087761081682520473600,
             -8824053728865840192000, 6049545098753157120000,
             -2774765838662800128000, 763806510427609497600, -95382575704033754400],
            [19835652870, -5395297580640, 364182586693200, -10763618673376800,
             174908803442373000, -1763080738699119840, 11851820521255194480,
             -55630994283442749600, 187754605706619279900, -463591619028689580000,
             843736746632215035600, -1129631016152221783200, 1098252376814660067000,
             -753830033789944188000, 346146444087219270000, -95382575704033754400,
             11922821963004219300]
            ])
        assert_array_equal(invhilbert(17, exact=True), invh17)
        assert_allclose(invhilbert(17), invh17.astype(float), rtol=1e-12)

    def test_inverse(self):
        for n in xrange(1, 10):
            a = hilbert(n)
            b = invhilbert(n)
            # The Hilbert matrix is increasingly badly conditioned,
            # so take that into account in the test
            c = cond(a)
            assert_allclose(a.dot(b), eye(n), atol=1e-15*c, rtol=1e-15*c)


class TestPascal(TestCase):

    cases = [
        (1, array([[1]]), array([[1]])),
        (2, array([[1, 1],
                   [1, 2]]),
            array([[1, 0],
                   [1, 1]])),
        (3, array([[1, 1, 1],
                   [1, 2, 3],
                   [1, 3, 6]]),
            array([[1, 0, 0],
                   [1, 1, 0],
                   [1, 2, 1]])),
        (4, array([[1,  1,  1,  1],
                   [1,  2,  3,  4],
                   [1,  3,  6, 10],
                   [1,  4, 10, 20]]),
            array([[1, 0, 0, 0],
                   [1, 1, 0, 0],
                   [1, 2, 1, 0],
                   [1, 3, 3, 1]])),
    ]

    def check_case(self, n, sym, low):
        assert_array_equal(pascal(n), sym)
        assert_array_equal(pascal(n, kind='lower'), low)
        assert_array_equal(pascal(n, kind='upper'), low.T)
        assert_array_almost_equal(pascal(n, exact=False), sym)
        assert_array_almost_equal(pascal(n, exact=False, kind='lower'), low)
        assert_array_almost_equal(pascal(n, exact=False, kind='upper'), low.T)

    def test_cases(self):
        for n, sym, low in self.cases:
            self.check_case(n, sym, low)

    def test_big(self):
        p = pascal(50)
        assert_equal(p[-1, -1], comb(98, 49, exact=True))


if __name__ == "__main__":
    run_module_suite()
