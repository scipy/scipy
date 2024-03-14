#include "specfun.h"
#include "special/specfun/specfun.h"

extern "C" {

void specfun_airyzo(int nt, int kf, double *xa, double *xb, double *xc, double *xd) {
    special::specfun::airyzo(nt, kf, xa, xb, xc, xd);
}

void specfun_aswfa(double x, int m, int n, double c, int kd, double cv, double *s1f, double *s1d) {
    special::specfun::aswfa(x, m, n, c, kd, cv, s1f, s1d);
}

void specfun_bernob(int n, double *bn) {
    special::specfun::bernob(n, bn);
}

double _Complex specfun_cerror(double _Complex z)  {
    std::complex<double> res = special::specfun::cerror(*reinterpret_cast<std::complex<double> *>(&z));
    return *reinterpret_cast<double _Complex *>(&res);
}

void specfun_cerzo(int nt, double _Complex *zo) {
    special::specfun::cerzo(nt, reinterpret_cast<std::complex<double> *>(zo));
}

double _Complex specfun_cchg(double a, double b, double _Complex z) {
    std::complex<double> res = special::specfun::cchg(a, b, *reinterpret_cast<std::complex<double> *>(&z));
    return *reinterpret_cast<double _Complex *>(&res);
}

void specfun_cfc(double _Complex z, double _Complex *zf, double _Complex *zd) {
    special::specfun::cfc(*reinterpret_cast<std::complex<double> *>(&z), reinterpret_cast<std::complex<double> *>(zf), reinterpret_cast<std::complex<double> *>(zd));
}

void specfun_cfs(double _Complex z, double _Complex *zf, double _Complex *zd) {
    special::specfun::cfs(*reinterpret_cast<std::complex<double> *>(&z), reinterpret_cast<std::complex<double> *>(zf), reinterpret_cast<std::complex<double> *>(zd));
}

double _Complex specfun_cgama(double _Complex z, int kf) {
    std::complex<double> res = special::specfun::cgama(*reinterpret_cast<std::complex<double> *>(&z), kf);
    return *reinterpret_cast<double _Complex *>(&res);
}

double specfun_chgm(double x, double a, double b) {
    return special::specfun::chgm(x, a, b);
}

double specfun_chgu(double x, double a, double b, int *md, int *isfer) {
    return special::specfun::chgu(x, a, b, md, isfer);
}

void specfun_clpmn(double _Complex z, int m, int n, int ntype, double _Complex *cpm, double _Complex *cpd) {
    special::specfun::clpmn(*reinterpret_cast<std::complex<double> *>(&z), m, n, ntype, reinterpret_cast<std::complex<double> *>(cpm), reinterpret_cast<std::complex<double> *>(cpd));
}

void specfun_clpn(int n, double _Complex z, double _Complex *cpn, double _Complex *cpd) {
    special::specfun::clpn(n, *reinterpret_cast<std::complex<double> *>(&z), reinterpret_cast<std::complex<double> *>(cpn), reinterpret_cast<std::complex<double> *>(cpd));
}

void specfun_clqmn(double _Complex z, int m, int n, double _Complex *cqm, double _Complex *cqd) {
    special::specfun::clqmn(*reinterpret_cast<std::complex<double> *>(&z), m, n, reinterpret_cast<std::complex<double> *>(cqm), reinterpret_cast<std::complex<double> *>(cqd));
}

void specfun_clqn(int n, double _Complex z, double _Complex *cqn, double _Complex *cqd) {
    special::specfun::clqn(n, *reinterpret_cast<std::complex<double> *>(&z), reinterpret_cast<std::complex<double> *>(cqn), reinterpret_cast<std::complex<double> *>(cqd));
}

void specfun_cpbdn(int n, double _Complex z, double _Complex *cpb, double _Complex *cpd) {
    special::specfun::cpbdn(n, *reinterpret_cast<std::complex<double> *>(&z), reinterpret_cast<std::complex<double> *>(cpb), reinterpret_cast<std::complex<double> *>(cpd));
}

double specfun_cva2(int kd, int m, double q) {
    return special::specfun::cva2(kd, m, q);
}

void specfun_cyzo(int nt, int kf, int kc, double _Complex *zo, double _Complex *zv) {
    special::specfun::cyzo(nt, kf, kc, reinterpret_cast<std::complex<double> *>(zo), reinterpret_cast<std::complex<double> *>(zv));
}

double specfun_eix(double x) {
    return special::specfun::eix(x);
}

double specfun_e1xb(double x) {
    return special::specfun::e1xb(x);
}

double _Complex specfun_eixz(double _Complex z) {
    std::complex<double> res = special::specfun::eixz(*reinterpret_cast<std::complex<double> *>(&z));
    return *reinterpret_cast<double _Complex *>(&res);
}

double _Complex specfun_e1z(double _Complex z) {
    std::complex<double> res = special::specfun::e1z(*reinterpret_cast<std::complex<double> *>(&z));
    return *reinterpret_cast<double _Complex *>(&res);
}

void specfun_eulerb(int n, double *en) {
    special::specfun::eulerb(n, en);
}

void specfun_fcoef(int kd, int m, double q, double a, double *fc) {
    special::specfun::fcoef(kd, m, q, a, fc);
}

void specfun_fcszo(int kf, int nt, double _Complex *zo) {
    special::specfun::fcszo(kf, nt, reinterpret_cast<std::complex<double> *>(zo));
}

void specfun_ffk(int ks, double x, double *fr, double *fi, double *fm, double *fa,
         double *gr, double *gi, double *gm, double *ga) {
    special::specfun::ffk(ks, x, fr, fi, fm, fa, gr, gi, gm, ga);
}

double _Complex specfun_hygfz(double a, double b, double c, double _Complex z, int *isfer) {
    std::complex<double> res = special::specfun::hygfz(a, b, c, *reinterpret_cast<std::complex<double> *>(&z), isfer);
    return *reinterpret_cast<double _Complex *>(&res);
}

void specfun_itairy(double x, double *apt, double *bpt, double *ant, double *bnt) {
    special::specfun::itairy(x, apt, bpt, ant, bnt);
}

void specfun_itika(double x, double *ti, double *tk) {
    special::specfun::itika(x, ti, tk);
}

void specfun_itjya(double x, double *tj, double *ty) {
    special::specfun::itjya(x, tj, ty);
}

double specfun_itsh0(double x) {
    return special::specfun::itsh0(x);
}

double specfun_itsl0(double x) {
    return special::specfun::itsl0(x);
}

double specfun_itth0(double x) {
    return special::specfun::itth0(x);
}

void specfun_ittika(double x, double *tti, double *ttk) {
    special::specfun::ittika(x, tti, ttk);
}

void specfun_ittjya(double x, double *ttj, double *tty) {
    special::specfun::ittjya(x, ttj, tty);
}

void specfun_jdzo(int nt, double *zo, int *n, int *m, int *p) {
    special::specfun::jdzo(nt, zo, n, m, p);
}

void specfun_jyzo(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1) {
    special::specfun::jyzo(n, nt, rj0, rj1, ry0, ry1);
}

void specfun_klvna(double x, double *ber, double *bei, double *ger, double *gei,
                   double *der, double *dei, double *her, double *hei) {
    special::specfun::klvna(x, ber, bei, ger, gei, der, dei, her, hei);
}

void specfun_klvnzo(int nt, int kd, double *zo) {
    special::specfun::klvnzo(nt, kd, zo);
}

void specfun_lamn(int n, double x, int *nm, double *bl, double *dl) {
    special::specfun::lamn(n, x, nm, bl, dl);
}

void specfun_lamv(double v, double x, double *vm, double *vl, double *dl)  {
    special::specfun::lamv(v, x, vm, vl, dl);
}

void specfun_lpmn(int m, int n, double x, double *pm, double *pd) { 
    special::specfun::lpmn(m, n, x, pm, pd);
}

double specfun_lpmv(double x, int m, double v) {
    return special::specfun::lpmv(x, m, v);
}

void specfun_lpn(int n, double x, double *pn, double *pd) {
    special::specfun::lpn(n, x, pn, pd);
}

void specfun_lqmn(double x, int m, int n, double *qm, double *qd) {
    special::specfun::lqmn(x, m, n, qm, qd);
}

void specfun_lqnb(int n, double x, double* qn, double* qd) {
    special::specfun::lqnb(n, x, qn, qd);
}

void specfun_mtu0(int kf, int m, double q, double x, double *csf, double *csd) {
    special::specfun::mtu0(kf, m, q, x, csf, csd);
}

void specfun_mtu12(int kf, int kc, int m, double q, double x, double *f1r, double *d1r, double *f2r, double *d2r) {
    special::specfun::mtu12(kf, kc, m, q, x, f1r, d1r, f2r, d2r);
}

void specfun_pbdv(double x, double v, double *dv, double *dp, double *pdf, double *pdd) {
    special::specfun::pbdv(x, v, dv, dp, pdf, pdd);
}

void specfun_pbvv(double x, double v, double *vv, double *vp, double *pvf, double *pvd) {
    special::specfun::pbvv(x, v, vv, vp, pvf, pvd);
}

void specfun_pbwa(double a, double x, double *w1f, double *w1d, double *w2f, double *w2d) {
    special::specfun::pbwa(a, x, w1f, w1d, w2f, w2d);
}

void specfun_rctj(int n, double x, int *nm, double *rj, double *dj) {
    special::specfun::rctj(n, x, nm, rj, dj);
}

void specfun_rcty(int n, double x, int *nm, double *ry, double *dy) {
    special::specfun::rcty(n, x, nm, ry, dy);
}

void specfun_rswfp(int m, int n, double c, double x, double cv, int kf, double *r1f, double *r1d, double *r2f, double *r2d) {
    special::specfun::rswfp(m, n, c, x, cv, kf, r1f, r1d, r2f, r2d);
}

void specfun_rswfo(int m, int n, double c, double x, double cv, int kf, double *r1f, double *r1d, double *r2f, double *r2d) {
    special::specfun::rswfo(m, n, c, x, cv, kf, r1f, r1d, r2f, r2d);
}

void specfun_sdmn(int m, int n, double c, double cv, int kd, double *df) {
    special::specfun::sdmn(m, n, c, cv, kd, df);
}

void specfun_segv(int m, int n, double c, int kd, double *cv, double *eg) {
    special::specfun::segv(m, n, c, kd, cv, eg);
}

}
