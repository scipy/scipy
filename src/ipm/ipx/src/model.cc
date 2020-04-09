// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "model.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "utils.h"

namespace ipx {

void Model::Load(const Control& control, Int num_constr, Int num_var,
                 const Int* Ap, const Int* Ai, const double* Ax,
                 const double* rhs, const char* constr_type, const double* obj,
                 const double* lbuser, const double* ubuser, Info* info) {
    clear();
    CopyInput(num_constr, num_var, Ap, Ai, Ax, rhs, constr_type, obj, lbuser,
              ubuser, info);
    if (info->errflag)
        return;
    control.Log()
        << "Input\n"
        << Textline("Number of variables:") << num_var_ << '\n'
        << Textline("Number of free variables:") << num_free_var_ << '\n'
        << Textline("Number of constraints:") << num_constr_ << '\n'
        << Textline("Number of equality constraints:") << num_eqconstr_ << '\n'
        << Textline("Number of matrix entries:") << num_entries_ << '\n';
    PrintCoefficientRange(control);
    ScaleModel(control);

    // Make an automatic decision for dualization if not specified by user.
    Int dualize = control.dualize();
    if (dualize < 0)
        dualize = num_constr > 2*num_var;
    if (dualize)
        LoadDual();
    else
        LoadPrimal();

    A_.clear();
    AIt_ = Transpose(AI_);
    assert(AI_.begin(num_cols_ + num_rows_) == AIt_.begin(num_rows_));
    FindDenseColumns();
    norm_c_ = Infnorm(c_);
    norm_bounds_ = Infnorm(b_);
    for (double x : lb_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
    for (double x : ub_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
    PrintPreprocessingLog(control);
    WriteInfo(info);
}

void Model::clear() {
    // clear computational form model
    dualized_ = false;
    num_rows_ = 0;
    num_cols_ = 0;
    num_dense_cols_ = 0;
    nz_dense_ = 0;
    AI_.clear();
    AIt_.clear();
    b_.resize(0);
    c_.resize(0);
    lb_.resize(0);
    ub_.resize(0);
    norm_bounds_ = 0.0;
    norm_c_ = 0.0;

    // clear user model
    num_constr_ = 0;
    num_eqconstr_ = 0;
    num_var_ = 0;
    num_free_var_ = 0;
    num_entries_ = 0;
    boxed_vars_.clear();
    constr_type_.clear();
    norm_obj_ = 0.0;
    norm_rhs_ = 0.0;
    scaled_obj_.resize(0);
    scaled_rhs_.resize(0);
    scaled_lbuser_.resize(0);
    scaled_ubuser_.resize(0);
    A_.clear();

    flipped_vars_.clear();
    colscale_.resize(0);
    rowscale_.resize(0);
}

void Model::PresolveStartingPoint(const double* x_user,
                                  const double* slack_user,
                                  const double* y_user,
                                  const double* z_user,
                                  Vector& x_solver,
                                  Vector& y_solver,
                                  Vector& z_solver) const {
    const Int m = rows();
    const Int n = cols();
    assert((int)x_solver.size() == n+m);
    assert((int)y_solver.size() == m);
    assert((int)z_solver.size() == n+m);

    Vector x_temp(num_var_);
    Vector slack_temp(num_constr_);
    Vector y_temp(num_constr_);
    Vector z_temp(num_var_);
    if (x_user)
        std::copy_n(x_user, num_var_, std::begin(x_temp));
    if (slack_user)
        std::copy_n(slack_user, num_constr_, std::begin(slack_temp));
    if (y_user)
        std::copy_n(y_user, num_constr_, std::begin(y_temp));
    if (z_user)
        std::copy_n(z_user, num_var_, std::begin(z_temp));
    ScaleBasicSolution(x_temp, slack_temp, y_temp, z_temp);
    DualizeBasicSolution(x_temp, slack_temp, y_temp, z_temp,
                         x_solver, y_solver, z_solver);
    (void)(n);
    (void)(m);
}

void Model::PostsolveInteriorSolution(const Vector& x_solver,
                                      const Vector& xl_solver,
                                      const Vector& xu_solver,
                                      const Vector& y_solver,
                                      const Vector& zl_solver,
                                      const Vector& zu_solver,
                                      double* x_user,
                                      double* xl_user, double* xu_user,
                                      double* slack_user,
                                      double* y_user,
                                      double* zl_user, double* zu_user) const {
    const Int m = rows();
    const Int n = cols();
    assert((int)x_solver.size() == n+m);
    assert((int)xl_solver.size() == n+m);
    assert((int)xu_solver.size() == n+m);
    assert((int)y_solver.size() == m);
    assert((int)zl_solver.size() == n+m);
    assert((int)zu_solver.size() == n+m);

    Vector x_temp(num_var_);
    Vector xl_temp(num_var_);
    Vector xu_temp(num_var_);
    Vector slack_temp(num_constr_);
    Vector y_temp(num_constr_);
    Vector zl_temp(num_var_);
    Vector zu_temp(num_var_);
    DualizeBackInteriorSolution(x_solver, xl_solver, xu_solver, y_solver,
                                zl_solver, zu_solver, x_temp, xl_temp, xu_temp,
                                slack_temp, y_temp, zl_temp, zu_temp);
    ScaleBackInteriorSolution(x_temp, xl_temp, xu_temp, slack_temp, y_temp,
                              zl_temp, zu_temp);
    if (x_user)
        std::copy(std::begin(x_temp), std::end(x_temp), x_user);
    if (xl_user)
        std::copy(std::begin(xl_temp), std::end(xl_temp), xl_user);
    if (xu_user)
        std::copy(std::begin(xu_temp), std::end(xu_temp), xu_user);
    if (slack_user)
        std::copy(std::begin(slack_temp), std::end(slack_temp), slack_user);
    if (y_user)
        std::copy(std::begin(y_temp), std::end(y_temp), y_user);
    if (zl_user)
        std::copy(std::begin(zl_temp), std::end(zl_temp), zl_user);
    if (zu_user)
        std::copy(std::begin(zu_temp), std::end(zu_temp), zu_user);
    (void)(n);
    (void)(m);
}

void Model::EvaluateInteriorSolution(const Vector& x_solver,
                                     const Vector& xl_solver,
                                     const Vector& xu_solver,
                                     const Vector& y_solver,
                                     const Vector& zl_solver,
                                     const Vector& zu_solver,
                                     Info* info) const {
    const Int m = rows();
    const Int n = cols();
    assert((int)x_solver.size() == n+m);
    assert((int)xl_solver.size() == n+m);
    assert((int)xu_solver.size() == n+m);
    assert((int)y_solver.size() == m);
    assert((int)zl_solver.size() == n+m);
    assert((int)zu_solver.size() == n+m);

    // Build solution to scaled user model.
    Vector x(num_var_);
    Vector xl(num_var_);
    Vector xu(num_var_);
    Vector slack(num_constr_);
    Vector y(num_constr_);
    Vector zl(num_var_);
    Vector zu(num_var_);
    DualizeBackInteriorSolution(x_solver, xl_solver, xu_solver, y_solver,
                                zl_solver, zu_solver, x, xl, xu, slack, y, zl,
                                zu);

    // Build residuals for scaled model.
    // rl = lb-x+xl
    Vector rl(num_var_);
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(scaled_lbuser_[j]))
            rl[j] = scaled_lbuser_[j] - x[j] + xl[j];
        else
            assert(xl[j] == INFINITY);
    }
    // ru = ub-x-xu
    Vector ru(num_var_);
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(scaled_ubuser_[j]))
            ru[j] = scaled_ubuser_[j] - x[j] - xu[j];
        else
            assert(xu[j] == INFINITY);
    }
    // rb = rhs-slack-A*x
    Vector rb = scaled_rhs_ - slack;
    assert((int)scaled_rhs_.size() == num_constr_);
    assert((int)rb.size() == num_constr_);
    MultiplyWithScaledMatrix(x, -1.0, rb, 'N');
    // rc = obj-zl+zu-A'y
    Vector rc = scaled_obj_ - zl + zu;
    MultiplyWithScaledMatrix(y, -1.0, rc, 'T');
    
    ScaleBackResiduals(rb, rc, rl, ru);
    double presidual = Infnorm(rb);
    presidual = std::max(presidual, Infnorm(rl));
    presidual = std::max(presidual, Infnorm(ru));
    double dresidual = Infnorm(rc);

    double pobjective = Dot(scaled_obj_, x);
    double dobjective = Dot(scaled_rhs_, y);
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(scaled_lbuser_[j]))
            dobjective += scaled_lbuser_[j] * zl[j];
        if (std::isfinite(scaled_ubuser_[j]))
            dobjective -= scaled_ubuser_[j] * zu[j];
    }
    assert(std::isfinite(dobjective));
    double objective_gap = (pobjective-dobjective) /
        (1.0 + 0.5*std::abs(pobjective+dobjective));

    double complementarity = 0.0;
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(scaled_lbuser_[j]))
            complementarity += xl[j]*zl[j];
        if (std::isfinite(scaled_ubuser_[j]))
            complementarity += xu[j]*zu[j];
    }
    for (Int i = 0; i < num_constr_; i++)
        complementarity -= y[i]*slack[i];

    // For computing the norms of the user variables, we have to scale back.
    ScaleBackInteriorSolution(x, xl, xu, slack, y, zl, zu);

    info->abs_presidual = presidual;
    info->abs_dresidual = dresidual;
    info->rel_presidual = presidual/(1.0+norm_rhs_);
    info->rel_dresidual = dresidual/(1.0+norm_obj_);
    info->pobjval = pobjective;
    info->dobjval = dobjective;
    info->rel_objgap = objective_gap;
    info->complementarity = complementarity;
    info->normx = Infnorm(x);
    info->normy = Infnorm(y);
    info->normz = std::max(Infnorm(zl), Infnorm(zu));
    (void)(n);
    (void)(m);
}

void Model::PostsolveBasicSolution(const Vector& x_solver,
                                   const Vector& y_solver,
                                   const Vector& z_solver,
                                   const std::vector<Int>& basic_status_solver,
                                   double* x_user, double* slack_user,
                                   double* y_user, double* z_user) const {
    const Int m = rows();
    const Int n = cols();
    assert((int)x_solver.size() == n+m);
    assert((int)y_solver.size() == m);
    assert((int)z_solver.size() == n+m);
    assert((int)basic_status_solver.size() == n+m);

    Vector x_temp(num_var_);
    Vector slack_temp(num_constr_);
    Vector y_temp(num_constr_);
    Vector z_temp(num_var_);
    std::vector<Int> cbasis_temp(num_constr_);
    std::vector<Int> vbasis_temp(num_var_);
    DualizeBackBasicSolution(x_solver, y_solver, z_solver, x_temp, slack_temp,
                             y_temp, z_temp);
    DualizeBackBasis(basic_status_solver, cbasis_temp, vbasis_temp);
    CorrectScaledBasicSolution(x_temp, slack_temp, y_temp, z_temp, cbasis_temp,
                               vbasis_temp);
    ScaleBackBasicSolution(x_temp, slack_temp, y_temp, z_temp);
    if (x_user)
        std::copy(std::begin(x_temp), std::end(x_temp), x_user);
    if (slack_user)
        std::copy(std::begin(slack_temp), std::end(slack_temp), slack_user);
    if (y_user)
        std::copy(std::begin(y_temp), std::end(y_temp), y_user);
    if (z_user)
        std::copy(std::begin(z_temp), std::end(z_temp), z_user);
    (void)(n);
    (void)(m);
}

void Model::EvaluateBasicSolution(const Vector& x_solver,
                                  const Vector& y_solver,
                                  const Vector& z_solver,
                                  const std::vector<Int>& basic_status_solver,
                                  Info* info) const {
    const Int m = rows();
    const Int n = cols();
    assert((int)x_solver.size() == n+m);
    assert((int)y_solver.size() == m);
    assert((int)z_solver.size() == n+m);
    assert((int)basic_status_solver.size() == n+m);

    // Build basic solution to scaled user model.
    Vector x(num_var_);
    Vector slack(num_constr_);
    Vector y(num_constr_);
    Vector z(num_var_);
    std::vector<Int> cbasis(num_constr_);
    std::vector<Int> vbasis(num_var_);
    DualizeBackBasicSolution(x_solver, y_solver, z_solver, x, slack, y, z);
    DualizeBackBasis(basic_status_solver, cbasis, vbasis);
    CorrectScaledBasicSolution(x, slack, y, z, cbasis, vbasis);
    double pobj = Dot(scaled_obj_, x);

    // Build infeasibilities in scaled user model.
    Vector xinfeas(num_var_);
    Vector sinfeas(num_constr_);
    Vector yinfeas(num_constr_);
    Vector zinfeas(num_var_);
    for (Int j = 0; j < num_var_; j++) {
        if (x[j] < scaled_lbuser_[j])
            xinfeas[j] = x[j]-scaled_lbuser_[j];
        if (x[j] > scaled_ubuser_[j])
            xinfeas[j] = x[j]-scaled_ubuser_[j];
        if (vbasis[j] != IPX_nonbasic_lb && z[j] > 0.0)
            zinfeas[j] = z[j];
        if (vbasis[j] != IPX_nonbasic_ub && z[j] < 0.0)
            zinfeas[j] = z[j];
    }
    for (Int i = 0; i < num_constr_; i++) {
        if (constr_type_[i] == '<') {
            if (slack[i] < 0.0)
                sinfeas[i] = slack[i];
            if (y[i] > 0.0)
                yinfeas[i] = y[i];
        }
        if (constr_type_[i] == '>') {
            if (slack[i] > 0.0)
                sinfeas[i] = slack[i];
            if (y[i] < 0.0)
                yinfeas[i] = y[i];
        }
    }

    // Scale back basic solution and infeasibilities.
    ScaleBackBasicSolution(x, slack, y, z);
    ScaleBackBasicSolution(xinfeas, sinfeas, yinfeas, zinfeas);

    info->primal_infeas = std::max(Infnorm(xinfeas), Infnorm(sinfeas));
    info->dual_infeas = std::max(Infnorm(zinfeas), Infnorm(yinfeas));
    info->objval = pobj;
    (void)(n);
    (void)(m);
}

void Model::PostsolveBasis(const std::vector<Int>& basic_status_solver,
                           Int* cbasis_user, Int* vbasis_user) const {
    const Int m = rows();
    const Int n = cols();
    assert((int)basic_status_solver.size() == n+m);

    std::vector<Int> cbasis_temp(num_constr_);
    std::vector<Int> vbasis_temp(num_var_);
    DualizeBackBasis(basic_status_solver, cbasis_temp, vbasis_temp);
    ScaleBackBasis(cbasis_temp, vbasis_temp);
    if (cbasis_user)
        std::copy(std::begin(cbasis_temp), std::end(cbasis_temp), cbasis_user);
    if (vbasis_user)
        std::copy(std::begin(vbasis_temp), std::end(vbasis_temp), vbasis_user);
    (void)(m);
    (void)(n);
}

// Checks if the vectors are valid LP data vectors. Returns 0 if OK and a
// negative value if a vector is invalid.
static int CheckVectors(Int m, Int n, const double* rhs,const char* constr_type,
                        const double* obj, const double* lb, const double* ub) {
    for (Int i = 0; i < m; i++)
        if (!std::isfinite(rhs[i]))
            return -1;
    for (Int j = 0; j < n; j++)
        if (!std::isfinite(obj[j]))
            return -2;
    for (Int j = 0; j < n; j++) {
        if (!std::isfinite(lb[j]) && lb[j] != -INFINITY)
            return -3;
        if (!std::isfinite(ub[j]) && ub[j] != INFINITY)
            return -3;
        if (lb[j] > ub[j])
            return -3;
    }
    for (Int i = 0; i < m; i++)
        if (constr_type[i] != '=' && constr_type[i] != '<' &&
            constr_type[i] != '>')
            return -4;
    return 0;
}

// Checks if A is a valid m-by-n matrix in CSC format. Returns 0 if OK and a
// negative value otherwise.
Int CheckMatrix(Int m, Int n, const Int *Ap, const Int *Ai, const double *Ax) {
    if (Ap[0] != 0)
        return -5;
    for (Int j = 0; j < n; j++)
        if (Ap[j] > Ap[j+1])
            return -5;
    for (Int p = 0; p < Ap[n]; p++)
        if (!std::isfinite(Ax[p]))
            return -6;
    // Test for out of bound indices and duplicates.
    std::vector<Int> marked(m, -1);
    for (Int j = 0; j < n; j++) {
        for (Int p = Ap[j]; p < Ap[j+1]; p++) {
            Int i = Ai[p];
            if (i < 0 || i >= m)
                return -7;
            if (marked[i] == j)
                return -8;
            marked[i] = j;
        }
    }
    return 0;
}

void Model::CopyInput(Int num_constr, Int num_var, const Int* Ap, const Int* Ai,
                      const double* Ax, const double* rhs,
                      const char* constr_type, const double* obj,
                      const double* lbuser, const double* ubuser, Info* info) {
    if (!(Ap && Ai && Ax && rhs && constr_type && obj && lbuser && ubuser)) {
        info->errflag = IPX_ERROR_argument_null;
        return;
    }
    if (num_constr < 0 || num_var <= 0) {
        info->errflag = IPX_ERROR_invalid_dimension;
        return;
    }
    if (CheckVectors(num_constr, num_var, rhs, constr_type, obj, lbuser, ubuser)
        != 0) {
        info->errflag = IPX_ERROR_invalid_vector;
        return;
    }
    if (CheckMatrix(num_constr, num_var, Ap, Ai, Ax) != 0) {
        info->errflag = IPX_ERROR_invalid_matrix;
        return;
    }
    info->errflag = 0;
    num_constr_ = num_constr;
    num_eqconstr_ = std::count(constr_type, constr_type+num_constr, '=');
    num_var_ = num_var;
    num_entries_ = Ap[num_var];
    num_free_var_ = 0;
    boxed_vars_.clear();
    for (Int j = 0; j < num_var; j++) {
        if (std::isinf(lbuser[j]) && std::isinf(ubuser[j]))
            num_free_var_++;
        if (std::isfinite(lbuser[j]) && std::isfinite(ubuser[j]))
            boxed_vars_.push_back(j);
    }
    constr_type_ = std::vector<char>(constr_type, constr_type+num_constr);
    scaled_obj_ = Vector(obj, num_var);
    scaled_rhs_ = Vector(rhs, num_constr);
    scaled_lbuser_ = Vector(lbuser, num_var);
    scaled_ubuser_ = Vector(ubuser, num_var);
    A_.LoadFromArrays(num_constr, num_var, Ap, Ap+1, Ai, Ax);
    norm_obj_ = Infnorm(scaled_obj_);
    norm_rhs_ = Infnorm(scaled_rhs_);
    for (double x : scaled_lbuser_)
        if (std::isfinite(x))
            norm_rhs_ = std::max(norm_rhs_, std::abs(x));
    for (double x : scaled_ubuser_)
        if (std::isfinite(x))
            norm_rhs_ = std::max(norm_rhs_, std::abs(x));
}

void Model::ScaleModel(const Control& control) {
    flipped_vars_.clear();
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(scaled_ubuser_[j]) && std::isinf(scaled_lbuser_[j])) {
            scaled_lbuser_[j] = -scaled_ubuser_[j];
            scaled_ubuser_[j] = INFINITY;
            ScaleColumn(A_, j, -1.0);
            scaled_obj_[j] *= -1.0;
            flipped_vars_.push_back(j);
        }
    }
    colscale_.resize(0);
    rowscale_.resize(0);

    // Choose scaling method.
    if (control.scale() >= 1)
        EquilibrateMatrix();

    // Apply scaling to vectors.
    if (colscale_.size() > 0) {
        assert((int)colscale_.size() == num_var_);
        scaled_obj_ *= colscale_;
        scaled_lbuser_ /= colscale_;
        scaled_ubuser_ /= colscale_;
    }
    if (rowscale_.size() > 0) {
        assert((int)rowscale_.size() == num_constr_);
        scaled_rhs_ *= rowscale_;
    }
}

void Model::LoadPrimal() {
    num_rows_ = num_constr_;
    num_cols_ = num_var_;
    dualized_ = false;

    // Copy A and append identity matrix.
    AI_ = A_;
    for (Int i = 0; i < num_constr_; i++) {
        AI_.push_back(i, 1.0);
        AI_.add_column();
    }
    assert(AI_.cols() == num_var_+num_constr_);

    // Copy vectors and set bounds on slack variables.
    b_ = scaled_rhs_;
    c_.resize(num_var_+num_constr_);
    c_ = 0.0;
    std::copy_n(std::begin(scaled_obj_), num_var_, std::begin(c_));
    lb_.resize(num_rows_+num_cols_);
    std::copy_n(std::begin(scaled_lbuser_), num_var_, std::begin(lb_));
    ub_.resize(num_rows_+num_cols_);
    std::copy_n(std::begin(scaled_ubuser_), num_var_, std::begin(ub_));
    for (Int i = 0; i < num_constr_; i++) {
        switch(constr_type_[i]) {
        case '=':
            lb_[num_var_+i] = 0.0;
            ub_[num_var_+i] = 0.0;
            break;
        case '<':
            lb_[num_var_+i] = 0.0;
            ub_[num_var_+i] = INFINITY;
            break;
        case '>':
            lb_[num_var_+i] = -INFINITY;
            ub_[num_var_+i] = 0.0;
            break;
        }
    }
}

void Model::LoadDual() {
    num_rows_ = num_var_;
    num_cols_ = num_constr_ + boxed_vars_.size();
    dualized_ = true;

    // Check that every variable with finite scaled_ubuser_ has finite
    // scaled_lbuser_ (must be the case after scaling).
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(scaled_ubuser_[j]))
            assert(std::isfinite(scaled_lbuser_[j]));
    }

    // Build AI.
    AI_ = Transpose(A_);
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(scaled_ubuser_[j])) {
            AI_.push_back(j, -1.0);
            AI_.add_column();
        }
    }
    assert(AI_.cols() == num_cols_);
    for (Int i = 0; i < num_rows_; i++) {
        AI_.push_back(i, 1.0);
        AI_.add_column();
    }

    // Build vectors.
    b_ = scaled_obj_;
    c_.resize(num_cols_+num_rows_);
    Int put = 0;
    for (double x : scaled_rhs_)
        c_[put++] = -x;
    for (double x : scaled_ubuser_)
        if (std::isfinite(x))
            c_[put++] = x;
    assert(put == num_cols_);
    for (double x : scaled_lbuser_)
        // If x is negative infinity, then the variable will be fixed and we can
        // give it any (finite) cost.
        c_[put++] = std::isfinite(x) ? -x : 0.0;
    lb_.resize(num_cols_+num_rows_);
    ub_.resize(num_cols_+num_rows_);
    for (Int i = 0; i < num_constr_; i++)
        switch(constr_type_[i]) {
        case '=':
            lb_[i] = -INFINITY;
            ub_[i] = INFINITY;
            break;
        case '<':
            lb_[i] = -INFINITY;
            ub_[i] = 0.0;
            break;
        case '>':
            lb_[i] = 0.0;
            ub_[i] = INFINITY;
            break;
        }
    for (Int j = num_constr_; j < num_cols_; j++) {
        lb_[j] = 0.0;
        ub_[j] = INFINITY;
    }
    for (Int j = 0; j < num_var_; j++) {
        lb_[num_cols_+j] = 0.0;
        ub_[num_cols_+j] = std::isfinite(scaled_lbuser_[j]) ? INFINITY : 0.0;
    }
}

// Returns a power-of-2 factor s such that s*2^exp becomes closer to the
// interval [2^expmin, 2^expmax].
static double EquilibrationFactor(int expmin, int expmax, int exp) {
    if (exp < expmin)
        return std::ldexp(1.0, (expmin-exp+1)/2);
    if (exp > expmax)
        return std::ldexp(1.0, -((exp-expmax+1)/2));
    return 1.0;
}

void Model::EquilibrateMatrix() {
    const Int m = A_.rows();
    const Int n = A_.cols();
    const Int* Ap = A_.colptr();
    const Int* Ai = A_.rowidx();
    double* Ax = A_.values();

    colscale_.resize(0);
    rowscale_.resize(0);

    // The absolute value of each nonzero entry of AI can be written as x*2^exp
    // with x in the range [0.5,1). We consider AI well scaled if each entry is
    // such that expmin <= exp <= expmax for parameters expmin and expmax.
    // For example,
    //
    //  expmin  expmax  allowed range
    //  -----------------------------
    //       0       2    [0.5,    4)
    //      -1       3    [0.25,   8)
    //       1       7    [1.0,  128)
    //
    // If AI is not well scaled, a recursive row and column equilibration is
    // applied. In each iteration, the scaling factors are roughly 1/sqrt(max),
    // where max is the maximum row or column entry. However, the factors are
    // truncated to powers of 2, so that no round-off errors occur.

    constexpr int expmin = 0;
    constexpr int expmax = 3;
    constexpr Int maxround = 10;

    // Quick return if entries are within the target range.
    bool out_of_range = false;
    for (Int p = 0; p < Ap[n]; p++) {
        int exp;
        std::frexp(std::abs(Ax[p]), &exp);
        if (exp < expmin || exp > expmax) {
            out_of_range = true;
            break;
        }
    }
    if (!out_of_range)
        return;

    colscale_.resize(n);
    rowscale_.resize(m);
    colscale_ = 1.0;
    rowscale_ = 1.0;
    Vector colmax(n), rowmax(m);

    for (Int round = 0; round < maxround; round++) {
        // Compute infinity norm of each row and column.
        rowmax = 0.0;
        for (Int j = 0; j < n; j++) {
            colmax[j] = 0.0;
            for (Int p = Ap[j]; p < Ap[j+1]; p++) {
                Int i = Ai[p];
                double xa = std::abs(Ax[p]);
                colmax[j] = std::max(colmax[j], xa);
                rowmax[i] = std::max(rowmax[i], xa);
            }
        }
        // Replace rowmax and colmax entries by scaling factors from this round.
        bool out_of_range = false;
        for (Int i = 0; i < m; i++) {
            int exp;
            std::frexp(rowmax[i], &exp);
            rowmax[i] = EquilibrationFactor(expmin, expmax, exp);
            if (rowmax[i] != 1.0) {
                out_of_range = true;
                rowscale_[i] *= rowmax[i];
            }
        }
        for (Int j = 0; j < n; j++) {
            int exp;
            std::frexp(colmax[j], &exp);
            colmax[j] = EquilibrationFactor(expmin, expmax, exp);
            if (colmax[j] != 1.0) {
                out_of_range = true;
                colscale_[j] *= colmax[j];
            }
        }
        if (!out_of_range)
            break;
        // Rescale A.
        for (Int j = 0; j < n; j++) {
            for (Int p = Ap[j]; p < Ap[j+1]; p++) {
                Ax[p] *= colmax[j];     // column scaling
                Ax[p] *= rowmax[Ai[p]]; // row scaling
            }
        }
    }
}

void Model::FindDenseColumns() {
    num_dense_cols_ = 0;
    nz_dense_ = rows() + 1;

    std::vector<Int> colcount(num_cols_);
    for (Int j = 0; j < num_cols_; j++)
        colcount[j] = AI_.end(j)-AI_.begin(j);
    std::sort(colcount.begin(), colcount.end());

    for (Int j = 1; j < num_cols_; j++) {
        // TODO: Magic numbers 401, 101 here:
        if (colcount[j] > std::max((Int)401, (Int)(101*colcount[j-1]))) {
            // j is the first dense column
            num_dense_cols_ = num_cols_ - j;
            nz_dense_ = colcount[j];
            break;
        }
    }

    if (num_dense_cols_ > 1000) {
        num_dense_cols_ = 0;
        nz_dense_ = rows() + 1;
    }
}

void Model::PrintCoefficientRange(const Control& control) const {
    double amin = INFINITY;
    double amax = 0.0;
    for (Int j = 0; j < A_.cols(); j++) {
        for (Int p = A_.begin(j); p < A_.end(j); p++) {
            double x = std::abs(A_.value(p));
            if (x != 0.0) {
                amin = std::min(amin, x);
                amax = std::max(amax, x);
            }
        }
    }
    if (amin == INFINITY)       // no nonzero entries in A_
        amin = 0.0;
    control.Log()
        << Textline("Matrix range:")
        << "[" << Scientific(amin, 5, 0) << ", "
        << Scientific(amax, 5, 0) << "]\n";

    double rhsmin = INFINITY;
    double rhsmax = 0.0;
    for (double x : scaled_rhs_) {
        if (x != 0.0) {
            rhsmin = std::min(rhsmin, std::abs(x));
            rhsmax = std::max(rhsmax, std::abs(x));
        }
    }
    if (rhsmin == INFINITY)     // no nonzero entries in rhs
        rhsmin = 0.0;
    control.Log()
        << Textline("RHS range:")
        << "[" << Scientific(rhsmin, 5, 0) << ", "
        << Scientific(rhsmax, 5, 0) << "]\n";

    double objmin = INFINITY;
    double objmax = 0.0;
    for (double x : scaled_obj_) {
        if (x != 0.0) {
            objmin = std::min(objmin, std::abs(x));
            objmax = std::max(objmax, std::abs(x));
        }
    }
    if (objmin == INFINITY)     // no nonzero entries in obj
        objmin = 0.0;
    control.Log()
        << Textline("Objective range:")
        << "[" << Scientific(objmin, 5, 0) << ", "
        << Scientific(objmax, 5, 0) << "]\n";

    double boundmin = INFINITY;
    double boundmax = 0.0;
    for (double x : scaled_lbuser_) {
        if (x != 0.0 && std::isfinite(x)) {
            boundmin = std::min(boundmin, std::abs(x));
            boundmax = std::max(boundmax, std::abs(x));
        }
    }
    for (double x : scaled_ubuser_) {
        if (x != 0.0 && std::isfinite(x)) {
            boundmin = std::min(boundmin, std::abs(x));
            boundmax = std::max(boundmax, std::abs(x));
        }
    }
    if (boundmin == INFINITY)   // no finite nonzeros entries in bounds
        boundmin = 0.0;
    control.Log()
        << Textline("Bounds range:")
        << "[" << Scientific(boundmin, 5, 0) << ", "
        << Scientific(boundmax, 5, 0) << "]\n";
}

void Model::PrintPreprocessingLog(const Control& control) const {
    // Find the minimum and maximum scaling factor.
    double minscale = INFINITY;
    double maxscale = 0.0;
    if (colscale_.size() > 0) {
        auto minmax = std::minmax_element(std::begin(colscale_),
                                          std::end(colscale_));
        minscale = std::min(minscale, *minmax.first);
        maxscale = std::max(maxscale, *minmax.second);
    }
    if (rowscale_.size() > 0) {
        auto minmax = std::minmax_element(std::begin(rowscale_),
                                          std::end(rowscale_));
        minscale = std::min(minscale, *minmax.first);
        maxscale = std::max(maxscale, *minmax.second);
    }
    if (minscale == INFINITY)
        minscale = 1.0;
    if (maxscale == 0.0)
        maxscale = 1.0;

    control.Log()
        << "Preprocessing\n"
        << Textline("Dualized model:") << (dualized() ? "yes" : "no") << '\n'
        << Textline("Number of dense columns:") << num_dense_cols() << '\n';
    if (control.scale() > 0) {
        control.Log()
            << Textline("Range of scaling factors:") << "["
            << Scientific(minscale, 8, 2) << ", "
            << Scientific(maxscale, 8, 2) << "]\n";
    }
}

void Model::WriteInfo(Info *info) const {
    info->num_var = num_var_;
    info->num_constr = num_constr_;
    info->num_entries = num_entries_;
    info->num_rows_solver = num_rows_;
    info->num_cols_solver = num_cols_ + num_rows_; // including slack columns
    info->num_entries_solver = AI_.entries();
    info->dualized = dualized_;
    info->dense_cols = num_dense_cols();
}

void Model::ScaleBasicSolution(Vector& x, Vector& slack, Vector& y, Vector& z)
    const {
    if (colscale_.size() > 0) {
        x /= colscale_;
        z *= colscale_;
    }
    if (rowscale_.size() > 0) {
        y /= rowscale_;
        slack *= rowscale_;
    }
    for (Int j : flipped_vars_) {
        x[j] *= -1.0;
        z[j] *= -1.0;
    }
}

void Model::ScaleBackInteriorSolution(Vector& x, Vector& xl, Vector& xu,
                                      Vector& slack, Vector& y, Vector& zl,
                                      Vector& zu) const {
    if (colscale_.size() > 0) {
        x *= colscale_;
        xl *= colscale_;
        xu *= colscale_;
        zl /= colscale_;
        zu /= colscale_;
    }
    if (rowscale_.size() > 0) {
        y *= rowscale_;
        slack /= rowscale_;
    }
    for (Int j : flipped_vars_) {
        assert(std::isfinite(xl[j]));
        assert(std::isinf(xu[j]));
        assert(zu[j] == 0.0);
        x[j] *= -1.0;
        xu[j] = xl[j];
        xl[j] = INFINITY;
        zu[j] = zl[j];
        zl[j] = 0.0;
    }
}

void Model::ScaleBackResiduals(Vector& rb, Vector& rc, Vector& rl,
                               Vector& ru) const {
    if (colscale_.size() > 0) {
        rc /= colscale_;
        rl *= colscale_;
        ru *= colscale_;
    }
    if (rowscale_.size() > 0)
        rb /= rowscale_;
    for (Int j : flipped_vars_) {
        rc[j] *= -1.0;
        assert(ru[j] == 0.0);
        ru[j] = -rl[j];
        rl[j] = 0.0;
    }
}

void Model::ScaleBackBasicSolution(Vector& x, Vector& slack, Vector& y,
                                   Vector& z) const {
    if (colscale_.size() > 0) {
        x *= colscale_;
        z /= colscale_;
    }
    if (rowscale_.size() > 0) {
        y *= rowscale_;
        slack /= rowscale_;
    }
    for (Int j : flipped_vars_) {
        x[j] *= -1.0;
        z[j] *= -1.0;
    }
}

void Model::ScaleBackBasis(std::vector<Int>& cbasis,
                           std::vector<Int>& vbasis) const {
    for (Int j : flipped_vars_) {
        assert(vbasis[j] != IPX_nonbasic_ub);
        if (vbasis[j] == IPX_nonbasic_lb)
            vbasis[j] = IPX_nonbasic_ub;
    }
}

void Model::DualizeBasicSolution(const Vector& x_user,
                                 const Vector& slack_user,
                                 const Vector& y_user,
                                 const Vector& z_user,
                                 Vector& x_solver,
                                 Vector& y_solver,
                                 Vector& z_solver) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + (int)boxed_vars_.size() == n);

        // Build dual solver variables from primal user variables.
        y_solver = -x_user;
        for (Int i = 0; i < num_constr_; i++)
            z_solver[i] = -slack_user[i];
        for (Int k = 0; k < (Int) boxed_vars_.size(); k++) {
            Int j = boxed_vars_[k];
            z_solver[num_constr_+k] = c(num_constr_+k) + y_solver[j];
        }
        for (Int i = 0; i < m; i++)
            z_solver[n+i] = c(n+i)-y_solver[i];

        // Build primal solver variables from dual user variables.
        std::copy_n(std::begin(y_user), num_constr_, std::begin(x_solver));
        std::copy_n(std::begin(z_user), num_var_, std::begin(x_solver) + n);
        for (Int k = 0; k < (Int) boxed_vars_.size(); k++) {
            Int j = boxed_vars_[k];
            if (x_solver[n+j] < 0.0) {
                // j is a boxed variable and z_user[j] < 0
                x_solver[num_constr_+k] = -x_solver[n+j];
                x_solver[n+j] = 0.0;
            } else {
                x_solver[num_constr_+k] = 0.0;
            }
        }
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        std::copy_n(std::begin(x_user), n, std::begin(x_solver));
        std::copy_n(std::begin(slack_user), m, std::begin(x_solver) + n);
        std::copy_n(std::begin(y_user), m, std::begin(y_solver));
        std::copy_n(std::begin(z_user), n, std::begin(z_solver));
        for (Int i = 0; i < m; i++)
            z_solver[n+i] = c(n+i)-y_solver[i];
    }
}

void Model::DualizeBackInteriorSolution(const Vector& x_solver,
                                        const Vector& xl_solver,
                                        const Vector& xu_solver,
                                        const Vector& y_solver,
                                        const Vector& zl_solver,
                                        const Vector& zu_solver,
                                        Vector& x_user,
                                        Vector& xl_user,
                                        Vector& xu_user,
                                        Vector& slack_user,
                                        Vector& y_user,
                                        Vector& zl_user,
                                        Vector& zu_user) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + (int)boxed_vars_.size() == n);
        x_user = -y_solver;

        // If the solution from the solver would be exact, we could copy the
        // first num_constr_ entries from x_solver into y_user. However, to
        // satisfy the sign condition on y_user even if the solution is not
        // exact, we have to use the xl_solver and xu_solver entries for
        // inequality constraints.
        for (Int i = 0; i < num_constr_; i++) {
            switch (constr_type_[i]) {
            case '=':
                y_user[i] = x_solver[i];
                break;
            case '<':
                y_user[i] = -xu_solver[i];
                break;
            case '>':
                y_user[i] = xl_solver[i];
                break;
            }
            assert(std::isfinite(y_user[i]));
        }

        // Dual variables associated with lbuser <= x in the scaled user model
        // are the slack variables from the solver. For an exact solution we
        // would have x_solver[n+1:n+m] == xl_solver[n+1:n+m]. Using xl_solver
        // guarantees that zl_user >= 0 in any case. If variable j has no lower
        // bound in the scaled user model (i.e. is free), then the j-th slack
        // variable was fixed at zero in the solver model, but the IPM solution
        // may not satisfy this. Hence we must set zl_user[j] = 0 explicitly.
        std::copy_n(std::begin(xl_solver) + n, num_var_, std::begin(zl_user));
        for (Int j = 0; j < num_var_; j++)
            if (!std::isfinite(scaled_lbuser_[j]))
                zl_user[j] = 0.0;

        // Dual variables associated with x <= ubuser in the scaled user model
        // are the primal variables that were added for boxed variables in the
        // solver model.
        zu_user = 0.0;
        Int k = num_constr_;
        for (Int j : boxed_vars_)
            zu_user[j] = xl_solver[k++];
        assert(k == n);

        // xl in the scaled user model is zl[n+1:n+m] in the solver model or
        // infinity.
        for (Int i = 0; i < m; i++) {
            if (std::isfinite(scaled_lbuser_[i]))
                xl_user[i] = zl_solver[n+i];
            else
                xl_user[i] = INFINITY;
        }

        // xu in the scaled user model are the entries in zl for columns of the
        // negative identity matrix (that were added for boxed variables).
        xu_user = INFINITY;
        k = num_constr_;
        for (Int j : boxed_vars_)
            xu_user[j] = zl_solver[k++];
        assert(k == n);
        
        for (Int i = 0; i < num_constr_; i++) {
            switch (constr_type_[i]) {
            case '=':
                slack_user[i] = 0.0;
                break;
            case '<':
                slack_user[i] = zu_solver[i];
                break;
            case '>':
                slack_user[i] = -zl_solver[i];
                break;
            }
        }
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        std::copy_n(std::begin(x_solver), num_var_, std::begin(x_user));

        // Instead of copying y_solver into y_user, we use the entries from
        // zl_solver and zu_solver for inequality constraints, so that the sign
        // condition on y_user is satisfied.
        for (Int i = 0; i < m; i++) {
            assert(lb_[n+i] == 0.0 || lb_[n+i] == -INFINITY);
            assert(ub_[n+i] == 0.0 || ub_[n+i] ==  INFINITY);
            assert(lb_[n+i] == 0.0 || ub_[n+i] == 0.0);
            switch (constr_type_[i]) {
            case '=':
                y_user[i] = y_solver[i];
                break;
            case '<':
                y_user[i] = -zl_solver[n+i];
                break;
            case '>':
                y_user[i] = zu_solver[n+i];
                break;
            }
            assert(std::isfinite(y_user[i]));
        }
        std::copy_n(std::begin(zl_solver), num_var_, std::begin(zl_user));
        std::copy_n(std::begin(zu_solver), num_var_, std::begin(zu_user));
        std::copy_n(std::begin(xl_solver), num_var_, std::begin(xl_user));
        std::copy_n(std::begin(xu_solver), num_var_, std::begin(xu_user));

        // If the solution would be exact, slack_user were given by the entries
        // of x_solver corresponding to slack columns. To satisfy the sign
        // condition in any case, we build the slack for inequality constraints
        // from xl_solver and xu_solver and set the slack for equality
        // constraints to zero.
        for (Int i = 0; i < m; i++) {
            switch (constr_type_[i]) {
            case '=':
                slack_user[i] = 0.0;
                break;
            case '<':
                slack_user[i] = xl_solver[n+i];
                break;
            case '>':
                slack_user[i] = -xu_solver[n+i];
                break;
            }
            assert(std::isfinite(slack_user[i]));
        }
    }
}

void Model::DualizeBackBasicSolution(const Vector& x_solver,
                                     const Vector& y_solver,
                                     const Vector& z_solver,
                                     Vector& x_user,
                                     Vector& slack_user,
                                     Vector& y_user,
                                     Vector& z_user) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + (int)boxed_vars_.size() == n);
        x_user = -y_solver;
        for (Int i = 0; i < num_constr_; i++)
            slack_user[i] = -z_solver[i];
        std::copy_n(std::begin(x_solver), num_constr_, std::begin(y_user));
        std::copy_n(std::begin(x_solver) + n, num_var_, std::begin(z_user));
        Int k = num_constr_;
        for (Int j : boxed_vars_)
            z_user[j] -= x_solver[k++];
        assert(k == n);
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        std::copy_n(std::begin(x_solver), num_var_, std::begin(x_user));
        std::copy_n(std::begin(x_solver) + n, num_constr_,
                    std::begin(slack_user));
        std::copy_n(std::begin(y_solver), num_constr_, std::begin(y_user));
        std::copy_n(std::begin(z_solver), num_var_, std::begin(z_user));
    }
    (void)(m);
}

void Model::DualizeBackBasis(const std::vector<Int>& basic_status_solver,
                             std::vector<Int>& cbasis_user,
                             std::vector<Int>& vbasis_user) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + (int)boxed_vars_.size() == n);
        for (Int i = 0; i < num_constr_; i++) {
            if (basic_status_solver[i] == IPX_basic)
                cbasis_user[i] = IPX_nonbasic;
            else
                cbasis_user[i] = IPX_basic;
        }
        for (Int j = 0; j < num_var_; j++) {
            // slack cannot be superbasic
            assert(basic_status_solver[n+j] != IPX_superbasic);
            if (basic_status_solver[n+j] == 0)
                vbasis_user[j] = std::isfinite(scaled_lbuser_[j]) ?
                    IPX_nonbasic_lb : IPX_superbasic;
            else
                vbasis_user[j] = IPX_basic;
        }
        Int k = num_constr_;
        for (Int j : boxed_vars_)
            if (basic_status_solver[k++] == IPX_basic) {
                assert(vbasis_user[j] == IPX_basic);
                vbasis_user[j] = IPX_nonbasic_ub;
            }
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        for (Int i = 0; i < num_constr_; i++) {
            // slack cannot be superbasic
            assert(basic_status_solver[n+i] != IPX_superbasic);
            if (basic_status_solver[n+i] == IPX_basic)
                cbasis_user[i] = IPX_basic;
            else
                cbasis_user[i] = IPX_nonbasic;
        }
        for (Int j = 0; j < num_var_; j++)
            vbasis_user[j] = basic_status_solver[j];
    }
    (void)(m);
}

void Model::CorrectScaledBasicSolution(Vector& x, Vector& slack, Vector& y,
                                       Vector& z,
                                       const std::vector<Int> cbasis,
                                       const std::vector<Int> vbasis) const {
    for (Int j = 0; j < num_var_; j++) {
        if (vbasis[j] == IPX_nonbasic_lb)
            x[j] = scaled_lbuser_[j];
        if (vbasis[j] == IPX_nonbasic_ub)
            x[j] = scaled_ubuser_[j];
        if (vbasis[j] == IPX_basic)
            z[j] = 0.0;
    }
    for (Int i = 0; i < num_constr_; i++) {
        if (cbasis[i] == IPX_nonbasic)
            slack[i] = 0.0;
        if (cbasis[i] == IPX_basic)
            y[i] = 0.0;
    }
}

void Model::MultiplyWithScaledMatrix(const Vector& rhs, double alpha,
                                     Vector& lhs, char trans) const {
    if (trans == 't' || trans == 'T') {
        assert((int)rhs.size() == num_constr_);
        assert((int)lhs.size() == num_var_);
        if (dualized())
            for (Int i = 0; i < num_constr_; i++)
                ScatterColumn(AI_, i, alpha*rhs[i], lhs);
        else
            for (Int j = 0; j < num_var_; j++)
                lhs[j] += alpha * DotColumn(AI_, j, rhs);
    }
    else {
        assert((int)rhs.size() == num_var_);
        assert((int)lhs.size() == num_constr_);
        if (dualized())
            for (Int i = 0; i < num_constr_; i++)
                lhs[i] += alpha * DotColumn(AI_, i, rhs);
        else
            for (Int j = 0; j < num_var_; j++)
                ScatterColumn(AI_, j, alpha*rhs[j], lhs);
    }
}

double PrimalInfeasibility(const Model& model, const Vector& x) {
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    assert(x.size() == lb.size());

    double infeas = 0.0;
    for (Int j = 0; j < (Int) x.size(); j++) {
        infeas = std::max(infeas, lb[j]-x[j]);
        infeas = std::max(infeas, x[j]-ub[j]);
    }
    return infeas;
}

double DualInfeasibility(const Model& model, const Vector& x,
                                const Vector& z) {
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    assert(x.size() == lb.size());
    assert(z.size() == lb.size());

    double infeas = 0.0;
    for (Int j = 0; j < (Int) x.size(); j++) {
        if (x[j] > lb[j])
            infeas = std::max(infeas, z[j]);
        if (x[j] < ub[j])
            infeas = std::max(infeas, -z[j]);
    }
    return infeas;
}

double PrimalResidual(const Model& model, const Vector& x) {
    const SparseMatrix& AIt = model.AIt();
    const Vector& b = model.b();
    assert((int)x.size() == AIt.rows());

    double res = 0.0;
    for (Int i = 0; i < (Int) b.size(); i++) {
        double r = b[i] - DotColumn(AIt, i, x);
        res = std::max(res, std::abs(r));
    }
    return res;
}

double DualResidual(const Model& model, const Vector& y, const Vector& z) {
    const SparseMatrix& AI = model.AI();
    const Vector& c = model.c();
    assert((int)y.size() == AI.rows());
    assert((int)z.size() == AI.cols());

    double res = 0.0;
    for (Int j = 0; j < (Int) c.size(); j++) {
        double r = c[j] - z[j] - DotColumn(AI, j, y);
        res = std::max(res, std::abs(r));
    }
    return res;
}

}  // namespace ipx
