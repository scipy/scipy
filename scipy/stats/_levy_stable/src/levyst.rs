//! Implements special functions for stable distribution calculations.
//!
//! A function g appears in the integrand in Nolan's method for calculating
//! stable densities and distribution functions. It takes a different form for
//! alpha = 1 vs alpha ≠ 1. See [NO] for more info.
//! 
//! Notes on translation to Rust:
//! The original C struct stores a `g` function pointer selected at
//! construction time (`g_alpha_ne_one` vs `g_alpha_eq_one`); here that's a
//! bool flag dispatched in `g()` instead, which is the more idiomatic Rust
//! shape and avoids the C version's raw pointer + manual malloc/free.
//!
//! References
//! [NO] John P. Nolan (1997) Numerical calculation of stable densities and
//! distribution functions.

use pyo3::prelude::*;
use std::f64::consts::{FRAC_1_PI, FRAC_2_PI, FRAC_PI_2};

#[pyclass]
struct Nolan {
    alpha: f64,
    zeta: f64,
    xi: f64,
    zeta_prefactor: f64,
    alpha_exp: f64,
    alpha_xi: f64,
    zeta_offset: f64,
    two_beta_div_pi: f64,
    pi_div_two_beta: f64,
    x0_div_term: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    alpha_eq_one: bool,
}

fn g_alpha_ne_one(sp: &Nolan, theta: f64) -> f64 {
    if theta == -sp.xi {
        return if sp.alpha < 1.0 { 0.0 } else { f64::INFINITY };
    }
    if theta == FRAC_PI_2 {
        return if sp.alpha < 1.0 { f64::INFINITY } else { 0.0 };
    }

    let cos_theta = theta.cos();
    sp.zeta_prefactor
        * (cos_theta / (sp.alpha_xi + sp.alpha * theta).sin() * sp.zeta_offset)
            .powf(sp.alpha_exp)
        * (sp.alpha_xi + (sp.alpha - 1.0) * theta).cos()
        / cos_theta
}

fn g_alpha_eq_one(sp: &Nolan, theta: f64) -> f64 {
    if theta == -sp.xi {
        return 0.0;
    }
    if theta == FRAC_PI_2 {
        return f64::INFINITY;
    }

    (1.0 + theta * sp.two_beta_div_pi)
        * ((sp.pi_div_two_beta + theta) * theta.tan() - sp.x0_div_term).exp()
        / theta.cos()
}

#[pymethods]
impl Nolan {
    #[new]
    fn new(alpha: f64, beta: f64, x0: f64) -> Self {
        // Stores results of intermediate computations so they need not be
        // recomputed when g is called many times during numerical
        // integration through QUADPACK.
        let zeta = -beta * (FRAC_PI_2 * alpha).tan();

        if alpha != 1.0 {
            let xi = (-zeta).atan() / alpha;
            let zeta_prefactor = (zeta.powi(2) + 1.0).powf(-1.0 / (2.0 * (alpha - 1.0)));
            let alpha_exp = alpha / (alpha - 1.0);
            let alpha_xi = (-zeta).atan();
            let zeta_offset = x0 - zeta;
            let (c1, c3) = if alpha < 1.0 {
                (0.5 - xi * FRAC_1_PI, FRAC_1_PI)
            } else {
                (1.0, -FRAC_1_PI)
            };
            let c2 = alpha * FRAC_1_PI / (alpha - 1.0).abs() / (x0 - zeta);

            Nolan {
                alpha,
                zeta,
                xi,
                zeta_prefactor,
                alpha_exp,
                alpha_xi,
                zeta_offset,
                two_beta_div_pi: 0.0,
                pi_div_two_beta: 0.0,
                x0_div_term: 0.0,
                c1,
                c2,
                c3,
                alpha_eq_one: false,
            }
        } else {
            let xi = FRAC_PI_2;
            let two_beta_div_pi = beta * FRAC_2_PI;
            let pi_div_two_beta = FRAC_PI_2 / beta;
            let x0_div_term = x0 / two_beta_div_pi;
            let c2 = 0.5 / beta.abs();

            Nolan {
                alpha,
                zeta,
                xi,
                zeta_prefactor: 0.0,
                alpha_exp: 0.0,
                alpha_xi: 0.0,
                zeta_offset: 0.0,
                two_beta_div_pi,
                pi_div_two_beta,
                x0_div_term,
                c1: 0.0,
                c2,
                c3: FRAC_1_PI,
                alpha_eq_one: true,
            }
        }
    }

    fn g(&self, theta: f64) -> f64 {
        if self.alpha_eq_one {
            g_alpha_eq_one(self, theta)
        } else {
            g_alpha_ne_one(self, theta)
        }
    }

    #[getter]
    fn zeta(&self) -> f64 {
        self.zeta
    }

    #[getter]
    fn xi(&self) -> f64 {
        self.xi
    }

    #[getter]
    fn c1(&self) -> f64 {
        self.c1
    }

    #[getter]
    fn c2(&self) -> f64 {
        self.c2
    }

    #[getter]
    fn c3(&self) -> f64 {
        self.c3
    }
}

#[pymodule]
fn levyst(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Nolan>()?;
    Ok(())
}
