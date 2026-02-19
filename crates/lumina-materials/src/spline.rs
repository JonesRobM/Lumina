//! Cubic spline interpolation for smooth material property curves.
//!
//! Tabulated material data (e.g. Johnson & Christy) is provided at discrete
//! wavelengths. Cubic spline interpolation provides a smooth, continuous
//! $\epsilon(\lambda)$ between data points, which is essential for computing
//! fine spectral features.

/// A natural cubic spline interpolator for real-valued data.
///
/// Given $n$ data points $(x_i, y_i)$, constructs piecewise cubic polynomials
/// with continuous first and second derivatives.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    /// Sorted x values (knots).
    xs: Vec<f64>,
    /// Corresponding y values.
    ys: Vec<f64>,
    /// Second derivatives at each knot (computed during construction).
    y2s: Vec<f64>,
}

impl CubicSpline {
    /// Construct a natural cubic spline from data points.
    ///
    /// # Arguments
    /// * `xs` - Strictly increasing x values.
    /// * `ys` - Corresponding y values (same length as `xs`).
    ///
    /// # Panics
    /// Panics if `xs` and `ys` have different lengths, or if `xs` is not
    /// strictly increasing, or if fewer than 2 points are provided.
    pub fn new(xs: Vec<f64>, ys: Vec<f64>) -> Self {
        assert_eq!(xs.len(), ys.len(), "xs and ys must have equal length");
        assert!(xs.len() >= 2, "Need at least 2 data points");
        for i in 1..xs.len() {
            assert!(
                xs[i] > xs[i - 1],
                "xs must be strictly increasing at index {}",
                i
            );
        }

        let n = xs.len();
        let mut y2s = vec![0.0; n];
        let mut u = vec![0.0; n - 1];

        // Forward sweep (tridiagonal system for natural spline)
        for i in 1..n - 1 {
            let sig = (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1]);
            let p = sig * y2s[i - 1] + 2.0;
            y2s[i] = (sig - 1.0) / p;
            u[i] = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
                - (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
            u[i] = (6.0 * u[i] / (xs[i + 1] - xs[i - 1]) - sig * u[i - 1]) / p;
        }

        // Back substitution
        for k in (0..n - 2).rev() {
            y2s[k + 1] = y2s[k + 1] * y2s[k + 2] + u[k + 1];
        }

        Self { xs, ys, y2s }
    }

    /// Evaluate the spline at a given x value.
    ///
    /// Extrapolation beyond the data range uses the boundary polynomial.
    pub fn evaluate(&self, x: f64) -> f64 {
        let n = self.xs.len();

        // Binary search for the enclosing interval
        let mut lo = 0;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if self.xs[mid] > x {
                hi = mid;
            } else {
                lo = mid;
            }
        }

        let h = self.xs[hi] - self.xs[lo];
        let a = (self.xs[hi] - x) / h;
        let b = (x - self.xs[lo]) / h;

        a * self.ys[lo]
            + b * self.ys[hi]
            + ((a * a * a - a) * self.y2s[lo] + (b * b * b - b) * self.y2s[hi]) * h * h / 6.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spline_passes_through_data_points() {
        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = vec![2.0, 3.0, 5.0, 4.0, 1.0];
        let spline = CubicSpline::new(xs.clone(), ys.clone());

        for (x, y) in xs.iter().zip(ys.iter()) {
            let result = spline.evaluate(*x);
            assert!(
                (result - y).abs() < 1e-10,
                "Spline({}) = {} but expected {}",
                x,
                result,
                y
            );
        }
    }
}
