//! Affine transformations for geometry manipulation.
//!
//! Provides scale, rotate, translate, and mirror operations that can be
//! applied to sets of points or to primitives. These are used both in the
//! TOML configuration (to position objects in a scene) and interactively
//! in the GUI.

use nalgebra::{Matrix3, Vector3};

/// An affine transformation: rotation/scale matrix + translation.
#[derive(Debug, Clone)]
pub struct Transform {
    /// 3x3 rotation/scale matrix.
    pub matrix: Matrix3<f64>,
    /// Translation vector (nm).
    pub translation: Vector3<f64>,
}

impl Default for Transform {
    fn default() -> Self {
        Self {
            matrix: Matrix3::identity(),
            translation: Vector3::zeros(),
        }
    }
}

impl Transform {
    /// Create a pure translation.
    pub fn translation(dx: f64, dy: f64, dz: f64) -> Self {
        Self {
            matrix: Matrix3::identity(),
            translation: Vector3::new(dx, dy, dz),
        }
    }

    /// Create a uniform scale about the origin.
    pub fn uniform_scale(factor: f64) -> Self {
        Self {
            matrix: Matrix3::identity() * factor,
            translation: Vector3::zeros(),
        }
    }

    /// Create a non-uniform scale.
    pub fn scale(sx: f64, sy: f64, sz: f64) -> Self {
        Self {
            matrix: Matrix3::from_diagonal(&Vector3::new(sx, sy, sz)),
            translation: Vector3::zeros(),
        }
    }

    /// Apply this transformation to a 3D point.
    pub fn apply(&self, point: &[f64; 3]) -> [f64; 3] {
        let v = Vector3::new(point[0], point[1], point[2]);
        let result = self.matrix * v + self.translation;
        [result.x, result.y, result.z]
    }

    /// Compose two transforms: self followed by other.
    pub fn then(&self, other: &Transform) -> Transform {
        Transform {
            matrix: other.matrix * self.matrix,
            translation: other.matrix * self.translation + other.translation,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_transform() {
        let t = Transform::default();
        let p = [1.0, 2.0, 3.0];
        let result = t.apply(&p);
        assert!((result[0] - 1.0).abs() < 1e-12);
        assert!((result[1] - 2.0).abs() < 1e-12);
        assert!((result[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_scale_and_translate() {
        let t = Transform::uniform_scale(2.0).then(&Transform::translation(1.0, 0.0, 0.0));
        let p = [1.0, 1.0, 1.0];
        let result = t.apply(&p);
        assert!((result[0] - 3.0).abs() < 1e-12);
        assert!((result[1] - 2.0).abs() < 1e-12);
        assert!((result[2] - 2.0).abs() < 1e-12);
    }
}
