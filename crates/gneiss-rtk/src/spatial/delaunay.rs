//! 2D Delaunay Triangulation and Barycentric Spatial Interpolation.
//!
//! Implements incremental Bowyer-Watson Delaunay triangulation with
//! super-triangle initialization, point-in-triangle location, barycentric
//! interpolation, and robust out-of-hull inverse distance fallback.

use std::fmt;
use nalgebra::Vector2;

/// Error type for spatial Delaunay operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DelaunayError {
    InsufficientPoints(usize),
    DegenerateMesh(String),
    InterpolationFailed(String),
}

impl fmt::Display for DelaunayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InsufficientPoints(n) => write!(f, "Insufficient points for triangulation: {n} < 3"),
            Self::DegenerateMesh(msg) => write!(f, "Degenerate mesh: {msg}"),
            Self::InterpolationFailed(msg) => write!(f, "Interpolation failed: {msg}"),
        }
    }
}

impl std::error::Error for DelaunayError {}

/// Alias for compatibility with project interface contract.
pub type EngineError = DelaunayError;

/// 2D point for spatial triangulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

impl Point2D {
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    #[inline]
    pub fn dist_sq(&self, other: &Point2D) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }
}

impl From<Vector2<f64>> for Point2D {
    fn from(v: Vector2<f64>) -> Self {
        Self { x: v.x, y: v.y }
    }
}

impl From<Point2D> for Vector2<f64> {
    fn from(p: Point2D) -> Self {
        Vector2::new(p.x, p.y)
    }
}

/// Triangle represented by indices of its three vertices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Triangle {
    pub vertices: [usize; 3],
}

impl Triangle {
    #[inline]
    pub fn new(a: usize, b: usize, c: usize) -> Self {
        Self { vertices: [a, b, c] }
    }

    #[inline]
    pub fn contains_vertex(&self, v: usize) -> bool {
        self.vertices[0] == v || self.vertices[1] == v || self.vertices[2] == v
    }

    #[inline]
    pub fn edges(&self) -> [(usize, usize); 3] {
        [
            (self.vertices[0], self.vertices[1]),
            (self.vertices[1], self.vertices[2]),
            (self.vertices[2], self.vertices[0]),
        ]
    }
}

/// 2D Delaunay triangulation mesh.
#[derive(Debug, Clone)]
pub struct Delaunay2D {
    pub points: Vec<Point2D>,
    pub triangles: Vec<Triangle>,
}

pub type DelaunayMesh = Delaunay2D;

impl Delaunay2D {
    /// Construct Delaunay triangulation from slice of Vector2 points.
    pub fn new(points: &[Vector2<f64>]) -> Result<Self, EngineError> {
        let pts: Vec<Point2D> = points.iter().map(|p| Point2D::from(*p)).collect();
        Self::from_points(&pts)
    }

    /// Construct Delaunay triangulation from slice of Point2D points.
    pub fn from_points(points: &[Point2D]) -> Result<Self, EngineError> {
        let deduped = deduplicate_points(points);
        if deduped.len() < 3 {
            return Err(EngineError::InsufficientPoints(deduped.len()));
        }
        if are_points_collinear(&deduped) {
            return Err(EngineError::DegenerateMesh("All points are collinear".into()));
        }
        let triangles = triangulate_bowyer_watson(&deduped)?;
        Ok(Self { points: deduped, triangles })
    }

    /// Locate triangle enclosing query point, returning triangle index and barycentric weights.
    pub fn locate_triangle(&self, p: Point2D) -> Option<(usize, [f64; 3])> {
        for (idx, tri) in self.triangles.iter().enumerate() {
            let a = self.points[tri.vertices[0]];
            let b = self.points[tri.vertices[1]];
            let c = self.points[tri.vertices[2]];
            if let Some(weights) = compute_barycentric(a, b, c, p) {
                return Some((idx, weights));
            }
        }
        None
    }

    /// Interpolate values at query point using barycentric weights or IDW fallback.
    pub fn interpolate_barycentric(&self, query: &Vector2<f64>, values: &[f64]) -> Option<f64> {
        if values.len() != self.points.len() || self.points.is_empty() {
            return None;
        }
        let p = Point2D::from(*query);
        Some(self.interpolate(values, p))
    }

    /// Interpolate scalar field across mesh with automatic out-of-hull fallback.
    pub fn interpolate(&self, values: &[f64], p: Point2D) -> f64 {
        if let Some((idx, w)) = self.locate_triangle(p) {
            let tri = &self.triangles[idx];
            return w[0] * values[tri.vertices[0]]
                + w[1] * values[tri.vertices[1]]
                + w[2] * values[tri.vertices[2]];
        }
        fallback_inverse_distance(&self.points, values, p)
    }
}

/// Deduplicate points within spatial tolerance.
fn deduplicate_points(points: &[Point2D]) -> Vec<Point2D> {
    let mut result: Vec<Point2D> = Vec::with_capacity(points.len());
    for p in points {
        let exists = result.iter().any(|existing| existing.dist_sq(p) < 1e-12);
        if !exists {
            result.push(*p);
        }
    }
    result
}

/// Check if all points lie along a single 2D line.
fn are_points_collinear(points: &[Point2D]) -> bool {
    if points.len() < 3 {
        return true;
    }
    let p0 = points[0];
    let p1 = points[1];
    for p2 in &points[2..] {
        let area = orient2d(p0, p1, *p2).abs();
        if area > 1e-9 {
            return false;
        }
    }
    true
}

/// 2D cross product: positive if A->B->C is counter-clockwise.
#[inline]
pub fn orient2d(a: Point2D, b: Point2D, c: Point2D) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Bowyer-Watson incremental Delaunay triangulation.
fn triangulate_bowyer_watson(points: &[Point2D]) -> Result<Vec<Triangle>, EngineError> {
    let n = points.len();
    let (super_pts, super_tri) = make_super_triangle(points);
    let mut all_points = points.to_vec();
    all_points.extend_from_slice(&super_pts);

    let mut triangles = vec![super_tri];
    for (pt_idx, p) in points.iter().enumerate() {
        insert_point_bw(&mut triangles, &all_points, pt_idx, *p);
    }
    triangles.retain(|tri| !tri.contains_vertex(n) && !tri.contains_vertex(n + 1) && !tri.contains_vertex(n + 2));
    if triangles.is_empty() {
        return Err(EngineError::DegenerateMesh("Triangulation produced no valid triangles".into()));
    }
    Ok(triangles)
}

/// Insert point into mesh using Bowyer-Watson cavity formation.
fn insert_point_bw(triangles: &mut Vec<Triangle>, pts: &[Point2D], pt_idx: usize, p: Point2D) {
    let mut bad_triangles = Vec::new();
    for (i, tri) in triangles.iter().enumerate() {
        if in_circumcircle(pts[tri.vertices[0]], pts[tri.vertices[1]], pts[tri.vertices[2]], p) {
            bad_triangles.push(i);
        }
    }
    let boundary = extract_cavity_boundary(triangles, &bad_triangles);
    bad_triangles.sort_unstable();
    for &idx in bad_triangles.iter().rev() {
        triangles.swap_remove(idx);
    }
    for (u, v) in boundary {
        let tri = oriented_triangle(pts, u, v, pt_idx);
        triangles.push(tri);
    }
}

/// Build oriented triangle ensuring counter-clockwise ordering.
fn oriented_triangle(pts: &[Point2D], u: usize, v: usize, p: usize) -> Triangle {
    if orient2d(pts[u], pts[v], pts[p]) > 0.0 {
        Triangle::new(u, v, p)
    } else {
        Triangle::new(v, u, p)
    }
}

/// Extract boundary edges belonging to exactly one bad triangle.
fn extract_cavity_boundary(triangles: &[Triangle], bad_indices: &[usize]) -> Vec<(usize, usize)> {
    let mut edges = Vec::new();
    for &idx in bad_indices {
        for edge in triangles[idx].edges() {
            edges.push(edge);
        }
    }
    let mut boundary = Vec::new();
    for (i, &(u1, v1)) in edges.iter().enumerate() {
        let shared = edges.iter().enumerate().any(|(j, &(u2, v2))| {
            i != j && ((u1 == u2 && v1 == v2) || (u1 == v2 && v1 == u2))
        });
        if !shared {
            boundary.push((u1, v1));
        }
    }
    boundary
}

/// Circumcircle test: true if point P lies inside circumcircle of ABC (ordered CCW).
fn in_circumcircle(a: Point2D, b: Point2D, c: Point2D, p: Point2D) -> bool {
    let ccw = orient2d(a, b, c) > 0.0;
    let (v0, v1, v2) = if ccw { (a, b, c) } else { (a, c, b) };
    let d0x = v0.x - p.x;
    let d0y = v0.y - p.y;
    let d1x = v1.x - p.x;
    let d1y = v1.y - p.y;
    let d2x = v2.x - p.x;
    let d2y = v2.y - p.y;

    let d0_sq = d0x * d0x + d0y * d0y;
    let d1_sq = d1x * d1x + d1y * d1y;
    let d2_sq = d2x * d2x + d2y * d2y;

    let det = d0x * (d1y * d2_sq - d1_sq * d2y)
        - d0y * (d1x * d2_sq - d1_sq * d2x)
        + d0_sq * (d1x * d2y - d1y * d2x);
    det > 1e-11
}

/// Create super-triangle enclosing point cloud bounding box.
fn make_super_triangle(points: &[Point2D]) -> ([Point2D; 3], Triangle) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in points {
        min_x = min_x.min(p.x);
        max_x = max_x.max(p.x);
        min_y = min_y.min(p.y);
        max_y = max_y.max(p.y);
    }
    let dx = (max_x - min_x).max(1.0);
    let dy = (max_y - min_y).max(1.0);
    let dmax = dx.max(dy);
    let x_mid = (min_x + max_x) * 0.5;
    let y_mid = (min_y + max_y) * 0.5;

    let s0 = Point2D::new(x_mid - 25.0 * dmax, y_mid - dmax);
    let s1 = Point2D::new(x_mid, y_mid + 25.0 * dmax);
    let s2 = Point2D::new(x_mid + 25.0 * dmax, y_mid - dmax);
    let n = points.len();
    ([s0, s1, s2], Triangle::new(n, n + 1, n + 2))
}

/// Compute barycentric coordinates for point P inside triangle ABC.
fn compute_barycentric(a: Point2D, b: Point2D, c: Point2D, p: Point2D) -> Option<[f64; 3]> {
    let total_area = orient2d(a, b, c);
    if total_area.abs() < 1e-12 {
        return None;
    }
    let sa = orient2d(p, b, c);
    let sb = orient2d(a, p, c);
    let sc = orient2d(a, b, p);

    let tol = -1e-6 * total_area.abs();
    let is_ccw = total_area > 0.0;
    let valid = if is_ccw {
        sa >= tol && sb >= tol && sc >= tol
    } else {
        sa <= -tol && sb <= -tol && sc <= -tol
    };
    if !valid {
        return None;
    }
    let wa = (sa / total_area).max(0.0);
    let wb = (sb / total_area).max(0.0);
    let wc = (sc / total_area).max(0.0);
    let sum = wa + wb + wc;
    if sum < 1e-12 {
        return None;
    }
    Some([wa / sum, wb / sum, wc / sum])
}

/// Fallback inverse distance weighting (IDW) interpolation.
fn fallback_inverse_distance(points: &[Point2D], values: &[f64], p: Point2D) -> f64 {
    let mut weight_sum = 0.0;
    let mut val_sum = 0.0;
    for (pt, &val) in points.iter().zip(values) {
        let d = pt.dist_sq(&p).sqrt();
        if d < 1e-7 {
            return val;
        }
        let w = 1.0 / (d * d);
        weight_sum += w;
        val_sum += w * val;
    }
    if weight_sum > 0.0 {
        val_sum / weight_sum
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delaunay_regular_square() {
        let pts = vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(10.0, 0.0),
            Vector2::new(10.0, 10.0),
            Vector2::new(0.0, 10.0),
        ];
        let mesh = Delaunay2D::new(&pts).expect("triangulation failed");
        assert_eq!(mesh.triangles.len(), 2);
        let vals = vec![0.0, 10.0, 20.0, 10.0];
        let query = Vector2::new(5.0, 5.0);
        let interp = mesh.interpolate_barycentric(&query, &vals).expect("interpolation failed");
        assert!((interp - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_delaunay_collinear_error() {
        let pts = vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 1.0),
            Vector2::new(2.0, 2.0),
            Vector2::new(3.0, 3.0),
        ];
        let res = Delaunay2D::new(&pts);
        assert!(res.is_err());
    }

    #[test]
    fn test_delaunay_duplicate_points() {
        let pts = vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(10.0, 0.0),
            Vector2::new(10.0, 0.0), // duplicate
            Vector2::new(0.0, 10.0),
        ];
        let mesh = Delaunay2D::new(&pts).expect("triangulation with duplicates");
        assert_eq!(mesh.points.len(), 3);
        assert_eq!(mesh.triangles.len(), 1);
    }

    #[test]
    fn test_delaunay_insufficient_points() {
        let pts = vec![Vector2::new(0.0, 0.0), Vector2::new(1.0, 0.0)];
        assert!(Delaunay2D::new(&pts).is_err());
    }

    #[test]
    fn test_delaunay_out_of_hull_fallback() {
        let pts = vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(10.0, 0.0),
            Vector2::new(0.0, 10.0),
        ];
        let mesh = Delaunay2D::new(&pts).expect("triangulation failed");
        let vals = vec![1.0, 2.0, 3.0];
        let query = Vector2::new(50.0, 50.0);
        let interp = mesh.interpolate_barycentric(&query, &vals).expect("fallback");
        assert!(interp > 0.0 && interp < 4.0);
    }
}
