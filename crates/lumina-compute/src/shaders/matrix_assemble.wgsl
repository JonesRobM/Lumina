// GPU matrix assembly for the free-space dyadic Green's function.
//
// Each thread computes one element (row, col) of the 3N×3N interaction matrix.
// Diagonal blocks (i == j) are written as zero — patched separately with α⁻¹.
// Off-diagonal blocks get −G_αβ(r_i, r_j) in f32 precision (~5×10⁻⁵ tolerance).
//
// Complex numbers stored as vec2<f32>: .x = real, .y = imag.

@group(0) @binding(0) var<storage, read>        positions : array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write>  matrix    : array<vec2<f32>>;

struct Params {
    dim       : u32,   // 3 * n_dipoles
    n_dipoles : u32,
    k         : f32,   // wavenumber (nm⁻¹), always real
    _pad      : f32,
}
@group(0) @binding(2) var<uniform> params : Params;

// Complex multiply: (a + ib)(c + id)
fn cmul(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
    return vec2<f32>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

// Compute G_alpha_beta for displacement vector (rx, ry, rz) at wavenumber k.
fn greens_component(rx: f32, ry: f32, rz: f32, k: f32, alpha: u32, beta: u32) -> vec2<f32> {
    let r_sq = rx * rx + ry * ry + rz * rz;
    if r_sq < 1e-20 {
        return vec2<f32>(0.0, 0.0);
    }
    let r  = sqrt(r_sq);
    let kr = k * r;
    let kr_sq = kr * kr;

    // exp(ikr) = cos(kr) + i·sin(kr)
    let exp_re = cos(kr);
    let exp_im = sin(kr);

    // prefactor = k²·exp(ikr) / (4π·r)
    let scale  = k * k / (4.0 * 3.14159265358979 * r);
    let pf_re  = scale * exp_re;
    let pf_im  = scale * exp_im;

    // a = 1 + (ikr − 1) / kr²  = (1 − 1/kr²) + i/kr
    let a_re = 1.0 - 1.0 / kr_sq;
    let a_im = 1.0 / kr;

    // b = (3 − 3ikr − kr²) / kr²  = (3−kr²)/kr² − 3i/kr
    let b_re = (3.0 - kr_sq) / kr_sq;
    let b_im = -3.0 / kr;

    // r̂ components
    var rhat : array<f32, 3>;
    rhat[0] = rx / r;
    rhat[1] = ry / r;
    rhat[2] = rz / r;

    let ri = rhat[alpha];
    let rj = rhat[beta];

    // δ_αβ
    var delta : f32 = 0.0;
    if alpha == beta { delta = 1.0; }

    // val = a·δ_αβ + b·r̂_α·r̂_β
    let rr     = ri * rj;
    let val_re = a_re * delta + b_re * rr;
    let val_im = a_im * delta + b_im * rr;

    // result = (pf_re + i·pf_im) · (val_re + i·val_im)
    let res_re = pf_re * val_re - pf_im * val_im;
    let res_im = pf_re * val_im + pf_im * val_re;

    return vec2<f32>(res_re, res_im);
}

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let row = gid.x;
    let col = gid.y;
    let dim = params.dim;

    if row >= dim || col >= dim {
        return;
    }

    let i     = row / 3u;
    let j     = col / 3u;
    let alpha = row % 3u;
    let beta  = col % 3u;

    if i == j {
        // Diagonal block: written as zero; patched with α⁻¹ on CPU after dispatch.
        matrix[row * dim + col] = vec2<f32>(0.0, 0.0);
        return;
    }

    let pi = positions[i];
    let pj = positions[j];
    let rx = pi.x - pj.x;
    let ry = pi.y - pj.y;
    let rz = pi.z - pj.z;

    // Interaction matrix: A_ij = α⁻¹ δ_ij − G_ij
    // Off-diagonal entry = −G_αβ(r_i, r_j)
    let g = greens_component(rx, ry, rz, params.k, alpha, beta);
    matrix[row * dim + col] = vec2<f32>(-g.x, -g.y);
}
