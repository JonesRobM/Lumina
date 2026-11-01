// GPU matrix assembly for the free-space dyadic Green's function.
//
// Each thread computes one element (row, col) of the 3Nx3N interaction matrix.
// Diagonal blocks (i == j) are written as zero -- patched separately with alpha^-1.
// Off-diagonal blocks get -G_alphabeta(r_i, r_j) in f32 precision (~5x10^-^5 tolerance).
//
// Complex numbers stored as vec2<f32>: .x = real, .y = imag.

@group(0) @binding(0) var<storage, read>        positions : array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write>  matrix    : array<vec2<f32>>;

struct Params {
    dim       : u32,   // 3 * n_dipoles
    n_dipoles : u32,
    k         : f32,   // wavenumber (nm^-1), always real
    _pad      : f32,
}
@group(0) @binding(2) var<uniform> params : Params;

// Compute G_alpha_beta for displacement vector (rx, ry, rz) at wavenumber k.
fn greens_component(rx: f32, ry: f32, rz: f32, k: f32, alpha: u32, beta: u32) -> vec2<f32> {
    let r_sq = rx * rx + ry * ry + rz * rz;
    if r_sq < 1e-20 {
        return vec2<f32>(0.0, 0.0);
    }
    let r  = sqrt(r_sq);
    let kr = k * r;
    let kr_sq = kr * kr;

    // exp(ikr) = cos(kr) + i*sin(kr)
    let exp_re = cos(kr);
    let exp_im = sin(kr);

    // prefactor = k^2*exp(ikr) / (4pi*r)
    let scale  = k * k / (4.0 * 3.14159265358979 * r);
    let pf_re  = scale * exp_re;
    let pf_im  = scale * exp_im;

    // a = 1 + (ikr - 1) / kr^2  = (1 - 1/kr^2) + i/kr
    let a_re = 1.0 - 1.0 / kr_sq;
    let a_im = 1.0 / kr;

    // b = (3 - 3ikr - kr^2) / kr^2  = (3-kr^2)/kr^2 - 3i/kr
    let b_re = (3.0 - kr_sq) / kr_sq;
    let b_im = -3.0 / kr;

    // r_hat components
    var rhat : array<f32, 3>;
    rhat[0] = rx / r;
    rhat[1] = ry / r;
    rhat[2] = rz / r;

    let ri = rhat[alpha];
    let rj = rhat[beta];

    // delta_alphabeta
    var delta : f32 = 0.0;
    if alpha == beta { delta = 1.0; }

    // val = a*delta_alphabeta + b*r_hat_alpha*r_hat_beta
    let rr     = ri * rj;
    let val_re = a_re * delta + b_re * rr;
    let val_im = a_im * delta + b_im * rr;

    // result = (pf_re + i*pf_im) * (val_re + i*val_im)
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
        // Diagonal block: written as zero; patched with alpha^-1 on CPU after dispatch.
        matrix[row * dim + col] = vec2<f32>(0.0, 0.0);
        return;
    }

    let pi = positions[i];
    let pj = positions[j];
    let rx = pi.x - pj.x;
    let ry = pi.y - pj.y;
    let rz = pi.z - pj.z;

    // Interaction matrix: A_ij = alpha^-1 delta_ij - G_ij
    // Off-diagonal entry = -G_alphabeta(r_i, r_j)
    let g = greens_component(rx, ry, rz, params.k, alpha, beta);
    matrix[row * dim + col] = vec2<f32>(-g.x, -g.y);
}
