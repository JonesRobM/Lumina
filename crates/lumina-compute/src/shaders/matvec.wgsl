// Complex matrix-vector product: y = A * x
//
// Complex numbers are stored as vec2<f32> where .x = real, .y = imag.
// The matrix is stored in row-major order as a flat array of vec2<f32>.
// Each workgroup thread computes one row of the output vector.

@group(0) @binding(0) var<storage, read> matrix: array<vec2<f32>>;
@group(0) @binding(1) var<storage, read> input_vec: array<vec2<f32>>;
@group(0) @binding(2) var<storage, read_write> output_vec: array<vec2<f32>>;

struct Params {
    dim: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}
@group(0) @binding(3) var<uniform> params: Params;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let row = gid.x;
    let dim = params.dim;
    if (row >= dim) {
        return;
    }

    var sum = vec2<f32>(0.0, 0.0);
    let row_offset = row * dim;

    for (var col = 0u; col < dim; col++) {
        let a = matrix[row_offset + col];
        let b = input_vec[col];
        // Complex multiply: (a.x + i*a.y) * (b.x + i*b.y)
        //   = (a.x*b.x - a.y*b.y) + i*(a.x*b.y + a.y*b.x)
        sum += vec2<f32>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
    }

    output_vec[row] = sum;
}
