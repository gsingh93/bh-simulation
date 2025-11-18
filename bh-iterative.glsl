#iChannel0 "file:///home/gsingh2011/Downloads-chromeos/stars2.jpg"

// ================== parameters ==================
const float PI = 3.14159265359;
const int MAX_STEPS = 96;      // integration steps
const float R_FAR = 50.0;      // stop once ray is this far from BH
const float STEP_MIN = 0.01;   // min step length (world units along ray)
const float STEP_MAX = 0.5;    // max step length
const float STEP_SCALE = 0.10; // ds ≈ STEP_SCALE * r
const float CAPTURE_R = 1.05;  // capture if r < CAPTURE_R * rs
const float BEND_GAIN = 1.0;   // 1.0–2.0: calibrates weak-field strength

struct BH {
  vec3 pos;
  float rs;
}; // rs = Schwarzschild radius

// ================== camera utilities ==================
void cameraBasis(in vec3 ro, in vec3 ta, out vec3 f, out vec3 r, out vec3 u) {
  f = normalize(ta - ro);
  vec3 up = (abs(f.y) > 0.99) ? vec3(0, 0, 1) : vec3(0, 1, 0);
  r = normalize(cross(f, up));
  u = cross(r, f);
}

vec3 makeRay(vec2 p, vec3 ro, vec3 ta, float fov) {
  vec3 f, r, u;
  cameraBasis(ro, ta, f, r, u);
  float h = tan(fov * 0.5);
  return normalize(r * (p.x * h) + u * (p.y * h) + f);
}

// ================== spherical sky sampling ==================
vec3 sampleSkySpherical(vec3 dir) {
  dir = normalize(dir);
  float phi = atan(dir.z, dir.x);              // longitude
  float theta = asin(clamp(dir.y, -1.0, 1.0)); // latitude
  vec2 uv = vec2(phi / (2.0 * PI) + 0.5, 0.5 - theta / PI);
  uv.x = fract(uv.x);
  uv.y = clamp(uv.y, 0.0, 1.0);
  return texture(iChannel0, uv).rgb;
}

// ================== GR weak-field acceleration ==================
// dn/ds ≈ K * rs * (I - n n^T) * (R / r^3),   R = bh.pos - pos
// The transverse projector (I - n n^T) keeps changes orthogonal to n.
vec3 photonAccel(vec3 pos, vec3 n, BH bh) {
  vec3 R = bh.pos - pos;
  float r2 = dot(R, R);
  float r = sqrt(max(r2, 1e-12));
  vec3 grad = R / (r2 * r);                  // ~ R / r^3
  vec3 transverse = grad - n * dot(n, grad); // remove component along n
  return BEND_GAIN * bh.rs * transverse;
}

// ================== one RK4 step on (pos, n) ==================
void rk4Step(inout vec3 pos, inout vec3 n, BH bh, float ds) {
  // k1
  vec3 k1_p = n;
  vec3 k1_n = photonAccel(pos, n, bh);

  // k2
  vec3 n2 = normalize(n + 0.5 * ds * k1_n);
  vec3 k2_p = n2;
  vec3 k2_n = photonAccel(pos + 0.5 * ds * k1_p, n2, bh);

  // k3
  vec3 n3 = normalize(n + 0.5 * ds * k2_n);
  vec3 k3_p = n3;
  vec3 k3_n = photonAccel(pos + 0.5 * ds * k2_p, n3, bh);

  // k4
  vec3 n4 = normalize(n + ds * k3_n);
  vec3 k4_p = n4;
  vec3 k4_n = photonAccel(pos + ds * k3_p, n4, bh);

  pos += (ds / 6.0) * (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p);
  n = normalize(n + (ds / 6.0) * (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n));
}

// ================== iterative bending integrator ==================
vec3 bendRayIter(vec3 ro, vec3 rd, BH bh) {
  vec3 pos = ro;
  vec3 n = normalize(rd);

  // quick capture pre-test (optional but cheap)
  vec3 toBH0 = bh.pos - ro;
  float b0 = length(cross(toBH0, n)); // impact parameter of initial ray
  float bcrit = 2.598076211 * bh.rs;  // photon-sphere criterion
  if (dot(n, normalize(toBH0)) > 0.0 &&
      b0 < bcrit)     // aimed toward BH and sub-critical b
    return vec3(0.0); // captured

  for (int i = 0; i < MAX_STEPS; ++i) {
    vec3 R = bh.pos - pos;
    float r = length(R);

    // capture during integration
    if (r < CAPTURE_R * bh.rs)
      return vec3(0.0);

    // far enough that further bending is negligible
    if (r > R_FAR)
      break;

    // adaptive step: smaller near the BH, larger when far
    float ds = clamp(STEP_SCALE * r, STEP_MIN, STEP_MAX);

    rk4Step(pos, n, bh, ds);
  }

  return normalize(n);
}

// ================== render ==================
vec3 render(vec2 p) {
  vec3 ro = vec3(0.0, 0.0, -4.0);
  vec3 ta = vec3(0.0, 0.0, 0.0);
  float fov = radians(50.0);

  vec3 rd = makeRay(p, ro, ta, fov);

  BH bh;
  bh.pos = vec3(0.0);
  bh.rs = 0.02; // tweak rs

  vec3 rdb = bendRayIter(ro, rd, bh);
  bool captured = (length(rdb) < 1e-6);

  vec3 col = captured ? vec3(0.0) : sampleSkySpherical(rdb);

  // For LDR sky textures (PNG/JPG), skip tonemapping/gamma.
  // If using HDR (.hdr/.exr), enable:
  // col = col / (1.0 + col);
  // col = pow(col, vec3(1.0/2.2));

  return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;
  fragColor = vec4(render(p), 1.0);
}
