#define vscode 0

#if vscode == 0
// WebGL / browser uniforms
uniform sampler2D iChannel0; // environment (equirectangular) texture
uniform sampler2D iChannel1; // yaw buffer (if used in your host)
uniform vec3 iResolution;    // (width, height, depth)
uniform vec4 iMouse;         // xy = current, zw = click
uniform float radius;        // black hole Schwarzschild radius (world units)
uniform float fov;           // vertical field of view in degrees
#else
// Shadertoy-style setup for local testing in VSCode
#iChannel0 "file://./blue_nebulae_2.png"
#iChannel1 "file://./buf.glsl"
#iKeyboard
#iUniform float radius = 0.1 in{0.00, 0.4 }
#iUniform float fov = 50.0 in{10, 80 }
#endif

// ----------------------------------------------------------
// Black-hole lensing with a spherical (equirectangular) sky
// ----------------------------------------------------------

const float PI = 3.14159265359;

// Black hole parameters (single BH)
struct BH {
  vec3 pos; // position in world space
  float rs; // Schwarzschild radius (controls lensing strength)
};

// ------------------------ Camera -------------------------

// Build an orthonormal camera basis (forward, right, up)
void cameraBasis(in vec3 ro, in vec3 ta, out vec3 f, out vec3 r, out vec3 u) {
  f = normalize(ta - ro);
  vec3 up = (abs(f.y) > 0.99) ? vec3(0, 0, 1) : vec3(0, 1, 0);
  r = normalize(cross(f, up));
  u = cross(r, f);
}

// Construct a ray direction from NDC p in [-1,1]^2
vec3 makeRay(vec2 p, vec3 ro, vec3 ta, float fovRad) {
  vec3 f, r, u;
  cameraBasis(ro, ta, f, r, u);
  float h = tan(fovRad * 0.5);
  return normalize(r * (p.x * h) + u * (p.y * h) + f);
}

// -------------------- Lensing (single kick) --------------------

// Impact parameter b: shortest distance from ray to BH center
float impactParameter(vec3 ro, vec3 rd, vec3 p) {
  // rd is assumed normalized
  return length(cross(p - ro, rd));
}

// Weak-field deflection angle α(b, rs) in radians
float deflectionAngle(float b, float rs) {
  float x = rs / max(b, 1e-6);
  // Leading GR term + next-order correction
  return 2.0 * x + (15.0 * PI / 16.0) * x * x;
}

// Critical impact parameter at photon sphere
float bCrit(float rs) {
  // (3*sqrt(3)/2)*rs
  return 2.598076211 * rs;
}

// Apply a single deflection at closest approach (thin-lens approximation)
vec3 bendRay(vec3 ro, vec3 rd, BH bh) {
  vec3 toBH = bh.pos - ro;
  float b = impactParameter(ro, rd, bh.pos);

  // Capture: ray is heading toward BH and b is inside photon-sphere threshold
  if (dot(rd, normalize(toBH)) > 0.0 && b < bCrit(bh.rs))
    return vec3(0.0); // sentinel: captured

  // Stable rotation axis; avoid degeneracy when rd ‖ toBH
  vec3 axis = cross(rd, toBH);
  if (dot(axis, axis) < 1e-16)
    axis =
        normalize(cross(rd, abs(rd.y) < 0.9 ? vec3(0, 1, 0) : vec3(1, 0, 0)));
  else
    axis = normalize(axis);

  float alpha = deflectionAngle(b, bh.rs);
  return normalize(rotateAround(rd, axis, alpha));
}

// ----------------- Environment sampling (equirect) -----------------

// Equirectangular (lat-long) mapping:
//   - dir normalized
//   - U wraps longitude φ in (−π..π]
//   - V maps latitude θ in [−π/2..+π/2] to [0..1]
vec3 sampleSkySpherical(vec3 dir) {
  dir = normalize(dir);
  float phi = atan(dir.z, dir.x);              // longitude
  float theta = asin(clamp(dir.y, -1.0, 1.0)); // latitude
  vec2 uv = vec2(
      phi / (2.0 * PI) + 0.5, // U: [-π,π] → [0,1]
      0.5 + theta /
                PI // V: [-π/2,π/2] → [0,1] (this variant expects "up" this way)
  );
  uv.x = fract(uv.x);           // wrap U
  uv.y = clamp(uv.y, 0.0, 1.0); // clamp V
  return texture(iChannel0, uv).rgb;
}

// ---------------------------- Render -----------------------------

vec3 render(vec2 p) {
  // Camera orbit around BH in the XZ-plane
  vec3 ro;
  float mx = iMouse.x;
  float sx = mx / iResolution.x; // [0,1]

#if vscode == 0
  // In browser / host mode: yaw directly from mouse X
  float yaw = sx * 2.0 * PI;
#else
  // In Shadertoy/VSCode mode: yaw stored in buffer (e.g., via a separate pass)
  float yaw = texture(iChannel1, vec2(0.5, 0.5)).r + 1.0;
#endif

  float R = 4.0;         // camera distance from BH
  vec3 center = vec3(0); // BH at origin

  ro.x = center.x + R * cos(yaw);
  ro.y = center.y; // fixed height; pitch can be added later if desired
  ro.z = center.z + R * sin(yaw);

  vec3 ta = center; // look at BH

  // Convert FOV uniform (degrees) to radians for ray construction
  float fovRad = radians(fov);

  // Primary ray in world space
  vec3 rd = makeRay(p, ro, ta, fovRad);

  // Black hole parameters
  BH bh;
  bh.pos = vec3(0.0);
  bh.rs = radius;

  // Single-kick bending
  vec3 rdb = bendRay(ro, rd, bh);
  bool captured = (length(rdb) < 1e-5);

  // Sample environment along bent direction; black hole appears as a dark disk
  vec3 col;
  if (captured) {
    col = vec3(0.0);
  } else {
    col = sampleSkySpherical(rdb);
  }

  return col;
}

// --------------------------- Entry point -------------------------

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // NDC p in [-1,1]^2 with aspect correction on X
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;

  fragColor = vec4(render(p), 1.0);
}
