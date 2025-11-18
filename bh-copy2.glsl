#iChannel0 "file:///home/gsingh2011/Downloads-chromeos/stars.jpg"
/* ============================================================
   Black-hole lensing over a flat 2D texture (billboard plane)
   - Works with any 2D image bound to iChannel0
   - Covers the screen frustum in BOTH axes
   - Safe fallbacks + less aggressive capture
   ============================================================ */

// Tunables
#define EPS 1e-6
#define CAPTURE_ON 1        // set to 0 to temporarily disable capture
#define CAPTURE_FACTOR 0.95 // b < CAPTURE_FACTOR * bCrit(rs)
const float PI = 3.14159265359;

/* ----------------------------- Types ----------------------------- */
struct BH {
  vec3 pos; // black hole center (world)
  float rs; // Schwarzschild radius (lensing strength)
};

/* --------------------------- Utilities --------------------------- */

/* Orthonormal camera basis */
void cameraBasis(in vec3 ro, in vec3 ta, out vec3 f, out vec3 r, out vec3 u) {
  f = normalize(ta - ro);
  // Fallback if f nearly parallel to world up
  vec3 up = abs(f.y) > 0.99 ? vec3(0, 0, 1) : vec3(0, 1, 0);
  r = normalize(cross(f, up));
  u = cross(r, f);
}

/* Pinhole primary ray for NDC p in [-1,1]^2 (aspect already applied) */
vec3 makeRay(vec2 p, vec3 ro, vec3 ta, float fov) {
  vec3 f, r, u;
  cameraBasis(ro, ta, f, r, u);
  float h = tan(fov * 0.5);
  return normalize(r * (p.x * h) + u * (p.y * h) + f);
}

/* Rodrigues rotation */
vec3 rotateAround(vec3 v, vec3 axis, float ang) {
  float c = cos(ang), s = sin(ang);
  return v * c + cross(axis, v) * s + axis * dot(axis, v) * (1.0 - c);
}

/* -------------------- Lensing (Schwarzschild) -------------------- */
float impactParameter(vec3 ro, vec3 rd, vec3 p) {
  return length(cross(p - ro, rd)); // rd assumed normalized
}

float deflectionAngle(float b, float rs) {
  float x = rs / max(b, 1e-6);
  return 2.0 * x + (15.0 * PI / 16.0) * x * x;
}

float bCrit(float rs) {
  return 2.598076211 * rs; // (3*sqrt(3)/2)*rs
}

/* Bend ray toward BH by α. Returns vec3(0) as a "captured" sentinel. */
vec3 bendRay(vec3 ro, vec3 rd, BH bh) {
  vec3 toBH = bh.pos - ro;
  float b = impactParameter(ro, rd, bh.pos);

#if CAPTURE_ON
  // Gentler capture threshold to avoid over-capture while tuning
  if (dot(rd, normalize(toBH)) > 0.0 && b < CAPTURE_FACTOR * bCrit(bh.rs))
    return vec3(0.0); // captured
#endif

  vec3 axis = cross(rd, toBH);
  float ax2 = dot(axis, axis);
  if (ax2 < 1e-16)
    return rd; // colinear: no well-defined bend

  float alpha = deflectionAngle(b, bh.rs);
  return normalize(rotateAround(rd, normalize(axis), alpha));
}

/* ---------------------- Billboard (flat image) ------------------- */

/* Fallback-friendly image aspect ratio:
   - Try textureSize(iChannel0, 0)
   - Else try iChannelResolution[0]
   - Else fall back to screen aspect */
float imageAspectFallback(float ARscreen) {
  float ARi = ARscreen;
  // Some environments don’t support textureSize; guard it.
  // If your environment supports it, this will work:
  ivec2 ts = textureSize(iChannel0, 0);
  if (ts.x > 0 && ts.y > 0) {
    ARi = float(ts.x) / float(ts.y);
  } else if (iChannelResolution[0].x > 0.0 && iChannelResolution[0].y > 0.0) {
    ARi = iChannelResolution[0].x / iChannelResolution[0].y;
  }
  return ARi;
}

/* Compute plane half-sizes (W,H) at distance D that COVER the screen frustum */
void coverPlaneSizes(float fov, float D, float ARscreen, float ARimage,
                     out float W, out float H) {
  float Hs = tan(fov * 0.5) * D;
  float Ws = Hs * ARscreen;
  float s = max(Hs, Ws / max(ARimage, 1e-6)); // guard
  s *= 5.0;                                   // small padding
  H = s;
  W = s * ARimage;
}

/* Ray-plane intersection: plane (C, n). Returns false if no hit. */
bool rayPlane(vec3 ro, vec3 rd, vec3 C, vec3 n, out float tHit) {
  float denom = dot(rd, n);
  if (denom <= EPS)
    return false; // looking away/parallel
  float t = dot(C - ro, n) / denom;
  if (t <= 0.0)
    return false; // behind eye
  tHit = t;
  return true;
}

/* Map a point P on the plane to [0,1] UV within rectangle [-W..W]×[-H..H] */
vec2 planePointToUV(vec3 P, vec3 C, vec3 r, vec3 u, float W, float H) {
  vec3 L = P - C;
  float x = dot(L, r);                                   // in [-W..W]
  float y = dot(L, u);                                   // in [-H..H]
  return vec2(x / (2.0 * W) + 0.5, 0.5 - y / (2.0 * H)); // flip V
}

// -------------------------------------------------------
// Map a direction to screen-space UV using the camera pinhole
// Vertical FOV = fov, screen aspect = ARs
// Returns (uv, ok); ok=false if the ray points behind the camera.
// -------------------------------------------------------
bool directionToUV(vec3 dir, vec3 r, vec3 u, vec3 f, float fov, float ARs,
                   out vec2 uv) {
  // dir in camera basis
  float dx = dot(dir, r);
  float dy = dot(dir, u);
  float dz = dot(dir, f);
  if (dz <= 1e-6)
    return false; // pointing away

  // perspective divide onto plane z=1 in camera space
  float px = dx / dz;
  float py = dy / dz;

  // normalize by plane half-sizes: tan(fov/2) vertically, *ARs horizontally
  float h = tan(fov * 0.5);
  float x_ndc = px / (h * ARs); // [-1,1] when inside FOV
  float y_ndc = py / h;

  // NDC -> UV
  uv = vec2(x_ndc * 0.5 + 0.5, 0.5 - y_ndc * 0.5);
  return true;
}

// -------------------------------------------------------
// Render: flat 2D texture projected from camera, with GR bending
// -------------------------------------------------------
vec3 render(vec2 p) {
  // Camera
  vec3 ro = vec3(0.0, 0.0, -4.0);
  vec3 ta = vec3(0.0, 0.0, 0.0);
  float fov = radians(50.0);
  float ARs = iResolution.x / iResolution.y;

  // Camera basis
  vec3 f, r, u;
  cameraBasis(ro, ta, f, r, u);

  // Primary ray and bending
  vec3 rd = makeRay(p, ro, ta, fov);
  BH bh;
  bh.pos = vec3(0.0);
  bh.rs = 0.1; // tweak rs
  vec3 rdb = bendRay(ro, rd, bh);
  bool captured = (length(rdb) < 1e-5);

  // Project bent direction to UV (no finite plane)
  vec2 uv;
  bool ok = directionToUV(rdb, r, u, f, fov, ARs, uv);

  // If the bent ray leaves the FOV, show a neutral gray instead of black
  bool inside =
      ok && (uv.x >= 0.0 && uv.x <= 1.0 && uv.y >= 0.0 && uv.y <= 1.0);
  vec3 img = inside ? texture(iChannel0, uv).rgb : vec3(0.2);

  // Compose: captured rays = black hole; else sample image
  vec3 col = captured ? vec3(0.0) : img;

  return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;
  fragColor = vec4(render(p), 1.0);
}
