#iChannel0 "file:///home/gsingh2011/Downloads-chromeos/stars.png"

/* ============================================================
   Black-hole lensing over a flat 2D texture on a billboard
   - Uses a plane that covers the screen frustum at distance D
   - Works with any regular 2D image in iChannel0
   - Clean, modular, and commented
   ============================================================ */

// Shadertoy: iChannel0 bound to any flat 2D PNG/JPG
// Example: stars, a photo, etc. (NOT an equirectangular pano)

#define EPS 1e-6
const float PI = 3.14159265359;

/* ----------------------------- Types ----------------------------- */
struct BH {
  vec3 pos; // black hole center (world)
  float rs; // Schwarzschild radius (controls lens strength)
};

/* --------------------------- Utilities --------------------------- */
// float saturate(float x) { return clamp(x, 0.0, 1.0); }

/* Orthonormal camera basis from eye (ro) to target (ta) */
void cameraBasis(in vec3 ro, in vec3 ta, out vec3 f, out vec3 r, out vec3 u) {
  f = normalize(ta - ro);                       // forward
  r = normalize(cross(f, vec3(0.0, 1.0, 0.0))); // right (use a fixed world-up)
  u = cross(r, f);                              // up
}

/* Pinhole camera primary ray dir for NDC p in [-1,1]^2 (with aspect applied) */
vec3 makeRay(vec2 p, vec3 ro, vec3 ta, float fov) {
  vec3 f, r, u;
  cameraBasis(ro, ta, f, r, u);
  float h = tan(fov * 0.5); // half-height of image plane @ unit distance
  return normalize(r * (p.x * h) + u * (p.y * h) + f);
}

/* Rodrigues rotation of vector v around 'axis' by angle 'ang' (radians) */
vec3 rotateAround(vec3 v, vec3 axis, float ang) {
  float c = cos(ang), s = sin(ang);
  return v * c + cross(axis, v) * s + axis * dot(axis, v) * (1.0 - c);
}

/* -------------------- Lensing (Schwarzschild) -------------------- */
float impactParameter(vec3 ro, vec3 rd, vec3 p) {
  // shortest distance between infinite line (ro + t*rd) and point p
  return length(cross(p - ro, rd)); // rd assumed normalized
}

float deflectionAngle(float b, float rs) {
  // α(b) ≈ 2*rs/b + (15π/16)*(rs/b)^2 (radians); guard b
  float x = rs / max(b, 1e-6);
  return 2.0 * x + (15.0 * PI / 16.0) * x * x;
}

float bCrit(float rs) {
  // photon-sphere impact parameter: (3√3/2)*rs
  return 2.598076211 * rs;
}

/* Bend a ray toward BH by α around the plane normal to (rd, toBH) */
vec3 bendRay(vec3 ro, vec3 rd, BH bh) {
  vec3 toBH = bh.pos - ro;
  float b = impactParameter(ro, rd, bh.pos);

  // Capture if heading toward BH and below photon-sphere threshold
  if (dot(rd, normalize(toBH)) > 0.0 && b < bCrit(bh.rs))
    return vec3(0.0); // sentinel: "captured"

  vec3 axis = cross(rd, toBH);
  float ax2 = dot(axis, axis);
  if (ax2 < 1e-16)
    return rd; // colinear (no well-defined bend)

  float alpha = deflectionAngle(b, bh.rs);
  return normalize(rotateAround(rd, normalize(axis), alpha));
}

/* ---------------------- Billboard (flat image) ------------------- */
/* Compute plane half-sizes (W,H) at distance D that COVER the screen frustum,
   preserving the image aspect ARi. Adds a small padding. */
void coverPlaneSizes(float fov, float D, float ARscreen, float ARimage,
                     out float W, out float H) {
  float Hs = tan(fov * 0.5) * D; // screen half-height at distance D
  float Ws = Hs * ARscreen;      // screen half-width  at distance D

  // scale so image rectangle (W,H) fully covers (Ws,Hs) without stretching
  float s = max(Hs, Ws / ARimage);

  // small pad so bent rays near edges still hit
  s *= 1.02;

  H = s;
  W = s * ARimage;
}

/* Ray-plane intersection: plane point C, normal n (assume normalized-ish) */
bool rayPlane(vec3 ro, vec3 rd, vec3 C, vec3 n, out float tHit) {
  float denom = dot(rd, n);
  if (denom <= EPS)
    return false; // looking away or parallel
  float t = dot(C - ro, n) / denom;
  if (t <= 0.0)
    return false; // behind eye
  tHit = t;
  return true;
}

/* Map a point P on the plane to [0,1] UV in the billboard rectangle
 * [-W..W]×[-H..H] */
vec2 planePointToUV(vec3 P, vec3 C, vec3 r, vec3 u, float W, float H) {
  vec3 L = P - C;
  float x = dot(L, r); // in [-W..W]
  float y = dot(L, u); // in [-H..H]
  // Flip V so "up" is v=0 (conventional textures)
  return vec2(x / (2.0 * W) + 0.5, 0.5 - y / (2.0 * H));
}

/* ---------------------------- Render ----------------------------- */
vec3 render(vec2 p) {
  // Camera
  vec3 ro = vec3(0.0, 0.0, -4.0);
  vec3 ta = vec3(0.0, 0.0, 0.0);
  float fov = radians(50.0);

  vec3 rd = makeRay(p, ro, ta, fov);

  // Camera basis for plane placement and UV mapping
  vec3 f, r, u;
  cameraBasis(ro, ta, f, r, u);

  // Black hole (set rs=0.0 to disable lensing)
  BH bh;
  bh.pos = vec3(0.0);
  bh.rs = 0.1;

  // Bend primary ray
  vec3 rdb = bendRay(ro, rd, bh);
  bool captured = (length(rdb) < 1e-5);

  // Billboard plane at distance D, sized to cover the screen in BOTH axes
  float D = 3.0;
  vec3 C = ro + f * D; // plane center
  vec3 n = f;          // plane normal (faces camera)

  // Screen & image aspects
  float ARs = iResolution.x / iResolution.y;
  ivec2 ts = textureSize(iChannel0, 0);
  float ARi = float(ts.x) / float(ts.y);

  float W, H;
  coverPlaneSizes(fov, D, ARs, ARi, W, H);

  // Intersect bent ray with plane
  float tHit;
  if (!rayPlane(ro, rdb, C, n, tHit))
    return vec3(0.0);

  vec3 P = ro + rdb * tHit; // hit point on plane
  vec2 uv = planePointToUV(P, C, r, u, W, H);

  // Avoid clamping (which smears). Mask outside instead.
  bool inside = (uv.x >= 0.0 && uv.x <= 1.0 && uv.y >= 0.0 && uv.y <= 1.0);

  // captured = false;
  // inside = false;
  vec3 col =
      captured ? vec3(0.0) : (inside ? texture(iChannel0, uv).rgb : vec3(0.0));

  // col = texture(iChannel0, uv).rgb;
  // col = vec3(uv.x * 20.0, uv.y, 1.0);
  //  Simple tonemap + gamma
  col = col / (1.0 + col);
  col = pow(col, vec3(1.0 / 2.2));
  return col;
}

/* --------------------------- Entry point ------------------------- */
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // NDC in [-1,1]^2 with aspect correction
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;

  vec3 col = render(p);
  fragColor = vec4(col, 1.0);
}
