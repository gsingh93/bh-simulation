#iChannel0 "file:///home/gsingh2011/code/bh-simulation/blue_nebulae_2.png"

/* ============================================================
   Black-hole lensing with a spherical (equirectangular) sky map
   - Bind an equirectangular sky texture to iChannel0 (U wraps)
   - No billboard/cover plane; sample by bent ray direction
   ============================================================ */

const float PI = 3.14159265359;

/* ----------------------------- Types ----------------------------- */
struct BH {
  vec3 pos; // black hole position (world)
  float rs; // Schwarzschild radius (controls lensing strength)
};

/* --------------------------- Utilities --------------------------- */
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
  return 2.0 * x + (15.0 * PI / 16.0) * x * x; // weak-field approx (radians)
}

float bCrit(float rs) {
  return 2.598076211 * rs; // (3*sqrt(3)/2)*rs (photon sphere)
}

vec3 bendRay(vec3 ro, vec3 rd, BH bh) {
  vec3 toBH = bh.pos - ro;
  float b = impactParameter(ro, rd, bh.pos);

  // Capture: heading toward BH and below photon-sphere threshold
  if (dot(rd, normalize(toBH)) > 0.0 && b < bCrit(bh.rs))
    return vec3(0.0); // sentinel: captured

  // Stable bend axis (avoid degeneracy when cross≈0)
  vec3 axis = cross(rd, toBH);
  if (dot(axis, axis) < 1e-16)
    axis =
        normalize(cross(rd, abs(rd.y) < 0.9 ? vec3(0, 1, 0) : vec3(1, 0, 0)));
  else
    axis = normalize(axis);

  float alpha = deflectionAngle(b, bh.rs);
  return normalize(rotateAround(rd, axis, alpha));
}

/* --------------------- Spherical sky sampling -------------------- */
/* Equirectangular (lat-long) mapping:
   - U wraps longitude φ in (−π..π]
   - V maps latitude θ in [−π/2..+π/2] to [1..0] (V=0 is up) */
vec3 sampleSkySpherical(vec3 dir) {
  dir = normalize(dir);
  float phi = atan(dir.z, dir.x);              // longitude
  float theta = asin(clamp(dir.y, -1.0, 1.0)); // latitude
  vec2 uv = vec2(phi / (2.0 * PI) + 0.5, 0.5 - theta / PI);
  uv.x = fract(uv.x);           // wrap U
  uv.y = clamp(uv.y, 0.0, 1.0); // clamp V
  return texture(iChannel0, uv).rgb;
}

/* ---------------------------- Render ----------------------------- */
vec3 render(vec2 p) {
  // Camera
  vec3 ro = vec3(0.0, 4.0, -7.0);
  vec3 ta = vec3(0.0, 0.0, 0.0);
  float fov = radians(50.0);

  // Primary ray
  vec3 rd = makeRay(p, ro, ta, fov);

  // Black hole (set rs=0.0 to disable lensing)
  BH bh;
  float mouseX = 3.0 * (iMouse.x / iResolution.x) * 2.0 - 1.0;
  float mouseY = 2.0 * (iMouse.y / iResolution.y) * 2.0 - 1.0;
  float mouseRs = 0.02 + 0.1 * (iMouse.y / iResolution.y);
  bh.pos = vec3(-1.0 * mouseX, 0.0 /*mouseY*/, 0.0);
  bh.rs = mouseRs;

  // Bend ray
  vec3 rdb = bendRay(ro, rd, bh);
  bool captured = (length(rdb) < 1e-5);

  // Sample spherical sky with bent direction
  vec3 col = captured ? vec3(0.0) : sampleSkySpherical(rdb);

  return col;
}

/* --------------------------- Entry point ------------------------- */
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // NDC p in [-1,1]^2 with aspect correction on X
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;

  fragColor = vec4(render(p), 1.0);
}
