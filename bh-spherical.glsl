/*
#iChannel0 "file://./blue_nebulae_2.png"
#iChannel1 "file://./buf.glsl"
#iKeyboard
#iUniform float radius = 0.1 in{0.00, 0.4 }
#iUniform float fov = 50 in{10, 80 }
*/

uniform sampler2D iChannel0;
uniform sampler2D iChannel1;
uniform vec3 iResolution;
uniform vec4 iMouse;
uniform float radius;
uniform float fov;

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
/*
vec3 debugPattern(vec3 dir) {
  dir = normalize(dir);
  float phi = atan(dir.z, dir.x);              // [-π,π]
  float theta = asin(clamp(dir.y, -1.0, 1.0)); // [-π/2,π/2]

  // Normalize to [0,1]
  float u = phi / (2.0 * PI) + 0.5;
  float v = 0.5 - theta / PI;

  // Big colored quadrants for orientation
  // left/right: u < 0.5 vs u ≥ 0.5
  // top/bottom: v < 0.5 vs v ≥ 0.5
  vec3 base;
  if (u < 0.5 && v < 0.5)
    base = vec3(1, 0, 0); // red
  else if (u >= 0.5 && v < 0.5)
    base = vec3(0, 1, 0); // green
  else if (u < 0.5 && v >= 0.5)
    base = vec3(0, 0, 1); // blue
  else
    base = vec3(1, 1, 0); // yellow

  // Optional: add grid lines in u,v
  float gridU = step(0.95, fract(u * 8.0)); // 8 vertical bands
  float gridV = step(0.95, fract(v * 4.0)); // 4 horizontal bands
  float grid = max(gridU, gridV);

  // Blend grid lines in white
  return mix(base, vec3(1.0), grid * 0.6);
}

vec3 debugPattern2(vec3 dir) {
  dir = normalize(dir);

  float phi = atan(dir.z, dir.x);
  float theta = asin(clamp(dir.y, -1.0, 1.0));

  float u = phi / (2.0 * PI) + 0.5;
  float v = 0.5 - theta / PI;

  float ar = iResolution.y / iResolution.x;

  float cu = floor(u * 64.0 * ar);
  float cv = floor(v * 64.0);
  float idx = mod(cu + cv, 4.0);

  if (idx < 1.0)
    return vec3(1, 0, 0); // red
  if (idx < 2.0)
    return vec3(0, 1, 0); // green
  if (idx < 3.0)
    return vec3(0, 0, 1); // blue
  return vec3(1, 1, 0);   // yellow
}

/*float checkerAA(vec2 p) {
  vec2 q = sin(PI * p * vec2(20, 10));
  float m = q.x * q.y;
  return .5 - m / fwidth(m);
}*/
/* ---------------------------- Render ----------------------------- */
vec3 render(vec2 p) {
  // Camera
  vec3 ro;
  float mx = iMouse.x;           // mouse X in pixels
  float sx = mx / iResolution.x; // 0..1
  float yaw = sx * 2.0 * PI;     // 0..2π
  // float yaw = texture(iChannel1, vec2(0.5, 0.5)).r + 1.0;
  float R = 4.0;         // distance from BH
  vec3 center = vec3(0); // usually vec3(0.0)
  ro.x = center.x + R * cos(yaw);
  ro.y = center.y; // keep level; or add pitch later
  ro.z = center.z + R * sin(yaw);

  vec3 ta = center; // always look at the black hole
  // vec3 ro = vec3(0.0, 4.0, -7.0);
  // vec3 ta = vec3(0.0, 0.0, 0.0);
  // float fov = 50.0; // delete me
  float fov = radians(fov);

  // Primary ray
  vec3 rd = makeRay(p, ro, ta, fov);

  // Black hole (set rs=0.0 to disable lensing)
  BH bh;

  float mouseX = 1.0 * (iMouse.x / iResolution.x) * 2.0 - 1.0;
  float mouseY = 1.0 * (iMouse.y / iResolution.y) * 2.0 - 1.0;
  float mouseRs = 0.02 + 0.2 * (iMouse.y / iResolution.y);

  /*
    float radius = mouseRs;

    if (isKeyPressed(Key_RightArrow))
      radius += 0.01; // right arrow
    if (isKeyPressed(Key_LeftArrow))
      radius -= 0.01; // left arrow

    radius = clamp(radius, 0.01, 5.0);
    */

  // float radius = 0.1; // delete me
  //  bh.pos = vec3(-1.0 * mouseX, 0.0 /*mouseY*/, 0.0);
  bh.pos = vec3(0.0);
  bh.rs = radius;

  // Bend ray
  vec3 rdb = bendRay(ro, rd, bh);
  bool captured = (length(rdb) < 1e-5);

  // Sample spherical sky with bent direction
  // vec3 col = captured ? vec3(0.0) : sampleSkySpherical(rdb);
  vec3 col;
  if (captured) {
    col = vec3(0.0);
  } else {
    col = sampleSkySpherical(rdb);
    // float ty = 0.5 + 0.5 * rdb.y;
    // float tx = 0.5 + 0.5 * rdb.x;
    // col = vec3(0.0, tx, ty);
    // col = debugPattern2(rdb);
    // col = );
    // col = mix(vec3(0.1, 0.1, 0.1), vec3(0.0, 0.0, 1.0), t);
  }

  return col;
}

/* --------------------------- Entry point ------------------------- */
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // NDC p in [-1,1]^2 with aspect correction on X
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;

  fragColor = vec4(render(p), 1.0);
}
