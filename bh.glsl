#iChannel0 "file:///home/gsingh2011/Downloads-chromeos/stars.png"

const float PI = 3.14;

struct BH {
  vec3 pos;
  float rs; // Fix radius instead of calculating based on mass (2GM/c^2)
};

float impactParameter(vec3 ro, vec3 rd, vec3 p) {
  return length(cross(p - ro, rd));
}

float deflectionAngle(float b, float rs) {
  float x = rs / max(b, 1e-6);
  return 2.0 * x + (15.0 * PI / 16.0) * x * x;
}

float bCrit(float rs) { return 2.598076211 * rs; }

vec3 rotateAround(vec3 v, vec3 axis, float ang) {
  float c = cos(ang), s = sin(ang);
  return v * c + cross(axis, v) * s + axis * dot(axis, v) * (1.0 - c);
}

vec3 bendRay(vec3 ro, vec3 rd, BH bh) {
  vec3 toBH = bh.pos - ro;
  float b = impactParameter(ro, rd, bh.pos);
  if (dot(rd, normalize(toBH)) > 0.0 && b < 0.8 * bCrit(bh.rs))
    return vec3(0.0); // captured
  vec3 axis = cross(rd, toBH);
  if (dot(axis, axis) < 1e-16)
    return rd;
  float alpha = deflectionAngle(b, bh.rs);
  // alpha *= 2.0;
  return normalize(rotateAround(rd, normalize(axis), alpha));
}

vec3 makeRay(vec2 p, vec3 ro, vec3 ta, float fov) {
  vec3 f = normalize(ta - ro);
  vec3 r = normalize(cross(f, vec3(0, 1, 0)));
  vec3 u = cross(r, f);

  float h = tan(fov * 0.5) * 15.0;
  vec3 rd = normalize(r * (p.x * h) + u * (p.y * h) + f);
  return rd;
}

vec3 sampleSky(vec3 rd) {
  // Simple sky gradient driven by rd.y
  vec3 skyTop = vec3(0.1, 0.1, 0.1);
  vec3 skyBot = vec3(0.00, 0.00, 0.00);
  float t = clamp(0.5 + 0.5 * rd.y, 0.0, 1.0); // rd.y in [-1,1] → [0,1]
  vec3 col = mix(skyBot, skyTop, t);
  return col;
}

vec3 render(vec2 p, vec2 fragCoord) {
  vec3 ro = vec3(0, 0, -2), ta = vec3(0, 0, 0);
  float fov = radians(50.0);
  vec3 rd = makeRay(p, ro, ta, fov);

  BH bh;
  bh.pos = vec3(0);
  bh.rs = 1.1;

  vec3 rdb = bendRay(ro, rd, bh);
  bool captured = (length(rdb) < 1e-5);

  // vec4 texel = texture(iChannel0, fragCoord / iResolution.xy);

  // assumes an equirectangular env in iChannel0 (sampler2D)
  // +Z is forward, +X right, +Y up
  const float PI = 3.14159265359;

  vec3 d = normalize(rdb);                   // bent ray direction
  float phi = atan(d.z, d.x);                // azimuth in [-pi,pi]
  float theta = asin(clamp(d.y, -1.0, 1.0)); // elevation in [-pi/2,pi/2]

  // map to [0,1]^2; flip V so up is at v=0
  vec2 uv = vec2(phi / (2.0 * PI) + 0.5, 0.5 - theta / PI);

  // wrap U to avoid seam issues
  uv.x = fract(uv.x);

  // sample
  vec4 texel = texture(iChannel0, uv);

  vec3 col = captured ? vec3(0.0) : texel.rgb; // sampleSky(rdb);
  float b = impactParameter(ro, rd, bh.pos);
  float dim = smoothstep(bh.rs * 1.2, bh.rs * 3.5, b);
  col *= mix(0.35, 1.0, dim);

  col = col / (1.0 + col);
  col = pow(col, vec3(1.0 / 2.2));
  float bVis = smoothstep(0.0, 4.0, b);
  // col = vec3(bVis);
  // col = rdb;
  float alpha = deflectionAngle(b, bh.rs);
  // col = vec3(smoothstep(0.0, 4.0, alpha));
  // col = vec3(alpha);
  return col;
}

/*
vec3 render(vec2 p, vec2 fragCoord) {
  // --- Camera ---
  vec3 ro = vec3(0, 0, -4), ta = vec3(0, 0, 0);
  float fov = radians(50.0);
  vec3 rd = makeRay(p, ro, ta, fov);

  // Camera basis (to place the billboard plane)
  vec3 f = normalize(ta - ro);
  vec3 r = normalize(cross(f, vec3(0, 1, 0)));
  vec3 u = cross(r, f);

  // --- Black hole ---
  BH bh;
  bh.pos = vec3(0);
  bh.rs = 1.1; // tweak rs to taste

  // Bend primary ray
  vec3 rdb = bendRay(ro, rd, bh);
  bool captured = (length(rdb) < 1e-5);

  // --- Billboard plane matching the FOV ---
  // Plane located at distance D in front of camera, normal = f
  float D = 3.0;
  vec3 C = ro + f * D; // plane center
  vec3 n = f;          // plane normal (faces camera)

  // float pad = 1.02;

  // Frustum half-height at distance D
  float Hscreen = tan(fov * 0.5) * D;
  float Wscreen = Hscreen * (iResolution.x / iResolution.y);

  // Choose image half-sizes to COVER the screen without stretching:
  float Ai =
      float(textureSize(iChannel0, 0).x) / float(textureSize(iChannel0, 0).y);

  // Find the minimal scale s so that (Wimg>=Wscreen) and (Himg>=Hscreen)
  float s = max(Hscreen, Wscreen / Ai);

  // Optional padding so bent rays still hit the plane
  s *= 1.05;

  // Final plane half-sizes that preserve image aspect and cover the screen
  float H = s;
  float W = s * Ai;

  // float H = tan(fov * 0.5) * D; // half-height on plane
  float AR = iResolution.x / iResolution.y;
  // float W = H * AR; // half-width

  // Ray-plane intersection
  float denom = dot(rdb, n);
  if (denom <= 1e-6) {
    // Ray not hitting the forward-facing plane; show black (or background
    // color)
    return vec3(0.0);
  }
  float t = dot(C - ro, n) / denom;
  t = max(t, 0.0);
  vec3 P = ro + rdb * t;

  // Map hit point on plane to [0,1] UVs for the flat 2D texture
  vec3 L = P - C;
  float x = dot(L, r); // in [-W..W]
  float y = dot(L, u); // in [-H..H]
  vec2 uv =
      vec2(x / (2.0 * W) + 0.5, 0.5 - y / (2.0 * H)); // flip V so up is v=0

  /*bvec2 in01 = bvec2(uv.x >= 0.0 && uv.x <= 1.0, uv.y >= 0.0 && uv.y <= 1.0);
  bool inside = all(in01);

  vec3 col =
      captured ? vec3(0.0) : (inside ? texture(iChannel0, uv).rgb : vec3(0.0));

  // Sample the texture; clamp to avoid sampling outside
  vec3 col = captured ? vec3(0.0) : texture(iChannel0, clamp(uv, 0.0, 1.0)).rgb;

  // (Optional) simple tonemap + gamma
  col = col / (1.0 + col);
  col = pow(col, vec3(1.0 / 2.2));
  return col;
  }


  // Keep your existing makeRay(...) and bendRay(...)
  bool hitSphere(vec3 ro, vec3 rd, vec3 c, float r, out float tHit) {
    vec3 oc = ro - c;
    float b = dot(oc, rd);
    float c2 = dot(oc, oc) - r * r;
    float h = b * b - c2;
    if (h < 0.0)
      return false;
    h = sqrt(h);
    float t0 = -b - h;
    float t1 = -b + h;
    tHit = (t0 > 0.0) ? t0 : ((t1 > 0.0) ? t1 : -1.0);
    return tHit > 0.0;
  }

  /*
  vec3 render(vec2 p, vec2 fragCoord) {
    // --- Camera ---
    vec3 ro = vec3(0, 0, -4);
    vec3 ta = vec3(0, 0, 0);
    float fov = radians(50.0);

    // Primary ray
    vec3 rd = makeRay(p, ro, ta, fov);

    // Camera basis (for billboard plane)
    vec3 f = normalize(ta - ro);
    vec3 r = normalize(cross(f, vec3(0, 1, 0)));
    vec3 u = cross(r, f);

    // --- Black hole ---
    BH bh;
    bh.pos = vec3(0);
    bh.rs = 0.0; // set to 0.0 to disable lensing

    // Bend ray
    vec3 rdb = bendRay(ro, rd, bh);
    bool captured = (length(rdb) < 1e-5);

    // --- Billboard that COVERS the screen frustum at distance D ---
    float D = 3.0;
    float ARs = iResolution.x / iResolution.y;
    float Hs = tan(fov * 0.5) * D; // screen half-height at D
    float Ws = Hs * ARs;           // screen half-width  at D

    ivec2 ts = textureSize(iChannel0, 0);
    float ARi = float(ts.x) / float(ts.y); // image aspect

    // Scale s so image (W,H) covers (Ws,Hs) without stretching; small pad
    float s = max(Hs, Ws / ARi) * 1.02;
    float H = s;
    float W = s * ARi;

    // Plane center/normal
    vec3 C = ro + f * D;
    vec3 n = f;

    // Ray-plane intersection
    float denom = dot(rdb, n);
    if (denom <= 1e-6)
      return vec3(0.0);
    float t = dot(C - ro, n) / denom;
    if (t <= 0.0)
      return vec3(0.0);

    // Build bend plane + angle again (you already have these pieces)
    vec3 toBH = bh.pos - ro;
    vec3 axis = cross(rd, toBH);
    float axisLen2 = dot(axis, axis);
    float alpha = (axisLen2 < 1e-16)
                      ? 0.0
                      : deflectionAngle(impactParameter(ro, rd, bh.pos), bh.rs);
    axis = (axisLen2 < 1e-16) ? vec3(1, 0, 0) : normalize(axis);

    // Choose where the debug path ends: at the billboard hit distance
    float T = t;                           // path length (to your plane)
    int N = 24;                            // spheres along the path
    float rad = 0.02 * tan(fov * 0.5) * D; // thickness scales with view

    // Accumulate a bright overlay where the camera ray hits any debug sphere
    vec3 overlay = vec3(0.0);
    for (int i = 0; i < N; ++i) {
      float s = float(i) / float(N - 1); // 0..1 along the path
      // Rotate direction progressively by s*alpha (circular-arc approximation)
      vec3 dirS = normalize(rotateAround(rd, axis, s * alpha));
      vec3 posS = ro + dirS * (s * T); // from camera to plane distance
      float tHit;
      if (hitSphere(ro, rd, posS, rad, tHit)) {
        // Add a color (e.g., cyan) that doesn’t tonemap away too much
        overlay += vec3(0.1, 0.9, 1.0);
      }
    }

    // ... compute your usual 'col' for the scene ...
    // then blend the overlay on top (additive is fine for a debug line)

    vec3 P = ro + rdb * t;

    // Map hit point on plane to UV in [0,1] (no clamp → use mask)
    vec3 L = P - C;
    float x = dot(L, r);
    float y = dot(L, u);
    vec2 uv = vec2(x / (2.0 * W) + 0.5, 0.5 - y / (2.0 * H));

    bool inside = (uv.x >= 0.0 && uv.x <= 1.0 && uv.y >= 0.0 && uv.y <= 1.0);
    vec3 img = inside ? texture(iChannel0, uv).rgb : vec3(0.0);

    vec3 col = captured ? vec3(0.0) : img;
    col = min(col + overlay, 1.0);

    // Tonemap + gamma
    col = col / (1.0 + col);
    col = pow(col, vec3(1.0 / 2.2));
    return col;
  }
  */
/*
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;
  fragColor = vec4(render(p, fragCoord), 1.0);
}
*/

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 p = (fragCoord / iResolution.xy) * 2.0 - 1.0;
  p.x *= iResolution.x / iResolution.y;
  vec4 texel = texture(iChannel0, fragCoord / iResolution.xy);
  fragColor = vec4(render(p, fragCoord), 1.0);
  /// fragColor = texel;
}
