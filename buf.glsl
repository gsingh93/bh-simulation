#define vscode 1

#if vscode == 0

uniform sampler2D iChannel0;
uniform int iFrame;
uniform vec4 iMouse;
uniform vec3 iResolution;

#else

#iChannel0 "self"

#endif

const float PI = 3.14159265359;

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // Previous state stored in center texel:
  //   R = yaw
  //   G = lastDx (horizontal drag reference)
  //   B = pitch
  //   A = lastDy (vertical drag reference)
  vec4 prev = (iFrame == 0) ? vec4(0.0, 0.0, 0.0,
                                   0.0) // yaw=0, pitch=0, lastDx=0, lastDy=0
                            : texture(iChannel0, vec2(0.5, 0.5));

  float yaw = prev.r;
  float lastDx = prev.g;
  float pitch = prev.b;
  float lastDy = prev.a;

  // iMouse.zw = click position; z >= 0 means mouse is down
  bool down = (iMouse.z >= 0.0);

  float dx = 0.0;
  float dy = 0.0;

  if (down) {
    // Drag offsets from press point, normalized to [-1,1]
    dx = (iMouse.x - iMouse.z) / iResolution.x;
    dy = (iMouse.y - iMouse.w) / iResolution.y;

    // Incremental changes since last frame
    float deltaDx = dx - lastDx;
    float deltaDy = dy - lastDy;

    // Horizontal drag → yaw (full-width drag → 2π rotation)
    yaw += deltaDx * 2.0 * PI;

    // Vertical drag → pitch
    // Map full-height drag to a usable pitch range, then clamp
    float pitchScale = -PI; // adjust if you want more/less sensitivity
    pitch += deltaDy * pitchScale;
    pitch = clamp(pitch, -0.9, 0.9); // limit pitch to avoid flipping over

    // Remember current drag offsets for next frame
    lastDx = dx;
    lastDy = dy;
  } else {
    // Not dragging: keep yaw/pitch, reset drag references
    lastDx = 0.0;
    lastDy = 0.0;
  }

  // Store back: yaw, lastDx, pitch, lastDy
  fragColor = vec4(yaw, lastDx, pitch, lastDy);
}
