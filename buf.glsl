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
  // Read previous state (center texel)
  vec4 prev = (iFrame == 0)
                  ? vec4(0.0, 0.0, 0.0, 1.0) // yaw=0, lastDx=0 on first frame
                  : texture(iChannel0, vec2(0.5, 0.5));

  float yaw = prev.r;    // current yaw angle
  float lastDx = prev.g; // last drag offset (normalized)

  // Shadertoy: iMouse.z >= 0 means mouse is down
  bool down = (iMouse.z >= 0.0);

  float dx = 0.0;
  if (down) {
    // Drag offset from press point, normalized to [-1,1]
    dx = (iMouse.x - iMouse.z) / iResolution.x;

    // Incremental change since last frame
    float deltaDx = dx - lastDx;

    // Add only the incremental drag to yaw
    // Full-width drag → 2π rotation
    yaw += deltaDx * 2.0 * PI;

    // Remember current drag offset for next frame
    lastDx = dx;
  } else {
    // Not dragging: keep yaw, reset drag reference
    lastDx = 0.0;
  }

  // Store back: yaw in R, lastDx in G
  fragColor = vec4(yaw, lastDx, 0.0, 1.0);
}
