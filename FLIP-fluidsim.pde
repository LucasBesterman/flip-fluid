final int SIZE = 80;
final int NUM_PARTICLES = 15_000;
float dx, dy;

final float RADIUS = 0.3;
float radiusSq;

final float DT = 1;
final int SUBSTEPS = 2;
final float REPEL_STRENGTH = 0;
final boolean SEP_PARTICLES = true;

final int PROJECTION_STEPS = 30;
final float DRIFT_CORRECTION_COEFF = 0.05;
final float TARGET_DENSITY = 5;

final float FLIP_AMOUNT = 0.9;

final PVector GRAVITY = new PVector(0, 0.01);

final float OBSTACLE_RADIUS = 2;

final float FOAM_THRESHOLD = 0.7;

final boolean DISPLAY_GRID = false;
final boolean DISPLAY_PARTICLES = true;

float[][] density;
float[][] vx;
float[][] vy;
boolean[][] s;

PVector[] pos;
PVector[] vel;
float[] foam;

int mapSize;
int[][][] map;

PVector obstacle;
int lastMoveFrame;

void setup() {
  size(800, 800);
  noStroke();

  dx = width / float(SIZE);
  dy = height / float(SIZE);
  radiusSq = sq(RADIUS);

  density = new float[SIZE][SIZE];
  vx = new float[SIZE][SIZE];
  vy = new float[SIZE][SIZE];
  s = new boolean[SIZE][SIZE];

  for (int i = 0; i < SIZE; i++) {
    s[0][i] = true;
    s[SIZE-1][i] = true;
    s[i][0] = true;
    s[i][SIZE-1] = true;
  }

  pos = new PVector[NUM_PARTICLES];
  vel = new PVector[NUM_PARTICLES];
  foam = new float[NUM_PARTICLES];

  for (int i = 0; i < NUM_PARTICLES; i++) {
    pos[i] = new PVector(random(SIZE/2, SIZE-RADIUS-1), random(RADIUS+1, SIZE-RADIUS-1));
    vel[i] = new PVector();
  }

  mapSize = floor(SIZE / (RADIUS * 2));

  obstacle = new PVector(SIZE/2, SIZE/2);
  lastMoveFrame = 0;
}

void draw() {
  // Construct hashmap
  map = new int[mapSize][mapSize][];

  for (int x = 0; x < mapSize; x++) for (int y = 0; y < mapSize; y++) map[x][y] = new int[] {};

  for (int i = 0; i < NUM_PARTICLES; i++) {
    int mx = constrain(floor(pos[i].x / (RADIUS * 2)), 0, mapSize-1);
    int my = constrain(floor(pos[i].y / (RADIUS * 2)), 0, mapSize-1);
    map[mx][my] = append(map[mx][my], i);
  }

  // Handle particle collisions
  if (SEP_PARTICLES) {
    for (int n = 0; n < SUBSTEPS; n++) {
      for (int x = 0; x < mapSize; x++) {
        for (int y = 0; y < mapSize; y++) {
          for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
              int ax = x+dx, ay = y+dy;
              if (ax < 0 || ax > mapSize-1 || ay < 0 || ay > mapSize-1) continue;
              handle_collisions(x, y, x+dx, y+dy);
            }
          }
        }
      }
    }
  }
  
  // Update partices
  boolean[][] water = new boolean[SIZE][SIZE];

  PVector obstacleVel = new PVector();

  if (mousePressed) {
    PVector prev = obstacle.copy();
    obstacle.x = mouseX / dx;
    obstacle.y = mouseY / dy;
    obstacleVel = PVector.sub(obstacle, prev).mult(DT / (frameCount - lastMoveFrame));
    lastMoveFrame = frameCount;
  }

  for (int i = 0; i < NUM_PARTICLES; i++) {
    vel[i].add(PVector.mult(GRAVITY, DT));
    pos[i].add(PVector.mult(vel[i], DT));

    PVector sep = PVector.sub(pos[i], obstacle);
    float dist = sep.mag();

    float sumRadius = OBSTACLE_RADIUS + RADIUS;

    if (dist < sumRadius) {
      PVector normal = PVector.div(sep, dist);
      float overlap = sumRadius - dist;
      PVector correction = PVector.mult(normal, overlap);
      pos[i].add(correction);

      PVector relativeVel = PVector.sub(vel[i], obstacleVel);
      float dot = PVector.dot(relativeVel, normal);

      PVector normalVel = PVector.mult(normal, dot);
      vel[i].sub(normalVel);
    }

    if (pos[i].x < RADIUS+1 || pos[i].x > SIZE-RADIUS-1) {
      vel[i].x *= -0.25;
      pos[i].x = constrain(pos[i].x, RADIUS+1, SIZE-RADIUS-1);
    }

    if (pos[i].y < RADIUS+1 || pos[i].y > SIZE-RADIUS-1) {
      vel[i].y *= -0.25;
      pos[i].y = constrain(pos[i].y, RADIUS+1, SIZE-RADIUS-1);
    }

    int cx = constrain(floor(pos[i].x), 1, SIZE-2);
    int cy = constrain(floor(pos[i].y), 1, SIZE-2);
    water[cx][cy] = true;
  }

  for (int x = 1; x < SIZE-1; x++) for (int y = 1; y < SIZE-1; y++) s[x][y] = !water[x][y];

  // Transfer particle velocities to grid
  density = new float[SIZE][SIZE];
  vx = new float[SIZE][SIZE];
  vy = new float[SIZE][SIZE];

  for (int i = 0; i < NUM_PARTICLES; i++) {
    int cx = floor(pos[i].x);
    int cy = floor(pos[i].y);

    float u = pos[i].x - cx;
    float v = pos[i].y - cy;

    int i0 = constrain(cx, 1, SIZE-2);
    int i1 = constrain(cx+1, 2, SIZE-1);
    int j0 = constrain(cy, 1, SIZE-2);
    int j1 = constrain(cy+1, 2, SIZE-1);

    float w00 = (s[i1][j0]? 1 : (1-u)) * (s[i0][j1]? 1 : (1-v));
    float w01 = (s[i1][j1]? 1 : (1-u)) * (s[i0][j0]? 1 : v);
    float w10 = (s[i0][j0]? 1 : u)     * (s[i1][j1]? 1 : (1-v));
    float w11 = (s[i0][j1]? 1 : u)     * (s[i1][j0]? 1 : v);

    if (!s[i0][j0]) transfer_velocity(i0, j0, vel[i], w00);
    if (!s[i0][j1]) transfer_velocity(i0, j1, vel[i], w01);
    if (!s[i1][j0]) transfer_velocity(i1, j0, vel[i], w10);
    if (!s[i1][j1]) transfer_velocity(i1, j1, vel[i], w11);
  }

  for (int x = 0; x < SIZE; x++) {
    for (int y = 0; y < SIZE; y++) {
      if (s[x][y]) {
        vx[x][y] = 0;
        vy[x][y] = 0;
      }
      else if (density[x][y] != 0) {
        vx[x][y] /= density[x][y];
        vy[x][y] /= density[x][y];
      }
    }
  }

  // Make grid velocities incompressible
  float[][] vxn = vx;
  float[][] vyn = vy;
  project(vxn, vyn, PROJECTION_STEPS);

  // Transfer grid velocities to particles
  for (int i = 0; i < NUM_PARTICLES; i++) {
    int cx = floor(pos[i].x);
    int cy = floor(pos[i].y);

    float u = pos[i].x - cx;
    float v = pos[i].y - cy;

    int i0 = constrain(cx, 1, SIZE-2);
    int i1 = constrain(cx+1, 2, SIZE-1);
    int j0 = constrain(cy, 1, SIZE-2);
    int j1 = constrain(cy+1, 2, SIZE-1);

    float w00 = s[i0][j0]? 0 : (1-u) * (1-v);
    float w01 = s[i0][j1]? 0 : (1-u) * v;
    float w10 = s[i1][j0]? 0 : u * (1-v);
    float w11 = s[i1][j1]? 0 : u * v;
    
    float sum = w00 + w01 + w10 + w11;
    if (sum == 0) continue;
    float invSum = 1 / sum;

    float velX = (vx[i0][j0] * w00 + vx[i0][j1] * w01 + vx[i1][j0] * w10 + vx[i1][j1] * w11) * invSum;
    float velY = (vy[i0][j0] * w00 + vy[i0][j1] * w01 + vy[i1][j0] * w10 + vy[i1][j1] * w11) * invSum;

    float velXn = (vxn[i0][j0] * w00 + vxn[i0][j1] * w01 + vxn[i1][j0] * w10 + vxn[i1][j1] * w11) * invSum;
    float velYn = (vyn[i0][j0] * w00 + vyn[i0][j1] * w01 + vyn[i1][j0] * w10 + vyn[i1][j1] * w11) * invSum;

    PVector velOld = new PVector(velX, velY);
    PVector velNew = new PVector(velXn, velYn);

    PVector change = PVector.sub(velNew, velOld);
    vel[i] = PVector.lerp(velNew, PVector.add(vel[i], change), FLIP_AMOUNT);

    foam[i] = max(foam[i] - 0.01, 0);
    if (density[i0][j0] < TARGET_DENSITY * FOAM_THRESHOLD) foam[i] = 1;
  }

  // Display grid and particles
  if (DISPLAY_GRID) {
    colorMode(HSB);

    loadPixels();
    for (int x = 0; x < width; x++) {
      for (int y = 0; y < height; y++) {
        int px = x * SIZE / width;
        int py = y * SIZE / height;

        float theta = atan2(vx[px][py], vy[px][py]) + PI;
        float speed = mag(vx[px][py], vy[px][py]);
        float dens = density[px][py];
        pixels[x + y*width] = color(theta * 255 / TAU, speed * 255, dens * 20);
      }
    }
    updatePixels();
  } else background(0);

  if (DISPLAY_PARTICLES) {
    colorMode(RGB);

    for (int i = 0; i < NUM_PARTICLES; i++) {
      float px = pos[i].x * dx;
      float py = pos[i].y * dy;
      float diam = RADIUS * dx * 2;

      fill(lerpColor(#0000ff, #ffffff, foam[i]));
      circle(px, py, diam);
    }
  }

  fill(#ff0000);
  circle(obstacle.x * dx, obstacle.y * dy, OBSTACLE_RADIUS * dx * 2);
}

void transfer_velocity(int x, int y, PVector v, float weight) {
  vx[x][y] += v.x * weight;
  vy[x][y] += v.y * weight;
  density[x][y] += weight;
}

// Handle collisions between two hashmap cells
void handle_collisions(int x1, int y1, int x2, int y2) {
  for (int i : map[x1][y1]) {
    for (int j : map[x2][y2]) {
      if (i >= j) continue;

      PVector sep = PVector.sub(pos[i], pos[j]);
      float distSq = sep.magSq();

      if (distSq >= radiusSq * 4 || distSq == 0) continue;

      float dist = sep.mag();

      float overlap = RADIUS * 2 - dist;
      PVector correction = PVector.mult(sep, overlap / dist / 2 / SUBSTEPS);

      pos[i].add(correction);
      pos[j].sub(correction);

      float force = overlap / RADIUS * REPEL_STRENGTH; // Restoring force
      sep.mult(force / dist);

      vel[i].add(sep);
      vel[j].sub(sep);
    }
  }
}

// Make fluid incompressible by eliminating divergences
void project(float[][] lx, float[][] ly, int iter) {
  float[][] div = new float[SIZE][SIZE];

  for (int i = 0; i < iter; i++) {
    for (int x = 1; x < SIZE-1; x++) {
      for (int y = 1; y < SIZE-1; y++) {
        if (s[x][y]) continue;

        int sum = int(s[x+1][y]) + int(s[x-1][y]) + int(s[x][y+1]) + int(s[x][y-1]);
        div[x][y] = (lx[x+1][y] - lx[x-1][y] + ly[x][y+1] - ly[x][y-1]) / max(4 - sum, 1);
        div[x][y] -= max(DRIFT_CORRECTION_COEFF * (density[x][y] - TARGET_DENSITY), 0); // Drift correction
      }
    }

    for (int x = 1; x < SIZE-1; x++) {
      for (int y = 1; y < SIZE-1; y++) {
        if (s[x][y]) continue;

        if (!s[x+1][y]) lx[x+1][y] -= div[x][y];
        if (!s[x-1][y]) lx[x-1][y] += div[x][y];
        if (!s[x][y+1]) ly[x][y+1] -= div[x][y];
        if (!s[x][y-1]) ly[x][y-1] += div[x][y];
      }
    }
  }
}
