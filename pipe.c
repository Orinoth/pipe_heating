#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"                 
#include "navier-stokes/conserving.h"  
#include "tension.h"
#include "contact.h"
#include "curvature.h"

#define L_pipe 0.005
#define R_pipe 0.005

vector h[];
h.t[top] = contact_angle (60. * pi / 180.);

//边界条件 
uf.n[left] = 0.; 
uf.n[right] = 0.; 
uf.n[top] = 0.;

int main() {
  size (L_pipe);
  origin (0, 0);
  
  init_grid (256); 
  
  rho1 = 1000.0; mu1 = 0.08;      
  rho2 = 1.2;    mu2 = 0.000018;   
  f.sigma = 0.072;
  
  TOLERANCE = 1e-4;

  f.height = h;
  run();
}

event init (t = 0) {
  fraction (f, L_pipe/2.0 - x);
}

// 实时监测模块
event log_progress (i += 10) {
  double max_u = statsf(u.x).max;
  fprintf(stderr, "计算步数: %d | 物理时间: %.6f 秒 | dt: %.2e | 最大流速: %.4f m/s\n", 
          i, t, dt, max_u);
}

// 智能刹车
event auto_stop (t += 0.01) {
  double current_max_u = statsf(u.x).max;
  
  if (current_max_u < 0.001 && t > 0.05) {
    fprintf(stderr, "\n大哥，液面已经完全静止，平衡态已达到！耗时: %g 秒。准备下班！\n", t);
    return 1; 
  }
}

event end_backup (t = 2.0) {
    fprintf(stderr, "时间达到 2 秒保底限制，强制结束。\n");
}

// 追踪界面位置
event track_interface (t += 0.003) {
  double x_center = 0.0, x_wall = 0.0;
  foreach() {
    if (f[] > 0.01 && f[] < 0.99) {
      if (y < R_pipe * 0.1) x_center = x; 
      if (y > R_pipe * 0.9) x_wall = x;   
    }
  }
  fprintf(stderr, ">>> 追踪: 中心X=%.5f, 管壁X=%.5f, 液面弯曲差=%.5f\n", 
          x_center, x_wall, x_wall - x_center);
}

// 输出相图
event movie (t += 0.001) {
  char name[80];
  sprintf(name, "f-%06.4f.png", t); 
  output_ppm (f, file = name, n = 512, 
              box = {{0,0},{L_pipe, R_pipe}}, 
              min = 0, max = 1, map = cool_warm);
}
