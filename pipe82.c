#include "axi.h"
#include "ecipriano/navier-stokes/velocity-jump.h" 
#include "two-phase.h"                  
#include "tension.h"
#include "contact.h"
#include "curvature.h"
#include "ecipriano/boiling.h" 
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define L_pipe 0.005
#define R_pipe 0.005

#define T_HOT 398.15      
#define T_SAT 373.15
#define T_COLD 372.15

vector h[];
h.t[top] = contact_angle (60. * pi / 180.);

T[left] = dirichlet(T_HOT);   
T[right] = dirichlet(T_COLD); 
T[top] = neumann(0.);         

uf.n[left] = 0.; 
uf.n[right] = 0.; 
uf.n[top] = 0.;

char outdir[128];

int main() {
  size (L_pipe);
  origin (0, 0);
  init_grid (128); 
  
  rho1 = 1000.0; mu1 = 0.08;      
  rho2 = 1.2;    mu2 = 0.000018;   
  
  f.height = h;
  f.sigma = 0.072; 
  
  cp1 = 4182.0;        
  cp2 = 1005.0;        
  lambda1 = 6.0;     
  lambda2 = 0.025;   
  dhev = 2.26e6;     
  
  double Tsat = T_SAT;  
  TL0 = Tsat;        
  TG0 = Tsat;        
  TIntVal = Tsat;    
  
  TOLERANCE = 1e-4;
  nv = 2; 

  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  sprintf(outdir, "Result_Thot%.2f_Tcold%.2f_%04d%02d%02d_%02d%02d", 
          T_HOT, T_COLD, 
          timeinfo->tm_year + 1900, 
          timeinfo->tm_mon + 1,     
          timeinfo->tm_mday,        
          timeinfo->tm_hour,        
          timeinfo->tm_min);        

  mkdir(outdir, 0777); 

  run();
}

event init (t = 0) {
  //：在右侧管壁处（x > 4.8mm）人为初始化一层 0.2mm 厚的液膜，提供冷凝界面
  fraction (f, max(L_pipe/2.0 - x, x - 0.0048));
  
  scalar TL = liq->T, TG = gas->T;
  
  foreach() {
    // 将右侧冷凝液膜及其附近的气体直接初始化为冷端温度，加速冷凝启动
    if (x > 0.0045) {
        TL[] = f[] * T_COLD;
        TG[] = (1. - f[]) * T_COLD;
    } else {
        TL[] = f[] * TL0;
        TG[] = (1. - f[]) * TG0;
    }
    T[]  = TL[] + TG[]; 
  }
  
  copy_bcs ({TL, TG}, T);
}

event log_progress (i += 10) {
  double max_u = statsf(u.x).max;
  double max_T = statsf(T).max;
  fprintf(stderr, "计算步数: %d | 物理时间: %.6f s | dt: %.2e | 最大流速: %.4f m/s | 最高温度: %.2f K\n", 
          i, t, dt, max_u, max_T);
}

event track_interface (t += 0.003) {
  double x_center = 0.0, x_wall = 0.0;
  foreach() {
    if (f[] > 0.01 && f[] < 0.99) {
      if (y < R_pipe * 0.1) x_center = x; 
      if (y > R_pipe * 0.9) x_wall = x;   
    }
  }
}

event movie (t += 0.05) { 
  char name_f[256], name_T[256], name_u[256];
  
  sprintf(name_f, "%s/f-%06.4f.png", outdir, t); 
  output_ppm (f, file = name_f, n = 512, 
              box = {{0,0},{L_pipe, R_pipe}}, 
              min = 0, max = 1, map = cool_warm);
              
  sprintf(name_T, "%s/T-%06.4f.png", outdir, t);
  output_ppm (T, file = name_T, n = 512, 
              box = {{0,0},{L_pipe, R_pipe}}, 
              min = T_COLD, max = T_HOT, map = jet); 
              
  scalar un[];
  foreach() un[] = norm(u);
  sprintf(name_u, "%s/u-%06.4f.png", outdir, t);
  output_ppm (un, file = name_u, n = 512, 
              box = {{0,0},{L_pipe, R_pipe}}, 
              min = 0, max = 0.5, map = jet); 
}

double last_x_center = 0.0;

event check_steady (t += 0.05) { 
  double current_x_center = 0.0;
  double n_count = 0.0;
  
  foreach() {
    // 这里如果存在两端液面，追踪逻辑可能会混淆，暂时仅作为参考输出
    if (f[] > 0.01 && f[] < 0.99 && y < R_pipe * 0.1 && x < L_pipe/2.0 + 0.001) {
       current_x_center += x;
       n_count += 1.0;
    }
  }
  if (n_count > 0.0) current_x_center /= n_count;
  
  double dx = fabs(current_x_center - last_x_center);
  double current_max_u = statsf(u.x).max;
  
  if (t > 0.1 && dx < 1e-7 && current_max_u < 1e-4) { 
     fprintf(stderr, ">>> 【收敛报告】系统已达到稳态！\n");
     return 1; 
  }
  
  last_x_center = current_x_center;
}

event end_backup (t = 10.0) {
    fprintf(stderr, "时间达到保底限制，强制结束。\n");
}
