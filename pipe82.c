#include "axi.h"
#include "ecipriano/navier-stokes/velocity-jump.h" 
#include "two-phase.h"                  
#include "tension.h"
#include "contact.h"
#include "curvature.h"
#include "ecipriano/boiling.h" 

#define L_pipe 0.005
#define R_pipe 0.005

#define T_HOT 413.15      
#define T_SAT 373.15     

vector h[];
h.t[top] = contact_angle (60. * pi / 180.);

// 设定 T 的边界条件：左侧高温，右侧维持在饱和温度
T[left] = dirichlet(T_HOT);   
T[right] = dirichlet(T_SAT); 
T[top] = neumann(0.);         

uf.n[left] = 0.; 
uf.n[right] = 0.; 
uf.n[top] = 0.;

int main() {
  size (L_pipe);
  origin (0, 0);
  init_grid (256); 
  
  rho1 = 1000.0; mu1 = 0.08;      
  rho2 = 1.2;    mu2 = 0.000018;   
  
  f.height = h;
  f.sigma = 0.072; 
  
  cp1 = 4182.0;       
  cp2 = 1005.0;       
  lambda1 = 30.0;     // 液相导热率
  lambda2 = 0.025;    
  dhev = 2.26e6; //2.26e6     
  
  double Tsat = T_SAT;  
  TL0 = Tsat;        // 液体初始处于饱和预热态
  TG0 = Tsat;        // 气体初始处于饱和预热态     
  TIntVal = Tsat;    
  
  TOLERANCE = 1e-4;
  nv = 2; 

  run();
}

event init (t = 0) {
  fraction (f, L_pipe/2.0 - x);
  
  // 获取系统底层真实的液相和气相温度指针
  scalar TL = liq->T, TG = gas->T;
  
  foreach() {
    TL[] = f[] * TL0;
    TG[] = (1. - f[]) * TG0;
    T[]  = TL[] + TG[]; 
  }
  
  // 将外面 T 的恒温边界条件，拷贝给底层真正参与计算的 TL 和 TG
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
  fprintf(stderr, ">>> 追踪: 中心X=%.5f, 管壁X=%.5f\n", x_center, x_wall);
}

event movie (t += 0.001) {
  char name_f[80], name_T[80], name_u[80];
  
  sprintf(name_f, "f-%06.4f.png", t); 
  output_ppm (f, file = name_f, n = 512, 
              box = {{0,0},{L_pipe, R_pipe}}, 
              min = 0, max = 1, map = cool_warm);
              
  sprintf(name_T, "T-%06.4f.png", t);
  output_ppm (T, file = name_T, n = 512, 
              box = {{0,0},{L_pipe, R_pipe}}, 
              min = T_SAT, max = T_HOT, map = jet); 
              
  scalar un[];
  foreach() un[] = norm(u);
  sprintf(name_u, "u-%06.4f.png", t);
  output_ppm (un, file = name_u, n = 512, 
              box = {{0,0},{L_pipe, R_pipe}}, 
              min = 0, max = 0.5, map = jet); 
}

// 定义全局变量记录上一次的液面中心位置
double last_x_center = 0.0;

// 独立的收敛检查模块 (每 0.01 秒检查一次)
event check_steady (t += 0.01) {
  double current_x_center = 0.0;
  double n_count = 0.0;
  
  // 寻找中心轴附近的液面X坐标
  foreach() {
    if (f[] > 0.01 && f[] < 0.99 && y < R_pipe * 0.1) {
       current_x_center += x;
       n_count += 1.0;
    }
  }
  if (n_count > 0.0) current_x_center /= n_count;
  
  // 计算0.01秒内的位移差
  double dx = fabs(current_x_center - last_x_center);
  double current_max_u = statsf(u.x).max;
  
  // 收敛判定条件：物理时间超过0.1秒，且位移极小，且最大流速极低
  if (t > 0.1 && dx < 1e-6 && current_max_u < 1e-3) {
     fprintf(stderr, "\n======================================================\n");
     fprintf(stderr, ">>> 【收敛报告】系统已达到热力学稳态！\n");
     fprintf(stderr, ">>> 最终耗时: %.4f 秒\n", t);
     fprintf(stderr, ">>> 最终液面中心位置: %.6f m\n", current_x_center);
     fprintf(stderr, ">>> 内部残余最大流速: %.6f m/s\n", current_max_u);
     fprintf(stderr, "======================================================\n\n");
     return 1; // 触发此条件时，提前终止整个 run()，完成计算
  }
  
  last_x_center = current_x_center;
}

event end_backup (t = 2.0) {
    fprintf(stderr, "时间达到保底限制，强制结束。\n");
}
