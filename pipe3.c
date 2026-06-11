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

#define L_pipe 0.04     // 长度 40mm
#define R_inner 0.005   // 内部流体管径 5mm

#define T_HOT 403.15      
#define T_SAT 373.15
#define T_COLD 373.15

// 接触角验证宏
vector h[];
h.t[top] = contact_angle (60. * pi / 180.);

// 四面墙壁彻底焊死，完全封闭边界
T[left] = neumann(0.);   
T[right] = neumann(0.); 
uf.n[left] = 0.; 
uf.n[right] = 0.; 
uf.n[top] = 0.;

scalar is_solid[];
char outdir[128];

int main() {
  size (L_pipe); 
  origin (0, 0);
  init_grid (256); // 高精度网格，为了捕捉薄液膜
  
  rho1 = 1000.0; mu1 = 0.08;      
  rho2 = 1.2;    mu2 = 0.000018;   
  
  f.height = h;
  f.sigma = 0.072; 
  cp1 = 4182.0; cp2 = 1005.0;        
  lambda1 = 6.0; lambda2 = 0.025;   
  dhev = 2.26e6;     
  
  TL0 = T_SAT; TG0 = T_SAT; TIntVal = T_SAT;    
  TOLERANCE = 1e-4; nv = 2; 

  time_t rawtime; struct tm * timeinfo;
  time (&rawtime); timeinfo = localtime (&rawtime);
  
  // 将 T_HOT 和 T_COLD 拼接到文件夹命名中，保留到小数点后两位
  sprintf(outdir, "THot%.2f_TCold%.2f_%04d%02d%02d_%02d%02d", 
          T_HOT, T_COLD,
          timeinfo->tm_year + 1900, timeinfo->tm_mon + 1,     
          timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min);        
  
  mkdir(outdir, 0777);
  
  // 补全缺失的启动命令和闭合大括号
  run();
}

event init (t = 0) {
  scalar TL = liq->T, TG = gas->T;
  
  foreach() {
    is_solid[] = (y >= R_inner) ? 1. : 0.;
    
    if (is_solid[] > 0.) {
        f[] = 1.0; 
        // 初始时刻，直接让三块砖头满载温度
        double T_local = T_SAT;
        if (x < L_pipe / 6.0) {
            T_local = T_HOT;
        } else if (x > 5.0 * L_pipe / 6.0) {
            T_local = T_COLD;
        }
        TL[] = T_local; TG[] = T_local; T[] = T_local;
        
    } else {
        // 内部流体通道排布
        if (x <= L_pipe / 2.0) {
            f[] = 1.0; 
        } else if (x > 5.0 * L_pipe / 6.0 && y > R_inner - 0.001) {
            f[] = 1.0; // 注意：0.001 是 1mm 厚
        } else {
            f[] = 0.0;
        }
        TL[] = f[] * T_SAT;
        TG[] = (1. - f[]) * T_SAT;
        T[]  = TL[] + TG[];
    }
  }
  copy_bcs ({TL, TG}, T);
}

// ==========================================
// 核心魔法：全体积强制锁温！
// ==========================================
event solid_physics (i++) {
  scalar TL = liq->T;
  foreach() {
    if (is_solid[] > 0.) {
      u.x[] = 0.; 
      u.y[] = 0.;
      
      if (x < L_pipe / 6.0) {
          TL[] = T_HOT; 
          T[] = T_HOT;
      } 
      else if (x > 5.0 * L_pipe / 6.0) {
          TL[] = T_COLD; 
          T[] = T_COLD;
      } 
      else {
          TL[] = T_SAT; 
          T[] = T_SAT; 
      }
    }
  }
}

// ==========================================
// 综合监控探针：实时排查物理与数值异常
// ==========================================
event log_progress (i += 10) {
  double max_u = statsf(u.x).max;
  double max_p = statsf(p).max;
  double min_p = statsf(p).min;
  
  double evap_expansion = 0.; 
  double cond_shrinkage = 0.; 
  double total_liq_vol = 0.;  

  foreach(reduction(+:evap_expansion) reduction(+:cond_shrinkage) reduction(+:total_liq_vol)) {
    // 1. 监控流体区域（y < 5mm）的真实液体总体积，排查液膜是否凭空消失
    if (y < R_inner) {
        total_liq_vol += f[] * dv();
    }
    
    // 2. 监控系统体积平衡（散度）
    double div = 0.;
    foreach_dimension() div += (uf.x[1] - uf.x[])/Delta;
    
    if (div > 1e-6) {
        evap_expansion += div * dv();
    } else if (div < -1e-6) {
        cond_shrinkage += div * dv();
    }
  }
  
  // 净体积误差：正数代表整体在膨胀（憋气），负数代表整体在收缩（抽真空）
  double net_error = evap_expansion + cond_shrinkage;

  if (t > 0) {
      fprintf(stderr, "步数:%5d | 时间:%.5fs | 液体总体积:%.4e | 净体积误差:% .4e | 压强:[% .1f, % .1f] | 最大流速:%.2f\n", 
              i, t, total_liq_vol, net_error, min_p, max_p, max_u);
  }
}

event movie (t += 0.005) { 
  char name_T[256], name_water[256], name_u[256];
  
  // 提取纯水图像，剥离固体红色伪装
  scalar pure_water[];
  foreach() pure_water[] = (y < R_inner) ? f[] : 0.0; 
  sprintf(name_water, "%s/f_PureWater-%06.4f.png", outdir, t); 
  output_ppm (pure_water, file = name_water, n = 1024, 
              box = {{0,0},{L_pipe, 0.010}}, 
              min = 0, max = 1, map = cool_warm);
              
  sprintf(name_T, "%s/T-%06.4f.png", outdir, t);
  output_ppm (T, file = name_T, n = 1024, 
              box = {{0,0},{L_pipe, 0.010}}, 
              min = T_COLD, max = T_HOT, map = jet); 

  scalar un[];
  foreach() un[] = norm(u);
  sprintf(name_u, "%s/u-%06.4f.png", outdir, t);
  output_ppm (un, file = name_u, n = 1024, 
              box = {{0,0},{L_pipe, 0.010}}, 
              min = 0, max = 0.5, map = jet); 
}

event end_backup (t = 10.0) {
    fprintf(stderr, "保底结束。\n");
}
