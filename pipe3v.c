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

vector h[];
h.t[top] = contact_angle (60. * pi / 180.);

// 四面墙壁绝对封闭边界
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
  init_grid (256); 
  
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
  
  sprintf(outdir, "MassCheck_THot%.2f_TCold%.2f_%04d%02d%02d_%02d%02d", 
          T_HOT, T_COLD,
          timeinfo->tm_year + 1900, timeinfo->tm_mon + 1,     
          timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min);        
  
  mkdir(outdir, 0777);
  run();
}

event init (t = 0) {
  scalar TL = liq->T, TG = gas->T;
  
  foreach() {
    is_solid[] = (y >= R_inner) ? 1. : 0.;
    
    if (is_solid[] > 0.) {
        f[] = 1.0; 
        double T_local = T_SAT;
        if (x < L_pipe / 6.0) { T_local = T_HOT; } 
        else if (x > 5.0 * L_pipe / 6.0) { T_local = T_COLD; }
        TL[] = T_local; TG[] = T_local; T[] = T_local;
    } else {
        if (x <= L_pipe / 2.0) {
            f[] = 1.0; 
        } else if (x > 5.0 * L_pipe / 6.0 && y > R_inner - 0.001) {
            f[] = 1.0; 
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

event solid_physics (i++) {
  scalar TL = liq->T;
  foreach() {
    if (is_solid[] > 0.) {
      u.x[] = 0.; u.y[] = 0.;
      if (x < L_pipe / 6.0) { TL[] = T_HOT; T[] = T_HOT; } 
      else if (x > 5.0 * L_pipe / 6.0) { TL[] = T_COLD; T[] = T_COLD; } 
      else { TL[] = T_SAT; T[] = T_SAT; }
    }
  }
}

// ==========================================
// 核心对账事件：每一维度的质量与流向审计
// ==========================================
event log_progress (i += 10) {
  double min_p = statsf(p).min;
  double max_p = statsf(p).max;
  
  // 初始化流体区（y < 5mm）的各项物理质量账本
  double fluid_liquid_mass = 0.; 
  double fluid_gas_mass = 0.;    
  
  // 初始化固体区（y >= 5mm）的质量账本
  double solid_water_mass = 0.;

  foreach(reduction(+:fluid_liquid_mass) reduction(+:fluid_gas_mass) reduction(+:solid_water_mass)) {
    double cell_vol = dv();
    
    if (y < R_inner) {
        // 物理步骤 1：基于 VOF 体积分数 f 计算当前网格的液相与气相质量
        fluid_liquid_mass += f[] * rho1 * cell_vol;
        fluid_gas_mass    += (1.0 - f[]) * rho2 * cell_vol;
    } else {
        // 物理步骤 2：监控上方伪固体内的质量，看其是否越界
        solid_water_mass  += f[] * rho1 * cell_vol;
    }
  }
  
  // 计算流体计算域内的总质量
  double total_fluid_mass = fluid_liquid_mass + fluid_gas_mass;

  if (t > 0) {
      fprintf(stderr, "\n--- [质量对账单 t = %.5f s] ---\n", t);
      fprintf(stderr, " 1. 固体区总质量: %.6e (理论上应恒定不变)\n", solid_water_mass);
      fprintf(stderr, " 2. 流体区液相质量: %.6e (水变蒸汽时此值应减小)\n", fluid_liquid_mass);
      fprintf(stderr, " 3. 流体区气相质量: %.6e (水变蒸汽时此值应增大)\n", fluid_gas_mass);
      fprintf(stderr, " 4. 流体区全相总质量: %.6e (两相转换的最终基石，若波动则不守恒)\n", total_fluid_mass);
      fprintf(stderr, " 5. 泊松极值压强:  [% .1f, % .1f]\n", min_p, max_p);
  }
}

event movie (t += 0.005) { 
  char name_T[256], name_water[256], name_uy[256];
  
  scalar pure_water[];
  foreach() pure_water[] = (y < R_inner) ? f[] : 0.0; 
  sprintf(name_water, "%s/f_PureWater-%06.4f.png", outdir, t); 
  output_ppm (pure_water, file = name_water, n = 1024, box = {{0,0},{L_pipe, 0.010}}, min = 0, max = 1, map = cool_warm);
              
  sprintf(name_T, "%s/T-%06.4f.png", outdir, t);
  output_ppm (T, file = name_T, n = 1024, box = {{0,0},{L_pipe, 0.010}}, min = T_COLD, max = T_HOT, map = jet); 

  // 物理步骤 3：输出垂直法向速度场，直接监控交界面上的源项物理流动方向
  scalar uy[];
  foreach() uy[] = u.y[];
  sprintf(name_uy, "%s/uy-%06.4f.png", outdir, t);
  output_ppm (uy, file = name_uy, n = 1024, box = {{0,0},{L_pipe, 0.010}}, min = -0.1, max = 0.1, map = cool_warm); 
}

event end_backup (t = 10.0) {
    fprintf(stderr, "保底结束。\n");
}
