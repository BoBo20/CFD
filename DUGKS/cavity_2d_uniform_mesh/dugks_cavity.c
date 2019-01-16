#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define L 1.0
#define M0 30
#define N0 30
#define M (M0+2)
#define N (N0+2)
#define M1 (M+1)
#define N1 (N+1)
#define Qx 3
#define Qy 3
#define imin 1
#define imax (N-2)
#define jmin 1
#define jmax (M-2)
#define Re 1000.0
#define CFL 0.1
#define sgn(x) ((x)==0.0?0.0:((x)<0.0?-1.0:1.0))
#define fmin(x,y) ((x)<(y)?(x):(y))
#define u_wall 0.1

void gks_ini(void);
void SlopeY(int i);
void SlopeX(int j);
void InterpX(int j);
void InterpY(int i);
void boundary(void);
void Evol(void);
double averageRHO(void);
double feq(int kx, int ky, double Ux, double Uy, double RHO);
void  datadeal(void);
void output(void);
double dx,dy, dt,RT,w,wb,tau,nu;
double f[M][N][Qx][Qy]; // f~ at the cell center
double f_plus[M][N][Qx][Qy]; // f_bar^+ at the cell center
double rho[M][N],ux[M][N],uy[M][N]; //cell center
double xc[N],yc[M];
double xf_face[N1][Qx][Qy], yf_face[M1][Qx][Qy];      // cell interface
//double minmod(double a, double b, double c);
double ex[Qx]={0., 1., -1.};
double ey[Qy]={0., 1., -1.};
double tpx[Qx]={2.0/3.0, 1.0/6.0, 1.0/6.0};
double tpy[Qy]={2.0/3.0, 1.0/6.0, 1.0/6.0};
int re[Qx]={0,2,1};

int main()
{
  int     m;
  int     readdata;
  int     mmax;

  double  err;
  double  u_old;
  double  t;

  clock_t time_start, time_end;

  dx = L / (imax - imin + 1);
  dy = L / (jmax - jmin + 1);
  RT = 1.0 / 3;
  nu = L * u_wall / Re;

  for(m = 0; m < Qx; m++)
  {
    ex[m] = sqrt(3*RT) * ex[m];
    ey[m] = sqrt(3*RT) * ey[m];
  }

  for(m = imin - 1; m <= imax + 1; m++)
  { 
    xc[m] = yc[m] = 0.5 * dx + (m-1) * dx;
  }


  dt = CFL * dx / sqrt(6*RT);
  tau = nu / RT;
  w = 1.5 * dt / (2*tau + dt);
  wb = dt / (2*tau + 0.5*dt); //for interface

  printf("w=%e \n",w);
  gks_ini();
  u_old = ux[M/2][N/2];
  m = 0;
  time_start = clock();
AA:
  printf("input mmax:\n");
  scanf("%d",&mmax);
  mmax += m;
  printf("\ndt=%lf mmax=%d\n", dt, mmax);
  err = 1.0;
  printf("start iteration...\n\n");
  while((m<mmax) && err>1.0E-6)
  {
    m++;
    Evol();
    if(m%1000 == 0)
    {
      err = fabs(ux[M/2][N/2] - u_old) / u_wall;
      u_old = ux[M/2][N/2];
      time_end = clock();
      printf("m=%d  err:%e  u:%e  average_rho:%e duration:%.1fs\n", 
              m,
              err,
              u_old,
              averageRHO(),
              (double)(time_end-time_start)/CLOCKS_PER_SEC);
    }  
  }
  datadeal();
  printf("Continue? (yes=1 no=0)\n");
  scanf("%d",&readdata);
  if(readdata)
  { 
    goto AA;
  }

  output();
  return 0;
}

void gks_ini()
{
  int i,j, kx, ky;
  for(j = jmin; j <= jmax; j++)
  for(i = imin; i <= imax; i++)
  {
      ux[j][i] = 0.0; 
      uy[j][i] = 0.0;
      rho[j][i] = 1.0;
      for(kx=0;kx<Qx;kx++) 
      for(ky=0;ky<Qy;ky++) 
         f[j][i][kx][ky]=feq(kx,ky,ux[j][i],uy[j][i],rho[j][i]);
  }
}
#pragma acc routine seq
inline double feq(int kx, int ky, double Ux, double Uy, double RHO)
{
  double uv,eu,x;

  eu = (ex[kx]*Ux + ey[ky]*Uy) / RT;
  uv = (Ux*Ux + Uy*Uy) / RT;
  x = tpx[kx]*tpy[ky]*RHO*(1.0+eu+0.5*(eu*eu-uv));
  //x=tpx[kx]*tpy[ky]*RHO*(1.0+3*eu+4.5*eu*eu-1.5*uv);
  //x=x*(1.0+tau*(ex[kx]-Ux)*force/RT);
  return x;
}
/*
double minmod(double a, double b, double c)
{
  double sa, sb, sc;
  sa=sgn(a); sb=sgn(b);sc=sgn(c);
  if(sa==sb && sb==sc)
     return sa*fmin(fabs(a),fmin(fabs(b),fabs(c)));
  else return 0;
}
*/

void boundary()
{
  int i,j, kx, ky;
  double rho_b;

//left & right walls
  for(j=jmin;j<=jmax;j++)
  for(ky=0;ky<Qy;ky++)
  for(kx=0;kx<Qx;kx++)
  {
    f_plus[j][imin-1][kx][ky] = 2.0*f_plus[j][imin][kx][ky]-f_plus[j][imin+1][kx][ky];
    f_plus[j][imax+1][kx][ky] = 2.0*f_plus[j][imax][kx][ky]-f_plus[j][imax-1][kx][ky];
  }      

  for(j=jmin;j<=jmax;j++) //left wall
  for(kx=0;kx<Qx;kx++) //bounce back
  {
    if(ex[kx]>0) 
      for(ky=0;ky<Qy;ky++)
      {
        f_plus[j][imin-1][kx][ky]=f_plus[j][imin-1][re[kx]][re[ky]]+f_plus[j][imin][re[kx]][re[ky]]-f_plus[j][imin][kx][ky];
      }
  }  
  for(j=jmin;j<=jmax;j++) //right wall
  for(kx=0;kx<Qx;kx++) //bounce back
  {
    if(ex[kx]<0) 
      for(ky=0;ky<Qy;ky++)
      {
        f_plus[j][imax+1][kx][ky]=f_plus[j][imax+1][re[kx]][re[ky]]+f_plus[j][imax][re[kx]][re[ky]]-f_plus[j][imax][kx][ky];
      }
  }

// top & bottom
  for(i=imin-1;i<=imax+1;i++) 
  for(ky=0;ky<Qy;ky++) 
  for(kx=0;kx<Qx;kx++)
  {
    f_plus[jmin-1][i][kx][ky]= 2*f_plus[jmin][i][kx][ky]-f_plus[jmin+1][i][kx][ky];
    f_plus[jmax+1][i][kx][ky]= 2*f_plus[jmax][i][kx][ky]-f_plus[jmax-1][i][kx][ky];
  }      

  for(i=imin-1;i<=imax+1;i++) //bottom wall
    for(ky=0;ky<Qy;ky++) //bounce back
    {
      if(ey[ky]>0) 
        for(kx=0;kx<Qx;kx++)
        {
          f_plus[jmin-1][i][kx][ky]=f_plus[jmin-1][i][re[kx]][re[ky]]+f_plus[jmin][i][re[kx]][re[ky]]-f_plus[jmin][i][kx][ky];
        }
    }

  for(i=imin-1;i<=imax+1;i++)
  {
    rho_b=0.0;
    for(ky=0;ky<Qy;ky++) //bounce back
    {
      if(ey[ky]==0) 
        for(kx=0;kx<Qx;kx++) 
          rho_b+=0.5*(f_plus[jmax+1][i][kx][ky]+f_plus[jmax][i][kx][ky]);
      else if(ey[ky]>0) 
        for(kx=0;kx<Qx;kx++) 
          rho_b+=(f_plus[jmax+1][i][kx][ky]+f_plus[jmax][i][kx][ky]);
    }
    for(ky=0;ky<Qy;ky++)
    {
      if(ey[ky]<0) 
        for(kx=0;kx<Qx;kx++)
          f_plus[jmax+1][i][kx][ky]=f_plus[jmax+1][i][re[kx]][re[ky]]+f_plus[jmax][i][re[kx]][re[ky]]-f_plus[jmax][i][kx][ky]
                                + 4*rho_b*tpx[kx]*tpy[ky]*(ex[kx]*u_wall)/RT;//why 4, thesis is 2
    }
  }
}

void InterpX(int j)   // f at cell interface: X-direction
{
  int i, kx, ky, iL,jL,jR;
  double x, y, fc, dfx, dfy, ux_face, uy_face, rho_face;
  double hR,hL,AL,AR,AC;

  jL=j-1;  
  jR=j+1;
  hL = dy;
  hR = dy;
  AL= 1/(2*dy); 
  AR= 1/(2*dy);
  for(i=imin;i<=imax+1;i++) // inner nodes
  {
   iL=i-1;
   
	 for(ky=0;ky<Qy;ky++) 
   for(kx=0;kx<Qx;kx++)
    {
      fc=0.5*(f_plus[j][i][kx][ky]+f_plus[j][iL][kx][ky]);

      dfx=(f_plus[j][i][kx][ky]-f_plus[j][iL][kx][ky])/dx;

      dfy=0.5*(AR*(f_plus[jR][iL][kx][ky]+f_plus[jR][i][kx][ky])
             - AL*(f_plus[jL][iL][kx][ky]+f_plus[jL][i][kx][ky]));

      x=0.5*ex[kx]*dt; 
      y=0.5*ey[ky]*dt;//half time step
      xf_face[i][kx][ky] = fc-x*dfx-y*dfy;
    }
  }

//the original f at interface
  for(i=imin;i<=imax+1;i++)
  {
    ux_face = 0.0;
    uy_face = 0.0;
    rho_face= 0.0;
    for(kx=0;kx<Qx;kx++) 
    for(ky=0;ky<Qy;ky++)
    {
      rho_face += xf_face[i][kx][ky];
      ux_face += ex[kx]*xf_face[i][kx][ky];
      uy_face += ey[ky]*xf_face[i][kx][ky];
    }   
    ux_face /= rho_face; 
    uy_face /= rho_face;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
        xf_face[i][kx][ky]=(1.0-0.5*wb)*xf_face[i][kx][ky]+0.5*wb*feq(kx,ky, ux_face, uy_face,rho_face);
  }
}

void InterpY(int i)   // f at cell interface
{
  int j, kx, ky, iL, iR, jL;
  double fc, x, y, ux_face, uy_face, rho_face, dfx, dfy;
  double hR,hL,AL,AR,AC;

  iL = i-1; 
  iR = i+1;

  hL = dx;
  hR = dx;

  AC = 0; 
  AL = 1.0/(2.0*dx);
  AR = 1.0/(2.0*dx);

// y-direction: no-slip
  for(j=jmin;j<=jmax+1;j++)
  {
    jL=j-1;
    
    for(ky=0;ky<Qy;ky++) 
    for(kx=0;kx<Qx;kx++)
    {
      fc = 0.5*(f_plus[j][i][kx][ky]+f_plus[jL][i][kx][ky]);

      dfy = (f_plus[j][i][kx][ky]-f_plus[jL][i][kx][ky])/dy;
      dfx = 0.5*(AR*(f_plus[jL][iR][kx][ky]+f_plus[j][iR][kx][ky])
               - AL*(f_plus[jL][iL][kx][ky]+f_plus[j][iL][kx][ky]));

      y=0.5*ey[ky]*dt; 
      x=0.5*ex[kx]*dt;//half time step
      yf_face[j][kx][ky]=fc-x*dfx-y*dfy;
    }
  }

//origional DFs
  for(j=jmin;j<=jmax+1;j++)
  {
    ux_face = 0.0;
    uy_face = 0.0;
    rho_face= 0.0;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    {
      ux_face += ex[kx] * yf_face[j][kx][ky];
      uy_face += ey[ky] * yf_face[j][kx][ky];

      rho_face += yf_face[j][kx][ky];
    }
    ux_face/=rho_face;
    uy_face/=rho_face;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
      yf_face[j][kx][ky] = (1.0-0.5*wb) * yf_face[j][kx][ky] 
                         + 0.5 * wb * feq(kx,ky, ux_face, uy_face,rho_face);
  }
}

void Evol()
{
  int i,j, kx, ky;
  double FM;

  // f_plus in each cell
  for(j=jmin;j<=jmax;j++)
  for(i=imin;i<=imax;i++)
  for(kx=0;kx<Qx;kx++)
  for(ky=0;ky<Qy;ky++)
  {
    FM = feq(kx, ky, ux[j][i], uy[j][i], rho[j][i]);
    f_plus[j][i][kx][ky] = f[j][i][kx][ky] - w*(f[j][i][kx][ky] - FM);
  }

  boundary();

  //update f: X-direction
  for(j=jmin;j<=jmax;j++)
  {
	  InterpX(j);
    for(i=imin;i<=imax;i++)
    {
      for(kx=0;kx<Qx;kx++)
      for(ky=0;ky<Qy;ky++)
      {
       f[j][i][kx][ky] = (4.0 * f_plus[j][i][kx][ky] - f[j][i][kx][ky]) / 3.0
		                  + ex[kx] * dt / dx * (xf_face[i][kx][ky] - xf_face[i+1][kx][ky]);
      }
    }
  }

  //update f: Y-direction
  for(i=imin;i<=imax;i++)
  {  
	  InterpY(i);
    for(j=jmin;j<=jmax;j++)
    {
      
      for(kx=0;kx<Qx;kx++) 
      for(ky=0;ky<Qy;ky++)
      {
       f[j][i][kx][ky] += ey[ky] * dt / dy * (yf_face[j][kx][ky] - yf_face[j+1][kx][ky]);
      }
    }
  }
  
  //update macroscopic variables in each cell
  for(j=jmin;j<=jmax;j++) 
  for(i=imin;i<=imax;i++)
  {
    rho[j][i] = 0.0;
    ux[j][i] = 0.0;
    uy[j][i] = 0.0;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    {
      rho[j][i]+=f[j][i][kx][ky];
      ux[j][i]+=ex[kx]*f[j][i][kx][ky];
      uy[j][i]+=ey[ky]*f[j][i][kx][ky];
    }
    ux[j][i]/=rho[j][i]; 
    uy[j][i]/=rho[j][i];
  }
  
}

double averageRHO()
{
  double RHO_average;
  int i,j;
  for(j=jmin;j<=jmax;j++) 
  for(i=imin;i<=imax;i++)
  {
    RHO_average += rho[j][i];
  }
  return RHO_average/(M0*N0);
}

void  datadeal()
   {
     int i,j;
     FILE *fp;
     fp=fopen("./output/u_center_line.dat","w");
     fprintf(fp,"y coord,u\n");
     for(j=jmin;j<=jmax;j++)  {
      fprintf(fp,"%e,%e",yc[j],ux[j][N0/2]);
      fprintf(fp,"\n");
      }
     fclose(fp);

     fp=fopen("./output/v_center_line.dat","w");
     fprintf(fp,"x coord, v\n");
     for (i=imin; i<=imax; i++)
     {
      fprintf(fp,"%e,%e",xc[i],uy[M0/2][i]);
      fprintf(fp,"\n");
    }
     fclose(fp);
}


void output()
{
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>output to file\n");
  char file_name[100];
  sprintf(file_name, "./output/Re%d_N%d.csv", (int)Re, N0);
  FILE* fpt = NULL;
  fpt = fopen(file_name, "w");
  fprintf(fpt, "x coord,y coord,z coord,rho,u,v,w\n");
  for (int j = jmin; j <= jmax; j++)
    for (int i = imin; i <= imax; i++) 
      fprintf(fpt, "%f,%f,%f,%f,%f,%f,%f\n", xc[i], yc[j], 0.0, rho[j][i], ux[j][i], uy[j][i],0.0);
  fclose(fpt);
}