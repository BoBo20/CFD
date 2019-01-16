#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define L 1.0
#define M0 30
#define N0 30
#define K0 30
#define M (M0+2)
#define N (N0+2)
#define K (K0+2)
#define M1 (M+1)
#define N1 (N+1)
#define K1 (K+1)
#define Qx 3
#define Qy 3
#define Qz 3
#define imin 1
#define imax (N-2)
#define jmin 1
#define jmax (M-2)
#define kmin 1
#define kmax (K-2)
#define Re 400.0
#define CFL 0.5
#define sgn(x) ((x)==0.0?0.0:((x)<0.0?-1.0:1.0))
#define fmin(x,y) ((x)<(y)?(x):(y))
#define u_wall 0.1

void gks_ini(double****** f,double*** ux, double ***uy, double ***uz, double ***rho);

void InterpX(int k, int j,double****** f_plus, double**** xf_face);
void InterpY(int k, int i,double****** f_plus, double**** yf_face);
void InterpZ(int j, int i,double****** f_plus, double**** zf_face);
void boundary(double****** f_plus);
void Evol(double****** f_plus, double****** f, double*** ux, double*** uy, double*** uz, double*** rho,
          double**** xf_face, double**** yf_face, double**** zf_face);
double averageRHO(double*** rho);
double feq(int kx, int ky, int kz, double Ux, double Uy, double Uz, double RHO);

void output_csv(double*** rho,double*** ux,double*** uy,double*** uz,
                double* xc,double* yc,double* zc);
double dx,dy,dz,dt,RT,w,wb,tau,nu;
//double f[M][N][K][Qx][Qy][Qz]; // f~ at the cell center
//double f_plus[M][N][K][Qx][Qy][Qz]; // f_bar^+ at the cell center
//double rho[M][N][K],ux[M][N][K],uy[M][N][K],uz[M][N][K]; //cell center
//double xc[N],yc[M],zc[K];
//double xf_face[N1][Qx][Qy][Qz], yf_face[M1][Qx][Qy][Qz], zf_face[K1][Qx][Qy][Qz];      // cell interface
//double minmod(double a, double b, double c);
double ex[Qx]={0., 1., -1.};
double ey[Qy]={0., 1., -1.};
double ez[Qz]={0., 1., -1.};

double tpx[Qx]={2.0/3.0, 1.0/6.0, 1.0/6.0};
double tpy[Qy]={2.0/3.0, 1.0/6.0, 1.0/6.0};
double tpz[Qz]={2.0/3.0, 1.0/6.0, 1.0/6.0};

int re[Qx]={0,2,1};





int main()
{
  int     m;
  int     readdata;
  int     mmax;

  double****** f;
  double****** f_plus;
  double*** rho;
  double*** ux;
  double*** uy;
  double*** uz;
  double**** xf_face;
  double**** yf_face;
  double**** zf_face;
  double* xc;
  double* yc;
  double* zc;
{
  //assign memory f[K][M][N][Qx][Qy][Qz]
  f = (double******)malloc(sizeof(double*****)*K);
  for(int k=0;k<K;k++)
  {
    f[k] = (double*****)malloc(sizeof(double****)*M);
    for(int m=0;m<M;m++)
    {
      f[k][m] = (double****)malloc(sizeof(double***)*N);
      for(int n=0;n<N;n++)
      {
        f[k][m][n] = (double***)malloc(sizeof(double**)*Qx);
        for(int qx=0;qx<Qx;qx++)
        {
          f[k][m][n][qx] = (double**)malloc(sizeof(double*)*Qy);
          for(int qy=0;qy<Qy;qy++)
          {
            f[k][m][n][qx][qy] = (double*)malloc(sizeof(double)*Qz);
          }
        }
      }
    }
  }

  //assign memory f_plus[K][M][N][Qx][Qy][Qz]
  f_plus = (double******)malloc(sizeof(double*****)*K);
  for(int k=0;k<K;k++)
  {
    f_plus[k] = (double*****)malloc(sizeof(double****)*M);
    for(int m=0;m<M;m++)
    {
      f_plus[k][m] = (double****)malloc(sizeof(double***)*N);
      for(int n=0;n<N;n++)
      {
        f_plus[k][m][n] = (double***)malloc(sizeof(double**)*Qx);
        for(int qx=0;qx<Qx;qx++)
        {
          f_plus[k][m][n][qx] = (double**)malloc(sizeof(double*)*Qy);
          for(int qy=0;qy<Qy;qy++)
          {
            f_plus[k][m][n][qx][qy] = (double*)malloc(sizeof(double)*Qz);
          }
        }
      }
    }
  }

  //assign memory rho[K][M][N]
  rho = (double***)malloc(sizeof(double**)*K);
  for(int k=0;k<K;k++)
  {
    rho[k] = (double**)malloc(sizeof(double*)*M);
    for(int m=0;m<M;m++)
    {
      rho[k][m] = (double*)malloc(sizeof(double)*N);
    }
  }
  //assign memory ux[K][M][N]
  ux = (double***)malloc(sizeof(double**)*K);
  for(int k=0;k<K;k++)
  {
    ux[k] = (double**)malloc(sizeof(double*)*M);
    for(int m=0;m<M;m++)
    {
      ux[k][m] = (double*)malloc(sizeof(double)*N);
    }
  }
  //assign memory uy[K][M][N]
  uy = (double***)malloc(sizeof(double**)*K);
  for(int k=0;k<K;k++)
  {
    uy[k] = (double**)malloc(sizeof(double*)*M);
    for(int m=0;m<M;m++)
    {
      uy[k][m] = (double*)malloc(sizeof(double)*N);
    }
  }
  //assign memory uz[K][M][N]
  uz = (double***)malloc(sizeof(double**)*K);
  for(int k=0;k<K;k++)
  {
    uz[k] = (double**)malloc(sizeof(double*)*M);
    for(int m=0;m<M;m++)
    {
      uz[k][m] = (double*)malloc(sizeof(double)*N);
    }
  }

  //assign memory xf_face
  xf_face = (double****)malloc(sizeof(double***)*N);
  for(int n=0;n<N;n++)
  {
    xf_face[n] = (double***)malloc(sizeof(double**)*Qx);
    for(int qx=0;qx<Qx;qx++)
    {
      xf_face[n][qx] = (double**)malloc(sizeof(double*)*Qy);
      for(int qy=0;qy<Qy;qy++)
      {
        xf_face[n][qx][qy] = (double*)malloc(sizeof(double)*Qz);
      }
    }
  }
  //assign memory yf_face
  yf_face = (double****)malloc(sizeof(double***)*N);
  for(int n=0;n<N;n++)
  {
    yf_face[n] = (double***)malloc(sizeof(double**)*Qx);
    for(int qx=0;qx<Qx;qx++)
    {
      yf_face[n][qx] = (double**)malloc(sizeof(double*)*Qy);
      for(int qy=0;qy<Qy;qy++)
      {
        yf_face[n][qx][qy] = (double*)malloc(sizeof(double)*Qz);
      }
    }
  }
  //assign memory zf_face
  zf_face = (double****)malloc(sizeof(double***)*N);
  for(int n=0;n<N;n++)
  {
    zf_face[n] = (double***)malloc(sizeof(double**)*Qx);
    for(int qx=0;qx<Qx;qx++)
    {
      zf_face[n][qx] = (double**)malloc(sizeof(double*)*Qy);
      for(int qy=0;qy<Qy;qy++)
      {
        zf_face[n][qx][qy] = (double*)malloc(sizeof(double)*Qz);
      }
    }
  }

  //assign memory xc[M],yc[N],zc[K]
  xc = (double*)malloc(sizeof(double)*M);
  yc = (double*)malloc(sizeof(double)*M);
  zc = (double*)malloc(sizeof(double)*M);

}
  double  err;
  double  u_old;
  double  t;

  clock_t time_start, time_end;

  dx = L / (imax - imin + 1);
  dy = L / (jmax - jmin + 1);
  dz = L / (kmax - kmin + 1);
  RT = 1.0 / 3;
  nu = L * u_wall / Re;

  for(m = 0; m < Qx; m++)
  {
    ex[m] = sqrt(3*RT) * ex[m];
    ey[m] = sqrt(3*RT) * ey[m];
    ez[m] = sqrt(3*RT) * ez[m];
  }

  for(m = imin - 1; m <= imax + 1; m++)
  { 
    xc[m] = 0.5 * dx + (m-1) * dx;
    yc[m] = 0.5 * dy + (m-1) * dy;
    zc[m] = 0.5 * dz + (m-1) * dz;
  }


  dt = CFL * dx / sqrt(6*RT);
  tau = nu / RT;
  w = 1.5 * dt / (2*tau + dt);
  wb = dt / (2*tau + 0.5*dt); //for interface

  printf("w=%e \n",w);
  gks_ini(f,ux, uy, uz, rho);
  u_old = ux[M/2][N/2][K/2];
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
    Evol(f_plus,f,ux, uy, uz, rho, xf_face, yf_face,zf_face);
    if(m%1000 == 0)
    {
      err = fabs(ux[K/2][M/2][N/2] - u_old) / u_wall;
      u_old = ux[M/2][N/2][K/2];
      time_end = clock();
      printf("m=%d  err:%e  u:%e  average_rho:%f duration:%.1fs\n", 
              m,
              err,
              u_old,
              averageRHO(rho),
              (double)(time_end-time_start)/CLOCKS_PER_SEC);
    }  
  }
  printf("Continue? (yes=1 no=0)\n");
  scanf("%d",&readdata);
  if(readdata)
  { 
    goto AA;
  }

  output_csv(rho,ux, uy, uz, xc,yc,zc);


  return 0;
}

void gks_ini(double****** f,double*** ux, double ***uy, double ***uz, double ***rho)
{
  int i, j, k, kx, ky, kz;
  for(k = kmin; k <= kmax; k++)
  for(j = jmin; j <= jmax; j++)
  for(i = imin; i <= imax; i++)
  {
      ux[k][j][i] = 0.0; 
      uy[k][j][i] = 0.0;
      uz[k][j][i] = 0.0;
      rho[k][j][i] = 1.0;
      for(kx = 0; kx < Qx; kx++) 
      for(ky = 0; ky < Qy; ky++) 
      for(kz = 0; kz < Qz; kz++)
         f[k][j][i][kx][ky][kz]=feq(kx,ky,kz,ux[k][j][i],uy[k][j][i],uz[k][j][i],rho[k][j][i]);
  }
}

double feq(int kx, int ky, int kz, double Ux, double Uy, double Uz, double RHO)
{
  double uvw,eu,x;

  eu = (ex[kx]*Ux + ey[ky]*Uy + ez[kz]*Uz) / RT;
  uvw = (Ux*Ux + Uy*Uy + Uz*Uz) / RT;
  x = tpx[kx]*tpy[ky]*tpy[kz]*RHO*(1.0+eu+0.5*(eu*eu-uvw));
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

void boundary(double****** f_plus)
{
  int i, j, k, kx, ky, kz;
  double rho_b;

// front & back walls
  for(k=kmin;k<=kmax;k++)
  for(i=imin;i<=imax;i++)
  for(kx=0;kx<Qx;kx++)
  for(ky=0;ky<Qy;ky++)
  for(kz=0;kz<Qz;kz++)
  {

    f_plus[k][jmin-1][i][kx][ky][kz] = 2.0*f_plus[k][jmin][i][kx][ky][kz]-f_plus[k][jmin+1][i][kx][ky][kz];

    //printf("f_plus[%d][jmax+1][%d][%d][%d][%d]=%f\n", k,i,kx,ky,kz,f_plus[k][jmax+1][i][kx][ky][kz]);

    f_plus[k][jmax+1][i][kx][ky][kz] = 2.0*f_plus[k][jmax][i][kx][ky][kz]-f_plus[k][jmax-1][i][kx][ky][kz];      
  }


  for(k=kmin;k<=kmax;k++)
  for(i=imin;i<=imax;i++) //front wall
  for(ky=0;ky<Qy;ky++) //bounce back
  {
    if(ey[ky]>0) 
      for(kx=0;kx<Qx;kx++)
      for(kz=0;kz<Qz;kz++)
      {
        f_plus[k][jmin-1][i][kx][ky][kz]=f_plus[k][jmin-1][i][re[kx]][re[ky]][re[kz]]
        +f_plus[k][jmin][i][re[kx]][re[ky]][re[kz]]-f_plus[k][jmin][i][kx][ky][kz];
      }
  }

  for(k=kmin;k<=kmax;k++)
  for(i=imin;i<=imax;i++) //back wall
  for(ky=0;ky<Qy;ky++) //bounce back
  {
    if(ey[ky]<0) 
      for(kx=0;kx<Qx;kx++)
      for(kz=0;kz<Qz;kz++)
      {
        f_plus[k][jmax+1][i][kx][ky][kz]=f_plus[k][jmax+1][i][re[kx]][re[ky]][re[kz]]
        +f_plus[k][jmax][i][re[kx]][re[ky]][re[kz]]-f_plus[k][jmax][i][kx][ky][kz];
      }
  }
  //left & right walls
    for(k=kmin;k<=kmax;k++)
    for(j=jmin;j<=jmax;j++)
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qz;kz++)
    {
      f_plus[k][j][imax+1][kx][ky][kz] = 2.0*f_plus[k][j][imax][kx][ky][kz]-f_plus[k][j][imax-1][kx][ky][kz];
      f_plus[k][j][imin-1][kx][ky][kz] = 2.0*f_plus[k][j][imin][kx][ky][kz]-f_plus[k][j][imin+1][kx][ky][kz];
    }      

    for(k=kmin;k<=kmax;k++)
    for(j=jmin;j<=jmax;j++) //left wall
    for(kx=0;kx<Qx;kx++) //bounce back
    {
      if(ex[kx]>0) 
        for(ky=0;ky<Qy;ky++)
        for(kz=0;kz<Qz;kz++)
        {
          f_plus[k][j][imin-1][kx][ky][kz]=f_plus[k][j][imin-1][re[kx]][re[ky]][re[kz]]
          +f_plus[k][j][imin][re[kx]][re[ky]][re[kz]]-f_plus[k][j][imin][kx][ky][kz];
        }
    }

    for(k=kmin;k<=kmax;k++)  
    for(j=jmin;j<=jmax;j++) //right wall
    for(kx=0;kx<Qx;kx++) //bounce back
    {
      if(ex[kx]<0) 
        for(ky=0;ky<Qy;ky++)
        for(kz=0;kz<Qz;kz++)
        {
          f_plus[k][j][imax+1][kx][ky][kz]=f_plus[k][j][imax+1][re[kx]][re[ky]][re[kz]]
          +f_plus[k][j][imax][re[kx]][re[ky]][re[kz]]-f_plus[k][j][imax][kx][ky][kz];
        }
    }

// top & bottom
  for(j=jmin-1;j<=jmax+1;j++)
  for(i=imin-1;i<=imax+1;i++) 
  for(kx=0;kx<Qx;kx++)
  for(ky=0;ky<Qy;ky++)
  for(kz=0;kz<Qz;kz++) 
  {
    f_plus[kmin-1][j][i][kx][ky][kz]= 2*f_plus[kmin][j][i][kx][ky][kz]-f_plus[kmin+1][j][i][kx][ky][kz];
    f_plus[kmax+1][j][i][kx][ky][kz]= 2*f_plus[kmax][j][i][kx][ky][kz]-f_plus[kmax-1][j][i][kx][ky][kz];
  }      

  for(j=jmin-1;j<=jmax+1;j++)
  for(i=imin-1;i<=imax+1;i++) //bottom wall
    for(kz=0;kz<Qz;kz++) //bounce back
    {
      if(ez[kz]>0) 
        for(kx=0;kx<Qx;kx++)
        for(ky=0;ky<Qy;ky++)
        {
          f_plus[kmin-1][j][i][kx][ky][kz]=f_plus[kmin-1][j][i][re[kx]][re[ky]][re[kz]]
          +f_plus[kmin][j][i][re[kx]][re[ky]][re[kz]]-f_plus[kmin][j][i][kx][ky][kz];
        }
    }

  for(j=jmin-1;j<=jmax+1;j++)
  for(i=imin-1;i<=imax+1;i++)
  {
    rho_b=0.0;
    for(kz=0;kz<Qz;kz++) //bounce back
    {
      if(ez[kz]==0) 
        for(kx=0;kx<Qx;kx++)
        for(ky=0;ky<Qy;ky++) 
          rho_b+=0.5*(f_plus[kmax+1][j][i][kx][ky][kz]+f_plus[kmax][j][i][kx][ky][kz]);
      else if(ez[kz]>0) 
        for(kx=0;kx<Qx;kx++) 
        for(ky=0;ky<Qy;ky++) 
          rho_b+=(f_plus[kmax+1][j][i][kx][ky][kz]+f_plus[kmax][j][i][kx][ky][kz]);
    }
    for(kz=0;kz<Qz;kz++)
    {
      if(ez[kz]<0) 
        for(kx=0;kx<Qx;kx++)
        for(ky=0;ky<Qy;ky++)
          f_plus[kmax+1][j][i][kx][ky][kz]=
          f_plus[kmax+1][j][i][re[kx]][re[ky]][re[kz]]
          +f_plus[kmax][j][i][re[kx]][re[ky]][re[kz]]-f_plus[kmax][j][i][kx][ky][kz]
          + 4*rho_b*tpx[kx]*tpy[ky]*tpz[kz]*(ex[kx]*u_wall)/RT;//why 4, thesis is 2
    }
  }
}

void InterpX(int k, int j, double****** f_plus, double**** xf_face)   // f at cell interface: X-direction
{
  int i, kx, ky, kz, iL,jL,jR;
  double x, y, z, fc, dfx, dfy, dfz, ux_face, uy_face, uz_face, rho_face;
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
   
   for(kx=0;kx<Qx;kx++)
   for(ky=0;ky<Qy;ky++)
   for(kz=0;kz<Qz;kz++) 
    {
      fc=0.5*(f_plus[k][j][i][kx][ky][kz]+f_plus[k][j][iL][kx][ky][kz]);

      dfx=(f_plus[k][j][i][kx][ky][kz]-f_plus[k][j][iL][kx][ky][kz])/dx;

      dfy=0.5*AL*((f_plus[k][jR][iL][kx][ky][kz]+f_plus[k][jR][i][kx][ky][kz])
                - (f_plus[k][jL][iL][kx][ky][kz]+f_plus[k][jL][i][kx][ky][kz]));

      dfz=0.5*AL*((f_plus[k+1][j][iL][kx][ky][kz]+f_plus[k+1][j][i][kx][ky][kz])
                - (f_plus[k-1][j][iL][kx][ky][kz]+f_plus[k-1][j][i][kx][ky][kz]));

      x=0.5*ex[kx]*dt; 
      y=0.5*ey[ky]*dt;//half time step
      z=0.5*ez[kz]*dt; 

      xf_face[i][kx][ky][kz] = fc-x*dfx-y*dfy-z*dfz;
    }
  }

//the original f at interface
  for(i=imin;i<=imax+1;i++)
  {
    ux_face = 0.0;
    uy_face = 0.0;
    uz_face = 0.0;
    rho_face= 0.0;
    for(kx=0;kx<Qx;kx++) 
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qz;kz++)
    {
      rho_face += xf_face[i][kx][ky][kz];
      ux_face += ex[kx]*xf_face[i][kx][ky][kz];
      uy_face += ey[ky]*xf_face[i][kx][ky][kz];
      uz_face += ez[kz]*xf_face[i][kx][ky][kz];
    }   
    ux_face /= rho_face; 
    uy_face /= rho_face;
    uz_face /= rho_face;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qy;kz++)
        xf_face[i][kx][ky][kz]=(1.0-0.5*wb)*xf_face[i][kx][ky][kz]+0.5*wb*feq(kx,ky,kz, ux_face, uy_face,uz_face, rho_face);
  }
}

void InterpY(int k, int i, double****** f_plus, double**** yf_face)   // f at cell interface
{
  int j, kx, ky, kz, iL, iR, jL;
  double fc, x, y,z,  ux_face, uy_face, uz_face, rho_face, dfx, dfy, dfz;
  double hR,hL,AL,AR;

  iL = i-1; 
  iR = i+1;

  hL = dx;
  hR = dx;

  AL = 1.0/(2.0*dx);
  AR = 1.0/(2.0*dx);

// y-direction: no-slip
  for(j=jmin;j<=jmax+1;j++)
  {
    jL=j-1;
    
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++) 
    for(kz=0;kz<Qz;kz++)
    {
      fc = 0.5*(f_plus[k][j][i][kx][ky][kz]+f_plus[k][jL][i][kx][ky][kz]);

      dfy = (f_plus[k][j][i][kx][ky][kz]-f_plus[k][jL][i][kx][ky][kz])/dy;

      dfx = 0.5*AL*((f_plus[k][jL][iR][kx][ky][kz]+f_plus[k][j][iR][kx][ky][kz])
                  - (f_plus[k][jL][iL][kx][ky][kz]+f_plus[k][j][iL][kx][ky][kz]));

      dfz = 0.5*AL*((f_plus[k+1][jL][i][kx][ky][kz]+f_plus[k+1][j][i][kx][ky][kz])
                  - (f_plus[k-1][jL][i][kx][ky][kz]+f_plus[k-1][j][i][kx][ky][kz]));

      x=0.5*ex[kx]*dt;
      y=0.5*ey[ky]*dt; 
      z=0.5*ez[kz]*dt;//half time step

      yf_face[j][kx][ky][kz]=fc-x*dfx-y*dfy-z*dfz;
    }
  }

//origional DFs
  for(j=jmin;j<=jmax+1;j++)
  {
    ux_face = 0.0;
    uy_face = 0.0;
    uz_face = 0.0;
    rho_face= 0.0;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qz;kz++)
    {
      ux_face += ex[kx] * yf_face[j][kx][ky][kz];
      uy_face += ey[ky] * yf_face[j][kx][ky][kz];
      uz_face += ez[kz] * yf_face[j][kx][ky][kz];

      rho_face += yf_face[j][kx][ky][kz];
    }
    ux_face/=rho_face;
    uy_face/=rho_face;
    uz_face/=rho_face;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qz;kz++)
      yf_face[j][kx][ky][kz] = (1.0-0.5*wb) * yf_face[j][kx][ky][kz] 
                         + 0.5 * wb * feq(kx,ky,kz, ux_face, uy_face,uz_face, rho_face);
  }
}

void InterpZ(int j, int i, double****** f_plus, double**** zf_face)   // f at cell interface
{
  int k, kx, ky, kz, iR, iL, jR, jL;
  double fc, x, y,z,  ux_face, uy_face, uz_face, rho_face, dfx, dfy, dfz;
  double hR,hL,AL,AR;

  iL = i-1; 
  iR = i+1;
  jL = j-1; 
  jR = j+1; 

  AL = 1.0/(2.0*dx);
  AR = 1.0/(2.0*dx);

// y-direction: no-slip
  for(k=kmin;k<=kmax+1;k++)
  {
    
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++) 
    for(kz=0;kz<Qz;kz++)
    {
      fc = 0.5*(f_plus[k-1][j][i][kx][ky][kz]+f_plus[k][j][i][kx][ky][kz]);

      dfz = (f_plus[k][j][i][kx][ky][kz]-f_plus[k-1][j][i][kx][ky][kz])/dz;

      dfx = 0.5*AL*((f_plus[k-1][j][iR][kx][ky][kz]+f_plus[k][j][iR][kx][ky][kz])
                  - (f_plus[k-1][j][iL][kx][ky][kz]+f_plus[k][j][iL][kx][ky][kz]));

      dfy = 0.5*AL*((f_plus[k-1][jR][i][kx][ky][kz]+f_plus[k][jR][i][kx][ky][kz])
                  - (f_plus[k-1][jL][i][kx][ky][kz]+f_plus[k][jL][i][kx][ky][kz]));

      x=0.5*ex[kx]*dt;
      y=0.5*ey[ky]*dt; 
      z=0.5*ez[kz]*dt;//half time step

      zf_face[k][kx][ky][kz]=fc-x*dfx-y*dfy-z*dfz;
    }
  }

//origional DFs
  for(k=jmin;k<=jmax+1;k++)
  {
    ux_face = 0.0;
    uy_face = 0.0;
    uz_face = 0.0;
    rho_face= 0.0;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qz;kz++)
    {
      uy_face += ey[ky] * zf_face[k][kx][ky][kz];
      uz_face += ez[kz] * zf_face[k][kx][ky][kz];
      ux_face += ex[kx] * zf_face[k][kx][ky][kz];

      rho_face += zf_face[k][kx][ky][kz];
    }
    ux_face/=rho_face;
    uy_face/=rho_face;
    uz_face/=rho_face;
    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qz;kz++)
      zf_face[k][kx][ky][kz] = (1.0-0.5*wb) * zf_face[k][kx][ky][kz] 
                         + 0.5 * wb * feq(kx,ky,kz, ux_face, uy_face,uz_face, rho_face);
  }
}

void Evol(double****** f_plus, double****** f, double*** ux, double*** uy, double*** uz, double*** rho,
          double**** xf_face, double**** yf_face, double**** zf_face)
{
  int i,j,k, kx, ky,kz;
  double FM;

  // f_plus in each cell
  for(k=kmin;k<=kmax;k++)
  for(j=jmin;j<=jmax;j++)
  for(i=imin;i<=imax;i++)
  for(kx=0;kx<Qx;kx++)
  for(ky=0;ky<Qy;ky++)
  for(kz=0;kz<Qz;kz++)
  {
    FM = feq(kx, ky,kz, ux[k][j][i], uy[k][j][i],uz[k][j][i], rho[k][j][i]);
    f_plus[k][j][i][kx][ky][kz] = f[k][j][i][kx][ky][kz] - w*(f[k][j][i][kx][ky][kz] - FM);
  }

  boundary(f_plus);

  //update f: X-direction
  for(k=jmin;k<=jmax;k++)
  for(j=jmin;j<=jmax;j++)
  {
	  InterpX(k,j,f_plus,xf_face);
    for(i=imin;i<=imax;i++)
    {
      for(kx=0;kx<Qx;kx++)
      for(ky=0;ky<Qy;ky++)
      for(kz=0;kz<Qz;kz++)
      {
       f[k][j][i][kx][ky][kz] = (4.0 * f_plus[k][j][i][kx][ky][kz] - f[k][j][i][kx][ky][kz]) / 3.0
		                  + ex[kx] * dt / dx * (xf_face[i][kx][ky][kz] - xf_face[i+1][kx][ky][kz]);
      }
    }
  }

  //update f: Y-direction
  for(k=imin;k<=imax;k++)
  for(i=imin;i<=imax;i++)
  {  
	  InterpY(k,i,f_plus, yf_face);
    for(j=jmin;j<=jmax;j++)
    {
      
      for(kx=0;kx<Qx;kx++) 
      for(ky=0;ky<Qy;ky++)
      for(kz=0;kz<Qz;kz++)
      {
       f[k][j][i][kx][ky][kz] += ey[ky] * dt / dy * (yf_face[j][kx][ky][kz] - yf_face[j+1][kx][ky][kz]);
      }
    }
  }
  //update f: Z-direction
  for(j=jmin;j<=jmax;j++)
  for(i=imin;i<=imax;i++)
  {  
    InterpZ(j,i,f_plus, zf_face);
    for(k=kmin;k<=kmax;k++)
    {
      
      for(kx=0;kx<Qx;kx++) 
      for(ky=0;ky<Qy;ky++)
      for(kz=0;kz<Qz;kz++)
      {
       f[k][j][i][kx][ky][kz] += ez[kz] * dt / dz * (zf_face[k][kx][ky][kz] - zf_face[k+1][kx][ky][kz]);
      }
    }
  }
  
  //update macroscopic variables in each cell
  for(k=kmin;k<=kmax;k++) 
  for(j=jmin;j<=jmax;j++) 
  for(i=imin;i<=imax;i++)
  {
    rho[k][j][i] = 0.0;
    ux[k][j][i] = 0.0;
    uy[k][j][i] = 0.0;
    uz[k][j][i] = 0.0;

    for(kx=0;kx<Qx;kx++)
    for(ky=0;ky<Qy;ky++)
    for(kz=0;kz<Qz;kz++)
    {
      rho[k][j][i]+=f[k][j][i][kx][ky][kz];
      ux[k][j][i]+=ex[kx]*f[k][j][i][kx][ky][kz];
      uy[k][j][i]+=ey[ky]*f[k][j][i][kx][ky][kz];
      uz[k][j][i]+=ez[kz]*f[k][j][i][kx][ky][kz];
    }
    ux[k][j][i]/=rho[k][j][i]; 
    uy[k][j][i]/=rho[k][j][i];
    uz[k][j][i]/=rho[k][j][i];
  }
  
}

double averageRHO(double*** rho)
{
  double RHO_average;
  int i,j,k;
  for(k=kmin;k<=kmax;k++) 
  for(j=jmin;j<=jmax;j++) 
  for(i=imin;i<=imax;i++)
  {
    RHO_average += rho[k][j][i];
  }
  return RHO_average/(M0*N0*K0);
}



void output_csv(double*** rho,double*** ux,double*** uy,double*** uz,
                double* xc,double* yc,double* zc)
{
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>output to file\n");
  char file_name[100];
  sprintf(file_name, "./output/Re%d_N%d.csv", (int)Re, N1);
  FILE* fpt = NULL;
  fpt = fopen(file_name, "w");

  fprintf(fpt, "x coord,y coord,z coord,rho,U,V,W\n");

  for (int k = kmin; k <= kmax; k++)
  for (int j = jmin; j <= jmax; j++)
  for (int i = imin; i <= imax; i++) 
      fprintf(fpt, "%f,%f,%f,%f,%f,%f,%f\n", xc[i], yc[j],zc[k], rho[k][j][i], ux[k][j][i], uy[k][j][i],uz[k][j][i]);
  fclose(fpt);
}
