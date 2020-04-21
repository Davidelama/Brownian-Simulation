// compile with: g++ Br_Sim.cpp nrutil.cpp random_mars.cpp -o soft_discs_Br -O2 -w

#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <sys/time.h>
#include "Br_Sim.h"
#include "nrutil.h"
#include "random_mars.h"


#define PI 3.1415926

//this program simulates a gas of hard spheres with 3D MD and PBC

int main() {
  int i,j, total, per, ETA, hr, min, sec;
  double Ek, U, P;
  
  //time stuff
  struct timeval tp, tpend;
  long int ms, msend;
  //close time stuff
  
  mu=.001; //mobility coefficient
  T = 1.0; //temperature in units of epsilon/kB.
  dt = 0.005; //time-step length
  rho=0.6; //rho=density
  N = 64; //216
  sigma=pow(2,1./6.); //soft diameter
  Nblock=10;
  Nlevel=5;
  Nmeas=rint(pow(Nblock,3));
  Neq=rint(pow(Nblock,Nlevel-1));
  Nsim=rint(pow(Nblock,Nlevel+1));  //total time has to be Nblock*Nblock^Nlevel
  dim=3;
  
  L=pow((double)N/rho,1/(double)dim);
  
  random_mars = new RanMars(time(NULL));
  
  x = dvector(0,N-1);
  y = dvector(0,N-1);
  z = dvector(0,N-1);
  
  vx = dvector(0,N-1);
  vy = dvector(0,N-1);
  vz = dvector(0,N-1);
  
  fx = dvector(0,N-1);
  fy = dvector(0,N-1);
  fz = dvector(0,N-1);
  
  double **vx_blk, **vy_blk, **vz_blk, **Sumx_blk, **Sumy_blk, **Sumz_blk; //contains v for particle i, position in block j, level k
  double *Zx, *Zy, *Zz, *Drx, *Dry, *Drz;  //contains particle average of Z, (Nblock-1)*Nlevel elements
  init_block(vx_blk, vy_blk, vz_blk, Zx, Zy, Zz); //we pass to the function the value x by reference &x. If x is of type type*, we pass directly x
  init_block(Sumx_blk, Sumy_blk, Sumz_blk, Drx, Dry, Drz);
  
  
  FILE *fout, *gout, *sfout, *zout, *drout;
  init_output(fout, gout, sfout, zout, drout);  //output files
  
  latt_init_pos();
  //rand_init_pos();
  
  init_vel();
  init_cell_list();
  
  print_cell(0);
  
  print_pos(0);
  print_single_g(0);
  
  //equilibration
  for (i=0; i<Neq; i++)
  {
    if (i%Nmeas==0)
    {
      rescale_vel();  //in this way after the particles go into "realistic" positions, the temperature is still the one we wanted 
      Ek=calculate_kin();
      U=calculate_pot();
      P=calculate_pressure();
      printf("Ek=%f, U=%f, E=%f, P=%f (stabilization)\n",Ek, U,Ek+U, P);
    }
    
    //md_step();
    
    bd_step();
    
  }
  
  //start of simulation

  gettimeofday(&tp, NULL);
  ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;

  for (i=0; i<Nsim; i++)
  {
    if (i%Nmeas==0) //actual measurements
    {
      Ek=calculate_kin();
      U=calculate_pot();
      P=calculate_pressure();
      print_g(gout);
      print_sf(sfout);
      fprintf(fout,"%f %f %f\n", Ek, U, P);
    }

    if (i%(Nsim/100)==0 && i!=0)  //Just because I like ETAs
    {
      per=i/(Nsim/100);
      gettimeofday(&tpend, NULL);
      msend = tpend.tv_sec * 1000 + tpend.tv_usec / 1000;
      ETA=int((double)(msend-ms)/per*(100-per)/1000);
      hr=ETA/3600;
      min=(ETA%3600)/60;
      sec=(ETA%60);
      printf("%d%%, ETA: %dh:%dm:%ds, Ek=%f, U=%f, E=%f, P=%f \n",per,hr,min,sec,Ek, U,Ek+U, P);
      
    }
    //md_step();  //normal newtonian dynamics
    bd_step();  //brownian dynamics
    block_corr(i, vx_blk, vy_blk, vz_blk, Zx, Zy, Zz, zout);
    block_intg(i, Sumx_blk, Sumy_blk, Sumz_blk, Drx, Dry, Drz, drout);
  }
  
  close_output(fout, gout, sfout, zout, drout);
  print_pos(Nsim-1);
  print_single_g(Nsim-1);
  print_block(vx_blk, vy_blk, vz_blk);
  
  return 0;
}

/*********************************************************************************************************************************/

void rand_init_pos()  //random initial position
{
  int i, j;
  
  for (i=0; i<N; i++)
  {
    x[i]=(random_mars->uniform()-0.5)*L;
    y[i]=(random_mars->uniform()-0.5)*L;
    z[i]=(random_mars->uniform()-0.5)*L;
    
    //for (j=0; j<i; j++)
    //{
    //  if (dist(x[i],y[i],x[j],y[j])<sigma)
    //  {
    //    i--;
    //    break;
    //  }
    //}
  }
}

/*********************************************************************************************************************************/

void latt_init_pos()  //packed initial potition
{
  int i, a;
  double dpart;
  
  a=ceil(pow((double)N,1/((double)3))); //this version doesn't work with hard shells
  dpart=L/((double)a);
  
  //a=floor(L/sigma);
  
  for (i=0; i<N; i++)
  {
    x[i]=(i%a)*dpart-L*0.5;
    y[i]=((i%(a*a))/a)*dpart-L*0.5;
    z[i]=(i/(a*a))*dpart-L*0.5;
    if (z[i]>=L*0.5-dpart*0.5)
    {
      z[i]-=L*rint(z[i]/L)-dpart*.5;
      y[i]+=dpart*.5;
      x[i]+=dpart*.5;
    }
  }
}

/*********************************************************************************************************************************/

double dist(double x1, double y1, double z1, double x2, double y2, double z2)  //determines overlap of two particles
{
  double dx, dy, dz;
  
  dx=x1-x2-L*rint((x1-x2)/L); //cool way to do it
  dy=y1-y2-L*rint((y1-y2)/L);
  dz=z1-z2-L*rint((z1-z2)/L);
  
  return sqrt(dx*dx+dy*dy+dz*dz);
}

/*********************************************************************************************************************************/

double calculate_pot()
{
  double dx, dy, dz, r2, r2i, r6i, frj, sigma2, U;
  int i, p, n, j, nindex;
  
  sigma2=sigma*sigma;
  U=0;
  
  for (i=0; i<N; i++)
  {
    for (n=0; n<Nneigh; n++)
    {
      nindex=neighbor[cell_list_index[i]][n];
      for (p=0; p<Nmax; p++)
      {
        j=cell_list[nindex][p];
        if (j==-1)
          break;
        if (j!=i)
        {
          dx=x[i]-x[j]-L*rint((x[i]-x[j])/L);
          dy=y[i]-y[j]-L*rint((y[i]-y[j])/L);
          dz=z[i]-z[j]-L*rint((z[i]-z[j])/L);
          
          r2=dx*dx+dy*dy+dz*dz;
          if (r2<sigma2)
          {
            r2i=1./r2;
            r6i=r2i*r2i*r2i;
            U+=4.*(r6i*r6i-r6i)+1.;
          }
        }
      }
    }      
  }
  return U/2;
}

/*********************************************************************************************************************************/

double calculate_pressure() 
{
  double dx, dy, dz, P, r2, r2i, r6i, sigma2;
  int i, p, n, j, nindex;
  P=0.0;
  
  sigma2=sigma*sigma;   //serious way to do it
  
  for (i=0; i<N; i++) 
  {
    for (n=0; n<Nneigh; n++)
    {
      nindex=neighbor[cell_list_index[i]][n];
      for (p=0; p<Nmax; p++)
      {
        j=cell_list[nindex][p];
        if (j==-1)
          break;
        if (j!=i)
        {
          dx=x[i]-x[j]-L*rint((x[i]-x[j])/L);
          dy=y[i]-y[j]-L*rint((y[i]-y[j])/L);
          dz=z[i]-z[j]-L*rint((z[i]-z[j])/L);
          
          r2=dx*dx+dy*dy+dz*dz;
          if (r2<sigma2)
          {
            r2i=1./r2;
            r6i=r2i*r2i*r2i;
            P+=12.*(2.*r6i*r6i-r6i); // divided by 2 because we take i!=j and not i<j
          }
        }
      }
    }
  }
  
  //for (i=0;i<N;i++) //this should work but it doesn't
  //{
  //  P+=x[i]*fx[i]+y[i]*fy[i];
  //}
  return rho*T+P/(dim*pow(L,dim));
}

/*********************************************************************************************************************************/

double calculate_force()
{
  double dx, dy, dz, r2, r2i, r6i, frj, sigma2, fr;
  int i, p, n, j, nindex;
  
  sigma2=sigma*sigma;
  
  for (i=0; i<N; i++)
  {
    fx[i]=0;
    fy[i]=0;
    fz[i]=0;
    for (n=0; n<Nneigh; n++)
    {
      nindex=neighbor[cell_list_index[i]][n];
      for (p=0; p<Nmax; p++)
      {
        j=cell_list[nindex][p];
        if (j==-1)
          break;
        if (j!=i)
        {
          dx=x[i]-x[j]-L*rint((x[i]-x[j])/L);
          dy=y[i]-y[j]-L*rint((y[i]-y[j])/L);
          dz=z[i]-z[j]-L*rint((z[i]-z[j])/L);
          
          r2=dx*dx+dy*dy+dz*dz;
          if (r2<sigma2)
          {
            r2i=1./r2;
            r6i=r2i*r2i*r2i;
            fr=24.*(2.*r6i*r6i-r6i)*r2i;
            fx[i]+=fr*dx;
            fy[i]+=fr*dy;
            fz[i]+=fr*dz;
          }
        }
      }
    }      
  }
}

/*********************************************************************************************************************************/

void md_step()
{
  int i;
  
  for(i=0; i<N; i++)
  {
    vx[i]+=.5*dt*fx[i];
    vy[i]+=.5*dt*fy[i];
    vz[i]+=.5*dt*fz[i];
    
    x[i]+=dt*vx[i];
    y[i]+=dt*vy[i];
    z[i]+=dt*vz[i];
    
    x[i]-=L*rint(x[i]/L);
    y[i]-=L*rint(y[i]/L);
    z[i]-=L*rint(z[i]/L);
  }
  
  reinit_cell_list();
  
  calculate_force();
  
  for(i=0; i<N; i++)
  {
    vx[i]+=.5*dt*fx[i];
    vy[i]+=.5*dt*fy[i];
    vz[i]+=.5*dt*fz[i];
  }
}

/*********************************************************************************************************************************/

void bd_step()
{
  double dx,dy,dz,B;
  int i;
  
  B=sqrt(2.0*T*mu*dt);
  
  for(i=0; i<N; i++)
  {

    dx=mu*dt*fx[i]+B*random_mars->gaussian();
    dy=mu*dt*fy[i]+B*random_mars->gaussian();
    dz=mu*dt*fz[i]+B*random_mars->gaussian();
    
    vx[i]=dx/dt; //at least these velocities kind of make sense
    vy[i]=dy/dt;
    vz[i]=dz/dt;
    
    x[i]+=dx;
    y[i]+=dy;
    z[i]+=dz;
    
    x[i]-=L*rint(x[i]/L);
    y[i]-=L*rint(y[i]/L);
    z[i]-=L*rint(z[i]/L);
  }
  
  reinit_cell_list();
  
  calculate_force();
 
}

double calculate_kin()
{
  double Vsum;
  int i;
  
  Vsum=0;
  
  for (i=0; i<N; i++)
  {
    Vsum+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];   
  }
  Vsum*=0.5;
  
  return Vsum;
}

/*********************************************************************************************************************************/

void init_vel()
{
  int i;
  double Vx, Vy, Vz, Tnow;
  
  Vx=0;
  Vy=0;
  Vz=0;
  
  for (i=0;i<N;i++)
  {
    vx[i]= random_mars->gaussian();
    vy[i]= random_mars->gaussian();
    vz[i]= random_mars->gaussian();
    
    Vx+=vx[i];
    Vy+=vy[i];
    Vz+=vz[i];
  }
  
  //remove center of mass of moentum
  for (i=0;i<N;i++)
  {
    vx[i]-=Vx/((double)N);
    vy[i]-=Vy/((double)N);
    vz[i]-=Vz/((double)N);
  }
  
  //rescale temperature
  for (i=0;i<N;i++)
  {
    Tnow+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  Tnow/=(double)dim*N;
  
  for (i=0;i<N;i++)
  {
    vx[i]*=sqrt(T/Tnow);
    vy[i]*=sqrt(T/Tnow);
    vz[i]*=sqrt(T/Tnow);
  }
  
}

/*********************************************************************************************************************************/

void rescale_vel()
{
  int i;
  double Tnow;
  
  Tnow=0;
  
  for (i=0;i<N;i++)
  {
    Tnow+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  Tnow/=(double)dim*N;
  
  for (i=0;i<N;i++)
  {
    vx[i]*=sqrt(T/Tnow);
    vy[i]*=sqrt(T/Tnow);
    vz[i]*=sqrt(T/Tnow);
  }
  
}

/*********************************************************************************************************************************/

void print_pos(int step)  //prints the positions of the particles
{
  int i;
  
  std::stringstream sstr;
  sstr << "positions_step" << step << ".dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();
  
  FILE * out = fopen(cstr,"w");
  
  fprintf(out,"%f %f %f\n",L,sigma,rc);
  
  for (i=0; i<N; i++)
  {
    fprintf(out,"%f %f %f\n",x[i],y[i], z[i]);
  }
  fclose(out);
}

/*********************************************************************************************************************************/

void init_cell_list()
{
  int i, j, k, n, Ncelld, cellnow, ix, iy, iz, jx, jy, jz, kx, ky, kz;
  Nneigh=pow(3,dim);
  Ncell = floor(L/sigma);  //cells per side
  rc=L/Ncell;
  Ncelld=rint(pow(Ncell,dim));
  printf("rc=%f, Ncell=%d, sigma=%f\n",rc,Ncell,sigma);
  Nmax = 9; //maximum amount of particles we allow in a cell
  cell_list_index = ivector(0,N-1); //vector that links each particle to the index of its cell
  cell_list = imatrix(0,Ncelld-1,0,Nmax-1);  //matrix that for each cell tells the particles contained
  neighbor = imatrix(0,Ncelld-1,0,Nneigh-1);  //list of cells neighboring a certain cell, same cell included. For other purposes only half cells may be required
  
  //we initialize the cell_list matrix by filling each cell of Nmax "-1". In this way we will know to stop counting when we hit the "-1"
  
  for (i=0; i<Ncelld; i++)
  {
    for (j=0; j<Nmax; j++)
    {
      cell_list[i][j]=-1;
    }
  }
  
  
  //we fill out this matrix with real particles
  
  for (n=0; n<N; n++)
  {
    ix=floor((x[n]+0.5*L)/rc);
    iy=floor((y[n]+0.5*L)/rc);
    iz=floor((z[n]+0.5*L)/rc);
    cellnow=ix+iy*Ncell+iz*Ncell*Ncell;
    k=0;
    while (cell_list[cellnow][k]!=-1)
    {
      k++;
      if (k==Nmax)
      {
        printf("Too many particles in cell %d!\n", cellnow);
        break;
      }
    }
    cell_list[cellnow][k]=n;
    cell_list_index[n]=cellnow;
  }
  
  //we generate the list of the cells that neighbor each other
  
  for (ix=0; ix<Ncell; ix++)
  {
    for (iy=0; iy<Ncell; iy++)
    {
      for (iz=0; iz<Ncell; iz++)
      {
        n=0;
        for (jx=ix-1; jx<ix+2; jx++)
        {
          for (jy=iy-1; jy<iy+2; jy++)
          {
            for (jz=iz-1; jz<iz+2; jz++)
            {
              kx=jx;
              ky=jy;
              kz=jz;
              if (jx<0)
                kx=Ncell-1;
              if (jx==Ncell)
                kx=0;
              if (jy<0)
                ky=Ncell-1;
              if (jy==Ncell)
                ky=0;
              if (jz<0)
                kz=Ncell-1;
              if (jz==Ncell)
                kz=0;
              neighbor[ix+Ncell*iy+Ncell*Ncell*iz][n]=kx+ky*Ncell+kz*Ncell*Ncell;
              n++;
            }
          }
        }
      }
    }
  }
}

/*********************************************************************************************************************************/

void reinit_cell_list()
{
  int i, j, k, n, o, cellnow, Ncelld;
  
  Ncelld=rint(pow(Ncell,dim));
  //we initialize the cell_list matrix by filling each cell of Nmax "-1". In this way we will know to stop counting when we hit the "-1"
  
  for (i=0; i<Ncelld; i++)
  {
    for (j=0; j<Nmax; j++)
    {
      cell_list[i][j]=-1;
    }
  }
  
  //we fill out this matrix with real particles
  
  for (n=0; n<N; n++)
  {
    i=floor((x[n]+0.5*L)/rc);
    j=floor((y[n]+0.5*L)/rc);
    o=floor((z[n]+0.5*L)/rc);
    cellnow=i+j*Ncell+o*Ncell*Ncell;
    k=0;
    while (cell_list[cellnow][k]!=-1)
    {
      k++;
      if (k==Nmax)
      {
        printf("%d\n",n);
        printf("Too many particles in cell %d", cellnow);
        break;
      }
    }
    cell_list[cellnow][k]=n;
    cell_list_index[n]=cellnow;
  }
}

/*********************************************************************************************************************************/

void print_single_g(int step)
{
  double Lmax, dhist, A, a;
  int i, j, hidx;

  Lmax=L;
  dhist=sigma*0.5;
  
  for (i=0;i<Nhist;i++)
  {
    hist[i]=0;
  }
  
  //printf("lol");
  for (i=0; i<N; i++)
  {
    for (j=0; j<i; j++)
    {
      hidx=(int)(dist(x[i],y[i],z[i],x[j],y[j],z[j])/dhist);
      if (hidx<Nhist)
        hist[hidx]+=2;
    }
  }
  std::stringstream sstr;
  sstr << "g_step" << step << ".dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();

  FILE * out = fopen(cstr,"w");

  fprintf(out,"%f %f %f\n",L,sigma,rc);
  
  a=4./3.*PI*pow(dhist,dim);
  
  for (i=0; i<Nhist; i++) 
  {
    A = a* (pow(i+1,dim)-pow(i,dim));
    fprintf(out,"%f %f\n",i*dhist,hist[i]/((double) N * N *A));
  }
  
  fclose(out);
  
}

/*********************************************************************************************************************************/

void print_g(FILE* &gout)
{
  double Lmax, dhist, A, a;
  int i, j, hidx;

  Lmax=L;
  dhist=sigma*0.5;
  
  for (i=0;i<Nhist;i++)
  {
    hist[i]=0;
  }
  
  //printf("lol");
  for (i=0; i<N; i++)
  {
    for (j=0; j<i; j++)
    {
      hidx=(int)(dist(x[i],y[i],z[i],x[j],y[j],z[j])/dhist);
      if (hidx<Nhist)
        hist[hidx]+=2;
    }
  }
  
  a=4./3.*PI*pow(dhist,dim);
  
  for (i=0; i<Nhist; i++) 
  {
    A = a* (pow(i+1,dim)-pow(i,dim));
    fprintf(gout,"%f ",hist[i]/((double) N * N *A));
  }
  fprintf(gout,"\n");
  
}

/*********************************************************************************************************************************/

void print_sf(FILE* &sfout)
{
  double SF, q, N2;
  int i, j, k, Tot;
  
  Tot=20;
  N2=2./(3.*(double)N);
  
  for (k=1; k<=Tot; k++) //k=0 is pointless (gives the number of particles
  {
    SF=0.;
    q=2.*PI/L*k;
    
    for (i=0;i<N;i++)
    {
      for (j=0;j<i;j++)
      {
        SF+=cos(q*((x[i]-x[j])-L*rint((x[i]-x[j])/L)));
        SF+=cos(q*((y[i]-y[j])-L*rint((y[i]-y[j])/L)));
        SF+=cos(q*((z[i]-z[j])-L*rint((z[i]-z[j])/L)));
      }
    }
    fprintf(sfout,"%f ",1.+SF*N2);
  }
 
  fprintf(sfout,"\n");
  
}

/*********************************************************************************************************************************/

void init_block(double** &x_blk, double** &y_blk, double** &z_blk, double* &Zx, double* &Zy, double* &Zz)
{
  
  x_blk = dmatrix(0,N-1,0,Nblock*Nlevel-1);
  y_blk = dmatrix(0,N-1,0,Nblock*Nlevel-1);
  z_blk = dmatrix(0,N-1,0,Nblock*Nlevel-1);
  
  Zx = dvector(0,(Nblock-1)*Nlevel-1);
  Zy = dvector(0,(Nblock-1)*Nlevel-1);
  Zz = dvector(0,(Nblock-1)*Nlevel-1);
}

/*********************************************************************************************************************************/

void block_intg(int step, double** &x_blk, double** &y_blk, double** &z_blk, double* &Zx, double* &Zy, double* &Zz, FILE* &drout)
{
  int i, ib, il, q, index, index1, index2, low_level_start, low_level_end, il1, Zindex, Zindex0, c;
  double dt2N;
  bool BREAK;
  
  dt2N=dt*dt/((double)N);
  
  il=Nlevel;
  q=pow(Nblock,il);
  index1=q;
  index=(step+1)%q;
  while(index!=0&&q!=1) 
  {
    index1=index;
    il-=1;
    q/=Nblock;
    index=index%q;
  }
  ib=index1/q-1;
  if (ib<0)
    ib+=Nblock;
  //printf("step=%d, block=%d, level=%d\n",step,ib,il);
  
  
  index=il*Nblock+ib;
  
  
  if (il==0)  //simple velocity for il=0
  {
    for (i=0; i<N; i++)
    {
      x_blk[i][index]=vx[i];  //assign quantity to right slot
      y_blk[i][index]=vy[i];
      z_blk[i][index]=vz[i]; 
    }
  }
  else  //sum of lower velocities for il>0
  {
    BREAK=false;
    for (il1=0;il1<=il;il1++)
    {
      if(il1==il && (x_blk[0][(il1+1)*Nblock-1]==0||il1<3||il1==Nlevel))  //the last condition, paired with il<Nlevel with BREAK afterwards, is FUNDAMENTAL: avoids last level from messing up
      {
        BREAK=true;
        break;
      }
      
      low_level_end=(il1+1)*Nblock-1;
      low_level_start=il1*Nblock;
      
      Zindex0=il1*(Nblock-1);
      if (x_blk[0][(il1+1)*Nblock-1]==0||il1<3)
      {
        if (il1==0)
        {
          for (i=0;i<N;i++)
          {
            x_blk[i][low_level_end]=vx[i];  //assign "missed" value in the 0th level
            y_blk[i][low_level_end]=vy[i];
            z_blk[i][low_level_end]=vz[i];
            
            Zx[Zindex0]+=x_blk[i][low_level_end]*x_blk[i][low_level_end];  //calculate correlation for last element of 0th level. No need for next levels
            Zy[Zindex0]+=y_blk[i][low_level_end]*y_blk[i][low_level_end];
            Zz[Zindex0]+=z_blk[i][low_level_end]*z_blk[i][low_level_end];
          }
        }
        else
        {
          for (i=0;i<N;i++)
          {
            x_blk[i][low_level_end]=x_blk[i][low_level_start-Nblock];  //assign "missed" value in the other levels
            y_blk[i][low_level_end]=y_blk[i][low_level_start-Nblock];
            z_blk[i][low_level_end]=z_blk[i][low_level_start-Nblock];
          }
        }

        for (index1=low_level_end-1;index1>=low_level_start;index1--) //perform partial sum
        {
          
          Zindex=Zindex0+(low_level_end-index1);
          for (i=0;i<N;i++)
          {
            x_blk[i][index1]+=x_blk[i][index1+1];  
            y_blk[i][index1]+=y_blk[i][index1+1];
            z_blk[i][index1]+=z_blk[i][index1+1];
            
            
            Zx[Zindex]+=x_blk[i][index1]*x_blk[i][index1]; //calculate correlation for the other elements 
            Zy[Zindex]+=y_blk[i][index1]*y_blk[i][index1];  
            Zz[Zindex]+=z_blk[i][index1]*z_blk[i][index1];  
          }
        }
 
      }            
      else
      {
        low_level_end=(il1+1)*Nblock-1;
        low_level_start=il1*Nblock;
      
        Zindex0=il1*(Nblock-1);
        //printf("%d, %d\n",step, il1);
        for (i=0;i<N;i++)
        {
          for (index1=low_level_end-Nblock+1;index1<low_level_end;index1++) //for higher levels we only delete one term
          {
            x_blk[i][index1]=x_blk[i][index1+1]+x_blk[i][low_level_start-Nblock];
            y_blk[i][index1]=y_blk[i][index1+1]+y_blk[i][low_level_start-Nblock];
            z_blk[i][index1]=z_blk[i][index1+1]+z_blk[i][low_level_start-Nblock];
            
            Zindex=Zindex0+(low_level_end-index1);
            Zx[Zindex]+=x_blk[i][index1]*x_blk[i][index1];
            Zy[Zindex]+=y_blk[i][index1]*y_blk[i][index1];
            Zz[Zindex]+=z_blk[i][index1]*z_blk[i][index1];

          }
          
          x_blk[i][low_level_end]=x_blk[i][low_level_start-Nblock];
          y_blk[i][low_level_end]=y_blk[i][low_level_start-Nblock];
          z_blk[i][low_level_end]=z_blk[i][low_level_start-Nblock];

        }
      }

      
      if (il1==0)
      {
        fprintf(drout,"%d %.8f %.8f %.8f\n",Zindex0,Zx[Zindex0]*dt2N,Zy[Zindex0]*dt2N,Zz[Zindex0]*dt2N);
        Zx[Zindex0]=0.;
        Zy[Zindex0]=0.;
        Zz[Zindex0]=0.;
      }

      if (il1==Nlevel)
        c=1;
      else
        c=0;
      
      for (Zindex=Zindex0+1;Zindex<Zindex0+Nblock-c;Zindex++) //print correlation
      {
        fprintf(drout,"%d %.8f %.8f %.8f\n",Zindex,Zx[Zindex]*dt2N,Zy[Zindex]*dt2N,Zz[Zindex]*dt2N);
        Zx[Zindex]=0.;
        Zy[Zindex]=0.;
        Zz[Zindex]=0.;
      }
      
    }
    if (BREAK && il<Nlevel)
      for (i=0;i<N;i++)
      {
        x_blk[i][index]=x_blk[i][low_level_start];  //first element in last level is the new element of the current level
        y_blk[i][index]=y_blk[i][low_level_start];
        z_blk[i][index]=z_blk[i][low_level_start];
      }
  }
}

/*********************************************************************************************************************************/

void block_corr(int step, double** &x_blk, double** &y_blk, double** &z_blk, double* &Zx, double* &Zy, double* &Zz, FILE* &zout)
{
  int i, ib, il, q, index, index1, index2;
  double corrx0,corry0,corrz0;
  
  il=Nlevel;
  q=pow(Nblock,il);
  index=(step+1)%q;
  while(index!=0&&q!=1) 
  {
    index1=index;
    il-=1;
    q/=Nblock;
    index=index%q;
  }
  ib=index1/q-1;
  if (ib<0)
    ib+=Nblock;
  //printf("step=%d, block=%d, level=%d\n",step,ib,il);
  
  //if (il==Nlevel)
  //  il--;
  
  index=il*Nblock+ib;
  if (index<Nblock*Nlevel)
  {
    if (x_blk[0][(il+1)*Nblock-1]==0||il<3) //for lower levels we delete the whole level
    {
      for (i=0; i<N; i++)
      {
        x_blk[i][index]=vx[i];  //assign quantity to right slot
        y_blk[i][index]=vy[i];
        z_blk[i][index]=vz[i];
      }
    }
    else
    {
      index=(il+1)*Nblock-1;
      for (i=0; i<N; i++)
      {
        for (index1=index-Nblock+1;index1<index;index1++) //for higher levels we only delete one term
        {
          x_blk[i][index1]=x_blk[i][index1+1];
          y_blk[i][index1]=y_blk[i][index1+1];
          z_blk[i][index1]=z_blk[i][index1+1];
        }

      }
      il++; //will be useful for next section
    }
  }
  
  while (il>0)
  {
    index=il*Nblock-1;  //in this way we will be working from last block of old level, instead of working from first block on new level (it's better)
    
    corrx0=0.;
    corry0=0.;
    corrz0=0.;
    
    for (i=0; i<N; i++)
    {
      x_blk[i][index]=vx[i];  //assign "missed" quantity to first slot of new level when reaching new level (since ib=0 is never truly reached)
      y_blk[i][index]=vy[i];
      z_blk[i][index]=vz[i];
      
      corrx0+=x_blk[i][index-Nblock+1]*x_blk[i][index-Nblock+1];
      corry0+=y_blk[i][index-Nblock+1]*y_blk[i][index-Nblock+1];
      corrz0+=z_blk[i][index-Nblock+1]*z_blk[i][index-Nblock+1];
      
      for (index1=index-Nblock+2; index1<=index; index1++) //correlate over last level
      {
        index2=index1-il; //il still refers to new level even if index refers to old level
        Zx[index2]+=x_blk[i][index-Nblock+1]*x_blk[i][index1]; //The correlation functions contains one element less than Nblock per level bc we don't correlate a value with itself
        Zy[index2]+=y_blk[i][index-Nblock+1]*y_blk[i][index1];
        Zz[index2]+=z_blk[i][index-Nblock+1]*z_blk[i][index1];
      }
    }

    for (index1=index-Nblock-il+2; index1<=index-il; index1++) //print the new values of Z
    {
      fprintf(zout,"%d %f %f %f\n",index1,Zx[index1]/corrx0,Zy[index1]/corry0,Zz[index1]/corrz0);  //Z is a dumb array I will remove it
      Zx[index1]=0;
      Zy[index1]=0;
      Zz[index1]=0;
    }
    il--;
  }
  
}

/*********************************************************************************************************************************/

void init_output(FILE* &fout,FILE* &gout,FILE* &sfout, FILE* &zout, FILE* &drout) //here I pass by reference (&) a FILE* type (i.e. pointer to FILE), that allows me to treat the FILE* variable as a global one 
{
 
  fout = fopen("output.dat","w");

  gout = fopen("g_total.dat","w");
  
  Nhist=ceil(L/sigma); //for the g function histogram
  hist=ivector(0,Nhist-1);

  fprintf(gout,"%f\n", sigma*0.5);
  
  sfout = fopen("sf_total.dat","w");
  
  fprintf(sfout,"%f\n", 2.*PI/L);
  
  zout = fopen("Z_corr.dat","w");
  
  fprintf(zout,"%d %d %f\n", Nblock, Nlevel, dt);
  
  drout = fopen("Dr.dat","w");
  
  fprintf(drout,"%d %d %f\n", Nblock, Nlevel, dt);
  
  //print parameters
  FILE * parameters = fopen("info.dat","w");
  fprintf(parameters,"L N T dt sigma rho dim\n"); 
  fprintf(parameters,"%f %d %f %f %f %f %d", L, N, T, dt, sigma, rho, dim);
  fclose(parameters); 
}

/*********************************************************************************************************************************/

void close_output(FILE* &gout, FILE* &fout, FILE* &sfout, FILE* &zout, FILE* &drout)
{
  fclose(gout);
  fclose(fout);
  fclose(sfout);
  fclose(zout);
  fclose(drout);
}

/*********************************************************************************************************************************/

void print_cell(int step)
{
  int i, j, Ncelld;
  std::stringstream sstr;
  sstr << "cell_step" << step << ".dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();
  
  FILE * out = fopen(cstr,"w");
  
  fprintf(out,"%f %f %f\n",L,sigma,rc);
  
  Ncelld=rint(pow(Ncell,dim));
  for (i=0;i<Ncelld;i++)
  {
    fprintf(out,"%d ",i);
    for (j=0;j<Nmax;j++)
    {
      fprintf(out,"%d %",cell_list[i][j]);
    }
    fprintf(out,"\n");
  }
  
  for (i=0;i<Ncelld;i++)
  {
    fprintf(out,"%d ",i);
    for (j=0;j<Nneigh;j++)
    {
      fprintf(out,"%d %",neighbor[i][j]);
    }
    fprintf(out,"\n");
  }
  
  fprintf(out,"\n");
  for (i=0;i<N;i++)
    fprintf(out,"%d, %d \n%",i, cell_list_index[i]);
  fclose(out);
}

/*********************************************************************************************************************************/

void test_cell(int step)
{
  int j, k, i, Ncelld;
  bool trigger;
  Ncelld=rint(pow(Ncell,dim));
  for (i=0;i<N;i++)
  {
    trigger=false;
    for (j=0; j<Ncelld;j++)
      for (k=0; k<Nmax;k++)
        if (cell_list[j][k]==i)
          trigger=true;
    
    if (!trigger)
      printf("WTF! Step %d, cell %d, particle %d\n",step, i, j);
  }
}

/*********************************************************************************************************************************/

void print_block(double** &x_blk, double** &y_blk, double** &z_blk)
{
  int i, j, k, index;
  
  FILE* bout = fopen("block.dat", "w");
  
  
  for (i=0; i<N; i++)
  {
    for (j=0; j<Nlevel; j++)
    {
      for (k=0; k<Nblock; k++)
      {
        index=j*Nblock+k;
        fprintf(bout, "%f %f %f %d %d %d\n", x_blk[i][index], y_blk[i][index], z_blk[i][index], i, j, k);
      }
    }
  }
  
  fclose(bout);
}
  
/*********************************************************************************************************************************/