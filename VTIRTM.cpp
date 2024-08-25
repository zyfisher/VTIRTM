#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<windows.h>
#define pi  3.1415926
#define f0  15.0
#define DX  10
#define DZ  10
#define NX  441 
#define NZ  241
//#define MX  175
//#define MZ  20
//#define MZ1 20
#define PML 20  /*PML层数*/
#define SNAPTIME 1.5
#define SNAPTIME1 0
#define DT 0.001
//#define NT int(SNAPTIME/DT)
//#define NT1 int(SNAPTIME1/DT)
#define NT 3000
#define NT1 0
#define R        0.000001/*理论反射系数*/
#define Niterations     50/*迭代次数*/
#define MX1       0/*第一炮炮点*/
#define MX2       8/*炮点距*/
int tomo(int a1)
{
  return int(10.0*sin(2*pi*(a1-25-PML)/100.0)+10.0+PML);/*起伏地表函数可修改*/
}
float **creat2Ddata(int a1,int b)
{
    int I;
    float **d;
    d=new float*[a1]; 
	for(I=0;I<a1;I++)
		d[I]=new float[b];
    return d;
}
float ***creat3Ddata(int a1,int b)
{  
    int J,K;
    float ***g;
    g=new float**[2];
    for(J=0;J<2;J++)
		g[J]=new float*[a1];
     for(J=0;J<2;J++)
	{
		for(K=0;K<a1;K++)
			g[J][K]=new float[b];
	}
    return g;
}
void fft(float *xr,float *xi,int nr,float sign) 
{
	float tr,ti,treal,timag;
	int i,p,n,m,k,v,r,w,q,n1,N1;
	int j=0,u=0;
	
    N1=pow(2,nr);  
	for(i=1;i<(N1-1);i++)
	{
	   for(p=1;p<=nr;p++)
	    {
		 k=(int)(N1/pow(2,p));
		 if(j<k)
		 {
		  j=j+k;
		  break;
		 }
		 else
		 {
		 j=j-k;
		 }
	   }
	   if(i<j)
	   {
	    tr=xr[i];
	    ti=xi[i];
	    xr[i]=xr[j];
	    xi[i]=xi[j];
	    xr[j]=tr;
	    xi[j]=ti;
	    }
	}
	for(i=nr;i>=1;i--)
	{
	  m=0;
	  v=(int)pow(2,i-1);
	  w=(int)(N1/(v*2));
	  r=(int)pow(2,u);
	  q=2*r;
	 for(j=1;j<=v;j++)
	  {
	   for(n=1;n<=w;n++)
	   {
	   treal=xr[m];
	   timag=xi[m];
	   xr[m]=treal+xr[m+r]*cos(2.0*pi*m/q)+sign*xi[m+r]*sin(2.0*pi*m/q);
	   xi[m]=timag-sign*xr[m+r]*sin(2.0*pi*m/q)+xi[m+r]*cos(2.0*pi*m/q);
	   xr[m+r]=2.0*treal-xr[m];
	   xi[m+r]=2.0*timag-xi[m];
	   m=m+1;
	   }
	   m=m+r;
	  } 
	 u=u+1;
	}
	
	for(i=0;i<N1;i++)
	{
		if(sign<0.0)
	    {
		xr[i]=xr[i]*(1.0/(N1));
		xi[i]=xi[i]*(1.0/(N1));
		}
		else
		{
		xr[i]=xr[i];
		xi[i]=xi[i];	
	    }
	}
}
void  filter(float *wavelet,int nt)
{
    float W;
	int I,J,K,j,k1,k2,nfft1,nfft2,NK,L;
    float *xr,*xi;
	/*计算fft点数*/
	L=nt+1;
	k2=log(L)/log(2);
	if(L>pow(2,k2))k2=k2+1;
	nfft2=pow(2,k2);
	
	//printf("nfft2=%d k2=%d\n",nfft2,k2);
	
    xr=new float[nfft2]; 
    xi=new float[nfft2];
    
	for(K=0;K<nfft2;K++)
	   {
	    xr[K]=0.0;
	    xi[K]=0.0;
	   }	
	for(K=0;K<=nt;K++)
	   {
	    xr[K]=wavelet[K];
	    xi[K]=0.0;
	   }
	fft(xr,xi,k2,1.0);
	for(K=0;K<nfft2;K++)
	   {
		   W=(K+0.1)*(2.0*pi/(DT*nfft2));
		   xr[K]=xr[K]/(W*W);
		   xi[K]=xi[K]/(W*W);
	   }
    fft(xr,xi,k2,-1.0);
    for(K=0;K<=nt;K++)
	    {
	     wavelet[K]=xr[K];
	    }
	 delete []xr;
     delete []xi;
}
/*void  smooth(float **g,float **g1, int nx, int nz, int pml, int length)
{
  int I,J,K,j,k,nx1,nz1;
  float **g2,**sum1,m,**sum2;
  nx1=nx-2*pml;
  nz1=nz-2*pml;
  m=float(length);
  
    g2=new float*[nx1]; 
	for(I=0;I<nx1;I++)
		g2[I]=new float[nz1];

    sum1=new float*[nx1]; 
	for(I=0;I<nx1;I++)
		sum1[I]=new float[nz1];

    sum2=new float*[nx1]; 
	for(I=0;I<nx1;I++)
		sum2[I]=new float[nz1];
    
    for(J=0;J<nx1;J++) 
      for(K=0;K<nz1;K++)
	  {
	   g2[J][K]=g[J+pml][K+pml];
	  }

   for(I=0;I<1;I++)
   {
     for(J=0;J<nx1;J++)
		for(K=0;K<nz1;K++)
		{
           sum1[J][K]=0.0;
           if((K>=length/2)&&(K<=nz1-(length/2)))
		   {
			  for(k=-length/2;k<length/2;k++)
			  {
               sum1[J][K]=sum1[J][K]+g2[J][K+k];
			  }
             g2[J][K]=sum1[J][K]/m;
		   } 
           else
		   {
             g2[J][K]=g[J][K];
		   }
       }
      for(J=0;J<nx1;J++)
		for(K=0;K<nz1;K++)
		{
           sum2[J][K]=0.0;
           if((J>=length/2)&&(J<=nx1-(length/2)))
		   {
			  for(j=-length/2;j<length/2;j++)
			  {
               sum2[J][K]=sum2[J][K]+g2[J+j][K];
			  }
             g2[J][K]=sum2[J][K]/m;
		   }
           else
		   {
             g2[J][K]=g2[J][K]; 
		   }  
       }
	}
	for(J=0;J<nx1;J++) 
      for(K=0;K<nz1;K++)
	  {
	   g1[J+pml][K+pml]=g2[J][K];
	  }
}*/
void  smooth(float **g,float **g1, int nx, int nz, int pml, int length)
{
  int I,J,K,j,k,nx1,nz1;
  float **g2,**sum1,m,**sum2;
  //nx1=nx-2*pml;
  //nz1=nz-2*pml;
  nx1=nx;
  nz1=nz;
  m=float(length);
  
    g2=new float*[nx1]; 
	for(I=0;I<nx1;I++)
		g2[I]=new float[nz1];

    sum1=new float*[nx1]; 
	for(I=0;I<nx1;I++)
		sum1[I]=new float[nz1];

    sum2=new float*[nx1]; 
	for(I=0;I<nx1;I++)
		sum2[I]=new float[nz1];
    
    for(J=0;J<nx1;J++) 
      for(K=0;K<nz1;K++)
	  {
	   g2[J][K]=g[J][K];
	  }

   for(I=0;I<1;I++)
   {
     for(J=0;J<nx1;J++)
		for(K=0;K<nz1;K++)
		{
           sum1[J][K]=0.0;
           if((K>=length/2)&&(K<=nz1-(length/2)))
		   {
			  for(k=-length/2;k<length/2;k++)
			  {
               sum1[J][K]=sum1[J][K]+g2[J][K+k];
			  }
             g2[J][K]=sum1[J][K]/m;
		   } 
           else
		   {
             g2[J][K]=g[J][K];
		   }
       }
      for(J=0;J<nx1;J++)
		for(K=0;K<nz1;K++)
		{
           sum2[J][K]=0.0;
           if((J>=length/2)&&(J<=nx1-(length/2)))
		   {
			  for(j=-length/2;j<length/2;j++)
			  {
               sum2[J][K]=sum2[J][K]+g2[J+j][K];
			  }
             g2[J][K]=sum2[J][K]/m;
		   }
           else
		   {
             g2[J][K]=g2[J][K]; 
		   }  
       }
	}
	for(J=0;J<nx1;J++) 
      for(K=0;K<nz1;K++)
	  {
	   g1[J][K]=g2[J][K];
	  }
}
void main()
{
    int I,J,K,i,j,k,MX,MZ,MZ1,N;
	/*差分系数*/
	float A1=1.22134,A2=-0.0969315,A3=0.0174477,A4=-0.00296729,A5=0.000359005,A6=-0.0000218478;	
	float Sum1,Sum2,Sum21,Sum22,Sum41,Sum42,Sum51,Sum52,Sum61,Sum62,g2max,g3max,max,max1,max2,max3,min1,min2,min3,C1max,C2max,C3max,MAX1,MAX2,MAX3,step,Sum1g,Sum0g;
	float A,B,C,D,E,M,step1,step2,step3,L2,L4,L6,beta,Mpmax,Mpmin,Msmax,Msmin,tao;
	
    float ***U,***UX,***UZ,***Up,***UpX,***UpZ,***Us,***UsX,***UsZ,***V,***VX,***VZ,***Vp,***VpX,***VpZ,***Vs,***VsX,***VsZ,***P,***PX,***PZ,***Q,***QX,***QZ,***S,***SX,***SZ,***tH,***tHX,***tHZ,***tV,***tVX,***tVZ;
    float ***U2,***U2X,***U2Z,***U2p,***U2s,***V2,***V2X,***V2Z,***V2p,***V2s,***P2,***P2X,***P2Z,***Q2,***Q2X,***Q2Z,***S2,***S2X,***S2Z,***tH2,***tH2X,***tH2Z,***tV2,***tV2X,***tV2Z;
    float ***U4,***U4X,***U4Z,***U4p,***U4pX,***U4pZ,***U4s,***U4sX,***U4sZ,***V4,***V4X,***V4Z,***V4p,***V4pX,***V4pZ,***V4s,***V4sX,***V4sZ,***P4,***P4X,***P4Z,***Q4,***Q4X,***Q4Z,***S4,***S4X,***S4Z,***tH4,***tH4X,***tH4Z,***tV4,***tV4X,***tV4Z;
    float **sum1,**sum2,**sum3,**sum4,**sum5,**sum6,**sum7,**sum8,**sum9,**sum10,**SUM1,**SUM2,**SUM3,**SUM4,**SUM5,**SUM6,**SUM7,**SUM8,**SUM9,**SUM10;
    float **Image1,**Image2,**Image3,**Image4,**Image5,**Image6,**Image7,**Image8,**Image9,**Image10;
    float **I1,**I2,**I3,**I4,**I5,**I6,**I7,**I8,**I9,**I10,**W1,**W2,**W3,**W4,**W5,**W6,**W7,**W8,**W9,**W10,**filter,**filter1;
    float **tsum_U2,**tsum_U2p,**tsum_U2s,**tsum_V2,**tsum_V2p,**tsum_V2s,**tsum_U4,**tsum_U4p,**tsum_U4s,**tsum_V4,**tsum_V4p,**tsum_V4s,**U2p1,**V2p1,**U4p1,**U4s1,**V4p1,**V4s1;
    
    float **UsX_obs,**UsZ_obs;
    	
    float **VP0,**VS0,**VP,**VS,**c11,**c33,**c13,**c44,**C11,**C33,**C13,**C44,**C55,**C15,**C35,**DEN,**e,**r,**det,**dx,**dz,**e0,**det0,**DEN0;
    float ***UU2p,***UU2s,***VV2p,***VV2s,***UU2,***VV2,***UU4p,***UU4s,***VV4p,***VV4s,***UU4,***VV4;
    float **D11p,**D12p,**D13p,**D21p,**D22p,**D23p,**D11s,**D12s,**D13s,**D21s,**D22s,**D23s;
    
    
    float **UsX1,**UsX2,**UsX3,**UsX4,**UsX5,**UsX6,**UsX7,**UsX8,**UsX9,**UsX10,**UsX11,**UsX12,**UsX13,**UsX14,**UsX15,**UsX16,**UsX17,**UsX18,**UsX19,**UsX20,**UsX21,**UsX22,**UsX23,**UsX24;
    float **UsZ1,**UsZ2,**UsZ3,**UsZ4,**UsZ5,**UsZ6,**UsZ7,**UsZ8,**UsZ9,**UsZ10,**UsZ11,**UsZ12,**UsZ13,**UsZ14,**UsZ15,**UsZ16,**UsZ17,**UsZ18,**UsZ19,**UsZ20,**UsZ21,**UsZ22,**UsZ23,**UsZ24;
    float wavelet[NT+1],w1,W,POW;
    float misfit[Niterations];
	
    void PML_coefficient(float **dx,float **dz);
    void parameters(float **VP);
    void laplace_filter(float **I,float **I1);
    
    void foward(int I,int MX,int MZ,float A1,float A2,float A3,float A4,float A5,float A6,float **C11,float **C13,float **C33,float **C44,float **c13,float **DEN,float **dx, float **dz,float ***P,float ***PX,float ***PZ,float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,float ***tH,float ***tHX,float ***tHZ,float ***tV,float ***tVX,float ***tVZ,float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs);
    void reconstruct(float A1,float A2,float A3,float A4,float A5,float A6,float **C11,float **C13,float **C33,float **C44,float **c13,float **DEN,float ***P,float ***Q,float ***S,float ***tH,float ***tV,float ***U,float ***Up,float ***Us,float ***V,float ***Vp,float ***Vs);
    void adjoint(float A1,float A2,float A3,float A4,float A5,float A6,float **C11,float **C13,float **C33,float **C44,float **c13,float **DEN,float **dx, float **dz,float ***P,float ***PX,float ***PZ,float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,float ***tH,float ***tHX,float ***tHZ,float ***tV,float ***tVX,float ***tVZ,float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs);
    
    void foward1(int I,int MX,int MZ,float A1,float A2,float A3,float A4,float A5,float A6,
             float **C11,float **C13,float **C33,float **C44,float **D11p,float **D12p,float **D13p,float **D21p,float **D22p,float **D23p,
             float **D11s,float **D12s,float **D13s,float **D21s,float **D22s,float **D23s,float **DEN,float **dx, float **dz,
             float ***P,float ***PX,float ***PZ,float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,
             float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***UsX,float ***UsZ,
             float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs,float ***VsX,float ***VsZ);
             
    void reconstruct1(float A1,float A2,float A3,float A4,float A5,float A6,
                  float **C11,float **C13,float **C33,float **C44,
                  float **D11p,float **D12p,float **D13p,float **D21p,float **D22p,float **D23p,
                  float **D11s,float **D12s,float **D13s,float **D21s,float **D22s,float **D23s,
                  float **DEN,
                  float ***P,float ***Q,float ***S,float ***U,float ***Up,float ***Us,float ***V,float ***Vp,float ***Vs);
                  
    void adjoint1(float A1,float A2,float A3,float A4,float A5,float A6,
               float **C11,float **C13,float **C33,float **C44,
               float **D11p,float **D12p,float **D13p,float **D21p,float **D22p,float **D23p,
               float **D11s,float **D12s,float **D13s,float **D21s,float **D22s,float **D23s,float **DEN,
               float **dx, float **dz,float ***P,float ***PX,float ***PZ,
               float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,
               float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***UsX,float ***UsZ,
               float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs,float ***VsX,float ***VsZ);
                                      
   /*开辟动态数组*/
    U=creat3Ddata(NX,NZ);
    UX=creat3Ddata(NX,NZ);
    UZ=creat3Ddata(NX,NZ);
    
    Up=creat3Ddata(NX,NZ);
    UpX=creat3Ddata(NX,NZ);
    UpZ=creat3Ddata(NX,NZ);
    Us=creat3Ddata(NX,NZ);
    UsX=creat3Ddata(NX,NZ);
    UsZ=creat3Ddata(NX,NZ);
    
    U2=creat3Ddata(NX,NZ);
    U2X=creat3Ddata(NX,NZ);
    U2Z=creat3Ddata(NX,NZ);
    
	U2p=creat3Ddata(NX,NZ);
    U2s=creat3Ddata(NX,NZ);

	U4=creat3Ddata(NX,NZ);
    U4X=creat3Ddata(NX,NZ);
    U4Z=creat3Ddata(NX,NZ);
    
    U4p=creat3Ddata(NX,NZ);
    U4pX=creat3Ddata(NX,NZ);
    U4pZ=creat3Ddata(NX,NZ);
    U4s=creat3Ddata(NX,NZ);
    U4sX=creat3Ddata(NX,NZ);
    U4sZ=creat3Ddata(NX,NZ);
    
    UU2=creat3Ddata(NX,NZ);
    VV2=creat3Ddata(NX,NZ);
    UU2p=creat3Ddata(NX,NZ);
    UU2s=creat3Ddata(NX,NZ);
    VV2p=creat3Ddata(NX,NZ);
    VV2s=creat3Ddata(NX,NZ);
    
    V=creat3Ddata(NX,NZ);
    VX=creat3Ddata(NX,NZ);
    VZ=creat3Ddata(NX,NZ);
    
    Vp=creat3Ddata(NX,NZ);
    VpX=creat3Ddata(NX,NZ);
    VpZ=creat3Ddata(NX,NZ);
    Vs=creat3Ddata(NX,NZ);
    VsX=creat3Ddata(NX,NZ);
    VsZ=creat3Ddata(NX,NZ);
    
    V2=creat3Ddata(NX,NZ);
    V2X=creat3Ddata(NX,NZ);
    V2Z=creat3Ddata(NX,NZ);
    
	V2p=creat3Ddata(NX,NZ);
    V2s=creat3Ddata(NX,NZ);

	V4=creat3Ddata(NX,NZ);
    V4X=creat3Ddata(NX,NZ);
    V4Z=creat3Ddata(NX,NZ);
    
    V4p=creat3Ddata(NX,NZ);
    V4pX=creat3Ddata(NX,NZ);
    V4pZ=creat3Ddata(NX,NZ);
    V4s=creat3Ddata(NX,NZ);
    V4sX=creat3Ddata(NX,NZ);
    V4sZ=creat3Ddata(NX,NZ);
    UU4=creat3Ddata(NX,NZ);
    VV4=creat3Ddata(NX,NZ);
    UU4p=creat3Ddata(NX,NZ);
    UU4s=creat3Ddata(NX,NZ);
    VV4p=creat3Ddata(NX,NZ);
    VV4s=creat3Ddata(NX,NZ);
    
    P=creat3Ddata(NX,NZ);
    PX=creat3Ddata(NX,NZ);
    PZ=creat3Ddata(NX,NZ);

    P2=creat3Ddata(NX,NZ);
    P2X=creat3Ddata(NX,NZ);
    P2Z=creat3Ddata(NX,NZ);

    P4=creat3Ddata(NX,NZ);
    P4X=creat3Ddata(NX,NZ);
    P4Z=creat3Ddata(NX,NZ);
    
    Q=creat3Ddata(NX,NZ);
    QX=creat3Ddata(NX,NZ);
    QZ=creat3Ddata(NX,NZ);

    Q2=creat3Ddata(NX,NZ);
    Q2X=creat3Ddata(NX,NZ);
    Q2Z=creat3Ddata(NX,NZ);

    Q4=creat3Ddata(NX,NZ);
    Q4X=creat3Ddata(NX,NZ);
    Q4Z=creat3Ddata(NX,NZ);
    
    S=creat3Ddata(NX,NZ);
    SX=creat3Ddata(NX,NZ);
    SZ=creat3Ddata(NX,NZ);

    S2=creat3Ddata(NX,NZ);
    S2X=creat3Ddata(NX,NZ);
    S2Z=creat3Ddata(NX,NZ);

    S4=creat3Ddata(NX,NZ);
    S4X=creat3Ddata(NX,NZ);
    S4Z=creat3Ddata(NX,NZ);
    
    tH=creat3Ddata(NX,NZ);
    tHX=creat3Ddata(NX,NZ);
    tHZ=creat3Ddata(NX,NZ);

    tH2=creat3Ddata(NX,NZ);
    tH2X=creat3Ddata(NX,NZ);
    tH2Z=creat3Ddata(NX,NZ);

    tH4=creat3Ddata(NX,NZ);
    tH4X=creat3Ddata(NX,NZ);
    tH4Z=creat3Ddata(NX,NZ);
    
    tV=creat3Ddata(NX,NZ);
    tVX=creat3Ddata(NX,NZ);
    tVZ=creat3Ddata(NX,NZ);

    tV2=creat3Ddata(NX,NZ);
    tV2X=creat3Ddata(NX,NZ);
    tV2Z=creat3Ddata(NX,NZ);

    tV4=creat3Ddata(NX,NZ);
    tV4X=creat3Ddata(NX,NZ);
    tV4Z=creat3Ddata(NX,NZ);
    
    
    U2p1=creat2Ddata(NX,NZ);
    V2p1=creat2Ddata(NX,NZ);
    U4p1=creat2Ddata(NX,NZ);
    U4s1=creat2Ddata(NX,NZ);
    V4p1=creat2Ddata(NX,NZ);
    V4s1=creat2Ddata(NX,NZ);
    
    sum1=creat2Ddata(NX,NZ);
    sum2=creat2Ddata(NX,NZ);
    sum3=creat2Ddata(NX,NZ);
    sum4=creat2Ddata(NX,NZ);
    sum5=creat2Ddata(NX,NZ);
    sum6=creat2Ddata(NX,NZ);
    sum7=creat2Ddata(NX,NZ);
    sum8=creat2Ddata(NX,NZ);
    sum9=creat2Ddata(NX,NZ);
    sum10=creat2Ddata(NX,NZ);
    
    SUM1=creat2Ddata(NX,NZ);
	SUM2=creat2Ddata(NX,NZ);
	SUM3=creat2Ddata(NX,NZ);
	SUM4=creat2Ddata(NX,NZ);
	SUM5=creat2Ddata(NX,NZ);
	SUM6=creat2Ddata(NX,NZ);
	SUM7=creat2Ddata(NX,NZ);
	SUM8=creat2Ddata(NX,NZ);
	SUM9=creat2Ddata(NX,NZ);
	SUM10=creat2Ddata(NX,NZ);
	
    W1=creat2Ddata(NX,NZ);
    W2=creat2Ddata(NX,NZ);
    W3=creat2Ddata(NX,NZ);
    W4=creat2Ddata(NX,NZ);
    W5=creat2Ddata(NX,NZ);
    W6=creat2Ddata(NX,NZ);
    W7=creat2Ddata(NX,NZ);
    W8=creat2Ddata(NX,NZ);
    W9=creat2Ddata(NX,NZ);
    W10=creat2Ddata(NX,NZ);
    
    filter=creat2Ddata(NX,NZ);
    filter1=creat2Ddata(NX,NZ);
    
    Image1=creat2Ddata(NX,NZ);
    Image2=creat2Ddata(NX,NZ);
    Image3=creat2Ddata(NX,NZ);
    Image4=creat2Ddata(NX,NZ);
    Image5=creat2Ddata(NX,NZ);
    Image6=creat2Ddata(NX,NZ);
    Image7=creat2Ddata(NX,NZ);
    Image8=creat2Ddata(NX,NZ);
    Image9=creat2Ddata(NX,NZ);
    Image10=creat2Ddata(NX,NZ);
    
    I1=creat2Ddata(NX,NZ);
    I2=creat2Ddata(NX,NZ);
    I3=creat2Ddata(NX,NZ);
    I4=creat2Ddata(NX,NZ);
    I5=creat2Ddata(NX,NZ);
    I6=creat2Ddata(NX,NZ);
    I7=creat2Ddata(NX,NZ);
    I8=creat2Ddata(NX,NZ);
    I9=creat2Ddata(NX,NZ);
    I10=creat2Ddata(NX,NZ);
    
    VP0=creat2Ddata(NX,NZ);
    VS0=creat2Ddata(NX,NZ);
    VP=creat2Ddata(NX,NZ);
    VS=creat2Ddata(NX,NZ);
    DEN=creat2Ddata(NX,NZ);
    DEN0=creat2Ddata(NX,NZ);
    
    c11=creat2Ddata(NX,NZ);
    c33=creat2Ddata(NX,NZ);
    c13=creat2Ddata(NX,NZ);
    c44=creat2Ddata(NX,NZ);
    C11=creat2Ddata(NX,NZ);
    C33=creat2Ddata(NX,NZ);
    C13=creat2Ddata(NX,NZ);
    C44=creat2Ddata(NX,NZ);
    C55=creat2Ddata(NX,NZ);
    C15=creat2Ddata(NX,NZ);
    C35=creat2Ddata(NX,NZ);
    e=creat2Ddata(NX,NZ);
    r=creat2Ddata(NX,NZ);
    det=creat2Ddata(NX,NZ);
    e0=creat2Ddata(NX,NZ);
    det0=creat2Ddata(NX,NZ);
    
    dx=creat2Ddata(NX,NZ);
    dz=creat2Ddata(NX,NZ);
    D11p=creat2Ddata(NX,NZ);
    D12p=creat2Ddata(NX,NZ);
    D13p=creat2Ddata(NX,NZ);
    D21p=creat2Ddata(NX,NZ);
    D22p=creat2Ddata(NX,NZ);
    D23p=creat2Ddata(NX,NZ);
    
    D11s=creat2Ddata(NX,NZ);
    D12s=creat2Ddata(NX,NZ);
    D13s=creat2Ddata(NX,NZ);
    D21s=creat2Ddata(NX,NZ);
    D22s=creat2Ddata(NX,NZ);
    D23s=creat2Ddata(NX,NZ);
    
  
    UsX1=creat2Ddata(NX,NT+1);
    UsX2=creat2Ddata(NX,NT+1);
    UsX3=creat2Ddata(NX,NT+1);
    UsX4=creat2Ddata(NX,NT+1);
    UsX5=creat2Ddata(NX,NT+1);
    UsX6=creat2Ddata(NX,NT+1);
    UsX7=creat2Ddata(NX,NT+1);
    UsX8=creat2Ddata(NX,NT+1); 
	UsX9=creat2Ddata(NX,NT+1);
    UsX10=creat2Ddata(NX,NT+1);
    UsX11=creat2Ddata(NX,NT+1);
    UsX12=creat2Ddata(NX,NT+1);
    
    UsX13=creat2Ddata(NZ,NT+1);
    UsX14=creat2Ddata(NZ,NT+1);
    UsX15=creat2Ddata(NZ,NT+1);
    UsX16=creat2Ddata(NZ,NT+1);
    UsX17=creat2Ddata(NZ,NT+1);
    UsX18=creat2Ddata(NZ,NT+1);
    UsX19=creat2Ddata(NZ,NT+1);
    UsX20=creat2Ddata(NZ,NT+1);
    UsX21=creat2Ddata(NZ,NT+1);
    UsX22=creat2Ddata(NZ,NT+1);
    UsX23=creat2Ddata(NZ,NT+1);
    UsX24=creat2Ddata(NZ,NT+1);
    
    
    UsZ1=creat2Ddata(NX,NT+1);
    UsZ2=creat2Ddata(NX,NT+1);
    UsZ3=creat2Ddata(NX,NT+1);
    UsZ4=creat2Ddata(NX,NT+1);
    UsZ5=creat2Ddata(NX,NT+1);
    UsZ6=creat2Ddata(NX,NT+1);
    UsZ7=creat2Ddata(NX,NT+1);
    UsZ8=creat2Ddata(NX,NT+1); 
	UsZ9=creat2Ddata(NX,NT+1);
    UsZ10=creat2Ddata(NX,NT+1);
    UsZ11=creat2Ddata(NX,NT+1);
    UsZ12=creat2Ddata(NX,NT+1);
    
    
    UsZ13=creat2Ddata(NZ,NT+1);
    UsZ14=creat2Ddata(NZ,NT+1);
    UsZ15=creat2Ddata(NZ,NT+1);
    UsZ16=creat2Ddata(NZ,NT+1);
    UsZ17=creat2Ddata(NZ,NT+1);
    UsZ18=creat2Ddata(NZ,NT+1);
    UsZ19=creat2Ddata(NZ,NT+1);
    UsZ20=creat2Ddata(NZ,NT+1);
    UsZ21=creat2Ddata(NZ,NT+1);
    UsZ22=creat2Ddata(NZ,NT+1);
    UsZ23=creat2Ddata(NZ,NT+1);
    UsZ24=creat2Ddata(NZ,NT+1);
    
    UsX_obs=creat2Ddata(NX,NT+1);
    UsZ_obs=creat2Ddata(NX,NT+1);
    
    tsum_U2=creat2Ddata(NX,NZ);
    tsum_U2p=creat2Ddata(NX,NZ);
    tsum_U2s=creat2Ddata(NX,NZ);
    tsum_V2=creat2Ddata(NX,NZ);
    tsum_V2p=creat2Ddata(NX,NZ);
    tsum_V2s=creat2Ddata(NX,NZ);
    
    tsum_U4=creat2Ddata(NX,NZ);
    tsum_U4p=creat2Ddata(NX,NZ);
    tsum_U4s=creat2Ddata(NX,NZ);
    tsum_V4=creat2Ddata(NX,NZ);
    tsum_V4p=creat2Ddata(NX,NZ);
    tsum_V4s=creat2Ddata(NX,NZ);

   
	float LS1, LS2,LS3,LS4,LS5, LS6,LS7,LS8,LS9,LS10,LS11,LS12,LS13,LS14,LS15,LS16,LS17,LS18,LS19,LS20;
	FILE *f1,*f2,*f3,*f4,*f5;
	FILE *f6[NX],*f7[NX],*f8[NX],*f9[NX],*f10[NX];
	FILE *f11[NT+1],*f12[NT+1];
    char filename1[1024],filename2[1024],filename3[1024],filename4[1024],filename5[1024],filename6[1024],filename7[1024],filename8[1024],filename9[1024],filename10[1024],filename11[1024],filename12[1024];
	for(N=0;N<Niterations;N++)
	{
      misfit[N]=0.0;
	}
   for(I=0;I<=NT;I++)
	{
      wavelet[I]=0.0;
	}
    for(J=0;J<NX;J++)
		for(K=0;K<NZ;K++)
		{
         VP0[J][K]=0.0;
		 VS0[J][K]=0.0;
         VP[J][K]=0.0;
		 VS[J][K]=0.0;
		 DEN[J][K]=0.0;
		 det0[J][K]=0.0;
		 e0[J][K]=0.0;
         det[J][K]=0.0;
		 e[J][K]=0.0;
         DEN0[J][K]=0.0;
		}
         	  
    /*读取计算区域模型参数*/
	  
	  f1=fopen("VP2.dat","rb");
      f2=fopen("VS2.dat","rb");
	  for(J=PML;J<NX-PML;J++)
		 for(K=PML;K<NZ-PML;K++)
		 {
           fread(&LS1,sizeof(float),1,f1);
           fread(&LS2,sizeof(float),1,f2);
	       VP0[J][K]=LS1;
           VS0[J][K]=LS2;
		 }
	    fclose(f1);
        fclose(f2);
      f1=fopen("delta2.dat","rb");
      f2=fopen("epsilon2.dat","rb");
	  for(J=PML;J<NX-PML;J++)
		 for(K=PML;K<NZ-PML;K++)
		 {
           fread(&LS1,sizeof(float),1,f1);
           fread(&LS2,sizeof(float),1,f2);
	       det0[J][K]=LS1;
           e0[J][K]=LS2;
		 }
	    fclose(f1);
        fclose(f2);
	f1=fopen("DEN2.dat","rb");
	  for(J=PML;J<NX-PML;J++)
		 for(K=PML;K<NZ-PML;K++)
		 {
           fread(&LS1,sizeof(float),1,f1);
	       DEN0[J][K]=LS1;
		 }
	    fclose(f1);
	
     parameters(VP0);
     parameters(VS0);
     parameters(det0);
     parameters(e0);
	 parameters(DEN0);
     
	 /*计算区域背景模型参数*/
	smooth(VP0,VP,NX,NZ,PML,10);
    smooth(VS0,VS,NX,NZ,PML,10);
    smooth(det0,det,NX,NZ,PML,10);
    smooth(e0,e,NX,NZ,PML,10); 
    smooth(DEN0,DEN,NX,NZ,PML,10); 
    
	/*for(J=PML;J<NX-PML;J++)
		for(K=PML;K<NZ-PML;K++)
		 {
	       DEN[J][K]=1500.0;
		 }*/
	PML_coefficient(dx,dz);
    parameters(VP);
    parameters(VS);
    parameters(DEN);
    parameters(det);
    parameters(e);

    tao=0.75; 	 
     
	for(J=0;J<NX;J++)
		 for(K=0;K<NZ;K++)
		  {
			filter[J][K]=0.0;
		 }
	for(J=PML;J<NX-PML;J++)
		for(K=PML;K<NZ-PML;K++)
		 {
	       if(K>=int(10*sin(2*pi*(J-25-PML)/100)+10+PML))
	       {
	        filter[J][K]=1.0;
	       }
	       else
	       {
	        filter[J][K]=0.0;
	       }
		 }
    for(J=0;J<=NX-1;J++)	
         for(K=0;K<=NZ-1;K++)
         {
           //det[J][K]=det0[J][K];
           //e[J][K]=e0[J][K];
           
           C33[J][K]=VP[J][K]*VP[J][K]*DEN[J][K];
           C44[J][K]=VS[J][K]*VS[J][K]*DEN[J][K];
           C11[J][K]=(2*e[J][K]+1)*C33[J][K];
           C13[J][K]=sqrt(((2*det[J][K]+1)*C33[J][K]-C44[J][K])*(C33[J][K]-C44[J][K]))-C44[J][K];
           c13[J][K]=sqrt(2*det[J][K]+1)*C33[J][K];
           //c13[J][K]=sqrt(((2*det[J][K]+1)*C33[J][K]-(fabs(e[J][K]-det[J][K])/tao)*C33[J][K])*(C33[J][K]-(fabs(e[J][K]-det[J][K])/tao)*C33[J][K]))-(fabs(e[J][K]-det[J][K])/tao)*C33[J][K];
           
           A=(C13[J][K]+2.0*C44[J][K])/C33[J][K];
           B=-2.0;
           M=C13[J][K]*C13[J][K]-C11[J][K]*C33[J][K];
           
           D11p[J][K]=(C13[J][K]*A*C33[J][K]-C11[J][K]*C33[J][K])/M;
           D12p[J][K]=(C13[J][K]*C11[J][K]-C11[J][K]*A*C33[J][K])/M;
           D13p[J][K]=0.0;
           
           D21p[J][K]=(C13[J][K]*C33[J][K]-C33[J][K]*A*C33[J][K])/M;
           D22p[J][K]=(C13[J][K]*A*C33[J][K]-C11[J][K]*C33[J][K])/M;
           D23p[J][K]=0.0;
           
           D11s[J][K]=(C13[J][K]*B*C44[J][K])/M;
           D12s[J][K]=(-C11[J][K]*B*C44[J][K])/M;
           D13s[J][K]=1.0;
           
           D21s[J][K]=(-C33[J][K]*B*C44[J][K])/M;
           D22s[J][K]=(C13[J][K]*B*C44[J][K])/M;
           D23s[J][K]=1.0;
         }     	       
              	 	 	 
    /*for(J=PML;J<NX-PML;J++)
		for(K=PML;K<NZ-PML;K++)
		 {
	       if(K>=int(30.0*sin(2*pi*(J-50)/200.0)+60.0))
	       {
	        filter[J][K]=1.0;
	       }
	       else
	       {
	        filter[J][K]=0.0;
	       }
		 }*/
	
	//parameters(filter);	
	
    /*for(I=0;I<=NT;I++)
	{
		W=pi*f0*DT*(I-100);
	    POW=(1.0-2*W*W)*exp(-W*W);
		wavelet[I]=POW; 
	} */
     DWORD start,end;
     start=GetTickCount();
	 
	 Sum1=0.0;
     Sum2=0.0;   
 
    for(J=0;J<NX;J++)
		for(K=0;K<NZ;K++)
		{
           I1[J][K]=0.0;
           I2[J][K]=0.0;
           I3[J][K]=0.0;
           I4[J][K]=0.0;
           I5[J][K]=0.0;
           I6[J][K]=0.0;
           I7[J][K]=0.0;
           I8[J][K]=0.0;
           I9[J][K]=0.0;
           I10[J][K]=0.0;
           
           Image1[J][K]=0.0;
           Image2[J][K]=0.0;
           Image3[J][K]=0.0;
           Image4[J][K]=0.0;
           Image5[J][K]=0.0;
           Image6[J][K]=0.0;
           Image7[J][K]=0.0;
           Image8[J][K]=0.0;
           Image9[J][K]=0.0;
           Image10[J][K]=0.0;
           
           W1[J][K]=0.0; 
           W2[J][K]=0.0;
           W3[J][K]=0.0; 
           W4[J][K]=0.0;
           W5[J][K]=0.0; 
           W6[J][K]=0.0;
           W7[J][K]=0.0; 
           W8[J][K]=0.0;
           W9[J][K]=0.0; 
           W10[J][K]=0.0;
		}
    for(MX=PML+MX1;MX<=NX-PML-1;MX=MX+MX2)
    {
	  MZ=tomo(MX);
      printf("MX=%d\n",MX);
      printf("MZ=%d\n",MZ);
      
	  for(J=0;J<NX;J++)
		for(K=0;K<NZ;K++)
		{
            sum1[J][K]=0.0;
            sum2[J][K]=0.0;
            sum3[J][K]=0.0;
            sum4[J][K]=0.0;
            sum5[J][K]=0.0;
            sum6[J][K]=0.0;
            sum7[J][K]=0.0;
            sum8[J][K]=0.0;
            sum9[J][K]=0.0;
            sum10[J][K]=0.0;
            
            SUM1[J][K]=0.0;
			SUM2[J][K]=0.0;
			SUM3[J][K]=0.0;
			SUM4[J][K]=0.0;
			SUM5[J][K]=0.0;
			SUM6[J][K]=0.0;
			SUM7[J][K]=0.0;
			SUM8[J][K]=0.0;
			SUM9[J][K]=0.0;
			SUM10[J][K]=0.0;
		} 
		/*变量初始化*/
      for(I=0;I<2;I++)
        for(J=0;J<NX;J++)
		 for(K=0;K<NZ;K++)
		  {
			P[I][J][K]=0.0;
            PX[I][J][K]=0.0;
            PZ[I][J][K]=0.0;
            Q[I][J][K]=0.0;
            QX[I][J][K]=0.0;
            QZ[I][J][K]=0.0;
            S[I][J][K]=0.0;
            SX[I][J][K]=0.0;
            SZ[I][J][K]=0.0;
            
            tH[I][J][K]=0.0;
			tHX[I][J][K]=0.0;
			tHZ[I][J][K]=0.0;
			tV[I][J][K]=0.0;
			tVX[I][J][K]=0.0;
			tVZ[I][J][K]=0.0;
            
            U[I][J][K]=0.0;
			UX[I][J][K]=0.0;
            UZ[I][J][K]=0.0;
            Up[I][J][K]=0.0;
			UpX[I][J][K]=0.0;
            UpZ[I][J][K]=0.0;
            Us[I][J][K]=0.0;
            UsX[I][J][K]=0.0;
            UsZ[I][J][K]=0.0;
            
            V[I][J][K]=0.0;
            VX[I][J][K]=0.0;
            VZ[I][J][K]=0.0;
            Vp[I][J][K]=0.0;
            VpX[I][J][K]=0.0;
            VpZ[I][J][K]=0.0;
            Vs[I][J][K]=0.0;
            VsX[I][J][K]=0.0;
            VsZ[I][J][K]=0.0;
		  }
      for(J=0;J<NX;J++) 
	    for(I=0;I<=NT;I++)
	    {
		 UsX1[J][I]=0.0;
         UsX2[J][I]=0.0;
         UsX3[J][I]=0.0;
         UsX4[J][I]=0.0;
         UsX5[J][I]=0.0;
         UsX6[J][I]=0.0;
         UsX7[J][I]=0.0;
         UsX8[J][I]=0.0;
         UsX9[J][I]=0.0;
         UsX10[J][I]=0.0;
         UsX11[J][I]=0.0;
         UsX12[J][I]=0.0;
        
	     UsZ1[J][I]=0.0;
         UsZ2[J][I]=0.0;
         UsZ3[J][I]=0.0;
         UsZ4[J][I]=0.0;
         UsZ5[J][I]=0.0;
         UsZ6[J][I]=0.0;
         UsZ7[J][I]=0.0;
         UsZ8[J][I]=0.0;
         UsZ9[J][I]=0.0;
         UsZ10[J][I]=0.0;
         UsZ11[J][I]=0.0;
         UsZ12[J][I]=0.0;
	    }
      for(K=0;K<NZ;K++) 
	    for(I=0;I<=NT;I++)
	    {
	     
         UsX13[K][I]=0.0;
         UsX14[K][I]=0.0;
         UsX15[K][I]=0.0;
         UsX16[K][I]=0.0;
         UsX17[K][I]=0.0;
         UsX18[K][I]=0.0;
         UsX19[K][I]=0.0;
         UsX20[K][I]=0.0;
         UsX21[K][I]=0.0;
         UsX22[K][I]=0.0;
         UsX23[K][I]=0.0;
         UsX24[K][I]=0.0;
         
         UsZ13[K][I]=0.0;
         UsZ14[K][I]=0.0;
         UsZ15[K][I]=0.0;
         UsZ16[K][I]=0.0;
         UsZ17[K][I]=0.0;
         UsZ18[K][I]=0.0;
         UsZ19[K][I]=0.0;
         UsZ20[K][I]=0.0;
         UsZ21[K][I]=0.0;
         UsZ22[K][I]=0.0;
         UsZ23[K][I]=0.0;
         UsZ24[K][I]=0.0;
	    }
  for(I=0;I<=NT;I++)
   {
	  //printf("%d\n",I);
	 
      //U[0][MX][MZ]=U[0][MX][MZ]+w1;
      //V[0][MX][MZ]=V[0][MX][MZ]+w1;
	 //foward(I,MX,MZ,A1,A2,A3,A4,A5,A6,C11,C13,C33,C44,c13,DEN,dx,dz,P,PX,PZ,Q,QX,QZ,S,SX,SZ,tH,tHX,tHZ,tV,tVX,tVZ,U,UX,UZ,Up,UpX,UpZ,Us,V,VX,VZ,Vp,VpX,VpZ,Vs);
	 foward1(I,MX,MZ,A1,A2,A3,A4,A5,A6,C11,C13,C33,C44,D11p,D12p,D13p,D21p,D22p,D23p,D11s,D12s,D13s,D21s,D22s,D23s,DEN,dx,dz,
	         P,PX,PZ,Q,QX,QZ,S,SX,SZ,U,UX,UZ,Up,UpX,UpZ,Us,UsX,UsZ,V,VX,VZ,Vp,VpX,VpZ,Vs,VsX,VsZ);
	 for(J=0;J<=NX-1;J++)
	  {
        UsX1[J][I]=U[1][J][PML-1];
        UsX2[J][I]=U[1][J][PML-2];
        UsX3[J][I]=U[1][J][PML-3];
        UsX4[J][I]=U[1][J][PML-4];
        UsX5[J][I]=U[1][J][PML-5];
        UsX6[J][I]=U[1][J][PML-6];
	    
        UsX7[J][I]=U[1][J][NZ-PML+4];
        UsX8[J][I]=U[1][J][NZ-PML+3];
        UsX9[J][I]=U[1][J][NZ-PML+2];
        UsX10[J][I]=U[1][J][NZ-PML+1];
        UsX11[J][I]=U[1][J][NZ-PML];
        UsX12[J][I]=U[1][J][NZ-PML-1];
        
	    UsZ1[J][I]=V[1][J][PML];
	    UsZ2[J][I]=V[1][J][PML-1];
        UsZ3[J][I]=V[1][J][PML-2];
        UsZ4[J][I]=V[1][J][PML-3];
        UsZ5[J][I]=V[1][J][PML-4];
        UsZ6[J][I]=V[1][J][PML-5];
	    
        UsZ7[J][I]=V[1][J][NZ-PML+5];
        UsZ8[J][I]=V[1][J][NZ-PML+4];
        UsZ9[J][I]=V[1][J][NZ-PML+3];
        UsZ10[J][I]=V[1][J][NZ-PML+2];
        UsZ11[J][I]=V[1][J][NZ-PML+1];
        UsZ12[J][I]=V[1][J][NZ-PML];
	  }
      for( K=0;K<=NZ-1;K++)    
	  {
        UsX13[K][I]=U[1][PML][K];
		UsX14[K][I]=U[1][PML-1][K];
		UsX15[K][I]=U[1][PML-2][K];
		UsX16[K][I]=U[1][PML-3][K];
		UsX17[K][I]=U[1][PML-4][K];
		UsX18[K][I]=U[1][PML-5][K];
		
		UsX19[K][I]=U[1][NX-PML+5][K];
		UsX20[K][I]=U[1][NX-PML+4][K];
		UsX21[K][I]=U[1][NX-PML+3][K];
		UsX22[K][I]=U[1][NX-PML+2][K];
		UsX23[K][I]=U[1][NX-PML+1][K];
		UsX24[K][I]=U[1][NX-PML][K];
		
	    UsZ13[K][I]=V[1][PML-1][K];
	    UsZ14[K][I]=V[1][PML-2][K];
		UsZ15[K][I]=V[1][PML-3][K];
		UsZ16[K][I]=V[1][PML-4][K];
		UsZ17[K][I]=V[1][PML-5][K];
		UsZ18[K][I]=V[1][PML-6][K];
		
		UsZ19[K][I]=V[1][NX-PML+4][K];
		UsZ20[K][I]=V[1][NX-PML+3][K];
		UsZ21[K][I]=V[1][NX-PML+2][K];
		UsZ22[K][I]=V[1][NX-PML+1][K];
		UsZ23[K][I]=V[1][NX-PML][K];
		UsZ24[K][I]=V[1][NX-PML-1][K];
	  }
/*.................将当前的（1时刻）波场值赋给0时刻作为下一时刻的0时刻值...............................*/	  
	for(J=0;J<=NX-1;J++)
		for(K=0;K<=NZ-1;K++)
		{
		U[0][J][K]=U[1][J][K];
        UX[0][J][K]=UX[1][J][K];
        UZ[0][J][K]=UZ[1][J][K];
        Up[0][J][K]=Up[1][J][K];
        UpX[0][J][K]=UpX[1][J][K];
        UpZ[0][J][K]=UpZ[1][J][K];
        Us[0][J][K]=Us[1][J][K];
        UsX[0][J][K]=UsX[1][J][K];
        UsZ[0][J][K]=UsZ[1][J][K];
        
		V[0][J][K]=V[1][J][K];
        VX[0][J][K]=VX[1][J][K];
        VZ[0][J][K]=VZ[1][J][K];
        Vp[0][J][K]=Vp[1][J][K];
        VpX[0][J][K]=VpX[1][J][K];
        VpZ[0][J][K]=VpZ[1][J][K];
        Vs[0][J][K]=Vs[1][J][K];
        VsX[0][J][K]=VsX[1][J][K];
        VsZ[0][J][K]=VsZ[1][J][K];
        
		P[0][J][K]=P[1][J][K];
        PX[0][J][K]=PX[1][J][K];
        PZ[0][J][K]=PZ[1][J][K];
        Q[0][J][K]=Q[1][J][K];
        QX[0][J][K]=QX[1][J][K];
        QZ[0][J][K]=QZ[1][J][K];
        S[0][J][K]=S[1][J][K];
        SX[0][J][K]=SX[1][J][K];
        SZ[0][J][K]=SZ[1][J][K];
        
        tH[0][J][K]=tH[1][J][K];
        tHX[0][J][K]=tHX[1][J][K];
        tHZ[0][J][K]=tHZ[1][J][K];
        tV[0][J][K]=tV[1][J][K];
		tVX[0][J][K]=tVX[1][J][K];
		tVZ[0][J][K]=tVZ[1][J][K];
		}	
    }
     for(J=0;J<NX;J++) 
	    for(I=0;I<=NT;I++)
	    {
		  UsX_obs[J][I]=0.0;
		  UsZ_obs[J][I]=0.0;
	    }

/*..................................读取数据...........................*/ 
       sprintf(filename6,"outputsX(起伏地表,可导,各向异性) X=%d.dat",MX);
       sprintf(filename7,"outputsZ(起伏地表,可导,各向异性) X=%d.dat",MX);
       f6[MX]=fopen(filename6,"rb");
       f7[MX]=fopen(filename7,"rb");
		 for(J=PML;J<NX-PML;J++)
           for(I=0;I<=NT;I++)
		   {
           fread(&LS6,sizeof(float),1,f6[MX]);
		   fread(&LS7,sizeof(float),1,f7[MX]);
	       UsX_obs[J][I]=LS6;
           UsZ_obs[J][I]=LS7;
		   }
         fclose(f6[MX]);
         fclose(f7[MX]);
	   /*sprintf(filename6,"outputsX(X=22) X=%d.dat",MX);
       sprintf(filename7,"outputsZ(X=22) X=%d.dat",MX);
       f6[MX]=fopen(filename6,"rb");
       f7[MX]=fopen(filename7,"rb");
		
		 for(J=PML;J<NZ-PML;J++)
           for(I=0;I<=NT;I++)
		   {
           fread(&LS6,sizeof(float),1,f6[MX]);
		   fread(&LS7,sizeof(float),1,f7[MX]);
	       USX[J][I]=LS6;
           USZ[J][I]=LS7;
		   }
         fclose(f6[MX]);
         fclose(f7[MX]);
       sprintf(filename6,"outputsX(X=171)1_down X=%d.dat",MX);
       sprintf(filename7,"outputsZ(X=171)1_down X=%d.dat",MX);
       f6[MX]=fopen(filename6,"rb");
       f7[MX]=fopen(filename7,"rb");
		
		 for(J=PML;J<NZ-PML;J++)
           for(I=0;I<=NT;I++)
		   {
           fread(&LS6,sizeof(float),1,f6[MX]);
		   fread(&LS7,sizeof(float),1,f7[MX]);
	       USX1[J][I]=LS6;
           USZ1[J][I]=LS7;
		   }
         fclose(f6[MX]);
         fclose(f7[MX]);*/

     //printf("%30.22f\n",UsX[80][600]);
     //printf("%30.22f\n",UsZ[80][600]); 
 	/*正演波场（震源）重构*/       
	 for(I=0;I<2;I++)
        for(J=0;J<NX;J++)
		 for(K=0;K<NZ;K++)
		  {
            P2[I][J][K]=0.0;
            Q2[I][J][K]=0.0;
            S2[I][J][K]=0.0;
            
            tH2[I][J][K]=0.0;
			tV2[I][J][K]=0.0;
			
            U2[I][J][K]=0.0;
            V2[I][J][K]=0.0;
            U2p[I][J][K]=0.0;
            U2s[I][J][K]=0.0;
            V2p[I][J][K]=0.0;
            V2s[I][J][K]=0.0;
            
            P4[I][J][K]=0.0;
            P4X[I][J][K]=0.0;
            P4Z[I][J][K]=0.0;
            
            Q4[I][J][K]=0.0;
            Q4X[I][J][K]=0.0;
            Q4Z[I][J][K]=0.0;
            
            S4[I][J][K]=0.0;
            S4X[I][J][K]=0.0;
            S4Z[I][J][K]=0.0;
            
            tH4[I][J][K]=0.0;
			tH4X[I][J][K]=0.0;
			tH4Z[I][J][K]=0.0;
			tV4[I][J][K]=0.0;
			tV4X[I][J][K]=0.0;
			tV4Z[I][J][K]=0.0;
			
            U4[I][J][K]=0.0;
			U4X[I][J][K]=0.0;
            U4Z[I][J][K]=0.0;
            
            U4p[I][J][K]=0.0;
			U4pX[I][J][K]=0.0;
            U4pZ[I][J][K]=0.0;
            U4s[I][J][K]=0.0;
            U4sX[I][J][K]=0.0;
            U4sZ[I][J][K]=0.0;
            
            V4[I][J][K]=0.0;
            V4X[I][J][K]=0.0;
            V4Z[I][J][K]=0.0;
            V4p[I][J][K]=0.0;
            V4pX[I][J][K]=0.0;
            V4pZ[I][J][K]=0.0;
            V4s[I][J][K]=0.0;
            V4sX[I][J][K]=0.0;
            V4sZ[I][J][K]=0.0;
		  }
	
	 for(J=0;J<NX;J++)
		for(K=0;K<NZ;K++)
		{
 	        P2[1][J][K]=P[0][J][K];
 	        Q2[1][J][K]=Q[0][J][K];
 	        S2[1][J][K]=S[0][J][K];
 	        tH2[1][J][K]=tH[0][J][K];
 	        tV2[1][J][K]=tV[0][J][K];
 	        
            U2[1][J][K]=U[0][J][K];
            V2[1][J][K]=V[0][J][K];
            U2p[1][J][K]=Up[0][J][K];
            V2p[1][J][K]=Vp[0][J][K];
            U2s[1][J][K]=Us[0][J][K];
            V2s[1][J][K]=Vs[0][J][K];
        }
 /*用于正演波场重构*/
/*用于波场快照的逆时外推初始条件*/
   for(I=NT;I>=NT1;I--)	
    {   
	   for(J=PML;J<=NX-PML-1;J++)
	  {
        U2[1][J][PML-1]=UsX1[J][I];
        U2[1][J][PML-2]=UsX2[J][I];
        U2[1][J][PML-3]=UsX3[J][I];
        U2[1][J][PML-4]=UsX4[J][I];
        U2[1][J][PML-5]=UsX5[J][I];
        U2[1][J][PML-6]=UsX6[J][I];
        
        U2[1][J][NZ-PML+4]=UsX7[J][I];
        U2[1][J][NZ-PML+3]=UsX8[J][I];
        U2[1][J][NZ-PML+2]=UsX9[J][I];
        U2[1][J][NZ-PML+1]=UsX10[J][I];
        U2[1][J][NZ-PML]=UsX11[J][I];
        U2[1][J][NZ-PML-1]=UsX12[J][I];
        
        
	    V2[1][J][PML]=UsZ1[J][I];
	    V2[1][J][PML-1]=UsZ2[J][I];
        V2[1][J][PML-2]=UsZ3[J][I];
        V2[1][J][PML-3]=UsZ4[J][I];
        V2[1][J][PML-4]=UsZ5[J][I];
        V2[1][J][PML-5]=UsZ6[J][I];
	    
        V2[1][J][NZ-PML+5]=UsZ7[J][I];
        V2[1][J][NZ-PML+4]=UsZ8[J][I];
        V2[1][J][NZ-PML+3]=UsZ9[J][I];
        V2[1][J][NZ-PML+2]=UsZ10[J][I];
        V2[1][J][NZ-PML+1]=UsZ11[J][I];
        V2[1][J][NZ-PML]=UsZ12[J][I];
	  }
       for(K=PML;K<=NZ-PML-1;K++)   
	  {
        U2[1][PML][K]=UsX13[K][I];
		U2[1][PML-1][K]=UsX14[K][I];
		U2[1][PML-2][K]=UsX15[K][I];
		U2[1][PML-3][K]=UsX16[K][I];
		U2[1][PML-4][K]=UsX17[K][I];
		U2[1][PML-5][K]=UsX18[K][I];
		
		U2[1][NX-PML+5][K]=UsX19[K][I];
		U2[1][NX-PML+4][K]=UsX20[K][I];
		U2[1][NX-PML+3][K]=UsX21[K][I];
		U2[1][NX-PML+2][K]=UsX22[K][I];
		U2[1][NX-PML+1][K]=UsX23[K][I];
		U2[1][NX-PML][K]=UsX24[K][I];
		
	    V2[1][PML-1][K]=UsZ13[K][I];
	    V2[1][PML-2][K]=UsZ14[K][I];
		V2[1][PML-3][K]=UsZ15[K][I];
		V2[1][PML-4][K]=UsZ16[K][I];
		V2[1][PML-5][K]=UsZ17[K][I];
		V2[1][PML-6][K]=UsZ18[K][I];
		
		V2[1][NX-PML+4][K]=UsZ19[K][I];
		V2[1][NX-PML+3][K]=UsZ20[K][I];
		V2[1][NX-PML+2][K]=UsZ21[K][I];
		V2[1][NX-PML+1][K]=UsZ22[K][I];
		V2[1][NX-PML][K]=UsZ23[K][I];
		V2[1][NX-PML-1][K]=UsZ24[K][I];
		
	  }
	//reconstruct(A1,A2,A3,A4,A5,A6,C11,C13,C33,C44,c13,DEN,P2,Q2,S2,tH2,tV2,U2,U2p,U2s,V2,V2p,V2s);
	reconstruct1(A1,A2,A3,A4,A5,A6,C11,C13,C33,C44,D11p,D12p,D13p,D21p,D22p,D23p,D11s,D12s,D13s,D21s,D22s,D23s,DEN,P2,Q2,S2,U2,U2p,U2s,V2,V2p,V2s);
    /*.................将当前的（1时刻）波场值赋给0时刻作为下一时刻的0时刻值...............................*/ 
	
	for(J=PML;J<=NX-PML-1;J++)
	  for(K=PML;K<=NZ-PML-1;K++)
		{
		 P2[1][J][K]=P2[0][J][K];
		 Q2[1][J][K]=Q2[0][J][K];
		 S2[1][J][K]=S2[0][J][K];
		 tH2[1][J][K]=tH2[0][J][K];
		 tV2[1][J][K]=tV2[0][J][K];
		 U2[1][J][K]=U2[0][J][K];
		 V2[1][J][K]=V2[0][J][K];
		 
		 U2p[1][J][K]=U2p[0][J][K];
		 V2p[1][J][K]=V2p[0][J][K];
		 U2s[1][J][K]=U2s[0][J][K];
		 V2s[1][J][K]=V2s[0][J][K];
		}	  	 
    /*产生Ur-------4波场值（代表接收点波场下行波）*/      
	   for(J=PML;J<NX-PML;J++)
      {
        
            MZ1=tomo(J);
	        U4[1][J][MZ1]=UsX_obs[J][I];
	        V4[1][J][MZ1]=UsZ_obs[J][I];
	  }
     //adjoint(A1,A2,A3,A4,A5,A6,C11,C13,C33,C44,c13,DEN,dx,dz,P4,P4X,P4Z,Q4,Q4X,Q4Z,S4,S4X,S4Z,tH4,tH4X,tH4Z,tV4,tV4X,tV4Z,U4,U4X,U4Z,U4p,U4pX,U4pZ,U4s,V4,V4X,V4Z,V4p,V4pX,V4pZ,V4s);
       adjoint1(A1,A2,A3,A4,A5,A6,C11,C13,C33,C44,D11p,D12p,D13p,D21p,D22p,D23p,D11s,D12s,D13s,D21s,D22s,D23s,DEN,dx,dz,
                P4,P4X,P4Z,Q4,Q4X,Q4Z,S4,S4X,S4Z,U4,U4X,U4Z,U4p,U4pX,U4pZ,U4s,U4sX,U4sZ,V4,V4X,V4Z,V4p,V4pX,V4pZ,V4s,V4sX,V4sZ);
/*.................将当前的（1时刻）波场值赋给0时刻作为下一时刻的0时刻值...............................*/
	 
	 for(J=0;J<=NX-1;J++)
	  for(K=0;K<=NZ-1;K++)
	  {
		U4[1][J][K]=U4[0][J][K];
        U4X[1][J][K]=U4X[0][J][K];
        U4Z[1][J][K]=U4Z[0][J][K];
		V4[1][J][K]=V4[0][J][K];
        V4X[1][J][K]=V4X[0][J][K];
        V4Z[1][J][K]=V4Z[0][J][K];
        
        U4p[1][J][K]=U4p[0][J][K];
        U4pX[1][J][K]=U4pX[0][J][K];
        U4pZ[1][J][K]=U4pZ[0][J][K];
        U4s[1][J][K]=U4s[0][J][K];
        U4sX[1][J][K]=U4sX[0][J][K];
        U4sZ[1][J][K]=U4sZ[0][J][K];
        
        V4p[1][J][K]=V4p[0][J][K];
        V4pX[1][J][K]=V4pX[0][J][K];
        V4pZ[1][J][K]=V4pZ[0][J][K];
        V4s[1][J][K]=V4s[0][J][K];
        V4sX[1][J][K]=V4sX[0][J][K];
        V4sZ[1][J][K]=V4sZ[0][J][K];
        
		P4[1][J][K]=P4[0][J][K];
        P4X[1][J][K]=P4X[0][J][K];
        P4Z[1][J][K]=P4Z[0][J][K];

		Q4[1][J][K]=Q4[0][J][K];
        Q4X[1][J][K]=Q4X[0][J][K];
        Q4Z[1][J][K]=Q4Z[0][J][K];
        
        S4[1][J][K]=S4[0][J][K];
        S4X[1][J][K]=S4X[0][J][K];
        S4Z[1][J][K]=S4Z[0][J][K];
        
        tH4[1][J][K]=tH4[0][J][K];
        tH4X[1][J][K]=tH4X[0][J][K];
        tH4Z[1][J][K]=tH4Z[0][J][K];
        
        tV4[1][J][K]=tV4[0][J][K];
        tV4X[1][J][K]=tV4X[0][J][K];
        tV4Z[1][J][K]=tV4Z[0][J][K];
	  }
	 /*for(J=0;J<=NX-1;J++)
	   for(K=0;K<=NZ-1;K++)
	   {
	    U2p1[J][K]=filter[J][K]*U2p[0][J][K];
	    V2p1[J][K]=filter[J][K]*V2p[0][J][K];
	    U4p1[J][K]=filter[J][K]*U4p[0][J][K];
	    V4p1[J][K]=filter[J][K]*V4p[0][J][K];
	    U4s1[J][K]=filter[J][K]*U4s[0][J][K];
	    V4s1[J][K]=filter[J][K]*V4s[0][J][K];
	   } */
	 for(J=PML;J<NX-PML;J++)
		for(K=PML;K<NZ-PML;K++)
		 {
	   
	     //sum7[J][K]+=U2p1[J][K]*U4p1[J][K]+V2p1[J][K]*V4p1[J][K];  
	     //sum8[J][K]+=U2p1[J][K]*U4s1[J][K]+V2p1[J][K]*V4s1[J][K];
	    
	     //SUM7[J][K]+=U2p1[J][K]*U2p1[J][K]+V2p1[J][K]*V2p1[J][K];
	     sum7[J][K]+=U2p[0][J][K]*U4p[0][J][K]+V2p[0][J][K]*V4p[0][J][K];  
	     sum8[J][K]+=U2p[0][J][K]*U4s[0][J][K]+V2p[0][J][K]*V4s[0][J][K];
	     SUM7[J][K]+=U2p[0][J][K]*U2p[0][J][K]+V2p[0][J][K]*V2p[0][J][K];
		} 
	  }
   printf("%30.22f\n",U2s[0][80][80]);
   printf("%30.22f\n",U4s[0][80][80]);
   //printf("%30.22f\n",UsX[80][600]);
   //printf("%30.22f\n",UsZ[80][600]); 
   //printf("%30.22f\n",sum3[100][100]);	
   printf("MX=%d\n",MX);
	for(J=PML;J<NX-PML;J++)
		for(K=PML;K<NZ-PML;K++)
		{
		   /*Image1[J][K]+=(sum1[J][K]/SUM1[J][K]);
           Image2[J][K]+=(sum2[J][K]/SUM2[J][K]);
           Image3[J][K]+=(sum3[J][K]/SUM3[J][K]);
           Image4[J][K]+=(sum4[J][K]/SUM4[J][K]);
           Image5[J][K]+=(sum5[J][K]/SUM5[J][K]);
           Image6[J][K]+=(sum6[J][K]/SUM6[J][K]);
           Image7[J][K]+=(sum7[J][K]/SUM7[J][K]);
           Image8[J][K]+=(sum8[J][K]/SUM8[J][K]);
           Image9[J][K]+=(sum9[J][K]/SUM9[J][K]);
           Image10[J][K]+=(sum10[J][K]/SUM10[J][K]);*/

           //Image1[J][K]+=sum1[J][K];
           //Image2[J][K]+=sum2[J][K];
           //Image3[J][K]+=sum3[J][K];
           //Image4[J][K]+=sum4[J][K];
           //Image5[J][K]+=sum5[J][K];
           //Image6[J][K]+=sum6[J][K];
           Image7[J][K]+=sum7[J][K];
           Image8[J][K]+=sum8[J][K];
           //Image9[J][K]+=sum9[J][K];
           //Image10[J][K]+=sum10[J][K];
           
           //W1[J][K]+=SUM1[J][K]; 
           //W2[J][K]+=SUM2[J][K];
           //W3[J][K]+=SUM3[J][K]; 
           //W4[J][K]+=SUM4[J][K];
           //W5[J][K]+=SUM5[J][K]; 
           //W6[J][K]+=SUM6[J][K];
           W7[J][K]+=SUM7[J][K]; 
           //W8[J][K]+=SUM8[J][K];
           //W9[J][K]+=SUM9[J][K]; 
           //W10[J][K]+=SUM10[J][K];
		}
	}
	for(J=PML;J<NX-PML;J++)
		for(K=PML;K<NZ-PML;K++)
		{
		   /*I1[J][K]=Image1[J][K];
           I2[J][K]=Image2[J][K];
           I3[J][K]=Image3[J][K];
           I4[J][K]=Image4[J][K];
           I5[J][K]=Image5[J][K];
           I6[J][K]=Image6[J][K];
           I7[J][K]=Image7[J][K];
           I8[J][K]=Image8[J][K];
           I9[J][K]=Image9[J][K];
           I10[J][K]=Image10[J][K];*/


           //I1[J][K]=Image1[J][K]/W1[J][K];
           //I2[J][K]=Image2[J][K]/W2[J][K];
           //I3[J][K]=Image3[J][K]/W3[J][K];
           //I4[J][K]=Image4[J][K]/W4[J][K];
           //I5[J][K]=Image5[J][K]/W5[J][K];
           //I6[J][K]=Image6[J][K]/W6[J][K];
           I7[J][K]=Image7[J][K]/W7[J][K];
           I8[J][K]=Image8[J][K]/W7[J][K];
           //I9[J][K]=Image9[J][K]/W9[J][K];
           //I10[J][K]=Image10[J][K]/W10[J][K];
		}
	    
	    f2=fopen("PP剖面(起伏地表,可导,分解方法2,滤波前).dat","wb");
	    f3=fopen("PS剖面(起伏地表,可导,分解方法2,滤波前).dat","wb");
	    
	   for(J=PML;J<NX-PML;J++)
	     for(K=PML;K<NZ-PML;K++)
		 {
		 LS2=float(I7[J][K]);
		 LS3=float(I8[J][K]);
		 fwrite(&LS2,sizeof(float),1,f2);
		 fwrite(&LS3,sizeof(float),1,f3);
		 }  
	    
	    fclose(f2);
	    fclose(f3);
	   
	
	laplace_filter(I7,Image7);
	laplace_filter(I8,Image8);
	
	    f2=fopen("PP剖面(起伏地表,可导,分解方法2,滤波后).dat","wb");
	    f3=fopen("PS剖面(起伏地表,可导,分解方法2,滤波后).dat","wb");
	
	   for(J=PML;J<NX-PML;J++)
	     for(K=PML;K<NZ-PML;K++)
		 {
		 LS2=float(Image7[J][K]);
		 LS3=float(Image8[J][K]);
		 fwrite(&LS2,sizeof(float),1,f2);
		 fwrite(&LS3,sizeof(float),1,f3);
		 }  
	   
	    fclose(f2);
	    fclose(f3);
	   
	 for(J=PML;J<NX-PML;J++)
	     for(K=PML;K<NZ-PML;K++)
	   {
	    Image7[J][K]=filter[J][K]*Image7[J][K];
	    Image8[J][K]=filter[J][K]*Image8[J][K];
	   } 
	   
	    f2=fopen("PP剖面(起伏地表,可导,分解方法2,最终成像结果).dat","wb");
	    f3=fopen("PS剖面(起伏地表,可导,分解方法2,最终成像结果).dat","wb");
	   
	    
	   for(J=PML;J<NX-PML;J++)
	     for(K=PML;K<NZ-PML;K++)
		 {
		
		 LS2=float(Image7[J][K]);
		 LS3=float(Image8[J][K]);
		
		 
		
		 fwrite(&LS2,sizeof(float),1,f2);
		 fwrite(&LS3,sizeof(float),1,f3);
		
		 }  
	   
	    fclose(f2);
	    fclose(f3);
	    
	       
    end=GetTickCount()-start;
	printf("time=%d N=%d\n",end,N);    
}
void laplace_filter(float **I,float **I1)
{
 int J,K;

for(J=0;J<NX;J++) 
      for(K=0;K<NZ;K++)
	  {
	      I1[J][K]=0;
	  }
/*for(J=PML;J<NX-PML;J++)
	for(K=PML;K<NZ-PML;K++)
	{
	    if((J>PML)&&(J<NX-1-PML)&&(K>PML)&&(K<NZ-1-PML))
		{ 
		I1[J][K]=(1.0/(DX*DX))*(I[J+1][K]-2*I[J][K]+I[J-1][K])+(1.0/(DZ*DZ))*(I[J][K+1]-2*I[J][K]+I[J][K-1]);
		}
	    if((J==PML)&&(K!=PML)&&(K!=NZ-1-PML))
		{  
		I1[J][K]=(1.0/DX)*(I[J+1][K]-I[J][K])+(1.0/(DZ*DZ))*(I[J][K+1]-2*I[J][K]+I[J][K-1]);
		}
        if((J==NX-1-PML)&&(K!=PML)&&(K!=NZ-1-PML))
		{ 
	    I1[J][K]=(1.0/DX)*(I[J][K]-I[J-1][K])+(1.0/(DZ*DZ))*(I[J][K+1]-2*I[J][K]+I[J][K-1]);
		}
        if((K==PML)&&(J!=PML)&&(J!=NX-1-PML))
	    { 
	    I1[J][K]=(1.0/(DX*DX))*(I[J+1][K]-2*I[J][K]+I[J-1][K])+(1.0/DZ)*(I[J][K+1]-I[J][K]);
		}
        if((K==NZ-1-PML)&&(J!=PML)&&(J!=NX-1-PML))
        { 
	    I1[J][K]=(1.0/(DX*DX))*(I[J+1][K]-2*I[J][K]+I[J-1][K])+(1.0/DZ)*(I[J][K]-I[J][K-1]);
		}
        if((J==PML)&&(K==PML))
	    {   
		 I1[J][K]=(1.0/DX)*(I[J+1][K]-I[J][K])+(1.0/DZ)*(I[J][K+1]-I[J][K]);
		}
        if((J==PML)&&(K==NZ-1-PML))
		{   
		 I1[J][K]=(1.0/DX)*(I[J+1][K]-I[J][K])+(1.0/DZ)*(I[J][K]-I[J][K-1]);
		}
        if((K==PML)&&(J==NX-1-PML))
	    { 
	     I1[J][K]=(1.0/DX)*(I[J][K]-I[J-1][K])+(1.0/DZ)*(I[J][K+1]-I[J][K]);
		}
        if((J==NX-1-PML)&&(K==NZ-1-PML))
		{ 
	     I1[J][K]=(1.0/DX)*(I[J][K]-I[J-1][K])+(1.0/DZ)*(I[J][K]-I[J][K-1]);
		}
	}*/
for(J=PML;J<NX-PML;J++)
	for(K=PML;K<NZ-PML;K++)
	{
	 I1[J][K]=(I[J+1][K]-2*I[J][K]+I[J-1][K])+(I[J][K+1]-2*I[J][K]+I[J][K-1]);
	}
}
/*计算衰减系数*/
void PML_coefficient(float **dx,float **dz)
{
    int J,K;
    float x,z,plx,plz,a;
    a=1.79;
    
    plx=PML*DX;      
    plz=PML*DZ;
    for(J=0;J<=NX-1;J++)	
      for(K=0;K<=NZ-1;K++)
      {
	    if((J>=0)&&(J<=PML-1)&&(K>=0)&&(K<=PML-1))	 
	   {
		 x=(PML-J)*DX;z=(PML-K)*DZ;
	     dx[J][K]=2*pi*f0*a*(float(x)/plx)*(float(x)/plx);
	     dz[J][K]=2*pi*f0*a*(float(z)/plz)*(float(z)/plz);
	   }
       else if((J>=0)&&(J<=PML-1)&&(K>=NZ-PML)&&(K<=NZ-1))
	   {
	     x=(PML-J)*DX;z=(K-(NZ-1-PML))*DZ;
         dx[J][K]=2*pi*f0*a*(float(x)/plx)*(float(x)/plx);
	     dz[J][K]=2*pi*f0*a*(float(z)/plz)*(float(z)/plz);
	   }
       else if ((J>=NX-PML)&&(J<=NX-1)&&(K>=0)&&(K<=PML-1))
	    {    
		 x=(J-(NX-1-PML))*DX;z=(PML-K)*DZ;
         dx[J][K]=2*pi*f0*a*(float(x)/plx)*(float(x)/plx);
	     dz[J][K]=2*pi*f0*a*(float(z)/plz)*(float(z)/plz);
	    }
        else if((J>=NX-PML)&&(J<=NX-1)&&(K>=NZ-PML)&&(K<=NZ-1))
	   {
		 x=(J-(NX-1-PML))*DX;z=(K-(NZ-1-PML))*DZ;
         dx[J][K]=2*pi*f0*a*(float(x)/plx)*(float(x)/plx);
	     dz[J][K]=2*pi*f0*a*(float(z)/plz)*(float(z)/plz);
	   }
       else if((J>=0)&&(J<=PML-1)&&(K>=PML)&&(K<NZ-PML))
	   {
		 x=(PML-J)*DX;z=0;
         dx[J][K]=2*pi*f0*a*(float(x)/plx)*(float(x)/plx);
	     dz[J][K]=0;
	   }
       else if((J>=NX-PML)&&(J<=NX-1)&&(K>=PML)&&(K<NZ-PML))
	   {      
	     x=(J-(NX-1-PML))*DX;z=0;
         dx[J][K]=2*pi*f0*a*(float(x)/plx)*(float(x)/plx);
	     dz[J][K]=0;
	   }
       else if((J>=PML)&&(J<=NX-1-PML)&&(K>=0)&&(K<=PML-1))
	   {    
	     x=0;z=(PML-K)*DZ;
	     dx[J][K]=0;
	     dz[J][K]=2*pi*f0*a*(float(z)/plz)*(float(z)/plz);
	   }
       else if((J>=PML)&&(J<=NX-1-PML)&&(K>=NZ-PML)&&(K<=NZ-1))
	   {     
		 x=0;z=(K-(NZ-1-PML))*DZ;
         dx[J][K]=0;
	     dz[J][K]=2*pi*f0*a*(float(z)/plz)*(float(z)/plz);
	   }
       else if((J>=PML)&&(J<=NX-1-PML)&&(K>=PML)&&(K<=NZ-1-PML))
	   { 
		  dx[J][K]=0;
		  dz[J][K]=0;
	   }
    }
}
 /*求边界处模型参数*/
void  parameters(float **VP)
{
  int I,J,K;
   for(J=0;J<=NX-1;J++)	
     for(K=0;K<=NZ-1;K++)
    {
	   if((J>=0)&&(J<=PML-1)&&(K>=0)&&(K<=PML-1))	 
	   {
		VP[J][K]=VP[PML][PML];
	   }
       if((J>=0)&&(J<=PML-1)&&(K>=NZ-PML)&&(K<=NZ-1))
	   {
	    VP[J][K]=VP[PML][NZ-PML-1];
	   }
       if ((J>=NX-PML)&&(J<=NX-1)&&(K>=0)&&(K<=PML-1))
	   {    
		VP[J][K]=VP[NX-PML-1][PML];
	   }
       if((J>=NX-PML)&&(J<=NX-1)&&(K>=NZ-PML)&&(K<=NZ-1))
	   {
	    VP[J][K]=VP[NX-PML-1][NZ-PML-1];
	   }
       if((J>=0)&&(J<=PML-1)&&(K>=PML)&&(K<NZ-PML))
	   {
		VP[J][K]=VP[PML][K];
	   }
       if((J>=NX-PML)&&(J<=NX-1)&&(K>=PML)&&(K<NZ-PML))
	   {   
		VP[J][K]=VP[NX-PML-1][K];
	   }
       if((J>=PML)&&(J<=NX-1-PML)&&(K>=0)&&(K<=PML-1))
	   {    
	    VP[J][K]=VP[J][PML];
	   }
       if((J>=PML)&&(J<=NX-1-PML)&&(K>=NZ-PML)&&(K<=NZ-1))
	   {     
		VP[J][K]=VP[J][NZ-1-PML];
	   }
       if((J>=PML)&&(J<=NX-1-PML)&&(K>=PML)&&(K<=NZ-1-PML))
	   { 
		VP[J][K]=VP[J][K];
	   }
    }
 }
void foward(int I,int MX,int MZ,float A1,float A2,float A3,float A4,float A5,float A6,float **C11,float **C13,float **C33,float **C44,float **c13,float **DEN,float **dx, float **dz,float ***P,float ***PX,float ***PZ,float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,float ***tH,float ***tHX,float ***tHZ,float ***tV,float ***tVX,float ***tVZ,float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs)
 {
        int J,K;
        float   W,POW,DX_P,DZ_P,DX_Q,DZ_Q,DX_S,DZ_S,DX_U,DZ_U,DX_V,DZ_V,DX_tH,DZ_tV; 
         
         W=pi*f0*DT*(I-100);
	     POW=(1.0-2*W*W)*exp(-W*W);
         P[0][MX][MZ]=P[0][MX][MZ]+POW;
         Q[0][MX][MZ]=Q[0][MX][MZ]+POW;   
           DX_P=0.0;
	       DZ_P=0.0;
	       DX_Q=0.0;
	       DZ_Q=0.0;
		   DX_S=0.0;
		   DZ_S=0.0;
		   DX_U=0.0;
	       DZ_U=0.0;
	       DX_V=0.0;
	       DZ_V=0.0;
	       DX_tH=0.0;
	       DZ_tV=0.0; 
/* .......................计算速度分量..............................*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
		{
	       DX_P=A1*(P[0][J][K]-P[0][J-1][K])+A2*(P[0][J+1][K]-P[0][J-2][K])+A3*(P[0][J+2][K]-P[0][J-3][K])+A4*(P[0][J+3][K]-P[0][J-4][K])+A5*(P[0][J+4][K]-P[0][J-5][K])+A6*(P[0][J+5][K]-P[0][J-6][K]);
           DZ_P=A1*(P[0][J][K]-P[0][J][K-1])+A2*(P[0][J][K+1]-P[0][J][K-2])+A3*(P[0][J][K+2]-P[0][J][K-3])+A4*(P[0][J][K+3]-P[0][J][K-4])+A5*(P[0][J][K+4]-P[0][J][K-5])+A6*(P[0][J][K+5]-P[0][J][K-6]);
           DX_Q=A1*(Q[0][J][K]-Q[0][J-1][K])+A2*(Q[0][J+1][K]-Q[0][J-2][K])+A3*(Q[0][J+2][K]-Q[0][J-3][K])+A4*(Q[0][J+3][K]-Q[0][J-4][K])+A5*(Q[0][J+4][K]-Q[0][J-5][K])+A6*(Q[0][J+5][K]-Q[0][J-6][K]);
           DZ_Q=A1*(Q[0][J][K]-Q[0][J][K-1])+A2*(Q[0][J][K+1]-Q[0][J][K-2])+A3*(Q[0][J][K+2]-Q[0][J][K-3])+A4*(Q[0][J][K+3]-Q[0][J][K-4])+A5*(Q[0][J][K+4]-Q[0][J][K-5])+A6*(Q[0][J][K+5]-Q[0][J][K-6]);
           
           DX_tH=A1*(tH[0][J][K]-tH[0][J-1][K])+A2*(tH[0][J+1][K]-tH[0][J-2][K])+A3*(tH[0][J+2][K]-tH[0][J-3][K])+A4*(tH[0][J+3][K]-tH[0][J-4][K])+A5*(tH[0][J+4][K]-tH[0][J-5][K])+A6*(tH[0][J+5][K]-tH[0][J-6][K]);
           DZ_tV=A1*(tV[0][J][K]-tV[0][J][K-1])+A2*(tV[0][J][K+1]-tV[0][J][K-2])+A3*(tV[0][J][K+2]-tV[0][J][K-3])+A4*(tV[0][J][K+3]-tV[0][J][K-4])+A5*(tV[0][J][K+4]-tV[0][J][K-5])+A6*(tV[0][J][K+5]-tV[0][J][K-6]);
		   
		   DX_S=A1*(S[0][J+1][K]-S[0][J][K])+A2*(S[0][J+2][K]-S[0][J-1][K])+A3*(S[0][J+3][K]-S[0][J-2][K])+A4*(S[0][J+4][K]-S[0][J-3][K])+A5*(S[0][J+5][K]-S[0][J-4][K])+A6*(S[0][J+6][K]-S[0][J-5][K]);
           DZ_S=A1*(S[0][J][K+1]-S[0][J][K])+A2*(S[0][J][K+2]-S[0][J][K-1])+A3*(S[0][J][K+3]-S[0][J][K-2])+A4*(S[0][J][K+4]-S[0][J][K-3])+A5*(S[0][J][K+5]-S[0][J][K-4])+A6*(S[0][J][K+6]-S[0][J][K-5]);
	       
	       U[1][J][K]=U[0][J][K]+(1.0/DEN[J][K])*(DT/DX)*DX_P+(1.0/DEN[J][K])*(DT/DZ)*DZ_S;
           V[1][J][K]=V[0][J][K]+(1.0/DEN[J][K])*(DT/DX)*DX_S+(1.0/DEN[J][K])*(DT/DZ)*DZ_Q;
           
           Up[1][J][K]=Up[0][J][K]+(1.0/DEN[J][K])*(DT/DX)*DX_tH;
           Vp[1][J][K]=Vp[0][J][K]+(1.0/DEN[J][K])*(DT/DZ)*DZ_tV;
	    }
/*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
          
	    }

  
/*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }    
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }    
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }   
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {     
	      UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	   for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }
 /*.....................四周降阶处理............................................*/
    if((J=0)&&(K=0))
		{
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
        }
    if((J=0)&&(K=NZ-1))
		{
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
		}
    if((J=NX-1)&&(K=NZ-1))
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K]; 
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
		}
    if((J=NX-1)&&(K=0))
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
		}

    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
          
		}
        if(K=NZ-1)
	    { 
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
		}
        if(J=NX-1)
	    { 
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+(tH[0][J][K]-tH[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];

          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+(tV[0][J][K]-tV[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K];
	    }
	}
    for(J=0;J<=NX-1;J++)
		for(K=0;K<=NZ-1;K++)
	    {
	     Us[1][J][K]=U[1][J][K]-Up[1][J][K];
	     Vs[1][J][K]=V[1][J][K]-Vp[1][J][K];
	    }
    
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
	    {
            DX_U=A1*(U[1][J+1][K]-U[1][J][K])+A2*(U[1][J+2][K]-U[1][J-1][K])+A3*(U[1][J+3][K]-U[1][J-2][K])+A4*(U[1][J+4][K]-U[1][J-3][K])+A5*(U[1][J+5][K]-U[1][J-4][K])+A6*(U[1][J+6][K]-U[1][J-5][K]);
            DZ_U=A1*(U[1][J][K]-U[1][J][K-1])+A2*(U[1][J][K+1]-U[1][J][K-2])+A3*(U[1][J][K+2]-U[1][J][K-3])+A4*(U[1][J][K+3]-U[1][J][K-4])+A5*(U[1][J][K+4]-U[1][J][K-5])+A6*(U[1][J][K+5]-U[1][J][K-6]);
            DX_V=A1*(V[1][J][K]-V[1][J-1][K])+A2*(V[1][J+1][K]-V[1][J-2][K])+A3*(V[1][J+2][K]-V[1][J-3][K])+A4*(V[1][J+3][K]-V[1][J-4][K])+A5*(V[1][J+4][K]-V[1][J-5][K])+A6*(V[1][J+5][K]-V[1][J-6][K]);
            DZ_V=A1*(V[1][J][K+1]-V[1][J][K])+A2*(V[1][J][K+2]-V[1][J][K-1])+A3*(V[1][J][K+3]-V[1][J][K-2])+A4*(V[1][J][K+4]-V[1][J][K-3])+A5*(V[1][J][K+5]-V[1][J][K-4])+A6*(V[1][J][K+6]-V[1][J][K-5]);
			
			P[1][J][K]=P[0][J][K]+(DT/DX)*C11[J][K]*DX_U+(DT/DZ)*C13[J][K]*DZ_V;			
	        
	        Q[1][J][K]=Q[0][J][K]+(DT/DX)*C13[J][K]*DX_U+(DT/DZ)*C33[J][K]*DZ_V;		
	        
	        S[1][J][K]=S[0][J][K]+(DT/DX)*C44[J][K]*DX_V+(DT/DZ)*C44[J][K]*DZ_U;
	        
	        tH[1][J][K]=tH[0][J][K]+(DT/DX)*C11[J][K]*DX_U+(DT/DZ)*c13[J][K]*DZ_V;
	        
	        tV[1][J][K]=tV[0][J][K]+(DT/DX)*c13[J][K]*DX_U+(DT/DZ)*C33[J][K]*DZ_V;
	    }
          
  
/*计算边界处各应力分量值*/
 /*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {  
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	    }
  
    /*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {

           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	    }   
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	    }   
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {	   
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    { 
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	    }  
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];

	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];

	    }
 /*.....................边界条件四周降阶处理............................................*/

    if((J=0)&&(K=0))
	 {	
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];                                      
       }
    if((J=0)&&(K=NZ-1))
      {
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
      }
    if((J=NX-1)&&(K=NZ-1))
      {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	     
	       tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
      }

    if((J=NX-1)&&(K=0))
      {
		   PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
      
      }
   
         
    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
		   PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
          
		}
       if(K=NZ-1)
	    { 
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
		         
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
		}
       if(J=NX-1)
	   { 
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
           
           tHX[1][J][K]=((1-0.5*DT*dx[J][K])*tHX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   tHZ[1][J][K]=((1-0.5*DT*dz[J][K])*tHZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[1][J][K]=tHX[1][J][K]+tHZ[1][J][K];

           tVX[1][J][K]=((1-0.5*DT*dx[J][K])*tVX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   tVZ[1][J][K]=((1-0.5*DT*dz[J][K])*tVZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[1][J][K]=tVX[1][J][K]+tVZ[1][J][K];
	   }
	}
   
}
void foward1(int I,int MX,int MZ,float A1,float A2,float A3,float A4,float A5,float A6,
             float **C11,float **C13,float **C33,float **C44,float **D11p,float **D12p,float **D13p,float **D21p,float **D22p,float **D23p,
             float **D11s,float **D12s,float **D13s,float **D21s,float **D22s,float **D23s,float **DEN,float **dx, float **dz,
             float ***P,float ***PX,float ***PZ,float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,
             float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***UsX,float ***UsZ,
             float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs,float ***VsX,float ***VsZ)
 {
        int J,K;
        float   W,POW,DX_P,DZ_P,DX_Q,DZ_Q,DX_S,DZ_S,DX_U,DZ_U,DX_V,DZ_V,DX_tH,DZ_tV; 
         
         W=pi*f0*DT*(I-100);
	     POW=(1.0-2*W*W)*exp(-W*W);
         P[0][MX][MZ]=P[0][MX][MZ]+POW;
         Q[0][MX][MZ]=Q[0][MX][MZ]+POW;   
           DX_P=0.0;
	       DZ_P=0.0;
	       DX_Q=0.0;
	       DZ_Q=0.0;
		   DX_S=0.0;
		   DZ_S=0.0;
		   DX_U=0.0;
	       DZ_U=0.0;
	       DX_V=0.0;
	       DZ_V=0.0;
	       DX_tH=0.0;
	       DZ_tV=0.0; 
/* .......................计算速度分量..............................*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
		{
	       DX_P=A1*(P[0][J][K]-P[0][J-1][K])+A2*(P[0][J+1][K]-P[0][J-2][K])+A3*(P[0][J+2][K]-P[0][J-3][K])+A4*(P[0][J+3][K]-P[0][J-4][K])+A5*(P[0][J+4][K]-P[0][J-5][K])+A6*(P[0][J+5][K]-P[0][J-6][K]);
           DZ_P=A1*(P[0][J][K]-P[0][J][K-1])+A2*(P[0][J][K+1]-P[0][J][K-2])+A3*(P[0][J][K+2]-P[0][J][K-3])+A4*(P[0][J][K+3]-P[0][J][K-4])+A5*(P[0][J][K+4]-P[0][J][K-5])+A6*(P[0][J][K+5]-P[0][J][K-6]);
           
           DX_Q=A1*(Q[0][J][K]-Q[0][J-1][K])+A2*(Q[0][J+1][K]-Q[0][J-2][K])+A3*(Q[0][J+2][K]-Q[0][J-3][K])+A4*(Q[0][J+3][K]-Q[0][J-4][K])+A5*(Q[0][J+4][K]-Q[0][J-5][K])+A6*(Q[0][J+5][K]-Q[0][J-6][K]);
           DZ_Q=A1*(Q[0][J][K]-Q[0][J][K-1])+A2*(Q[0][J][K+1]-Q[0][J][K-2])+A3*(Q[0][J][K+2]-Q[0][J][K-3])+A4*(Q[0][J][K+3]-Q[0][J][K-4])+A5*(Q[0][J][K+4]-Q[0][J][K-5])+A6*(Q[0][J][K+5]-Q[0][J][K-6]);
           
		   DX_S=A1*(S[0][J+1][K]-S[0][J][K])+A2*(S[0][J+2][K]-S[0][J-1][K])+A3*(S[0][J+3][K]-S[0][J-2][K])+A4*(S[0][J+4][K]-S[0][J-3][K])+A5*(S[0][J+5][K]-S[0][J-4][K])+A6*(S[0][J+6][K]-S[0][J-5][K]);
           DZ_S=A1*(S[0][J][K+1]-S[0][J][K])+A2*(S[0][J][K+2]-S[0][J][K-1])+A3*(S[0][J][K+3]-S[0][J][K-2])+A4*(S[0][J][K+4]-S[0][J][K-3])+A5*(S[0][J][K+5]-S[0][J][K-4])+A6*(S[0][J][K+6]-S[0][J][K-5]);
	       
	       U[1][J][K]=U[0][J][K]+(1.0/DEN[J][K])*(DT/DX)*DX_P+(1.0/DEN[J][K])*(DT/DZ)*DZ_S;
           V[1][J][K]=V[0][J][K]+(1.0/DEN[J][K])*(DT/DX)*DX_S+(1.0/DEN[J][K])*(DT/DZ)*DZ_Q;
           
           Up[1][J][K]=Up[0][J][K]+(D11p[J][K]/DEN[J][K])*(DT/DX)*DX_P+(D12p[J][K]/DEN[J][K])*(DT/DX)*DX_Q+(D13p[J][K]/DEN[J][K])*(DT/DZ)*DZ_S;
           Vp[1][J][K]=Vp[0][J][K]+(D21p[J][K]/DEN[J][K])*(DT/DZ)*DZ_P+(D22p[J][K]/DEN[J][K])*(DT/DZ)*DZ_Q+(D23p[J][K]/DEN[J][K])*(DT/DX)*DX_S;
           
           Us[1][J][K]=Us[0][J][K]+(D11s[J][K]/DEN[J][K])*(DT/DX)*DX_P+(D12s[J][K]/DEN[J][K])*(DT/DX)*DX_Q+(D13s[J][K]/DEN[J][K])*(DT/DZ)*DZ_S;
           Vs[1][J][K]=Vs[0][J][K]+(D21s[J][K]/DEN[J][K])*(DT/DZ)*DZ_P+(D22s[J][K]/DEN[J][K])*(DT/DZ)*DZ_Q+(D23s[J][K]/DEN[J][K])*(DT/DX)*DX_S;
	    }
/*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }

  
/*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }    
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }    
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }   
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {     
	      UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	   for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }
 /*.....................四周降阶处理............................................*/
    if((J=0)&&(K=0))
		{
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
        }
    if((J=0)&&(K=NZ-1))
		{
          UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
		}
    if((J=NX-1)&&(K=NZ-1))
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K]; 
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
		}
    if((J=NX-1)&&(K=0))
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
		}

    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
          
		}
        if(K=NZ-1)
	    { 
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K])/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K]+(S[0][J+1][K]-S[0][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
		}
        if(J=NX-1)
	    { 
		  UX[1][J][K]=((1-0.5*DT*dx[J][K])*UX[0][J][K]+(P[0][J][K]-P[0][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[1][J][K]=((1-0.5*DT*dz[J][K])*UZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[1][J][K]=UX[1][J][K]+UZ[1][J][K];

          VX[1][J][K]=((1-0.5*DT*dx[J][K])*VX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[1][J][K]=((1-0.5*DT*dz[J][K])*VZ[0][J][K]+(Q[0][J][K]-Q[0][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[1][J][K]=VX[1][J][K]+VZ[1][J][K];
          
          UpX[1][J][K]=((1-0.5*DT*dx[J][K])*UpX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[1][J][K]=((1-0.5*DT*dz[J][K])*UpZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[1][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[1][J][K]=((1-0.5*DT*dx[J][K])*VpX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[1][J][K]=((1-0.5*DT*dz[J][K])*VpZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[1][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[1][J][K]=((1-0.5*DT*dx[J][K])*UsX[0][J][K]+((P[0][J][K]-P[0][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[1][J][K]=((1-0.5*DT*dz[J][K])*UsZ[0][J][K]+(S[0][J][K+1]-S[0][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[1][J][K]=UsX[1][J][K]+UsZ[1][J][K];
          
          VsX[1][J][K]=((1-0.5*DT*dx[J][K])*VsX[0][J][K])/(1+0.5*DT*dx[J][K]);
		  VsZ[1][J][K]=((1-0.5*DT*dz[J][K])*VsZ[0][J][K]+((P[0][J][K]-P[0][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[0][J][K]-Q[0][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[1][J][K]=VsX[1][J][K]+VsZ[1][J][K]; 
	    }
	}
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
	    {
            DX_U=A1*(U[1][J+1][K]-U[1][J][K])+A2*(U[1][J+2][K]-U[1][J-1][K])+A3*(U[1][J+3][K]-U[1][J-2][K])+A4*(U[1][J+4][K]-U[1][J-3][K])+A5*(U[1][J+5][K]-U[1][J-4][K])+A6*(U[1][J+6][K]-U[1][J-5][K]);
            DZ_U=A1*(U[1][J][K]-U[1][J][K-1])+A2*(U[1][J][K+1]-U[1][J][K-2])+A3*(U[1][J][K+2]-U[1][J][K-3])+A4*(U[1][J][K+3]-U[1][J][K-4])+A5*(U[1][J][K+4]-U[1][J][K-5])+A6*(U[1][J][K+5]-U[1][J][K-6]);
            
            DX_V=A1*(V[1][J][K]-V[1][J-1][K])+A2*(V[1][J+1][K]-V[1][J-2][K])+A3*(V[1][J+2][K]-V[1][J-3][K])+A4*(V[1][J+3][K]-V[1][J-4][K])+A5*(V[1][J+4][K]-V[1][J-5][K])+A6*(V[1][J+5][K]-V[1][J-6][K]);
            DZ_V=A1*(V[1][J][K+1]-V[1][J][K])+A2*(V[1][J][K+2]-V[1][J][K-1])+A3*(V[1][J][K+3]-V[1][J][K-2])+A4*(V[1][J][K+4]-V[1][J][K-3])+A5*(V[1][J][K+5]-V[1][J][K-4])+A6*(V[1][J][K+6]-V[1][J][K-5]);
			
			P[1][J][K]=P[0][J][K]+(DT/DX)*C11[J][K]*DX_U+(DT/DZ)*C13[J][K]*DZ_V;			
	        Q[1][J][K]=Q[0][J][K]+(DT/DX)*C13[J][K]*DX_U+(DT/DZ)*C33[J][K]*DZ_V;		
	        S[1][J][K]=S[0][J][K]+(DT/DX)*C44[J][K]*DX_V+(DT/DZ)*C44[J][K]*DZ_U;
	    }
          
  
/*计算边界处各应力分量值*/
 /*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {  
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }
  
    /*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {

           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }   
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }   
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {	   
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    { 
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }  
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }
 /*.....................边界条件四周降阶处理............................................*/

    if((J=0)&&(K=0))
	 {	
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];                             
       }
    if((J=0)&&(K=NZ-1))
      {
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
      }
    if((J=NX-1)&&(K=NZ-1))
      {
           PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
      }

    if((J=NX-1)&&(K=0))
      {
		   PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
      }
   
         
    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
		   PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
		}
       if(K=NZ-1)
	    { 
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K])/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
		         
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K]+(U[1][J+1][K]-U[1][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
		}
       if(J=NX-1)
	   { 
	       PX[1][J][K]=((1-0.5*DT*dx[J][K])*PX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[1][J][K]=((1-0.5*DT*dz[J][K])*PZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[1][J][K]=PX[1][J][K]+PZ[1][J][K];

           QX[1][J][K]=((1-0.5*DT*dx[J][K])*QX[0][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[1][J][K]=((1-0.5*DT*dz[J][K])*QZ[0][J][K]+(V[1][J][K+1]-V[1][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[1][J][K]=QX[1][J][K]+QZ[1][J][K];

           SX[1][J][K]=((1-0.5*DT*dx[J][K])*SX[0][J][K]+(V[1][J][K]-V[1][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[1][J][K]=((1-0.5*DT*dz[J][K])*SZ[0][J][K]+(U[1][J][K]-U[1][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[1][J][K]=SX[1][J][K]+SZ[1][J][K];
	   }
	}
   
}
void reconstruct(float A1,float A2,float A3,float A4,float A5,float A6,float **C11,float **C13,float **C33,float **C44,float **c13,float **DEN,float ***P,float ***Q,float ***S,float ***tH,float ***tV,float ***U,float ***Up,float ***Us,float ***V,float ***Vp,float ***Vs)
{
    int J,K;
    float   DX_P,DZ_P,DX_Q,DZ_Q,DX_S,DZ_S,DX_U,DZ_U,DX_V,DZ_V,DX_tH,DZ_tV; 
           
           DX_P=0.0;
	       DZ_P=0.0;
	       DX_Q=0.0;
	       DZ_Q=0.0;
		   DX_S=0.0;
		   DZ_S=0.0;
		   DX_U=0.0;
	       DZ_U=0.0;
	       DX_V=0.0;
	       DZ_V=0.0;
	       DX_tH=0.0;
	       DZ_tV=0.0; 
	     
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
	    {
            DX_U=A1*(U[1][J+1][K]-U[1][J][K])+A2*(U[1][J+2][K]-U[1][J-1][K])+A3*(U[1][J+3][K]-U[1][J-2][K])+A4*(U[1][J+4][K]-U[1][J-3][K])+A5*(U[1][J+5][K]-U[1][J-4][K])+A6*(U[1][J+6][K]-U[1][J-5][K]);
            DZ_U=A1*(U[1][J][K]-U[1][J][K-1])+A2*(U[1][J][K+1]-U[1][J][K-2])+A3*(U[1][J][K+2]-U[1][J][K-3])+A4*(U[1][J][K+3]-U[1][J][K-4])+A5*(U[1][J][K+4]-U[1][J][K-5])+A6*(U[1][J][K+5]-U[1][J][K-6]);
            DX_V=A1*(V[1][J][K]-V[1][J-1][K])+A2*(V[1][J+1][K]-V[1][J-2][K])+A3*(V[1][J+2][K]-V[1][J-3][K])+A4*(V[1][J+3][K]-V[1][J-4][K])+A5*(V[1][J+4][K]-V[1][J-5][K])+A6*(V[1][J+5][K]-V[1][J-6][K]);
            DZ_V=A1*(V[1][J][K+1]-V[1][J][K])+A2*(V[1][J][K+2]-V[1][J][K-1])+A3*(V[1][J][K+3]-V[1][J][K-2])+A4*(V[1][J][K+4]-V[1][J][K-3])+A5*(V[1][J][K+5]-V[1][J][K-4])+A6*(V[1][J][K+6]-V[1][J][K-5]);
			
			P[0][J][K]=P[1][J][K]-(DT/DX)*C11[J][K]*DX_U-(DT/DZ)*C13[J][K]*DZ_V;			
	        
	        Q[0][J][K]=Q[1][J][K]-(DT/DX)*C13[J][K]*DX_U-(DT/DZ)*C33[J][K]*DZ_V;		
	        
	        S[0][J][K]=S[1][J][K]-(DT/DX)*C44[J][K]*DX_V-(DT/DZ)*C44[J][K]*DZ_U;
	        
	        tH[0][J][K]=tH[1][J][K]-(DT/DX)*C11[J][K]*DX_U-(DT/DZ)*c13[J][K]*DZ_V;
	        
	        tV[0][J][K]=tV[1][J][K]-(DT/DX)*c13[J][K]*DX_U-(DT/DZ)*C33[J][K]*DZ_V;
	    }
    
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
		{
	       DX_P=A1*(P[0][J][K]-P[0][J-1][K])+A2*(P[0][J+1][K]-P[0][J-2][K])+A3*(P[0][J+2][K]-P[0][J-3][K])+A4*(P[0][J+3][K]-P[0][J-4][K])+A5*(P[0][J+4][K]-P[0][J-5][K])+A6*(P[0][J+5][K]-P[0][J-6][K]);
           DZ_P=A1*(P[0][J][K]-P[0][J][K-1])+A2*(P[0][J][K+1]-P[0][J][K-2])+A3*(P[0][J][K+2]-P[0][J][K-3])+A4*(P[0][J][K+3]-P[0][J][K-4])+A5*(P[0][J][K+4]-P[0][J][K-5])+A6*(P[0][J][K+5]-P[0][J][K-6]);
           DX_Q=A1*(Q[0][J][K]-Q[0][J-1][K])+A2*(Q[0][J+1][K]-Q[0][J-2][K])+A3*(Q[0][J+2][K]-Q[0][J-3][K])+A4*(Q[0][J+3][K]-Q[0][J-4][K])+A5*(Q[0][J+4][K]-Q[0][J-5][K])+A6*(Q[0][J+5][K]-Q[0][J-6][K]);
           DZ_Q=A1*(Q[0][J][K]-Q[0][J][K-1])+A2*(Q[0][J][K+1]-Q[0][J][K-2])+A3*(Q[0][J][K+2]-Q[0][J][K-3])+A4*(Q[0][J][K+3]-Q[0][J][K-4])+A5*(Q[0][J][K+4]-Q[0][J][K-5])+A6*(Q[0][J][K+5]-Q[0][J][K-6]);
           
           DX_tH=A1*(tH[0][J][K]-tH[0][J-1][K])+A2*(tH[0][J+1][K]-tH[0][J-2][K])+A3*(tH[0][J+2][K]-tH[0][J-3][K])+A4*(tH[0][J+3][K]-tH[0][J-4][K])+A5*(tH[0][J+4][K]-tH[0][J-5][K])+A6*(tH[0][J+5][K]-tH[0][J-6][K]);
           DZ_tV=A1*(tV[0][J][K]-tV[0][J][K-1])+A2*(tV[0][J][K+1]-tV[0][J][K-2])+A3*(tV[0][J][K+2]-tV[0][J][K-3])+A4*(tV[0][J][K+3]-tV[0][J][K-4])+A5*(tV[0][J][K+4]-tV[0][J][K-5])+A6*(tV[0][J][K+5]-tV[0][J][K-6]);
		   
		   DX_S=A1*(S[0][J+1][K]-S[0][J][K])+A2*(S[0][J+2][K]-S[0][J-1][K])+A3*(S[0][J+3][K]-S[0][J-2][K])+A4*(S[0][J+4][K]-S[0][J-3][K])+A5*(S[0][J+5][K]-S[0][J-4][K])+A6*(S[0][J+6][K]-S[0][J-5][K]);
           DZ_S=A1*(S[0][J][K+1]-S[0][J][K])+A2*(S[0][J][K+2]-S[0][J][K-1])+A3*(S[0][J][K+3]-S[0][J][K-2])+A4*(S[0][J][K+4]-S[0][J][K-3])+A5*(S[0][J][K+5]-S[0][J][K-4])+A6*(S[0][J][K+6]-S[0][J][K-5]);
	       
	       U[0][J][K]=U[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_P-(1.0/DEN[J][K])*(DT/DZ)*DZ_S;
           V[0][J][K]=V[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_S-(1.0/DEN[J][K])*(DT/DZ)*DZ_Q;
           
           Up[0][J][K]=Up[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_tH;
           Vp[0][J][K]=Vp[1][J][K]-(1.0/DEN[J][K])*(DT/DZ)*DZ_tV;
           
           Us[0][J][K]=U[0][J][K]-Up[0][J][K];
		   Vs[0][J][K]=V[0][J][K]-Vp[0][J][K];
	    }
}
void reconstruct1(float A1,float A2,float A3,float A4,float A5,float A6,
                  float **C11,float **C13,float **C33,float **C44,
                  float **D11p,float **D12p,float **D13p,float **D21p,float **D22p,float **D23p,
                  float **D11s,float **D12s,float **D13s,float **D21s,float **D22s,float **D23s,
                  float **DEN,
                  float ***P,float ***Q,float ***S,float ***U,float ***Up,float ***Us,float ***V,float ***Vp,float ***Vs)
{
    int J,K;
    float   DX_P,DZ_P,DX_Q,DZ_Q,DX_S,DZ_S,DX_U,DZ_U,DX_V,DZ_V,DX_Up,DZ_Up,DX_Vp,DZ_Vp,DX_Us,DZ_Us,DX_Vs,DZ_Vs,DX_tH,DZ_tV; 
           
           DX_P=0.0;
	       DZ_P=0.0;
	       DX_Q=0.0;
	       DZ_Q=0.0;
		   DX_S=0.0;
		   DZ_S=0.0;
		   DX_U=0.0;
	       DZ_U=0.0;
	       DX_V=0.0;
	       DZ_V=0.0;
	       DX_Up=0.0;
	       DZ_Up=0.0;
	       DX_Vp=0.0;
	       DZ_Vp=0.0;
	       DX_Us=0.0;
	       DZ_Us=0.0;
	       DX_Vs=0.0;
	       DZ_Vs=0.0;
	       
	       DX_tH=0.0;
	       DZ_tV=0.0; 
	     
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
	    {
            DX_U=A1*(U[1][J+1][K]-U[1][J][K])+A2*(U[1][J+2][K]-U[1][J-1][K])+A3*(U[1][J+3][K]-U[1][J-2][K])+A4*(U[1][J+4][K]-U[1][J-3][K])+A5*(U[1][J+5][K]-U[1][J-4][K])+A6*(U[1][J+6][K]-U[1][J-5][K]);
            DZ_U=A1*(U[1][J][K]-U[1][J][K-1])+A2*(U[1][J][K+1]-U[1][J][K-2])+A3*(U[1][J][K+2]-U[1][J][K-3])+A4*(U[1][J][K+3]-U[1][J][K-4])+A5*(U[1][J][K+4]-U[1][J][K-5])+A6*(U[1][J][K+5]-U[1][J][K-6]);
            DX_V=A1*(V[1][J][K]-V[1][J-1][K])+A2*(V[1][J+1][K]-V[1][J-2][K])+A3*(V[1][J+2][K]-V[1][J-3][K])+A4*(V[1][J+3][K]-V[1][J-4][K])+A5*(V[1][J+4][K]-V[1][J-5][K])+A6*(V[1][J+5][K]-V[1][J-6][K]);
            DZ_V=A1*(V[1][J][K+1]-V[1][J][K])+A2*(V[1][J][K+2]-V[1][J][K-1])+A3*(V[1][J][K+3]-V[1][J][K-2])+A4*(V[1][J][K+4]-V[1][J][K-3])+A5*(V[1][J][K+5]-V[1][J][K-4])+A6*(V[1][J][K+6]-V[1][J][K-5]);
			
			P[0][J][K]=P[1][J][K]-(DT/DX)*C11[J][K]*DX_U-(DT/DZ)*C13[J][K]*DZ_V;			
	        Q[0][J][K]=Q[1][J][K]-(DT/DX)*C13[J][K]*DX_U-(DT/DZ)*C33[J][K]*DZ_V;		
	        S[0][J][K]=S[1][J][K]-(DT/DX)*C44[J][K]*DX_V-(DT/DZ)*C44[J][K]*DZ_U;
	    }
    
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
		{
	       DX_P=A1*(P[0][J][K]-P[0][J-1][K])+A2*(P[0][J+1][K]-P[0][J-2][K])+A3*(P[0][J+2][K]-P[0][J-3][K])+A4*(P[0][J+3][K]-P[0][J-4][K])+A5*(P[0][J+4][K]-P[0][J-5][K])+A6*(P[0][J+5][K]-P[0][J-6][K]);
           DZ_P=A1*(P[0][J][K]-P[0][J][K-1])+A2*(P[0][J][K+1]-P[0][J][K-2])+A3*(P[0][J][K+2]-P[0][J][K-3])+A4*(P[0][J][K+3]-P[0][J][K-4])+A5*(P[0][J][K+4]-P[0][J][K-5])+A6*(P[0][J][K+5]-P[0][J][K-6]);
           
           DX_Q=A1*(Q[0][J][K]-Q[0][J-1][K])+A2*(Q[0][J+1][K]-Q[0][J-2][K])+A3*(Q[0][J+2][K]-Q[0][J-3][K])+A4*(Q[0][J+3][K]-Q[0][J-4][K])+A5*(Q[0][J+4][K]-Q[0][J-5][K])+A6*(Q[0][J+5][K]-Q[0][J-6][K]);
           DZ_Q=A1*(Q[0][J][K]-Q[0][J][K-1])+A2*(Q[0][J][K+1]-Q[0][J][K-2])+A3*(Q[0][J][K+2]-Q[0][J][K-3])+A4*(Q[0][J][K+3]-Q[0][J][K-4])+A5*(Q[0][J][K+4]-Q[0][J][K-5])+A6*(Q[0][J][K+5]-Q[0][J][K-6]);
           
		   DX_S=A1*(S[0][J+1][K]-S[0][J][K])+A2*(S[0][J+2][K]-S[0][J-1][K])+A3*(S[0][J+3][K]-S[0][J-2][K])+A4*(S[0][J+4][K]-S[0][J-3][K])+A5*(S[0][J+5][K]-S[0][J-4][K])+A6*(S[0][J+6][K]-S[0][J-5][K]);
           DZ_S=A1*(S[0][J][K+1]-S[0][J][K])+A2*(S[0][J][K+2]-S[0][J][K-1])+A3*(S[0][J][K+3]-S[0][J][K-2])+A4*(S[0][J][K+4]-S[0][J][K-3])+A5*(S[0][J][K+5]-S[0][J][K-4])+A6*(S[0][J][K+6]-S[0][J][K-5]);
	       
	       U[0][J][K]=U[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_P-(1.0/DEN[J][K])*(DT/DZ)*DZ_S;
           V[0][J][K]=V[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_S-(1.0/DEN[J][K])*(DT/DZ)*DZ_Q;
           
           Up[0][J][K]=Up[1][J][K]-(D11p[J][K]/DEN[J][K])*(DT/DX)*DX_P-(D12p[J][K]/DEN[J][K])*(DT/DX)*DX_Q-(D13p[J][K]/DEN[J][K])*(DT/DZ)*DZ_S;
           Vp[0][J][K]=Vp[1][J][K]-(D21p[J][K]/DEN[J][K])*(DT/DZ)*DZ_P-(D22p[J][K]/DEN[J][K])*(DT/DZ)*DZ_Q-(D23p[J][K]/DEN[J][K])*(DT/DX)*DX_S;
           
           Us[0][J][K]=Us[1][J][K]-(D11s[J][K]/DEN[J][K])*(DT/DX)*DX_P-(D12s[J][K]/DEN[J][K])*(DT/DX)*DX_Q-(D13s[J][K]/DEN[J][K])*(DT/DZ)*DZ_S;
           Vs[0][J][K]=Vs[1][J][K]-(D21s[J][K]/DEN[J][K])*(DT/DZ)*DZ_P-(D22s[J][K]/DEN[J][K])*(DT/DZ)*DZ_Q-(D23s[J][K]/DEN[J][K])*(DT/DX)*DX_S;
	    }
}
void adjoint(float A1,float A2,float A3,float A4,float A5,float A6,float **C11,float **C13,float **C33,float **C44,float **c13,float **DEN,float **dx, float **dz,float ***P,float ***PX,float ***PZ,float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,float ***tH,float ***tHX,float ***tHZ,float ***tV,float ***tVX,float ***tVZ,float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs)
 {
        int J,K;
        float   DX_P,DZ_P,DX_Q,DZ_Q,DX_S,DZ_S,DX_U,DZ_U,DX_V,DZ_V,DX_tH,DZ_tV; 
           
           DX_P=0.0;
	       DZ_P=0.0;
	       DX_Q=0.0;
	       DZ_Q=0.0;
		   DX_S=0.0;
		   DZ_S=0.0;
		   DX_U=0.0;
	       DZ_U=0.0;
	       DX_V=0.0;
	       DZ_V=0.0;
	       DX_tH=0.0;
	       DZ_tV=0.0; 
/* .......................计算速度分量..............................*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
		{
	       DX_P=A1*(P[1][J][K]-P[1][J-1][K])+A2*(P[1][J+1][K]-P[1][J-2][K])+A3*(P[1][J+2][K]-P[1][J-3][K])+A4*(P[1][J+3][K]-P[1][J-4][K])+A5*(P[1][J+4][K]-P[1][J-5][K])+A6*(P[1][J+5][K]-P[1][J-6][K]);
           DZ_P=A1*(P[1][J][K]-P[1][J][K-1])+A2*(P[1][J][K+1]-P[1][J][K-2])+A3*(P[1][J][K+2]-P[1][J][K-3])+A4*(P[1][J][K+3]-P[1][J][K-4])+A5*(P[1][J][K+4]-P[1][J][K-5])+A6*(P[1][J][K+5]-P[1][J][K-6]);
           DX_Q=A1*(Q[1][J][K]-Q[1][J-1][K])+A2*(Q[1][J+1][K]-Q[1][J-2][K])+A3*(Q[1][J+2][K]-Q[1][J-3][K])+A4*(Q[1][J+3][K]-Q[1][J-4][K])+A5*(Q[1][J+4][K]-Q[1][J-5][K])+A6*(Q[1][J+5][K]-Q[1][J-6][K]);
           DZ_Q=A1*(Q[1][J][K]-Q[1][J][K-1])+A2*(Q[1][J][K+1]-Q[1][J][K-2])+A3*(Q[1][J][K+2]-Q[1][J][K-3])+A4*(Q[1][J][K+3]-Q[1][J][K-4])+A5*(Q[1][J][K+4]-Q[1][J][K-5])+A6*(Q[1][J][K+5]-Q[1][J][K-6]);
           
           DX_tH=A1*(tH[1][J][K]-tH[1][J-1][K])+A2*(tH[1][J+1][K]-tH[1][J-2][K])+A3*(tH[1][J+2][K]-tH[1][J-3][K])+A4*(tH[1][J+3][K]-tH[1][J-4][K])+A5*(tH[1][J+4][K]-tH[1][J-5][K])+A6*(tH[1][J+5][K]-tH[1][J-6][K]);
           DZ_tV=A1*(tV[1][J][K]-tV[1][J][K-1])+A2*(tV[1][J][K+1]-tV[1][J][K-2])+A3*(tV[1][J][K+2]-tV[1][J][K-3])+A4*(tV[1][J][K+3]-tV[1][J][K-4])+A5*(tV[1][J][K+4]-tV[1][J][K-5])+A6*(tV[1][J][K+5]-tV[1][J][K-6]);
		   
		   DX_S=A1*(S[1][J+1][K]-S[1][J][K])+A2*(S[1][J+2][K]-S[1][J-1][K])+A3*(S[1][J+3][K]-S[1][J-2][K])+A4*(S[1][J+4][K]-S[1][J-3][K])+A5*(S[1][J+5][K]-S[1][J-4][K])+A6*(S[1][J+6][K]-S[1][J-5][K]);
           DZ_S=A1*(S[1][J][K+1]-S[1][J][K])+A2*(S[1][J][K+2]-S[1][J][K-1])+A3*(S[1][J][K+3]-S[1][J][K-2])+A4*(S[1][J][K+4]-S[1][J][K-3])+A5*(S[1][J][K+5]-S[1][J][K-4])+A6*(S[1][J][K+6]-S[1][J][K-5]);
	       
	       U[0][J][K]=U[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_P-(1.0/DEN[J][K])*(DT/DZ)*DZ_S;
           V[0][J][K]=V[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_S-(1.0/DEN[J][K])*(DT/DZ)*DZ_Q;
           
           Up[0][J][K]=Up[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_tH;
           Vp[0][J][K]=Vp[1][J][K]-(1.0/DEN[J][K])*(DT/DZ)*DZ_tV;
	    }
 


/*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
          
	    }

  
/*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }    
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }    
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }   
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {     
	      UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	   for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }
 /*.....................四周降阶处理............................................*/
    if((J=0)&&(K=0))
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
        }
    if((J=0)&&(K=NZ-1))
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
		}
    if((J=NX-1)&&(K=NZ-1))
		{
          
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
		}
    if((J=NX-1)&&(K=0))
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
		}

    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
          
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
          
		}
        if(K=NZ-1)
	    { 
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
		}
        if(J=NX-1)
	    { 
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-(tH[1][J][K]-tH[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[0][J][K]+UpZ[0][J][K];

          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-(tV[1][J][K]-tV[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[0][J][K]+VpZ[0][J][K];
	    }
	}
    for(J=0;J<=NX-1;J++)
		for(K=0;K<=NZ-1;K++)
	    {
	     Us[0][J][K]=U[0][J][K]-Up[0][J][K];
	     Vs[0][J][K]=V[0][J][K]-Vp[0][J][K];
	    }
    
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
	    {
            DX_U=A1*(U[0][J+1][K]-U[0][J][K])+A2*(U[0][J+2][K]-U[0][J-1][K])+A3*(U[0][J+3][K]-U[0][J-2][K])+A4*(U[0][J+4][K]-U[0][J-3][K])+A5*(U[0][J+5][K]-U[0][J-4][K])+A6*(U[0][J+6][K]-U[0][J-5][K]);
            DZ_U=A1*(U[0][J][K]-U[0][J][K-1])+A2*(U[0][J][K+1]-U[0][J][K-2])+A3*(U[0][J][K+2]-U[0][J][K-3])+A4*(U[0][J][K+3]-U[0][J][K-4])+A5*(U[0][J][K+4]-U[0][J][K-5])+A6*(U[0][J][K+5]-U[0][J][K-6]);
            DX_V=A1*(V[0][J][K]-V[0][J-1][K])+A2*(V[0][J+1][K]-V[0][J-2][K])+A3*(V[0][J+2][K]-V[0][J-3][K])+A4*(V[0][J+3][K]-V[0][J-4][K])+A5*(V[0][J+4][K]-V[0][J-5][K])+A6*(V[0][J+5][K]-V[0][J-6][K]);
            DZ_V=A1*(V[0][J][K+1]-V[0][J][K])+A2*(V[0][J][K+2]-V[0][J][K-1])+A3*(V[0][J][K+3]-V[0][J][K-2])+A4*(V[0][J][K+4]-V[0][J][K-3])+A5*(V[0][J][K+5]-V[0][J][K-4])+A6*(V[0][J][K+6]-V[0][J][K-5]);
			
			P[0][J][K]=P[1][J][K]-(DT/DX)*C11[J][K]*DX_U-(DT/DZ)*C13[J][K]*DZ_V;			
	        
	        Q[0][J][K]=Q[1][J][K]-(DT/DX)*C13[J][K]*DX_U-(DT/DZ)*C33[J][K]*DZ_V;		
	        
	        S[0][J][K]=S[1][J][K]-(DT/DX)*C44[J][K]*DX_V-(DT/DZ)*C44[J][K]*DZ_U;
	        
	        tH[0][J][K]=tH[1][J][K]-(DT/DX)*C11[J][K]*DX_U-(DT/DZ)*c13[J][K]*DZ_V;
	        
	        tV[0][J][K]=tV[1][J][K]-(DT/DX)*c13[J][K]*DX_U-(DT/DZ)*C33[J][K]*DZ_V;
	    }
          
  
/*计算边界处各应力分量值*/
 /*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {  
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	    }
  
    /*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {

           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	    }   
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	    }   
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {	   
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    { 
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	    }  
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];

	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];

	    }
 /*.....................边界条件四周降阶处理............................................*/

    if((J=0)&&(K=0))
	 {	
           
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];                                      
       }
    if((J=0)&&(K=NZ-1))
      {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
      }
    if((J=NX-1)&&(K=NZ-1))
      {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
      }

    if((J=NX-1)&&(K=0))
      {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
      }
   
         
    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
          
		}
       if(K=NZ-1)
	    { 
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*c13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
		}
       if(J=NX-1)
	   { 
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
           
           tHX[0][J][K]=((1-0.5*DT*dx[J][K])*tHX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   tHZ[0][J][K]=((1-0.5*DT*dz[J][K])*tHZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*c13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tH[0][J][K]=tHX[0][J][K]+tHZ[0][J][K];

           tVX[0][J][K]=((1-0.5*DT*dx[J][K])*tVX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   tVZ[0][J][K]=((1-0.5*DT*dz[J][K])*tVZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           tV[0][J][K]=tVX[0][J][K]+tVZ[0][J][K];
	   }
	}
 }
 void adjoint1(float A1,float A2,float A3,float A4,float A5,float A6,
               float **C11,float **C13,float **C33,float **C44,
               float **D11p,float **D12p,float **D13p,float **D21p,float **D22p,float **D23p,
               float **D11s,float **D12s,float **D13s,float **D21s,float **D22s,float **D23s,float **DEN,
               float **dx, float **dz,float ***P,float ***PX,float ***PZ,
               float ***Q,float ***QX,float ***QZ,float ***S,float ***SX,float ***SZ,
               float ***U,float ***UX,float ***UZ,float ***Up,float ***UpX,float ***UpZ,float ***Us,float ***UsX,float ***UsZ,
               float ***V,float ***VX,float ***VZ,float ***Vp,float ***VpX,float ***VpZ,float ***Vs,float ***VsX,float ***VsZ)
 {
        int J,K;
        float   DX_P,DZ_P,DX_Q,DZ_Q,DX_S,DZ_S,DX_U,DZ_U,DX_V,DZ_V,DX_tH,DZ_tV; 
           
           DX_P=0.0;
	       DZ_P=0.0;
	       DX_Q=0.0;
	       DZ_Q=0.0;
		   DX_S=0.0;
		   DZ_S=0.0;
		   DX_U=0.0;
	       DZ_U=0.0;
	       DX_V=0.0;
	       DZ_V=0.0;
	       DX_tH=0.0;
	       DZ_tV=0.0; 
/* .......................计算速度分量..............................*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
		{
	       DX_P=A1*(P[1][J][K]-P[1][J-1][K])+A2*(P[1][J+1][K]-P[1][J-2][K])+A3*(P[1][J+2][K]-P[1][J-3][K])+A4*(P[1][J+3][K]-P[1][J-4][K])+A5*(P[1][J+4][K]-P[1][J-5][K])+A6*(P[1][J+5][K]-P[1][J-6][K]);
           DZ_P=A1*(P[1][J][K]-P[1][J][K-1])+A2*(P[1][J][K+1]-P[1][J][K-2])+A3*(P[1][J][K+2]-P[1][J][K-3])+A4*(P[1][J][K+3]-P[1][J][K-4])+A5*(P[1][J][K+4]-P[1][J][K-5])+A6*(P[1][J][K+5]-P[1][J][K-6]);
           
           DX_Q=A1*(Q[1][J][K]-Q[1][J-1][K])+A2*(Q[1][J+1][K]-Q[1][J-2][K])+A3*(Q[1][J+2][K]-Q[1][J-3][K])+A4*(Q[1][J+3][K]-Q[1][J-4][K])+A5*(Q[1][J+4][K]-Q[1][J-5][K])+A6*(Q[1][J+5][K]-Q[1][J-6][K]);
           DZ_Q=A1*(Q[1][J][K]-Q[1][J][K-1])+A2*(Q[1][J][K+1]-Q[1][J][K-2])+A3*(Q[1][J][K+2]-Q[1][J][K-3])+A4*(Q[1][J][K+3]-Q[1][J][K-4])+A5*(Q[1][J][K+4]-Q[1][J][K-5])+A6*(Q[1][J][K+5]-Q[1][J][K-6]);
           
		   DX_S=A1*(S[1][J+1][K]-S[1][J][K])+A2*(S[1][J+2][K]-S[1][J-1][K])+A3*(S[1][J+3][K]-S[1][J-2][K])+A4*(S[1][J+4][K]-S[1][J-3][K])+A5*(S[1][J+5][K]-S[1][J-4][K])+A6*(S[1][J+6][K]-S[1][J-5][K]);
           DZ_S=A1*(S[1][J][K+1]-S[1][J][K])+A2*(S[1][J][K+2]-S[1][J][K-1])+A3*(S[1][J][K+3]-S[1][J][K-2])+A4*(S[1][J][K+4]-S[1][J][K-3])+A5*(S[1][J][K+5]-S[1][J][K-4])+A6*(S[1][J][K+6]-S[1][J][K-5]);
	       
	       U[0][J][K]=U[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_P-(1.0/DEN[J][K])*(DT/DZ)*DZ_S;
           V[0][J][K]=V[1][J][K]-(1.0/DEN[J][K])*(DT/DX)*DX_S-(1.0/DEN[J][K])*(DT/DZ)*DZ_Q;
           
           Up[0][J][K]=Up[1][J][K]-(D11p[J][K]/DEN[J][K])*(DT/DX)*DX_P-(D12p[J][K]/DEN[J][K])*(DT/DX)*DX_Q-(D13p[J][K]/DEN[J][K])*(DT/DZ)*DZ_S;
           Vp[0][J][K]=Vp[1][J][K]-(D21p[J][K]/DEN[J][K])*(DT/DZ)*DZ_P-(D22p[J][K]/DEN[J][K])*(DT/DZ)*DZ_Q-(D23p[J][K]/DEN[J][K])*(DT/DX)*DX_S;
           
           Us[0][J][K]=Us[1][J][K]-(D11s[J][K]/DEN[J][K])*(DT/DX)*DX_P-(D12s[J][K]/DEN[J][K])*(DT/DX)*DX_Q-(D13s[J][K]/DEN[J][K])*(DT/DZ)*DZ_S;
           Vs[0][J][K]=Vs[1][J][K]-(D21s[J][K]/DEN[J][K])*(DT/DZ)*DZ_P-(D22s[J][K]/DEN[J][K])*(DT/DZ)*DZ_Q-(D23s[J][K]/DEN[J][K])*(DT/DX)*DX_S;
	    }
 


/*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }

  
/*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }    
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }    
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }   
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {     
	      UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	   for(K=NZ-PML;K<=NZ-2;K++)
	    {
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }
 /*.....................四周降阶处理............................................*/
    if((J=0)&&(K=0))
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
        }
    if((J=0)&&(K=NZ-1))
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
		}
    if((J=NX-1)&&(K=NZ-1))
		{
          
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
		}
    if((J=NX-1)&&(K=0))
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
		}

    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
          
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
          
		}
        if(K=NZ-1)
	    { 
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K])/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23p[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K]-(S[1][J+1][K]-S[1][J][K])*(D23s[J][K]/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
		}
        if(J=NX-1)
	    { 
          UX[0][J][K]=((1-0.5*DT*dx[J][K])*UX[1][J][K]-(P[1][J][K]-P[1][J-1][K])*(1.0/DEN[J][K])*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UZ[0][J][K]=((1-0.5*DT*dz[J][K])*UZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          U[0][J][K]=UX[0][J][K]+UZ[0][J][K];

          VX[0][J][K]=((1-0.5*DT*dx[J][K])*VX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VZ[0][J][K]=((1-0.5*DT*dz[J][K])*VZ[1][J][K]-(Q[1][J][K]-Q[1][J][K-1])*(1.0/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          V[0][J][K]=VX[0][J][K]+VZ[0][J][K];
          
          
          UpX[0][J][K]=((1-0.5*DT*dx[J][K])*UpX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12p[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UpZ[0][J][K]=((1-0.5*DT*dz[J][K])*UpZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13p[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Up[0][J][K]=UpX[1][J][K]+UpZ[1][J][K];
          
          VpX[0][J][K]=((1-0.5*DT*dx[J][K])*VpX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VpZ[0][J][K]=((1-0.5*DT*dz[J][K])*VpZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21p[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22p[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vp[0][J][K]=VpX[1][J][K]+VpZ[1][J][K]; 
          
          UsX[0][J][K]=((1-0.5*DT*dx[J][K])*UsX[1][J][K]-((P[1][J][K]-P[1][J-1][K])*(D11s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J-1][K])*(D12s[J][K]/DEN[J][K]))*DT/DX)/(1+0.5*DT*dx[J][K]);
		  UsZ[0][J][K]=((1-0.5*DT*dz[J][K])*UsZ[1][J][K]-(S[1][J][K+1]-S[1][J][K])*(D13s[J][K]/DEN[J][K])*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Us[0][J][K]=UsX[0][J][K]+UsZ[0][J][K];
          
          VsX[0][J][K]=((1-0.5*DT*dx[J][K])*VsX[1][J][K])/(1+0.5*DT*dx[J][K]);
		  VsZ[0][J][K]=((1-0.5*DT*dz[J][K])*VsZ[1][J][K]-((P[1][J][K]-P[1][J][K-1])*(D21s[J][K]/DEN[J][K])+(Q[1][J][K]-Q[1][J][K-1])*(D22s[J][K]/DEN[J][K]))*DT/DZ)/(1+0.5*DT*dz[J][K]);
          Vs[0][J][K]=VsX[0][J][K]+VsZ[0][J][K]; 
          
	    }
	}
    
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=PML;K<=NZ-PML-1;K++)
	    {
            DX_U=A1*(U[0][J+1][K]-U[0][J][K])+A2*(U[0][J+2][K]-U[0][J-1][K])+A3*(U[0][J+3][K]-U[0][J-2][K])+A4*(U[0][J+4][K]-U[0][J-3][K])+A5*(U[0][J+5][K]-U[0][J-4][K])+A6*(U[0][J+6][K]-U[0][J-5][K]);
            DZ_U=A1*(U[0][J][K]-U[0][J][K-1])+A2*(U[0][J][K+1]-U[0][J][K-2])+A3*(U[0][J][K+2]-U[0][J][K-3])+A4*(U[0][J][K+3]-U[0][J][K-4])+A5*(U[0][J][K+4]-U[0][J][K-5])+A6*(U[0][J][K+5]-U[0][J][K-6]);
            DX_V=A1*(V[0][J][K]-V[0][J-1][K])+A2*(V[0][J+1][K]-V[0][J-2][K])+A3*(V[0][J+2][K]-V[0][J-3][K])+A4*(V[0][J+3][K]-V[0][J-4][K])+A5*(V[0][J+4][K]-V[0][J-5][K])+A6*(V[0][J+5][K]-V[0][J-6][K]);
            DZ_V=A1*(V[0][J][K+1]-V[0][J][K])+A2*(V[0][J][K+2]-V[0][J][K-1])+A3*(V[0][J][K+3]-V[0][J][K-2])+A4*(V[0][J][K+4]-V[0][J][K-3])+A5*(V[0][J][K+5]-V[0][J][K-4])+A6*(V[0][J][K+6]-V[0][J][K-5]);
			
			P[0][J][K]=P[1][J][K]-(DT/DX)*C11[J][K]*DX_U-(DT/DZ)*C13[J][K]*DZ_V;			
	        Q[0][J][K]=Q[1][J][K]-(DT/DX)*C13[J][K]*DX_U-(DT/DZ)*C33[J][K]*DZ_V;		
	        S[0][J][K]=S[1][J][K]-(DT/DX)*C44[J][K]*DX_V-(DT/DZ)*C44[J][K]*DZ_U;
	    }
          
  
/*计算边界处各应力分量值*/
 /*左边界*/
    for(J=1;J<=PML-1;J++)	 
	   for(K=PML;K<=NZ-PML-1;K++)
	    {  
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }
  
    /*右边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
       for(K=PML;K<=NZ-PML-1;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }
/*上边界*/

    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {

           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }   
/*下边界*/
    for(J=PML;J<=NX-PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }   
/*左上角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=1;K<=PML-1;K++)
	    {	   
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }

/*右上边界*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=1;K<=PML-1;K++)
	    { 
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }  
    
/*左下角*/
    for(J=1;J<=PML-1;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }

/*右下角*/
    for(J=NX-PML;J<=NX-2;J++)	 
	    for(K=NZ-PML;K<=NZ-2;K++)
	    {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }
 /*.....................边界条件四周降阶处理............................................*/

    if((J=0)&&(K=0))
	 {	
           
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];                                  
       }
    if((J=0)&&(K=NZ-1))
      {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
      }
    if((J=NX-1)&&(K=NZ-1))
      {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
      }

    if((J=NX-1)&&(K=0))
      {
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
      }
   
         
    for(J=1;J<=NX-2;J++)
	{
		if(K=0)
		{
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
		}
       if(K=NZ-1)
	    { 
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K])/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	    }
	}
    for(K=1;K<=NZ-2;K++)
    {
		if(J=0)
		{
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C11[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K]-(U[0][J+1][K]-U[0][J][K])*C13[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K])/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
		}
       if(J=NX-1)
	   { 
           PX[0][J][K]=((1-0.5*DT*dx[J][K])*PX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   PZ[0][J][K]=((1-0.5*DT*dz[J][K])*PZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C13[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           P[0][J][K]=PX[0][J][K]+PZ[0][J][K];

           QX[0][J][K]=((1-0.5*DT*dx[J][K])*QX[1][J][K])/(1+0.5*DT*dx[J][K]);
		   QZ[0][J][K]=((1-0.5*DT*dz[J][K])*QZ[1][J][K]-(V[0][J][K+1]-V[0][J][K])*C33[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           Q[0][J][K]=QX[0][J][K]+QZ[0][J][K];

           SX[0][J][K]=((1-0.5*DT*dx[J][K])*SX[1][J][K]-(V[0][J][K]-V[0][J-1][K])*C44[J][K]*(DT/DX))/(1+0.5*DT*dx[J][K]);
	       SZ[0][J][K]=((1-0.5*DT*dz[J][K])*SZ[1][J][K]-(U[0][J][K]-U[0][J][K-1])*C44[J][K]*(DT/DZ))/(1+0.5*DT*dz[J][K]);
           S[0][J][K]=SX[0][J][K]+SZ[0][J][K];
	   }
	}
 }