/* This is John Lujano's research with Dr. Tausch for the paper entitled
An Optimization Method for Moving Interface Problems Goverend by the Heat Equation*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define ONETHRD 0.333333333333333333333333333333
#define M_SQPI 1.77245385090552
#define M_SQPINV 0.564189583547756287

/*
 * Compute the radius at a given time t
 * Parameters 
 *    t    time
 * Return
 *    g    radius
 */

double rKnown(double t){
  
  double g;
  
  g=t;
  
  return g;
}

/*
 * Compute the derivative of the radius with respect to time
 * Parameters 
 *    t    time
 * Return
 *    g    dr/dt
 */

double dRKnown(double t){
  
  double g;
  
  g=1;
  
  return g;
}

/*
 * Give the xo for a test u(x,t)
 * Parameters 
 * Return
 *    g    xo
 */

double xofun(){
  
  double g;
  
  g=2.0;
  
  return g;
}

/*
 * Give the to for a test u(x,t)
 * Parameters 
 * Return
 *    g    to
 */

double tofun(){
  
  double g;
  
  g=2.0;
  
  return g;
}

/*
 * Give the value for a test u(x,t)
 * Parameters 
 *    x    the radius at a given time
 *    t    the time
 * Return
 *    g    u
 */

double u(double x,double t){
  
  double g, xo=xofun() , to=tofun();
  
  g=exp(-(x-xo)*(x-xo)/(4.0*(t+to)))/sqrt(t+to);
  
  return g;
}

/*
 * Give the value for a test du(x,t)
 * Parameters 
 *    x    the radius at a given time
 *    t    the time
 * Return
 *    g    du
 */

double du(double x,double t){

  double g, xo=xofun() , to=tofun();

  g = exp(-(x-xo)*(x-xo)/(4.0*(t+to)))/sqrt(t+to);
  g *= -(x-xo)/(2.0*(t+to));

  return g;
}

/*
 * Give the value for a test uo(x,t)
 * Parameters 
 *    x    the radius at a given time
 * Return
 *    g    uo
 */

double uo(double x){
  
  double g;
  
  //g=u(x,0.0);
  //g=0;
  g=1.0/(1.0+x);
  
  return g;
}

/*
 * Give the value for a test duo(x,t)
 * Parameters 
 *    x    the radius at a given time
 * Return
 *    g    duo
 */

double duo(double x){
  
  double g;
  
  //g=du(x,0.0);
  g=-1.0/((1.0+x)*(1.0+x));
  
  return g;
}

/*
 * Give the value for a test q(t)
 * Parameters 
 *    t    the time
 * Return
 *    q    q
 */

double qKnown(double t){

  double r,r1,q;
  
  r  = rKnown(t);
  r1 = dRKnown(t);
  q = du(r, t) + u(r, t)*r1;

  return q;
}

/*
 * Give the value for kv-> part of the green's function
 * Parameters 
 *    x    radius at the given timestep
 *    xtau the radius spanning up to the given timestep
 *    t    the time at our timestep
 *    tau  the time spanning up to the given timestep
 * Return
 *    g    duo
 */

double kv(double x, double xtau,double t, double tau){
  double g;

  g=0.5*M_SQPINV*exp(-(x-xtau)*(x-xtau)/(4.0*(t-tau)));
  
  return g;
}

/*Used to test the single layer potential with trivial kv*/
/*
 * Give the value for kv-> part of the green's function
 * Parameters 
 *    x    radius at the given timestep
 *    xtau the radius spanning up to the given timestep
 *    t    the time at our timestep
 *    tau  the time spanning up to the given timestep
 * Return
 *    g    duo
 */

double kv2(double x,double xo,double t,double to){
  double g;
  
  g=1.0;
  
  return g;
}

/*
 * Give the value for kv at the current timestep which approaches the value
 * Parameters 
 * Return
 *    g    kvtntn
 */

double kvtntn(){
  double g;
  
  g=0.5*M_SQPINV;
  
  return g;
}

/*Used to test the single layer potential with trivial kv*/
/*
 * Parameters 
 * Return
 *    g    kvtntn
 */

double kvtntn2(){
  double g;
  
  g=1.0;
  
  return g;
}

/*
 * Give the value for kd-> part of the green's function
 * Parameters 
 *    x    radius at the given timestep
 *    xtau the radius spanning up to the given timestep
 *    t    the time at our timestep
 *    tau  the time spanning up to the given timestep
 * Return
 *    g    kd
 */

double kd(double x,double xtau,double t,double tau){
  double g;
  
  g=exp(-(x-xtau)*(x-xtau)/(4.0*(t-tau)))*(x-xtau)*0.25*M_SQPINV/((t-tau));
  
  return g;
}

/*
 * Give the value for kd-> part of the green's function
 * At the current timestep
 * Parameters 
 *    t    the time at our timestep
 * Return
 *    g    kdtntn
 */

double kdtntn(double dr, double t){
  double g;
  
  g=dr*0.25*M_SQPINV;
  
  return g;
}

/*
 * Algorithm for finding the single layer potential in the integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 *    r    radius time history
 *    q    history of q values
 * Return
 *         final value
 */

double sLay(int n,double h,double *r,double *q){
  double tn=n*h;
  int i;
  double sum=0.5*kv(r[n],r[0],tn,0.0)*q[0]/sqrt(tn-0.0);
  
  for (i=1; i<n; i++){
    sum += kv(r[n],r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
  }
  
  return sum*h;
}

/*Same as previous, but the trivial kv2*/

double sLay2(int n,double h,double *r,double *q){
  double tn=n*h;
  int i;
  double sum=0.5*kv2(r[n],r[0],tn,0.0)*q[0]/sqrt(tn-0.0);
  
  for (i=1; i<n; i++){
    sum += kv2(r[n],r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
  }
  
  return sum*h;
}

/*
 * Algorithm for finding the double layer potential in the integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 *    r    radius time history
 *    u    history of u values
 * Return
 *         final value
 */

double dLay(int n,double h,double *r,double *u){
  double tn=n*h;
  int i;
  double sum=0.5*kd(r[n],r[0],tn,0.0)*u[0]/sqrt(tn-0.0);

  for (i=1; i<n; i++){
    sum += kd(r[n],r[i],tn,i*h)*u[i]/sqrt(tn-i*h);
  }
  
  return sum*h;
}

/*
 * Algorithm for finding the mu value in the discretization
 * Parameters 
 *    n    timestep index
 *    h    dt
 * Return
 *         final value
 */

double mun(int n,double h){
  double tn=n*h;
  int i;
  double sum=0.5*h/sqrt(tn-0.0);
  
  for (i=1; i<n; i++){
    sum += h/sqrt(tn-i*h);
  }
  
  return 2.0*sqrt(tn)-sum;
}

/*
 * Algorithm for finding the initial potential value in the integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 *    r    history of r values
 * Return
 *         final value
 */

double initPot(int n, double h, double *r){
  int numInts=10000;
  double b=10.0, z, tn=(n)*h, s4t=sqrt(4.0*tn),a=r[n]/s4t;

  if (a<b){
    a=a;
  }
  else{
    a=b;
  }

  double delta=(b+a)/numInts;
  double sum=0.5*(exp(-a*a)*uo(r[n]-a*s4t)+exp(-b*b)*uo(r[n]+b*s4t));
  int i;
  
  for(i=1; i<numInts; i++){
    z = -a+delta*i;
    sum += exp(-z*z)*uo(r[n]+z*s4t);
  }
  
  return sum*delta*M_SQPINV;
}

/*
 * Used to test one of the integration schemes
 * Parameters 
 *    t    time
 * Return
 *         final value
 */

double qTestExact(double t){
  double g;
  
  g=exp(t)*sqrt(M_PI)*erf(sqrt(t));
  
  return g;
}

/*
 * Used to test one of the integration schemes
 * Parameters 
 *    t    time
 * Return
 *         final value
 */

double qEval(double t){
  double ans=(1+t*t+t*t*t*t);
  
  return ans;
}

/*
 * Used to test one of the integration schemes
 * Parameters 
 *    t    time
 * Return
 *         final value
 */

double q2Exact(double t){
  double ans=M_PI*(1+3.0*t*t/8.0+35.0*t*t*t*t/128.0);
  
  return ans;
}

/*
 * Algorithm for finding the single layer potential in the adjoint integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 *    r    radius time history
 *    q    history of q values
 * Return
 *         final value
 */

double sLayAdj(int n, double h, double *r, double *qt){
  double ktntjsum=0.0 , tn=n*h ;
  int i;
  
  for (i=1;i<n;i++){
    ktntjsum+=kv(r[n],r[i],tn,i*h)*qt[i]/(sqrt(tn-i*h));
  }
  
  return h*ktntjsum;
}

/*Same as previous, but the trivial kv2*/

double sLayAdj2(int n, double h, double *r, double *qt){
  double ktntjsum=0.0,tn=n*h;
  int i;
  
  for (i=1;i<n;i++){
    ktntjsum+=kv2(r[n],r[i],tn,i*h)*qt[i]/(sqrt(tn-i*h));
  }
  
  return h*ktntjsum;
}

/*
 * Algorithm for finding the coefficient of q[0] in the adjoint
 * Parameters 
 *    n    timestep index
 *    h    dt
 * Return
 *         final value
 */

double qto(int n, double h){
  double sum=M_PI/2.0,tn=n*h;
  int i;
  
  for (i=1;i<n;i++){
    sum-=h*sqrt(tn-i*h)/(tn*sqrt(i*h));
  }
  
  return sum;
}

/*
 * Algorithm for finding the mu value in the adjoint integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 * Return
 *         final value
 */

double munAdj(int n, double h){
  double tn=n*h , sum=0.5*M_PI*tn;
  int i;
  
  for (i=1;i<n;i++){
    sum-=h*sqrt(i*h)/(sqrt(tn-i*h));
  }
  
  return sum/sqrt(tn);
}

/*
 * Algorithm for solving the state equation
 * Parameters 
 *    n    number of intervals
 *    tmax maximum times
 *    r    history of radius values
 * Return
 *    qs   the state q values
 */

double * solveState(int n, double tmax, double *r){
  int i , j;
  double *qs , *us;
  double dro , duoVal , dt,dr;
  
  qs=(double*)malloc( (n+1)*sizeof(double) );
  us=(double*)malloc( (n+1)*sizeof(double) );
  
  dt=tmax/(n);
  
  for (i=0;i<=n;i++){
    //us[i]=u(r[i],i*dt);
    us[i]=1.0;
  }
  
  //first order
  //dro=(r[1]-r[0])/dt;
  //third order
  dro=(-11.0*r[0]/6.0+3.0*r[1]-3.0*r[2]/2.0+1.0*r[3]/3.0)/dt;
  duoVal=duo(r[0]);
  
  //qs[0]=qKnown(0);
  qs[0]=dro+duoVal;
  //qs[0]=0;
  
  for (j=1;j<=n;j++){
    if (j!=n){
      dr=(r[j+1]-r[j-1])/(2.0*dt);
    }
    else{
      dr=(r[j]-r[j-1])/dt;
    }
    qs[j]=(-0.5*us[j]-sLay(j,dt,r,qs)+dLay(j,dt,r,us)+kdtntn(dr,j*dt)*mun(j,dt)*us[j]+initPot(j,dt,r))/(mun(j,dt)*kvtntn());
    //printf("first %f slay %f dlay %f mun %f init %f\n",0.5*us[j],sLay(j,dt,r,qs),dLay(j,dt,r,us),kdtntn(dr,j*dt)*mun(j,dt)*us[j],initPot(j,dt,r));
  }
  
  return qs;
}

/*
 * Algorithm for solving the adjoint equation
 * Parameters 
 *    n    timestep index
 *    tmax maximum times
 *    r    history of radius values
 *    qs   state q values
 * Return
 *    qa   the adjoint q values
 */

double * solveAdj(int n, double tmax, double *r ,double *qs){
  int i , j,k;
  double *qa , *us, *rrev, *qarev;
  double dt,dr;
  
  qa=(double*)malloc( (n+1)*sizeof(double) );
  qarev=(double*)malloc( (n+1)*sizeof(double) );
  us=(double*)malloc( (n+1)*sizeof(double) );
  rrev=(double*)malloc( (n+1)*sizeof(double) );
  
  dt=tmax/(n);
  
  for (i=0;i<=n;i++){
    //us[i]=u(r[i],i*dt);
    us[i]=qs[n-i];
    rrev[i]=r[n-i];
  }  
  
  //qa[0]=sqrt(M_PI)*(u(r[0],0)-uo(0));
  //qa[0]=M_SQPINV*qs[n-1];
  //qa[0]=0;
  qa[0]=1.0*M_SQPINV*qs[n];
  
  for (j=1;j<=n;j++){
    if (j!=n){
      dr=(rrev[j+1]-rrev[j-1])/(2.0*dt);
    }
    else{
      dr=(rrev[j]-rrev[j-1])/dt;
    }
    qa[j]=(-0.5*us[j]-sLayAdj(j,dt,rrev,qa)-qto(j,dt)*kv(rrev[j],rrev[0],j*dt,0.0)*qa[0]+dLay(j,dt,rrev,us)+kdtntn(dr,j*dt)*mun(j,dt)*us[j])/(munAdj(j,dt)*kvtntn());
  }
  
  for (k=0;k<=n;k++){
    qarev[k]=qa[n-k];
  }
  
  return qarev;
}

/*
 * legendre recursive algorithm
 * Parameters 
 *    n    degree of legendre polynomial
 *    x    x value
 * Return
 *         final value
 */

double leg(double n, double x){
  double g;
  if (n==0){
    g=1.0;
  }
  else if(n==1){
    g=x;
  }
  else {
    g=1.0/n*((2.0*(n-1)+1.0)*x*leg(n-1,x)-(n-1.0)*leg(n-2,x));
  }
return g;
}

/*
 * derivative of legendre recursive algorithm
 * Parameters 
 *    n    degree of legendre polynomial
 *    x    x value
 * Return
 *         final value
 */

double dleg(double n, double x){
  double g;
  int i;
  
  if (n==0){
    g=0.0;
  }
  else if(n==1){
    g=1.0;
  }
  else {
    g=(2.0*(n-1.0)+1.0)*leg(n-1,x);
    for (i=3;i<n+1;i++){
      if (i%2==0){
        g=g;
      }
      else{
        g+=(2.0*((n-1.0)-(i-1.0))+1.0)*leg((n-1)-(i-1),x);
      }
    }
  }
return g;
}

/*
 * compute the boundary curve based on the Lagrange coefficients
 * Parameters 
 *    N    degree of Legendre polynomial
 *    M    number of time steps
 *    a    Legendre Polynomial Coefficients
 *    tmax Maximum Time
 * Return
 *    r    The front location as a function of time
 */

double * bdry_curve(int N, int M, double *a, double tmax){
  double *r;
  double dt, xi, ti;
  int i, j;
  
  dt=tmax/M;
  
  r=(double*)calloc(M+1, sizeof(double));
  
  for (i=0;i<M+1;i++){
    ti=i*dt;
    xi=2.0/tmax*(ti-tmax/2.0);
    for (j=0;j<N+1;j++){
      r[i]+=ti*a[j]*leg(j,xi);
    }
  }
  
  return r;
}

/*
 * returns the functional and the modified Neumann data on the fixed boundary
 *
 * Parameters
 *   M     number of time steps
 *   tMax  maximal time
 *   r     boundary curve
 *   q     modified Neumann data on the boundary curve (return)
 * Return
 *   fun   Value of the error functional
 */

double compute_functional(int M, double tMax, double *q){
  double fun,dt;
  int i;
  
  dt=tMax/(M);
  
  fun=(0.5*q[0]*q[0]+0.5*q[M]*q[M]);
  
  for (i=1;i<M;i++){
    fun+=q[i]*q[i];
  }
  
  return 0.5*dt*fun;
}

/*
 * compute the shape gradient using the adjoint method.
 *
 * Parameters:
 *   N     degree of Legendre polynomial        
 *   M     number of time steps
 *   a     legendre polynomial coefficients
 *   tMax  maximal time
 *   qS    modified Neumann data computed by compute_functional()
 *   dr    r'
 *   a2    a values defined by t*pi(z(t))
 *   
 * Return:
 *   g     gradient
 */

double * compute_gradient(int N, int M, double tMax, double *a, double *qS, double *r){
  
  double *qa,*g,*dr;
  double dt,tk,xk,tj,xj,dra,totmax,totmaxinv;
  int i,j,k,l;
  
  g=(double*)calloc(N+1, sizeof(double));
  dr=(double*)calloc(M+1, sizeof(double));
  
  qa=solveAdj(M,tMax,r,qS);
  
  dt=tMax/M;
  totmax=2.0/tMax;
  totmaxinv=1.0/totmax;
  
  for (k=0;k<M+1;k++){
    tk=k*dt;
    xk=totmax*(tk-totmaxinv);
    for (l=0;l<N+1;l++){
      dr[k]+=a[l]*leg(l,xk)+tk*totmax*a[l]*dleg(l,xk);
    }
  }
  
  for (i=0;i<N+1;i++){
    for (j=0;j<M+1;j++){
      tj=j*dt;
      xj=totmax*(tj-totmaxinv);
      if (j==0){
        //g[i]+=0.5*(qS[j]*(totmax*tj*dleg(i,xj)+leg(i,xj))-(qa[j]+qS[j]*dr[j])*tj*leg(i,xj)*(qS[j]-dr[j]));
        g[i]+=0.5*(qS[j]*(totmax*tj*dleg(i,xj)+leg(i,xj))-(qa[j]+qS[j]*dr[j])*tj*leg(i,xj)*(qS[j]-dr[j])-qa[j]*(qS[j]-dr[j])*tj*leg(i,xj));
      }
      else if(j==M){
        g[i]+=0.5*(qS[j]*(totmax*tj*dleg(i,xj)+leg(i,xj))-(qa[j]+qS[j]*dr[j])*tj*leg(i,xj)*(qS[j]-dr[j]));
      }
      else{
        //g[i]+=qS[j]*(totmax*tj*dleg(i,xj)+leg(i,xj))-(qa[j]+qS[j]*dr[j])*tj*leg(i,xj)*(qS[j]-dr[j]);
        g[i]+=qS[j]*(totmax*tj*dleg(i,xj)+leg(i,xj))-(qa[j]+qS[j]*dr[j])*tj*leg(i,xj)*(qS[j]-dr[j]);
      }
    }
    //g[i]*=dt;
    //g[i]=g[i]*dt+mun(M,dt)*qa[M]*(qS[M]-dr[M])*tj*leg(i,xj);
    g[i]=g[i]*dt+mun(M,dt)*qa[M]*(qS[M]-dr[M])*tj*leg(i,xj)*0;
  }
  return g;
}

/*
 * compute the shape gradient using the adjoint method.
 *
 * Parameters:
 *   N     degree of Legendre polynomial        
 *   M     number of time steps
 *   tMax  maximal time
 *   
 * Return:
 *   NONE
 */
 
void test_gradient(int N, int M, double tMax) {
  int l, j;
  double fcnl, fcnl1, fcnl2;
  double *q, *a, *r, *g, Gd;
  double eps = 1e-4;

  q  = (double*)calloc(M+1,sizeof(double));
  a  = (double*)calloc(N+1,sizeof(double));
  g  = (double*)calloc(N+1,sizeof(double));


  /* test curve and exact gradient */
  for (j=0; j<=N; j++ ){
    if (j==0||j==1){
      a[j]=1.0;
    }
    else{
      a[j]=0;
    }
  }
  r = bdry_curve(N, M, a, tMax);
  q = solveState(M,tMax,r);
  fcnl = compute_functional(M, tMax, q);
  g = compute_gradient(N, M, tMax, a, q, r);


  for (l=0; l<=N; l++) {
    /* first approximation using eps/2 */
    a[l] += eps;
    r=bdry_curve(N, M, a, tMax);
    a[l] -= eps;
    q=solveState(M,tMax,r);
    fcnl1 = compute_functional(M, tMax, q);
    
    a[l] -= eps;
    r=bdry_curve(N, M, a, tMax);
    a[l] += eps;
    q=solveState(M,tMax,r);
    fcnl2 = compute_functional(M, tMax, q);
    
    Gd = (fcnl1 - fcnl2)/(2.0*eps);
    printf("%d grad=%lf fd=%lf ", l, g[l], Gd);

    /* second approximation using eps/4 */
    a[l] += 0.25*eps;
    r=bdry_curve(N, M, a, tMax);
    a[l] -= 0.25*eps;
    q=solveState(M,tMax,r);
    fcnl1 = compute_functional(M, tMax, q);
    
    a[l] -= 0.25*eps;
    r=bdry_curve(N, M, a, tMax);
    a[l] += 0.25*eps;
    q=solveState(M,tMax,r);
    fcnl2 = compute_functional(M, tMax, q);
    
    Gd = 4.0*(fcnl1 - fcnl2)/(2.0*eps);
    printf("fd=%lf ", Gd);
    
    /* third approximation using eps/8 */
    a[l] += 0.125*eps;
    r=bdry_curve(N, M, a, tMax);
    a[l] -= 0.125*eps;
    q=solveState(M,tMax,r);
    fcnl1 = compute_functional(M, tMax, q);
    
    a[l] -= 0.125*eps;
    r=bdry_curve(N, M, a, tMax);
    a[l] += 0.125*eps;
    q=solveState(M,tMax,r);
    fcnl2 = compute_functional(M, tMax, q);
    
    Gd = 8.0*(fcnl1 - fcnl2)/(2.0*eps);
    printf("fd=%lf  ", Gd);
    printf("diff = %lf\n", Gd-g[l]);
  }
  exit(1);

}

/*
 * determines the Quasi-Newton-direction
 *
 * Parameters:
 *   a     edge coefficients        
 *   b     gradients
 *   x     direction to be used
 *   m     number of iterations
 *   n     number of unknowns
 *   
 * Return:
 *   y     direction determined
 */

double *quasi_newton(double *a,double *b,double *x,unsigned int m,unsigned int n){
unsigned int	i;
double		c, d, e, *y;

if (m == 0)
{  y = (double*) malloc(n*sizeof(double));
   memcpy(y,x,n*sizeof(double));
   }
else
{  c = d = 0;
   for (i=0; i<n; i++) {  
      e = a[m*n+i]-a[(m-1)*n+i];
      c += e*(b[m*n+i]-b[(m-1)*n+i]);
      d += e*x[i];
      }
   if (c == 0) 	/* Then this was just a gradient step */{
      y = quasi_newton(a,b,x,m-1,n);
      }
   else		/* Then this was a quasi-Newton step */{
      for (i=0; i<n; i++){
        x[i] -= (b[m*n+i]-b[(m-1)*n+i])*d/c;
      }
      e = 0;
      y = quasi_newton(a,b,x,m-1,n);
      for (i=0; i<n; i++){ 
        e += (b[m*n+i]-b[(m-1)*n+i])*y[i];
        
      }
      for (i=0; i<n; i++){
        y[i] += (a[m*n+i]-a[(m-1)*n+i])*(d-e)/c;
      }
      }
   }
return(y);
}

int main(int nargs, char *argv[]){
  double tmax=1.0,dt,fcnl,eps=1e-4;
  int i,j,k,x,l,numInts=512,verbose=0,deg=20,M,job=1,inner_itermax=10,outer_itermax=10;
  double *r , *qKnownList  , *qs,*qa,*grad;
  double *qTest , *uPasser , *sLays;
  double *dLays,*munValsD,*initPote,*denom,*qtList;
  double *sLaysAdj,*dLaysAdj,*munValsAdj,*initPoteAdj;
  double *denomAdj,*qAdj;
  double *qEvalList, *qtoList, *qExact;
  int     line_search_itermax=5;   /* steps in line search        */
  double  *a, *a_new;	             /* current domains             */
  double  *g;		                   /* gradient                    */
  double  *A, *G;		               /* list of domains, gradients  */
  double  df, f, f_new;	           /* values of functionals       */
  double  grad_norm;	             /* norm of gradients           */
  double  h,h_0=0.3;		               /* stepsize                    */
  double  *d;		                   /* search direction            */
  double  *q;	                     /* modified Neumann data       */
  double  error;		               /* error of q                  */
  FILE	  *file;	   	             /* protocol file               */
  time_t  t1, t2, t3;              /* time variables              */
 
  /* parse the command line */
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
      case 'v': verbose = atoi( argv[i]+3 );
        break;
      case 'T': tmax = atof( argv[i]+3 );
        break;
      case 'N': numInts = atoi( argv[i]+3 );
        break;
      case 'M': deg = atoi( argv[i]+3 );  
        break;        
      case 'j': job = atoi( argv[i]+3 );
        break;
      }

  dt=tmax/numInts;

  qa=(double*)malloc( (numInts+1)*sizeof(double) );
  r=(double*)malloc( (numInts+1)*sizeof(double) );
  grad=(double*)malloc( (numInts+1)*sizeof(double) );
  qKnownList=(double*)malloc( (numInts+1)*sizeof(double) );
  uPasser=(double*)malloc( (numInts+1)*sizeof(double) );
  //For the SLay Test
  qTest=(double*)malloc( (numInts+1)*sizeof(double) );
  //These are the qs for the second integration strategy testing module
  qEvalList=(double*)malloc( (numInts+1)*sizeof(double) );
  qExact=(double*)malloc( (numInts+1)*sizeof(double) );
  time(&t1);
  d = (double*)calloc( deg+1,sizeof(double) );
  a = (double*)calloc( deg+1,sizeof(double) );
  a_new = (double*)calloc( deg+1,sizeof(double) );
  g = (double*)calloc(deg+1,sizeof(double));
  A = (double*)calloc( inner_itermax*(deg+1),sizeof(double) );
  G = (double*)calloc( inner_itermax*(deg+1),sizeof(double) );
  q = (double*)calloc(numInts+1,sizeof(double*));  
  
  for(i=0;i<numInts+1;i++){
    qKnownList[i]=qKnown(i*dt);
    uPasser[i]=u(r[i],i*dt);
    qTest[i]=exp(i*dt);
    qEvalList[i]=qEval(i*dt)/sqrt(i*dt);
    qExact[i]=qKnown(i*dt);
  }
#if 0
// for the theoretical initial value in solve state
  for (j=0; j<=deg; j++ ){
    if (j==0){
      a[j]=0.74614;
    }
    else if(j==1){
      a[j]=-0.099985;
    }
    else if(j==2){
      a[j]=0.0426055;
    }
    else if(j==3){
      a[j]=-0.0235205;
    }
    else if(j==4){
      a[j]=0.0146;
    }
    else if(j==5){
      a[j]=-0.00960941;
    }
    else if(j==6){
      a[j]=0.00646947;
    }
    else if(j==7){
      a[j]=-0.00432055;
    }
    else if(j==8){
      a[j]=0.00273173;
    }
    else if(j==9){
      a[j]=-0.00148395;
    }
    else if(j==10){
      a[j]=0.000463357;
    }    
  }
#endif
#if 0

// for the 0 initial value in solve state  
  for (j=0; j<=deg; j++ ){
    if (j==0){
      a[j]=0.746019;
    }
    else if(j==1){
      a[j]=-0.0997427;
    }
    else if(j==2){
      a[j]=0.0422516;
    }
    else if(j==3){
      a[j]=-0.0231323;
    }
    else if(j==4){
      a[j]=0.0141652;
    }
    else if(j==5){
      a[j]=-0.00918195;
    }
    else if(j==6){
      a[j]=0.00607715;
    }
    else if(j==7){
      a[j]=-0.00397111;
    }
    else if(j==8){
      a[j]=0.00245463;
    }
    else if(j==9){
      a[j]=-0.00130293;
    }
    else if(j==10){
      a[j]=0.000412902;
    }    
  }  
  
r=bdry_curve(deg,numInts,a,tmax);
q=solveState(numInts,tmax,r);
qa=solveAdj(numInts,tmax,r,q);

for (k=0;k<numInts+1;k++){
  printf("%f\n",q[k]);
}
#endif

  if ( job==0 ){
    test_gradient(deg,numInts,tmax);
  }
  /* initial guess */
  a[0] = 1.0;
  r=bdry_curve(deg, numInts, a, tmax);
  q=solveState(numInts,tmax,r);
  f_new = f = compute_functional(numInts, tmax, q);

  file = fopen("protocol.m","w+");
  if (file == NULL) {  
    printf("ERROR: cannot open file \n");
  }
  else {  
    k = outer_itermax*inner_itermax+1;
    fprintf(file,"f = zeros(%d,1);\ng = zeros(%d,1);\na = zeros(%d,%d);\n\n", k,k,deg+1,k);
    fprintf(file,"f(1)=%g;\n",f);
    fprintf(file,"\na(:,1)=[");
    for (i=0; i<=deg; i++) fprintf(file,"\n%g",a[i]);
    fprintf(file,"\n];\n");
    fclose(file);
  }

  /* outer loop */
  for (i=0; i<outer_itermax; i++) {     
    /* inner loop */
    for (j=0; j<inner_itermax; j++) {  
      /* compute the current gradient */
      time(&t2);
      r=bdry_curve(deg, numInts, a, tmax);
      q=solveState(numInts,tmax,r);
      g=compute_gradient(deg, numInts, tmax, a, q, r);
      for ( grad_norm=0, k=0; k<=deg; k++ ) {  
        grad_norm += g[k]*g[k];
      }
      grad_norm = sqrt(grad_norm); 
      for ( error=0, k=0; k<=deg; k++ ) {  
        error += q[k]*q[k];
      }
      error = sqrt(error);
       
      
      printf("Step %d | Functional %g | Norm of Gradient %g | Error q %g\n",
             i*inner_itermax+j,f,grad_norm,error);
	 
      /* determine the search direction and initial step size */
      memcpy(&A[j*(deg+1)],a,(deg+1)*sizeof(double));
      memcpy(&G[j*(deg+1)],g,(deg+1)*sizeof(double));
      d = quasi_newton(A,G,g,j,deg+1);
      h = 0;
      for (k=0; k<=deg; k++) h += d[k]*d[k];
      if (h > 1) {  
        h = sqrt(h);
        for (k=0; k<=deg; k++) d[k] /= h;
      }
   
      /* Line search */
      printf("Line-search: ");
      h = h_0;   	/* initial step length */
      for (k=0; k<line_search_itermax; k++) {   
        
        /* boundary update */
        printf("%d->%g(%g) | ",k,h,f_new);
        for (l=0; l<=deg+1; l++) a_new[l] = a[l]-h*d[l];
        r=bdry_curve(deg, numInts, a_new, tmax);
        q=solveState(numInts,tmax,r);
        f_new = compute_functional(numInts,tmax,q);
       
        /* if necessary new step length */
        if ((f_new >= f+eps) && (h != 0))
          {  df = 0;
            for (l=0; l<=deg; l++)  df += h*g[l]*d[l];
            if (df/(df+f_new-f) > 0) h *= 0.5*df/(df+f_new-f); else h = 0;
          }
        else break;
      }

      /* print stuff */
      time(&t3);
      printf("%g sec\n",difftime(t3,t2));
      file = fopen("protocol.m","a+");
      if (file == NULL) {  
        printf("ERROR: cannot open file \n");
      }
      else { 
        fprintf(file,"g(%d)=%g;\n\n",i*inner_itermax+j+1,grad_norm);
        fprintf(file,"f(%d)=%g;\n",i*inner_itermax+j+2,f_new);
        fprintf(file,"a(:,%d)=[",i*inner_itermax+j+2);
        for (l=0; l<=deg; l++) fprintf(file,"\n%g",a_new[l]);
        fprintf(file,"\n];\n");
        fclose(file);
      }
         
      /* final boundary update */
      f = f_new;
      memcpy(a,a_new,(deg+1)*sizeof(double));
    } /* inner loop */
    printf("outer jump \n");
  } /* outer loop */

  time(&t2);	/* Endzeit */
  printf("OVERALL CPU-TIME: %g sec\n",difftime(t2,t1));

  return(0);
}