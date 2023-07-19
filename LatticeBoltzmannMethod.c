#include <stdio.h>
#include <math.h>
#include <string.h>

#define D 16
#define LX (D*10)
#define LY (D*4)

// D2Q9
#define F 9
#define cs_sq (1./3.)
#define wC (4./9.)
#define wS (1./9.)
#define wL (1./36.)

// LES turbulence model
#define C_LES (.19)

enum LatVel {_C,_N,_S,_W,_E,_NE,_NW,_SE,_SW};
enum BC {BC_IN, BC_OUT, FLUID, SOLID};

double cx[F] = {0,0, 0,-1,1,1,-1, 1,-1};
double cy[F] = {0,1,-1, 0,0,1, 1,-1,-1};
double w[F] = {wC,wS,wS,wS,wS,wL,wL,wL,wL};

// distribution functions
double cells[LX*LY*F];
double tmp_cells[LX*LY*F];
#define AR(x,y,f) (f+F*(x)+F*LX*(y))

// SOLID/FLUID cell type
int map[LX*LY];
double Rho[LX*LY];
double Ux[LX*LY];
double Uy[LX*LY];
#define ARM(x,y) ((x)+LX*(y))

double tau;
double uinLB, vinLB;

double vel_max = 0;

double feq(int i, double rho, double u, double v) {
  double cu = (cx[i]*u + cy[i]*v)/cs_sq;
  return w[i]*rho*(1 + cu + 0.5*cu*cu - 0.5*(u*u+v*v)/cs_sq);
}

void collide_stream(double *c, double *tc) {
  double rho,ux,uy;
  double fstar[F];
  int i,x,y,xp,xm,ym,yp;
  double fc,fe,fn,fw,fs,fne,fnw,fsw,fse;
  double k3,k4,k5,k6,k7,k8,kxxyy_at;

  vel_max = 0;

  for (y = 0; y < LY; y++)
    for (x = 0; x < LX; x++) {
      int ctype = map[ARM(x,y)];
      if (ctype != SOLID) {
      /* FLUID CELLS*/
        if (ctype == FLUID) {
          /* INSIDE OF THE COMP DOMAIN */
          rho = 0; ux = 0; uy = 0;
          // MAKRO
          for (i = 0; i < F; i++) {
            rho += c[AR(x,y,i)];
		        ux += cx[i]*c[AR(x,y,i)];
		        uy += cy[i]*c[AR(x,y,i)];
	        }
	        ux /= rho; uy /= rho;
	        Rho[ARM(x,y)] = rho;
	        Ux[ARM(x,y)] = ux;
	        Uy[ARM(x,y)] = uy;
          vel_max = sqrt(ux*ux + uy*uy) > vel_max ? sqrt(ux*ux + uy*uy) : vel_max;

          fc = c[AR(x,y,_C)];
          fe = c[AR(x,y,_E)];
          fn = c[AR(x,y,_N)];
          fw = c[AR(x,y,_W)];
          fs = c[AR(x,y,_S)];
          fne = c[AR(x,y,_NE)];
          fnw = c[AR(x,y,_NW)];
          fsw = c[AR(x,y,_SW)];
          fse = c[AR(x,y,_SE)];
		  
		  k3 = 1./12.*(rho*(ux*ux + uy*uy) - fe - fn -fs - fw 
			   -2.*(fse + fsw + fne + fnw - rho/3.));
		  k4 = .25/tau*(fn + fs - fe - fw + rho*(ux*ux - uy*uy));
		  k5 = .25/tau*(fne + fsw - fnw - fse - ux*uy*rho);

		  kxxyy_at = (6*k3 + 2*k4 - rho*ux*ux + fe + fw + fne + fnw + fse + fsw) *
					 (6*k3 - 2*k4 - rho*uy*uy + fn + fs + fne + fnw + fse + fsw);
					/* rho/9 */

		  k6 = -((fse + fsw - fne - fnw - 2*ux*ux*uy*rho + uy*(rho - fn - fs - fc))*.25
				+ .5*ux*(fne - fnw - fse + fsw) - .5*uy*(-3*k3 - k4) - 2*ux*k5);

		  k7 = -((fsw + fnw - fse - fne - 2*uy*uy*ux*rho + ux*(rho - fw - fe - fc))*.25
				+ .5*uy*(fne + fsw - fse - fnw) - .5*ux*(-3*k3 + k4) - 2*uy*k5);

		  k8 = .25*(kxxyy_at - fne - fnw - fse - fsw + 
				2*(ux*(fne - fnw + fse - fsw) + uy*(fne + fnw - fse -fsw)) +
				4*ux*uy*(fnw - fne + fse - fsw) -
				ux*ux*(fn + fne + fnw + fs + fse + fsw) +
				uy*uy*(3*ux*ux*rho - fe - fne - fnw - fse - fsw - fw)) -
				2*k3 - 2*ux*k7 - 2*uy*k6 + 4*ux*uy*k5 +
				(-1.5*k3 + .5*k4)*ux*ux +
				(-1.5*k3 - .5*k4)*uy*uy;
	/*
    NP = .25f * (rho * (kxxyy) - Ne - Nw - Se - Sw - 8 * P
    + 2 * (v.x * (Ne - Nw + Se - Sw - 4 * RIGHT) + v.y * (Ne + Nw - Se - Sw - 4 * UP))
    + 4 * v.x * v.y * (-Ne + Nw + Se - Sw + 4 * V)
    + v.x * v.x * (-N - Ne - Nw - S - Se - Sw + 2 * NE - 6 * P)
    + v.y * v.y * ((-E - Ne - Nw - Se - Sw - W - 2 * NE - 6 * P) + 3 * v.x * v.x * rho));
	*/
			fstar[_C] = fc + 4*(-k3 + k8);
	    	fstar[_W] = fw - k3 - 2*k8 + k4 - 2*k7;
		    fstar[_E] = fe - k3 - 2*k8 + k4 + 2*k7;
		    fstar[_N] = fn - k3 - 2*k8 - k4 + 2*k6;
		    fstar[_S] = fs - k3 - 2*k8 - k4 - 2*k6;
		    fstar[_NW] = fnw + 2*k3 + k8 + k5 - k6 + k7;
		    fstar[_NE] = fne + 2*k3 + k8 - k5 - k6 - k7;
		    fstar[_SW] = fsw + 2*k3 + k8 - k5 + k6 + k7;
		    fstar[_SE] = fse + 2*k3 + k8 + k5 + k6 - k7;


		  			
          /*
          for (i = 0; i < F; i++) {
            if (fstar[i]<0)
              fstar[i] = feq(i,rho, ux, uy);
          }
          */
        }
        if (ctype == BC_OUT) { //do-nothing BC ~ f_i(x-1,y,t) copy
          /* OUTLET */
          /*c[AR(x,y,_NW)] = c[AR(x-1,y,_NW)];
          c[AR(x,y,_W)] = c[AR(x-1,y,_W)];
          c[AR(x,y,_SW)] = c[AR(x-1,y,_SW)];
          rho = 0; ux = 0; uy = 0;
          // MAKRO
          for (i = 0; i < F; i++) {
            rho += c[AR(x,y,i)];
		        ux += cx[i]*c[AR(x,y,i)];
		        uy += cy[i]*c[AR(x,y,i)];
	        }
          ux /= rho; uy /= rho;
          Rho[ARM(x,y)] = rho;
          Ux[ARM(x,y)] = ux;
          Uy[ARM(x,y)] = uy;
          // COLLISION
          for (i = 0; i < F; i++)
       	    fstar[i] = (1-1/tau)*c[AR(x,y,i)] + feq(i, rho, ux, uy)/tau;
          */
          rho = 1.;
          double uoutLB = (Ux[ARM(x-1,y)]+Ux[ARM(x-1,y+1)]+Ux[ARM(x-1,y-1)])/3.;
          double voutLB = 0;
          for (i = 0; i < F; i++)
       	    fstar[i] = feq(i, rho, uoutLB, voutLB);
          Rho[ARM(x,y)] = rho;
          Ux[ARM(x,y)] = uoutLB;
          Uy[ARM(x,y)] = voutLB;
        }
        if (ctype == BC_IN) { // velocity inlet (rho=1) ~ f_i^eq(rho,u_inlet)
        /* INLET */
          //rho = (Rho[ARM(x+1,y)]+Rho[ARM(x+1,y+1)]+Rho[ARM(x+1,y-1)])/3.;
          rho = 1.0;
          for (i = 0; i < F; i++)
       	    fstar[i] = feq(i, rho, uinLB, vinLB);

	        Rho[ARM(x,y)] = rho;
	        Ux[ARM(x,y)] = uinLB;
	        Uy[ARM(x,y)] = vinLB;
        }
      } else {
        /* SOLID CELLS */
		    //fstar[_C] = c[AR(x,y,_C)];
		    fstar[_W] = c[AR(x,y,_E)];
		    fstar[_E] = c[AR(x,y,_W)];
		    fstar[_N] = c[AR(x,y,_S)];
		    fstar[_S] = c[AR(x,y,_N)];
		    fstar[_NW] = c[AR(x,y,_SE)];
		    fstar[_NE] = c[AR(x,y,_SW)];
		    fstar[_SW] = c[AR(x,y,_NE)];
		    fstar[_SE] = c[AR(x,y,_NW)];
      }
      // PERIODIC BC
      xp = ( x == LX - 1 ) ? 0 : x + 1;
      yp = ( y == LY - 1 ) ? 0 : y + 1;
      xm = ( x == 0 ) ? LX - 1 : x - 1;
      ym = ( y == 0 ) ? LY - 1 : y - 1;

      // STREAMING
      tc[AR(x,y,_C)] = fstar[_C];
      tc[AR(x,yp,_N)] = fstar[_N];
      tc[AR(x,ym,_S)] = fstar[_S];
      tc[AR(xm,y,_W)] = fstar[_W];
      tc[AR(xp,y,_E)] = fstar[_E];
      tc[AR(xp,yp,_NE)] = fstar[_NE];
      tc[AR(xp,ym,_SE)] = fstar[_SE];
      tc[AR(xm,yp,_NW)] = fstar[_NW];
      tc[AR(xm,ym,_SW)] = fstar[_SW];
    }
}

void setIC(double rho0, double u0, double v0) {
  int x,y,i;

  for (y = 0; y < LY; y++)
    for (x = 0; x < LX; x++)
      for (i = 0; i < F; i++) {
        cells[AR(x,y,i)] = feq(i, rho0, u0, v0);
        Rho[ARM(x,y)] = rho0;
        Ux[ARM(x,y)] = u0;
        Uy[ARM(x,y)] = v0;
        map[ARM(x,y)] = FLUID;
        if (x == 0 && y > 3*D) map[ARM(x,y)] = BC_IN;
        if (x == LX-1 && (y > 3*D || y < D)) map[ARM(x,y)] = BC_OUT;
        if (y == 0 || y == LY-1) {
          map[ARM(x,y)]=SOLID;
          //Rho[ARM(x,y)]=0;
          //Ux[ARM(x,y)] = 0;
          //Uy[ARM(x,y)] = 0;
        }
	    }

  // SCIANA PIONOWA 0
  for ( y = 3*D - 1; y > 0; y--)
    map[ARM(0,y)] = SOLID;

  // SCIANA PIONOWA 1
  for ( y = 3*D; y > 2*D; y--)
    map[ARM(D,y)] = SOLID;

  // SCIANA PIONOWA 2
  for ( y = 4*D - 1; y > 2*D; y--)
    map[ARM(2*D,y)] = SOLID;

  // SCIANA PIONOWA 3
  for ( y = 3*D; y > 0; y--)
    map[ARM(4*D,y)] = SOLID;

  // SCIANA PIONOWA 4
  for ( y = 4*D - 1; y > D; y--)
    map[ARM(6*D,y)] = SOLID;

  // SCIANA PIONOWA 5
  for ( y = 2*D; y > 0; y--)
    map[ARM(9*D,y)] = SOLID;

  // SCIANA PIONOWA 6
  for ( y = 3*D; y > D; y--)
    map[ARM(10*D-1,y)] = SOLID;

  // SCIANA POZIOMA 0
  for ( x = 1; x < D; x++)
    map[ARM(x,3*D)] = SOLID;

  // SCIANA POZIOMA 1
  for ( x = 8*D; x < 10*D-1; x++)
    map[ARM(x,3*D)] = SOLID;

  // SCIANA POZIOMA 2
  for ( x = 8*D; x < 9*D; x++)
    map[ARM(x,2*D)] = SOLID;
}
// ParaView
// VisIt

void dumpStateVTK(char *fname) {
 FILE *fp;
 int x,y;

 fp = fopen(fname, "w");
 fprintf(fp,"# vtk DataFile Version 2.0\n");
 fprintf(fp,"2D-ADE data file \n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET RECTILINEAR_GRID\n");
 fprintf(fp,"DIMENSIONS %d %d %d\n", LX, LY, 1);
 fprintf(fp,"X_COORDINATES %d int\n", LX);
 for (x = 0; x < LX; x++)
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Y_COORDINATES %d int\n", LY);
 for (x = 0; x < LY; x++)
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Z_COORDINATES 1 int\n");
 fprintf(fp, "0\n");
 fprintf(fp,"POINT_DATA %d \n", LX*LY);
 fprintf(fp,"SCALARS density double 1\n");
 fprintf(fp,"LOOKUP_TABLE default\n");
 for (y = 0; y < LY; y++)
   for (x = 0; x < LX; x++)
    fprintf(fp, "%e \n", Rho[ARM(x,y)]); 
 fprintf(fp,"VECTORS velocity double\n");
 for (y = 0; y < LY; y++)
   for (x = 0; x < LX; x++)
    fprintf(fp, "%e %e 0.0\n", Ux[ARM(x,y)], Uy[ARM(x,y)]); 
 fprintf(fp,"\n");
 fclose(fp);
}

int main() {
int iter = 0;
int ITERMAX = 1e5;
double rho0LB,u0LB,v0LB,nuLB;
double Re, Ma, N;
char fname[256];

  Ma = 0.1;		// Ma = U/c_s
  Re = 10;	// Re = UL/nu
  N = LX;

  // IC setup
  rho0LB = 1.;
  u0LB = Ma * sqrt(cs_sq);
  v0LB = 0;

  // BC setup
  uinLB = u0LB;
  vinLB = 0;

  setIC(rho0LB, u0LB, v0LB);
  dumpStateVTK("state0.vtk");
 
  // PHYS setup
  nuLB = u0LB * N / Re;

  // LBM setup
  tau = nuLB/cs_sq + 0.5;
  printf(" nu_LB: %lf\n tau: %lf\n u0LB: %lf\n", nuLB, tau, u0LB);
  do {
    if (iter%2)
      collide_stream(tmp_cells, cells);
    else
      collide_stream(cells, tmp_cells);
    iter++;
	if ( !(iter%100) ) {
    printf("# vel_max = %g\n", vel_max);
    sprintf(fname,"state%d.vtk",iter/100);
    dumpStateVTK(fname);
  }
  } while ( iter < ITERMAX);
   //dumpStateVTK("state.vtk");
   return 0;
}