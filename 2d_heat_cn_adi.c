/*
Solving 2D Heat equation using the Crank-Nicolson with Alternate Direction Implicit method
*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define xmax 100
#define ymax 100
#define itr 50000

double analytical(double, double, double);
void f_analytical(double**, double**, double, double***);
void tri(int, double**,double**,double**,double**,double**,double**);
void update_old(double***,double***);
double error(double***,double***);

double L=0.5;
double tmax=10;
double alpha=0.02;
double tmp0=100;

//analytical function
double analytical(double x, double y, double t){
	double tmp;
	tmp = exp(-2*L*L*t)*cos(L*x)*sin(L*y);
	return tmp;
}

//analytical solution
void f_analytical(double** x, double** y, double t, double*** tmp_analytical){
	int i,j;

	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){
			(*tmp_analytical)[i][j] = analytical((*x)[i],(*y)[j],t);
		}
	}
}

//solve tridiagonal system
void tri(int max, double** a, double** b, double** c, double** d, double** r, double** s){
	int i;
	for(i=1;i<max;i++){
		if(i==1){
			(*r)[1]=(*c)[1]/(*b)[1];
			(*s)[1]=(*d)[1]/(*b)[1];
		}else{
			(*r)[i]=(*c)[i]/((*b)[i]-((*a)[i]*(*r)[i-1]));
			(*r)[max-1]=0;
			(*s)[i]=((*d)[i]-((*a)[i]*(*s)[i-1]))/((*b)[i]-((*a)[i]*(*r)[i-1]));
		}
	}		
}

//update old values for next timestep
void update_old(double*** tmp_old,double*** tmp_new){
	int i,j;
	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){
			(*tmp_old)[i][j]=(*tmp_new)[i][j];
		}
	}
}

//calculate error
double error(double*** tmp_cn, double*** tmp_analytic){
	double err=0.0,sum=0.0;
	int i,j;
	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){
			sum = sum + (((*tmp_cn)[i][j] - (*tmp_analytic)[i][j]) * ((*tmp_cn)[i][j] - (*tmp_analytic)[i][j]));			
		}
	}
	err = (sqrt(sum))/(((double)xmax+1)*((double)ymax+1));
	return err;
}

int main(){
	int i,j,k;
	double **tmp_old, **tmp_new, **tmp_analytical;
	double *x, *y;
	double err;
	double zx,zy;
	double x_nd, y_nd, t_nd, dx, dy, dt, cfl;		
	double *a, *b, *c, *d, *r, *s;

	//calculating non-dimentional range of quantities	
	x_nd = L/L;
	y_nd = L/L;
	t_nd = (alpha*tmax)/(L*L);

	//calculating increments in time and space
	dx=x_nd/((double)xmax);
	dy=y_nd/((double)ymax);
	dt=t_nd/(double)itr;

	//calculating cfl numbers in x and y direction
	zx = (dt)/(dx*dx);
	zy = (dt)/(dy*dy);

	cfl=zx+zy;

	//memory allocation
	x=(double*)malloc((xmax+1)*sizeof(double));
	y=(double*)malloc((ymax+1)*sizeof(double));

	tmp_old=(double**)malloc((xmax+1)*sizeof(double*));
	for(i=0;i<=xmax;i++){
	tmp_old[i]=(double*)malloc((ymax+1)*sizeof(double));
	}
	tmp_new=(double**)malloc((xmax+1)*sizeof(double*));
	for(i=0;i<=xmax;i++){
	tmp_new[i]=(double*)malloc((ymax+1)*sizeof(double));
	}
	tmp_analytical=(double**)malloc((xmax+1)*sizeof(double*));
	for(i=0;i<=xmax;i++){
	tmp_analytical[i]=(double*)malloc((ymax+1)*sizeof(double));
	}

	//generate grid
	for(i=0;i<=xmax;i++){
	x[i] = dx*i;
	}
	for(j=0;j<=ymax;j++){
	y[j] = dy*j;
	}	

	//initialize temperature at t=0
	f_analytical(&x, &y, 0.0, &tmp_old);

	// Crank-Nicolson ADI
	//loop for incrementing time-steps	
	for(k=0;k<=itr;k++){
		//x-sweep
		//time at half point
		double tx= (k+0.5)*dt;

		//allocate memory for abcdrs arrays	
		a=(double*)malloc((xmax+1)*sizeof(double));
		b=(double*)malloc((xmax+1)*sizeof(double));
		c=(double*)malloc((xmax+1)*sizeof(double));
		d=(double*)malloc((xmax+1)*sizeof(double));
		r=(double*)malloc((xmax+1)*sizeof(double));
		s=(double*)malloc((xmax+1)*sizeof(double));

		//initialize unused locations of abcdrs arrays
		a[0]=0, b[0]=0, c[0]=0, d[0]=0, r[0]=0, s[0]=0;
		a[xmax]=0, b[xmax]=0, c[xmax]=0, d[xmax]=0, r[xmax]=0, s[xmax]=0;

		//calculate boundary points from analytical solution
		for(i=0;i<=xmax;i++){
			tmp_new[i][0] = analytical(x[i],y[0],tx);
			tmp_new[i][ymax] = analytical(x[i],y[ymax],tx);
		}
		for(j=0;j<=ymax;j++){
			tmp_new[0][j] = analytical(x[0],y[j],tx);
			tmp_new[xmax][j] = analytical(x[xmax],y[j],tx);
		}

		for(j=1;j<ymax;j++){
			//populate abcd arrays of tridiagonal system
			for(i=1;i<xmax;i++){
				a[i]= -1 * zx;
				b[i]= 2 + 2*zx;
				c[i]= -1 * zx;
			}
			a[1]=0;
			c[xmax-1]=0;

			for(i=2;i<xmax-1;i++){
				d[i] = 2*tmp_old[i][j] + zy*(tmp_old[i][j+1]-2*tmp_old[i][j]+tmp_old[i][j-1]);
			}	
			d[1] = 2*tmp_old[1][j]+ zy*(tmp_old[1][j+1]-2*tmp_old[1][j]+tmp_old[1][j-1]) + zx*tmp_new[0][j];
			d[xmax-1]= 2*tmp_old[xmax-1][j] + zy*(tmp_old[xmax-1][j+1]-2*tmp_old[xmax-1][j]+tmp_old[xmax-1][j-1]) + zx*tmp_new[xmax][j];

			//solving the tridiagonal system of eqn
			tri(xmax, &a, &b, &c, &d, &r, &s);

			//populate tmp_new
			tmp_new[xmax-1][j]=s[xmax-1];
			for(i=xmax-2;i>=1;i--){
				tmp_new[i][j]=s[i]-r[i]*tmp_new[i+1][j];
			}
		}

		//update tmp_old
		update_old(&tmp_old, &tmp_new);

		free(a);
		free(b);
		free(c);
		free(d);
		free(r);
		free(s);	

		// y sweep
		//calculate time at n+1
		double ty = tx+(0.5*dt);

		//allocate memory for abcdrs arrays
		a=(double*)malloc((ymax+1)*sizeof(double));
		b=(double*)malloc((ymax+1)*sizeof(double));
		c=(double*)malloc((ymax+1)*sizeof(double));
		d=(double*)malloc((ymax+1)*sizeof(double));
		r=(double*)malloc((ymax+1)*sizeof(double));
		s=(double*)malloc((ymax+1)*sizeof(double));

		//initialize unused locations of abcdrs arrays
		a[0]=0, b[0]=0, c[0]=0, d[0]=0, r[0]=0, s[0]=0;
		a[ymax]=0, b[ymax]=0, c[ymax]=0, d[ymax]=0, r[ymax]=0, s[ymax]=0;

		//calculate boundary points from analytical solution
		for(i=0;i<=xmax;i++){
			tmp_new[i][0] = analytical(x[i],y[0],ty);
			tmp_new[i][ymax] = analytical(x[i],y[ymax],ty);
		}
		for(j=0;j<=ymax;j++){
			tmp_new[0][j] = analytical(x[0],y[j],ty);
			tmp_new[xmax][j] = analytical(x[xmax],y[j],ty);
		}

		for(i=1;i<xmax;i++){
			//populate abcd arrays
			for(j=1;j<ymax;j++){
				a[j]= -1 * zy;
				b[j]= 2 + 2*zy;
				c[j]= -1 * zy;
			}
			a[1]=0;
			c[ymax-1]=0;

			for(j=2;j<ymax-1;j++){
				d[j]= 2*tmp_old[i][j] + zx*(tmp_old[i+1][j]-2*tmp_old[i][j]+tmp_old[i-1][j]);
			}	
			d[1]= 2*tmp_old[i][1] + zx*(tmp_old[i+1][1]-2*tmp_old[i][1]+tmp_old[i-1][1]) + zy*tmp_new[i][0];
			d[ymax-1]= 2*tmp_old[i][ymax-1] + zx*(tmp_old[i+1][ymax-1]-2*tmp_old[i][ymax-1]+tmp_old[i-1][ymax-1]) + zy*tmp_new[i][ymax];

			//solve the tridiagonal system of equations
			tri(ymax, &a, &b, &c, &d, &r, &s);

			//populate tmp_mew		
			tmp_new[i][ymax-1]=s[ymax-1];

			for(j=ymax-2;j>=1;j--){
				tmp_new[i][j]=s[j]-(r[j]*tmp_new[i][j+1]);
			}
		}

		//update tmp_old	
		update_old(&tmp_old, &tmp_new);	

		free(a);
		free(b);
		free(c);
		free(d);
		free(r);
		free(s);	

	}	

	return 0;
}
