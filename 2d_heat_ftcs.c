/*
Solving 2D Heat equation using the Forward Time Central Space explicit method and the Crank-Nicolson with Alternate Direction Implicit method
*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double analytical(double, double, double);
void f_analytical(double**, double**, double, double***);
void ftcs_diffusion_2d(double**,double**,double,double,double,double,double***,double***);
void update_old(double***,double***);
double error(double***,double***);

double L=0.5;
double tmax=10;
double alpha=0.02;
double tmp0=100;
int xmax=100;
int ymax=100;
int itr=50000;

//function for calculating ftcs solution
void ftcs_diffusion_2d(double** x, double** y, double t, double dx, double dy, double dt, double*** tmp_old, double*** tmp_new){
	int i,j;
	double zx, zy;

	//calculate z1 and z2
	zx = dt/(dx*dx);
	zy = dt/(dy*dy);

	//ftcs to find tmp_new
	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){

			//update boundary points using analytical solution
			if(i==0 || i==xmax || j==0 || j==ymax){
				(*tmp_new)[i][j]=analytical((*x)[i],(*y)[j],t);
			}

			else{	
				(*tmp_new)[i][j] = (*tmp_old)[i][j] + zx*((*tmp_old)[i+1][j]-2*(*tmp_old)[i][j]+(*tmp_old)[i-1][j]) + zy*((*tmp_old)[i][j+1]-2*(*tmp_old)[i][j]+(*tmp_old)[i][j-1]);
			}
		}
	}

}

//update old for next timestep
void update_old(double*** tmp_old,double*** tmp_new){
	int i,j;

	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){
			(*tmp_old)[i][j]=(*tmp_new)[i][j];
		}
	}
}

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

//calculate error
double error(double*** tmp_ftcs, double*** tmp_analytic){
	double err=0.0,sum=0.0;
	int i,j;
	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){
			sum += ((*tmp_ftcs)[i][j]- (*tmp_analytic)[i][j]) * ((*tmp_ftcs)[i][j]- (*tmp_analytic)[i][j]);			
		}
	}
	err = (sqrt(sum))/((xmax+1)*(ymax+1));
	return err;
}

int main(){
	int i,j,k;
	double x_nd, y_nd, t_nd;
	double dx,dy,dt,t,cfl;
	double **tmp_old, **tmp_new, **tmp_analytical;
	double *x, *y;
	double err;

	//finding non-dimensional equivalents
	x_nd = L/L;
	y_nd = L/L;
	t_nd = (alpha*tmax)/(L*L);

	//calculating increments in space and time
	dx=x_nd/xmax;
	dy=y_nd/ymax;
	dt=t_nd/itr;

	cfl=(dt/(dx*dx))+(dt/(dy*dy));

	//check for stability
	if(cfl>0.5){
		printf("check cfl\n");
		exit(0);
	}

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

	//initialize using the analytical equation
	f_analytical(&x, &y, 0.0, &tmp_old);

	//ftcs calculations
	for(k=1;k<=itr;k++){
		t=dt*k;
		ftcs_diffusion_2d(&x,&y,t,dx,dy,dt,&tmp_old,&tmp_new);
		update_old(&tmp_old,&tmp_new);
	}

	//free memory
	for(i=0;i<=xmax;i++){
		free(tmp_old[i]);
	}
	free(tmp_old);
	for(i=0;i<=xmax;i++){
		free(tmp_new[i]);
	}
	free(tmp_new);
	for(i=0;i<=xmax;i++){
		free(tmp_analytical[i]);
	}
	free(tmp_analytical);
	free(x);
	free(y);

	return 0;
}



