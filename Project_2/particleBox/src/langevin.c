#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

// Modules


////////////////
#define PI 3.14159265359f

#define N 1000  // Number of Particles
#define L 100.0 // Box Length

#define SIGMA 1.0 // Particle size
#define EPSILON 1.0 // WCA Adjustment

int main(){

    float clock_start, clock_end;
    clock_start = clock();
	
	int t=0,tmax=10000;

	int i, j;
	
	double avg_x=0, avg_x2=0;
	double avg_y=0, avg_y2=0;
    double averagePressure=0.0; // Average pressure over the last 10 steps
    int pressWindow = 100; // Number of steps to average pressure over
	
	double r;
	
	double dt=0.01;   /* time step size, energy, drag, mass */
	double kT=1.0;
	double g=1.0;
	// double m=1.0;

    double Fx=0.0, Fy=0.0;
    double r_cutoff = pow(2.0, 1.0/6.0)*SIGMA;
    double dx,dy,F_ij;
	
	
	double x[N],vx[N],vhx[N],ax[N];
	double y[N],vy[N],vhy[N],ay[N];
    double m[N];
	
	double ran,ran1,ran2,ran3;
	double gr1,gr2;
	double pref1,pref2;
    int IS_placed,IS_touching;
    double try_x,try_y, theta;

    // Seed random
    srand(time(NULL));
	
	pref1=g;
	pref2=sqrt(2*kT*g/dt);
	
	FILE *outt,*outstats,*outpressure;
	outt=fopen("langevin2D_traj.xyz","w");
    outpressure=fopen("langevin2D_pressure.dat","w");
	outstats=fopen("langevin2D_stats.txt","w");

	/* initialize the pasticle positions and velocities */
	avg_x=0;
	avg_x2=0;
	// for(i=0;i<N;i++){
	// 	x[i]=0;
	// 	y[i]=0;
	// 	ran1=(float)rand()/RAND_MAX;
	// 	ran2=(float)rand()/RAND_MAX;
	// 	gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*PI*ran2);
	// 	gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*PI*ran2);
	// 	vx[i]=gr1;
	// 	vy[i]=gr2;
	// 	ran1=(float)rand()/RAND_MAX;
	// 	ran2=(float)rand()/RAND_MAX;
	// 	gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*PI*ran2);
	// 	vy[i]=gr1;		
	// 	ax[i]=0;
	// 	ay[i]=0;
	// }
    printf("Placing Particles... \n");
    for (i = 0; i < N; i++) {
        IS_placed = 0;
        m[i] = 1.0;
        while (!IS_placed) {
            // Attempt to place particle at random position
            try_x = (float)rand() / (float)RAND_MAX * L - L/2;
            try_y = (float)rand() / (float)RAND_MAX * L - L/2;

            // Test placement attempt for overlap
            IS_touching = 0;
            for (j = 0; j < i; j++) {
                dx = try_x - x[j];
                dy = try_y - y[j];
                r = sqrt(dx * dx + dy * dy);
                
                // If particles are touching, break and try again
                if (r < r_cutoff) {
                    IS_touching = 1;
                    break;
                }
            }
            // If no overlap, place particle
            if (!IS_touching) {
                x[i] = try_x;
                y[i] = try_y;
                IS_placed = 1;
            }
        }

        // Randomize velocity
        theta = (float)rand() / (float)RAND_MAX * 2 * PI;
        vx[i] = sqrt(2 * kT / m[i]) * cos(theta);
        vy[i] = sqrt(2 * kT / m[i]) * sin(theta);

        m[i] = 1.0;

        // Clean Arrays
        ax[i] = 0;
        ay[i] = 0;
        vhx[i] = 0;
        vhy[i] = 0;
    }
		
	fprintf(outstats,"%E %E %E %E\n",avg_x,avg_x2,avg_y,avg_y2);

    // Main Loop //
    printf("Starting Simulation... \n");
	while(t<tmax){
		avg_x=0;
		avg_x2=0;
		avg_y=0;
		avg_y2=0;
		
		if(t%10==0){
            // Write to traj header
			fprintf(outt,"%i\n",N);
			fprintf(outt,"title\n");
		}

        if(t%pressWindow == 0){
            // Write pressure data
            fprintf(outpressure,"%d %E\n",t,averagePressure/(t * 4.0 * L)); // Pressure = Force / Area, Area = 4*L for 2D box with 4 walls
            // averagePressure = 0.0; // Reset pressure accumulator after writing
        }
		/* run the velocity-Verlet algorithm */
		for(i=0;i<N;i++){
			if(t%10==0){
				// Write to traj
                fprintf(outt,"a%i %f %f 0.0\n",i,x[i],y[i]);
            }
            // Half Step Velocity Update
			vhx[i]=vx[i]+0.5*ax[i]*dt;
			vhy[i]=vy[i]+0.5*ay[i]*dt;
			
            // Position Update
			x[i]=x[i]+vhx[i]*dt;
			y[i]=y[i]+vhy[i]*dt;
        }

        // Interaction Loop
        for (i = 0; i < N; i++) {

            // Calculate net force on particle [i]
            Fx = 0.0;
            Fy = 0.0;
            // Get force from all the other particles
            for (j = 0; j < N; j++) {
                if (i != j) {
                    // Calculate distance between particles
                    dx = x[j] - x[i];
                    dy = y[j] - y[i];
                    r = sqrt(dx * dx + dy * dy);

                    if (r < r_cutoff) { // If the particles are within interaction distance (2^(1/6) * sigma)

                        // Calculate force
                        F_ij = -(4.0 * EPSILON) * (1.0/r) * (6.0 * pow(SIGMA/r, 6.0) - 12.0 * pow(SIGMA/r, 12.0));
                        // F_ij = 0.0;

                        if (F_ij > 1000) {printf("Force over limit at t = %d! F = %f r = %f\n", t, F_ij, r);} // Flag large forces
                        if (F_ij < 0) {printf("Negative Force at t = %f! F = %f r = %f\n", t, F_ij, r);} // Flag negative (attractive) forces

                        Fx -= F_ij * dx / r; // x-component of force
                        Fy -= F_ij * dy / r; // y-component of force
                    }
                }
            }



            // Acceleration Update
			ax[i]=0;
			ay[i]=0;
				
			ran1=(double)rand()/RAND_MAX-0.5;
			ran2=(double)rand()/RAND_MAX-0.5;
			
			ax[i]+=-pref1*vhx[i]+pref2*ran1 + Fx;
			ay[i]+=-pref1*vhy[i]+pref2*ran2 + Fy;

            ax[i] /= m[i];
            ay[i] /= m[i];
			
            // Full Step Velocity Update
			vx[i]=vhx[i]+0.5*ax[i]*dt;
			vy[i]=vhy[i]+0.5*ay[i]*dt;

                        // Apply Reflecting Boundary Conditions
            /* When a particle crosses the boundary, its new position is the wall position minus the distance it overshot */
            /* x_reflected = wall_position - (x[i] - wall_position) = wall_position*2 - x[i] */

            if(x[i]<-L/2.0){ // Left Wall
                x[i]= -L - x[i];
                vx[i]=-vx[i];
                averagePressure += 2*m[i]*fabs(vhx[i])/dt;
            }
            if(x[i]>L/2.0){  // Right Wall
                x[i]= L - x[i];
                vx[i]=-vx[i];
                averagePressure += 2*m[i]*fabs(vhx[i])/dt;
            }
            if(y[i]<-L/2.0){ // Bottom Wall
                y[i]= -L - y[i];
                vy[i]=-vy[i];
                averagePressure += 2*m[i]*fabs(vhy[i])/dt;
            }
            if(y[i]>L/2.0){  // Top Wall
                y[i]= L - y[i];
                vy[i]=-vy[i];
                averagePressure += 2*m[i]*fabs(vhy[i])/dt;
            }
            ////////////////////////////////////////


			// Record Data
			avg_x+=x[i];
			avg_x2+=x[i]*x[i];
			avg_y+=y[i];
			avg_y2+=y[i]*y[i];
		}
		
		avg_x=avg_x/(double)N;
		avg_x2=avg_x2/(double)N;
		avg_y=avg_y/(double)N;
		avg_y2=avg_y2/(double)N;
		
		fprintf(outstats,"%E %E %E %E\n",avg_x,avg_x2,avg_y,avg_y2);
		
		t++;
	}

    clock_end = clock();
    double elapsed_time = (clock_end - clock_start) / CLOCKS_PER_SEC;
    printf("Elapsed Time: %f seconds\n", elapsed_time);
}
			
			