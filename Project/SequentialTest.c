/*
 * SequentialTest.c
 *
 *  Created on: 1 giu 2020
 *      Author: Fulvio
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct { int x, y, z, vx, vy, vz, mass; } Body;

void bodyForce(Body *p, int dt, int n) {

  for (int i = 0; i < n; i++) {
     int Fx = 0; int Fy = 0; int Fz = 0;

    for (int j = 0; j < n; j++) {
      int dx = p[j].x - p[i].x;
      int dy = p[j].y - p[i].y;
      int dz = p[j].z - p[i].z;
      int distSqr = dx*dx + dy*dy + dz*dz;
      int invDist = 3;//1.0f / sqrtf(distSqr);
      int invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
    }

    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
}

int main(const int argc, const char** argv) {

	  int j, nBodies = 1000;
	  int x, y, z, vx, vy, vz, mass;
	  FILE *fp;
	  time_t t;

	  const int dt = 1; // time step
	  const int nIters = 10;  // simulation iterations

	  int bytes = nBodies*sizeof(Body);
	  float *buf = (float*)malloc(bytes);
	  Body *p = (Body*)buf;
	  fp = fopen("TestValues.txt", "w");
	  srand((unsigned) time(&t));

	  for(j = 0; j < nBodies; j++){
		  x = rand() % 50;
		  y = rand() % 50;
		  z = rand() % 50;
		  vx = rand() % 50;
		  vy = rand() % 50;
		  vz = rand() % 50;
		  mass = rand() % 50;

		  fprintf(fp, "%d %d %d %d %d %d %d\n", x, y, z, vx, vy, vz, mass);

		  p[j].x = x;
		  p[j].y = y;
		  p[j].z = z;
		  p[j].vx = vx;
		  p[j].vy = vy;
		  p[j].vz = vz;
		  p[j].mass = mass;
	  }

	  fclose(fp);
	  fp = fopen("TestLog.txt", "w");

	  for (int iter = 1; iter <= nIters; iter++) {

		bodyForce(p, dt, nBodies); // compute interbody forces

		for (int i = 0 ; i < nBodies; i++) { // integrate position
		  p[i].x += p[i].vx*dt;
		  p[i].y += p[i].vy*dt;
		  p[i].z += p[i].vz*dt;
		}

		fprintf(fp, "**ITERATION %d**\n", iter);
		for(j = 0; j < nBodies; j++)
			fprintf(fp, "%d) %d %d %d %d %d %d %d\n", j, p[j].x, p[j].y, p[j].z, p[j].vx, p[j].vy, p[j].vz, p[j].mass);
	 }

	  fclose(fp);
	  free(buf);

	  return 0;
}
