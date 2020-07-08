/*
 * ParallelTest.c
 *
 *  Created on: 1 giu 2020
 *      Author: Fulvio
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"

typedef struct{
	int x, y, z, vx, vy, vz, mass;
}Particle;

void calcForce(int particle, Particle *part, int bodies){
	int prtcls;
	int fx = 0, fy = 0, fz = 0, interval = 1;

	for(prtcls = 0; prtcls < bodies; prtcls++){
		if(prtcls != particle){
			int dx = part[prtcls].x - part[particle].x;
			int dy = part[prtcls].y - part[particle].y;
			int dz = part[prtcls].z - part[particle].z;
			int dist = dx * dx + dy * dy + dz * dz;
			int invDist = 3;
			int invDist3 = invDist * invDist * invDist;
			fx += dx * invDist3;
			fy += dy * invDist3;
			fz += dz * invDist3;
		}
	}

	part[particle].vx += interval * fx;
	part[particle].vy += interval * fy;
	part[particle].vz += interval * fz;
}

int main(int argc, char *argv[]){
	int my_rank, p, nBodies, i, tag = 0, iterations;
	int interval = 1;
	FILE *log;
	MPI_Status status;
	time_t t;
	srand((unsigned) time(&t));
	nBodies = 1000;
	int space = nBodies * sizeof(Particle);
	float *mem = (float*)malloc(space);
	Particle *pr = (Particle*)mem;
	clock_t start;

	iterations = 10;

	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	/*Create MPI datatype for struct*/
	const int nItems = 7;
	int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype types[7] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype mpi_particle;
	MPI_Aint     offsets[7];

	offsets[0] = offsetof(Particle, x);
	offsets[1] = offsetof(Particle, y);
	offsets[2] = offsetof(Particle, z);
	offsets[3] = offsetof(Particle, vx);
	offsets[4] = offsetof(Particle, vy);
	offsets[5] = offsetof(Particle, vz);
	offsets[6] = offsetof(Particle, mass);

	MPI_Type_create_struct(nItems, blocklengths, offsets, types, &mpi_particle);
	MPI_Type_commit(&mpi_particle);

	if(my_rank == 0){
		int j, values;
		FILE *initValues;
		initValues = fopen("TestValues.txt", "r");

		start = clock();

		for(i = 0; i < nBodies; i++){
			fscanf(initValues, "%d", &pr[i].x);
			fscanf(initValues, "%d", &pr[i].y);
			fscanf(initValues, "%d", &pr[i].z);
			fscanf(initValues, "%d", &pr[i].vx);
			fscanf(initValues, "%d", &pr[i].vy);
			fscanf(initValues, "%d", &pr[i].vz);
			fscanf(initValues, "%d", &pr[i].mass);
		}

		fclose(initValues);

		printf("Number of iterations: %d\n", iterations);
		printf("Number of bodies: %d\n", nBodies);

		//Send the initial values to the Workers
		for(j = 1; j < p; j++){
			for(values = 0; values < nBodies; values++)
				MPI_Send(&pr[values], 1, mpi_particle, j, tag, MPI_COMM_WORLD);
		}

		//Calculate the force applied on the interval of particles given and change velocity
		for(i = 0; i < iterations; i++){
			for(j = 0; j < (nBodies + p - 1) / p; j++){
				calcForce(j, pr, nBodies);
			}

			//Modify the position according to the new velocity
			//Done in a separate for because values must be the same for every particle
			for(j = 0; j < (nBodies + p - 1) / p; j++){
				pr[j].x += pr[j].vx * interval;
				pr[j].y += pr[j].vy * interval;
				pr[j].z += pr[j].vz * interval;
			}

			//Receive all the new values
			for(j = 1; j < p; j++){
				for(values = j * (nBodies + p - 1) / p; values < (j * ((nBodies + p - 1) / p)) + ((nBodies + p - 1) / p) && values < nBodies; values++)
					MPI_Recv(&pr[values], 1, mpi_particle, j, tag, MPI_COMM_WORLD, &status);
			}

			log = fopen("log.txt", "a+");
			fprintf(log, "**ITERATION %d**\n", i);

			for(j = 0; j < nBodies; j++){
				fprintf(log, "%d) %d %d %d %d %d %d %d\n", j, pr[j].x, pr[j].y, pr[j].z, pr[j].vx, pr[j].vy, pr[j].vz, pr[j].mass);
			}

			//Send the updated particles to the Workers
			for(j = 1; j < p; j++){
				for(values = 0; values < nBodies; values++)
					MPI_Send(&pr[values], 1, mpi_particle, j, tag, MPI_COMM_WORLD);
			}

			fclose(log);
		}

		clock_t end = clock();
		double interval = (double) (end - start) / CLOCKS_PER_SEC;

		printf("Interval time: %f\n", interval);
	}else{
		int relBase, k, l;
		FILE *tst;
		tst = fopen("test.txt", "a+");

		for(k = 0; k < nBodies; k++)
			MPI_Recv(&pr[k], 1, mpi_particle, 0, tag, MPI_COMM_WORLD, &status);

		relBase = my_rank * ((nBodies + p - 1) / p);
		printf("**PROCESS %d** Interval of particles: %d - %d\n", my_rank, relBase, relBase + ((nBodies + p - 1) / p));

		for(k = 0; k < iterations; k++){
			for(relBase = my_rank * ((nBodies + p - 1) / p); relBase <  my_rank * ((nBodies + p - 1) / p) + ((nBodies + p - 1) / p) || relBase < nBodies; relBase++){
				calcForce(relBase, pr, nBodies);
			}

			//Modify the position according to the new velocity
			//Done in a separate for because values must be the same for every particle
			for(relBase = my_rank * ((nBodies + p - 1) / p); relBase <  my_rank * ((nBodies + p - 1) / p) + ((nBodies + p - 1) / p); relBase++){
				pr[relBase].x += pr[relBase].vx * interval;
			    pr[relBase].y += pr[relBase].vy * interval;
			    pr[relBase].z += pr[relBase].vz * interval;
			}

			fprintf(tst, "**ITERATION %d**\n", k);

			for(l = 500; l < nBodies; l++){
				fprintf(tst, "%d) %d %d %d %d %d %d %d\n", l, pr[l].x, pr[l].y, pr[l].z, pr[l].vx, pr[l].vy, pr[l].vz, pr[l].mass);
			}

			//Send all the new values
			for(relBase = my_rank * (nBodies + p - 1) / p; relBase < (my_rank * (nBodies + p - 1) / p) + ((nBodies + p - 1) / p) && relBase < nBodies; relBase++)
				MPI_Send(&pr[relBase], 1, mpi_particle, 0, tag, MPI_COMM_WORLD);

			//Receive the updated particles from the Master
			for(relBase = 0; relBase < nBodies; relBase++)
				MPI_Recv(&pr[relBase], 1, mpi_particle, 0, tag, MPI_COMM_WORLD, &status);
		}

		fprintf(tst, "It worked! Iterations: %d, bodies: %d, rank: %d\n", iterations, nBodies, my_rank);

		fclose(tst);
	}

	free(pr);

	if(my_rank == 0){
		clock_t end = clock();
		double interval = (double) (end - start) / CLOCKS_PER_SEC;

		printf("Interval time: %f\n", interval);
	}

	MPI_Type_free(&mpi_particle);
	MPI_Finalize();

	return 0;
}
