/*
 * Nbodytest.c
 *
 *  Created on: 8 lug 2020
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
	float x, y, z, vx, vy, vz, mass;
}Particle;

void calcForce(int particle, Particle *part, int bodies){
	int prtcls;
	float fx = 0, fy = 0, fz = 0, interval = 0.1f;

	for(prtcls = 0; prtcls < bodies; prtcls++){
		if(prtcls != particle){
			float dx = part[prtcls].x - part[particle].x;
			float dy = part[prtcls].y - part[particle].y;
			float dz = part[prtcls].z - part[particle].z;
			float dist = sqrtf((dx * dx) + (dy * dy) + (dz * dz));
			float invDist3 = (1.0f / dist) * (1.0f / dist) * (1.0f / dist);
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
	int my_rank, p, nBodies, i, iterations, it, m;
	float interval = 0.1f;
	FILE *log;
	MPI_Status status;
	time_t t;
	srand((unsigned) time(&t));
	nBodies = atoi(argv[2]);
	clock_t start;
	time_t begin;

	begin = time(NULL);

	iterations = atoi(argv[1]);

	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//Manage number of particles to round them according to processes
	int roundBodies, intervalRoundBodies;

	if(nBodies % p != 0)
		roundBodies = nBodies + (p - (nBodies % p));
	else
		roundBodies = nBodies;

	intervalRoundBodies = roundBodies / p;

	int space = roundBodies * sizeof(Particle);
	float *mem = (float*)malloc(space);
	Particle *pr = (Particle*)mem;

	/*Tmp variable for the computed values in every iteration*/
	float *tmp = (float*)malloc((int)intervalRoundBodies * sizeof(Particle));
	Particle *tmpPr = (Particle*)tmp;

	/*Tmp variable for the values gathered in every iteration*/
	float *tmpRecv = (float*)malloc(roundBodies * sizeof(Particle));
	Particle *tmpRecvPr = (Particle*)tmpRecv;

	/*Create MPI datatype for struct*/
	const int nItems = 7;
	int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype types[7] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
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
		float x, y, z, vx, vy, vz, mass;
		float a = 3;
		FILE *initValues;
		initValues = fopen("InitValues.txt", "w");

		start = clock();

		//Generate random values for every particle and save it to file and struct Particle
		for(i = 0; i < roundBodies; i++){

			x = (float)rand() / (float)(RAND_MAX / a) * (rand() % 500 + 1);
			y = (float)rand() / (float)(RAND_MAX / a) * (rand() % 500 + 1);
			z = (float)rand() / (float)(RAND_MAX / a) * (rand() % 500 + 1);
			vx = (float)rand() / (float)(RAND_MAX / a) * (rand() % 500 + 1);
			vy = (float)rand() / (float)(RAND_MAX / a) * (rand() % 500 + 1);
			vz = (float)rand() / (float)(RAND_MAX / a) * (rand() % 500 + 1);
			mass = (float)rand() / (float)(RAND_MAX / a) * (rand() % 500 + 1);

			fprintf(initValues, "%d)Values: x=%f, y=%f, z=%f, vx=%f, vy=%f, vz=%f, mass=%f\n", i, x, y, z, vx, vy, vz, mass);

			pr[i].x = x;
			pr[i].y = y;
			pr[i].z = z;
			pr[i].vx = vx;
			pr[i].vy = vy;
			pr[i].vz = vz;
			pr[i].mass = mass;
		}

		fclose(initValues);
	}

	for(it = 0; it < iterations; it++){

		if(my_rank == 0){
			if(it != 0){
				int s = intervalRoundBodies;//(nBodies + p - 1) / p;

				//Copy the updated values gathered to the main struct for new broadcast
				for(m = intervalRoundBodies; m < roundBodies; m++){
					pr[m].x = tmpRecvPr[s].x;
					pr[m].y = tmpRecvPr[s].y;
					pr[m].z = tmpRecvPr[s].z;
					pr[m].vx = tmpRecvPr[s].vx;
					pr[m].vy = tmpRecvPr[s].vy;
					pr[m].vz = tmpRecvPr[s].vz;
					pr[m].mass = tmpRecvPr[s].mass;
					s++;
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Bcast(pr, roundBodies, mpi_particle, 0, MPI_COMM_WORLD);

		if(my_rank == 0){
			int j;

			printf("Number of iterations: %d\n", iterations);
			printf("Number of bodies: %d\n", nBodies);

			for(j = 0; j < intervalRoundBodies; j++){
				calcForce(j, pr, nBodies);
			}

			//Modify the position according to the new velocity
			//Done in a separate for because values must be the same for every particle
			for(j = 0; j < intervalRoundBodies; j++){
				pr[j].x += pr[j].vx * interval;
				pr[j].y += pr[j].vy * interval;
				pr[j].z += pr[j].vz * interval;
			}

			for(j = 0; j < intervalRoundBodies; j++){
				tmpPr[j].x = 0;
				tmpPr[j].y = 0;
				tmpPr[j].z = 0;
				tmpPr[j].vx = 0;
				tmpPr[j].vy = 0;
				tmpPr[j].vz = 0;
				tmpPr[j].mass = 0;
			}
		}else{
			int relBase;

			relBase = my_rank * intervalRoundBodies;
			printf("**PROCESS %d** Interval of particles: %d - %d\n", my_rank, relBase, relBase + intervalRoundBodies);

			for(relBase = my_rank * intervalRoundBodies; relBase < my_rank * intervalRoundBodies + intervalRoundBodies && relBase < nBodies; relBase++){
				calcForce(relBase, pr, nBodies);
			}

			printf("**PROCESS %d** Force calculated\n", my_rank);

			//Modify the position according to the new velocity
			//Done in a separate for because values must be the same for every particle
			for(relBase = my_rank * intervalRoundBodies; relBase <  my_rank * intervalRoundBodies + intervalRoundBodies && relBase < nBodies; relBase++){
				pr[relBase].x += pr[relBase].vx * interval;
				pr[relBase].y += pr[relBase].vy * interval;
				pr[relBase].z += pr[relBase].vz * interval;
			}

			printf("**PROCESS %d** Position modified\n", my_rank);

			//Save the updated data to tmp in order to send it in gather
			int n = 0;

			for(relBase = my_rank * intervalRoundBodies; relBase < my_rank * intervalRoundBodies + intervalRoundBodies && relBase < nBodies; relBase++){
				tmpPr[n].x = pr[relBase].x;
				tmpPr[n].y = pr[relBase].y;
				tmpPr[n].z = pr[relBase].z;
				tmpPr[n].vx = pr[relBase].vx;
				tmpPr[n].vy = pr[relBase].vy;
				tmpPr[n].vz = pr[relBase].vz;
				tmpPr[n].mass = pr[relBase].mass;
				n++;

				if(relBase == (my_rank * (nBodies + p - 1) / p) + ((nBodies + p - 1) / p) - 1)
					printf("**PROCESS %d** Transferring particle n.%d [%d]\n", my_rank, relBase, n);
			}

			printf("**PROCESS %d** Values copied in tmp\n", my_rank);

			if(my_rank == p - 1 && relBase < roundBodies){
				while(relBase < roundBodies){
					tmpPr[n].x = 0;
					tmpPr[n].y = 0;
					tmpPr[n].z = 0;
					tmpPr[n].vx = 0;
					tmpPr[n].vy = 0;
					tmpPr[n].vz = 0;
					tmpPr[n].mass = 0;
					n++;
					relBase++;
				}

				printf("**PROCESS %d** Added padding struct nodes [%d]\n", my_rank, n);
			}
		}

		FILE *tst;
		tst = fopen("test.txt", "a+");

		fprintf(tst, "It worked! Iterations: %d, bodies: %d, rank: %d\n", iterations, nBodies, my_rank);

		fclose(tst);

		printf("**PROCESS %d** Waiting for the other processes\n", my_rank);

		MPI_Barrier(MPI_COMM_WORLD);

		printf("**PROCESS %d** Starting gather operation\n", my_rank);

		MPI_Gather(tmpPr, intervalRoundBodies, mpi_particle, tmpRecvPr, intervalRoundBodies, mpi_particle, 0, MPI_COMM_WORLD);

		if(my_rank == 0){
			int l;
			log = fopen("log.txt", "a+");
			fprintf(log, "**ITERATION %d**\n", it);

			for(l = 0; l < nBodies; l++){
				fprintf(log, "%d)Values: x=%f, y=%f, z=%f, vx=%f, vy=%f, vz=%f, mass=%f\n", l, pr[l].x, pr[l].y, pr[l].z, pr[l].vx, pr[l].vy, pr[l].vz, pr[l].mass);
			}

			fclose(log);
		}
	}

	free(mem);
	free(tmpPr);
	free(tmpRecvPr);

	if(my_rank == 0){
		clock_t end = clock();
		double interval = (double) (end - start) / CLOCKS_PER_SEC;
		time_t finish = time(NULL);
		float difference = finish - begin;

		printf("Interval time: %f (clock), %f (seconds)\n", interval, difference);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Type_free(&mpi_particle);
	MPI_Finalize();

	return 0;
}
