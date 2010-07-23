#ifndef BUCKET_H
#define BUCKET_H

#include "Particle.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#define THERE_IS_FLUID 1
#define THERE_NO_FLUID 0

class Bucket {
public:
	vector<Particle> part;
	bool status; // THERE_IS_FLUID or THERE_NO_FLUID
	double bucketX; // beginning х-coordinate
	double bucketY; // beginning y-coordinate
	double bucketZ; // beginning z-coordinate
	// constructor
	Bucket();
	Bucket(double bucketX, double bucketY, double bucketZ);

	void BucketInFile(ofstream &out);
	void BucketParticlesInFile (ofstream &out);

	void FluInFile(ofstream &out);
	void FluidParticlesInFile (ofstream &out);

	double bucDen(Particle& a, double re);
	void bucForse(Particle& a, double re, double res[6]);
	void bucFnorm(Particle& a,double re,double res[3]);
	double bucFsurf(Particle& a,double re);

	double bucRiw(Particle& a);
	void bucPosition(Particle& a, double timestep);
};

#endif
