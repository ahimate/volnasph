#include "Bucket.h"
#include "constants.h"

Bucket::Bucket() {
	bucketX = bucketY = bucketZ = 0;
	status = THERE_NO_FLUID;
}
Bucket::Bucket(double bucketX, double bucketY, double bucketZ) {
	this->bucketX = bucketX;
	this->bucketY = bucketY;
	this->bucketZ = bucketZ;
	status = THERE_NO_FLUID;//there are no fluid particles
}
void Bucket::BucketParticlesInFile(ofstream &out) {
	int numberOfParticles = particlesVector.size();
	for (int i = 0; i < numberOfParticles; i++) {
		out << particlesVector[i].x << "\t";
		out << particlesVector[i].y << "\t";
		out << particlesVector[i].z << "\t";
		out << particlesVector[i].cU << "\t";
		out << particlesVector[i].cV << "\t";
		out << particlesVector[i].cW << "\t";
		out << particlesVector[i].cRo << "\t";
		out << particlesVector[i].cP << "\t";
		out << endl;
	}
}

void Bucket::FluidParticlesInFile(ofstream &out) {
	int numberOfParticles = particlesVector.size();
	for (int i = 0; i < numberOfParticles; i++) {
		if (particlesVector[i].type == PARTICLE_TYPE_FLUID) {
			out << particlesVector[i].x << "\t";
			out << particlesVector[i].y << "\t";
			out << particlesVector[i].z << "\t";
			out << particlesVector[i].cU << "\t";
			out << particlesVector[i].cV << "\t";
			out << particlesVector[i].cW << "\t";
			out << particlesVector[i].cRo << "\t";
			out << particlesVector[i].cP << "\t";
			out << endl;
		}
	}
}

double Bucket::bucRiw(Particle& a)//
{
	if (a.type == PARTICLE_TYPE_FLUID)//if а - fluid
	{
		double riw = a.getDistance(a, particlesVector[0]);
		int numberOfParticles = particlesVector.size();
		for (int i = 1; i < numberOfParticles; i++) {
			if (particlesVector[i].type == PARTICLE_TYPE_WALL)//slope
			{
				if (a.getDistance(a, particlesVector[i]) < riw) {
					riw = a.getDistance(a, particlesVector[i]);
				}
			}
		}
		//		cout << "Riw=" << riw << "  ";
		return riw;
	} else {
		return 0;
	}
}
double Bucket::bucDen(Particle& a, double re) {
	double ro = 0;
	//if (a.type == PARTICLE_TYPE_FLUID) { //change density for all particles, not only for fluid
	//		if (status == THERE_IS_FLUID) // there are fluid particles // all particles influence on density
	//		{
	int numberOfParticles = particlesVector.size();
	for (int i = 0; i < numberOfParticles; i++) {
		//	if (particlesVector[i].type == PARTICLE_TYPE_FLUID)//a fluid flow // this is being checking in fluRo() function
		//	{
		//if (a.getDistance(a, particlesVector[i]) > epsilon) { //do not need this check for density
		ro += a.fluRo(a, particlesVector[i], re);
		//}
		//	} else {
		//		wallm = particlesVector[i].wallM; //boundary by
		//	}
	}
	//		}
	//}
	return ro;
}

void Bucket::bucForse(Particle& a, double re, double result[6])//производная частная пока
{
	int numberOfParticles = particlesVector.size();
	//	if ((a.type == PARTICLE_TYPE_FLUID) && (status == THERE_IS_FLUID))//там есть частицы склона
	//	{
	for (int i = 0; i < numberOfParticles; i++) {
		//		if ((particlesVector[i].type == PARTICLE_TYPE_FLUID)
		//			&& (a.getDistance(a, particlesVector[i]) > epsilon)) {
		if ((a.getDistance(a, particlesVector[i]) > epsilon) && (a.getDistance(
				a, particlesVector[i]) < re)) {
			result[0] += a.fluFpress(a, particlesVector[i], re, a.x
					- particlesVector[i].x);//presx
			result[1] += a.fluFpress(a, particlesVector[i], re, a.y
					- particlesVector[i].y);//presy
			result[2] += a.fluFpress(a, particlesVector[i], re, a.z
					- particlesVector[i].z);//presz
			result[3] += a.fluFvis(a, particlesVector[i], re,
					particlesVector[i].cU - a.cU);//visx
			result[4] += a.fluFvis(a, particlesVector[i], re,
					particlesVector[i].cV - a.cV);//visy
			result[5] += a.fluFvis(a, particlesVector[i], re,
					particlesVector[i].cW - a.cW);//visz
		} else // collision or the same particle
		{

		}
	}
	//	}
}

void Bucket::bucFnorm(Particle& a, double re, double result[3]) {
	int numberOfParticles = particlesVector.size();
	//	if (status == THERE_IS_FLUID) {
	for (int i = 0; i < numberOfParticles; i++) {
		result[0] += a.fluNorm(a, particlesVector[i], re, a.x
				- particlesVector[i].x);
		result[1] += a.fluNorm(a, particlesVector[i], re, a.y
				- particlesVector[i].y);
		result[2] += a.fluNorm(a, particlesVector[i], re, a.z
				- particlesVector[i].z);
	}
	//	}
}

double Bucket::bucFsurf(Particle& a, double re) {
	int numberOfParticles = particlesVector.size();
	double result = 0;
	//	if (status == THERE_IS_FLUID) {
	for (int i = 0; i < numberOfParticles; i++) {
		result += a.fluFsur(a, particlesVector[i], re);
		//		}
	}
	return result;
}

void Bucket::bucPosition(Particle& a, double timestep) {
	if (a.type == PARTICLE_TYPE_FLUID) {
		a.x = a.x + a.cU * timestep;
		a.y = a.y + a.cV * timestep;
		a.z = a.z + a.cW * timestep;
	}
}
