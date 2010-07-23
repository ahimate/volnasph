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
	status = THERE_IS_FLUID;//there are fluid particles
}
void Bucket::BucketInFile(ofstream &out) {
	int numberOfParticles = part.size();
	for (int i = 0; i < numberOfParticles; i++) {
		out << part[i].x << "  ";
		out << part[i].y << "  ";
		out << part[i].z << "  ";
		out << endl;
	}
}
void Bucket::BucketParticlesInFile(ofstream &out) {
	int numberOfParticles = part.size();
	out << "x" << "  ";
	out << "y" << "  ";
	out << "z" << "  ";
	out << "U" << "  ";
	out << "V" << "  ";
	out << "W" << "  ";
	out << "Ro" << "  ";
	out << "P" << "  ";
	out << endl;
	for (int i = 0; i < numberOfParticles; i++) {
		out << part[i].x << "  ";
		out << part[i].y << "  ";
		out << part[i].z << "  ";
		out << part[i].cU << "  ";
		out << part[i].cV << "  ";
		out << part[i].cW << "  ";
		out << part[i].cRo << "  ";
		out << part[i].cP << "  ";
		out << endl;
	}
}

void Bucket::FluInFile(ofstream &out) {
	int numberOfParticles = part.size();
	for (int i = 0; i < numberOfParticles; i++) {
		if (part[i].type == PARTICLE_TYPE_FLUID) {
			out << part[i].x << "  ";
			out << part[i].y << "  ";
			out << part[i].z << "  ";
			out << endl;
		}
	}
}
void Bucket::FluidParticlesInFile(ofstream &out) {
	int numberOfParticles = part.size();
	out << "x" << "  ";
	out << "y" << "  ";
	out << "z" << "  ";
	out << "U" << "  ";
	out << "V" << "  ";
	out << "W" << "  ";
	out << "Ro" << "  ";
	out << "P" << "  ";
	out << endl;
	for (int i = 0; i < numberOfParticles; i++) {
		if (part[i].type == PARTICLE_TYPE_FLUID) {
			out << part[i].x << "  ";
			out << part[i].y << "  ";
			out << part[i].z << "  ";
			out << part[i].cU << "  ";
			out << part[i].cV << "  ";
			out << part[i].cW << "  ";
			out << part[i].cRo << "  ";
			out << part[i].cP << "  ";
			out << endl;
		}
	}
}

double Bucket::bucRiw(Particle& a)//
{
	if (a.type == PARTICLE_TYPE_FLUID)//if а - fluid
	{
		double riw = a.getDistance(a, part[0]);
		int numberOfParticles = part.size();
		for (int i = 1; i < numberOfParticles; i++) {
			if (part[i].type == PARTICLE_TYPE_WALL)//slope
			{
				if (a.getDistance(a, part[i]) < riw) {
					riw = a.getDistance(a, part[i]);
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
	if (a.type == PARTICLE_TYPE_FLUID) {
		if (status == THERE_IS_FLUID) // there are fluid particles
		{
			int numberOfParticles = part.size();
			for (int i = 0; i < numberOfParticles; i++) {
//				if (part[i].type == PARTICLE_TYPE_FLUID)//a fluid flow
//				{
					if (a.getDistance(a, part[i]) > epsilon) {
						ro += a.fluRo(a, part[i], re);
					}
//				} else {
//					wallm=part[i].wallM; //boundary by Harada
//				}
			}
		}
		//boundary by Harada
		/*
		 if(wallm!=0)
		 {
		 ro+=a.fluRoWall(bucRiw(a),wallm,re);
		 }
		 a.cRo=ro;
		 */
	}
	return ro;
}

void Bucket::bucForse(Particle& a, double re, double result[6])//производная частная пока
{
	int numberOfParticles = part.size();
	if ((a.type == PARTICLE_TYPE_FLUID) && (status == THERE_IS_FLUID))//там есть частицы склона
	{
		for (int i = 0; i < numberOfParticles; i++) {
			if ((part[i].type == PARTICLE_TYPE_FLUID) && (a.getDistance(a, part[i]) > epsilon)) {
				result[0] += a.fluFpress(a, part[i], re, a.x - part[i].x);//presx
				result[1] += a.fluFpress(a, part[i], re, a.y - part[i].y);//presy
				result[2] += a.fluFpress(a, part[i], re, a.z - part[i].z);//presz
				result[3] += a.fluFvis(a, part[i], re, part[i].cU - a.cU);//visx
				result[4] += a.fluFvis(a, part[i], re, part[i].cV - a.cV);//visy
				result[5] += a.fluFvis(a, part[i], re, part[i].cW - a.cW);//visz
			} else // boundary condition (?)
			{

			}
		}
	}
}

void Bucket::bucFnorm(Particle& a,double re,double result[3])
{
	int numberOfParticles=part.size();
//	if(status==0)
//	{
		for(int i=0;i<numberOfParticles;i++)
		{
			result[0]+=a.fluNorm(a,part[i],re,a.x-part[i].x);
			result[1]+=a.fluNorm(a,part[i],re,a.y-part[i].y);
			result[2]+=a.fluNorm(a,part[i],re,a.z-part[i].z);
		}
//	}
}

double Bucket::bucFsurf(Particle& a,double re)
{
	int numberOfParticles=part.size();
	double result=0;
	if(status==THERE_IS_FLUID)
	{
		for(int i=0;i<numberOfParticles;i++)
		{
			result+=a.fluFsur(a,part[i],re);
		}
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
