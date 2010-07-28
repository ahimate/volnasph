#include "Particle.h"
#include "constants.h"

double koefWroo = 315/(64*pi); //1.566681480; //
double koefWpres = 14.32394496; // 45/pi
double koefWvis = 14.32394496; // 45/pi
double koefWnorm = 9.400088875; // 945/(32*pi)
double koefWsur = 9.400088875; // 945/(32*pi)

double Particle::fluM;// = 0.02; // fluid particle mass
double Particle::fluNu;// = 3.5; // viscosity coefficient
double Particle::fluROO;// = 998.29; // fluid initial density
double Particle::fluK;// = 3; // state equation coefficient (gas stiffness)
double Particle::wallM;// = 0.04; // wall particle mass
double Particle::wallROO;
double Particle::surfaceTension; // = 0.0728 // surface tension coefficient

Particle::Particle(double _x, double _y, double _z) {
	this->x = _x;
	this->y = _y;
	this->z = _z;
	cU = 0;
	cV = 0;
	cW = 0;
	cRo = wallROO;
	cP = 0;
	type = PARTICLE_TYPE_WALL;
}
Particle::Particle(double _x, double _y, double _z, double _cU, double _cV,
		double _cW) {
	this->x = _x;
	this->y = _y;
	this->z = _z;
	cU = _cU;
	cV = _cV;
	cW = _cW;
	cRo = fluROO;
	cP = 0;
	type = PARTICLE_TYPE_FLUID;
}

double Particle::getDistance(Particle& a, Particle& b) {
	return sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z
			- a.z) * (b.z - a.z));
}

double Particle::fluFpress(Particle& a, Particle& b, double re, double ba) {
	double distance = Particle::getDistance(a, b);
	if (distance > re) {
		return 0;
	} else {
		double weightCoeff = (-1) * koefWpres / (pow(re, 6));
		double Wpress = weightCoeff * pow((re - distance), 2) * ba / distance;
		double pressureCoeff = (a.cP / (pow(a.cRo, 2))) + (b.cP
				/ (pow(b.cRo, 2)));
		if (b.type == PARTICLE_TYPE_FLUID) {
			return (-1) * b.fluM * Wpress * pressureCoeff;
		} else {
			return (-1) * b.wallM * Wpress * pressureCoeff;
		}
	}
}

double Particle::fluFvis(Particle& a, Particle& b, double re, double ba) {
	double distance = Particle::getDistance(a, b);
	if ((distance > re) || (distance == 0)) {
		return 0;
	} else {
		double Wvis = koefWvis * (re - distance) / (pow(re, 6));
		if (b.type == PARTICLE_TYPE_FLUID) {
			return b.fluM * Wvis * ba / b.cRo;
		} else {
			return b.wallM * Wvis * ba / b.cRo;
		}
		return 0;
	}
}

double Particle::fluRo(Particle& a, Particle& b, double re) {
	double distance = a.getDistance(a, b);
	if (distance > re) {
		return 0;
	} else {
		double W = koefWroo * (pow(re * re - distance * distance, 3)) / pow(re,
				9);
		if (b.type == PARTICLE_TYPE_FLUID) {
			return b.fluM * W;
		} else {
			return b.wallM * W;
		}
	}
}

double Particle::fluRoWall(double riw, double wallMass, double re)// check up
{
	double weightCoeff = 315 / (64 * pi * pow(re, 9));
	double W = weightCoeff * (pow(abs(re * re - riw * riw), 3));
	return wallMass * W;
}

double Particle::fluNorm(Particle& a, Particle& b, double re, double ba) {
	double distance = a.getDistance(a, b);
	if (distance > re) {
		return 0;
	} else {
		double weightCoeff = (-1) * koefWnorm / pow(re, 9);
		double W = weightCoeff * (re * re - distance * distance) * (3 * re * re
				- 7 * distance * distance) * ba;
		if (b.type == PARTICLE_TYPE_FLUID) {
			return b.fluM * W / b.cRo;
		} else {
			return b.wallM * W / b.cRo;
		}
	}
}

double Particle::fluFsur(Particle& a, Particle& b, double re) {
	double distance = a.getDistance(a, b);
	if (distance > re) {
		return 0;
	} else {
		double weightCoeff = (-1) * koefWsur / pow(re, 9);
		double W = weightCoeff * re * (re * re - distance * distance) * (re * re - distance * distance);
		if (b.type == PARTICLE_TYPE_FLUID) {
			return b.fluM * W / b.cRo;
		} else {
			return b.wallM * W / b.cRo;
		}
	}
}
