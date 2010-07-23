#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <iostream>
#include <limits>


using namespace std;

#define PARTICLE_TYPE_FLUID 1
#define PARTICLE_TYPE_WALL 0

class Particle {
public:
	double x; // x-coordinate
	double y; // y-coordinate
	double z; // z-coordinate
	double cU; // velosity along х direction
	double cV; // velosity along у direction
	double cW; // velosity along z direction
	double cRo;// density
	double cP;// pressure
	double cU5;
	double cV5;
	double cW5;

	bool type;// PARTICLE_TYPE_WALL (0) - wall ,PARTICLE_TYPE_FLUID (1) - fluid (moving particles)

	// constructor
	Particle(double _x, double _y, double _z);
	Particle(double _x, double _y, double _z, double _cU, double _cV, double cW);

	double static getDistance(Particle& a, Particle& b); // calculate distance between 2 particles

	double static fluFpress(Particle& a, Particle& b, double re, double ba); // calculate pressure force, ba=b.x(y,z)-a.x(y,z), a=i,b=j
	double static fluFvis(Particle& a, Particle& b, double re, double ba); // calculate viscosity ba=b.cU(v,w)-a.cU(v,w), a=i,b=j
	double static fluRo(Particle& a, Particle& b, double re); // calculate density: a - calculating for, b - a neighbor

	double static fluRoWall(double riw, double wallMass, double re);

	double static fluNorm(Particle& a,Particle& b, double re,double ba);
	double static fluFsur(Particle& a,Particle& b, double re);

	static double wallM; // wall particle mass
	static double fluM; // fluid particle mass
	static double fluNu; // viscosity coefficient
	static double fluROO; // fluid initial density
	static double wallROO; // wall initial density
	static double fluK; // state equation coefficient (gas stiffness)
	static double surfaceTension; // surface tension coefficient

	/*double getX();
	double getY();
	double getZ();
	double getU();
	double getV();
	double getW();
	double getRo();
	double getP();

	void setX(double _x);
	void setY(double _y);
	void setZ(double _z);
	void setU(double _cU);
	void setV(double _cV);
	void setW(double _cW);
	void setRo(double density);
	void setP(double pressure);
	*/
};
#endif
