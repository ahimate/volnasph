#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <iostream>
#include <limits>
using namespace std;

class Particle {
public:
	double x; // х-координата
	double y; // y-координата
	double z; // z-координата
	double cU; // скорость частицы по х
	double cV; // скорость частицы по у
	double cW; // скорость частицы по у
	double cRo;//плотность
	double cP;//давление

	bool type;// 0-склон-мыс,1-поток-движ частицы
	// constructor
	Particle(double x, double y, double z);
	Particle(double x, double y, double z, double u, double v, double w);

	double static getDistance(Particle& a, Particle& b); //расчет расстояния между 2 частицами

	double static fluFpress(Particle& a, Particle& b, double re, double ba); //функция давления ba=b.x(y,z)-a.x(y,z, a=i,b=j
	double static fluFvis(Particle& a, Particle& b, double re, double ba); //функция для вязкости ba=b.cU(v,w)-a.cU(v,w), a=i,b=j
	double static fluRoWall(double riw, double wallMass, double re);
	double static fluRo(Particle& a, Particle& b, double re); //функция для расчета плотности a-искомая,b-из ячейки

	static double wallM; // масса
	static double fluM; // масса
	static double fluNu; //коэф вязкости
	static double fluROO; //начальная плотность
	static double fluK; //коэф в уравн состояния
};
#endif
