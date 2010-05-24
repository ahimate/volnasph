#include <iostream>
#include <cmath>
using namespace std;
#ifndef PARTICLE_H
#define PARTICLE_H 1111
class Particle
{
public:
	double cX; // х-координата
	double cY; // y-координата
	double cZ; // z-координата
	double cU; // скорость частицы по х
	double cV; // скорость частицы по у
	double cW; // скорость частицы по у
	double cRo;//плотность
	double cP;//давление

	bool type;// 0-склон-мыс,1-поток-движ частицы
	// constructor
	Particle(double cpX, double cpY, double cpZ); 
	Particle(double cpX, double cpY, double cpZ, double cpU, double cpV,double cpW); 
	
	double static parDistanse(Particle& a,Particle& b); //расчет расстояния между 2 частицами

	double static fluFpress(Particle& a, Particle& b, double re,double ba); //функция давления ba=b.cX(y,z)-a.cX(y,z, a=i,b=j
	double static fluFvis(Particle& a, Particle& b, double re,double ba); //функция для вязкости ba=b.cU(v,w)-a.cU(v,w), a=i,b=j
	double static fluRoWall(double riw,double wallmass,double re);
	double static fluRo(Particle& a,Particle& b, double re); //функция для расчета плотности a-искомая,b-из ячейки
	
	static double wallM; // масса
	static double fluM; // масса
	static double fluNu; //коэф вязкости
	static double fluROO; //начальная плотность
	static double fluK; //коэф в уравн состояния
};
#endif

/*double getX(); // получить координату по x
double getY(); // получить координату по у
double getZ();// получить координату по z
double getU(); // получить скорость по x
double getV(); // получить скорость по у
double getW(); // получить скорость по z
double getRo();// получить плотность
double getP();// получить давление

void setX(double Spe); // установить координату по х
void setY(double Spe); // установить координату по у
void setZ(double Spe); // установить координату по z
void setU(double Spe); // установить скорость по х
void setV(double Spe); // установить скорость по у
void setW(double Spe); // установить скорость по z
void setRo(double Den); // установить плотность
void setP(double press); // установить давление
*/