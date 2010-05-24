#include <iostream>
#include <cmath>
using namespace std;
#ifndef PARTICLE_H
#define PARTICLE_H 1111
class Particle
{
public:
	double cX; // �-����������
	double cY; // y-����������
	double cZ; // z-����������
	double cU; // �������� ������� �� �
	double cV; // �������� ������� �� �
	double cW; // �������� ������� �� �
	double cRo;//���������
	double cP;//��������

	bool type;// 0-�����-���,1-�����-���� �������
	// constructor
	Particle(double cpX, double cpY, double cpZ); 
	Particle(double cpX, double cpY, double cpZ, double cpU, double cpV,double cpW); 
	
	double static parDistanse(Particle& a,Particle& b); //������ ���������� ����� 2 ���������

	double static fluFpress(Particle& a, Particle& b, double re,double ba); //������� �������� ba=b.cX(y,z)-a.cX(y,z, a=i,b=j
	double static fluFvis(Particle& a, Particle& b, double re,double ba); //������� ��� �������� ba=b.cU(v,w)-a.cU(v,w), a=i,b=j
	double static fluRoWall(double riw,double wallmass,double re);
	double static fluRo(Particle& a,Particle& b, double re); //������� ��� ������� ��������� a-�������,b-�� ������
	
	static double wallM; // �����
	static double fluM; // �����
	static double fluNu; //���� ��������
	static double fluROO; //��������� ���������
	static double fluK; //���� � ����� ���������
};
#endif

/*double getX(); // �������� ���������� �� x
double getY(); // �������� ���������� �� �
double getZ();// �������� ���������� �� z
double getU(); // �������� �������� �� x
double getV(); // �������� �������� �� �
double getW(); // �������� �������� �� z
double getRo();// �������� ���������
double getP();// �������� ��������

void setX(double Spe); // ���������� ���������� �� �
void setY(double Spe); // ���������� ���������� �� �
void setZ(double Spe); // ���������� ���������� �� z
void setU(double Spe); // ���������� �������� �� �
void setV(double Spe); // ���������� �������� �� �
void setW(double Spe); // ���������� �������� �� z
void setRo(double Den); // ���������� ���������
void setP(double press); // ���������� ��������
*/