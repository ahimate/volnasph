#include "Bucket.h"

using namespace std;
double eps = 1e-20; //погрешность
double g = 9.80665;
Bucket::Bucket(double x, double y, double z) {
	bX = x;
	bY = y;
	bZ = z;
	status = 1;//пустой
}
void Bucket::BucketInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		out << part[tmp2].cX << "  ";
		out << part[tmp2].cY << "  ";
		out << part[tmp2].cZ << "  ";
		out << endl;
	}
}
void Bucket::BucketVelInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		out << part[tmp2].cU << "  ";
		out << part[tmp2].cV << "  ";
		out << part[tmp2].cW << "  ";
		out << endl;
	}
}
void Bucket::BucketRoInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		out << part[tmp2].cRo << "  ";
		out << endl;
	}
}
void Bucket::BucketPInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		out << part[tmp2].cP << "  ";
		out << endl;
	}
}

void Bucket::FluInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		if (part[tmp2].type == 1) {
			out << part[tmp2].cX << "  ";
			out << part[tmp2].cY << "  ";
			out << part[tmp2].cZ << "  ";
			out << endl;
		}
	}
}
void Bucket::FluVelInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		if (part[tmp2].type == 1) {
			out << part[tmp2].cU << "  ";
			out << part[tmp2].cV << "  ";
			out << part[tmp2].cW << "  ";
			out << endl;
		}
	}
}
void Bucket::FluRoInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		if (part[tmp2].type == 1) {
			out << part[tmp2].cRo << "  ";
			out << endl;
		}
	}
}
void Bucket::FluPInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		if (part[tmp2].type == 1) {
			out << part[tmp2].cP << "  ";
			out << endl;
		}
	}
}

double Bucket::bucRiw(Particle& a)///////////////////////////////////
{
	if (a.type == 1)//а- флюид
	{
		double riw = a.parDistanse(a, part[0]);
		int part_size = part.size();
		for (int i = 1; i < part_size; i++) {
			if (part[i].type == 0)//склон
			{
				if (a.parDistanse(a, part[i]) < riw) {
					riw = a.parDistanse(a, part[i]);
				}
			}
		}
		cout << riw << "  ";
		return riw;
	} else {
		return 0;
	}
}
double Bucket::bucDen(Particle& a, double re) {
	double ro = 0;
	if (a.type == 1) {
		if (status == 0)//там есть частицы склона
		{
			int part_size = part.size();
			for (int i = 0; i < part_size; i++) {
				if (part[i].type == 1)//поток
				{
					if ((abs(a.cX - part[i].cX) > eps) && (abs(a.cY
							- part[i].cY) > eps) && (abs(a.cZ - part[i].cZ)
							> eps)) {
						ro += a.fluRo(a, part[i], re);
					}
				} else {
					//	wallm=part[i].wallM;
				}
			}
		}
		///граница
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
void Bucket::bucForse(Particle& a, double re, double res[6])//производная частная пока
{
	int part_size = part.size();
	if ((a.type == 1) && (status == 0))//там есть частицы склона
	{
		for (int i = 0; i < part_size; i++) {
			if ((part[i].type == 1) && (abs(a.cX - part[i].cX) > eps) && (abs(
					a.cY - part[i].cY) > eps) && (abs(a.cZ - part[i].cZ) > eps)) {
				res[0] += a.fluFpress(a, part[i], re, a.cX - part[i].cX);//presx
				res[1] += a.fluFpress(a, part[i], re, a.cY - part[i].cY);//presy
				res[2] += a.fluFpress(a, part[i], re, a.cZ - part[i].cZ);//presz
				res[3] += a.fluFvis(a, part[i], re, part[i].cU - a.cU);//visx
				res[4] += a.fluFvis(a, part[i], re, part[i].cV - a.cV);//visy
				res[5] += a.fluFvis(a, part[i], re, part[i].cW - a.cW);//visz
			} else//граничное--
			{

			}
		}
	}
}

void Bucket::bucPosition(Particle& a, double timestep) {
	if (a.type == 1) {
		a.cX = a.cX + a.cU * timestep;
		a.cY = a.cY + a.cV * timestep;
		a.cZ = a.cZ + a.cW * timestep;
	}
}
