#include "Bucket.h"

using namespace std;
double eps = 1e-20; //погрешность
double g = 9.80665;
Bucket::Bucket() {
	bucketX = bucketY = bucketZ = 0;
}
Bucket::Bucket(double bucketX, double bucketY, double bucketZ) {
	this->bucketX = bucketX;
	this->bucketY = bucketY;
	this->bucketZ = bucketZ;
	status = 1;//не пустой
}
void Bucket::BucketInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		out << part[tmp2].x << "  ";
		out << part[tmp2].y << "  ";
		out << part[tmp2].z << "  ";
		out << endl;
	}
}
void Bucket::BucketParticlesInFile(ofstream &out) {
	int tmp1 = part.size();
	out << "x" << "  ";
	out << "y" << "  ";
	out << "z" << "  ";
	out << "U" << "  ";
	out << "V" << "  ";
	out << "W" << "  ";
	out << "Ro" << "  ";
	out << "P" << "  ";
	out << endl;
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		out << part[tmp2].x << "  ";
		out << part[tmp2].y << "  ";
		out << part[tmp2].z << "  ";
		out << part[tmp2].cU << "  ";
		out << part[tmp2].cV << "  ";
		out << part[tmp2].cW << "  ";
		out << part[tmp2].cRo << "  ";
		out << part[tmp2].cP << "  ";
		out << endl;
	}
}
/*
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
 */
void Bucket::FluInFile(ofstream &out) {
	int tmp1 = part.size();
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		if (part[tmp2].type == 1) {
			out << part[tmp2].x << "  ";
			out << part[tmp2].y << "  ";
			out << part[tmp2].z << "  ";
			out << endl;
		}
	}
}
void Bucket::FluidParticlesInFile(ofstream &out) {
	int tmp1 = part.size();
	out << "x" << "  ";
	out << "y" << "  ";
	out << "z" << "  ";
	out << "U" << "  ";
	out << "V" << "  ";
	out << "W" << "  ";
	out << "Ro" << "  ";
	out << "P" << "  ";
	out << endl;
	for (int tmp2 = 0; tmp2 < tmp1; tmp2++) {
		if (part[tmp2].type == 1) {
			out << part[tmp2].x << "  ";
			out << part[tmp2].y << "  ";
			out << part[tmp2].z << "  ";
			out << part[tmp2].cU << "  ";
			out << part[tmp2].cV << "  ";
			out << part[tmp2].cW << "  ";
			out << part[tmp2].cRo << "  ";
			out << part[tmp2].cP << "  ";
			out << endl;
		}
	}
}
/*
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
*/

double Bucket::bucRiw(Particle& a)///////////////////////////////////
{
	if (a.type == 1)//а- флюид
	{
		double riw = a.getDistance(a, part[0]);
		int part_size = part.size();
		for (int i = 1; i < part_size; i++) {
			if (part[i].type == 0)//склон
			{
				if (a.getDistance(a, part[i]) < riw) {
					riw = a.getDistance(a, part[i]);
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
					if ((abs(a.x - part[i].x) > eps) && (abs(a.y - part[i].y)
							> eps) && (abs(a.z - part[i].z) > eps)) {
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
			if ((part[i].type == 1) && (abs(a.x - part[i].x) > eps) && (abs(a.y
					- part[i].y) > eps) && (abs(a.z - part[i].z) > eps)) {
				res[0] += a.fluFpress(a, part[i], re, a.x - part[i].x);//presx
				res[1] += a.fluFpress(a, part[i], re, a.y - part[i].y);//presy
				res[2] += a.fluFpress(a, part[i], re, a.z - part[i].z);//presz
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
		a.x = a.x + a.cU * timestep;
		a.y = a.y + a.cV * timestep;
		a.z = a.z + a.cW * timestep;
	}
}
