//попытка дописать 29 января
#include <iostream>
#include <cmath>
#include <vector>
#include "Particle.h"
#include "Bucket.h"
#include <string>
#include <fstream> 
#include <map>
using namespace std;

double grdXmin,grdXmax,grdYmin,grdYmax,grdZmin,grdZmax;
int parXcount, parYcount,parZcount; //число точек по x,y,z
double fluXStep, fluYStep, fluZStep; //шаги по x,y,z
int wallXcount, wallYcount;
double wallXmin,wallXmax,wallYmin,wallYmax,wallZmin,wallZmax;

double epsil=0.00001; //погрешность
int gridXstep, gridYstep;
double** wall;//под батиметрию-высоты частиц склона/////////////
double** flu;//////////////////////////////////////
int filename_count=5;
double timestep=0.01;
int kolvotimestep=1;

///Bucket
double bucRadius; // радиус ячейки
map <int, Bucket> mapa;
map <int, Bucket> mapa_new;
//const int bucCountMax=3000;//ограничение на количество ячеек
int kolvoParN;//колво частиц всего
int kolvoParInRadX=40;

/////////flu
double fluROOmain=998.29; //начальная плотность
double fluMmain=0.06;//масса частички
int fluPmain=101325;/////////////////////////////////

/////////wall
double wallRoo=1000;//начальная плотность

int kodir(double x,double y,double z)//c нуля,а количество с 1
{
	if((x<wallXmin)||(x>(wallXmax+bucRadius))||(y<wallYmin)||(y>(wallYmax+bucRadius)))
	{
		cout<<"error kodir: point outside the area("<<x<<", "<<y<<", "<<z<<" )"<<endl;
		int tmp1;
		cin>>tmp1;
		return 0;
	}
	else
	{
		int xx=(int)((x-wallXmin)/(2*bucRadius));
		int yy=(int)((y-wallYmin)/(2*bucRadius));
		int zz=(int)((z-wallZmin)/(2*bucRadius));
		if((xx<999)&&(yy<999)&&(z<999))
		{
			return xx*1000000+yy*1000+zz;
		}
		else 
		{
			cout<<"error kodir: count bucket is big"<<endl;
			int tmp1;
			cin>>tmp1;
			return 0;
		}
	}
}


void outputData(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		(*it).second.BucketInFile(out);
	}
	out.close();
	cout<<"create "<<outname<<endl;
}

void outputDataVel(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		(*it).second.BucketVelInFile(out);
	}
	out.close();
	cout<<"create "<<outname<<endl;
}
void outputDataRo(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		(*it).second.BucketRoInFile(out);
	}
	out.close();
	cout<<"create "<<outname<<endl;
}
void outputDataP(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		(*it).second.BucketPInFile(out);
	}
	out.close();
	cout<<"create "<<outname<<endl;
}

void outputFluData(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		if((*it).second.status==0)
		{
			(*it).second.FluInFile(out);
		}
	}
	out.close();
	cout<<"create "<<outname<<endl;
}

void outputFluDataVel(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		if((*it).second.status==0)
		{
			(*it).second.FluVelInFile(out);
		}
	}
	out.close();
	cout<<"create "<<outname<<endl;
}
void outputFluDataRo(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		if((*it).second.status==0)
		{
			(*it).second.FluRoInFile(out);
		}
	}
	out.close();
	cout<<"create "<<outname<<endl;
}
void outputFluDataP(map<int, Bucket> mapama,char* outname)
{
	ofstream out(outname);
	for(map<int, Bucket>::iterator it = mapama.begin(); it!=mapama.end();it++)
	{
		if((*it).second.status==0)
		{
			(*it).second.FluPInFile(out);
		}
	}
	out.close();
	cout<<"create "<<outname<<endl;
}


void inputData(char* fluidSurf,char* wallSurf,char* outnm) // ввод начальных данных из файла
{
	cout<<endl<<"INPUT DATA "<<fluidSurf<<" , "<<wallSurf<<endl;
	ofstream out(outnm);
	ifstream influ(fluidSurf);
	ifstream inwall(wallSurf);

	if(influ.fail())
	{
		cout<<"Error reading "<<fluidSurf<<"'."<<endl;
		int tmp1;
		cin>>tmp1;
	}
	if(inwall.fail())
	{
		cout<<"Error reading "<<wallSurf<<"'."<<endl;
		int tmp1;
		cin>>tmp1;
	}
	wall = new double*[gridYstep];
	for(int i=0;i<gridYstep;i++)
		wall[i]	= new double [gridXstep];
	flu = new double*[gridYstep];
	for(int i=0;i<gridYstep;i++)
		flu[i] = new double [gridXstep];

	string tmp2;
	inwall>>tmp2>>wallXcount>>wallYcount>>wallXmin>>wallXmax>>wallYmin>>wallYmax>>wallZmin>>wallZmax;
	influ>>tmp2>>parXcount>>parYcount>>grdXmin>>grdXmax>>grdYmin>>grdYmax>>grdZmin>>grdZmax;
	
	if(wallXcount!=parXcount)
	{
		cout<<"Error size X in file"<<endl;
		int tmp1;
		cin>>tmp1;
	}
	if(wallYcount!=parYcount)
	{
		cout<<"Error size Y in file"<<endl;
		int tmp1;
		cin>>tmp1;
	}

	fluXStep=(grdXmax-grdXmin)/(parXcount-1);
	fluYStep=(grdYmax-grdYmin)/(parYcount-1);
	fluZStep=fluXStep;
	
	double tmp3=0;
	for(int j=0;j<parYcount;j++)
	{
		for(int i=0;i<parXcount;i++)
		{
			influ>>flu[i][j];
			inwall>>wall[i][j];
			tmp3+=flu[i][j];
		}
	}
	double Vall=tmp3*fluXStep*fluYStep;
	kolvoParN=(int)((fluROOmain*Vall)/fluMmain);
	double Vpar=Vall/kolvoParN;
	double parStep=pow(Vpar,0.333333333);
	tmp3=(3*Vall*kolvoParInRadX)/(4*3.1415926358*kolvoParN);
	bucRadius=pow(tmp3,0.33333333333);
	out<<"kolvoParN   "<<kolvoParN<<endl;
	out<<"parStep   "<<parStep<<endl;
	out<<"bucRadius   "<<bucRadius<<endl;
	if(bucRadius<parStep)
	{
		cout<<"Error kolvoParN and input data"<<endl;
		int tmp1;
		cin>>tmp1;
	}
	
	//склон разбиваем на частицы
	for(int i=0;i<parXcount;i++)
	{
		for(int j=0;j<parYcount;j++)
		{
			double parx=i*fluXStep+wallXmin;
			double pary=j*fluYStep+wallYmin;
			double parz=wall[i][j];
			
			//bool tmp4=0;
			double tmp3=parz-bucRadius;
			for(double k=tmp3;k<parz;k+=parStep)
			{
				for(double tmpx=0;tmpx<fluXStep;tmpx+=parStep)
				{
					for(double tmpy=0;tmpy<fluYStep;tmpy+=parStep)
					{
						Particle a(parx+tmpx,pary+tmpy,k);
						/*if(k==parz)
						{
							tmp4=1;
						}*/
						if(((parx+tmpx)<wallXmin)||((parx+tmpx)>(wallXmax+bucRadius))||((pary+tmpy)<wallYmin)||((pary+tmpy)>(wallYmax+bucRadius)))
						{
							//послед раз добавленные частицы-лишние отметаются
						}
						else
						{
							if(mapa.find(kodir(parx+tmpx,pary+tmpy,k)) == mapa.end())
							{
								Bucket b(int((parx+tmpx-wallXmin)/(2*bucRadius))*(2*bucRadius)+wallXmin,
							     int((pary+tmpy-wallYmin)/(2*bucRadius))*(2*bucRadius)+wallYmin, 
								 int((k-wallZmin)/(2*bucRadius))*(2*bucRadius)+wallZmin);
								mapa[kodir(parx+tmpx,pary+tmpy,k)]=b;
							}	
							mapa[kodir(parx+tmpx,pary+tmpy,k)].part.push_back(a);
							mapa[kodir(parx+tmpx,pary+tmpy,k)].status=1;
						}
					}
				}
			}
		/*	if(tmp4==0)
			{
				double k=parz;
				Particle a(parx,pary,k);
			}
		*/
		}
	}

	//поток разбиваем на частицы
	for(int j=0;j<parYcount;j++)
	{
		for(int i=0;i<parXcount;i++)
		{
			double parz=flu[i][j]+wall[i][j];
			double parx=i*fluXStep+grdXmin;
			double pary=j*fluYStep+grdYmin;

			for(double k=wall[i][j];k<=parz;k+=parStep)
			{
				for(double tmpx=0;tmpx<fluXStep;tmpx+=parStep)
				{
					for(double tmpy=0;tmpy<fluYStep;tmpy+=parStep)
					{
						Particle a(parx+tmpx,pary+tmpy,k,0,0,0);
						/*if(k==parz)
						{
						tmp4=1;
						}*/
						if(((parx+tmpx)<wallXmin)||((parx+tmpx)>(wallXmax+bucRadius))||((pary+tmpy)<wallYmin)||((pary+tmpy)>(wallYmax+bucRadius)))
						{
							//послед раз добавленные частицы-лишние отметаются
						}
						else
						{
							if(mapa.find(kodir(parx+tmpx,pary+tmpy,k)) == mapa.end())
							{
								Bucket b(int((parx+tmpx-wallXmin)/(2*bucRadius))*(2*bucRadius)+wallXmin,
									int((pary+tmpy-wallYmin)/(2*bucRadius))*(2*bucRadius)+wallYmin, 
									int((k-wallZmin)/(2*bucRadius))*(2*bucRadius)+wallZmin);
								mapa[kodir(parx+tmpx,pary+tmpy,k)]=b;
							}	
							mapa[kodir(parx+tmpx,pary+tmpy,k)].part.push_back(a);
							mapa[kodir(parx+tmpx,pary+tmpy,k)].status=0;
						}
					}
				}
			}
		}
	}
	cout<<"Finish INPUT DATA "<<endl;
}

/*
//бежит по мапа_нью и в пустой мапа добавляет ячейки и частицы
void updateBucket()
{
	mapa.clear();
	int count=0;
	for(map<int, Bucket>::iterator it = mapa_new.begin(); it!=mapa_new.end();it++)
	{
		int part_size=(*it).second.part.size();
		for(int i=0;i<part_size;i++)
		{
			double tmpx=(*it).second.part[i].cX;
			double tmpy=(*it).second.part[i].cY;
			double tmpz=(*it).second.part[i].cZ;
			
			if((tmpx<wallXmin)||(tmpx>wallXmax)||(tmpy<wallYmin)||(tmpy>wallYmax))
			{
				//точку вне области не считаем
				count++;
			}
			else
			{
				if(mapa.find(kodir(tmpx,tmpy,tmpz)) == mapa.end())
				{
					Bucket a(int((tmpx-grdXmin)/bucDiameter)*bucDiameter+grdXmin,
					int((tmpy-grdYmin)/bucDiameter)*bucDiameter+grdYmin, 
					int((tmpz-grdZmin)/bucDiameter)*bucDiameter+grdZmin);
					mapa[kodir(tmpx,tmpy,tmpz)]=a;
				}	
				mapa[kodir(tmpx,tmpy,tmpz)].part.push_back((*it).second.part[i]);
				mapa[kodir(tmpx,tmpy,tmpz)].status=0;
			}
		}
	}
	cout<<"Free particles "<<count<<endl;
}
*/
void BucketStep(double t)
{
	//создаем массив на новом временном шаге
	mapa_new=mapa;
	//расчет плотности
	for(map<int, Bucket>::iterator it = mapa_new.begin(); it!=mapa_new.end();it++)
	{
		int part_size=(*it).second.part.size();
		for(int i=0;i<part_size;i++)
		{
			//double tmp1=(*it).second.bucDen((*it).second.part[i],(bucDiameter/2));
			double tmp1=0;
			for(int x1=-1;x1<2;x1++)
			{
				for(int x2=-1;x2<2;x2++)
				{
					for(int x3=-1;x3<2;x3++)
					{
						if(mapa_new.find((*it).first+x1*1000000+x2*1000+x3)==mapa_new.end())
						{
							//соседа нету
						}
						else
						{
						tmp1+=mapa[((*it).first+x1*1000000+x2*1000+x3)].bucDen((*it).second.part[i],bucRadius);
						}
					}
				}
			}
			(*it).second.part[i].cRo=tmp1;
			cout<<tmp1<<" ";
			if(tmp1<=epsil)
			{
				cout<<"Ro==0: "<<i;
				int tmp11;
				cin>>tmp11;
			}
			//cout<<tmp1<<"  ";
			if((*it).second.part[i].type==1)
			{
				(*it).second.part[i].cP=fluPmain+(*it).second.part[i].fluK*((*it).second.part[i].cRo-(*it).second.part[i].fluROO);
			//	cout<<(*it).second.part[i].cP<<endl;
			}
		}		
	}
	cout<<"Density and press";
	//расчет сил-проверить все
	for(map<int, Bucket>::iterator it = mapa_new.begin(); it!=mapa_new.end();it++)
	{
		int part_size=(*it).second.part.size();
		for(int i=0;i<part_size;i++)
		{
			double forse[6]={0,0,0,0,0,0};
			//////////
			for(int x1=-1;x1<2;x1++)
			{
				for(int x2=-1;x2<2;x2++)
				{
					for(int x3=-1;x3<2;x3++)
					{
						if(mapa_new.find((*it).first+x1*1000000+x2*1000+x3)==mapa_new.end())
						{
							//соседа нету
						}
						else
						{
							double f[6]={0,0,0,0,0,0};
							mapa_new[((*it).first+x1*1000000+x2*1000+x3)].bucForse((*it).second.part[i],bucRadius,f);
							for(int k=0;k<6;k++){forse[k]+=f[k];}
						}
					}
				}
			}
			double acceleration[3]={0,0,0};
			for(int k=0;k<3;k++)
			{
				double tmp8=(*it).second.part[i].cRo*forse[k]+(*it).second.part[i].fluNu*forse[k+3];
				if(k==2)
				{
					tmp8+=(*it).second.part[i].cRo*9.80665;//g
				}
				acceleration[k]=tmp8/((*it).second.part[i].cRo);
			//	cout<<acceleration[k]<<"  ";
			}
			//cout<<endl;
		}	
	}
	//		(*it).second.bucVel((*it).second.part[i],bucDiameter,timestep);
	//		(*it).second.bucPosition((*it).second.part[i],timestep );
	cout<<"Computation "<<endl;
}
void inputParams(char* input)
{
	ifstream in(input);
	if(in.fail())
	{
		cout<<"Error reading "<<input<<"'."<<endl;
		int tmp1;
		cin>>tmp1;
	}
	in>>Particle::fluM>>Particle::fluNu>>Particle::fluROO>>Particle::fluK>>Particle::wallM>>gridXstep>>gridYstep;
}
int main(int argc, char** argv)
{
	if(argc!=4)
	{
		cout<<"Error argv"<<endl;
		return 1;
	}
	inputParams(argv[3]);
	inputData(argv[1], argv[2],"Data.txt");//	inputData("f.grd","w.grd");
//	inputParams("input.txt");
//	inputData("fluid.grd","wall.grd","Data.txt");//	inputData("f.grd","w.grd");

//	outputData(mapa,"000.txt");
//	outputDataVel(mapa,"vel000.txt");
//	outputDataRo(mapa,"ro000.txt");
	
	outputFluData(mapa,"Flu_000.txt");
//	outputFluDataVel(mapa,"Flu_vel000.txt");
//	outputFluDataRo(mapa,"Flu_ro000.txt");
	for(int t=1;t<=kolvotimestep;t++)
	{
		cout<<"t   "<<t<<endl;
		BucketStep(t*timestep);
	//	updateBucket();

//		outputData(mapa_new,"001.txt");
//		outputDataVel(mapa_new,"vel001.txt");
//		outputDataRo(mapa_new,"ro001.txt");
//		outputDataP(mapa_new,"p001.txt");
		
//		outputFluData(mapa_new,"Flu_001.txt");
//		outputFluDataVel(mapa_new,"Flu_vel001.txt");
		outputFluDataRo(mapa_new,"Flu_ro001.txt");
		outputFluDataP(mapa_new,"Flu_p001.txt");
	}
	cout<<"Enter"<<endl;
	int i;
	cin>>i; 
//	return 0;
}
//outputData
/*string filename="./movie/";
char name[12];
int iter=t,
int len=1;
sprintf_s(name, "%d.txt", t);
while(iter/=10) len++; // определяем количество цифр в номере итерации
for(int i=0;i<filename_count-len;i++) filename+="0"; // добавляем нули в имени файла
filename+=string(name); //добавляем номер итерации в имя файла
outputData(mapa,filename);*/
