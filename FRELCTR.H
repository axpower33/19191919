const int Npart=300;

struct RPoint
		  {double X;
		   double Y;
		   double Z;};
struct
		   Particle{
		   double X,Y,Z;
		   double Vx,Vy,Vz;
		   double Fx,Fy,Fz;
		   double Nf;
		   double mass;
		   double q;
		   double R;
		   int	  N;
		   int	  agr;
		   bool stop;
		   Particle* next;
		   Particle* pred;};

struct Agregat{
		 double X,Y,Z;
		 double Vx,Vy,Vz;
		 double Jx,Jy,Jz;
		 double Mx,My,Mz;
	      };
struct Electron{
		double X,Y,Z;
		double Vx,Vy,Vz;
		double Fx,Fy,Fz;
		double tr;
		Electron* Next;
		Electron* Pred;
		};

/* extern Particle *PNp,*TempP,*Pi,*Pj,*Pk,*FirstPat;
extern Agregat dF[Npart/2] ;
extern int Pagregat[Npart/2];
extern Agregat CMass[Npart/2];
extern int ConPat[Npart+1];
extern unsigned short LdSv;
extern unsigned short PatInt;
extern bool needCharge;
extern bool needCulon;
extern bool needVDVaalse;
extern bool needLoad;
extern bool needSave;
extern unsigned char VDn;
extern int N;
extern int nn;
extern int s;
extern double dt;
extern	double Rmin;
extern	double Rmax;
extern	double Rmid;
extern	double Tmshft;
extern	double DensAg;
extern	double Xmin;
extern	double Xmax;
extern	double Ymin;
extern	double Ymax;
extern	double Zmax;
extern	int Tk;
extern	int MaxQ;
extern	int Kmax;
extern	int ScR;
extern	double Koeff;
extern	double RNG;
extern	double Xmn;
extern	double Xmx;
extern	double del;
extern  int countCall;
extern	double Lmin;
extern	double Lmax;
*/