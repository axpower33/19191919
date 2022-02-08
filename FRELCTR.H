const int Npart=300;

struct RPoint
		  {float X;
		   float Y;
		   float Z;};
struct
		   Particle{
		   float X,Y,Z;
		   float Vx,Vy,Vz;
		   float Fx,Fy,Fz;
		   float Nf;
		   float mass;
		   float q;
		   float R;
		   int	  N;
		   int	  agr;
		   bool stop;
		   Particle* next;
		   Particle* pred;};

struct Agregat{
		 float X,Y,Z;
		 float Vx,Vy,Vz;
		 float Jx,Jy,Jz;
		 float Mx,My,Mz;
	      };
struct Electron{
		float X,Y,Z;
		float Vx,Vy,Vz;
		float Fx,Fy,Fz;
		float tr;
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
extern float dt;
extern	float Rmin;
extern	float Rmax;
extern	float Rmid;
extern	float Tmshft;
extern	float DensAg;
extern	float Xmin;
extern	float Xmax;
extern	float Ymin;
extern	float Ymax;
extern	float Zmax;
extern	int Tk;
extern	int MaxQ;
extern	int Kmax;
extern	int ScR;
extern	float Koeff;
extern	float RNG;
extern	float Xmn;
extern	float Xmx;
extern	float del;
extern  int countCall;
extern	float Lmin;
extern	float Lmax;
*/