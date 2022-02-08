// fRSH2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <stdlib.h>
#include <fcntl.h>
#include <conio.h>
#include <io.h>
#include <sys\stat.h>
#include <math.h>
#include "frelctr.h"
#include <dos.h>
#include <stdio.h>
#include <string.h>

#define M_PI 3.14
#define Rmx 35
//#define Smx=25
//#define Vmin=0
//#define GetmaxX=640;
//#define GetmixX=640;
//#define GetmaxY=480;
//#define GetmixY=480;
//#define GetmaxZ=640;
//#define GetmixZ=640;

bool WorkShow = false;
float DX;
float DY;
int pg = 0;

struct v3
{
    float x, y, z;
};

//-----------€§¬Ґ­пҐ¬лҐ Їа®Ја ¬­лҐ Ї а ¬Ґвал------------//
bool needCharge = true;
bool needLoad = false;
bool needSave = false;
bool needRand = true;

//---------- ”Ё§ЁзҐбЄЁҐ ўҐ«ЁзЁ­л -----------//
int N = 50;
float dt = 1e-6;                //c        Ј Ї® ўаҐ¬Ґ­Ё
float dte = 5e-11;              //c        Ј Ї® ўаҐ¬Ґ­Ё
float Rmax = 13e-7;             //c¬      Њ ЄбЁ¬ «м­л© а ¤Ёгб г з бвЁжл
float Rmin = 13e-7;             //c¬      ЊЁ­Ё¬ «м­л© а ¤Ёгб г з бвЁжл
float Rmid = 100e-7;            //c¬      ‘аҐ¤­ҐҐ а ббв®п­ЁҐ ¬Ґ¦¤г з бвЁж ¬Ё
float Tmshft = 10e-3;           //б       ‚аҐ¬п бў®Ў®¤­®Ј® Їа®ЎҐЈ  з бвЁж
float DensAg = 10.5;            //Ј/c¬^3  Џ«®в­®бвм бҐаҐЎа 
int Tk = 300;                   //K       ’Ґ¬ЇҐа вга 
int MaxQ = 0;                   //Љг«®­   Њ ЄбЁ¬ «м­л© § ап¤
float A = 3.4e-19;              //„¦      ђ Ў®в  ўле®¤  н«ҐЄва®­  Ё§ бҐаҐЎа 
float L = 4e-7;                 //¬       „«Ё­  ў®«­л Ї ¤ ойҐЈ® бўҐв 

//------------ќ«ҐЄваЁзҐбЄЁҐ Є®­бв ­вл-----------------//
float kk = 1.38E-23;             //„¦/K    Џ®бв®п­­ п Ѓ®«мж¬ ­ 
float E0 = 8.854E-14;            //”/c¬    ќ«ҐЄваЁзҐбЄ п Ї®бв®п­­ п
int E = 100;                     //        „Ён«ҐЄваЁзҐбЄ п Їа®­Ёж Ґ¬®бвм
float C_SI = 1 / (4 * M_PI * E0 * E);    //        Љ®нддЁжЁҐ­в ЇҐаҐе®¤  Ё§ ‘ѓ‘ ў ‘€
float eQulon = -1.6E-19;         //Љ«      ‡ ап¤ н«ҐЄва®­ 
float h = 6e-34;                 //„¦*б    Џ®бв®п­­ п Џ« ­Є 
float c = 3e8;                  //¬/б     ‘Є®а®бвм бўҐв 
float Me = 1e-27;                //Ј       Њ бб  н«ҐЄва®­ 

Particle* PNp;
Agregat dF[Npart / 2];
int Pagregat[Npart / 2];
Agregat CMass[Npart / 2];
int ConPat[Npart + 1];
float TmPat[Npart + 1];
int s; float t; float te;
unsigned char VDn = 1;
float OldRX; float OldRY;
Particle* FirstPat = NULL;
Particle* LastPat = NULL;
Particle* TempP;
Particle* Pi;
Particle* Pj;
Particle** Mp = NULL;
float Xmax = 5e-5;
float Ymax = 5e-5 * 480 / 640;
float Zmax = 5e-5;
unsigned short PatInt;
char path[128];
char ch;
bool Se = false;

void MakeArray(int N)
{
    FirstPat = new Particle;
    FirstPat->pred = NULL;
    FirstPat->N = 1;
    LastPat = FirstPat;
    for (int i = 2; i != N + 1; i++)
    {
        LastPat->next = TempP = new Particle;
        TempP->pred = LastPat;
        TempP->next = NULL;
        TempP->N = i;
        LastPat = TempP;
    }
    Mp = new Particle * [N];
}


void DistroyArray()
{
    if (FirstPat == NULL) return;
    TempP = FirstPat->next;
    do
        delete FirstPat;
    while (FirstPat = TempP, TempP = TempP->next, FirstPat != NULL);
    FirstPat = NULL;

    delete[] Mp;
}


Particle* Npat(int Np)
{
    TempP = FirstPat;
    for (int i = 1; i < Np; i++) TempP = TempP->next;
    return (TempP);
}


void SizePatDistr()
{
    int k;
    float sizep[25];
    int Kol[25];
    char fdat[128];
    fdat[0] = 0;
    _fmode = O_TEXT;
    FILE* cl = fopen(fdat, "rt");
    if (cl != NULL)
    {
        for (k = 0; k < 25; k++)
        {
            if (!feof(cl))
            {
                fscanf(cl, "%e %d", &sizep[k], &Kol[k]);
                sizep[k] *= 1e-7;
            }
            else break;
        }
        fclose(cl);
        Pi = FirstPat;
        for (k = 0; k < 25; k++)
        {
            while (Kol[k])
            {
                Pi->R = sizep[k];
                Pi = Pi->next;
                Kol[k]--;
            }
        }
    }
    else
        for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
            if (Pi->N <= N * 0.2) Pi->R = Rmin;
            else if (Pi->N <= N * 0.8) Pi->R = (Rmin + Rmax) / 2;
            else Pi->R = Rmax;
}

void SpontDstr()
{
    float Vel, alpha, theta;
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    {
        Vel = 1.41 * sqrt(2 * kk * Tk / Pi->mass) * float(rand()) / RAND_MAX;
        alpha = 2 * M_PI * float(rand()) / RAND_MAX;
        theta = 2 * M_PI * float(rand()) / RAND_MAX;
        Pi->Vx = Vel * cos(theta) * sin(alpha);
        Pi->Vy = Vel * cos(theta) * cos(alpha);
        Pi->Vz = Vel * sin(theta);
    }
}


void zeroDstr()
{
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next) Pi->Vx = Pi->Vy = Pi->Vz = 0;
}

void MxwDstr()
{
    int Vminz = 0;
    float Vmax = 0.25;
    int kmax;
    if (N != 2) kmax = N / 3; else kmax = 1;
    int i, j;
    double F;
    float dN, V, Vel, dv, alpha, theta;
    int rdn;

    j = 0; i = 1;
    dv = (Vmax - Vminz) / kmax;

    Pj = FirstPat;

    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    {
        V = Vminz + dv * i;
        F = sqrt(Pi->mass / (2 * M_PI * kk * Tk));
        F = pow(F, 3);
        F = F * exp(-(Pi->mass * pow(V, 2)) / (2 * kk * Tk)) * 4 * M_PI * pow(V, 2);
        dN = N * F * dv;
        rdn = (int)(dN + 0.5);
        i++;
        if ((j + rdn) < N)
            while (rdn)
            {
                j++;
                Vel = V;
                alpha = 2 * M_PI * float(rand()) / RAND_MAX;
                theta = 2 * M_PI * float(rand()) / RAND_MAX;
                Pj->Vx = Vel * cos(theta) * sin(alpha);
                Pj->Vy = Vel * cos(theta) * cos(alpha);
                Pj->Vz = Vel * sin(theta);
                rdn--;
                Pj = Pj->next;
            }
    }
    int cnt = N - j;

    while (cnt--)
    {
        Vel = sqrt((8 * kk * Tk) / (M_PI * Pj->mass));
        alpha = 2 * M_PI * float(rand()) / RAND_MAX;
        theta = M_PI * float(rand()) / RAND_MAX;
        Pj->Vx = Vel * cos(theta) * sin(alpha);
        Pj->Vy = Vel * cos(theta) * cos(alpha);
        Pj->Vz = Vel * sin(theta);
        Pj = Pj->next;
    }
}

void InitParticle()
{
    int sign = 0;
    float sq, alpha;

    SizePatDistr();
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    {
        sign = 1 - sign;
    m1: Pi->X = -Xmax + 640 * float(rand()) / RAND_MAX;
        Pi->Y = -Ymax + 480 * float(rand()) / RAND_MAX;
        Pi->Z = -Xmax + 640 * float(rand()) / RAND_MAX;

        for (Pj = FirstPat; Pj != Pi; Pj = Pj->next)
        {
            float dist_ij = pow(Pi->X - Pj->X, 2) + pow(Pi->Y - Pj->Y, 2) + pow(Pi->Z - Pj->Z, 2);
            if (dist_ij <= pow(Pi->R + Pj->R, 2)) goto m1;
        }
        //if (sign) Pi->q=27e-19;else Pi->q=-27e-19;
        Pi->q = 0;
        Pi->mass = 4 / 3 * M_PI * pow(Pi->R, 3) / DensAg;
        //TmPat[Pi->N]=Tmshft*float(rand())/RAND_MAX;
        Pi->Nf = 0;
        Pi->agr = 0;
        Pi->stop = false;
        Pi->Fx = 0;
        Pi->Fy = 0;
        Pi->Fz = 0;
    }
    switch (VDn)
    {
    case 0: MxwDstr(); break;
    case 1: SpontDstr(); break;
    case 2: zeroDstr(); break;
    default:break;
    }
}

Electron* FirstEl = NULL, * LastEl = NULL, * TempEl, * El;

void InitEl()
{
    FirstEl = LastEl = NULL;
}

float Intense = 0.5; //‚в Њ®й­®бвм бўҐв 
float Ef = 1e19; //д®в®­®ў ­  1 „¦
float QuntExit = 1e4;
int Ie = 1;
void ElEmit()
{
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    {
        Pi->Nf += 2 * M_PI * Pi->R * Pi->R * dt * Intense * Ef;
        if (Pi->Nf >= QuntExit)
        {
            Pi->Nf = 0;
            Se = true;
            if (FirstEl == NULL)
            {
                FirstEl = LastEl = new Electron;
                LastEl->Pred = NULL;
                LastEl->Next = NULL;
            }
            else
            {
                LastEl->Next = new Electron;
                LastEl->Next->Pred = LastEl;
                LastEl = LastEl->Next;
                LastEl->Next = NULL;
            }
            float V = sqrt(Pi->Vx * Pi->Vx + Pi->Vy * Pi->Vy + Pi->Vz * Pi->Vz);

            LastEl->X = Pi->X + Pi->R * Pi->Vx / V;
            LastEl->Y = Pi->Y + Pi->R * Pi->Vy / V;
            LastEl->Z = Pi->Z + Pi->R * Pi->Vz / V;

            double Ve = sqrt((h * c / L - A) * 2 / Me);

            LastEl->Vx = Ve * Pi->Vx / V;
            LastEl->Vy = Ve * Pi->Vy / V;
            LastEl->Vz = Ve * Pi->Vz / V;
            LastEl->tr = 0;
            LastEl->Fx = LastEl->Fy = LastEl->Fz = 0;
            Pi->q -= eQulon;
            Ie++;
        }
    }
}

void ExclEl(Electron* ExEl)
{
    if (ExEl->Next != NULL)
        ExEl->Next->Pred = ExEl->Pred; else LastEl = ExEl->Pred;
    if (ExEl->Pred != NULL)
        ExEl->Pred->Next = ExEl->Next; else FirstEl = ExEl->Next;
    delete ExEl;
}

void ElAbsorbe()
{
    float dElPat;
    for (El = FirstEl; El != NULL; El = El->Next)
    {
        if (El->tr <= 2 * dt) continue;
        for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
        {
            dElPat = (Pi->X - El->X) * (Pi->X - El->X) + (Pi->Y - El->Y) * (Pi->Y - El->Y) +
                (Pi->Z - El->Z) * (Pi->Z - El->Z);
            if (dElPat <= Pi->R * Pi->R)
            {
                Pi->q += eQulon;
                ExclEl(El);
                Ie--;
            }
        }
    }
}

void DisposeEl()
{
    if (FirstEl == NULL) return;

    El = FirstEl->Next;
    do
        delete FirstEl;
    while (FirstEl = El, El = El->Next, FirstEl != NULL);
    FirstEl = NULL;
}


void ElMove()
{
    float dVx, dVy, dVz, dx, dy, dz;
    for (El = FirstEl; El != NULL; El = El->Next)
    {
        /*   dVx=dt*El->Fx/Me;
             dVy=dt*El->Fy/Me;
             dVz=dt*El->Fz/Me;

             El->Vx+=dVx;
             El->Vy+=dVy;
             El->Vz+=dVz;*/

        dx = dte * El->Vx;
        dy = dte * El->Vy;
        dz = dte * El->Vz;

        El->X += dx;
        El->Y += dy;
        El->Z += dz;
        if ((El->X > 640) || (El->X < 0)) { El->Vx = -El->Vx; El->X -= dx; }
        if ((El->Y > 480) || (El->Y < 0)) { El->Vy = -El->Vy; El->Y -= dy; }
        if ((El->Z > 640) || (El->Z < 0)) { El->Vz = -El->Vz; El->Z -= dz; }

        El->tr += dt;
    }
}

void InitAgr()
{
    int c;
    for (c = 0; c != N / 2; c++) Pagregat[c] = 0;
    for (c = 0; c != N + 1; c++) ConPat[c] = 0;
}

bool IsPatInAgr(Particle* Pi, Particle* Pj)
{
    return  (bool)((Pi->agr != 0) && (Pi->agr == Pj->agr));
}

void CulonForce()
{
    int Nagr;
    float Temp, dX, dY, dZ;
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    {
        Pi->Fx = 0;
        Pi->Fy = 0;
        Pi->Fz = 0;
    }
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    {
        for (Pj = Pi->next; Pj != NULL; Pj = Pj->next)
        {
            if (IsPatInAgr(Pi, Pj)) continue;
            dX = Pi->X - Pj->X;
            dY = Pi->Y - Pj->Y;
            dZ = Pi->Z - Pj->Z;
            Temp = Pi->q * Pj->q / pow(dX * dX + dY * dY + dZ * dZ, 1.5);
            Pj->Fx -= Temp * dX;
            Pj->Fy -= Temp * dY;
            Pj->Fz -= Temp * dZ;
            Pi->Fx += Temp * dX;
            Pi->Fy += Temp * dY;
            Pi->Fz += Temp * dZ;
        }
        Pi->Fx *= C_SI;
        Pi->Fy *= C_SI;
        Pi->Fz *= C_SI;
        //  float U=Pi->Fx*dX+Pi->Fy*dY+Pi->Fz*dZ;
        //  float Ek=Pi->mass*(pow(Pi->Vx,2)+pow(Pi->Vy,2)+pow(Pi->Vz,2))/2;
    }
}

void CulPatEl()
{
    float Temp, dX, dY, dZ;
    for (El = FirstEl; El != NULL; El = El->Next) El->Fx = El->Fy = El->Fz = 0;

    for (El = FirstEl; El != NULL; El = El->Next)
    {
        for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
        {
            dX = El->X - Pi->X;
            dY = El->Y - Pi->Y;
            dZ = El->Z - Pi->Z;
            Temp = eQulon * Pi->q / pow(dX * dX + dY * dY + dZ * dZ, 1.5);
            El->Fx += Temp * dX;
            El->Fy += Temp * dY;
            El->Fz += Temp * dZ;
        }
        El->Fx *= C_SI;
        El->Fy *= C_SI;
        El->Fz *= C_SI;
    }
}


#define Pprt(a) (*(Particle**)a)
int sortOfZ(const void* a, const void* b)
{
    if (Pprt(a)->Z > Pprt(b)->Z) return  1; else
        if (Pprt(a)->Z < Pprt(b)->Z) return -1; else
            return  0;
}

//void ShowPicture()
//{
//    char st[100];
//    float Dzx, Dzy;
//    if (WorkShow) return;
//    WorkShow = true;
//
//    for (El = FirstEl; El != NULL; El = El->Next)
//    {
//        Point(LIGHTCYAN, El->X + DX, El->Y + DY);
//    }
//    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
//    {
//        if (Pi->agr != 0)
//        {
//            Dzx = Pi->Z * (Pi->X - CMass[Pi->agr].X) / Zmax;
//            Dzy = Pi->Z * (Pi->Y - CMass[Pi->agr].Y) / Zmax;
//        }
//        else Dzx = Dzy = 0;
//
//        if (Pi->q == 0) Circle(MAGENTA, Pi->X + Dzx + DX, Pi->Y + Dzy + DY, Pi->R + Pi->Z * Pi->R / Zmax);
//        else
//            if (Pi->q < 0)
//            {
//                Line(BLUE, Pi->X + Dzx + DX - 2.0 / rangeX, Pi->Y + Dzy + DY, Pi->X + Dzx + 2.0 / rangeX + DX, Pi->Y + Dzy + DY);
//                Circle(BLUE, Pi->X + Dzx + DX, Pi->Y + Dzy + DY, Pi->R + Pi->Z * Pi->R / Zmax);
//            }
//            else
//            {
//                Line(RED, Pi->X + Dzx + DX - 2.0 / rangeX, Pi->Y + Dzy + DY, Pi->X + Dzx + 2.0 / rangeX + DX, Pi->Y + Dzy + DY);
//                Line(RED, Pi->X + Dzx + DX, Pi->Y - 2.0 / rangeY + Dzy + DY, Pi->X + Dzx + DX, Pi->Y + Dzy + 2.0 / rangeY + DY);
//                Circle(RED, Pi->X + Dzx + DX, Pi->Y + Dzy + DY, Pi->R + Pi->Z * Pi->R / Zmax);
//            }
//
//        /* Setcolor (WHITE);
//           itoa (Pi->N,st,10);
//           OutTextXY(Pi->X+DX, Pi->Y+DY, st);
//        */
//        if (!Pi->stop) Line(LIGHTGRAY, Pi->X + DX, Pi->Y + DY, Pi->X + DX + Pi->Fx * 1e10, Pi->Y + DY + Pi->Fy * 1e10);
//    }
//    for (int i = 1; i < s; i++)
//    {
//        Circle(GREEN, CMass[i].X + DX, CMass[i].Y + DY, 1.0 / rangeR);
//        Line(LIGHTGRAY, CMass[i].X + DX, CMass[i].Y + DY, CMass[i].X + DX + dF[i].X * 1e10, CMass[i].Y + DY + dF[i].Y * 1e10);
//    }
//
//    Rectangle(GREEN, GetminX + DX, GetminY + DY, 640 + DX, 480 + DY);
//    drwline(SET, RED, maxx / 2 - 2, maxy / 2, maxx / 2 + 2, maxy / 2);
//    drwline(SET, RED, maxx / 2, maxy / 2 - 2, maxx / 2, maxy / 2 + 2);
//    drwbox(SET, RED, maxx / 2 - 4, maxy / 2 + 4, maxx / 2 + 4, maxy / 2 - 4);
//
//    pagedisplay(0, 0, pg);
//    pg = 1 - pg;
//    pageactive(pg);
//    fillpage(0);
//    WorkShow = false;
//}
//
//bool int9Work = false;

//void far finishInt9h()
//{
//    asm{
//            in al,61h
//        mov ah,al
//            or al,80h
//            out 61h,al
//            xchg ah,al
//            out 61h,al
//            mov al,20h
//            out 20h,al
//    }
//}

//void interrupt ScPict(...)
//{
//    unsigned char  scanCode, keyCode;
//    unsigned char  fpuState[108];
//
//    scanCode = inport(0x60);
//    keyCode = (scanCode & 0x7F);
//
//    if (int9Work) { finishInt9h(); return; }
//    int9Work = true;

    //switch (keyCode)
    //{
    //case 78: case 74: case 75: case 77: case 80:
    //case 73: case 79: case 81: case 71: case 76:
    //case 72: case 1:
    //    if ((scanCode & 0x80))
    //    {
    //        finishInt9h();
    //        int9Work = false;
    //        return;
    //    } break;
    //default:
    //{ KeyInt();
    //int9Work = false;
    //return;
    //}
    //}

        //switch (keyCode)
        //{
        //case 78:
        //{ DX -= maxx / (2 * rangeX);
        //DY += maxy / (2 * rangeY);
        //rangeX *= 1.05;
        //rangeY *= 1.05;
        //DX += maxx / (2 * rangeX);
        //DY -= maxy / (2 * rangeY);
        //} break;
        //case 74:
        //{ DX -= maxx / (2 * rangeX);
        //DY += maxy / (2 * rangeY);
        //rangeX *= (1.0 / 1.05);
        //rangeY *= (1.0 / 1.05);
        //DX += maxx / (2 * rangeX);
        //DY -= maxy / (2 * rangeY);
        //} break;
        //case 76:
        //{ DX = 0; DY = 0;
        //rangeX = OldRX; rangeY = OldRY;
        //} break;
        //case 75: DX += GetScX / 150; break;
        //case 77: DX -= GetScX / 150; break;
        //case 80: DY += GetScY / 150; break;
        //case 72: DY -= GetScY / 150; break;
        //case 71: DX += GetScX / 5; break;
        //case 79: DX -= GetScX / 5; break;
        //case 81: DY += GetScY / 5; break;
        //case 73: DY -= GetScY / 5; break;
        //case 1: ch = 'q'; break;
        //}

    //   ShowPicture();
//}

//void RandWalk(float& TmP, float& Vlx, float& Vly, float& Vlz)
//{
//    int sign;
//    if (TmP < t)
//    {
//        sign = (2 * random(2) - 1);
//        Vlx = Vlx * sign;
//        sign = (2 * random(2) - 1);
//        Vly = Vly * sign;
//        sign = (2 * random(2) - 1);
//        Vlz = Vlz * sign;
//        TmP += Tmshft * float(rand()) / RAND_MAX;
//    }
//}

void MovePart()
{
    int i, j;
    char st[10];
    float dVx, dVy, dVz, dx, dy, dz;

    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
        if (Pi->stop == false)
        {
            WorkShow = true;
            dVx = dt * Pi->Fx / Pi->mass;
            dVy = dt * Pi->Fy / Pi->mass;
            dVz = dt * Pi->Fz / Pi->mass;
            Pi->Vx += dVx;
            Pi->Vy += dVy;
            Pi->Vz += dVz;
            // RandWalk(TmPat[Pi->N],Pi->Vx,Pi->Vy,Pi->Vz);

            Pi->X += dx = dt * Pi->Vx;
            Pi->Y += dy = dt * Pi->Vy;
            Pi->Z += dz = dt * Pi->Vz;

            if ((Pi->X > 640) || (Pi->X < 0)) { Pi->Vx = -Pi->Vx; Pi->X -= dx; }
            if ((Pi->Y > 480) || (Pi->Y < 0)) { Pi->Vy = -Pi->Vy; Pi->Y -= dy; }
            if ((Pi->Z > 640) || (Pi->Z < 0)) { Pi->Vz = -Pi->Vz; Pi->Z -= dz; }

            WorkShow = false;
        }
}

void AngleSpeed(int kk, RPoint& W)
{
    // Џ®«­л© ¬®¬Ґ­в бЁ«, ¤Ґ©бвўгойЁ© ­   ЈаҐЈ в б® бв®а®­л ®бв «м­ле з бвЁж
    float Nx = 0; float Ny = 0; float Nz = 0;

    int i = Pagregat[kk];
    while (i)
    {
        Nx += (Pi->Y - CMass[kk].Y) * Pi->Fz - (Pi->Z - CMass[kk].Z) * Pi->Fy;
        Ny += (Pi->Z - CMass[kk].Z) * Pi->Fx - (Pi->X - CMass[kk].X) * Pi->Fz;
        Nz += (Pi->X - CMass[kk].X) * Pi->Fy - (Pi->Y - CMass[kk].Y) * Pi->Fx;
        i = ConPat[i];
    }
    CMass[kk].Mx += dt * Nx;
    CMass[kk].My += dt * Ny;
    CMass[kk].Mz += dt * Nz;

    W.X = CMass[kk].Mx / CMass[kk].Jx;
    W.Y = CMass[kk].My / CMass[kk].Jy;
    W.Z = CMass[kk].Mz / CMass[kk].Jz;
}

// ‚лзЁб«пҐв жҐ­ва ¬ бб, Ї®«­л© ¬®¬Ґ­в, ¬®¬Ґ­в Ё­ҐажЁЁ  ЈаҐЈ в  Nagr
void ChangeCMass(int Nagr)
{
    int Np = Pagregat[Nagr];
    float Mass = 0;
    CMass[Nagr].X = 0;
    CMass[Nagr].Y = 0;
    CMass[Nagr].Z = 0;

    CMass[Nagr].Vx = 0;
    CMass[Nagr].Vy = 0;
    CMass[Nagr].Vz = 0;
    // ‚лзЁб«Ґ­ЁҐ жҐ­ва  ¬ бб  ЈаҐЈ в  Nagr

    while (Np)
    {
        PNp = Npat(Np);
        Mass += PNp->mass;
        CMass[Nagr].X += PNp->mass * PNp->X;
        CMass[Nagr].Y += PNp->mass * PNp->Y;
        CMass[Nagr].Z += PNp->mass * PNp->Z;

        CMass[Nagr].Vx += PNp->mass * PNp->Vx;
        CMass[Nagr].Vy += PNp->mass * PNp->Vy;
        CMass[Nagr].Vz += PNp->mass * PNp->Vz;

        Np = ConPat[Np];
    }

    CMass[Nagr].X /= Mass;
    CMass[Nagr].Y /= Mass;
    CMass[Nagr].Z /= Mass;

    CMass[Nagr].Vx /= Mass;
    CMass[Nagr].Vy /= Mass;
    CMass[Nagr].Vz /= Mass;

    // ‚лзЁб«Ґ­ЁҐ Ї®«­®Ј® ¬®¬Ґ­в  M Ё ¬®¬Ґ­в  Ё­ҐажЁЁ J  ЈаҐЈ в  Nagr

    Np = Pagregat[Nagr];

    CMass[Nagr].Jx = 0;
    CMass[Nagr].Mx = CMass[Nagr].My = CMass[Nagr].Mz = 0;

    while (Np)
    {
        PNp = Npat(Np);
        CMass[Nagr].Jx += PNp->mass * ((pow(PNp->Z - CMass[Nagr].Z, 2) + pow(PNp->Y - CMass[Nagr].Y, 2)) + 0.4 * pow(PNp->R, 2));
        CMass[Nagr].Jy += PNp->mass * ((pow(PNp->X - CMass[Nagr].X, 2) + pow(PNp->Z - CMass[Nagr].Z, 2)) + 0.4 * pow(PNp->R, 2));
        CMass[Nagr].Jz += PNp->mass * ((pow(PNp->X - CMass[Nagr].X, 2) + pow(PNp->Y - CMass[Nagr].Y, 2)) + 0.4 * pow(PNp->R, 2));

        CMass[Nagr].Mx += PNp->mass * ((PNp->Y - CMass[Nagr].Y) * (PNp->Vz - CMass[Nagr].Vz) - (PNp->Z - CMass[Nagr].Z) * (PNp->Vy - CMass[Nagr].Vy));
        CMass[Nagr].My += PNp->mass * ((PNp->Z - CMass[Nagr].Z) * (PNp->Vx - CMass[Nagr].Vx) - (PNp->X - CMass[Nagr].X) * (PNp->Vz - CMass[Nagr].Vz));
        CMass[Nagr].Mz += PNp->mass * ((PNp->X - CMass[Nagr].X) * (PNp->Vy - CMass[Nagr].Vy) - (PNp->Y - CMass[Nagr].Y) * (PNp->Vx - CMass[Nagr].Vx));

        Np = ConPat[Np];
    }
}

// Џа®жҐ¤гал ¤«п а Ў®вл б  ЈаҐЈ в ¬Ё
float GetMass(int i)
{
    float ret = 0;
    if (Pi->stop)
    {
        i = Pagregat[Pi->agr];
        do
        {
            ret += Pi->mass;
            i = ConPat[i];
        } while (i != 0);
    }
    else ret = Pi->mass;
    return ret;
}

// “бв ­®ўЄ  «Ё­Ґ©­ле бЄ®а®бвҐ© г з бвЁж ў  ЈаҐЈ вҐ б®¤Ґа¦ йЁе i
void SetAgrSpeed(Particle* Pi)
{
    RPoint W;
    if (Pi->stop)
    {
        //‚лзЁб«Ґ­ЁҐ «Ё­Ґ©­®© бЄ®а®бвЁ з бвЁж ў  ЈаҐЈ вҐ }
        int i = Pagregat[Pi->agr];
        AngleSpeed(Pi->agr, W);
        do
        {
            Pi->Vx = CMass[Pi->agr].Vx + W.Y * (Pi->Z - CMass[Pi->agr].Z) - W.Z * (Pi->Y - CMass[Pi->agr].Y);
            Pi->Vy = CMass[Pi->agr].Vy + W.Z * (Pi->X - CMass[Pi->agr].X) - W.X * (Pi->Z - CMass[Pi->agr].Z);
            Pi->Vz = CMass[Pi->agr].Vz + W.X * (Pi->Y - CMass[Pi->agr].Y) - W.Y * (Pi->X - CMass[Pi->agr].X);
            i = ConPat[i];
        } while (i != 0);
    }
}

//  „®Ў ў«пҐв j з бвЁжг Є  ЈаҐЈ вг ­ зЁ­ ойҐ¬гбп б i з бвЁжл
void AddPattoAgr(Particle* Pi, Particle* Pj)
{
    Pj->agr = Pi->agr;
    int i = Pi->N; int j = Pj->N;
    while (ConPat[i]) i = ConPat[i];

    ConPat[i] = j;
    if (ConPat[j] != 0)
    {
        j = ConPat[j];
        while (j)
        {
            Pj = Npat(j);
            Pj->agr = Pi->agr;
            j = ConPat[j];
        }
    }
}

//      ђa§¤ўЁЈ Ґв ¤ўҐ з бвЁжл
void PushAway(Particle* Pi, Particle* Pj)
{
    float dX, dY, dZ, R;

    R = Pi->R + Pj->R;
    dX = Pj->X - Pi->X;
    dY = Pj->Y - Pi->Y;
    dZ = Pj->Z - Pi->Z;

    float dist_ij = sqrt(dX * dX + dY * dY + dZ * dZ);

    dX *= (R - dist_ij) / dist_ij;
    dY *= (R - dist_ij) / dist_ij;
    dZ *= (R - dist_ij) / dist_ij;

    if (Pj->agr != 0)
    {
        CMass[Pj->agr].X += dX;
        CMass[Pj->agr].Y += dY;
        CMass[Pj->agr].Z += dZ;
        int j = Pagregat[Pj->agr];
        while (j)
        {
            Pj = Npat(j);
            Pj->X += dX;
            Pj->Y += dY;
            Pj->Z += dZ;

            j = ConPat[j];
        }
    }
    else
    {
        Pj->X += dX;
        Pj->Y += dY;
        Pj->Z += dZ;
    }
}

void UnitPaticle(Particle* Pi, Particle* Pj)
{
    int l, kk, m;
    int i = Pi->N;
    int j = Pj->N;

    Particle* Pm; Particle* Pl;

    if (IsPatInAgr(Pi, Pj)) return;


    PushAway(Pi, Pj);

    SetAgrSpeed(Pi);
    SetAgrSpeed(Pj);

    if (Pi->stop == false)
    {
        if (Pj->stop == false)
        {
            Pagregat[s] = Pi->N;
            Pi->agr = s;
            Pj->agr = s;
            s++;
            ConPat[i] = j;
        }
        else AddPattoAgr(Pj, Pi);

    }
    else
        if (Pj->stop == false) AddPattoAgr(Pi, Pj);
        else
        {
            kk = Pi->agr;
            l = Pagregat[kk];
            s--;
            Pagregat[kk] = Pagregat[s];
            Pagregat[s] = 0;
            CMass[kk] = CMass[s];
            CMass[s].X = CMass[s].Y = CMass[s].Z =
                CMass[s].Vx = CMass[s].Vy = CMass[s].Vz =
                CMass[s].Mx = CMass[s].My = CMass[s].Mz =
                CMass[s].Jx = CMass[s].Jy = CMass[s].Jz = 0;
            m = Pagregat[kk];
            while (m)
            {
                Pm = Npat(m);
                Pm->agr = kk;
                m = ConPat[m];
            }
            Pl = Npat(l);
            AddPattoAgr(Pj, Pl);
        }

    // Џ®¬ҐвЁвм i Ё j Є Є ЇаЁ­ ¤«Ґ¦ йЁҐ Є  ЈаҐЈ в ¬
    Pi->stop = true;
    Pj->stop = true;
    ChangeCMass(Pj->agr);
}

void AgrForces()
{
    int i, kk;
    float mass, dVx, dVy, dVz, dX, dY, dZ, M1;
    for (kk = 1; kk < s; kk++)
    {
        i = Pagregat[kk];
        dF[kk].X = dF[kk].Y = dF[kk].Z = mass = 0;

        while (i)
        {
            Pi = Npat(i);
            dF[kk].X += Pi->Fx;
            dF[kk].Y += Pi->Fy;
            dF[kk].Z += Pi->Fz;
            mass += Pi->mass;
            i = ConPat[i];
        }
        if ((CMass[kk].X > 640) || (CMass[kk].X < 0))  CMass[kk].Vx = -CMass[kk].Vx;
        if ((CMass[kk].Y > 480) || (CMass[kk].Y < 0))  CMass[kk].Vy = -CMass[kk].Vy;
        if ((CMass[kk].Z > 640) || (CMass[kk].Z < 0))  CMass[kk].Vz = -CMass[kk].Vz;

        dVx = dt * dF[kk].X / mass;
        dVy = dt * dF[kk].Y / mass;
        dVz = dt * dF[kk].Z / mass;
        CMass[kk].Vx += dVx;
        CMass[kk].Vy += dVy;
        CMass[kk].Vz += dVz;
        CMass[kk].X += dt * CMass[kk].Vx;
        CMass[kk].Y += dt * CMass[kk].Vy;
        CMass[kk].Z += dt * CMass[kk].Vz;

        //-------- ‚лзЁб«Ґ­ЁҐ гЈ«  Phi --------//
        RPoint phi;
        AngleSpeed(kk, phi);
        phi.X = dt * phi.X;
        phi.Y = dt * phi.Y;
        phi.Z = dt * phi.Z;
        //  ‚а йҐ­ЁҐ  ЈаҐЈ в  Ї®баҐ¤бвў®¬ ЇаҐ®Ўа §®ў ­Ёп Є®®а¤Ё­ в б Ї®¬®ймо---//
                             //---- ¬ ваЁжл ўа йҐ­Ёп ---//

        i = Pagregat[kk];
        while (i)
        {
            Pi = Npat(i);
            Pi->X += dt * (CMass[kk].Vx + dVx);
            Pi->Y += dt * (CMass[kk].Vy + dVy);
            Pi->Z += dt * (CMass[kk].Vz + dVz);
            dX = Pi->X - CMass[kk].X;
            dY = Pi->Y - CMass[kk].Y;
            Pi->X = CMass[kk].X + dX * cos(phi.Z) - dY * sin(phi.Z);
            Pi->Y = CMass[kk].Y + dY * cos(phi.Z) + dX * sin(phi.Z);
            dY = Pi->Y - CMass[kk].Y;
            dZ = Pi->Z - CMass[kk].Z;
            Pi->Y = CMass[kk].Y + dY * cos(phi.X) - dZ * sin(phi.X);
            Pi->Z = CMass[kk].Z + dZ * cos(phi.X) + dY * sin(phi.X);

            dX = Pi->X - CMass[kk].X;
            dZ = Pi->Z - CMass[kk].Z;
            Pi->X = CMass[kk].X + dX * cos(phi.Y) - dZ * sin(phi.Y);
            Pi->Z = CMass[kk].Z + dZ * cos(phi.Y) + dX * sin(phi.Y);
            i = ConPat[i];
        }
    }
}

// ‡ ЇЁб вм б®бв®п­ЁҐ  ЈаҐЈ в®ў ў д ©«

void Save()
{
    char fdat[128]; fdat[0] = 0;
    //strcpy(fdat, PathToDat());
    strcat(fdat, "fractal.dat");
    _fmode = O_BINARY;
    int f = creat(fdat, S_IWRITE);
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    write(f, Pi, sizeof(Particle));
    write(f, CMass, sizeof(CMass));
    write(f, Pagregat, sizeof(Pagregat));
    write(f, ConPat, sizeof(ConPat));
    write(f, &s, sizeof(s));
    write(f, &dt, sizeof(t));
    write(f, &t, sizeof(t));
    close(f);

    fdat[0] = 0;
    //strcpy(fdat, PathToDat());
    strcat(fdat, "fract.dat");
    FILE* F = fopen(fdat, "w+");
    for (int i = 1; i != s; i++)
    {
        int j = Pagregat[i];
        while (j != 0)
        {
            Pj = Npat(j);
            fprintf(F, "%+e    %+e    %+e\xD\n", Pj->X, Pj->Y, Pj->Z);
            j = ConPat[j];
        }
        fprintf(F, "\xD\n");
    }
    fclose(F);

    fdat[0] = 0;
    //strcpy(fdat, PathToDat());
    strcat(fdat, "frsp.dat");
    FILE* F1 = fopen(fdat, "w+");
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
        fprintf(F1, "%+e    %+e    %+e\xD\n", Pi->X, Pi->Y, Pi->Z);
    fclose(F1);
}

void Load()
{
    _fmode = O_BINARY;
    char fdat[128]; fdat[0] = 0;
    //strcpy(fdat, PathToDat());
    strcat(fdat, "fractal.dat");
    int f = open(fdat, S_IREAD);

    TempP = new Particle;
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
    {
        TempP->next = Pi->next;
        TempP->pred = Pi->pred;
        read(f, Pi, sizeof(Particle));
        Pi->next = TempP->next;
        Pi->pred = TempP->pred;
    }
    delete TempP;

    read(f, CMass, sizeof(CMass));
    read(f, Pagregat, sizeof(Pagregat));
    read(f, ConPat, sizeof(ConPat));
    read(f, &s, sizeof(s));
    read(f, &dt, sizeof(t));
    read(f, &t, sizeof(t));
    close(f);
}

int NumPatOutAgr()
{
    int i = 0;
    for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
        if (!Pi->stop) i++;
    return i;
}
int NumPatInAgr()
{
    return N - NumPatOutAgr();
}

float tpa;
int main()
{
    //unsigned long mb = coreleft();
    //struct ffblk ffblk;
    char fdat[128];
    fdat[0] = 0;
    //strcpy(fdat, PathToDat()); strcat(fdat, "numpat.dat");
    _fmode = O_TEXT;
    FILE* cl;
    ch = 'n';
    int ti = 0;
    DX = 0; DY = 0;
    te = 0; t = 0; s = 1;
    //PaletteData c, nm;
    //palget(nm, 0, 15);
    //for (int i = 0; i < 64; i++) { c[i + 16].r = ((i * i) / 64); c[i + 16].g = 0; c[i + 16].b = 0; }
    //for (i = 0; i < 64; i++) { c[i + 64 + 16].r = ((i * i) / 64); c[i + 64 + 16].g = 0; c[i + 64 + 16].b = ((i * i) / 64); }
    //for (i = 0; i < 64; i++) { c[i + 128 + 16].r = 0; c[i + 128 + 16].g = 0; c[i + 128 + 16].b = ((i * i) / 64); }
    //palcopy(nm, c, 0, 15);
    //InitGraph(-Xmax, -Ymax, Xmax, Ymax);
    //OldRX = rangeX; OldRY = rangeY;
    //KeyInt = getvect(9);
    //setvect(9, ScPict);
    //randomize();
    InitEl();
    InitAgr();
    MakeArray(N);
    InitParticle();

    do
    {
        if (!Se)
        {
            ElEmit();
            CulonForce();
            for (Pi = FirstPat; Pi != NULL; Pi = Pi->next)
                for (Pj = Pi->next; Pj != NULL; Pj = Pj->next)
                {
                    float dist_ij = pow(Pi->X - Pj->X, 2) + pow(Pi->Y - Pj->Y, 2) + pow(Pi->Z - Pj->Z, 2);
                    if (dist_ij <= pow((Pi->R + Pj->R), 2)) UnitPaticle(Pi, Pj);
                }
            AgrForces();
            MovePart();
            t += dt;
        }
        else
        {
            ElMove();
            te += dte;
            if ((te >= dt) || (Ie == 1)) { te = 0; Se = false; }
        }
        ElAbsorbe();
        //ShowPicture(); 
        break;
        if (kbhit()) ch = getch();
        if ((t > 1e-2) && (s == 2)) ch = 'q';
    } while (ch != 'q');

    getch();
    


    DistroyArray();
    DisposeEl();
    //gotoxy(20, 13);
    printf("Aggregation time = %.3e\n", t);
   // gotoxy(1, 20);
    //printf("memory lost %lu\n", mb - coreleft());
    printf("Number free electrons=%d", Ie--);
    getch();

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
