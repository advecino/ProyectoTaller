#include "../include/JPL_Eph_DE430.h"
#include "../include/Cheb3D.h"
#include <cmath>
#include <stdexcept>
#include <vector>
#include "../include/global.h"

/*
%--------------------------------------------------------------------------
%
% JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
%   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Notes: Light-time is already taken into account
%
% Last modified:   2018/01/11   M. Mahooti
%
%--------------------------------------------------------------------------*/


#define IDX(i,j) PC(i,j)

// Corta una subsecuencia de PCtemp en un vector<double>
// Cambia slice por esta versión:
static void slice(const Matrix& row, int start, int len, std::vector<double>& out) {
    // Aseguramos capacidad
    out.reserve(out.size() + len);
    // row es una matriz 1×N, usamos siempre fila 1 y columnas start..start+len-1
    for (int k = 0; k < len; ++k) {
        out.push_back(row(1, start + k));
    }
}


PlanetaryPositions JPL_Eph_DE430(double Mjd_TDB) {
    PlanetaryPositions pos;

    // --- Sección 1: buscar fila válida en PC ---
    double JD = Mjd_TDB + 2400000.5;
    int nRows = PC.getFilas(), nCols = PC.getColumnas();
    int i = -1;
    for (int r = 1; r <= nRows; ++r) {
        if (PC(r,1) <= JD && JD <= PC(r,2)) { i = r; break; }
    }
    if (i < 1) throw std::runtime_error("JD fuera de rango en PC");
    Matrix PCtemp = PC.getSubMatrix(i, i, 1, nCols);

    double t1 = PCtemp(1,1) - 2400000.5;
    double dt = Mjd_TDB - t1;

    // --- Sección 2: Tierra (Earth) ---
    // bloques de 13 coeficientes, dos subseqs de 16 días cada uno
    std::vector<double> Cx_E, Cy_E, Cz_E, tmp;
    int base = 231; // inicio MATLAB
    slice(PCtemp, base, 13, Cx_E);
    slice(PCtemp, base+13, 13, Cy_E);
    slice(PCtemp, base+26, 13, Cz_E);
    slice(PCtemp, base+39, 13, tmp);  Cx_E.insert(Cx_E.end(), tmp.begin(), tmp.end()); tmp.clear();
    slice(PCtemp, base+52, 13, tmp);  Cy_E.insert(Cy_E.end(), tmp.begin(), tmp.end()); tmp.clear();
    slice(PCtemp, base+65, 13, tmp);  Cz_E.insert(Cz_E.end(), tmp.begin(), tmp.end());

    int j= (dt>16 ? 1 : 0);
    double Mjd0 = t1 + 16*j;
    Matrix mCx_E(13*2,1), mCy_E(13*2,1), mCz_E(13*2,1);
    for(int k=0;k<13*2;++k){
        mCx_E(k+1,1)=Cx_E[k];
        mCy_E(k+1,1)=Cy_E[k];
        mCz_E(k+1,1)=Cz_E[k];
    }
    pos.r_Earth = 1e3 * Cheb3D(Mjd_TDB,13,Mjd0,Mjd0+16, mCx_E,mCy_E,mCz_E);

    // --- Sección 3: Luna (Moon) ---
    std::vector<double> Cx_M, Cy_M, Cz_M;
    base = 441;
    slice(PCtemp, base, 13, Cx_M);
    slice(PCtemp, base+13, 13, Cy_M);
    slice(PCtemp, base+26, 13, Cz_M);
    for(int rep=1;rep<8;++rep){
        slice(PCtemp, base+39*rep, 13, tmp);  Cx_M.insert(Cx_M.end(), tmp.begin(), tmp.end()); tmp.clear();
        slice(PCtemp, base+39*rep+13, 13, tmp); Cy_M.insert(Cy_M.end(), tmp.begin(), tmp.end()); tmp.clear();
        slice(PCtemp, base+39*rep+26, 13, tmp); Cz_M.insert(Cz_M.end(), tmp.begin(), tmp.end()); tmp.clear();
    }
    // Subintervalos de 4 días
    j = int(dt/4); j = std::min(j,7);
    Mjd0 = t1 + 4*j;
    Matrix mCx_M(13*8,1), mCy_M(13*8,1), mCz_M(13*8,1);
    for(int k=0;k<13*8;++k){
        mCx_M(k+1,1)=Cx_M[k];
        mCy_M(k+1,1)=Cy_M[k];
        mCz_M(k+1,1)=Cz_M[k];
    }
    pos.r_Moon = 1e3 * Cheb3D(Mjd_TDB,13,Mjd0,Mjd0+4,mCx_M,mCy_M,mCz_M);

    // --- Sección 4: Sol (Sun) ---
    std::vector<double> Cx_Su, Cy_Su, Cz_Su; tmp.clear();
    base = 753;
    slice(PCtemp, base, 11, Cx_Su);
    slice(PCtemp, base+11,11, Cy_Su);
    slice(PCtemp, base+22,11, Cz_Su);
    slice(PCtemp, base+33,11,tmp);  Cx_Su.insert(Cx_Su.end(), tmp.begin(), tmp.end()); tmp.clear();
    slice(PCtemp, base+44,11,tmp);  Cy_Su.insert(Cy_Su.end(), tmp.begin(), tmp.end()); tmp.clear();
    slice(PCtemp, base+55,11,tmp);  Cz_Su.insert(Cz_Su.end(), tmp.begin(), tmp.end());
    j = (dt>16 ? 1 : 0);
    Mjd0 = t1 + 16*j;
    Matrix mCx_S(11*2,1), mCy_S(11*2,1), mCz_S(11*2,1);
    for(int k=0;k<11*2;++k){
        mCx_S(k+1,1)=Cx_Su[k];
        mCy_S(k+1,1)=Cy_Su[k];
        mCz_S(k+1,1)=Cz_Su[k];
    }
    pos.r_Sun = 1e3 * Cheb3D(Mjd_TDB,11,Mjd0,Mjd0+16,mCx_S,mCy_S,mCz_S);

    // --- Secciones 5–11: Mercurio, Venus, Marte, Júpiter, Saturno, Urano, Neptuno, Plutón ---
    // Repite exactamente el mismo patrón:
    auto loadBody=[&](int bstart,int block,int reps,double span){
        std::vector<double> Cx, Cy, Cz,tmp2;
        slice(PCtemp, bstart, block, Cx);
        slice(PCtemp, bstart+block, block, Cy);
        slice(PCtemp, bstart+2*block, block, Cz);
        for(int r=1;r<reps;++r){
            slice(PCtemp,bstart+reps*block+ (r-1)*(block+block+block),block,tmp2);
            Cx.insert(Cx.end(), tmp2.begin(), tmp2.end()); tmp2.clear();
            slice(PCtemp,bstart+reps*block+block+(r-1)*(3*block),block,tmp2);
            Cy.insert(Cy.end(), tmp2.begin(), tmp2.end()); tmp2.clear();
            slice(PCtemp,bstart+reps*block+2*block+(r-1)*(3*block),block,tmp2);
            Cz.insert(Cz.end(), tmp2.begin(), tmp2.end()); tmp2.clear();
        }
        int jj = int(dt/span); jj = std::min(jj, reps-1);
        double M0 = t1 + span*jj;
        Matrix mCx(Cx.size(),1), mCy(Cy.size(),1), mCz(Cz.size(),1);
        for(int k=0;k<(int)Cx.size();++k){
            mCx(k+1,1)=Cx[k];
            mCy(k+1,1)=Cy[k];
            mCz(k+1,1)=Cz[k];
        }
        return 1e3 * Cheb3D(Mjd_TDB, block, M0, M0+span, mCx,mCy,mCz);
    };

    pos.r_Mercury = loadBody(   3,14,4, 8);
    pos.r_Venus   = loadBody( 171,10,2,16);
    pos.r_Mars    = loadBody( 309,11,1,32);
    pos.r_Jupiter = loadBody( 342, 8,1,32);
    pos.r_Saturn  = loadBody( 366, 7,1,32);
    pos.r_Uranus  = loadBody( 387, 6,1,32);
    pos.r_Neptune = loadBody( 405, 6,1,32);
    pos.r_Pluto   = loadBody( 423, 6,1,32);

    // --- Sección 12: Nutaciones y Libraciones ---
    // En JPL_Eph_DE430.cpp, dentro de la sección 12:

    // Nutaciones: offset=819, block=10, reps=4, span=8
    {
        std::vector<double> Cx, Cy, tmp;
        std::vector<double> Cz;  // zeros
        int base = 819;
        slice(PCtemp, base,     10, Cx);
        slice(PCtemp, base+10,  10, Cy);
        Cz.resize(10, 0.0);
        for (int rep = 1; rep < 4; ++rep) {
            slice(PCtemp, base + rep*20,    10, Cx);
            slice(PCtemp, base + rep*20+10, 10, Cy);
            // Cz permanece en cero
        }
        int jn = std::min(int(dt/8), 3);
        double M0n = t1 + 8*jn;
        Matrix mCx(40,1), mCy(40,1), mCz(40,1);
        for (int k = 0; k < 40; ++k) {
            mCx(k+1,1) = Cx[k];
            mCy(k+1,1) = Cy[k];
            mCz(k+1,1) = 0.0;
        }
        pos.Nutations = Cheb3D(Mjd_TDB, 10, M0n, M0n+8, mCx, mCy, mCz);
    }

    // Libraciones: offset=899, block=10, reps=4, span=8
    {
        std::vector<double> Cx, Cy, Cz, tmp;
        int base = 899;
        slice(PCtemp, base,     10, Cx);
        slice(PCtemp, base+10,  10, Cy);
        slice(PCtemp, base+20,  10, Cz);
        for (int rep = 1; rep < 4; ++rep) {
            slice(PCtemp, base + rep*30,    10, Cx);
            slice(PCtemp, base + rep*30+10, 10, Cy);
            slice(PCtemp, base + rep*30+20, 10, Cz);
        }
        int jl = std::min(int(dt/8), 3);
        double M0l = t1 + 8*jl;
        Matrix mCx(40,1), mCy(40,1), mCz(40,1);
        for (int k = 0; k < 40; ++k) {
            mCx(k+1,1) = Cx[k];
            mCy(k+1,1) = Cy[k];
            mCz(k+1,1) = Cz[k];
        }
        pos.Librations = Cheb3D(Mjd_TDB, 10, M0l, M0l+8, mCx, mCy, mCz);
    }


    // --- Correcciones finales ---
    const double EMRAT=81.30056907419062;
    pos.r_Earth = pos.r_Earth - (1.0/(1.0+EMRAT)) * pos.r_Moon;
    auto rel=[&](const Matrix& v){ return -1*pos.r_Earth + v; };
    pos.r_Mercury=rel(pos.r_Mercury);
    pos.r_Venus  =rel(pos.r_Venus);
    pos.r_Mars   =rel(pos.r_Mars);
    pos.r_Jupiter=rel(pos.r_Jupiter);
    pos.r_Saturn =rel(pos.r_Saturn);
    pos.r_Uranus =rel(pos.r_Uranus);
    pos.r_Neptune=rel(pos.r_Neptune);
    pos.r_Pluto  =rel(pos.r_Pluto);
    pos.r_Sun    =rel(pos.r_Sun);

    return pos;
}