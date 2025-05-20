#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./include/Sat_const.h"
#include "./include/Mjday.h"

//--------------------------------------------------------------------------
//
// Initial Orbit Determination using Gauss and Extended Kalman Filter methods
//
// References:
//   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
//   Applications", Springer Verlag, Heidelberg, 2000.
//   
//   D. Vallado, "Fundamentals of Astrodynamics and Applications", 
//   4th Edition, 2013.
//
//   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
//
// Last modified:   2020/03/16   Meysam Mahooti
//--------------------------------------------------------------------------
using namespace std;

int main(){
	//double Cnm[362][362], Snm[362][362], **eopdata;
    double  **eopdata;
    double** Cnm = (double**) malloc(362 * sizeof(double*));
    double** Snm = (double**) malloc(362 * sizeof(double*));
    for (int i = 0; i < 362; ++i) {
        Cnm[i] = (double*) malloc(362 * sizeof(double));
        Snm[i] = (double*) malloc(362 * sizeof(double));
    }
    FILE *fp;
    int f, c;
    double aux1, aux2;

	 int n_eqn;

    fp = fopen("./data/egm.txt", "r");
    if(fp == NULL){
		printf("Fail open GGM03S.txt file\n");
		exit(EXIT_FAILURE);
	}

    for(int n = 0; n <= 360; n++){
		for(int m = 0; m <= n; m++){
			fscanf(fp, "%d%d%lf%lf%lf%lf", &f, &c, &Cnm[n+1][m+1], &Snm[n+1][m+1], &aux1, &aux2);
		}
	}
	fclose(fp);
    
    cout << Cnm[45][23] << endl;
    
    //read Earth orientation parameters
    //  ----------------------------------------------------------------------------------------------------
    // |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    // |(0h UTC)           "         "          s          s          "        "          "         "     s
    //  ----------------------------------------------------------------------------------------------------

    fp = fopen("./data/eop19620101.txt", "r");
    if(fp == NULL){
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    eopdata = (double **) malloc(14 * sizeof(double *));
    if(eopdata == NULL){
        printf("eopdata: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i <= 14; i++){
        eopdata[i] = (double *) malloc(21413 * sizeof(double));
        if(eopdata[i] == NULL){
            printf("eopdata[i]: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }

    for(int i = 0; i < 21413; i++){
        fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eopdata[1][i], &eopdata[2][i], &eopdata[3][i], &eopdata[4][i], &eopdata[5][i], &eopdata[6][i], &eopdata[7][i], &eopdata[8][i], &eopdata[9][i], &eopdata[10][i], &eopdata[11][i], &eopdata[12][i], &eopdata[13][i]);
    }
    fclose(fp);

    cout << eopdata[2][12454] << endl;
    
    double obs[18][4];

    // read observations
    fp = fopen("./data/GEOS3.txt", "r");
    if(fp == NULL){
        printf("Fail open GEOS3.txt file\n");
        exit(EXIT_FAILURE);
    }

    char line[55], y[5], mo[3], d[3], h[3], mi[3], s[7], a[9], e[9], di[10];
    int Y, M, D, hh, mm, ss, az, el, Dist;
    
    for(int i = 0; i < 18; i++){
        fgets(line, sizeof(line)+2, fp);

        strncpy(y, &line[0], 4);
        y[4]='\0';
        Y=atoi(y);

        strncpy(mo, &line[5], 2);
        mo[2]='\0';
        M=atoi(mo);

        strncpy(d, &line[8],2);
        d[2]='\0';
        D=atoi(d);

        strncpy(h, &line[12], 2);
        h[2]='\0';
        hh=atoi(h);

        strncpy(mi, &line[15], 2);
        mi[2]='\0';
        mm=atoi(mi);

        strncpy(s, &line[18], 6);
        s[6]='\0';
        ss=atof(s);

        strncpy(a, &line[25], 8);
        a[8]='\0';
        az=atof(a);

        strncpy(e, &line[35], 8);
        e[8]='\0';
        el=atof(e);

        strncpy(di, &line[44], 9);
        di[9]='\0';
        Dist=atof(di);

        cout << Y << "  " << M <<  "  " << D <<  "  " << hh <<  "  " << mm <<  "  " << ss <<  "  " << az <<  "  " << el <<  "  " << Dist << endl;
      }
    fclose(fp);
    
    for(int i = 0; i <= 14; i++){
        free(eopdata[i]);
    }
    free(eopdata);
    for (int i = 0; i < 362; ++i) {
        free(Cnm[i]);
        free(Snm[i]);
    }
    free(Cnm);
    free(Snm);
    return 0;

}
	
