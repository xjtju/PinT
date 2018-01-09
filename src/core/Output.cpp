#include "Output.h"


void Output::var_outer(FILE *fp, double *p) {

}

/**
 * 3D
 * it is hard to express 3D data to an ASCII file, so
 * the func outputs all cross surfaces(section) orderly from one direction.   
 * direction: 0-X (the YZ surface), 1-Y, 2-Z
 */
void Output::var_outer_X(FILE *fp, double *p){
    int ind;
    for(int i=sx-1; i>=0; i--) {
        if( (i==nguard-1) || (i==sx-nguard-1)) fprintf(fp, "  ===========================  ===========================\n");
        fprintf(fp, " x global index : %d \n", idx + i - nguard);
        for(int k=sz-1; k>=0 ; k--) { 
            if( (k==nguard-1) || (k==sz-nguard-1)) fprintf(fp, "  ----------  \n");
            for(int j = 0; j < sy ; j++){
                ind = grid->getOuterIdx(i, j, k); 
                if( (j==nguard) || (j==sy-nguard)) fprintf(fp, " | ");
                fprintf (fp, "  %10.5f  ", p[ind]);
            }
            fprintf(fp,"\n"); 
        }
        fprintf(fp,"\n"); 
    }
}

void Output::var_outer_Y(FILE *fp, double *p){
    int ind;
    for(int j=sy-1; j>=0; j--) {
        if( (j==nguard-1) || (j==sy-nguard-1)) fprintf(fp, "  ===========================  ===========================\n");
        fprintf(fp, " y global index : %d \n", idy + j - nguard);
        for(int k=sz-1; k>=0 ; k--) { 
            if( (k==nguard-1) || (k==sz-nguard-1)) fprintf(fp, "  ----------  \n");
            for(int i = 0; i < sx ; i++){
                ind = grid->getOuterIdx(i, j, k); 
                if( (i==nguard) || (i==sx-nguard)) fprintf(fp, " | ");
                fprintf (fp, "  %10.5f  ", p[ind]);
            }
            fprintf(fp,"\n"); 
        }
        fprintf(fp,"\n"); 
    }
}

void Output::var_outer_Z(FILE *fp, double *p){
    int ind;
    for(int k=sz-1; k>=0; k--) {
        if( (k==nguard-1) || (k==sz-nguard-1)) fprintf(fp, "  ===========================  ===========================\n");
        fprintf(fp, " z global index : %d \n", idz + k - nguard);
        for(int j=sy-1; j>=0 ; j--) { 
            if( (j==nguard-1) || (j==sy-nguard-1)) fprintf(fp, "  ----------  \n");
            for(int i = 0; i < sx ; i++){
                ind = grid->getOuterIdx(i, j, k); 
                if( (i==nguard) || (i==sx-nguard)) fprintf(fp, " | ");
                fprintf (fp, "  %10.5f  ", p[ind]);
            }
            fprintf(fp,"\n"); 
        }
        fprintf(fp,"\n"); 
    }
}
