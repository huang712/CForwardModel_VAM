//
// Created by Feixiong Huang on 10/23/17.
//

//***************************************************************************
// This file is part of the CYGNSS E2ES.
//
// CYGNSS E2ES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CYGNSS E2ES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CYGNSS E2ES.  If not, see <http://www.gnu.org/licenses/>.
//
//---------------------------------------------------------------------------
//
// This file contains various vector, matrix, and math routines used
// throughout the E2ES
//
//****************************************************************************/

#include "gnssr.h"

double cot(double z){ return 1.0 / tan(z); }
double sec(double z){ return 1.0 / cos(z); }
double csc(double z){ return 1.0 / sin(z); }

//-----------------------------------------------------------------------------------
void bubble(int *a,int n) //array a and length n
{
    int i,j,temp;
    for(i=0;i<n-1;i++){
        for(j=i+1;j<n;j++){
            if(a[i]>a[j]) {
                temp=a[i];
                a[i]=a[j];
                a[j]=temp;
            }
        }
    }
}

void bilinear_interp(double *x_vec, double *y_vec, int size_x, int size_y, double x, double y, int *bi_index, double *bi_weight, double resolution){
    //find indices from bilinear interpolation
    //resolution = lat/lon resolution
    //x=lat, y=lon, index= lon * NLAT + lat;
    int ix1,ix2,iy1,iy2;
    double x1,x2,y1,y2;
    //find ix1 and ix2, the two index that are nearest to x

    for (int i=0; i<size_x; i++){
        if (fabs(x_vec[i]-x)<resolution){
            ix1=i; ix2=i+1;
            x1=x_vec[ix1]; x2=x_vec[ix2];
            break;
        }
    }

    for (int i=0; i<size_y; i++){
        if (fabs(y_vec[i]-y)<resolution){
            iy1=i; iy2=i+1;
            y1=y_vec[iy1]; y2=y_vec[iy2];
            break;
        }
    }
    //printf("x1 x2 y1 y2= %f %f %f %f\n",x1,x2,y1,y2);
    bi_index[0] = iy1*size_x+ix1; //index of Q11
    bi_index[1] = iy1*size_x+ix2; //index of Q21
    bi_index[2] = iy2*size_x+ix1; //index of Q12
    bi_index[3] = iy2*size_x+ix2;//index of Q22

    bi_weight[0] = (x2-x)*(y2-y)/((x2-x1)*(y2-y1)); //weight of Q11
    bi_weight[1] = (x-x1)*(y2-y)/((x2-x1)*(y2-y1)); //weight of Q21
    bi_weight[2] = (x2-x)*(y-y1)/((x2-x1)*(y2-y1)); //weight of Q12
    bi_weight[3] = (x-x1)*(y-y1)/((x2-x1)*(y2-y1)); //weight of Q22

}



int find_nearest(double *vec, int size, double value){  //by Feixiong
    //find the nearest value in vec[] and return the index
    double diff, temp_diff;
    int index;

    index = 0;
    diff = fabs(vec[0]-value);
    for (int i = 0;i < size; i++){
        temp_diff = fabs(vec[i] - value);
        if (temp_diff < diff){
            diff = temp_diff;
            index = i;
        }
    }
    return index;

}

double linear_interp( double a, double b, int direction, double time_01 ){
    if( direction == 0 )
        return (1-time_01) * a + time_01 * b;
    else
        return (1-time_01) * b + time_01 * a;
}

void cubic_interpolation( double f0, double f1, double df0, double df1,
                          double t, double *ft, double *dft, double timeInterval_s){
    // t must be between 0 and 1
    if( (t < 0) || (t > 1) ){
        printf("Time error in cubic interp!\n");
        exit(0);
    }

    // special cases when no interpolation needed
    if(t == 0){ *ft  = f0; *dft = df0; return; }
    if(t == 1){ *ft  = f1; *dft = df1; return; }

    df0 = df0 * timeInterval_s;
    df1 = df1 * timeInterval_s;

    // solve for coef of thrid deg polynom
    const double a = 2*f0 - 2*f1 + df0 + df1;
    const double b = -3*f0 + 3*f1 - 2*df0 - df1;
    const double c = df0;
    const double d = f0;

    // values of polynom and derivative at t
    *ft  = a*pow(t,3) + b*pow(t,2) + c*t + d;
    *dft = 3*a*pow(t,2) + 2*b*t + c;

    *dft = *dft * ( 1.0 / timeInterval_s );
}


void cubic_interpolation_3vector( double f0[3], double f1[3], double df0[3], double df1[3],
                                  double t, double ft[3], double dft[3], double timeInterval_s){

    // cubic interpolate a 3d vector
    cubic_interpolation(f0[0], f1[0], df0[0], df1[0], t, &(ft[0]), &(dft[0]), timeInterval_s);
    cubic_interpolation(f0[1], f1[1], df0[1], df1[1], t, &(ft[1]), &(dft[1]), timeInterval_s);
    cubic_interpolation(f0[2], f1[2], df0[2], df1[2], t, &(ft[2]), &(dft[2]), timeInterval_s);

}

int getTimeIndicies( int numTimeSamples, double timeStart_s, double timeInc_s,
                     double time_s, int *n0, int *n1 ){
    *n0 = floor((time_s - timeStart_s) / timeInc_s );
    *n1 = *n0 + 1;

    if( (*n0 < 0) || (*n0 >= numTimeSamples) )
        return 0; // time_s outside bounds
    else
        return 1;
}

//-----------------------------------------------------------------------------------

void vector_orthoNorm( double a[3], double b[3] ){
    // take vector b, remove component parallel to a
    // and then normalize length to 1
    double s = vector_dot_product(a,b);
    b[0] = b[0] - s*a[0];
    b[1] = b[1] - s*a[1];
    b[2] = b[2] - s*a[2];
    vector_unit(b,b);
}

void vector_scale( double a[3], double b[3], double s ){
    // b = a * s
    b[0] = a[0]*s;
    b[1] = a[1]*s;
    b[2] = a[2]*s;
}

double  vector_dot_product (double *a, double *b) {
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

void vector_cross_product(double a[3], double b[3],double c[3]){
    // c = a x b
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
}

double vector_norm(double v[3]) {
    return(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}

void vector_unit(double *x, double *x_unit) {
    // takes x and produces a unit vector in the
    // same direction (x and x_unit can be the same
    // when calling this function)
    double X = vector_norm(x);
    x_unit[0] = x[0] / X;
    x_unit[1] = x[1] / X;
    x_unit[2] = x[2] / X;
}

void vector_add(double a[3], double b[3],double c[3]) {
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

void vector_subtract(double a[3], double b[3],double c[3]) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

void vector_constrainToPlane( double a[3], double b[3], double c[3] ){
    double tempa[3], tempb[3];
    vector_unit(a, tempa);
    vector_scale(b, tempb, 1);
    vector_orthoNorm(tempa, tempb);
    double c1 = vector_dot_product(tempa, c);
    double c2 = vector_dot_product(tempb, c);

    c[0] = tempa[0] * c1 + tempb[0] * c2;
    c[1] = tempa[1] * c1 + tempb[1] * c2;
    c[2] = tempa[2] * c1 + tempb[2] * c2;
}

//-----------------------------------------------------------------------------------

void matrix_multiply(unsigned rowsA,unsigned columnsA,unsigned columnsB,double *A,double *B,double *C) {
    unsigned i, j, k;

    for(i = 0; i < rowsA; i++) {
        for(j = 0; j < columnsB; j++) {
            C[i*columnsB +j] = 0.0;
            for(k = 0; k < columnsA; k++) {
                C[i*columnsB +j] +=  A[i*columnsA+k]*B[k*columnsB+j];
            }
        }
    }
}

void matrix_transpose(unsigned rows,unsigned columns,double *A,double *B) {
    unsigned i, j;
    for(i = 0; i < columns; i++) {
        for(j = 0; j < rows; j++) {
            B[i*rows+j] = A[j*columns+i];
        }
    }
}

void matrixVectorMult3x3( double M[9], double x[3], double y[3] ){
    // y = M*x
    y[0] = M[0] * x[0] + M[1] * x[1] + M[2] * x[2];
    y[1] = M[3] * x[0] + M[4] * x[1] + M[5] * x[2];
    y[2] = M[6] * x[0] + M[7] * x[1] + M[8] * x[2];
}

void matrix_form3x3( double row1[3], double row2[3], double row3[3], double *M ){
    // form a matrix M by combining three row vectors
    for(int i=0;i<3;i++){
        M[i]   = row1[i];
        M[i+3] = row2[i];
        M[i+6] = row3[i];
    }
}

double matrix_det_3x3(double m[3][3]){
    return   m[0][0] * (m[1][1]*m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0]*m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0]*m[2][1] - m[1][1] * m[2][0]);
}

void matrix_scaleAdjoint_3x3(double a[3][3],double s,double m[3][3]){
    a[0][0] = (s) * (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
    a[1][0] = (s) * (m[1][2] * m[2][0] - m[1][0] * m[2][2]);
    a[2][0] = (s) * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    a[0][1] = (s) * (m[0][2] * m[2][1] - m[0][1] * m[2][2]);
    a[1][1] = (s) * (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
    a[2][1] = (s) * (m[0][1] * m[2][0] - m[0][0] * m[2][1]);
    a[0][2] = (s) * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
    a[1][2] = (s) * (m[0][2] * m[1][0] - m[0][0] * m[1][2]);
    a[2][2] = (s) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
}


void matrix_invert_3x3(double A[9], double invA[9]){

    double invM[3][3], M[3][3];
    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            M[row][col] = A[col + row*3];
        }
    }

    matrix_scaleAdjoint_3x3(invM, (1.0 / matrix_det_3x3(M)), M);

    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            invA[col + row*3] = invM[row][col];
        }
    }

}

/*
void testBilinear(void){
    GLuint texture;

    int width, height;
    char * data;
    FILE * file;
    int wrap = 0;

    // allocate buffer
    width = 256;
    height = 256;
    data = malloc( width * height * 3 );


    // allocate a texture name
    glGenTextures( 1, &texture );

    // select our current texture
    glBindTexture( GL_TEXTURE_2D, texture );

    // select modulate to mix texture with color for shading
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );

    // when texture area is small, bilinear filter the closest MIP map
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_NEAREST );
    // when texture area is large, bilinear filter the first MIP map
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    // if wrap is true, the texture wraps over at the edges (repeat)
    //       ... false, the texture ends at the edges (clamp)
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    wrap ? GL_REPEAT : GL_CLAMP );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
                    wrap ? GL_REPEAT : GL_CLAMP );

    // build our texture MIP maps
    gluBuild2DMipmaps( GL_TEXTURE_2D, 3, width,
                      height, GL_RGB, GL_UNSIGNED_BYTE, data );

    // free buffer
    free( data );

    //return texture;
}
*/

