//Functions created by
// G V V Sharma
// October 27, 2023
// Revised November 13, 2023

//Function declaration

double **Mateigval(double **a);//eigenvalues of a 2x2 matrix
double **Matquad(double a,double b, double c);//roots of a quadratic equation
double **Mateye(int m);//generate identity matrix
double **createMat(int m,int n);//create m x n matrix array
void printMat(double **p,int m,int n);//print matrix
double **loadMat(char *str,int m,int n);//load matrix from file
double Matnorm(double **a, int m);//norm of a vector
double Matdot(double **a, double ** b, int m);//inner product
double **Matsub(double **a, double **b, int m, int n);//subtract two matrices
double **Matadd(double **a, double **b, int m, int n);//add two matrices
double **Matscale(double **a, int m, int n, double k);//scale matrix
double **Matinv(double **mat, int m);//invert an m  x m matrix, m <=3
double **Matmul(double **a, double **b, int m, int n, int p);//multiply matrices a and b
double **Mathstack(double **a, double **b, int m, int n, int p);//horizontall stack matrices a and b of size mxn and mxp respectively
double **transposeMat(double **a,  int m, int n);//transpose of a
double **rotMat(double theta); //rotation matrix
double **normVec(double **a); //normal vector
double **circulantMat(double **a, int m);//Generating a circulant matrix from a vector
double **Matsec(double **a, double ** b, int m, double k);//section formula
double Matrace(double **a, int m);//trace of a matrix
double Matdet(double **a);//determinant of a 2x2 matrix
double **Matcol(double **a,int m, int n);//Extract nth column
double **Matrow(double **a,int m, int n);//Extract mth row
double **Matunit(double **a,int m);//Generate unit vector
double **Mateigvec(double **a, double **lam);
//End function declaration

//section formula
double **Matsec(double **a, double ** b, int m, double k){
	double **temp=createMat(m,1);
	temp = Matscale(Matadd(a,Matscale(b,m,1,k),m,1),m,1,1/(k+1));
	return temp;
}

//add matrices
double **Matadd(double **a,double **b, int m, int n){
int i, j;
double **c;
c = createMat(m,n);

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
c[i][j]= a[i][j]+b[i][j];
  }
 }
return c;
}

//scale matrix
double **Matscale(double **a, int m, int n, double k){
int i, j;
double **c;
c = createMat(m,n);

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
c[i][j]= k*a[i][j];
  }
 }
return c;
}

//Obtaining the normal vector
double **normVec(double **m){
	double **temp;
	temp = Matmul(rotMat(M_PI/2),m,2,2,1);
	return temp;
}

//Defining the function for matrix creation
double **createMat(int m,int n)
{
 int i;
 double **a;
 
 //Allocate memory to the pointer
a = (double **)malloc(m * sizeof( *a));
    for (i=0; i<m; i++)
         a[i] = (double *)malloc(n * sizeof( *a[i]));

 return a;
}
//End function for matrix creation

//Extract column
//
double **Matcol(double **a,int m, int n){
	int i = 0;
	double **b = createMat(m,1);//create column with m rows

//extract column vector
	for (i = 0; i < m; i++){
		b[i][0] = a[i][n];
	}
return b;
}

//Extract row
//
double **Matrow(double **a,int m, int n){
	int i = 0;
	double **b = createMat(n,1);//create column with n rows

//extract column vector
	for (i = 0; i < n; i++){
		b[i][0] = a[m][i];
	}
return b;
}


//Read  matrix from file
double **loadMat(char *str,int m,int n)
{
FILE *fp;
double **a;
int i,j;


a = createMat(m,n);
fp = fopen(str, "r");

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
   fscanf(fp,"%lf",&a[i][j]);
  }
 }
//End function for reading matrix from file

fclose(fp);
 return a;

}


//Defining the function for printing
void printMat(double **p, int m,int n)
{
 int i,j;

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  printf("%lf ",p[i][j]);
 printf("\n");
 }
}
//End function for printing

//Rotation matrix

double **rotMat(double theta){
double **temp=createMat(2,2);//creating the matrix
double c = cos(theta), s = sin(theta);
temp[0][0] = c;
temp[0][1] = -s;
temp[1][0] = s;
temp[1][1] = c;

return temp;
}
//inner product
double Matdot(double **a, double ** b, int m){
	double **temp= Matmul(transposeMat(a,m,1),b,1,m,1);
	return temp[0][0];
}
//Defining the function for norm

double Matnorm(double **a, int m){
	return sqrt(Matdot(a, a, m));
}
//Defining the function for difference of matrices

double **Matsub(double **a, double **b, int m, int n)
{
int i, j;
double **c;
c = createMat(m,n);

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
c[i][j]= a[i][j]-b[i][j];
  }
 }
return c;

}
//End function for difference of matrices

//Defining the function for inverse of a matrix
//code adapted from the internet


double **Matinv(double **a, int m)
{
double **c, det=0;
int i,j;
c = createMat(m,m);
printMat(c,m,m);
if (m==2){
det = a[0][0]*a[1][1]-a[0][1]*a[1][0];

c[0][0] = a[1][1]/det;
c[0][1] = -a[1][0]/det;
c[1][0] = -a[0][1]/det;
c[1][1] = a[0][0]/det;
}
else if(m==3){
for(i=0;i<m;i++)
      det += a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]);
 
   for(i=0;i<m;i++){
      for(j=0;j<m;j++)
	   c[i][j]=((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/det;
   }
}
else {
	printf("Invalid input \n");
	exit(0);
}
 
return transposeMat(c,m,m);//revised based on abhimany koushik's suggestion.  
//return c;
}
// End  function for inverse of 2x2 matrix


//Defining the function for product of matrices

double **Matmul(double **a, double **b, int m, int n, int p)
{
int i, j, k;
double **c, temp =0;
c = createMat(m,p);

 for(i=0;i<m;i++)
 {
  for(k=0;k<p;k++)
  {
    for(j=0;j<n;j++)
    {
	temp= temp+a[i][j]*b[j][k];
    }
	c[i][k]=temp;
	temp = 0;
  }
 }
return c;

}
//End function for difference of matrices

//Defining the function for transpose of matrix

double **transposeMat(double **a,  int m, int n)
{
int i, j;
double **c;
//printf("I am here");
c = createMat(n,m);

 for(i=0;i<n;i++)
 {
  for(j=0;j<m;j++)
  {
c[i][j]= a[j][i];
//  printf("%lf ",c[i][j]);
  }
 }
return c;

}
//End function for transpose of matrix
//trace of a matrix
double Matrace(double **a, int m){
	double c=0;
 for(int i=0;i<m;i++){
	 c +=  a[i][i];
 }
 return c;
}
//determinant of a 2x2 matrix
double Matdet(double **a){
return a[0][0]*a[1][1]-a[0][1]*a[1][0];
}
//eigenvalues of a 2x2 matrix
double **Mateigval(double **a){
	double b = -Matrace(a,2);
	double c = Matdet(a);
	return Matquad(1,b,c);
}

//generate identity matrix
double **Mateye(int m){

	int i,j;//dummy integers
double **I = createMat(m,m);//create m x n matrix array
 for(int i=0;i<m;i++){
 	for(int j=0;j<m;j++){
		if(i==j)
			 I[i][j] = 1;
		else
			 I[i][j] = 0;
 }
 }
return I;
}
//Generate unit vector
double **Matunit(double **a,int m){
	double **c=createMat(2,1);
	double temp = Matnorm(a,m);
	for(int i = 0; i < m; i++){
		c[i][0] = a[i][0]/temp; 
	}
	return c;
}
//roots of a quadratic equation
double **Matquad(double a,double b, double c){
	double **lam = createMat(2,1);
	double D = sqrt(pow(b,2.0)-4*a*c);
	double den =2.0*a; 
lam[0][0] = (-b+D)/den;
lam[1][0] = (-b-D)/den;
return lam;
}
//eigenvector matrix for a 2x2 matrix
//double **Mateigvec(double **a){
double **Mateigvec(double **a, double **lam){
	double **b1, **b2;
	double **p1, **p2;
	double **temp1, **temp2;
	double **omat=rotMat(M_PI/2);
	//A-lambda I
	b1 = Matadd(a,Matscale(Mateye(2),2,2,-lam[0][0]),2,2);
	b2 = Matadd(a,Matscale(Mateye(2),2,2,-lam[1][0]),2,2);
	//Extract 1st row
	temp1 = Matrow(b1, 0, 2);
	temp2 = Matrow(b2, 0, 2);
	//free the matrices
	free(b1);
	free(b2);
	//Generate unit vector
	b1 = Matunit(temp1, 2);
	b2 = Matunit(temp2, 2);
	//free temp vectors
	free(temp1);
	free(temp2);
	//Find eigen vector
	p1 = Matmul(omat, b1, 2, 2, 1);
	p2 = Matmul(omat, b2, 2, 2, 1);
	//free vectors
	free(b1);
	free(b2);
	return Mathstack(p1, p2, 2, 1, 1);
}
//horizontall stack matrices a and b of size mxn and mxp respectively
double **Mathstack(double **a, double **b, int m, int n, int p){
	double **c = createMat(m,n+p);//create m x n+p matrix array
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			c[i][j] = a[i][j];//a values in c
		for(j = n; j < n+p; j++){
			c[i][j] = b[i][j-n];//b values in c
		}
		}
	}
		return c;
}


//Generating a circulant matrix from a vector
double **circulantMat(double **a, int m){
    int i,j,k;
    printMat(a,4,1);
double **c = createMat(m,m);
    for(i=0; i < m; i++){
    	for(j=0; j < m; j++){
		if (i >=j)
	    c[i][j] = a[i-j][0];
		else
	    c[i][j] = a[m-j+i][0];
    }
    }
return c;
}
