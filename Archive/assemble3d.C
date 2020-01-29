#include <stdio.h>
#include <stdlib.h>

typedef int INT;

/*  0<=i1,i2,i3<=m-1, 0<=j1,j2,j3<=n-1  without specifying

    (i1*m+i2)*m+i3
    (i1,  i2,  i3  ) (j1, j2, j3 )

    ((i3-1)*m+i1)*m+i2,              i3>0
    (i1,  i2,  i3  ) (j1, j2, n  )
    (i1,  i2,  i3-1) (j1, j2, n+1)

    ((i2-1)*m+i1)*m+i3,              i2>0
    (i1,  i2-1,i3  ) (j1, n+1,j3 )
    (i1,  i2,  i3  ) (j1, n,  j3 )

    ((i1-1)*m+i2)*m+i3,              i1>0
    (i1-1,i2,  i3  ) (n+1,j2, j3 )
    (i1,  i2,  i3  ) (n,  j2, j3 )

    ((i2-1)*(m-1)+i3-1)*m+i1,        i2,i3>0
    (i1,  i2-1,i3-1) (j1, n+1,n+1)
    (i1,  i2-1,i3  ) (j1, n+1,n  )
    (i1,  i2,  i3-1) (j1, n,  n+1)
    (i1,  i2,  i3  ) (j1, n,  n  )

    ((i1-1)*(m-1)+i3-1)*m+i2,        i1,i3>0
    (i1,  i2,  i3  ) (n,  j2, n  )
    (i1,  i2,  i3-1) (n,  j2, n+1)
    (i1-1,i2,  i3  ) (n+1,j2, n  )
    (i1-1,i2,  i3-1) (n+1,j2, n+1)

    ((i1-1)*(m-1)+i2-1)*m+i3,        i1,i2>0
    (i1-1,i2-1,i3  ) (n+1,n+1,j3 )
    (i1-1,i2,  i3  ) (n+1,n,  j3 )
    (i1,  i2-1,i3  ) (n,  n+1,j3 )
    (i1,  i2,  i3  ) (n,  n,  j3 )

    ((i1-1)*(m-1)+i2-1)*(m-1)+i3-1,  i1,i2,i3>0
    (i1-1,i2-1,i3-1) (n+1,n+1,n+1)
    (i1-1,i2-1,i3  ) (n+1,n+1,n  )
    (i1-1,i2,  i3-1) (n+1,n,  n+1)
    (i1-1,i2,  i3  ) (n+1,n,  n  )
    (i1,  i2-1,i3-1) (n,  n+1,n+1)
    (i1,  i2-1,i3  ) (n,  n+1,n  )
    (i1,  i2,  i3-1) (n,  n,  n+1)
    (i1,  i2,  i3  ) (n,  n,  n  )

*/

INT l2gindex(INT m, INT i1, INT i2, INT i3, INT n, INT j1, INT j2, INT j3)
{	     
 INT  off[8], p=-1, m1=m-1;

 n--;

 off[0] = 0;
 off[1] = off[0] + m*m*m*n*n*n;
 off[2] = off[1] + m1*m*m*n*n;
 off[3] = off[2] + m1*m*m*n*n;
 off[4] = off[3] + m1*m*m*n*n;
 off[5] = off[4] + m1*m1*m*n;
 off[6] = off[5] + m1*m1*m*n;
 off[7] = off[6] + m1*m1*m*n;
 
/* local-to-gobal mapping*/
 if ( j1<n ) {
 	 	 	
 	 if ( j2<n ) {
 	 	 	
 	 	 if ( j3<n ) { 
	 	 	 
 	 	 	 p = off[0] + ((((i1*m+i2)*m+i3)*n+j1)*n+j2)*n+j3;    
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[1] + ((((i3-1)*m+i1)*m+i2)*n+j1)*n+j2;
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /* right */
 	 	 	 
 	 	 	 p = off[1] + (((i3*m+i1)*m+i2)*n+j1)*n+j2;
 	 	 	 
 	 	 }
 	 } else if ( j2==n && i2>0 ) { /* left */
    	 
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[2] + ((((i2-1)*m+i1)*m+i3)*n+j1)*n+j3;
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[4] + (((i2-1)*m1+(i3-1))*m+i1)*n+j1;
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[4] + (((i2-1)*m1+i3)*m+i1)*n+j1;
 	 	 	 
 	 	 }
 	 } else if ( j2>n && i2<m1 ) { /* left */
 	 	 	
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[2] + (((i2*m+i1)*m+i3)*n+j1)*n+j3;
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[4] + ((i2*m1+(i3-1))*m+i1)*n+j1;  
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[4] + ((i2*m1+i3)*m+i1)*n+j1;
 	 	 	 
 	 	 }
 	 }
 } else if ( j1==n && i1>0 ) { /* left */
 	 	
 	 if ( j2<n ) { 
 	 	 	
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[3] + ((((i1-1)*m+i2)*m+i3)*n+j2)*n+j3;
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[5] + (((i1-1)*m1+(i3-1))*m+i2)*n+j2; 

 	 	 } else if ( j3>n && i3<m1 ) {  /* right */
 	 	 	 
 	 	 	 p = off[5] + (((i1-1)*m1+i3)*m+i2)*n+j2; 
 	 	 	  
 	 	 }   
 	 } else if ( j2==n && i2>0 ) { /* left */

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + (((i1-1)*m1+(i2-1))*m+i3)*n+j3;
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+(i2-1))*m1+(i3-1);
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+(i2-1))*m1+i3;
 	 	 	 
 	 	 }
 	 } else if ( j2>n && i2<m1 ) { /* left */

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + (((i1-1)*m1+i2)*m+i3)*n+j3;
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+i2)*m1+(i3-1);
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+i2)*m1+i3;
 	 	 	 
 	 	 }   
 	 }
 } else if ( j1>n && i1<m1 ) {
 	 	
 	 if ( j2<n ) {
 	 	 	
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[3] + (((i1*m+i2)*m+i3)*n+j2)*n+j3;
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[5] + ((i1*m1+(i3-1))*m+i2)*n+j2;
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /* right */
 	 	 	 
 	 	 	 p = off[5] + ((i1*m1+i3)*m+i2)*n+j2;
 	 	 	 
 	 	 }   
 	 } else if ( j2==n && i2>0 ) { /* left */

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + ((i1*m1+(i2-1))*m+i3)*n+j3;  
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + (i1*m1+(i2-1))*m1+(i3-1);
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	
 	 	 	 p = off[7] + (i1*m1+(i2-1))*m1+i3;
 	 	 	 
 	 	 }   
 	 } else if ( j2>n && i2<m1 ) { /* left */

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + ((i1*m1+i2)*m+i3)*n+j3;
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + (i1*m1+i2)*m1+(i3-1);
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[7] + (i1*m1+i2)*m1+i3;
 	 	 	 
 	 	 }
 	 }
 }     
 
 return p;
}

INT l2gindeg(INT m, INT i1, INT i2, INT i3, INT n, INT j1, INT j2, INT j3, INT *deg)
{	     
 INT  off[8], nh, n1, x, y, z=0, p=-1, m1=m-1;

 n--;
 n1 = n-1;
 nh = n/2-1;
 	
 off[0] = 0;
 off[1] = off[0] + m*m*m*n*n*n;
 off[2] = off[1] + m1*m*m*n*n;
 off[3] = off[2] + m1*m*m*n*n;
 off[4] = off[3] + m1*m*m*n*n;
 off[5] = off[4] + m1*m1*m*n;
 off[6] = off[5] + m1*m1*m*n;
 off[7] = off[6] + m1*m1*m*n;
 
/* local-to-gobal mapping*/
 if ( j1<n ) {
 	 x = (j1==nh||j1==n1)?(2+(i1==0?0:1)+(i1==m1?0:1)):3;
 	 	 	
 	 if ( j2<n ) {
 	 	 y = x * ((j2==nh||j2==n1)?(2+(i2==0?0:1)+(i2==m1?0:1)):3);
 	 	 	
 	 	 if ( j3<n ) { 
	 	 	 
 	 	 	 p = off[0] + ((((i1*m+i2)*m+i3)*n+j1)*n+j2)*n+j3;    
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[1] + ((((i3-1)*m+i1)*m+i2)*n+j1)*n+j2;
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /* right */
 	 	 	 
 	 	 	 p = off[1] + (((i3*m+i1)*m+i2)*n+j1)*n+j2;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }
 	 } else if ( j2==n && i2>0 ) { /* left */
    	 y = x * (5+(i2==1?0:1)+(i2==m1?0:1));
    	 
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[2] + ((((i2-1)*m+i1)*m+i3)*n+j1)*n+j3;
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[4] + (((i2-1)*m1+(i3-1))*m+i1)*n+j1;
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[4] + (((i2-1)*m1+i3)*m+i1)*n+j1;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }
 	 } else if ( j2>n && i2<m1 ) { /* left */
 	 	 y = x * (5+(i2==0?0:1)+(i2==m1-1?0:1));
 	 	 	
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[2] + (((i2*m+i1)*m+i3)*n+j1)*n+j3;
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[4] + ((i2*m1+(i3-1))*m+i1)*n+j1;  
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[4] + ((i2*m1+i3)*m+i1)*n+j1;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }
 	 }
 } else if ( j1==n && i1>0 ) { /* left */
 	 x = (5+(i1==1?0:1)+(i1==m1?0:1));
 	 	
 	 if ( j2<n ) { 
 	 	 y = x * ((j2==nh||j2==n1)?(2+(i2==0?0:1)+(i2==m1?0:1)):3);
 	 	 	
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[3] + ((((i1-1)*m+i2)*m+i3)*n+j2)*n+j3;
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[5] + (((i1-1)*m1+(i3-1))*m+i2)*n+j2; 
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /* right */
 	 	 	 
 	 	 	 p = off[5] + (((i1-1)*m1+i3)*m+i2)*n+j2;  
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }   
 	 } else if ( j2==n && i2>0 ) { /* left */
    	 y = x * (5+(i2==1?0:1)+(i2==m1?0:1));

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + (((i1-1)*m1+(i2-1))*m+i3)*n+j3;
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+(i2-1))*m1+(i3-1);
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));;
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+(i2-1))*m1+i3;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }
 	 } else if ( j2>n && i2<m1 ) { /* left */
 	 	 y = x * (5+(i2==0?0:1)+(i2==m1-1?0:1));

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + (((i1-1)*m1+i2)*m+i3)*n+j3;
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+i2)*m1+(i3-1);
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[7] + ((i1-1)*m1+i2)*m1+i3;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }   
 	 }
 } else if ( j1>n && i1<m1 ) {
 	 x = (5+(i1==0?0:1)+(i1==m1-1?0:1));
 	 	
 	 if ( j2<n ) {
 	 	 y = x * ((j2==nh||j2==n1)?(2+(i2==0?0:1)+(i2==m1?0:1)):3); 
 	 	 	
 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[3] + (((i1*m+i2)*m+i3)*n+j2)*n+j3;
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[5] + ((i1*m1+(i3-1))*m+i2)*n+j2;
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /* right */
 	 	 	 
 	 	 	 p = off[5] + ((i1*m1+i3)*m+i2)*n+j2;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }   
 	 } else if ( j2==n && i2>0 ) { /* left */
    	 y = x * (5+(i2==1?0:1)+(i2==m1?0:1));

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + ((i1*m1+(i2-1))*m+i3)*n+j3;  
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + (i1*m1+(i2-1))*m1+(i3-1);
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	
 	 	 	 p = off[7] + (i1*m1+(i2-1))*m1+i3;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }   
 	 } else if ( j2>n && i2<m1 ) { /* left */
 	 	 y = x * (5+(i2==0?0:1)+(i2==m1-1?0:1));

 	 	 if ( j3<n ) {
 	 	 	 
 	 	 	 p = off[6] + ((i1*m1+i2)*m+i3)*n+j3;
 	 	 	 
 	 	 	 z = y * ((j3==nh||j3==n1)?(2+(i3==0?0:1)+(i3==m1?0:1)):3);
 	 	 	 
 	 	 } else if ( j3==n && i3>0 ) { /* left */
 	 	 	 
 	 	 	 p = off[7] + (i1*m1+i2)*m1+(i3-1);
 	 	 	 
 	 	 	 z = y * (5+(i3==1?0:1)+(i3==m1?0:1));
 	 	 	 
 	 	 } else if ( j3>n && i3<m1 ) {  /*right*/
 	 	 	 
 	 	 	 p = off[7] + (i1*m1+i2)*m1+i3;
 	 	 	 
 	 	 	 z = y * (5+(i3==0?0:1)+(i3==m1-1?0:1));
 	 	 }
 	 }
 }     
 
 *deg = z;
 return p;
}


int main(int argc, char *argv[])
{ INT m=3, n=3;
  INT i1, i2, i3, j1, j2, j3, g, c, *nz;

 if (argc>1) m = atoi(argv[1]);
 if (argc>2) n = atoi(argv[2]);
 nz = (INT *) calloc((m*n-1)*(m*n-1)*(m*n-1),sizeof(INT));
	
 for ( i1=0; i1<m; i1++ )
 for ( i2=0; i2<m; i2++ ) 
 for ( i3=0; i3<m; i3++ ) 
     for ( j1=0; j1<=n; j1++ )	
     for ( j2=0; j2<=n; j2++ )	
     for ( j3=0; j3<=n; j3++ ) {	
 	
         // g=l2gindex(m, i1, i2, i3, n, j1, j2, j3);
         g=l2gindeg(m, i1, i2, i3, n, j1, j2, j3, &c);
         
         if (g<0) continue;
         if (nz[g] && nz[g]-c ) 
         	printf(" [%d %d %d], [%d %d %d], [%d %d %d]\n",i1,i2,i3,j1,j2,j3,g,nz[g],c);
         nz[g] = c;
         
     }
     
  c = m*n-1;
  c = c*c*c;
  for ( g=0; g<c; g++) {
    if (!nz[g])	printf("%d\n",g);
  }  
  
  
 p = 0;
 for ( i1=0; i1<m; i1++ )
 for ( i2=0; i2<m; i2++ ) 
 for ( i3=0; i3<m; i3++ ) {
 	
 	/************************************/
     for ( j1=0; j1<n-1; j1++ )	
     for ( j2=0; j2<n-1; j2++ )	
     for ( j3=0; j3<n-1; j3++ ) {	
 	
         g = l2gindex(m, i1, i2, i3, n, j1, j2, j3);
         if (g<0) continue;
         
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ )
         for ( k2=ib[j2][0]; k2<ib[j2][1]; k2++ ) 
         for ( k3=ib[j3][0]; k3<ib[j3][1]; k3++ )  {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         
         for ( k3=ib[j3][1]; k3<ib[j3+1]; k3++ ) {
         if ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) ) 
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ )
         for ( k2=ib[j2][0]; k2<ib[j2][1]; k2++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         

         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ ) 
         if ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ )
         for ( k3=ib[j3][0]; k2<ib[j3][1]; k3++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
                  
         
         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ ) 
         if ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) )
         for ( k2=ib[j2][0]; k2<ib[j2][1]; k2++ )
         for ( k3=ib[j3][0]; k2<ib[j3][1]; k3++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         

         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ )
         if  ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k3=ib[j3][1]; k2<ib[j3+1]; k3++ ) 
         if  ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) )
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         
         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ )
         if  ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) ) 
         for ( k3=ib[j3][1]; k2<ib[j3+1]; k3++ )
         if  ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) ) 
         for ( k2=ib[j2][0]; k2<ib[j2][0]; k2++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }         

         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ )
         if  ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) ) 
         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ ) 
         if  ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k3=ib[j3][0]; k2<ib[j3][1]; k3++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }   
         
         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ ) 
         if  ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) )
         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ ) 
         if  ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k3=ib[j3][1]; k2<ib[j3+1]; k3++ )
         if  ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }   
         
     }
     
 	/************************************/
     j3 = n-1;
     if  ( i3>0 ) /* (i3, n-1) and (i3-1,n) */
     for ( j1=0; j1<n-1; j1++ )	
     for ( j2=0; j2<n-1; j2++ ) {	
 	
         g = l2gindex(m, i1, i2, i3, n, j1, j2, j3);
         if (g<0) continue;

         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ )
         for ( k2=ib[j2][0]; k2<ib[j2][1]; k2++ ) 
         for ( k3=ib[j3][0]; k3<ib[j3][1]; k3++ )  {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }         	
         
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ )
         for ( k2=ib[j2][0]; k2<ib[j2][1]; k2++ ) 
         for ( k3=ib[j3][0]; k3<ib[j3][1]; k3++ )  {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         
         for ( k3=ib[j3][1]; k3<ib[j3+1]; k3++ ) {
         if ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) ) 
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ )
         for ( k2=ib[j2][0]; k2<ib[j2][1]; k2++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         

         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ ) 
         if ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ )
         for ( k3=ib[j3][0]; k2<ib[j3][1]; k3++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
                  
         
         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ ) 
         if ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) )
         for ( k2=ib[j2][0]; k2<ib[j2][1]; k2++ )
         for ( k3=ib[j3][0]; k2<ib[j3][1]; k3++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         

         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ )
         if  ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k3=ib[j3][1]; k2<ib[j3+1]; k3++ ) 
         if  ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) )
         for ( k1=ib[j1][0]; k1<ib[j1][1]; k1++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }
         
         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ )
         if  ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) ) 
         for ( k3=ib[j3][1]; k2<ib[j3+1]; k3++ )
         if  ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) ) 
         for ( k2=ib[j2][0]; k2<ib[j2][0]; k2++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }         

         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ )
         if  ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) ) 
         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ ) 
         if  ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k3=ib[j3][0]; k2<ib[j3][1]; k3++ ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }   
         
         for ( k1=ib[j1][1]; k1<ib[j1+1]; k1++ ) 
         if  ( (jb[k1]==n-1&&i1>0) || (jb[k1]==n&&i1<m-1) )
         for ( k2=ib[j2][1]; k2<ib[j2+1]; k2++ ) 
         if  ( (jb[k2]==n2-1&&i2>0) || (jb[k2]==n2&&i2<m2-1) )
         for ( k3=ib[j3][1]; k2<ib[j3+1]; k3++ )
         if  ( (jb[k3]==n-1&&i3>0) || (jb[k3]==n&&i3<m-1) ) {
         	
         	 c = l2gindex(m, i1, i2, i3, n, jb[k1], jb[k2], jb[k3]);
         	 if (c<0) continue;
         	 
         	 js[p] = c;	
         	 vs[p] = vb[k1]*vb[k2]*vb[k3];
         	 p++;
         }   
         
     }
              
     }
     
     
        
        
         	
         	
         
     }     

 
 return 0;    
              
}






 n--;
 n1 = n-1;

 

 for ( j1=0; j1<n; j1++ ) {
     if ( j1 && n/2-j1  ) { /* j1-1 */
 		
 		x = b[n+j1];
 		
        for ( j2=0; j2<n; j2++ ) {
        	r = j1*n+j2;
        	
            if ( j2 && n/2-j2  ) { /* j2-1 */
        	
        	y = x*b[n+j2];
        	
        	c = (j1-1)*n+(j2-1);
            }
            
            y = x*b[j2];  /* j2 */
            
            c = (j1-1)*n+j2;
            
            if ( n/2-1-j2 && n-1-j2 ) {
        	y = x*b[n+j2+1];
        	
        	c = (j1-1)*n+(j2-1);
            }
        }
        
            
            
            
            
        	
        	
 	    
 		
 		
 		
 		
 	}
 	

void extend(INT n, INT *ia, INT *ja, INT **ie)
{
	INT i, k;
	
	*ie = (INT *) malloc(sizeof(INT)*3*n); 
	
	for ( i=0; i<n; i++ ) {
		  (*ie)[3*i] = ia[i];
		   
	    for ( k=ia[i]; k<ia[i+1] && k<n-2; k++ );
	    
	    (*ie)[3*i+1] = k;
	    
	    for ( ; k<ia[i+1] && k<n-1; k++ );
	    
	    (*ie)[3*i+2] = k;
	}
}
 	

/*  0<=i1,i2,i3<=m-1, 0<=j1,j2,j3<=n-1  without specifying

    (i1*m+i2)*m+i3,                (j1*n+j2)*n+j3
    (i1,  i2,  i3  ) (j1, j2, j3 )

    ((i3-1)*m+i1)*m+i2,            j1*n+j2           i3>0
    (i1,  i2,  i3  ) (j1, j2, n  )
    (i1,  i2,  i3-1) (j1, j2, n+1)

    ((i2-1)*m+i1)*m+i3,            j1*n+j3           i2>0
    (i1,  i2-1,i3  ) (j1, n+1,j3 )
    (i1,  i2,  i3  ) (j1, n,  j3 )

    ((i1-1)*m+i2)*m+i3,            j2*n+j3           i1>0
    (i1-1,i2,  i3  ) (n+1,j2, j3 )
    (i1,  i2,  i3  ) (n,  j2, j3 )

    ((i2-1)*(m-1)+i3-1)*m+i1,       j1               i2,i3>0
    (i1,  i2-1,i3-1) (j1, n+1,n+1)
    (i1,  i2-1,i3  ) (j1, n+1,n  )
    (i1,  i2,  i3-1) (j1, n,  n+1)
    (i1,  i2,  i3  ) (j1, n,  n  )

    ((i1-1)*(m-1)+i3-1)*m+i2,       j2               i1,i3>0
    (i1,  i2,  i3  ) (n,  j2, n  )
    (i1,  i2,  i3-1) (n,  j2, n+1)
    (i1-1,i2,  i3  ) (n+1,j2, n  )
    (i1-1,i2,  i3-1) (n+1,j2, n+1)

    ((i1-1)*(m-1)+i2-1)*m+i3,       j3               i1,i2>0
    (i1-1,i2-1,i3  ) (n+1,n+1,j3 )
    (i1-1,i2,  i3  ) (n+1,n,  j3 )
    (i1,  i2-1,i3  ) (n,  n+1,j3 )
    (i1,  i2,  i3  ) (n,  n,  j3 )

    ((i1-1)*(m-1)+i2-1)*(m-1)+i3-1,  0               i1,i2,i3>0
    (i1-1,i2-1,i3-1) (n+1,n+1,n+1)
    (i1-1,i2-1,i3  ) (n+1,n+1,n  )
    (i1-1,i2,  i3-1) (n+1,n,  n+1)
    (i1-1,i2,  i3  ) (n+1,n,  n  )
    (i1,  i2-1,i3-1) (n,  n+1,n+1)
    (i1,  i2-1,i3  ) (n,  n+1,n  )
    (i1,  i2,  i3-1) (n,  n,  n+1)
    (i1,  i2,  i3  ) (n,  n,  n  )

*/
void assemble3d(INT n, INT j1, INT j2, int c, int a)
{
		
  n1--;
  n2--;
  n3--;
  p = q = 0;
  
  /* (i1,i2,i3), (j1,j2,j3) */
  for ( j1=0; j1<n1; j1++ ) 
  for ( j2=0; j2<n2; j2++ ) 
  for ( j3=0; j3<n3; j3++ ) {
  	
  	  /* body (i1,i2,i3), (j1,j2,j3) */
  	  i  = (i1*m2+i2)*m3+i3; 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ )
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* faces, (i1,i2,i3), (j1,j2,n3) */
      i  = ((i3-1)*m1+i1)*m2+i2;
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3 )
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* faces, (i1,i2,i3), (j1,j2,n3+1) */
      i  = (i3*m1+i1)*m2+i2;
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* faces, (i1,i2,i3), (j1,n2,j3) */
      i  = ((i2-1)*m1+i1)*m3+i3;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n3 )      	
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* faces, (i1,i2,i3), (j1,n2+1,j3) */
      i  = (i2*m1+i1)*m3+i3;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* faces, (i1,i2,i3), (n1,j2,j3) */
      i  = ((i1-1)*m2+i2)*m3+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
            	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }  
                
      /* faces, (i1,i2,i3), (n1+1,j2,j3) */
      i  = (i1*m2+i2)*m3+i3;          
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;      
      }  
      
      /* edges, (i1,i2,i3), (j1,n2,n3) */
      i  = ((i2-1)*(m3-1)+(i3-1))*m1+i1;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
      	
          jv[q] = off[4] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* edges, (i1,i2,i3), (j1,n2,n3+1) */
      i  = ((i2-1)*(m3-1)+i3)*m1+i1;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
      	
          jv[q] = off[4] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* edges, (i1,i2,i3), (j1,n+1,n) */
      i  = (i2*(m3-1)+(i3-1))*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )      	          	                    	
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
      	
          jv[q] = off[4] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* edges, (i1,i2,i3), (j1,n2+1,n3+1) */
      i  = (i2*(m3-1)+i3)*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )      	          	                    	
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[4] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
             
      /* edges, (i1,i2,i3), (n1,j2,n3) */
      i  = ((i1-1)*(m3-1)+(i3-1))*m2+i2;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i3>0   && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )	
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
      	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* edges, (i1,i2,i3), (n1,j2,n3+1) */
      i  = ((i1-1)*(m3-1)+i3)*m2+i2;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )	
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
      	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* edges, (i1,i2,i3), (n1+1,j2,n3) */
      i  = (i1*(m3-1)+(i3-1))*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )	
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
      	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* edges, (i1,i2,i3), (n1+1,j2,n3+1) */
      i  = (i1*(m3-1)+i3)*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )	
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
      	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* edges, (i1,i2,i3), (n1,n2,j3) */
      i  = ((i1-1)*(m2-1)+(i2-1))*m3+i3;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )	
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* edges, (i1,i2,i3), (n1,n2+1,j3) */
      i  = ((i1-1)*(m2-1)+i2)*m3+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )	
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* edges, (i1,i2,i3), (n1+1,n2,j3) */
      i  = (i1*(m2-1)+(i2-1))*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )	
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* edges, (i1,i2,i3), (n1+1,n2+1,j3) */
      i  = (i1*(m2-1)+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )	
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* vertices (i1,i2,i3), (n1,n2,n3) */
      i  = ((i1-1)*(m2-1)+(i2-1))*(m3-1)+(i3-1);
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   )	
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )	
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* vertices (i1,i2,i3), (n,n,n+1) */
      i  = ((i1-1)*(m2-1)+(i2-1))*(m3-1)+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   )	
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )	
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* vertices (i1,i2,i3), (n,n+1,n) */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+(i3-1);
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   )	
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )	
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* vertices (i1,i2,i3), (n,n+1,n+1) */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   )	
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )	
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* vertices (i1,i2,i3), (n+1,n,n) */
      i  = (i1*(m2-1)+(i2-1))*(m3-1)+(i3-1);
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )	
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )	
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* vertices (i1,i2,i3), (n+1,n,n+1) */
      i  = (i1*(m2-1)+(i2-1))*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )	
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )	
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* vertices (i1,i2,i3), (n+1,n+1,n) */
      i  = (i1*(m2-1)+i2)*(m3-1)+(i3-1);
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )	
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )	
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* vertices (i1,i2,i3), (n+1,n+1,n+1) */
      i  = (i1*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )	
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )	
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )	{
                     	      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      iv[++p] = q;
  }
  
/***********************************************************************/  
  /*row: faces, kind (i1,i2,i3-1), (j1,j2,n+1) 
                   = (i1,i2,i3),   (j1,j2,n) */
  if  ( i3>0 )      	
  for ( j1=0; j1<n1; j1++ ) 
  for ( j2=0; j2<n2; j2++ ) {
      	
  	  /* B1  (i1,i2,i3-1), [(j1,j2,n+1), (k1,k2,k3)] */
  	  i  = ((i1*m2+i2)*m3+i3-1;
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
  	  /* B2  (i1,i2,i3),   [(j1,j2,n  ), (k1,k2,k3)] */
  	  i  = ((i1*m2+i2)*m3+i3;
      for ( k1=ia[3*j1];  k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*j2];  k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*n3];  k3<ic[3*n3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,i3-1), [(j1,j2,n+1), (k1,k2,n)]  */
      i  = ((i3-2)*m1+i1)*m2+i2;
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   )      	        	
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,i3-1), [(j1,j2,n+1), (k1,k2,n+1)]
            +(i1,i2,i3),   [(j1,j2,n  ), (k1,k2,n  )]  */
      i  = ((i3-1)*m1+i1)*m2+i2;
      t  = 0;
      if  ( (k3=ic[3*n3+5])<ic[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
      if  ( (k3=ic[3*n3+1])<ic[3*n3+2] && jc[k3]==n3   ) t += c[k3];
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*t;
          q++;
      }
      
      /* F3, (i1,i2,i3  ), [(j1,j2,n), (k1,k2,n+1)] */
      i  = (i3*m1+i1)*m2+i2;
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 )      	        	
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F2, (i1,i2,i3-1), [(j1,j2, n+1), (k1,n,  k3 )] */
      i  = ((i2-1)*m1+i1)*m3+i3-1;
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2 )      	        	
      for ( k1=ia[3*j1];  k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n+3]; k3<ic[3*n+4];  k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F2, (i1,i2,i3  ), [(j1,j2,n  ), (k1,n, k3 )] */
      i  = ((i2-1)*m1+i1)*m3+i3;
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F2, (i1,i2,i3-1), [(j1,j2, n+1), (k1,n+1,k3 )] */
      i  = (i2*m1+i1)*m3+i3-1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F2, (i1,i2,i3  ), [(j1,j2, n  ), (k1,n+1,k3 )] */
      i  = (i2*m1+i1)*m3+i3;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F1, (i1,i2,i3-1), [(j1,j2, n+1), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2)*m3+i3-1;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F1, (i1,i2,i3  ), [(j1,j2, n  ), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2)*m3+i3;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*j2];  k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*n3];  k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F1, (i1,i2,i3-1), [(j1, j2, n+1), (n+1,k2, k3 )] */
      i  = (i1*m2+i2)*m3+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F1, (i1,i2,i3  ), [(j1, j2, n  ), (n+1,k2, k3 )] */
      i  = (i1*m2+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )
      for ( k2=ib[3*j2];  k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*n3];  k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* E1,  (i1,i2,i3-1), [(j1, j2, n+1), (k1,  n,n  )] */
      i  = ((i2-1)*(m3-1)+(i3-2))*m1+i1;
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      if  ( i3>1   && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E1,  (i1,i2,i3-1), [(j1, j2, n+1), (k1,  n,n+1)]
             +(i1,i2,i3  ), [(j1, j2, n  ), (k1,  n,n  )]  */
      i  = ((i2-1)*(m3-1)+(i3-1))*m1+i1;
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
          t  = 0;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3]; 
      	
          for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
              	
              jv[q] = off[4] + i*n1+jb[k1];
              vv[q] = a[k1]*b[k2]*t;
              q++;
          }          
      }
      
      /* E1,  (i1,i2,i3), [(j1, j2, n  ), (k1, n, n+1)] */
      i  = ((i2-1)*(m3-1)+i3)*m1+i1;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ia[3*n3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[4] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
 
      /* E1,  (i1,i2,i3-1), [(j1, j2, n+1), (k1,  n+1,n  )] */
      i  = (i2*(m3-1)+(i3-2))*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      if  ( i3>1    && (k3=ic[3*n3+4])<ia[3*n3+5] && jc[k3]==n3   ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E1,  (i1,i2,i3-1), [(j1, j2, n+1), (k1,  n+1,n+1)]
             +(i1,i2,i3  ), [(j1, j2, n  ), (k1,  n+1,n  )]  */
      i = (i2*(m3-1)+(i3-1))*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
          t  = 0;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3]; 
      	
          for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
              	
              jv[q] = off[4] + i*n1+jb[k2];
              vv[q] = a[k1]*b[k2]*t;
              q++;
          }          
      }
      
      /* E1,  (i1,i2,i3), [(j1, j2, n  ), (k1, n+1, n+1)] */
      i = (i2*(m3-1)+i3)*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ia[3*n3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[4] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
              
      /* E2, (i1,i2,i3-1), [(j1, j2, n+1), (n,  k2, n  )] */
      i = ((i1-1)*(m3-1)+(i3-2))*m2+i2;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 ) 
      if  ( i3>1   && (k3=ic[3*n3+4])<ia[3*n3+5] && jc[k3]==n3 ) 
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2,   (i1,i2,i3-1), [(j1, j2, n+1), (n,  k2, n+1)]
              +(i1,i2,i3  ), [(j1, j2, n  ), (n,  k2, n  )]  */
      i  = ((i1-1)*(m3-1)+(i3-1))*m2+i2;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) {
          t  = 0;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3]; 
      	
          for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
              	
              jv[q] = off[5] + i*n2+jb[k2];
              vv[q] = a[k1]*b[k2]*t;
              q++;
          }          
      }
      
      /* E2,  (i1,i2,i3), [(j1, j2, n  ), (n,  k2, n+1)] */
      i  = ((i1-1)*(m3-1)+i3)*m2+i2;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ia[3*n3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
            
      /* E2,  (i1,i2,i3-1), [(j1, j2, n+1), (n+1,  k2, n  )] */
      i  = (i1*(m3-1)+(i3-2))*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3>1    && (k3=ic[3*n3+4])<ia[3*n3+5] && jc[k3]==n3   ) 
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2,  (i1,i2,i3-1), [(j1, j2, n+1), (n+1,  k2, n+1)]
             +(i1,i2,i3  ), [(j1, j2, n  ), (n+1,  k2, n  )]  */
      i  = (i1*(m3-1)+(i3-1))*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
          t  = 0;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3]; 
      	
          for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
              	
              jv[q] = off[5] + i*n2+jb[k2];
              vv[q] = a[k1]*b[k2]*t;
              q++;
          }          
      }
      
      /* E2,  (i1,i2,i3), [(j1, j2, n  ), (n+1,  k2, n+1)] */
      i  = (i1*(m3-1)+i3)*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ia[3*n3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*j2]; k2<ib[3*j2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      
      /* E3, (i1,i2,i3-1),  [(j1, j2, n+1), (n,  n,  k3 )] */
      i  = ((i1-1)*(m2-1)+(i2-1))*m3+i3-1;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 ) 
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2 ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E3,  (i1,i2,i3  ), [(j1, j2, n  ), (n,  n,  k3 )] */
      i  = ((i1-1)*(m2-1)+(i2-1))*m3+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2 ) 
      for ( k3=ic[3*n3]; k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
            
      /* E3, (i1,i2,i3-1),  [(j1, j2, n+1), (n,  n+1,  k3 )] */
      i  = ((i1-1)*(m2-1)+i2)*m3+i3-1;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E3,  (i1,i2,i3  ), [(j1, j2, n  ), (n,  n+1,  k3 )] */
      i  = ((i1-1)*(m2-1)+i2)*m3+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*n3]; k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
            
      /* E3, (i1,i2,i3-1),  [(j1, j2, n+1), (n+1,  n,  k3 )] */
      i  = (i1*(m2-1)+(i2-1))*m3+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2  ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E3,  (i1,i2,i3  ), [(j1, j2, n  ), (n+1,  n,  k3 )] */
      i  = (i1*(m2-1)+(i2-1))*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      for ( k3=ic[3*n3]; k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
            
            
      /* E3, (i1,i2,i3-1),  [(j1, j2, n+1), (n+1,  n+1,  k3 )] */
      i  = (i1*(m2-1)+i2)*m3+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E3,  (i1,i2,i3  ), [(j1, j2, n  ), (n+1,  n+1,  k3 )] */
      i  = (i1*(m2-1)+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*n3]; k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n,  n,  n  )] */
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3-2;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n,  n,  n+1)]
             +(i1,i2,i3  ), [(j1, j2, n  ), (n,  n,  n  )]  */                  
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3-1;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) {
      	  t = 0.;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3];
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*t;
          q++;          
      }
      
          
      /* V ,  (i1,i2,i3  ), [(j1, j2, n  ), (n,  n,  n+1  )] */
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n,  n+1,  n  )] */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+i3-2;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n,  n+1,  n+1)]
             +(i1,i2,i3  ), [(j1, j2, n  ), (n,  n+1,  n  )]  */                  
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3-1;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) {
      	  t = 0.;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3];
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*t;
          q++;          
      }
      
      /* V ,  (i1,i2,i3  ), [(j1, j2, n  ), (n,  n+1,  n+1  )] */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n+1,  n,  n  )] */
      i  = (i1*(m2-1)+i2-1)*(m3-1)+i3-2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n+1,  n,  n+1)]
             +(i1,i2,i3  ), [(j1, j2, n  ), (n+1,  n,  n  )]  */                  
      i  = (i1*(m2-1)+i2-1)*(m3-1)+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) {
      	  t = 0.;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3];
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*t;
          q++;          
      }
          
      /* V ,  (i1,i2,i3  ), [(j1, j2, n  ), (n+1,  n,  n+1  )] */
      i  = (i1*(m2-1)+i2-1)*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      

      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n+1,  n+1,  n  )] */
      i  = (i1*(m2-1)+i2)*(m3-1)+i3-2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V ,  (i1,i2,i3-1), [(j1, j2, n+1), (n+1,  n+1,  n+1)]
             +(i1,i2,i3  ), [(j1, j2, n  ), (n+1,  n+1,  n  )]  */                  
      i  = (i1*(m2-1)+i2)*(m3-1)+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) {
      	  t = 0.;
          if  ( (k3=ic[3*n3+5])<ia[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
          if  ( (k3=ic[3*n3+1])<ia[3*n3+2] && jc[k3]==n3   ) t += c[k3];
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*t;
          q++;          
      }
          
      /* V ,  (i1,i2,i3  ), [(j1, j2, n  ), (n+1,  n+1,  n+1  )] */
      i  = (i1*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      iv[++p] = q;
  }                 
/***************************************************************/         
            

/***********************************************************************/  
  /*row: faces, kind (i1,i2-1,i3), (j1,n+1,j3) 
                   = (i1,i2,  i3), (j1,n,  j3) */
  if  ( i2>0 )      	
  for ( j1=0; j1<n1; j1++ ) 
  for ( j3=0; j3<n3; j3++ ) {
      	
  	  /* B1  (i1,i2-1,i3), [(j1,n+1,j3), (k1,k2,k3)] */
  	  i  = (i1*m2+i2-1)*m3+i3-1;
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
  	  /* B1  (i1,i2,i3), [(j1,n,  j3), (k1,k2,k3)] */
  	  i  = (i1*m2+i2)*m3+i3-1;
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F3, (i1,i2-1,i3), [(j1,n+1,j3), (k1,k2,n)]  */
      i  = ((i3-1)*m1+i1)*m2+i2-1;
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,i3), [(j1,n,j3), (k1,k2,n)]  */
      i  = ((i3-1)*m1+i1)*m2+i2;
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F3, (i1,i2-1,i3), [(j1,n+1,j3), (k1,k2,n+1)]  */
      i  = (i3*m1+i1)*m2+i2-1;
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,i3), [(j1,n,j3), (k1,k2,n+1)]  */
      i  = (i3*m1+i1)*m2+i2;
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F2, (i1,i2-1,i3), [(j1,n+1,j2), (k1,n,  k3 )] */
      i  = ((i2-2)*m1+i1)*m3+i3;
      if  ( i2>1   && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F2, (i1,i2-1,i3), [(j1,n+1,j2), (k1,n+1,k3 )]
            +(i1,i2,  i3), [(j1,n  ,j2), (k1,n,  k3 )] */
      i  = ((i2-1)*m1+i1)*m3+i3;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];   	        	
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];   	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*t*c[k3];
          q++;
      }

      /* F2, (i1,i2,i3), [(j1,n,j2), (k1,n+1,  k3 )] */
      i  = (i2*m1+i1)*m3+i3;
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F1, (i1,i2-1,i3), [(j1,n+1,j3), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2-1)*m3+i3;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F1, (i1,i2,i3), [(j1,n,j3), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2)*m3+i3;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F1, (i1,i2-1,i3), [(j1,n+1,j3), (n+1, k2, k3 )] */
      i  = (i1*m2+i2-1)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )      	        	
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F1, (i1,i2,i3), [(j1,n,j3), (n, k2, k3 )] */
      i  = (i1*m2+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )      	        	
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
            
      /* E1,  (i1,i2-1,i3), [(j1,n+1,j3), (k1,  n,n  )] */
      i  = ((i2-2)*(m3-1)+i3-1)*m1+i1;
      if  ( i2>1   && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   ) 
      if  ( i3>0   && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E1,  (i1,i2-1,i3), [(j1,n+1,j3), (k1,  n,n+1)] */
      i  = ((i2-2)*(m3-1)+i3)*m1+i1;
      if  ( i2>1   && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   ) 
      if  ( i3<m-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      

      /* E1,  (i1,i2-1,i3), [(j1,n+1,j3), (k1,  n+1,n  )] 
             +(i1,i2,  i3), [(j1,n,  j3), (k1,  n,  n  )] */
      i  = ((i2-1)*(m3-1)+i3-1)*m1+i1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i3>0   && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }
      
      /* E1,  (i1,i2-1,i3), [(j1,n+1,j3), (k1,  n+1,n+1)] 
             +(i1,i2,  i3), [(j1,n,  j3), (k1,  n,  n+1)] */
      i  = ((i2-1)*(m3-1)+i3)*m1+i1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }      
      
      /* E1,  (i1,i2,i3), [(j1,n,j3), (k1,  n+1,n  )] */
      i  = (i2*(m3-1)+i3-1)*m1+i1;
      if  ( i2>1   && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 ) 
      if  ( i3>0   && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E1,  (i1,i2,i3), [(j1,n,j3), (k1,  n+1,n+1)] */
      i  = (i2*(m3-1)+i3)*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
              
      /* E2, (i1,i2-1,i3), [(j1, n+1, j3), (n,  k2, n  )] */
      i = ((i1-1)*(m3-1)+(i3-1))*m2+i2-1;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 ) 
      if  ( i3>0   && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3 ) 
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1,i2,i3), [(j1, n, j3), (n,  k2, n  )] */
      i = ((i1-1)*(m3-1)+(i3-1))*m2+i2;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 ) 
      if  ( i3>0   && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3 ) 
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1,i2-1,i3), [(j1, n+1, j3), (n,  k2, n+1  )] */
      i = ((i1-1)*(m3-1)+(i3-1))*m2+i2-1;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i3>0   && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1,i2,i3), [(j1, n, j3), (n,  k2, n+1  )] */
      i = ((i1-1)*(m3-1)+(i3-1))*m2+i2;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i3>0   && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E2, (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  k2, n  )] */
      i = (i1*(m3-1)+(i3-1))*m2+i2-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1,i2,i3), [(j1, n, j3), (n+1,  k2, n  )] */
      i = (i1*(m3-1)+(i3-1))*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3  ) 
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
      
      /* E2, (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  k2, n+1  )] */
      i = (i1*(m3-1)+(i3-1))*m2+i2-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1,i2,i3), [(j1, n, j3), (n+1,  k2, n+1  )] */
      i = (i1*(m3-1)+(i3-1))*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      

      /* E3, (i1,i2-1,i3), [(j1, n+1, j3), (n,  n,  k3 )] */
      i  = ((i1-1)*(m2-1)+(i2-2))*m3+i3;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 ) 
      if  ( i2>1   && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E3, (i1,i2-1,i3), [(j1, n+1, j3), (n,  n+1,  k3 )]
            +(i1,i2,  i3), [(j1, n,   j3), (n,  n,    k3 )] */
      i  = ((i1-1)*(m2-1)+(i2-1))*m3+i3;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }
            
      /* E3, (i1,i2,  i3), [(j1, n, j3), (n,  n+1,  k3 )] */
      i  = ((i1-1)*(m2-1)+i2)*m3+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      

      /* E3, (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  n,  k3 )] */
      i  = (i1*(m2-1)+(i2-2))*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* E3, (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  n+1,  k3 )]
            +(i1,i2,  i3), [(j1, n,   j3), (n+1,  n,    k3 )] */
      i  = (i1*(m2-1)+(i2-1))*m3+i3;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }
            
      /* E3, (i1,i2,  i3), [(j1, n, j3), (n+1,  n+1,  k3 )] */
      i  = (i1*(m2-1)+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
            
      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n,  n,  n  )] */
      i  = ((i1-1)*(m2-1)+i2-2)*(m3-1)+i3-1;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n,  n,  n+1  )] */
      i  = ((i1-1)*(m2-1)+i2-2)*(m3-1)+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n,  n+1,  n  )] 
            +(i1,i2,  i3), [(j1, n,   j3), (n,  n,    n  )] */
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3-1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }

      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n,  n+1,  n+1  )] 
            +(i1,i2,  i3), [(j1, n,   j3), (n,  n,    n+1  )] */
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3-1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }
      
      /* V,  (i1,i2,i3), [(j1, n, j3), (n,  n+1,  n  )] */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+i3-1;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1  )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2,i3), [(j1, n, j3), (n,  n+1,  n+1  )] */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1>0    && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1  )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  n,  n  )] */
      i  = (i1*(m2-1)+i2-2)*(m3-1)+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  n,  n+1  )] */
      i  = (i1*(m2-1)+i2-2)*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  n+1,  n  )] 
            +(i1,i2,  i3), [(j1, n,   j3), (n+1,  n,    n  )] */
      i  = (i1*(m2-1)+i2-1)*(m3-1)+i3-1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }

      /* V,  (i1,i2-1,i3), [(j1, n+1, j3), (n+1,  n+1,  n+1  )] 
            +(i1,i2,  i3), [(j1, n,   j3), (n+1,  n,    n+1  )] */
      i  = (i1*(m2-1)+i2-1)*(m3-1)+i3-1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*t*c[k3];
          q++;          
      }
      
      /* V,  (i1,i2,i3), [(j1, n, j3), (n+1,  n+1,  n  )] */
      i  = (i1*(m2-1)+i2)*(m3-1)+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2,i3), [(j1, n, j3), (n+1,  n+1,  n+1  )] */
      i  = (i1*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
                  
      iv[++p] = q;
  }  
                 
/***************************************************************/ 

/***********************************************************************/  
  /*row: faces, kind (i1-1,i2,i3), (n+1,j2,j3) 
                   = (i1,  i2,i3), (n,  j2,j3) */
  if  ( i1>0 )      	
  for ( j2=0; j2<n2; j1++ ) 
  for ( j3=0; j3<n3; j3++ ) {
      	
  	  /* B1  (i1-1,i2,i3), [(n+1,j2,j3), (k1,k2,k3)] */
  	  i  = ((i1-1)*m2+i2)*m3+i3;
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ ) 
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
  	  /* B1  (i1,i2,i3), [(n,  j2,j3), (k1,k2,k3)] */
  	  i  = (i1*m2+i2)*m3+i3;
      for ( k1=ia[3*n1];   k1<ia[3*n1+1]; k1++ ) 
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F3, (i1-1,i2,i3), [(n+1,j2,j3), (k1,k2,n)]  */
      i  = ((i3-1)*m1+i1-1)*m2+i2;
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )      	        	
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ )
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,i3), [(n,j2,j3), (k1,k2,n)]  */
      i  = ((i3-1)*m1+i1)*m2+i2;
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   )      	        	
      for ( k1=ia[3*n1];   k1<ia[3*n1+1]; k1++ )
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F3, (i1-1,i2,i3), [(n+1,j2,j3), (k1,k2,n+1)]  */
      i  = (i3*m1+i1-1)*m2+i2;
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )      	        	
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ )
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,i3), [(n,j2,j3), (k1,k2,n+1)]  */
      i  = (i3*m1+i1)*m2+i2;
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 )      	        	
      for ( k1=ia[3*n1];   k1<ia[3*n1+1]; k1++ )
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
            
      /* F2, (i1-1,i2,i3), [(n+1,j2,j3), (k1,n,  k3 )] */
      i  = ((i2-1)*m1+i1-1)*m3+i3;
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2 )      	        	
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F2, (i1,i2,i3), [(n,j2,j3), (k1,n,  k3 )] */
      i  = ((i2-1)*m1+i1)*m3+i3;
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2 )      	        	
      for ( k1=ia[3*n1];   k1<ia[3*n1+1]; k1++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F2, (i1-1,i2,i3), [(n+1,j2,j3), (k1,n+1,  k3 )] */
      i  = ((i2-1)*m1+i1-1)*m3+i3;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )      	        	
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F2, (i1,i2,i3), [(n,j2,j3), (k1,n+1,  k3 )] */
      i  = ((i2-1)*m1+i1)*m3+i3;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )      	        	
      for ( k1=ia[3*n1];   k1<ia[3*n1+1]; k1++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
                  
      /* F1, (i1-1,i2,i3), [(n+1,j2,j3), (n, k2, k3 )] */
      i  = ((i1-2)*m2+i2)*m3+i3;
      if  ( i1>1   && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F1, (i1-1,i2,i3), [(n+1,j2,j3), (n+1, k2, k3 )]
            +(i1,  i2,i3), [(n,  j2,j3), (n,   k2, k3 )] */
      i  = ((i1-1)*m2+i2)*m3+i3;
      t  = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = t*b[k2]*c[k3];
          q++;
      }

      /* F1, (i1-1,i2,i3), [(n+1,j2,j3), (n+1, k2, k3 )] */
      i  = (i1*m2+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 )      	        	
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) 
      for ( k3=ic[3*j3];   k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
            
      /* E1,  (i1-1,i2,i3), [(n+1, j2,j3), (k1,  n,n)] */
      i  = ((i2-1)*(m3-1)+i3-1)*m1+i1-1;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
            
      /* E1,  (i1,i2,i3), [(n, j2,j3), (k1,  n,n)] */
      i  = ((i2-1)*(m3-1)+i3-1)*m1+i1;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k1=ia[3*n1]; k1<ia[3*n1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      

      /* E1,  (i1-1,i2,i3), [(n+1, j2,j3), (k1,  n,n+1)] */
      i  = ((i2-1)*(m3-1)+i3)*m1+i1-1;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
            
      /* E1,  (i1-1,i2,i3), [(n, j2,j3), (k1,  n,n+1)] */
      i  = ((i2-1)*(m3-1)+i3)*m1+i1;
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*n1]; k1<ia[3*n1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
            
      /* E1,  (i1-1,i2,i3), [(n+1, j2,j3), (k1,  n+1,n)] */
      i  = (i2*(m3-1)+i3-1)*m1+i1-1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
            
      /* E1,  (i1,i2,i3), [(n, j2,j3), (k1,  n+1,n)] */
      i  = (i2*(m3-1)+i3-1)*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k1=ia[3*n1]; k1<ia[3*n1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      

      /* E1,  (i1-1,i2,i3), [(n+1, j2,j3), (k1,  n+1,n+1)] */
      i  = (i2*(m3-1)+i3)*m1+i1-1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*n1+3]; k1<ia[3*n1+4]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
            
      /* E1,  (i1,i2,i3), [(n, j2,j3), (k1,  n+1,n+1)] */
      i  = (i2*(m3-1)+i3)*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*n1]; k1<ia[3*n1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      

      /* E2, (i1-1,i2,i3), [(n+1, j2, j3), (n,  k2, n  )] */
      i = ((i1-2)*(m3-1)+(i3-1))*m2+i2;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3 ) 
      for ( k2=ib[3*j2+3]; k2<ib[3*j2+4]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1-1,i2,i3), [(n+1, j2, j3), (n,  k2, n+1  )] */
      i = ((i1-2)*(m3-1)+i3)*m2+i2;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1   ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
              
      /* E2, (i1-1,i2,i3), [(n+1,j2, j3), (n+1,k2, n  )]
            +(i1,  i2,i3), [(n,  j2, j3), (n,  k2, n  )] */
      i = ((i1-1)*(m3-1)+i3-1)*m2+i2;
      t = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k2=ib[3*j2+3]; k2<ib[3*j2+4]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1-1,i2,i3), [(n+1,j2, j3), (n+1,k2, n+1  )]
            +(i1,  i2,i3), [(n,  j2, j3), (n,  k2, n+1  )] */
      i = ((i1-1)*(m3-1)+i3)*m2+i2;
      t = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = t*b[k2]*c[k3];
          q++;          
      }     
      
      /* E2, (i1,i2,i3), [(n, j2, j3), (n+1,  k2, n  )] */
      i = (i1*(m3-1)+(i3-1))*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) 
      for ( k2=ib[3*j2+3]; k2<ib[3*j2+4]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }

      /* E2, (i1,i2,i3), [(n, j2, j3), (n+1,  k2, n+1  )] */
      i = (i1*(m3-1)+i3)*m2+i2;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) 
      for ( k2=ib[3*j2];   k2<ib[3*j2+1]; k2++ ) {
          	
          jv[q] = off[5] + i*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      

      /* E3, (i1-1,i2,  i3), [(n+1, j2, j3), (n,  n,  k3 )] */
      i  = (i1*(m2-2)+i2-1)*m3+i3;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1   ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
                  
      /* E3, (i1-1,i2,  i3), [(n+1,j2, j3), (n,  n+1,  k3 )] */
      i  = ((i1-2)*(m2-1)+i2)*m3+i3;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
         
      /* E3, (i1-1,i2,i3), [(n+1,j2, j3), (n+1, n,  k3 )]
            +(i1,  i2,i3), [(n,  j2, j3), (n,   n,  k3 )] */
      i  = ((i1-1)*(m2-1)+i2-1)*m3+i3;
      t  = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i2>0   && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = t*b[k2]*c[k3];
          q++;          
      }
                  
      /* E3, (i1-1,i2,i3), [(n+1,j2, j3), (n+1, n+1,  k3 )]
            +(i1,  i2,i3), [(n,  j2, j3), (n,   n+1,  k3 )] */
      i  = ((i1-1)*(m2-1)+i2)*m3+i3;
      t  = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = t*b[k2]*c[k3];
          q++;          
      }

      /* E3, (i1,i2,  i3), [(n, j2, j3), (n+1,  n,  k3 )] */
      i  = (i1*(m2-1)+i2-1)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
                  
      /* E3, (i1,  i2,  i3), [(n,  j2, j3), (n+1,  n+1,  k3 )] */
      i  = (i1*(m2-1)+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 ) 
      for ( k3=ic[3*j3]; k3<ic[3*j3+1]; k3++ ) {
          	
          jv[q] = off[6] + i*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n,  n,  n  )] */
      i  = ((i1-2)*(m2-1)+i2-1)*(m3-1)+i3-1;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1   ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n,  n,  n+1  )] */
      i  = ((i1-2)*(m2-1)+i2-1)*(m3-1)+i3;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1   ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n,  n+1,  n  )] */
      i  = (i1*(m2-2)+i2)*(m3-1)+i3-1;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n,  n+1,  n+1  )] */
      i  = ((i1-2)*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1>1    && (k1=ia[3*n1+4])<ia[3*n1+5] && ja[k1]==n1   ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
            
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n+1,n,  n  )]
            +(i1,  i2,i3), [(n,   j2, j3), (n,  n,  n  )]  */
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3-1;
      t  = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = t*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n+1,n,  n+1  )]
            +(i1,  i2,i3), [(n,   j2, j3), (n,  n,  n+1  )]  */
      i  = ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3;
      t  = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = t*b[k2]*c[k3];
          q++;          
      }
      
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n+1,n+1,  n  )]
            +(i1,  i2,i3), [(n,   j2, j3), (n,  n+1,  n  )]  */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+i3-1;
      t  = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = t*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1-1,i2,i3), [(n+1, j2, j3), (n+1,n+1,  n+1  )]
            +(i1,  i2,i3), [(n,   j2, j3), (n,  n+1,  n+1  )]  */
      i  = ((i1-1)*(m2-1)+i2)*(m3-1)+i3;
      t  = 0.;
      if  ( (k1=ia[3*n1+5])<ia[3*n1+6] && ja[k1]==n1+1 ) t += a[k1];
      if  ( (k1=ia[3*n1+1])<ia[3*n1+2] && ja[k1]==n1   ) t += a[k1];
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = t*b[k2]*c[k3];
          q++;          
      }

      /* V,  (i1,i2,i3), [(n, j2, j3), (n+1,  n,  n  )] */
      i  = (i1*(m2-1)+i2-1)*(m3-1)+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2,i3), [(n, j2, j3), (n+1,  n,  n+1  )] */
      i  = (i1*(m2-1)+i2-1)*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i2>0    && (k2=ib[3*j2+1])<ib[3*j2+2] && jb[k2]==n2   )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
      
      /* V,  (i1,i2,i3), [(n, j2, j3), (n+1,  n+1,  n  )] */
      i  = (i1*(m2-1)+i2)*(m3-1)+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3>0    && (k3=ic[3*j3+1])<ic[3*j3+2] && jc[k3]==n3   ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
          
      /* V,  (i1,i2,i3), [(n, j2, j3), (n+1,  n+1,  n+1  )] */
      i  = (i1*(m2-1)+i2)*(m3-1)+i3;
      if  ( i1<m1-1 && (k1=ia[3*n1+2])<ia[3*n1+3] && ja[k1]==n1+1 ) 
      if  ( i2<m2-1 && (k2=ib[3*j2+2])<ib[3*j2+3] && jb[k2]==n2+1 )
      if  ( i3<m3-1 && (k3=ic[3*j3+2])<ic[3*j3+3] && jc[k3]==n3+1 ) {
      	
          jv[q] = off[7] + i;
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }
                  
      iv[++p] = q;
  }  
                 
/***************************************************************/ 

/***********************************************************************/  
  /*row: egdes,  (i1,i2-1,i3-1), (j1,n+1,n+1) 
               = (i1,i2-1,i3  ), (j1,n+1,n  )
               = (i1,i2,  i3-1), (j1,n,  n+1)
               = (i1,i2,  i3  ), (j1,n,  n  ) */
  if  ( i2>0 && i3>0 )      	
  for ( j1=0; j1<n1; j1++ ) {
      	
  	  /* B1  (i1,i2-1,i3-1), [(j1,n2+1,n3+1), (k1,k2,k3)] */
  	  i  = (i1*m2+i2-1)*m3+i3-1;
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
  	  /* B1  (i1,i2-1,i3  ), [(j1,n2+1,n3  ), (k1,k2,k3)] */
  	  i  = (i1*m2+i2-1)*m3+i3;
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
  	  /* B1  (i1,i2,  i3-1), [(j1,n2,  n3+1), (k1,k2,k3)] */
  	  i  = (i1*m2+i2)*m3+i3-1;
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
  	  /* B1  (i1,i2,  i3), [(j1,n2,n  3  ), (k1,k2,k3)] */
  	  i  = (i1*m2+i2)*m3+i3;
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
      	
          jv[q] = off[0] + (((i*n1+ja[k1])*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }      
      
      /* F3, (i1,i2-1,i3-1), [(j1,n+1,n+1), (k1,k2,n)]  */
      i  = ((i3-2)*m1+i1)*m2+i2-1;
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,  i3-1), [(j1,n,  n+1), (k1,k2,n)]  */
      i  = ((i3-2)*m1+i1)*m2+i2;
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3   )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F3, (i1,i2-1,i3-1), [(j1,n+1,n+1), (k1,k2,n+1)] 
            +(i1,i2-1,i3  ), [(j1,n+1,n  ), (k1,k2,n  )] */
      i  = ((i3-1)*m1+i1)*m2+i2-1;
      t  = 0.;
      if  ( (k3=ic[3*n3+5])<ic[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];     	        	
      if  ( (k3=ic[3*n3+1])<ic[3*n3+2] && jc[k3]==n3   ) t += c[k3];     	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,  i3-1), [(j1,n,  n+1), (k1,k2,n+1)] 
            +(i1,i2,  i3  ), [(j1,n,  n  ), (k1,k2,n  )] */
      i  = ((i3-1)*m1+i1)*m2+i2;
      t  = 0.;
      if  ( (k3=ic[3*n3+5])<ic[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];     	        	
      if  ( (k3=ic[3*n3+1])<ic[3*n3+2] && jc[k3]==n3   ) t += c[k3];     	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }  
      
      /* F3, (i1,i2-1,i3  ), [(j1,n+1,n  ), (k1,k2,n+1)]  */
      i  = (i3*m1+i1)*m2+i2-1;
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F3, (i1,i2,  i3  ), [(j1,n,  n  ), (k1,k2,n+1)]  */
      i  = (i3*m1+i1)*m2+i2;
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ )
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) {          	
      	  
          jv[q] = off[1] + ((i*n1+ja[k1])*n2+jb[k2];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }      
      
      /* F2, (i1,i2-1,i3-1), [(j1,n+1,n+1), (k1,n,  k3 )] */
      i  = ((i2-2)*m1+i1)*m3+i3-1;
      if  ( i2>1   && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F2, (i1,i2-1,i3  ), [(j1,n+1,n  ), (k1,n,  k3 )] */
      i  = ((i2-2)*m1+i1)*m3+i3;
      if  ( i2>1   && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F2, (i1,i2-1,i3-1), [(j1,n+1,n+1), (k1,n+1,k3 )]
            +(i1,i2,  i3-1), [(j1,n,  n+1), (k1,n,  k3 )] */
      i  = ((i2-1)*m1+i1)*m3+i3-1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];    	        	
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];    	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F2, (i1,i2-1,i3  ), [(j1,n+1,n  ), (k1,n+1,k3 )]
            +(i1,i2,  i3  ), [(j1,n,  n  ), (k1,n,  k3 )] */
      i  = ((i2-1)*m1+i1)*m3+i3;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];    	        	
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];    	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }
      
      /* F2, (i1,i2,  i3-1), [(j1,n+1,n+1), (k1,n+1,k3 )] */
      i  = (i2*m1+i1)*m3+i3-1;
      if  ( i2<m2-1 && (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F2, (i1,i2,  i3  ), [(j1,n+1,n  ), (k1,n+1,k3 )] */
      i  = (i2*m1+i1)*m3+i3;
      if  ( i2<m2-1 && (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 )      	        	
      for ( k1=ia[3*j1];   k1<ia[3*j1+1]; k1++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[2] + ((i*n1+ja[k1])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }  
      
      /* F1, (i1,i2-1,i3-1), [(j1,n+1,n+1), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2-1)*m3+i3-1;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F1, (i1,i2-1,i3  ), [(j1,n+1,n  ), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2-1)*m3+i3;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }  
      
      /* F1, (i1,i2,  i3-1), [(j1,n,  n+1), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2)*m3+i3-1;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F1, (i1,i2,  i3  ), [(j1,n,  n  ), (n, k2, k3 )] */
      i  = ((i1-1)*m2+i2)*m3+i3;
      if  ( i1>0   && (k1=ia[3*j1+1])<ia[3*j1+2] && ja[k1]==n1 )      	        	
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }   
      
      /* F1, (i1,i2-1,i3-1), [(j1,n+1,n+1), (n+1, k2, k3 )] */
      i  = (i1*m2+i2-1)*m3+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )      	        	
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F1, (i1,i2-1,i3  ), [(j1,n+1,n  ), (n+1, k2, k3 )] */
      i  = (i1*m2+i2-1)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )      	        	
      for ( k2=ib[3*n2+3]; k2<ib[3*n2+4]; k2++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }  
      
      /* F1, (i1,i2,  i3-1), [(j1,n,  n+1), (n+1, k2, k3 )] */
      i  = (i1*m2+i2)*m3+i3-1;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )      	        	
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*n3+3]; k3<ic[3*n3+4]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }

      /* F1, (i1,i2,  i3  ), [(j1,n,  n  ), (n+1, k2, k3 )] */
      i  = (i1*m2+i2)*m3+i3;
      if  ( i1<m1-1 && (k1=ia[3*j1+2])<ia[3*j1+3] && ja[k1]==n1+1 )      	        	
      for ( k2=ib[3*n2];   k2<ib[3*n2+1]; k2++ ) 
      for ( k3=ic[3*n3];   k3<ic[3*n3+1]; k3++ ) {
          	
          jv[q] = off[3] + ((i*n2+jb[k2])*n3+jc[k3];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;
      }  

      /* E1,  (i1,i2-1,i3-1), [(j1, n+1,n+1), (k1,  n,n)] */
      i  = ((i2-2)*(m3-1)+i3-2)*m1+i1;
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2 ) 
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
      
      /* E1,  (i1,i2-1,i3-1), [(j1, n+1,n+1), (k1,  n,n+1)]
              (i1,i2-1,i3  ), [(j1, n+1,n  ), (k1,  n,n  )] */
      i  = ((i2-2)*(m3-1)+i3-1)*m1+i1;
      t  = 0.;
      if  ( (k3=ic[3*n3+5])<ic[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
      if  ( (k3=ic[3*n3+1])<ic[3*n3+2] && jc[k3]==n3   ) t += c[k3];
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }  
           
      /* E1,  (i1,i2-1,i3  ), [(j1, n+1,n  ), (k1,  n,n+1)] */
      i  = ((i2-2)*(m3-1)+i3)*m1+i1;
      if  ( i2>1    && (k2=ib[3*n2+4])<ib[3*n2+5] && jb[k2]==n2   ) 
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      } 
      
      /* E1,  (i1,i2-1,i3-1), [(j1, n+1,n+1), (k1,  n+1,n)] 
             +(i1,i2,  i3-1), [(j1, n,  n+1), (k1,  n+1,n)] */
      i  = ((i2-1)*(m3-1)+i3-2)*m1+i1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
      
      /* E1,  (i1,i2-1,i3-1), [(j1, n+1,n+1), (k1,  n+1,n+1)]
             +(i1,i2-1,i3  ), [(j1, n+1,n  ), (k1,  n+1,n  )]
             +(i1,i2,  i3-1), [(j1, n,  n+1), (k1,  n,  n+1)]
             +(i1,i2,  i3  ), [(j1, n,  n  ), (k1,  n,  n  )] */
      i  = ((i2-1)*(m3-1)+i3-1)*m1+i1;
      t  = 0.;
      if  ( (k3=ic[3*n3+5])<ic[3*n3+6] && jc[k3]==n3+1 ) t += c[k2];
      if  ( (k3=ic[3*n3+1])<ic[3*n3+2] && jc[k3]==n3   ) t += c[k2];
      s  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) s += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) s += b[k2];
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*s*t;
          q++;          
      }  
           
      /* E1,  (i1,i2-1,i3  ), [(j1, n+1,n  ), (k1,  n+1,n+1)] 
             +(i1,i2,  i3  ), [(j1, n,  n  ), (k1,  n+1,n+1)] */
      i  = ((i2-1)*(m3-1)+i3)*m1+i1;
      t  = 0.;
      if  ( (k2=ib[3*n2+5])<ib[3*n2+6] && jb[k2]==n2+1 ) t += b[k2];
      if  ( (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2   ) t += b[k2];
      if  ( i3<m3-1 && (k3=ic[3*n3+5])<ic[3*n3+6] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }   
      
      
      /* E1,  (i1,i2,i3-1), [(j1, n,n+1), (k1,  n+1,n)] */
      i  = (i2*(m3-1)+i3-2)*m1+i1;
      if  ( i2<m3-1 && (k2=ib[3*n2+1])<ib[3*n2+2] && jb[k2]==n2+1 ) 
      if  ( i3>1    && (k3=ic[3*n3+4])<ic[3*n3+5] && jc[k3]==n3  ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }      
      
      /* E1,  (i1,i2,  i3-1), [(j1, n,n+1),   (k1,  n+1,n+1)]
             +(i1,i2,  i3  ), [(j1, n,n  ),   (k1,  n+1,n  )] */
      i  = (i2*(m3-1)+i3-1)*m1+i1;
      t  = 0.;
      if  ( (k3=ic[3*n3+5])<ic[3*n3+6] && jc[k3]==n3+1 ) t += c[k3];
      if  ( (k3=ic[3*n3+1])<ic[3*n3+2] && jc[k3]==n3   ) t += c[k3];
      if  ( i2<m2-1  && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      }  
           
      /* E1,  (i1,i2,i3  ), [(j1, n,n  ), (k1,  n+1,n+1)] */
      i  = (i2*(m3-1)+i3)*m1+i1;
      if  ( i2<m2-1 && (k2=ib[3*n2+2])<ib[3*n2+3] && jb[k2]==n2+1 ) 
      if  ( i3<m3-1 && (k3=ic[3*n3+2])<ic[3*n3+3] && jc[k3]==n3+1 ) 
      for ( k1=ia[3*j1]; k1<ia[3*j1+1]; k1++ ) {
          	
          jv[q] = off[5] + i*n1+ja[k1];
          vv[q] = a[k1]*b[k2]*c[k3];
          q++;          
      } 
      
            	  
  }
}
  
  
   