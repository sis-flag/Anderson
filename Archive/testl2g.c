#include <stdio.h>
#include <stdlib.h>

typedef int INT;


/*  0<=i1,i2,i3<=m-1, 0<=j1,j2,j3<=n-1  without specifying

    (i1*m+i2)*m+i3,                (j1*n+j2)*n+j3
    (i1,  i2,  i3  ) (j1, j2, j3 )

    (i3*m+i1)*m+i2,                j1*n+j2           i3<m-1
    (i1,  i2,  i3  ) (j1, j2, n  )
    (i1,  i2,  i3+1) (j1, j2, n+1)

    (i2*m+i1)*m+i3,                j1*n+j3           i2<m-1
    (i1,  i2,  i3  ) (j1, n,  j3 )
    (i1,  i2+1,i3  ) (j1, n+1,j3 )

    (i1*m+i2)*m+i3,                j2*n+j3           i1<m-1
    (i1,  i2,  i3  ) (n,  j2, j3 )
    (i1+1,i2,  i3  ) (n+1,j2, j3 )

    (i2*(m-1)+i3)*m+i1,            j1               i2,i3<m-1
    (i1,  i2,  i3  ) (j1, n,  n  )
    (i1,  i2,  i3+1) (j1, n,  n+1)
    (i1,  i2+1,i3  ) (j1, n+1,n  )
    (i1,  i2+1,i3+1) (j1, n+1,n+1)

    (i1*(m-1)+i3)*m+i2,            j2               i1,i3<m-1
    (i1,  i2,  i3  ) (n,  j2, n  )
    (i1,  i2,  i3+1) (n,  j2, n+1)
    (i1+1,i2,  i3  ) (n+1,j2, n  )
    (i1+1,i2,  i3+1) (n+1,j2, n+1)

    (i1*(m-1)+i2)*m+i3,            j3               i1,i2<m-1
    (i1,  i2,  i3  ) (n,  n,  j3 )
    (i1,  i2+1,i3  ) (n,  n+1,j3 )
    (i1+1,i2,  i3  ) (n+1,n,  j3 )
    (i1+1,i2+1,i3  ) (n+1,n+1,j3 )

    (i1*(m-1)+i2)*(m-1)+i3,       0               i1,i2,i3<m-1
    (i1,  i2,  i3  ) (n,  n,  n  )
    (i1,  i2,  i3+1) (n,  n,  n+1)
    (i1,  i2+1,i3  ) (n,  n+1,n  )
    (i1,  i2+1,i3+1) (n,  n+1,n+1)
    (i1+1,i2,  i3  ) (n+1,n,  n  )
    (i1+1,i2,  i3+1) (n+1,n,  n+1)
    (i1+1,i2+1,i3  ) (n+1,n+1,n  )
    (i1+1,i2+1,i3+1) (n+1,n+1,n+1)
    
*/


INT l2gindex(INT m1, INT m2, INT m3, INT i1, INT i2, INT i3, 
             INT n1, INT n2, INT n3, INT j1, INT j2, INT j3)
{	     
 INT  off[8], p=-1;

 n1--;
 n2--;
 n3--;

 off[0] = 0;
 off[1] = off[0] + m1*m2*m3*n1*n2*n3;
 off[2] = off[1] + (m3-1)*m1*m2*n1*n2;
 off[3] = off[2] + (m2-1)*m1*m3*n1*n3;
 off[4] = off[3] + (m1-1)*m2*m3*n2*n3;
 off[5] = off[4] + (m2-1)*(m3-1)*m1*n1;
 off[6] = off[5] + (m1-1)*(m3-1)*m2*n2;
 off[7] = off[6] + (m1-1)*(m2-1)*m3*n3;

/* local-to-gobal mapping*/
 if ( j1<n1 ) {
 	 if ( j2<n2 ) {
 	 	 if ( j3<n3 ) { 
 	 	 	 p = off[0] + ((((i1*m2+i2)*m3+i3)*n1+j1)*n2+j2)*n3+j3;    
 	 	 } else if ( j3==n3   && i3<m3-1 ) { /* F3R */
 	 	 	 p = off[1] + (((i3*m1+i1)*m2+i2)*n1+j1)*n2+j2;
 	 	 } else if ( j3==n3+1 && i3>0    ) {  /* F3L */
 	 	 	 p = off[1] + ((((i3-1)*m1+i1)*m2+i2)*n1+j1)*n2+j2; 	 	 	 
 	 	 }
 	 } else if ( j2==n2 && i2<m2-1 ) { 
 	 	 if ( j3<n3 ) {                       /* F2R */
 	 	 	 p = off[2] + (((i2*m1+i1)*m3+i3)*n1+j1)*n3+j3;
 	 	 } else if ( j3==n3   && i3<m3-1 ) {  /* E1RR */
 	 	 	 p = off[4] + ((i2*(m3-1)+i3)*m1+i1)*n1+j1;
 	 	 } else if ( j3==n3+1 && i3>0    ) {  /* E1RL */
 	 	 	 p = off[4] + ((i2*(m3-1)+i3-1)*m1+i1)*n1+j1;
 	 	 }
 	 } else if ( j2==n2+1 && i2>0 ) {
 	 	 if ( j3<n3 ) {                       /* F2L */
 	 	 	 p = off[2] + ((((i2-1)*m1+i1)*m3+i3)*n1+j1)*n3+j3;
 	 	 } else if ( j3==n3   && i3<m3-1 ) {  /* E1LR */
 	 	 	 p = off[4] + (((i2-1)*(m3-1)+i3)*m1+i1)*n1+j1;  
 	 	 } else if ( j3==n3+1 && i3>0    ) {  /* E1LL */
 	 	 	 p = off[4] + (((i2-1)*(m3-1)+i3-1)*m1+i1)*n1+j1;
 	 	 }
 	 }
 } else if ( j1==n1 && i1<m1-1 ) { 
 	 if ( j2<n2 ) { 
 	 	 if ( j3<n3 ) {                     /* F1R */
 	 	 	 p = off[3] + (((i1*m2+i2)*m3+i3)*n2+j2)*n3+j3;
 	 	 } else if ( j3==n3 && i3<m3-1 ) {  /* E2RR */
 	 	 	 p = off[5] + ((i1*(m3-1)+i3)*m2+i2)*n2+j2;
 	 	 } else if ( j3==n3+1 && i3>0 ) {   /* E2RL */
 	 	 	 p = off[5] + ((i1*(m3-1)+i3-1)*m2+i2)*n2+j2; 
 	 	 }   
 	 } else if ( j2==n2 && i2<m2-1 ) { 
 	 	 if ( j3<n3 ) {                     /* E3RR */
 	 	 	 p = off[6] + ((i1*(m2-1)+i2)*m3+i3)*n3+j3;
 	 	 } else if ( j3==n3 && i3<m3-1 ) {  /* VRRR */
 	 	 	 p = off[7] + (i1*(m2-1)+i2)*(m3-1)+i3;
 	 	 } else if ( j3==n3+1 && i3>0 ) {   /* VRRL */
 	 	 	 p = off[7] + (i1*(m2-1)+i2)*(m3-1)+i3-1;
 	 	 }
 	 } else if ( j2==n2+1 && i2>0 ) { 
 	 	 if ( j3<n3 ) {                    /* E3RL */
 	 	 	 p = off[6] + ((i1*(m2-1)+i2-1)*m3+i3)*n3+j3;
 	 	 } else if ( j3==n3 && i3<m3-1 ) { /* VRLR */
 	 	 	 p = off[7] + (i1*(m2-1)+i2-1)*(m3-1)+i3;
 	 	 } else if ( j3==n3+1 && i3>0 ) {  /* VRLL */
 	 	 	 p = off[7] + (i1*(m2-1)+i2-1)*(m3-1)+i3-1;
 	 	 }   
 	 }
 } else if ( j1==n1+1 && i1>0 ) {
 	 if ( j2<n2 ) {
 	 	 if ( j3<n3 ) {                    /* F1L */
 	 	 	 p = off[3] + ((((i1-1)*m2+i2)*m3+i3)*n2+j2)*n3+j3;
 	 	 } else if ( j3==n3 && i3<m3-1 ) { /* E2LR */
 	 	 	 p = off[5] + (((i1-1)*(m3-1)+i3)*m2+i2)*n2+j2;
 	 	 } else if ( j3==n3+1 && i3>0 ) {  /* E2LL */
 	 	 	 p = off[5] + (((i1-1)*(m3-1)+i3-1)*m2+i2)*n2+j2;
 	 	 }   
 	 } else if ( j2==n2 && i2<m2-1 ) { 
 	 	 if ( j3<n3 ) { 	 	 	       /* E3LR */
 	 	 	 p = off[6] + (((i1-1)*(m2-1)+i2)*m3+i3)*n3+j3;  
 	 	 } else if ( j3==n3 && i3<m3-1 ) { /* VLRR */
 	 	 	 p = off[7] + ((i1-1)*(m2-1)+i2)*(m3-1)+i3;
 	 	 } else if ( j3==n3+1 && i3>0 ) {  /* VLRL */
 	 	 	 p = off[7] + ((i1-1)*(m2-1)+i2)*(m3-1)+i3-1;
 	 	 }   
 	 } else if ( j2=n2+1 && i2>0 ) { 
 	 	 if ( j3<n3 ) {                    /* E3LL */
 	 	 	 p = off[6] + (((i1-1)*(m2-1)+i2-1)*m3+i3)*n3+j3;
 	 	 } else if ( j3==n3 && i3<m3-1 ) { /* VLLR */
 	 	 	 p = off[7] + ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3;
 	 	 } else if ( j3==n3+1 && i3>0 ) {  /* VLLL */
 	 	 	 p = off[7] + ((i1-1)*(m2-1)+i2-1)*(m3-1)+i3-1;
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
{ INT m1=3, m2=3, m3=3, n1=3, n2=3, n3=3, p=0;
  INT i1, i2, i3, j1, j2, j3, g, q, c, *nz;

 if (argc>1) m1 = atoi(argv[1]);
 if (argc>2) m2 = atoi(argv[2]);
 if (argc>3) m3 = atoi(argv[3]);
 if (argc>4) n1 = atoi(argv[4]);
 if (argc>5) n2 = atoi(argv[5]);
 if (argc>6) n3 = atoi(argv[6]);
 nz = (INT *) calloc((m1*n1-1)*(m2*n2-1)*(m3*n3-1),sizeof(INT));

 //printf("testing I: %d...\n",p); 
 for ( i1=0; i1<m1; i1++ )
 for ( i2=0; i2<m2; i2++ ) 
 for ( i3=0; i3<m3; i3++ ) 
     for ( j1=0; j1<n1-1; j1++ )	
     for ( j2=0; j2<n2-1; j2++ )	
     for ( j3=0; j3<n3-1; j3++ ) {	
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
         
     }
     

 //printf("testing F3: %d...\n",p); 
 q  = p;   
 for ( j3=n3-1; j3<=n3; j3++ ) {
 p = q;
 for ( i3=j3-n3+1; i3<j3-n3+m3; i3++ ) 
 for ( i1=0; i1<m1; i1++ )
 for ( i2=0; i2<m2; i2++ ) 
     for ( j1=0; j1<n1-1; j1++ )	
     for ( j2=0; j2<n2-1; j2++ )	{
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
     } 
 }   
        
 //printf("testing F2: %d...\n",p); 
 q  = p;   
 for ( j2=n2-1; j2<=n2; j2++ ) {
 p = q;
 for ( i2=j2-n2+1; i2<j2-n2+m2; i2++ ) 
 for ( i1=0; i1<m1; i1++ )
 for ( i3=0; i3<m3; i3++ ) 
     for ( j1=0; j1<n1-1; j1++ )	
     for ( j3=0; j3<n3-1; j3++ )	{
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
     } 
 } 
 
 //printf("testing F1: %d...\n",p); 
 q  = p;   
 for ( j1=n1-1; j1<=n1; j1++ ) {
 p = q;
 for ( i1=j1-n1+1; i1<j1-n1+m1; i1++ ) 
 for ( i2=0; i2<m2; i2++ )
 for ( i3=0; i3<m3; i3++ ) 
     for ( j2=0; j2<n2-1; j2++ )	
     for ( j3=0; j3<n3-1; j3++ )	{
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
     } 
 }  
        
 //printf("testing E1: %d...\n",p); 
 q  = p;   
 for ( j2=n2-1; j2<=n2; j2++ ) 
 for ( j3=n3-1; j3<=n3; j3++ ) {
 p = q;
 for ( i2=j2-n2+1; i2<j2-n2+m2; i2++ ) 
 for ( i3=j3-n3+1; i3<j3-n3+m3; i3++ ) 
 for ( i1=0; i1<m1; i1++ )
     for ( j1=0; j1<n1-1; j1++ )	{
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
     } 
 }  
 
 //printf("testing E2: %d...\n",p); 
 q  = p;   
 for ( j1=n1-1; j1<=n1; j1++ ) 
 for ( j3=n3-1; j3<=n3; j3++ ) {
 p = q;
 for ( i1=j1-n1+1; i1<j1-n1+m1; i1++ ) 
 for ( i3=j3-n3+1; i3<j3-n3+m3; i3++ ) 
 for ( i2=0; i2<m2; i2++ )
     for ( j2=0; j2<n2-1; j2++ )	{
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
     }
 }   
     
 //printf("testing E3: %d...\n",p); 
 q  = p;   
 for ( j1=n1-1; j1<=n1; j1++ ) 
 for ( j2=n2-1; j2<=n2; j2++ ) {
 p = q;
 for ( i1=j1-n1+1; i1<j1-n1+m1; i1++ ) 
 for ( i2=j2-n2+1; i2<j2-n2+m2; i2++ ) 
 for ( i3=0; i3<m3; i3++ )
     for ( j3=0; j3<n3-1; j3++ )	{
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
     }
 } 
 
 //printf("testing V: %d...\n",p); 
 q  = p;   
 for ( j1=n1-1; j1<=n1; j1++ ) 
 for ( j2=n2-1; j2<=n2; j2++ ) 
 for ( j3=n3-1; j3<=n3; j3++ ) {
 p = q;
 for ( i1=j1-n1+1; i1<j1-n1+m1; i1++ ) 
 for ( i2=j2-n2+1; i2<j2-n2+m2; i2++ ) 
 for ( i3=j3-n3+1; i3<j3-n3+m3; i3++ ) {
 	
         g=l2gindex(m1,m2,m3, i1,i2,i3, n1,n2,n3, j1, j2, j3);
         if (g-p) 
         	//printf(" [%d %d %d], [%d %d %d], [%d %d]\n",i1,i2,i3,j1,j2,j3,g,p);
         p++;
     }
 }       
         
 free(nz);      
}