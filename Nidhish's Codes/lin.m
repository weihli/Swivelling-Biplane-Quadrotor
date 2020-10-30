omg1 = 1; omg2 = 5*omg1;

A = [0,1,0,0;
     -omg1^2, -2*1.2*omg1, 1, 0;
     0,0,0,1;
     0,0,-omg2^2, -2*1.2*omg2];
 
 [vec,val] = eig(A)

 