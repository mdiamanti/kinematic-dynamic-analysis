% Function for the extraction of DH matrix of a joint,
% given the classic DH parameters of that joint

function A = DH(alpha, a, d, theta)
    
a11 = cos(theta);
a12 = -sin(theta)*round(cos(alpha),2);
a13 = sin(theta)*round(sin(alpha),2);
a14 = a*cos(theta);

a21 = sin(theta);
a22 = cos(theta)*round(cos(alpha),2);
a23 = -cos(theta)*round(sin(alpha),2);
a24 = a*sin(theta);

a31 = 0;
a32 = round(sin(alpha),2);
a33 = round(cos(alpha),2);
a34 = d;

a41 = 0;
a42 = 0;
a43 = 0;
a44 = 1;


A = [a11 a12 a13 a14; 
     a21 a22 a23 a24;
     a31 a32 a33 a34;
     a41 a42 a43 a44];
 
end



