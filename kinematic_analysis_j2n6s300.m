% Forward kinematics and forward differential kinematics analysis
% of Kinova Jaco2 6 DOF

clear all

syms q1 q2 q3 q4 q5 q6 D1 D2 D3 D4 D5 D6 e2 aa d4b d5b d6b 'real';

%{
D1 = 0.2755;
D2 = 0.41;
D3 = 0.2073;
D4 = 0.0741;
D5 = 0.0741;
D6 = 0.16;
e2 = 0.0098;

aa = pi/6;
d4b = D3 + sin(aa)/sin(2*aa)*D4;
d5b = sin(aa)/sin(2*aa)*(D4+D5);
d6b = sin(aa)/sin(2*aa)*D5 + D6;
%}

%% Forward kinematics

% Classic Denavit Hartenberg parameters
alpha = [pi/2 pi pi/2 pi/3 pi/3 pi];
a = [0 D2 0 0 0 0];
d = [D1 0 -e2 -d4b -d5b -d6b];
theta = [-q1 q2+pi/2 q3-pi/2 q4 q5+pi q6-pi/2];

% Compute DH matrices
A = simplify(DH(alpha(1), a(1), d(1), theta(1)));
for i = 2:6
    A = [A simplify(DH(alpha(i), a(i), d(i), theta(i)))]; 
end

% Compute relative homogeneous transformation matrices
A1_0 = A(:,1:4);
A2_1 = A(:,5:8);
A3_2 = A(:,9:12);
A4_3 = A(:,13:16);
A5_4 = A(:,17:20);
A6_5 = A(:,21:24);

% Compute jaco's kinematic equation
T = A1_0 * A2_1 * A3_2 * A4_3 * A5_4 * A6_5; % homogeneous transformation matrix from frame 6 to frame 0
T = simplify(T);

% End effector's rotation matrix
r11 = T(1,1);
r21 = T(2,1);
r31 = T(3,1);

r12 = T(1,2);
r22 = T(2,2);
r32 = T(3,2);

r13 = T(1,3);
r23 = T(2,3);
r33 = T(3,3);

% End effector's position vector
px = T(1,4);
py = T(2,4);
pz = T(3,4);

%% Forward differential kinematics

% Compute the homogeneous transformation matrices from frame i to inertial frame 0
T1_0 = A1_0;
T2_0 = simplify(T1_0*A2_1);
T3_0 = simplify(T2_0*A3_2);
T4_0 = simplify(T3_0*A4_3);
T5_0 = simplify(T4_0*A5_4);
T6_0 = simplify(T5_0*A6_5);

% Extract position vectors from each homogeneous transformation matrix
p0 = [0 0 0]'; 
p1 = T1_0(1:3,4);
p2 = T2_0(1:3,4);
p3 = T3_0(1:3,4);
p4 = T4_0(1:3,4);
p5 = T5_0(1:3,4);
p6 = T6_0(1:3,4);

% Define vectors around which each link rotates in the precedent coordinate frame
z0 = [0 0 1]';
z1 = T1_0(1:3,3);
z2 = T2_0(1:3,3);
z3 = T3_0(1:3,3);
z4 = T4_0(1:3,3);
z5 = T5_0(1:3,3);
z6 = T6_0(1:3,3);

% Jacobian matrix
j11 = cross(z0, p6-p0);
j12 = cross(z1, p6-p1);
j13 = cross(z2, p6-p2);
j14 = cross(z3, p6-p3);
j15 = cross(z4, p6-p4);
j16 = cross(z5, p6-p5);

j21 = z0;
j22 = z1;
j23 = z2;
j24 = z3;
j25 = z4;
j26 = z5;

J = simplify([j11 j12 j13 j14 j15 j16;
              j21 j22 j23 j24 j25 j26]);  
          
J = vpa(J,2);
          
% Extract linear jacobian
Jl = J(1:3,:);

% Extract angular jacobian
Ja = J(4:6,:);
        