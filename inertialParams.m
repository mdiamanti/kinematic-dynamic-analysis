function params = inertialParams(enable)
% Returns a struct with inertial parameters of Kinova jaco2 6 DOF
% obtained from KInova's xacro file

%{ 
Inertia tensor:
ixx=0.083333 * mass * (3*radius*radius + height*height) ixy=0 ixz=0
iyy=0.083333 * mass * (3*radius*radius + height*height) iyz=0 
izz=0.5*mass*radius*radius
%}

% If enable == TRUE then numeric values are returned, else symbolic 
% values are returned

params = struct();

if enable
    
    % Mass of links
    m = [0.7477 0.99 0.6763 0.1785 0.1785 0.727]';
    params.m = cell(6,1);
    params.m = m;

    % Inertia tensor of links
    r = [0.04 0.04 0.03 0.04 0.04 0.04]';
    h = [0.14 0.35 0.15 0.02 0.02 0.03]';
    
    I = [0.083333 * m(1) * (3*r(1)*r(1) + h(1)*h(1)) 0 0; 
         0 0.083333 * m(1) * (3*r(1)*r(1) + h(1)*h(1)) 0; 
         0 0 0.5*m(1)*r(1)*r(1)];
     
    for i=2:6   
        temp = [0.083333 * m(i) * (3*r(i)*r(i) + h(i)*h(i)) 0 0; 
               0 0.083333 * m(i) * (3*r(i)*r(i) + h(i)*h(i)) 0; 
               0 0 0.5*m(i)*r(i)*r(i)];
        I = [I temp];
    end
    
    params.I = cell(3,18);
    params.I = I;
    
    % Center of mass (COM) of link
    com1 = [0 -0.002 -0.0605];
    com2 = [0 -0.2065 -0.01];
    com3 = [0 0.081 -0.0086];
    com4 = [0 -0.037 -0.0642];
    com5 = [0 -0.037 -0.0642];
    com6 = [0 0 -0.06];
    
    params.com = cell(1,18);
    params.com = [com1 com2 com3 com4 com5 com6];
    
else
    
    syms m1 m2 m3 m4 m5 m6 'real';
    syms Ix1 Ix2 Ix3 Ix4 Ix5 Ix6 Iy1 Iy2 Iy3 Iy4 Iy5 Iy6 Iz1 Iz2 Iz3 Iz4 Iz5 Iz6 'real'; 
    syms c1x c1y c1z c2x c2y c2z c3x c3y c3z c4x c4y c4z c5x c5y c5z c6x c6y c6z 'real';
    
    params.m = sym(zeros(6,1));
    params.m = [m1 m2 m3 m4 m5 m6]';
    
    params.I = sym(zeros(3,18));
    params.I(:,1:3) = [Ix1 0 0; 0 Iy1 0; 0 0 Iz1];
    params.I(:,4:6) = [Ix2 0 0; 0 Iy2 0; 0 0 Iz2];
    params.I(:,7:9) = [Ix3 0 0; 0 Iy3 0; 0 0 Iz3];
    params.I(:,10:12) = [Ix4 0 0; 0 Iy4 0; 0 0 Iz4];
    params.I(:,13:15) = [Ix5 0 0; 0 Iy5 0; 0 0 Iz5];
    params.I(:,16:18) = [Ix6 0 0; 0 Iy6 0; 0 0 Iz6];
    
    params.com = sym(zeros(1,18));
    params.com(1:3) = [c1x c1y c1z];
    params.com(4:6) = [c2x c2y c2z];
    params.com(7:9) = [c3x c3y c3z];
    params.com(10:12) = [c4x c4y c4z];
    params.com(13:15) = [c5x c5y c5z];
    params.com(16:18) = [c6x c6y c6z];

end