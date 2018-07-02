% Inverse dynamics analysis
% of Kinova Jaco2 6 DOF

run('kinematic_analysis_j2n6s300');

syms dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 'real';

q = [q1 q2 q3 q4 q5 q6]';
dq = [dq1 dq2 dq3 dq4 dq5 dq6]';
ddq = [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6]';

%% Inertial parameters

% Run inertialParams(1) to obtain numeric values or inertialParams(0) to
% obtain symbolic values

params = inertialParams(1);

% Mass of links
m = params.m;

% Inertia tensor of links
I1 = params.I(:,1:3);
I2 = params.I(:,4:6);
I3 = params.I(:,7:9);
I4 = params.I(:,10:12);
I5 = params.I(:,13:15);
I6 = params.I(:,16:18);

% Center of mass (COM) of links
com1 = params.com(1:3);
com2 = params.com(4:6);
com3 = params.com(7:9);
com4 = params.com(10:12);
com5 = params.com(13:15);
com6 = params.com(16:18);


%% Compute inertia matrix D

r0_c1 = [eye(3) zeros(3,1)] * T1_0 * [com1 1]';
Jl_1 = jacobian(r0_c1,q);
Ja_1 = [z0 zeros(3,5)];

r0_c2 = [eye(3) zeros(3,1)] * T2_0 * [com2 1]';
Jl_2 = jacobian(r0_c2,q);
Ja_2 = [z0 z1 0 0 0 0];

r0_c3 = [eye(3) zeros(3,1)] * T3_0 * [com3 1]';
Jl_3 = jacobian(r0_c3,q);
Ja_3 = [z0 z1 z2 0 0 0];

r0_c4 = [eye(3) zeros(3,1)] * T4_0 * [com4 1]';
Jl_4 = jacobian(r0_c4,q);
Ja_4 = [z0 z1 z2 z3 0 0];

r0_c5 = [eye(3) zeros(3,1)] * T5_0 * [com5 1]';
Jl_5 = jacobian(r0_c5,q);
Ja_5 = [z0 z1 z2 z3 z4 0];

r0_c6 = [eye(3) zeros(3,1)] * T6_0 * [com6 1]';
Jl_6 = jacobian(r0_c6,q);
Ja_6 = [z0 z1 z2 z3 z4 z5];

% Inertia matrix
D = simplify(m(1) * Jl_1' * Jl_1 + Ja_1' * I1 * Ja_1) + ...
    simplify(m(2) * Jl_2' * Jl_2 + Ja_2' * I2 * Ja_2) + ...
    simplify(m(3) * Jl_3' * Jl_3 + Ja_3' * I3 * Ja_3) + ...
    simplify(m(4) * Jl_4' * Jl_4 + Ja_4' * I4 * Ja_4) + ...
    simplify(m(5) * Jl_5' * Jl_5 + Ja_5' * I5 * Ja_5) + ...
    simplify(m(6) * Jl_6' * Jl_6 + Ja_6' * I6 * Ja_6);
    
D = vpa(D,2);
    
%% Compute gravity vector g

g0 = [0 0 -9.81]; % gravity acceleration vector with respect to frame 0
Jl = [Jl_1 Jl_2 Jl_3 Jl_4 Jl_5 Jl_6];

g = sym(zeros(6,1)); % gravity vector
for i=1:6
    for j=i:6
        g(i) = g(i) + m(j) * g0 * Jl(:,(j-1)*6+i);        
    end
end

g = vpa(g,2);


%% Compute coriolis and centrifugal terms h

h = sym(zeros(6,1));
parfor i=1:6
    for j=1:6
        for k=1:6
            temp1 = q(k);
            temp2 = q(i);
            h(i) = h(i) + (diff(D(i,j),temp1) - 0.5 * diff(D(j,k),temp2))*dq(j)*dq(k);
        end
    end
end

h = vpa(h,2);


%% Final dynamic equation

%tau = D * ddq + h + g;
%tau = vpa(tau,2);


%% Generating matlab scripts

matlabFunction(D, 'vars', {q}, 'file', 'inertiaMatrix', 'Optimize', false);
matlabFunction(g, 'vars', {q}, 'file', 'gravityVector', 'Optimize', false);
matlabFunction(h, 'vars', {q,dq}, 'file', 'cor_centriTerms', 'Optimize', false);

