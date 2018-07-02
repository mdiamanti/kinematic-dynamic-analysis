function dqq = equation_of_motion(tau,q,dq)
% Get joint accelerations dqq by solving forward dynamics equation

dqq  = inertiaMatrix(q)\(tau - cor_centriTerms(q,dq));

end

