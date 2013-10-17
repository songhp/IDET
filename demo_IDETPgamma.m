%%
% demo_IDETPgamma

% % options: 'gm'=1 'gk'=2 'bm'=3 'bk'=4
% % g: gaussian 
% % b: binary 
% % m: fixed measurements 
% % k: fixed sparsity level
% clear all
% k = 20;
% m = 80;
% n = 200;
% ensemble = 'USE';
% opts.sigma = 0;
% opts.gauss = 1;
% 
% % [x, y, A] = gen_signal(k, m, n, ensemble, opts);
% [x, y, A] = gen_signal(k, m, n, ensemble, opts);


load demo03


n = size(A,2); % signal length


itr = 0;
xt = zeros(n,1);  % initial x
r = y;
stop = 0;

support = [];
while ~stop
    % support detection    
    actfun = A' * r + x;
    T = setdiff([1:n], support);
    actfun = actfun(T);
    [activeset, val] = find(abs(actfun) > 0.5 * max(abs(actfun)));
    support = union(support, T(activeset));
    
    
    % signal estimation
    % IDE-x approach
    xt = zeros(n,1);
    xt(support) = pinv(A(:, support)) * y;
    r = y - A*xt;
    
% 	% IDE-s approach
%     x = zeros(n, 1);
%     Aa = A(:, act_set);
%     inact_set = setdiff([1:n], act_set);
%     Ai = A(:, inact_set);
%     P  = pinv(Ai*Ai');
%     xa = pinv(Aa' * P * Aa) * Aa' * y;
%     x(act_set) = xa;
%     x(inact_set) = Ai' * P * (y -  Aa * xa);

    
%     if norm(r) < Tolerance || norm(r)/norm(y) < Tolerance  || itr >100  
	if  norm(r)/norm(y) < 1e-6  || itr >100  
        stop = 1; % convergence tolerance
    else
        itr = itr+1;
    end
    
    SupportDetection(x, xt)

end



% 
% xt = OMP(A, y, k);
% SupportDetection(x, xt)
% SupportDetection(x, xbar) % 对比原始信号与重构信号

