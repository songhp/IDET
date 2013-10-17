%%
% demo_IDET

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

load demo01

n = size(A,2); % signal length


n = size(A,2); % signal length

% initialization
t = 1; % iteration number
xt = zeros(n,1);  % initial x
% xt = pinv(A)*y; % initial x
r = y;  % initial residue
support_size = k;  % initial sparsity level
support = [];
stop = 0; % not convergent

while ~stop

    % support detection
    [val idx] = sort(abs(xt +  A' *r), 'descend');
	support = sort(idx(1:support_size));
    
    % signal estimation
    xt = zeros(n,1);
    xt(support) = pinv(A(:, support)) * y;
    r = y - A*xt;
    
%     if norm(r) < Tolerance || norm(r)/norm(y) < Tolerance  || t >100  
	if  norm(r)/norm(y) < 1e-6 || t >100 
        stop = 1; % convergent
    else
        t = t+1;
    end
    SupportDetection(x, xt)

end



% 
% xt = OMP(A, y, k);
% SupportDetection(x, xt)
% SupportDetection(x, xbar) % 对比原始信号与重构信号

