function [xt out] = IDET(A, y, k, varargin)
% Iterative Detection Estimation with Thresholding
% Coded by Heping Song (hepingsong@gmail.com)
% 
%%%%%%%%%%%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT(compulsory):
%   A: m*n measurement matrix
%   y: m*1 measurements
%   k: sparsity level
%   
% INPUT(optional):
% 'At': the adjoint function handle that computes A^T*y for function handle
% 'MaxIt': the maximum gradient steps
% 'Tolerance': the convergence tolerance
%
%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   xt: reconstructed sparse signal
%   out.iter: number of iterations

path(path, './subfunctions');

% A is a matrix or function handle 
explicitA = ~(ischar(A) || isa(A, 'function_handle'));

maxiter = 100; % maximum number of iterations 
Tolerance = 1e-6; % convergence tolerance

%Read the optional inputs
if (rem(length(varargin),2)==1)
    error('Optional inputs must go by pairs!');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case upper('At')
                At = varargin{i+1};                
            case upper('MaxIt')
                MaxIt = varargin{i+1};
            case upper('Tolerance')
                Tolerance=varargin{i+1};
            otherwise
                error(['Unrecognized optional input: ''' varargin{i} '''']);
        end
    end
end



if (explicitA)
	n = size(A,2); % signal length
else
	n= length(At(y));
end

%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%
t = 1; % iteration number
xt = zeros(n,1);  % initial x
r = y;  % initial residue
support_size = k;  % initial sparsity level
stop = 0; % not convergent


while ~stop

    if (explicitA)
    % support detection
    [val idx] = sort(abs(xt + A' *r), 'descend');
	support = sort(idx(1:support_size));
    
    % signal estimation
    xt = zeros(n,1);
    xt(support) = pinv(A(:, support)) * y;
    r = y - A*xt;
    
    else
    % support detection
    [val idx] = sort(abs(xt + At(r)), 'descend');
	support = sort(idx(1:support_size));
    
    % signal estimation
    xt = zeros(n,1);
    [xt r] =MySubsetCG(y, xt, A, At, support, 1e-16, 0, MaxIt); 

% 	zz = zeros(n,1);
%     zz(support) = 1;
%     xt = zeros(n,1);
%     [xt r] =MySubsetCG(y, xt, A, At, find(zz~=0), 1e-9, 0, MaxIt); 

    r = y - A(xt);
    end
    

	if  norm(r)/norm(y) < Tolerance  || t > maxiter
        stop = 1; % convergent
    else
        t = t+1;
    end
    

end
% fprintf('\n'); 

% output
out.iter = t;