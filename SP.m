function [xt Out] = SP(A, y, k, varargin)
% SP: Subspace pursuit
%----- input 
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
%----- Output
% xt: the recovery signal
%
% Written by Heping Song

path(path, './subfunctions');

% A is a matrix or function handle 
explicitA = ~(ischar(A) || isa(A, 'function_handle'));

MaxIt = 100; % maximum number of iterations 
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


% initialization
t = 1; % iteration number
xt = zeros(n,1);  % initial x
% xt = pinv(A)*y; % initial x
r = y;  % initial residue
support_size = k;  % initial sparsity level
support = [];
stop = 0; % not convergent

while ~stop
    
    if (explicitA)
    % support detection
    [v, idx] = sort(abs( A' *r), 'descend');
	activeset = sort(idx(1:support_size));
    support =sort(union(activeset, support));
    
    % signal estimation
    At = A(:, support);
    z = abs(pinv(At) * y);
    [v, idx] = sort(z,'descend');
    support = support(idx(1:support_size));
    
    At = A(:, support);
    xt = zeros(n,1);
    xt(support) = pinv(At) * y;
    r = y - A*xt;
    
    else
    % support detection
    [v idx] = sort(abs(At(r)), 'descend');
    activeset = sort(idx(1:support_size));
    support =sort(union(activeset, support));

    
    % signal estimation
	xt = zeros(n,1);
	[xt r] =MySubsetCG(y, xt, A, At, support, 1e-16, 0, MaxIt); 
    [v, idx] = sort(abs(xt),'descend');
    support = idx(1:support_size); 

    [xt r] =MySubsetCG(y, xt, A, At, support, 1e-16, 0, MaxIt); 
	r = y - A(xt);
    end
    
    
	if ( norm(r)/norm(y) < Tolerance  || t > MaxIt )
        stop = 1; % convergent
	else
        t = t+1;
    
    end
end

Out.iter = t;

