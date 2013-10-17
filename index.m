%%
function index

% options: 'gm'=1 'gk'=2 'bm'=3 'bk'=4
% g: gaussian 
% b: binary 
% m: fixed measurements 
% k: fixed sparsity level
clear all
k = 100;
m = 500;
n = 2000;
ensemble = 'USE';
opts.sigma = 0;
opts.gauss = 1;

% [x, y, A] = gen_signal(k, m, n, ensemble, opts);
[x, y, A] = gen_signal(k, m, n, ensemble, opts);

% 收敛误差
tol = 1e-6;

% tic
% [xt out] = IDET(A, y, k);
% toc
% out.iter
% SupportDetection(x, xt)
% 
% 
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% tic
% [xt out] = IDET(H, y, k, 'At', Ht, 'MaxIt', 20);
% toc
% out.iter
% SupportDetection(x, xt)





% [xt out] = IDETbeta(A, y, 0.7);
% out.iter
% SupportDetection(x, xt)
% 
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt out] = IDETbeta(H, y, 0.7, 'At', Ht, 'MaxIt', 20);
% out.iter
% SupportDetection(x, xt)




% tic
% [xt out] = IDETgamma(A, y, 0.7);
% toc
% out.iter
% SupportDetection(x, xt)
% 
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% tic
% [xt out] = IDETgamma(H, y, 0.7, 'At', Ht, 'MaxIt', 20);
% toc
% out.iter
% SupportDetection(x, xt)





% [xt out] = SWGP(A, y, 0.7);
% out.iter
% SupportDetection(x, xt)
% 
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt out] = SWGP(H, y, 0.7, 'At', Ht, 'MaxIt', 20);
% out.iter
% SupportDetection(x, xt)




%  xt = GDE(A, y, k );

% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% z_init = zeros(n, 1);
% [xt,A_index,Count,GML_record] = AHT(H,Ht,y,z_init,1,1);
% Count
% SupportDetection(x, xt)

% tic
% Len_thresh=1;          
% z_init=zeros(m,1);
% [xt,A_index_ADORE,Count_ADORE,Golden_Iter,USS]=ADORE(A,[],y,'SearchLen',Len_thresh,'Thresh',tol,'IsHOrthonormal',0);
% SupportDetection(x, xt)
% toc
% 
% 
% tic
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% Len_thresh=1;          
% z_init=zeros(m,1);
% [xt,A_index_ADORE,Count_ADORE,Golden_Iter,USS]=ADORE(H,Ht,y,'SearchLen',Len_thresh,'Thresh',tol,'IsHOrthonormal',0);
% SupportDetection(x, xt)
% toc


% tic
% [xt, A_index,Count,loglikelihood,delta2] =  DORE(A, [], y, k, 'Thresh',tol,'visibility',0,'IsHOrthonormal',0);
% Count
% SupportDetection(x, xt)
% toc
% 
% tic
% invAAt      = inv(A*A');
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% invHHt      = @(z) invAAt*z;
% [xt, A_index,Count,loglikelihood,delta2] =  DORE(H, Ht, y, k, 'Thresh',tol,'Inverse',invHHt,'visibility',0,'IsHOrthonormal',0);
% Count
% SupportDetection(x, xt)
% toc
% 
% tic
% [xt, Out] = TEMP(A, y, 1, 'Tolerance', tol);
% Out.iter
% SupportDetection(x, xt)
% toc

% 
% tic
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt, Out] = TEMP(H, y, 1, 'At', Ht, 'MaxIt',100);
% Out.iter
% SupportDetection(x, xt)
% toc

% tic
% [xt, Out] = EMTP(A, y, k, 1, 'Tolerance',tol);
% Out.iter
% SupportDetection(x, xt)
% toc

% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt, Out] = EMTP(H, y, k, 1, 'At', Ht, 'MaxIt',10);
% Out.iter
% SupportDetection(x, xt)


% tic
% [xt, Out] = EMTPbeta(A, y, 0.5, 'Tolerance',tol);
% Out.iter
% SupportDetection(x, xt)
% toc

% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt, Out] = EMTPbeta(H, y, 0.5, 'At', Ht, 'MaxIt',100);
% Out.iter
% SupportDetection(x, xt)

% tic
% [xt, Out] = EMTPgamma(A, y, 0.8, 'Tolerance',tol);
% Out.iter
% SupportDetection(x, xt)
% toc

% SupportDetection(x, xt)
% % SupportDetection(x, xbar) % 对比原始信号与重构信号



