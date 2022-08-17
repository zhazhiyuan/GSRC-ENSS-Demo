function  [Opts, model]  =  Par_Set( nsig,c,deta,eta  )

Opts.step = 3;       
Opts.IteNum = 100;  
Opts.nSig      =   nsig/255;

if nsig <= 10
    load './model/Kodak_6x6_step4_win50_similar80_delta0.004_cls65.mat';
    Opts.c1 = c*2*sqrt(2);
    Opts.delta =deta;
    Opts.eta=eta;
elseif nsig<=20
    load './model/Kodak_6x6_step4_win50_similar80_delta0.004_cls65.mat';
    Opts.c1 = c*2*sqrt(2);
    Opts.delta = deta;
    Opts.eta= eta;
elseif nsig <=30
    load './model/Kodak_7x7_step4_win50_similar90_delta0.004_cls65.mat';    
    Opts.c1 = c*2*sqrt(2);
    Opts.delta = deta;
    Opts.eta=eta;
elseif nsig<=40
    load './model/Kodak_8x8_step4_win50_similar100_delta0.005_cls65.mat';
    Opts.c1 = c*2*sqrt(2);
    Opts.delta = deta; 
    Opts.eta=eta; 
elseif nsig<=50
    load './model/Kodak_8x8_step4_win50_similar100_delta0.008_cls65.mat';
    Opts.c1 = c*2*sqrt(2);
    Opts.delta = deta; 
    Opts.eta=eta;  
elseif nsig<=75
    load './model/Kodak_9x9_step4_win50_similar120_delta0.016_cls65.mat';
    Opts.c1 = c*2*sqrt(2);
    Opts.delta = deta;  
    Opts.eta=eta;
else
    load './model/Kodak_9x9_step4_win50_similar120_delta0.016_cls65.mat'; 
    Opts.c1 = c*2*sqrt(2); 
    Opts.delta = deta;
    Opts.eta=eta;   
end
Opts.ps = ps;        % patch size
Opts.nlsp = nlsp;  % number of non-local patches
Opts.Win = win;   % size of window around the patch
% dictionary and regularization parameter
for i = 1:size(GMM_D,2)
    Opts.D(:,:,i) = reshape(single(GMM_D(:, i)), size(GMM_S,1), size(GMM_S,1));
end

