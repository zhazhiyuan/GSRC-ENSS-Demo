% =========================================================================
% GSRC-ENSS Denoising for image denoising, Version 1.0
% Copyright(c) 2017 Zhiyuan Zha
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------

function  [IM_Out] = GSRC_ENSS_Denoising(Opts,model)

IM_Out = Opts.nim;

Opts.nSig0 = Opts.nSig;

[h,  w] = size(IM_Out);

Opts.maxr = h-Opts.ps+1;

Opts.maxc = w-Opts.ps+1;

Opts.maxrc = Opts.maxr * Opts.maxc;

Opts.h = h;

Opts.w = w;

r = 1:Opts.step:Opts.maxr;

Opts.r = [r r(end)+1:Opts.maxr];

c = 1:Opts.step:Opts.maxc;

Opts.c = [c c(end)+1:Opts.maxc];

Opts.lenr = length(Opts.r);

Opts.lenc = length(Opts.c);

Opts.lenrc = Opts.lenr*Opts.lenc;

Opts.ps2 = Opts.ps^2;

AllPSNR = zeros (1,Opts.IteNum);


for ite = 1 : Opts.IteNum

    IM_Out = IM_Out + Opts.delta*(Opts.nim - IM_Out);
    
    % estimation of noise variance
    if ite == 1
        
        Opts.nSig = Opts.nSig0;
        
    else
        
        dif = mean( mean( (Opts.nim - IM_Out).^2 ) ) ;
        
        Opts.nSig = sqrt( abs( Opts.nSig0^2 - dif ) )*Opts.eta;
        
    end
    
    % Generate groups from Noisy image
    [Nonlocal_group,Blk_arr,DC,Opts] = Group_selection( IM_Out, Opts);

    PYZ = zeros(model.nmodels,size(DC,2));
    
    sigma2I = Opts.nSig^2*eye(Opts.ps2);
    
        for i = 1:model.nmodels
            
            sigma = model.covs(:,:,i) + sigma2I;
            
            [R,~] = chol(sigma);
            
            Q = R'\Nonlocal_group;
            
            TempPYZ = - sum(log(diag(R))) - dot(Q,Q,1)/2;
            
            TempPYZ = reshape(TempPYZ,[Opts.nlsp size(DC,2)]);
            
            PYZ(i,:) = sum(TempPYZ);
        end
        % find the best Gaussian component for each group
        [~,dicidx] = max(PYZ);
        
        dicidx = dicidx';
        
        [idx,  s_idx] = sort(dicidx);
        
        idx2 = idx(1:end-1) - idx(2:end);
        
        seq = find(idx2);
        
        seg = [0; seq; length(dicidx)];

    X_re = zeros(Opts.ps2, Opts.maxrc,'single');
    
    W = zeros(Opts.ps2,Opts.maxrc,'single');

    for   j = 1:length(seg)-1
        
        idx =   s_idx(seg(j)+1:seg(j+1));
        
        cls =   dicidx(idx(1));
        
        D   =   Opts.D(:,:, cls);

        for i = 1:size(idx,1)
            
            %Generat Nonlocal similar patches from Noisy image
            Y = Nonlocal_group(:,(idx(i)-1)*Opts.nlsp+1:idx(i)*Opts.nlsp);
            %%
            %Externel Sparsity
            b = D'*Y;  
            
            %%
            %Internel Sparsity
            
            [U,~]  = getsvd(Y);
            
            a  = U'*Y;
            %%
            %Adaptive Regularization Parameterr
            
            s0                 =    a- b;

            s0                 =    mean (s0.^2,2);

            s0                 =    max  (0, s0-Opts.nSig^2);

            lam            =    repmat ( Opts.c1*Opts.nSig^2./(sqrt(s0)+eps),[1 Opts.nlsp]);
            
            %%          
            % Group Sparsity Residual Constraint (GSRC)

            alpha   =   soft( a-b, lam ) + b;
            
            % add DC components and aggregation
           
            X_re(:,Blk_arr(:,idx(i))) = X_re(:,Blk_arr(:,idx(i)))+bsxfun(@plus,U*alpha, DC(:,idx(i)));

            W(:,Blk_arr(:,idx(i))) = W(:,Blk_arr(:,idx(i)))+ones(Opts.ps2,Opts.nlsp);
        end
    end
    
    IM_Out = zeros(h,w,'single');
    
    im_wei = zeros(h,w,'single');
    
    r = 1:Opts.maxr;
    
    c = 1:Opts.maxc;
    
    k = 0;
    
    for i = 1:Opts.ps
        for j = 1:Opts.ps
            k = k+1;
            IM_Out(r-1+i,c-1+j)  =  IM_Out(r-1+i,c-1+j) + reshape( X_re(k,:)', [Opts.maxr Opts.maxc]);
            im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + reshape( W(k,:)', [Opts.maxr Opts.maxc]);
        end
    end
    
    IM_Out  =  IM_Out./im_wei;
    
    PSNR =   csnr( IM_Out*255, Opts.I*255, 0, 0 );
    
    AllPSNR(ite) = PSNR;
    
    SSIM      =  cal_ssim( IM_Out*255, Opts.I*255, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n',ite, PSNR,SSIM);
    
    if ite>1
        if AllPSNR(ite)-AllPSNR(ite-1)<0.005
            break;
        end
    end
    
end

IM_Out(IM_Out > 1) = 1;
IM_Out(IM_Out < 0) = 0;
return;