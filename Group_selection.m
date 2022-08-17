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
function       [Nonlocal_group,Blk_arr,DC,Opts] = Group_selection( IM, Opts)
% record the non-local patch set and the index of each patch in
% of seed patches in image
IM = single(IM);
X   = zeros(Opts.ps2, Opts.maxrc, 'single');
k   = 0;
for i = 1:Opts.ps
    for j = 1:Opts.ps
        k = k+1;
        blk = IM(i:end-Opts.ps+i,j:end-Opts.ps+j);
        X(k,:) = blk(:)';
    end
end
% index of each patch in image
Index    =   (1:Opts.maxrc);
Index    =   reshape(Index,Opts.maxr,Opts.maxc);
% record the indexs of patches similar to the seed patch
Blk_arr  =  zeros(Opts.nlsp, Opts.lenrc ,'single');
% Patch Group Means
DC = zeros(Opts.ps2,Opts.lenrc,'single');
% non-local patch groups
Nonlocal_group = zeros(Opts.ps2,Opts.lenrc*Opts.nlsp,'single');
for  i  =  1 :Opts.lenr
    for  j  =  1 : Opts.lenc
        row = Opts.r(i);
        col = Opts.c(j);
        off = (col-1)*Opts.maxr + row;
        off1 = (j-1)*Opts.lenr + i;
        % the range indexes of the window for searching the similar patches
        rmin = max( row - Opts.Win, 1 );
        rmax = min( row + Opts.Win, Opts.maxr );
        cmin = max( col - Opts.Win, 1 );
        cmax = min( col + Opts.Win, Opts.maxc );
        idx     =   Index(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        neighbor = X(:,idx); % the patches around the seed in X
        seed  = X(:,off);
        dis = sum(bsxfun(@minus,neighbor, seed).^2,1);
        [~,ind] = sort(dis);
        indc = idx( ind( 1:Opts.nlsp ) );
        Blk_arr(:,off1) = indc;
        temp = X( : , indc );
        DC(:,off1) = mean(temp,2);
        Nonlocal_group(:,(off1-1)*Opts.nlsp+1:off1*Opts.nlsp) = bsxfun(@minus,temp,DC(:,off1));
    end
end