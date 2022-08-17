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

function  IM_Out     =  Test_GSRC_ENSS (filename, sigma, c, deta, eta)

[Opts, model]  =  Par_Set( sigma,c,deta,eta );

OriName = [filename,'.tif'];

I                =     imread(OriName);

colorI           =      I;

[~, ~, kk]       =     size (I);

if kk==3
    
    I     = rgb2gray (I);
    
    x_yuv = rgb2ycbcr(colorI);
    
    x_denoising_yuv(:,:,2) = x_yuv(:,:,2); % Copy U Componet
    
    x_denoising_yuv(:,:,3) = x_yuv(:,:,3); % Copy V Componet
    
end

Opts.I  =  single(I)/255;

randn('seed',0);
Opts.nim =   Opts.I + Opts.nSig*randn(size(Opts.I));

fprintf('Initial value: PSNR = %2.4f, SSIM = %2.4f \n', csnr( Opts.nim*255, Opts.I*255, 0, 0 ),cal_ssim( Opts.nim*255, Opts.I*255, 0, 0 ));
% GSR_GMM denoising

[IM_Out]  =  GSRC_ENSS_Denoising(Opts,model);

PSNR_Final = csnr( IM_Out*255, Opts.I*255, 0, 0 );
SSIM_Final =  cal_ssim( IM_Out*255, Opts.I*255, 0, 0 );
fprintf('PSNR = %2.4f, SSIM = %2.4f \n', PSNR_Final, SSIM_Final );


Final_denoisng= strcat(filename,'_GSRC_ENSS_','_sigma_',num2str(sigma),'_PSNR_',num2str(PSNR_Final),'_SSIM_',num2str(SSIM_Final),'.png');

imwrite(uint8(IM_Out*255),strcat('./Results/',Final_denoisng));

if (kk==3)
    
x_denoising_yuv(:,:,1) = uint8(IM_Out*255);

x_inpaint_rgb = ycbcr2rgb(uint8(x_denoising_yuv));
        
Final_color_denoisng= strcat(filename,'_GSRC_ENSS_Color_','_sigma_',num2str(sigma),'_PSNR_',num2str(PSNR_Final),'_SSIM_',num2str(SSIM_Final),'.png');

imwrite(uint8(x_inpaint_rgb),strcat('./Results/',Final_color_denoisng));

end



end

