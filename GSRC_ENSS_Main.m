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

clc
clear
close all;

for i =1:14
     
ImgNo =i;

switch ImgNo
            case 1
                filename = 'House256';
            case 2
                filename = 'Leaves256';
            case 3
                filename = 'Monarch256';
            case 4
                filename = 'airplane256';    
            case 5
                filename = 'bridge256';
            case 6
                filename = 'elaine256';
            case 7
                filename = 'fingerprint256';
            case 8
                filename = 'J.Bean_gray';    
            case 9
                filename = 'Lake';
            case 10
                filename = 'Miss_gray';                 
            case 11
                filename = 'Parthenon_gray';
            case 12
                filename = 'starfish256';    
            case 13
                filename = 'straw';
            case 14
                filename = 'foreman256';                                
end
        
for j  =  2
           
    .........................................00
 SigNum=[10,20,30,40,50,75,100];  
 sigma=SigNum(j);
 
 if j==1 %10
    c=0.14; 
    deta=0.19; 
    eta=1.08;
 elseif j==2%20
    c=0.13; 
    deta=0.2; 
    eta=1.05;    
 elseif j==3%30
    c=0.12; 
    deta=0.21; 
    eta=1.05;     
 elseif j==4%40
     c=0.11; 
    deta=0.22; 
    eta=1.05;    
 elseif j==5%50
     c=0.1; 
    deta=0.23; 
    eta=1.05;    
 elseif j==6%75
     c=0.09; 
    deta=0.24; 
    eta=1;    
 else%100
    c=0.08; 
    deta=0.25; 
    eta=1;    
 end
 
 fprintf('........The i=%d image..........\n',i)
 fprintf(filename);
 fprintf('........Sigma=%d image..........\n',sigma)
 fprintf('.................\n');
   
  im_out     =  Test_GSRC_ENSS (filename, sigma, c, deta, eta);
 
  imshow(im_out,[]);
 
end

end


         