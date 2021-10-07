% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Dealias 3D(x,y,z) datasets in x and y by smoothing a
% fraction of the 2D Fourier coefficients to zero.  

function [output]=dealias(term,fraction)
  [numy,numx,numz]=size(term);
  output=zeros(numy,numx,numz);
  m1y=ceil(fraction*numy/2)+1;
  m2y=numy-m1y+2;
  m1x=ceil(fraction*numx/2)+1;
  m2x=numx-m1x+2;
  for zcount=1:numz
      termtemp=squeeze(term(:,:,zcount));
      fft2term=fft2(termtemp);
      fft2term(m1y:m2y,:)=zeros(m2y-m1y+1,numx);
      fft2term(:,m1x:m2x)=zeros(numy,m2x-m1x+1);
      output(:,:,zcount)=ifft2(fft2term);
  end

end

