% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Calculate second partial derivative in x or y of a 2D
% data set spectrally.

function deriv2=D2(term,derdim1,derdim2,scalex,scaley)
  if derdim1=='x'; scalefactor=scalex; elseif derdim1=='y'; scalefactor=scaley; end   
  deriv1=D1(term,derdim1,scalefactor);
  if derdim2=='x'; scalefactor=scalex; elseif derdim2=='y'; scalefactor=scaley; end   
  deriv2=D1(deriv1,derdim2,scalefactor);
end
