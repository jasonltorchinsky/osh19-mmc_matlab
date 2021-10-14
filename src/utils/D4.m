% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Calculate second partial derivative in x or y of a 2D
% data set spectrally.

function deriv4=D4(term,derdim1,scalefactor)
  deriv1=D1(term,derdim1,scalefactor);
  deriv2=D1(deriv1,derdim1,scalefactor);
  deriv3=D1(deriv2,derdim1,scalefactor);
  deriv4=D1(deriv3,derdim1,scalefactor);
end
