% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Calculate first partial derivative in x or y of a 2D
% data set spectrally.  

% NOTE:  ASSUMES Y IS IN THE UP-DOWN (ROWS) DIRECTION, X IS IN THE
% LEFT-RIGHT (COLUMNS) DIRECTION!!!

function deriv1=D1(term,derdim,scalefactor)
    [numy,numx]=size(term);
    if derdim=='y'; dim=1;
    elseif derdim=='x'; dim=2;
    end
    fftprod=fft(term,[],dim);
    if derdim=='y'
        wavenumbervec=transpose([0:numy/2-1 0 -numy/2+1:-1]);
        for kcount=1:numx
            fftprod(:,kcount)=fftprod(:,kcount).*1i.*wavenumbervec.*scalefactor;
        end
    elseif derdim=='x'
        wavenumbervec=[0:numx/2-1 0 -numx/2+1:-1];
        for lcount=1:numy
            fftprod(lcount,:)=fftprod(lcount,:).*1i.*wavenumbervec.*scalefactor;
        end
    end
    deriv1=real(ifft(fftprod,[],dim));
end
