% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Prepare a linear solution for use as initial
% condition in nonlinear script.  

% Steps that will need to be done:
% Get 1D linear solution in Fourier space. - Need to move this to the right
% folder.

disp(strcat(['Calculating linear mode #',int2str(num_damped),', k=',int2str(wave_number)]));

% disp(strcat(['tauvec=',num2str(tauvec)]));
% disp(strcat(['Bvsdim=',num2str(Bvsdim)]));
% disp(strcat([']));
% disp(strcat([' ]));
tauvec;

M=numy; % in linear report
J=numlevs; % in linear report

% Zonal wavenumber you want to study
kbar_int=wave_number;

% % Dimensional kbar
kbar = kbar_int*2*pi/Pe;

% Make Y matrix that you'll need
Y = diag(yy);

% Make F (called Fmat to avoid confusion with flux vector F in nonlinear
% model)
Fmat=fft(eye(numy))/sqrt(M);
Fmatinv=ifft(eye(numy))*sqrt(M);

% % TRYING THE SINE TRANSFORM HERE!!
% Fmat=dst(eye(numy));
% Fmatinv=idst(eye(numy));

% Make Dx and Dxx
Dhatx = 1i*kbar*eye(numy);
Dhatxx = -kbar^2*eye(numy);

% Make Dy and Dyy
wavenumbery = [0:numy/2 -numy/2+1:-1];
wndimy = wavenumbery*2*pi/(2*pY);
Dhaty=diag(1i*wndimy);
Dhaty(M/2+1,M/2+1)=0;
Dhatyy=diag(-wndimy.^2);
Dy = Fmatinv*Dhaty*Fmat;
Dyf=diag(Dy*beta*yy');
invLap=inv(Dhatxx+Dhatyy);
Lap=Dhatxx+Dhatyy;

% Make Pu, Pv
% NEED TO FIX THE LAST TERM IN EACH OF THESE - DY IS NOT CORRECT HERE!!
% Pu = -beta*inv(Dhatxx+Dhatyy)*(eye(numy)+Fmat*Y*Dy*Fmatinv);
% Pu = (-beta*inv(Dhatxx+Dhatyy)*(Fmat*Dy*Y*Fmatinv+Fmat*Y*Dy*Fmatinv));
% Pv = (beta*inv(Dhatxx+Dhatyy)*(1i*kbar*Fmat*Y*Fmatinv));
% Bu = -Fmat*Dy*Fmatinv*inv(Dhatxx+Dhatyy);
% Bv = 1i*kbar*inv(Dhatxx+Dhatyy);

% Create block matrix A
numeqs=(4*numlevs-3)*numy;
A=zeros(numeqs,numeqs);

% Populate A except for last two row and column blocks
for Jcounti=1:J-1
    for Jcountj=1:J-1

        % Populate the U parts of the matrix Aij
        Uuij=zeros(numy,numy);
        if Jcounti==Jcountj
            Uuij=Uuij+1/tau_u*eye(numy);
        end
        if Jcounti==Jcountj
            Uvij=-beta*Fmat*Y*Fmatinv;
        else
            Uvij=zeros(numy,numy);
        end
        if Jcounti>Jcountj
            Utij=g*dz*1i*kbar*Jcountj/(theta0*J)*eye(numy);
        else
            Utij=-g*dz*1i*kbar*(J-Jcountj)/(theta0*J)*eye(numy);
        end
        Uqij=zeros(numy,numy);
        
        % Populate the V parts of the matrix Aij
        if Jcounti==Jcountj
            Vuij=beta*Fmat*Y*Fmatinv;
        else
            Vuij=zeros(numy,numy);
        end
        Vvij=zeros(numy,numy);
        if Jcounti==Jcountj
            Vvij=Vvij+1/tau_u*eye(numy);
        end
        if Jcounti>Jcountj
            Vtij=g*dz*Jcountj/(theta0*J)*Dhaty;
        else
            Vtij=-g*dz*(J-Jcountj)/(J*theta0)*Dhaty;
        end
        Vqij=zeros(numy,numy);
        
        % Populate the Theta parts of the matrix Aij
        if Jcounti>=Jcountj
            Tuij=-1i*kbar*Bdim*dz*eye(numy);
        else
            Tuij=zeros(numy,numy);
        end
        if Jcounti>=Jcountj
            Tvij=-Bdim*dz*Dhaty;
        else
            Tvij=zeros(numy,numy);
        end
        Ttij=zeros(numy,numy);
        if Jcounti==Jcountj
            Ttij=Ttij+1/tau_theta*eye(numy);
        end
        mode2coeff=0.5;
        if use_2BCConvAdj==1
            Tqij=-sqrt(2)*Lv/cp*(1/tau_upper*(sqrt(3)/J*sin(zvec(Jcountj)*pi/H)-4*sqrt(3)/J*sin(2*zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)-mode2coeff*sin(2*zvec(Jcounti)*pi/H))...
                +1/tau_lower*(sqrt(3)/J*sin(zvec(Jcountj)*pi/H)+4*sqrt(3)/J*sin(2*zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)+mode2coeff*sin(2*zvec(Jcounti)*pi/H)))*eye(numy);
        elseif use_2BCConvAdj==2
            Tqij=-sqrt(2)*Lv/cp*(1/tau_upper*(sqrt(3)/J*sin(zvec(Jcountj)*pi/H)-4*sqrt(3)/J*sin(2*zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)-mode2coeff*sin(2*zvec(Jcounti)*pi/H))...
                +1/tau_lower*(2/J*sin(zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)))*eye(numy);
        else
            Tqij=zeros(numy,numy);
            if Jcounti==Jcountj
                Tqij=Tqij-Lv/cp/tauvec(Jcounti)*eye(numy);
            end
        end
        
        
        % Populate the Q parts of the matrix Aij
        if Jcounti>=Jcountj
            Quij=-1i*kbar*Bvsvec(Jcounti)*Fmat*diag(BvsYvec)*Fmatinv*dz*eye(numy);
        else  
            Quij=zeros(numy,numy);
        end
        
%         if Jcounti>=Jcountj
%             Qvij=-Bvsvec(Jcounti)*dz*Dhaty;
%         else
%             Qvij=zeros(numy,numy);
%         end
        if Jcounti>=Jcountj
            Qvij=-Bvsvec(Jcounti)*Fmat*diag(BvsYvec)*Fmatinv*dz*Dhaty;
        else
            Qvij=zeros(numy,numy);
        end
        if Jcounti==Jcountj
            Qvij=Qvij+1/2*Fmat*diag(Dy*QBGvec(Jcounti)*transpose(BvsYvec))*Fmatinv;
        end
        if Jcounti==Jcountj-1
            Qvij=Qvij+1/2*Fmat*diag(Dy*QBGvec(Jcounti)*transpose(BvsYvec))*Fmatinv;
        end
        if Jcounti==J-1
            Qvij=Qvij-1/2*Fmat*diag(Dy*QBGvec(Jcounti)*transpose(BvsYvec))*Fmatinv;
        end
        Qtij=zeros(numy,numy);
        if use_2BCConvAdj==1
            Qqij=sqrt(2)*(1/tau_upper*(sqrt(3)/J*sin(zvec(Jcountj)*pi/H)-4*sqrt(3)/J*sin(2*zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)-mode2coeff*sin(2*zvec(Jcounti)*pi/H))...
                +1/tau_lower*(sqrt(3)/J*sin(zvec(Jcountj)*pi/H)+4*sqrt(3)/J*sin(2*zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)+mode2coeff*sin(2*zvec(Jcounti)*pi/H)))*eye(numy);
            if Jcounti==Jcountj
                Qqij=Qqij+kbar^2*bvec(Jcounti)*eye(numy)-bvec(Jcounti)*Dhatyy+2*Dvvec(Jcounti)/dz^2*eye(numy);
            elseif Jcounti==Jcountj-1 || Jcounti==Jcountj+1
                Qqij=Qqij-Dvvec(Jcounti)/dz^2*eye(numy);
            end
        elseif use_2BCConvAdj==2
            Qqij=sqrt(2)*(1/tau_upper*(sqrt(3)/J*sin(zvec(Jcountj)*pi/H)-4*sqrt(3)/J*sin(2*zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)-0.5*sin(2*zvec(Jcounti)*pi/H))...
                +1/tau_lower*(2/J*sin(zvec(Jcountj)*pi/H))...
                *(sin(zvec(Jcounti)*pi/H)))*eye(numy);
            if Jcounti==Jcountj
                Qqij=Qqij+kbar^2*bvec(Jcounti)*eye(numy)-bvec(Jcounti)*Dhatyy+2*Dvvec(Jcounti)/dz^2*eye(numy);
            elseif Jcounti==Jcountj-1 || Jcounti==Jcountj+1
                Qqij=Qqij-Dvvec(Jcounti)/dz^2*eye(numy);
            end
        else
            if Jcounti==Jcountj
                Qqij=1/tauvec(Jcounti)*eye(numy)+kbar^2*bvec(Jcounti)*eye(numy)-bvec(Jcounti)*Dhatyy+2*Dvvec(Jcounti)/dz^2*eye(numy);
            elseif Jcounti==Jcountj-1 || Jcounti==Jcountj+1
                Qqij=-Dvvec(Jcounti)/dz^2*eye(numy);
            else
                Qqij=zeros(numy,numy);
            end 
        end
%         if Jcounti==Jcountj
%             Qqij=1/tauvec(Jcounti)*eye(numy)+kbar^2*bvec(Jcounti)*eye(numy)-bvec(Jcounti)*Dhatyy+2*Dvvec(Jcounti)/dz^2*eye(numy);
%         elseif Jcounti==Jcountj-1 || Jcounti==Jcountj+1
%             Qqij=-Dvvec(Jcounti)/dz^2*eye(numy);
%         else
%             Qqij=zeros(numy,numy);
%         end
        
        % Put all the parts together to make Aij, then populate A
        % All block rows and block columns except last two
        Aij = [Uuij Uvij Utij Uqij;...
               Vuij Vvij Vtij Vqij;...
               Tuij Tvij Ttij Tqij;...
               Quij Qvij Qtij Qqij];
        A(1+(Jcounti-1)*4*M:Jcounti*4*M,1+(Jcountj-1)*4*M:Jcountj*4*M)=Aij;
    end
end

% Last block row (not including last block column)
A(1+(J-1)*4*M:end,1:(J-1)*4*M)=zeros(M,4*M*(J-1));

% phatbar=invLap*Fmat*(Dyf*Dy*Fmatinv*invLap*Fmat-Y)*Fmatinv;
phatbar=invLap*Fmat*(Dyf*Dy*Fmatinv*invLap*Fmat+beta*Y)*Fmatinv;

% Last block column (not including last block row)
for Jcounti=1:J-1
%     Uzij=1i*kbar*(Pu*Bu+Pv*Bv);
%     Uzij=1i*kbar*phatbar;
    Uzij=zeros(M,M); 
%     Vzij=Dhaty*(Pu*Bu+Pv*Bv);
%     Vzij=Dhaty*phatbar;
    Vzij=zeros(M,M); 
    Tzij=zeros(M,M);
    Qzij=zeros(M,M);
    if Jcounti==J-1
        Qzij=Qzij+1/2*Fmat*diag(Dy*QBGvec(Jcounti)*transpose(BvsYvec))*J*Fmatinv*Dhatx*invLap;
    end
    Aij=[Uzij;...
         Vzij;...
         Tzij;...
         Qzij];
    A(1+(Jcounti-1)*4*M:Jcounti*4*M,1+(J-1)*4*M:end)=Aij;
end

% Last block row and column
% Calculate derivative spectrally
Zzij = 1i*kbar*Fmat*Dyf*Fmatinv*invLap+1/tau_u*eye(numy);
% Calculate derivative analytically 
% Zzij = 1i*kbar*beta*invLap;
% Zzij = Zzij-Dhu*Lap;

A(1+(J-1)*4*M:end,1+(J-1)*4*M:end)=Zzij;


% Find eigenvalues of A
[evecsFsp evalsmat]=eig(A);
% clear A Uuij Uvij Utij Uqij Vuij Vvij Vtij Vqij Tuij Tvij Ttij Tqij Quij Qvij Qtij Qqij;
% clear Uzij Vzij Tzij Qzij Zzij
evals=-1i*diag(evalsmat)*3600*24; % evals now omega in 1/day units

sortbyspeed=0; 
% Positive - eastward;  negative - westward

% % Sort the eigenvalues according to speed/least damped. 
if sortbyspeed==1 [~,Index_Sorted] = sort(real(evals),'descend');... % Fastest eastward speed first
else [~,Index_Sorted] = sort(imag(evals),'descend'); % largest growth rates first
end

evalssort=evals(Index_Sorted);

% evalssort(1:5)

% Set sorted growthrates
Growthrates = imag(evalssort);
disp('Growthrates');
disp(num2str(Growthrates(1:5)));
 
% Set sorted frequencies
Frequencies = real(evalssort)/(2*pi);

% Set sorted Right Evecs
evecsFspsort = evecsFsp(:,Index_Sorted);

% Now make a copy with the top-level winds replacing the barotropic
% relative vorticity
evecsphspYwinds=zeros(M*(4*J-3),1);

evectouse=num_damped;
% Inverse FFT to get back to physical space in y
for Jcounti=1:4*J-3
    evecsphspYwinds(1+(Jcounti-1)*M:Jcounti*M,1)=(Fmatinv*evecsFspsort(1+(Jcounti-1)*M:Jcounti*M,evectouse));
end

% Convert eigenvector to physical space in x
evectoshow=zeros(M*(4*J-3),numx);
for rowcount=1:M*(4*J-3)
    evectoshow(rowcount,:)=real(evecsphspYwinds(rowcount,1)*exp(1i*kbar*xx));
end

% Construct *baroclinic* solution in physical x-y space
umatini=zeros(numy,numx,numlevs+1);
vmatini=zeros(numy,numx,numlevs+1);
zmatini=zeros(numy,numx,1);
thetamatini=zeros(numy,numx,numlevs+1);
qmatini=zeros(numy,numx,numlevs+1);

for zcount=2:numlevs
    umatini(:,:,zcount)=evectoshow(1+4*M*(zcount-2):M+4*M*(zcount-2),:);
    vmatini(:,:,zcount)=evectoshow(M+1+4*M*(zcount-2):2*M+4*M*(zcount-2),:);
    thetamatini(:,:,zcount)=evectoshow(2*M+1+4*M*(zcount-2):3*M+4*M*(zcount-2),:);
    qmatini(:,:,zcount)=evectoshow(3*M+1+4*M*(zcount-2):4*M+4*M*(zcount-2),:);
end
zmatini(:,:,1)=evectoshow(M*4*(J-1)+1:M*4*(J-1)+M,:);

% disp(strcat(['Finished linear stability analysis with k=',int2str(wave_number)]));

clear evecsFsp evalsmat evals evecsFspsort evecsphspYwinds 

