function PUR_prog2PE(beta)
%calculates power curves PE CO=t LO and critical values CV

%load inputs PUR_prog1J
load INPUTS.mat stepsize hmaxPE stepsizebeta betamax alpha
if nargin<1
    beta=4;
end
if or(beta<0,beta>betamax); stop; end
if abs(beta/stepsizebeta-round(beta/stepsizebeta))>1e-8
    betanew=round(beta/stepsizebeta)*stepsizebeta;
    fprintf('WARNING beta projected on betagrid: %10.6f to %10.6f\n',beta,betanew)
    beta=betanew;
end

sizehgrid=round(hmaxPE/stepsize+2); % 0:stepsize:hmaxPE + extreme point
extr=125/hmaxPE;
betagridpoint=round(beta/stepsizebeta+1); 

%load and calculate J
load J.mat J
J0=1/2; 
J=J0+beta*J;

%calculating critical values of point optimal tests
cv=zeros(sizehgrid,1);
c=0;
for hh=1:sizehgrid
    if hh<sizehgrid
        h=(hh-1)*stepsize;
    else
        h=extr*hmaxPE;
    end
    while mean(normcdf(-c./sqrt((h^2+(h==0))*J)-.5*sqrt(h^2*J)))>alpha
        c=c+100;
    end
    while mean(normcdf(-c./sqrt((h^2+(h==0))*J)-.5*sqrt(h^2*J)))<alpha
        c=c-100;
    end
    level=mean(normcdf(-c./sqrt((h^2+(h==0))*J)-.5*sqrt(h^2*J)));
    dc=100;
    while or(dc>1e-8,level-alpha>1e-8)
        dc=dc/2;
        levelnew=mean(normcdf(-(c+dc)./sqrt((h^2+(h==0))*J)-.5*sqrt(h^2*J)));
        if levelnew>=alpha
            level=levelnew;
            c=c+dc;
        end
    end
    cv(hh)=c;
end

%calculating PE CO LO power curves
power=zeros(sizehgrid,3);
for hh=1:sizehgrid
    if hh<sizehgrid
        h=(hh-1)*stepsize;
    else
        h=extr*hmaxPE;
    end
    c=cv(hh);
    powenv=mean(normcdf(-c./sqrt((h^2+(h==0))*J)+.5*sqrt(h^2*J)));
    powttest=mean(normcdf(norminv(alpha)+sqrt(h^2*J)));
    c=cv(1);
    powlocopt=mean(normcdf(-c./sqrt(J)+sqrt(h^2*J)));
    power(hh,:)=[powenv powttest powlocopt]; 
end

load PECOLOCVbeta.mat PE CO LO CV
PE(betagridpoint,:)=[beta,power(:,1)'];
CO(betagridpoint,:)=[beta,power(:,2)'];
LO(betagridpoint,:)=[beta,power(:,3)'];
CV(betagridpoint,:)=[beta,cv'];
% fprintf('%10.6f',PE(betagridpoint,:))
% fprintf('\n')
save PECOLOCVbeta.mat PE CO LO CV

end