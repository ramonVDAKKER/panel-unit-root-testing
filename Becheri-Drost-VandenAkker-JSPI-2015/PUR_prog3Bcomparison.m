function PUR_prog3Bcomparison(beta,hmax,hmin,delta)
%calculates preliminary power curve PO and some top functions

%load inputs PUR_prog1J
load INPUTS.mat stepsize hmaxPE stepsizebeta betamax
if nargin<1
    beta=4;
end
if nargin<2
    hmax=hmaxPE/2;
end
if nargin<3
    hmin=0;
end
if nargin<4
    delta=10;
end
if or(beta<0,beta>betamax); stop; end
if abs(beta/stepsizebeta-round(beta/stepsizebeta))>1e-8
    betanew=round(beta/stepsizebeta)*stepsizebeta;
    fprintf('WARNING beta projected on betagrid: %10.6f to %10.6f\n',beta,betanew)
    beta=betanew;
end

hgridpointL=round(hmin/stepsize+1);
hgridpointR=round(hmax/stepsize+1);
sizehgrid=round(hmaxPE/stepsize+2); % 0:stepsize:hmaxPE + extreme point
extr=125/hmaxPE;
betagridpoint=round(beta/stepsizebeta+1); 

%load and calculate J
load J.mat J
J0=1/2; 
J=J0+beta*J;

%load power curves PE CO=t LO and CV
load PECOLOCVbeta.mat PE CV

%determine critical values of PO power curves to be investigated
cv=CV(betagridpoint,1+(hgridpointL:hgridpointR));

%power curve of PE
power=PE(betagridpoint,1+(1:sizehgrid))';

%evaluate all point optimal tests to be investigated and determine left/right tops
powpointoptuu=zeros(sizehgrid,1);
powpointopt=zeros(sizehgrid,1);
maxdiffL=zeros(sizehgrid-1,1);
maxdiffR=zeros(sizehgrid-1,1);
h00opt=-99;
diffh00opt=99;
for uu=1:delta:hgridpointR-hgridpointL+1
    h00=(hgridpointL-1+uu-1)*stepsize;
    for hh=1:sizehgrid
        if hh<sizehgrid
            h=(hh-1)*stepsize;
        else
            h=extr*hmaxPE;
        end
        powpointoptuu(hh)=mean(normcdf(-cv(uu)./sqrt((h00^2+(h00==0))*J)+(h-.5*h00)*sqrt(J)));
    end
    maxdiffuu=max(-powpointoptuu+power(:,1));
    maxdiffL(hgridpointL-1+uu)=PUR_max(-powpointoptuu(1:hgridpointL-1+uu)+power(1:hgridpointL-1+uu,1));
    maxdiffR(hgridpointL-1+uu)=PUR_max(-powpointoptuu(hgridpointL-1+uu:sizehgrid)+power(hgridpointL-1+uu:sizehgrid,1));
    if maxdiffuu<diffh00opt
        h00opt=h00;
        diffh00opt=maxdiffuu;
        powpointopt=powpointoptuu;
    end;
end

load PObeta2.mat PO TOPL TOPR h0
PO(betagridpoint,:)=[beta,powpointopt'];
TOPL(betagridpoint,:)=[beta,maxdiffL'];
TOPR(betagridpoint,:)=[beta,maxdiffR'];
h0(betagridpoint,:)=[beta,h00opt];
% fprintf('%10.4f',PO(betagridpoint,:))
% fprintf('\n')
save PObeta2.mat PO TOPL TOPR h0

if or(h00opt<=hmin+2*stepsize,h00opt>=hmax-2*stepsize)
    fprintf('***** WARNNG h0 near or on boundary *****: %10.4f %10.4f %10.4f\n',[hmin,h00opt,hmax])
end

end