function PUR_prog6plots(beta)
%calculates distance meaures and plots figures

%load inputs PUR_prog1J
load INPUTS.mat stepsize hmaxPE stepsizebeta betamax
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
betagridpoint=round(beta/stepsizebeta+1); 

%load power curves PE CO=t LO and CV
load PECOLOCVbeta.mat PE CO LO
%load power curve PO comparisons
load POastbeta.mat POast hast
%load power curve top information
load PObeta2.mat TOPL TOPR
TOPLhulp=TOPL; 
TOPRhulp=TOPR; 
load PObeta.mat TOPL TOPR
TOPL=max(TOPLhulp,TOPL); 
TOPR=max(TOPRhulp,TOPR); 

figure(1)
subplot(1,2,1)
hold on
hgrid=stepsize*(0:sizehgrid-2);
plot(hgrid,[POast(betagridpoint,1+(1:sizehgrid-1));CO(betagridpoint,1+(1:sizehgrid-1));LO(betagridpoint,1+(1:sizehgrid-1))])
plot(hgrid,PE(betagridpoint,1+(1:sizehgrid-1)),'-k')
plot([0 hmaxPE],[PUR_max(PE(betagridpoint,1+(1:sizehgrid))-POast(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-CO(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-LO(betagridpoint,1+(1:sizehgrid)))]*ones(1,2))
subplot(1,2,2)
hold on
plot(hgrid,ones(3,1)*PE(betagridpoint,1+(1:sizehgrid-1))-[POast(betagridpoint,1+(1:sizehgrid-1));CO(betagridpoint,1+(1:sizehgrid-1));LO(betagridpoint,1+(1:sizehgrid-1))])
plot([0 hmaxPE],[PUR_max(PE(betagridpoint,1+(1:sizehgrid))-POast(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-CO(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-LO(betagridpoint,1+(1:sizehgrid)))]*ones(1,2))
uu=max(TOPL(betagridpoint,1+(1:sizehgrid-1)),TOPR(betagridpoint,1+(1:sizehgrid-1)));
uugrid=hgrid(uu>0);
uu=uu(uu>0);
plot(uugrid,uu,'-k')
plot(hast(betagridpoint,2),0,'d')

end