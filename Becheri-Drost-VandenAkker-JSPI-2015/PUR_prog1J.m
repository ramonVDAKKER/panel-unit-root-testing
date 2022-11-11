function PUR_prog1J(NN,stepsize,hmaxPE,stepsizebeta,betamax,alpha)
%simulates J=\int W^2_t dt

if nargin<1
    NN=6;
end
if nargin<2
    stepsize=1/100;
end
if nargin<3
    hmaxPE=10;
end
if nargin<4
    stepsizebeta=1/10;
end
if nargin<5
    betamax=5;
end
if nargin<6
    alpha=0.05;
end
%save relevant inputs for PUR_prog2PE
save INPUTS.mat stepsize hmaxPE stepsizebeta betamax alpha

sizehgrid=round(hmaxPE/stepsize+2); % 0:stepsize:hmaxPE + extreme point
sizebetagrid=round(betamax/stepsizebeta+1); % 0:stepsizebeta:betamax
%initialize variables for PUR_prog2PE
PE=zeros(sizebetagrid,1+sizehgrid);
CO=zeros(sizebetagrid,1+sizehgrid);
LO=zeros(sizebetagrid,1+sizehgrid);
CV=zeros(sizebetagrid,1+sizehgrid);
save PECOLOCVbeta.mat PE CO LO CV
%initialize variables for PUR_prog3comparison
PO=zeros(sizebetagrid,1+sizehgrid);
TOPL=zeros(sizebetagrid,sizehgrid);
TOPR=zeros(sizebetagrid,sizehgrid);
h0=zeros(sizebetagrid,2);
save PObeta.mat PO TOPL TOPR h0
save PObeta2.mat PO TOPL TOPR h0
%initialize variables for PUR_prog4POast
POast=zeros(sizebetagrid,1+sizehgrid);
hast=zeros(sizebetagrid,2);
Fast=zeros(sizebetagrid,2);
save POastbeta.mat POast hast Fast

%simulating J
if NN<=0
    stop%not allowed
else
    stream=RandStream('mt19937ar', 'Seed', 20140124);
    REPS=10^NN;
    dt=REPS;
    J=zeros(REPS,1);
%     fprintf('\nPROCESS              cumlative time')
    for rep=1:REPS
        J(rep)=sum(cumsum(randn(stream,1,dt)).^2/dt^2);
%         if floor(rep/REPS*10)*REPS/10==rep
%             fprintf('\nsimulating J%6.0f%% %15.0f',rep/REPS*100,toc)
%         end
    end
    save J.mat J
end

end