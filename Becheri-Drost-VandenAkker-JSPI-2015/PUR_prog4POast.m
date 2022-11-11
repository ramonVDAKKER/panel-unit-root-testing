function PUR_prog4POast(beta,delta)
%improves preliminary power curve PO to POast using improved value of h0 from top function information

%load inputs PUR_prog1J
load INPUTS.mat stepsize hmaxPE stepsizebeta betamax alpha
if nargin<1
    beta=4;
end
if nargin<2
    delta=1;
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

load PObeta.mat TOPL TOPR h0

range=1+(1:delta:sizehgrid-1);

%linearly interpolate left/right top functions to determine crossing point
y=-abs(TOPR(betagridpoint,range)-TOPL(betagridpoint,range)); 
yL=length(y);
ymin=min(y);
ycorr=ymin-1;
y=y+ycorr*(y==0).*or([0 y(1:yL-1)]==0,[y(2:yL) 0]==0);

[ymax,ypos]=max(y); 
pos=1+(ypos-1)*delta+0*ymax;
if TOPR(betagridpoint,1+pos)-TOPL(betagridpoint,1+pos)<0; pos=pos-delta; end;
htemp=(pos-1)*stepsize;

FT=[TOPR(betagridpoint,1+(pos:delta:pos+delta)) TOPL(betagridpoint,1+(pos:delta:pos+delta))];  
if abs(h0(betagridpoint,2)-(FT(1)>FT(4))*delta*stepsize-htemp)>delta*stepsize/2
    fprintf('WARNING %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',[h0(betagridpoint,2) htemp FT])
end

hopt=htemp+delta*stepsize*(FT(1)-FT(3))/(FT(4)-FT(3)-FT(2)+FT(1));
Fopt=(FT(1)*FT(4)-FT(2)*FT(3))/(FT(4)-FT(3)-FT(2)+FT(1)); 

%load and calculate J
load J.mat J
J0=1/2; 
J=J0+beta*J;

%calculating critical value
    c=0;
    h=hopt;
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
    cv=c;

%evaluate point optimal test
    powpointopt=zeros(sizehgrid,1);
    h00=hopt;
    for hh=1:sizehgrid
        if hh<sizehgrid
            h=(hh-1)*stepsize;
        else
            h=extr*hmaxPE;
        end
        powpointopt(hh)=mean(normcdf(-cv./sqrt((h00^2+(h00==0))*J)+(h-.5*h00)*sqrt(J)));
    end

load POastbeta.mat hast Fast POast
hast(betagridpoint,:)=[beta hopt];
Fast(betagridpoint,:)=[beta Fopt];
POast(betagridpoint,:)=[beta powpointopt'];
save POastbeta.mat hast Fast POast

end
