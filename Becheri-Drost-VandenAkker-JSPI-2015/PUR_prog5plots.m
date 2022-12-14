function PUR_prog5plots(beta)
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
sizebetagrid=round(betamax/stepsizebeta+1);

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
if beta==4
% Create textbox
annotation('textbox',...
    [0.246815476190483 0.842612419700214 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{PE}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.277529761904769 0.888993576017133 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{PO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
% Create textbox
annotation('textbox',...
    [0.284077380952388 0.83438972162741 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{CO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0.498039215686275 0]);
% Create textbox
annotation('textbox',...
    [0.313005952380959 0.7769164882227 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{LO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[1 0 0]);
% Create textbox
annotation('textbox',...
    [0.696242559523819 0.170376486902083 0.0570907738095167 0.00127233536772472],...
    'String',{'\pi_{PE} - \pi_{PO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
% Create textbox
annotation('textbox',...
    [0.762313988095247 0.234616315595874 0.0570907738095167 0.00127233536772472],...
    'String',{'\pi_{PE} - \pi_{CO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0.498039215686275 0]);
% Create textbox
annotation('textbox',...
    [0.816718750000007 0.467678414096943 0.0570907738095167 0.00127233536772472],...
    'String',{'\pi_{PE} - \pi_{LO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[1 0 0]);
% Create textbox
annotation('textbox',...
    [0.462187500000026 0.135974304068522 0.0512499999999921 0.00400825401617627],...
    'String',{'d(\pi_{PE},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
% Create textbox
annotation('textbox',...
    [0.903125000000009 0.208779443254818 0.0512499999999921 0.00287325604430539],...
    'String',{'d(\pi_{PE},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
% Create textbox
annotation('textbox',...
    [0.903593750000009 0.245586799234035 0.0512499999999921 0.000665877425493436],...
    'String',{'d(\pi_{PE},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0.498039215686275 0]);
% Create textbox
annotation('textbox',...
    [0.462656250000026 0.161670235546039 0.0512499999999921 0.00242650152346432],...
    'String',{'d(\pi_{PE},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0.498039215686275 0]);
% Create textbox
annotation('textbox',...
    [0.463906250000027 0.267665952890792 0.0512499999999921 0.00302145100894517],...
    'String',{'d(\pi_{PE},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[1 0 0]);
% Create textbox
annotation('textbox',...
    [0.904843750000009 0.850107066381156 0.0512499999999921 0.00225434632912638],...
    'String',{'d(\pi_{PE},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[1 0 0]);
% Create textbox
annotation('textbox',...
    [0.700156250000011 0.338329764453961 0.0512499999999921 0.000846446999794348],...
    'String',{'d(\pi_{PE},\pi_h)'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textarrow
annotation('textarrow',[0.43671875 0.46484375],...
    [0.0794271335452657 0.0793258119593626],'TextEdgeColor','none',...
    'FontSize',12,...
    'String',{'- h'});
% Create textarrow
annotation('textarrow',[0.876361607142871 0.904486607142871],...
    [0.0808404097765292 0.0807390881906261],'TextEdgeColor','none',...
    'FontSize',12,...
    'String',{'- h'});
% Create textbox
annotation('textbox',...
    [0.644642857142869 0.110278372591007 0.0273809523809409 0.00321199143468951],...
    'String',{'- h_0^\ast'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
end
subplot(1,2,1)
hold on
hgrid=stepsize*(0:sizehgrid-2);
plot(hgrid,[POast(betagridpoint,1+(1:sizehgrid-1));CO(betagridpoint,1+(1:sizehgrid-1));LO(betagridpoint,1+(1:sizehgrid-1))])
plot(hgrid,PE(betagridpoint,1+(1:sizehgrid-1)),'-k')
plot([0 hmaxPE],[PUR_max(PE(betagridpoint,1+(1:sizehgrid))-POast(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-CO(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-LO(betagridpoint,1+(1:sizehgrid)))]*ones(1,2))
axis([0 betamax 0 1])
subplot(1,2,2)
hold on
plot(hgrid,ones(3,1)*PE(betagridpoint,1+(1:sizehgrid-1))-[POast(betagridpoint,1+(1:sizehgrid-1));CO(betagridpoint,1+(1:sizehgrid-1));LO(betagridpoint,1+(1:sizehgrid-1))])
plot([0 hmaxPE],[PUR_max(PE(betagridpoint,1+(1:sizehgrid))-POast(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-CO(betagridpoint,1+(1:sizehgrid)));PUR_max(PE(betagridpoint,1+(1:sizehgrid))-LO(betagridpoint,1+(1:sizehgrid)))]*ones(1,2))
uu=max(TOPL(betagridpoint,1+(1:sizehgrid-1)),TOPR(betagridpoint,1+(1:sizehgrid-1)));
uugrid=hgrid(uu>0);
uu=uu(uu>0);
plot(uugrid,uu,'-k')
plot(hast(betagridpoint,2),0,'d')
%axis([0 hmaxPE 0 .2])
axis([0 betamax 0 .2])

figure(2)
if beta==4
% Create textbox
annotation('textbox',...
    [0.246815476190483 0.842612419700214 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{PE}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.277529761904769 0.888993576017133 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{PO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.284077380952388 0.83438972162741 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{CO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.313005952380959 0.7769164882227 0.000803571428565081 0.000368091388495009],...
    'String',{'\pi_{LO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.696242559523819 0.170376486902083 0.0570907738095167 0.00127233536772472],...
    'String',{'\pi_{PE} - \pi_{PO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.762313988095247 0.234616315595874 0.0570907738095167 0.00127233536772472],...
    'String',{'\pi_{PE} - \pi_{CO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.816718750000007 0.467678414096943 0.0570907738095167 0.00127233536772472],...
    'String',{'\pi_{PE} - \pi_{LO}'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.462187500000026 0.135974304068522 0.0512499999999921 0.00400825401617627],...
    'String',{'d(\pi_{PE},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.903125000000009 0.208779443254818 0.0512499999999921 0.00287325604430539],...
    'String',{'d(\pi_{PE},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.903593750000009 0.245586799234035 0.0512499999999921 0.000665877425493436],...
    'String',{'d(\pi_{PE},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.462656250000026 0.161670235546039 0.0512499999999921 0.00242650152346432],...
    'String',{'d(\pi_{PE},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.463906250000027 0.267665952890792 0.0512499999999921 0.00302145100894517],...
    'String',{'d(\pi_{PE},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.904843750000009 0.850107066381156 0.0512499999999921 0.00225434632912638],...
    'String',{'d(\pi_{PE},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textbox
annotation('textbox',...
    [0.700156250000011 0.338329764453961 0.0512499999999921 0.000846446999794348],...
    'String',{'d(\pi_{PE},\pi_h)'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'BackgroundColor','none',...
    'EdgeColor','none');
% Create textarrow
annotation('textarrow',[0.43671875 0.46484375],...
    [0.0794271335452657 0.0793258119593626],'TextEdgeColor','none',...
    'FontSize',12,...
    'String',{'- h'});
% Create textarrow
annotation('textarrow',[0.876361607142871 0.904486607142871],...
    [0.0808404097765292 0.0807390881906261],'TextEdgeColor','none',...
    'FontSize',12,...
    'String',{'- h'});
% Create textbox
annotation('textbox',...
    [0.644642857142869 0.110278372591007 0.0273809523809409 0.00321199143468951],...
    'String',{'- h_0^\ast'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none',...
    'EdgeColor','none');
end
subplot(1,2,1)
hold on
plot(hgrid,POast(betagridpoint,1+(1:sizehgrid-1)),'--k')
plot(hgrid,CO(betagridpoint,1+(1:sizehgrid-1)),':k','LineWidth',2)
plot(hgrid,LO(betagridpoint,1+(1:sizehgrid-1)),'-.k')
plot(hgrid,PE(betagridpoint,1+(1:sizehgrid-1)),'-k')
plot([0 hmaxPE],PUR_max(PE(betagridpoint,1+(1:sizehgrid))-POast(betagridpoint,1+(1:sizehgrid)))*ones(1,2),'--k')
plot([0 hmaxPE],PUR_max(PE(betagridpoint,1+(1:sizehgrid))-CO(betagridpoint,1+(1:sizehgrid)))*ones(1,2),':k','LineWidth',2)
plot([0 hmaxPE],PUR_max(PE(betagridpoint,1+(1:sizehgrid))-LO(betagridpoint,1+(1:sizehgrid)))*ones(1,2),'-.k')
axis([0 betamax 0 1])
subplot(1,2,2)
hold on
plot(hgrid,PE(betagridpoint,1+(1:sizehgrid-1))-POast(betagridpoint,1+(1:sizehgrid-1)),'--k')
plot(hgrid,PE(betagridpoint,1+(1:sizehgrid-1))-CO(betagridpoint,1+(1:sizehgrid-1)),':k','LineWidth',2)
plot(hgrid,PE(betagridpoint,1+(1:sizehgrid-1))-LO(betagridpoint,1+(1:sizehgrid-1)),'-.k')
plot([0 hmaxPE],ones(1,2)*PUR_max(PE(betagridpoint,1+(1:sizehgrid))-POast(betagridpoint,1+(1:sizehgrid))),'--k')
plot([0 hmaxPE],ones(1,2)*PUR_max(PE(betagridpoint,1+(1:sizehgrid))-CO(betagridpoint,1+(1:sizehgrid))),':k','LineWidth',2)
plot([0 hmaxPE],ones(1,2)*PUR_max(PE(betagridpoint,1+(1:sizehgrid))-LO(betagridpoint,1+(1:sizehgrid))),'-.k')
plot(uugrid,uu,'-k')
plot(hast(betagridpoint,2),0,'dk')
axis([0 betamax 0 .2])

figure(3)
if beta==4;
% Create textarrow
annotation('textarrow',[0.435014880952401 0.464702380952401],...
    [0.079244792399183 0.080244792399183],'TextEdgeColor','none','FontSize',12,...
    'String',{'\beta'});
% Create textarrow
annotation('textarrow',[0.875543154761908 0.905230654761908],...
    [0.081641636863192 0.082641636863192],'TextEdgeColor','none','FontSize',12,...
    'String',{'\beta'});
% Create textbox
annotation('textbox',...
    [0.413028273809527 0.190105333514872 0.0709374999999989 0.00204021301032475],...
    'String',{'d(\pi_{PE},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
% Create textbox
annotation('textbox',...
    [0.413898809523812 0.27631856437653 0.0659374999999991 0.00229415238647609],...
    'String',{'d(\pi_{PE},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0.498039215686275 0]);
% Create textbox
annotation('textbox',...
    [0.413950892857149 0.84078177024669 0.0613095238095239 0.00214132762312366],...
    'String',{'d(\pi_{PE},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[1 0 0]);
% Create textbox
annotation('textbox',...
    [0.85396577380953 0.930238516666463 0.0659374999999974 0.000253473206019498],...
    'String',{'d(\pi_{PO},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
% Create textbox
annotation('textbox',...
    [0.850074404761907 0.163005780346821 0.063020833333332 0.0046242774566474],...
    'String',{'- d(\pi_{LO},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0 1]);
% Create textbox
annotation('textbox',...
    [0.852633928571436 0.803143914544942 0.0759374999999942 0.00137417178952553],...
    'String',{'d(\pi_{CO},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0.498039215686275 0]);
% Create textbox
annotation('textbox',...
    [0.849099702380959 0.208950254360015 0.0860937499999949 0.00028727562615825],...
    'String',{'- d(\pi_{LO},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[0 0.498039215686275 0]);
% Create textbox
annotation('textbox',...
    [0.854508928571434 0.28749109202623 0.0701562499999959 0.000688284389459412],...
    'String',{'d(\pi_{PO},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[1 0 0]);
% Create textbox
annotation('textbox',...
    [0.608005952380959 0.154773427733287 0.067589285714279 0.000587647826337897],...
    'String',{'- d(\pi_{CO},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'Color',[1 0 0]);
end
betagrid=stepsizebeta*(0:sizebetagrid-1);
subplot(1,2,1)
hold on
diffPELO=zeros(1,sizebetagrid);
diffPECO=zeros(1,sizebetagrid);
diffPEPOast=zeros(1,sizebetagrid);
diffCOLO=zeros(1,sizebetagrid);
diffLOCO=zeros(1,sizebetagrid);
diffCOPOast=zeros(1,sizebetagrid);
diffPOastCO=zeros(1,sizebetagrid);
diffPOastLO=zeros(1,sizebetagrid);
diffLOPOast=zeros(1,sizebetagrid);
for i=1:sizebetagrid
    diffPELO(i)=PUR_max(PE(i,1+(1:sizehgrid))-LO(i,1+(1:sizehgrid)));
    diffPECO(i)=PUR_max(PE(i,1+(1:sizehgrid))-CO(i,1+(1:sizehgrid)));
    diffPEPOast(i)=PUR_max(PE(i,1+(1:sizehgrid))-POast(i,1+(1:sizehgrid)));
    diffCOLO(i)=PUR_max(CO(i,1+(1:sizehgrid))-LO(i,1+(1:sizehgrid)));
    diffLOCO(i)=PUR_max(LO(i,1+(1:sizehgrid))-CO(i,1+(1:sizehgrid)));
    diffCOPOast(i)=PUR_max(CO(i,1+(1:sizehgrid))-POast(i,1+(1:sizehgrid)));
    diffPOastCO(i)=PUR_max(POast(i,1+(1:sizehgrid))-CO(i,1+(1:sizehgrid)));
    diffLOPOast(i)=PUR_max(LO(i,1+(1:sizehgrid))-POast(i,1+(1:sizehgrid)));
    diffPOastLO(i)=PUR_max(POast(i,1+(1:sizehgrid))-LO(i,1+(1:sizehgrid)));
end
plot(betagrid,[diffPEPOast;diffPECO;diffPELO])
axis([0 betamax 0 .2])
subplot(1,2,2)
hold on
plot(betagrid,[diffPOastLO;diffCOLO;diffPOastCO])
plot(betagrid,-[diffLOPOast;diffLOCO;diffCOPOast])
%axis([0 betamax 0 .02])

figure(4)
if beta==4;
% Create textarrow
annotation('textarrow',[0.435014880952401 0.464702380952401],...
    [0.079244792399183 0.080244792399183],'TextEdgeColor','none','FontSize',12,...
    'String',{'\beta'});
% Create textarrow
annotation('textarrow',[0.875543154761908 0.905230654761908],...
    [0.081641636863192 0.082641636863192],'TextEdgeColor','none','FontSize',12,...
    'String',{'\beta'});
% Create textbox
annotation('textbox',...
    [0.413028273809527 0.190105333514872 0.0709374999999989 0.00204021301032475],...
    'String',{'d(\pi_{PE},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.413898809523812 0.27631856437653 0.0659374999999991 0.00229415238647609],...
    'String',{'d(\pi_{PE},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.413950892857149 0.84078177024669 0.0613095238095239 0.00214132762312366],...
    'String',{'d(\pi_{PE},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.85396577380953 0.930238516666463 0.0659374999999974 0.000253473206019498],...
    'String',{'d(\pi_{PO},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.850074404761907 0.163005780346821 0.063020833333332 0.0046242774566474],...
    'String',{'- d(\pi_{LO},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.852633928571436 0.803143914544942 0.0759374999999942 0.00137417178952553],...
    'String',{'d(\pi_{CO},\pi_{LO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.849099702380959 0.208950254360015 0.0860937499999949 0.00028727562615825],...
    'String',{'- d(\pi_{LO},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.854508928571434 0.28749109202623 0.0701562499999959 0.000688284389459412],...
    'String',{'d(\pi_{PO},\pi_{CO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
% Create textbox
annotation('textbox',...
    [0.608005952380959 0.154773427733287 0.067589285714279 0.000587647826337897],...
    'String',{'- d(\pi_{CO},\pi_{PO})'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor','none');
end
subplot(1,2,1)
hold on
plot(betagrid,diffPEPOast,'--k')
plot(betagrid,diffPECO,':k','LineWidth',2)
plot(betagrid,diffPELO,'-.k')
axis([0 betamax 0 .2])
subplot(1,2,2)
hold on
plot(betagrid,diffPOastLO,'--k')
plot(betagrid,diffCOLO,':k','LineWidth',2)
plot(betagrid,diffPOastCO,'-k')
plot(betagrid,-diffLOPOast,'--k')
plot(betagrid,-diffLOCO,':k','LineWidth',2)
plot(betagrid,-diffCOPOast,'-k')
%axis([0 betamax 0 .02])
end