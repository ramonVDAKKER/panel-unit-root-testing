function PURtable(REPS)

TESTS=0;
Q=0;
h=0;
load PURDEF

alpha=1/Q;
DIMh=length(h);

DISTS=size(TESTS,1);
DIMNTchoice=size(TESTS,2);
Ntests=size(TESTS,4)/DIMh;

NTchoice=[10 50;25 50;      100 50;10 100;25 100;       100 100;10 250;25 250;       100 250]; %choice of (N,T) pairs

Power=zeros(DISTS,DIMNTchoice,Ntests*DIMh);
PowerMPPast=zeros(DISTS,DIMNTchoice,2*DIMh);
for test=1:Ntests*DIMh
    for ii=1:DIMNTchoice
        if test<=2*DIMh
            CV=(chi2inv(alpha,NTchoice(ii,1))-NTchoice(ii,1))/sqrt(2*NTchoice(ii,1));
        else
            CV=norminv(alpha);
        end
        for dist=1:DISTS
            Power(dist,ii,test)=sum(TESTS(dist,ii,:,test)<=CV,3);
            if test<=2*DIMh
                PowerMPPast(dist,ii,test)=sum(TESTS(dist,ii,:,(Ntests-2)*DIMh+test)<=CV,3);
            end
        end
    end
end
fprintf('dist    rep    N    T    infeasible PUR/BM/MPP/MPP*    feasible PUR/BM/MPP/MPP*      PE\n')
for dist=1:DISTS
    for ii=1:DIMNTchoice
        fprintf('%4.0f %6.0f %4.0f %4.0f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',[dist*ones(1,1) REPS*ones(1,1) NTchoice(ii,:) reshape(Power(dist,ii,1:2*DIMh),1,2*DIMh)/REPS ones(1,1)*(normcdf(norminv(alpha)-h'/sqrt(2)))]')
        for jj=2:Ntests/2
            fprintf('                         %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',reshape(Power(dist,ii,2*(jj-1)*DIMh+1:2*jj*DIMh),1,2*DIMh)'/REPS)
        end
        fprintf('                         %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',reshape(PowerMPPast(dist,ii,1:2*DIMh),1,2*DIMh)'/REPS)
    end
end

UU=zeros(Ntests,REPS);
Power=zeros(DISTS,DIMNTchoice,Ntests*DIMh);
for dist=1:DISTS
    for ii=1:DIMNTchoice
        for jj=1:Ntests
            UU(jj,:)=sort(reshape(TESTS(1,ii,:,(jj-1)*DIMh+1),1,REPS));
        end
        if REPS<Q
            CV=ones(REPS,1)*kron(UU(:,floor(REPS/Q)+1)',ones(1,DIMh));
        else
            CV=(1-(REPS+1)/Q+floor(REPS/Q))*ones(REPS,1)*kron(UU(:,floor(REPS/Q))',ones(1,DIMh))+((REPS+1)/Q-floor(REPS/Q))*ones(REPS,1)*kron(UU(:,floor(REPS/Q)+1)',ones(1,DIMh));
        end
        Power(dist,ii,:)=sum(reshape(TESTS(dist,ii,:,:),REPS,Ntests*DIMh)<=CV,1);
    end
end

fprintf('dist    rep    N    T    infeasible PUR/BM/MPP         feasible PUR/BM/MPP           PE [all powers size corrected]\n')
for dist=1:DISTS
    for ii=1:DIMNTchoice
        fprintf('%4.0f %6.0f %4.0f %4.0f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',[dist*ones(1,1) REPS*ones(1,1) NTchoice(ii,:) reshape(Power(dist,ii,1:2*DIMh),1,2*DIMh)/(floor((REPS+1)/Q)*Q) ones(1,1)*(normcdf(norminv(alpha)-h'/sqrt(2)))]')
        for jj=2:Ntests/2
            fprintf('                         %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',reshape(Power(dist,ii,2*(jj-1)*DIMh+1:2*jj*DIMh),1,2*DIMh)'/(floor((REPS+1)/Q)*Q))
        end
    end
end

end