function PUR(REPS)

tic
toc

stream=RandStream('mt19937ar', 'Seed', 20140117);

Q=20;
if REPS<Q-1 
    fprintf('WARNING: no size corrected powers possible') 
end
alpha=1/Q;

R=-2;            %use invariance of tests and take R=-2

Ntests=3*2;

DISTS=5;         %nunmber of distributions: U(0,2) chi_1^2/1 chi_2^2/2 chi_4^2/4 1

h=[0 -1 -2 -4]'; %points on power curve
DIMh=length(h);

NTchoice=[10 50;25 50;100 50;10 100;25 100;100 100;10 250;25 250;100 250]; %choice of (N,T) pairs
DIMNTchoice=size(NTchoice,1);
NMAX=max(NTchoice(:,1));
TMAX=max(NTchoice(:,2));

TESTS=zeros(DISTS,DIMNTchoice,REPS,Ntests*DIMh);

Power=zeros(DISTS,DIMNTchoice,(Ntests+2)*DIMh); %include MPPast 
for rep=1:REPS
    innovation=randn(stream,NMAX,TMAX); %*** NMAXxTMAX *** innovation fixed in each replication whatever choice of distribution or (N,T)
    HH=rand(stream,1,NMAX);             %***  1xNMAX   *** quantiles of distribution fixed in each replication whatever choice of (N,T)
    SIG2=.5+rand(stream,1,NMAX);        %***  1xNMAX   *** heterogeneous sigma's; in fact not relevant due to invariance, MPP choice
    B=.5+rand(stream,NMAX,1);           %***  NMAXx1   *** intercepts; in fact not relevant due to invariance
    
    for ii=1:DIMNTchoice
        N=NTchoice(ii,1); 
        T=NTchoice(ii,2);
        b=0;
        sig2=ones(1,N);
        if R>=-1.5                                     %simulated intercepts and sigmas, due to invariance irrelevant
            if R==0                                    %homogeneous sigma^2
                b=reshape(kron(B(1:N),ones(DIMh,T)),DIMh,N,T);
                sig2=ones(1,N);
            elseif R==-1                               %equally spaced sigma^2 on .5 to 1.5
                b=reshape(kron(B(1:N),ones(DIMh,T)),DIMh,N,T);
                sig2=.5+1/(N+1):1/(N+1):1.5-1/(N+1)/2; 
            elseif R==1                                %random sigma^2 from U(.5,1.5)
                b=reshape(kron(B(1:N),ones(DIMh,T)),DIMh,N,T);
                sig2=SIG2(1:N); 
            elseif R==999                              %extreme sigma^2: 1 to N
                b=reshape(kron(B(1:N),ones(DIMh,T)),DIMh,N,T);
                sig2=1:N;                              %***   1xN    ***
            end
        end

        for dist=1:DISTS
            if dist==1
                H=unifinv(HH(1:N),0,2);
            elseif dist==2
                H=chi2inv(HH(1:N),1)/1; 
            elseif dist==3
                H=chi2inv(HH(1:N),2)/2; 
            elseif dist==4
                H=chi2inv(HH(1:N),4)/4;                %***   1xN    ***
            else H=ones(1,N);
            end

            rho=ones(DIMh,N); rho=rho+h*H/(sqrt(N)*T); %***  DIMhxN  *** model parameter
            
            Observation=zeros(DIMh,N,T);               %*** DIMhxNxT *** generating the DGP
            Observation(:,:,1)=ones(DIMh,1)*(sqrt(sig2).*innovation(1:N,1)');
            for t=2:T
                Observation(:,:,t)=rho.*Observation(:,:,t-1)+ones(DIMh,1)*(sqrt(sig2).*innovation(1:N,t)');
            end
            Observation=Observation+b;
            
            Delta=Observation(:,:,2:T)-Observation(:,:,1:T-1);                         %*** DIMh*N*(T-1) *** difference of Observations
            sig2hat=sum(Delta.^2,3)/(T-1);                                             %***    DIMh*N    *** sigmasquare estimates
            Deltadivsig=Delta./reshape(kron(ones(DIMh,T-1),sqrt(sig2)),DIMh,N,T-1);    %*** DIMh*N*(T-1) *** infeasible
            Deltadivsighat=Delta./reshape(kron(ones(1,T-1),sqrt(sig2hat)),DIMh,N,T-1); %*** DIMh*N*(T-1) *** feasible
            cumsumDeltadivsig=cumsum(Deltadivsig,3);                                   %*** DIMh*N*(T-1) ***
            cumsumDeltadivsighat=cumsum(Deltadivsighat,3);                             %*** DIMh*N*(T-1) ***
            
            info=sum(sum(cumsumDeltadivsig(:,:,1:T-2).^2,3),2)'/(N*T^2);
            infohat=sum(sum(cumsumDeltadivsighat(:,:,1:T-2).^2,3),2)'/(N*T^2);

            PURtestinf=sqrt(2)*sum(sum(cumsumDeltadivsig(:,:,1:T-2).*Deltadivsig(:,:,2:T-1),3),2)'/(sqrt(N)*T);
            PURtest=sqrt(2)*sum(sum(cumsumDeltadivsighat(:,:,1:T-2).*Deltadivsighat(:,:,2:T-1),3),2)'/(sqrt(N)*T);
            PURtestinfE=sqrt(.5./info).*PURtestinf;
            PURtestE=sqrt(.5./infohat).*PURtest;

            MPPtestinf=2*sum(sum((Observation(:,:,1:T-1)./reshape(kron(ones(DIMh,T-1),sqrt(sig2)),DIMh,N,T-1)).*Deltadivsig,3),2)'/(sqrt(N)*(T-1))+sum(sum((Observation(:,:,1:T-1)./reshape(kron(ones(DIMh,T-1),sqrt(sig2)),DIMh,N,T-1)).^2,3),2)'/(N*(T-1)^2)-1/2;
            MPPtest=2*sum(sum((Observation(:,:,1:T-1)./reshape(kron(ones(1,T-1),sqrt(sig2hat)),DIMh,N,T-1)).*Deltadivsighat,3),2)'/(sqrt(N)*(T-1))+sum(sum((Observation(:,:,1:T-1)./reshape(kron(ones(1,T-1),sqrt(sig2hat)),DIMh,N,T-1)).^2,3),2)'/(N*(T-1)^2)-1/2;
            MPPtestinf=MPPtestinf+sum((Observation(:,:,1).^2-1/(1+1/(N*(T-1)))*(Observation(:,:,1)+1/(sqrt(N)*(T-1))*(Observation(:,:,T)-Observation(:,:,1))+1/(N*(T-1)^2)*sum(Observation(:,:,1:T-1),3)).^2)./reshape(kron(ones(DIMh,1),sig2),DIMh,N),2)';
            MPPtest=MPPtest+sum((Observation(:,:,1).^2-1/(1+1/(N*(T-1)))*(Observation(:,:,1)+1/(sqrt(N)*(T-1))*(Observation(:,:,T)-Observation(:,:,1))+1/(N*(T-1)^2)*sum(Observation(:,:,1:T-1),3)).^2)./sig2hat,2)';
            MPPtestinf=MPPtestinf/sqrt(2);
            MPPtest=MPPtest/sqrt(2);
            cv=[ones(1,2*DIMh)*(chi2inv(alpha,N)-N)/sqrt(2*N) ones(1,(Ntests-2)*DIMh)*norminv(alpha) ones(1,2*DIMh)*(chi2inv(alpha,N)-N)/sqrt(2*N)];
            Power(dist,ii,:)=Power(dist,ii,:)+reshape(([PURtestinf PURtest PURtestinfE PURtestE MPPtestinf MPPtest MPPtestinf MPPtest]<cv),1,1,(Ntests+2)*DIMh); %include MPPast
       
            TESTS(dist,ii,rep,:)=[PURtestinf PURtest PURtestinfE PURtestE MPPtestinf MPPtest];
        end
    end
    if floor(rep/10000)*10000==rep
        fprintf(' %5.0f %20.4f\n',rep,toc)
    end
end

fprintf('dist   rep    N    T    infeasible PUR/BM/MPP/MPP*    feasible PUR/BM/MPP/MPP*      PE\n')
for dist=1:DISTS
    for ii=1:DIMNTchoice
        fprintf('%4.0f %5.0f %4.0f %4.0f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',[dist*ones(1,1) rep*ones(1,1) NTchoice(ii,:) reshape(Power(dist,ii,1:2*DIMh),1,2*DIMh)/REPS ones(1,1)*(normcdf(norminv(alpha)-h'/sqrt(2)))]')
        for jj=2:(Ntests+2)/2 %include MPPast
            fprintf('                        %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',reshape(Power(dist,ii,2*(jj-1)*DIMh+1:2*jj*DIMh),1,2*DIMh)'/REPS)
        end
    end
end

Power=zeros(DISTS,DIMNTchoice,Ntests*DIMh); %exclude irrelevant MPPast 
UU=zeros(Ntests,REPS);
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

fprintf('dist   rep    N    T    infeasible PUR/BM/MPP         feasible PUR/BM/MPP           PE [all powers size corrected]\n')
for dist=1:DISTS
    for ii=1:DIMNTchoice
        fprintf('%4.0f %5.0f %4.0f %4.0f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',[dist*ones(1,1) rep*ones(1,1) NTchoice(ii,:) reshape(Power(dist,ii,1:2*DIMh),1,2*DIMh)/(floor((REPS+1)/Q)*Q) ones(1,1)*(normcdf(norminv(alpha)-h'/sqrt(2)))]')
        for jj=2:Ntests/2
            fprintf('                        %5.3f  %5.3f  %5.3f  %5.3f    %5.3f  %5.3f  %5.3f  %5.3f\n',reshape(Power(dist,ii,2*(jj-1)*DIMh+1:2*jj*DIMh),1,2*DIMh)'/(floor((REPS+1)/Q)*Q))
        end
    end
end

save('PURDEF.mat','TESTS','h','Q')
toc
end