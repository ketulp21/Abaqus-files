clear all
clc
format long
%Calculation of lam3
lam2=1:0.1:5;
alpha1=1.632;
mu1=0.053;
D1=0.001;
lam1=1;
    A1=4*mu1/3/alpha1;
	A2=A1/2;
	p1=(2*alpha1/3)-1;
	p2=(-alpha1/3)-1;
	B=2*lam1*lam2/D1;
    C=2/D1;
    x(1,:)=1./(lam1.*lam2)
for n=1:5
      %  x(n,:)=1./(lam1.*lam2)
	fx(n,:)=A1.*lam1.^p2.*lam2.^p2.*x(n,:).^p1-A2.*lam1.^p1*lam2.^p2.*x(n,:).^p2-A2.*lam1.^p2.*lam2.^p1.*x(n,:).^p2+B.*x(n,:)-C
	f_x(n,:)=A1.*(p1-1).*lam1.^p2.*lam2.^p2.*x(n,:).^(p1-1)-A2.*(p2-1).*lam1.^p1.*lam2.^p2.*x(n,:).^(p2-1)-A2.*(p2-1).*lam1.^p2.*lam2.^p1.*x(n,:).^(p2-1)+B
        x(n+1,:)=x(n,:)-(fx(n,:)./f_x(n,:))        
end
lam3=x(6,:)
%Calculation of Stress
    detF=lam1.*lam2.*lam3;
    wbar1=detF.^(-1/3).*lam1;
    wbar2=detF.^(-1/3).*lam2;
    wbar3=detF.^(-1/3).*lam3;

    a1=(2/3).*wbar1.^alpha1 - (1/3).*wbar2.^alpha1-(1/3).*wbar3.^alpha1;
    a2=-(1/3).*wbar1.^alpha1 + (2/3).*wbar2.^alpha1-(1/3).*wbar3.^alpha1;
   % a3=-(1/3)*wbar1^alpha1 - (1/3)*wbar2^alpha1+(2/3)*wbar3^alpha1;
    t1=2.*mu1./(detF.*alpha1);
    t2=(2.*(detF-1))./D1;

	  sigma1=t1.*a1 + t2;
	  sigma2=t1.*a2 + t2;

B=xlsread("Book1.xlsx",1);
S2=B(:,3);
lam=B(:,1);
S1=B(:,2);
figure()
plot(lam2,sigma2,'-o')
xlabel('Stretch')
ylabel('Stress')
xlim([1 5.5])
ylim([0 1.1])
hold on
%S2=movmean(S2,500);
plot(lam,S2,'-r')
legend('Calculated','Abaqus')
hold off

figure()
plot(lam2,sigma1,'-o')
xlabel('stretch')
ylabel('Stress1')
xlim([1 5.5])
ylim([0 0.1])
hold on
%S1=movmean(S1,500);
plot(lam,S1,'-r')
legend('Calculated','Abaqus')
hold off
