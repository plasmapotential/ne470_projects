clear all;
clc;
format long g;
prompt1 = 'How many nodes are desired: ';
prompt2 = 'Width of the slab in cm: ';
N = input(prompt1);
W = input(prompt2);
Sigtr = 3.62.*10.^(-2);
Sigabs = .1532;
vsig = .1570;
S = 10.^8;
D = 1./(3.*Sigtr);
dx = W./(N-1);
x = 0:dx:W;
A = zeros(N-1);
A(1,1) = (Sigabs./2 + D./dx.^2);
A(1,2) = (-D./dx.^2);
A(N-1,N-2) = (-D./dx.^2);
A(N-1,N-1) = (Sigabs + 2.*D./dx.^2); 
for i = 2:N-2
    A(i,i) = (Sigabs + 2.*D./dx.^2);
    A(i,i+1) = (-D./dx.^2);
    A(i,i-1) = (-D./dx.^2);
end
for i = 1:N-1
    if i==1
        B(i,i) = vsig./2;
    else
        B(i,i) = vsig;
    end
end
% Initial Guess for Flux, k-eff, and Source
for i = 1:N-1;
flux(i,1) = 1;
end
k = 1.00;
S=(1./k).*B*flux;
plot(flux);
flxdiff=1.0;
 
% Iteration
while(flxdiff>0.001)
oldflux=flux;
oldk=k;
oldS=S;
 
flux=inv(A)*oldS;
 
S=1/oldk*B*flux;
k=oldk*sum(S)/sum(oldS);
flxdiff=0.0;
for i=1:N-1
flxdiff = flxdiff + (flux(i)-oldflux(i))^2;
end
flxdiff = sqrt(flxdiff)
end
critlength = pi./2.*sqrt(D./(vsig-Sigabs))
% Normalize Flux (Integral of flux set to unity)
 
total = sum(flux)
flux = flux/total
 
% Add last node (zero flux condition)
  
flux(N)=0
plot(x,flux,'-*')
 
kvalue=num2str(k)
 
xlabel('Width of Slab (cm)')
ylabel('Flux')
title('Flux vs Width of the Slab')
text(10,0.025,'k-eff = ')
text(17,0.025,kvalue)
 

