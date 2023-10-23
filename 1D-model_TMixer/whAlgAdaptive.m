function [w,L]=whAlgAdaptive(mom)
n=length(mom)/2;

sigma=zeros(1,n*(n+1));
sigma(1:2*n)=mom;

a=zeros(1,n);

atmp=sigma(2)/sigma(1);
a(1)=atmp;

for i=1:2*n-2
    sigma(i+2*n)=sigma(i+2)-atmp*sigma(i+1);
end

b=zeros(1,n-1);
btmp = sigma(2*n+1)/sigma(1);
b(1)=btmp;

zeta=zeros(1,2*n-1);
zeta(1)=atmp;

zetatmp = btmp / atmp;
zeta(2)=zetatmp;

if zetatmp>0
    atmp=sigma(2*n+2)/sigma(2*n+1)-atmp;
    a(2)= atmp;

    zeta(3)=atmp-zetatmp;
end

for j=3:n
    if zetatmp<=0
        break
    end
    nCol=2*n-2*j+2;

    fb=2*n*j+j*(1-j)-2;
    fb_1=fb-nCol-2;
    fb_2=fb_1-nCol-4;

    for i=1:nCol
        sigma(i+fb)=sigma(i+2+fb_1)-atmp*sigma(i+1+fb_1)-btmp*sigma(i+2+fb_2);
    end

    atmp=sigma(fb+2)/sigma(fb+1)-sigma(fb_1+2)/sigma(fb_1+1);
    btmp=sigma(fb+1)/sigma(fb_1+1);

    a(j)=atmp;
    b(j-1)=btmp;

    zetatmp=btmp/zeta(2*j-3);
    zeta(2*j-2)=zetatmp;
    zeta(2*j-1)=atmp-zeta(2*j-2);
end

[bool_,rN]=nodeReduction(zeta,n);
aR=zeros(1,rN);
bR=zeros(1,rN-1);

for i=1:rN-1
    aR(i)=a(i);
    bR(i)=-power(b(i),0.5);
end

aR(rN)=a(rN);
aR=diag(aR);
aR=aR+diag(bR,-1);
aR=aR+diag(bR,+1);
[V,D]=eig(aR);

w=zeros(1,n);
L=zeros(1,n);

for i=1:rN
    w(i)=power(V(1,i),2)*mom(1);
    if w(i) < 0
        disp('Negative weight')
    elseif w(i) > 1e-15
        L(i)=D(i,i);
        if L(i)<0
            disp('Negative node')
        end
    else
        L(i)=1e-20;
    end
end

for i=rN+1:n
    w(i)=0;
    L(i)=1e-20;
end
end
