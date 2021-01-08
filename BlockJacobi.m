clear;
testcase=input('Please input testcase ','s');
matfilename=strcat('../GraphSparsification/',testcase,'.mat');

load(matfilename,'LG','b');
disp('Graph read finished');
[n,n]=size(LG);

LG=LG(1:n-1,1:n-1);
b=b(1:n-1,:);n=n-1;
nparts=2;
x1=LG\b;

[i,j,k]=find(tril(LG,-1));
[p,q1,q2]=mexmetis(i,j,k,[n],[nparts]);
LG=LG(p,p);b=b(p,1);

p2=1:n;
factors=cell(nparts,2);
for i=1:nparts
    subdomain=q1(i):q1(i+1)-1;
    B=LG(subdomain,subdomain);
    P=amd(B);
    R=chol(B(P,P));
    factors(i,:)={R,subdomain};
    p2(1,subdomain)=P+(q1(i)-1);
end
LG=LG(p2,p2);b=b(p2,1);
tic;[x,flag,relres,iter]=pcg(LG,b,1e-3,1000,@(x) BlockJacobiPre(x,factors,nparts));toc;
x2(p2,1)=x;x(p,1)=x2;
iter

function y=BlockJacobiPre(x,factors,nparts)
    y=zeros(size(x,1),1);
    for i=1:nparts
        R=factors{i,1};subdomain=factors{i,2};
        y(subdomain,1)=R\(R'\x(subdomain,1));
    end
end