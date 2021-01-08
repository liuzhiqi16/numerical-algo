clear;
% testcase: m14b
testcase=input('Please input testcase ','s');
matfilename=strcat('./',testcase,'.mat');

load(matfilename,'LG','b');
disp('Graph read finished');
[n,n]=size(LG);

nparts=2;

[i,j,k]=find(tril(LG,-1));

% p 是重排序向量，q1(i) 代表第i个子区域起始位置，q2(i) 代表第i个子区域内部节点的起始位置
[p,q1,q2]=mexmetis(i,j,k,[n],[nparts]);

LG=LG(p,p);b=b(p,1);

% 根据W得到E
W=-LG(q2(1):q1(2)-1,q2(2):q1(3)-1);
[m,mm]=size(W);
X1=zeros(m,1);
for i=1:m
    X1(i)=sqrt(norm(W(i,:)));
    W(i,:)=W(i,:)/X1(i);
end
X1=spdiags(X1,[0],m,m);
X2=W;
E=[zeros(q2(1)-q1(1),m);X1;zeros(q2(2)-q1(2),m);X2'];

% 单级方法，直接对对角块做分解
B=LG+E*E';
m1=q1(2)-1;
B11=B(1:m1,1:m1);B22=B(m1+1:end,m1+1:end);
p2=1:n;
P=1:size(B11,1);P(1,1:q2(1)-q1(1))=amd(B11(1:q2(1)-1,1:q2(1)-1));R1=chol(B11(P,P));p2(1,1:m1)=P;
P=1:size(B22,1);P(1,1:q2(2)-q1(2))=amd(B22(1:q2(2)-q1(2),1:q2(2)-q1(2)));R2=chol(B22(P,P));p2(1,m1+1:end)=P+m1;
LG=LG(p2,p2);b=b(p2,1);

% 低秩修正，rank=2
k=2;
[U,S,V]=svds(@(b,tflag)Afun(b,tflag,E,R1,R2,m1),size(E),k);
U=U*S;
V=V(:,1:k);
H=inv(eye(k)-U'*E*V);

% PCG求解方程
tic;[x,flag,relres,iter,RESVEC]=pcg(LG,b,1e-3,1000,@(x)MLRsolve(x,R1,R2,m1,U,H));toc;
disp(['MLR: ',num2str(iter)]);

function y=Afun(b,tflag,E,R1,R2,m1)
    if strcmp(tflag,'notransp')
        x=E*b;
        x1=x(1:m1,:);y1=R1\(R1'\x1);
        x2=x(m1+1:end,:);y2=R2\(R2'\x2);
        y=[y1;y2];
    else
        x1=b(1:m1,:);y1=R1\(R1'\x1);
        x2=b(m1+1:end,:);y2=R2\(R2'\x2);
        y=E'*[y1;y2];
    end
end

function y=MLRsolve(x,R1,R2,m1,U,H)
    y=zeros(size(x,1),1);
    y(1:m1,:)=R1\(R1'\x(1:m1,:));
    y(m1+1:end,:)=R2\(R2'\x(m1+1:end,:));
    y0=U*H*(U'*x);
    y=y+y0;
end




