% 2D Model problem
N=41;
disp(['N=',num2str(N)]);
nx=N-1;ny=N-1;n=nx*ny;hx=1/(nx+1);hy=1/(1+ny);
rhs=zeros(n,1);xacc=zeros(n,1);
ijk=zeros((ny-1)*nx+(nx-1)*ny+nx*ny,3);
edgeptr=0;
for i=1:nx
    for j=1:ny
       index=i+(j-1)*nx;
       x=i*hx;y=j*hy;
       rhs(index,1)=2*pi*pi*sin(pi*x)*sin(pi*y);
       edgeptr=edgeptr+1;
       ijk(edgeptr,:)=[index,index,2/(hx*hx)+2/(hy*hy)];
       left=index-1;right=index+1;up=index+nx;down=index-nx;
       if(i>1)
           edgeptr=edgeptr+1;
           ijk(edgeptr,:)=[index,left,-1/(hx*hx)];
       end
       if(i<nx)
           edgeptr=edgeptr+1;
           ijk(edgeptr,:)=[index,right,-1/(hx*hx)];
       end
       if(j>1)
           edgeptr=edgeptr+1;
           ijk(edgeptr,:)=[index,down,-1/(hy*hy)];
       end
       if(j<ny)
           edgeptr=edgeptr+1;
           ijk(edgeptr,:)=[index,up,-1/(hy*hy)];
       end
       xacc(index,1)=sin(pi*x)*sin(pi*y);
    end
end
if(edgeptr~=size(ijk,1))
    disp('ERROR!');
end
A=sparse(ijk(:,1),ijk(:,2),ijk(:,3),n,n);

nlevel=1;
m1=floor(0.5*ny)*nx;m2=n-m1;
ijk=zeros(nx,3);
for i=1:nx
    ijk(i,:)=[m1-i+1,nx-i+1,1];
end
E1=sparse(ijk(:,1),ijk(:,2),ijk(:,3),m1,nx);
for i=1:nx
    ijk(i,:)=[i,i,1];
end
E2=sparse(ijk(:,1),ijk(:,2),ijk(:,3),m2,nx);
A11=A(1:m1,1:m1);B1=A11+E1*E1';
A22=A(m1+1:n,m1+1:n);B2=A22+E2*E2';
L1=ichol(A11);L2=ichol(A22);
E=[E1;E2];
k=5;
[U,S,V]=svds(@(b,tflag)Afun(b,tflag,E,L1,L2,m1,n),size(E),k,'largest','SubspaceDimension',10*k,'Tolerance',1e-3);
U=U*S;
V=V(:,1:k);
H=inv(eye(k)-U'*E*V);
[x,flag,relres,iter]=pcg(A,rhs,1e-3,200,@(x)MLRsolve(x,L1,L2,m1,n,U,H));
disp(iter);

function y=Afun(b,tflag,E,L1,L2,m1,n)
    if strcmp(tflag,'notransp')
        x=E*b;
        x1=x(1:m1,:);y1=L1'\(L1\x1);
        x2=x(m1+1:n,:);y2=L2'\(L2\x2);
        y=[y1;y2];
    else
        x1=b(1:m1,:);y1=L1'\(L1\x1);
        x2=b(m1+1:n,:);y2=L2'\(L2\x2);
        y=E'*[y1;y2];
    end
end
function y=MLRsolve(x,L1,L2,m1,n,U,H)
    x1=x(1:m1,:);y1=L1'\(L1\x1);
    x2=x(m1+1:n,:);y2=L2'\(L2\x2);
    y=[y1;y2];
    y0=U*H*(U'*x);
    y=y+y0;
end
function y=BJacobisolve(x,L1,L2,m1,n,U,H)
    x1=x(1:m1,:);y1=L1'\(L1\x1);
    x2=x(m1+1:n,:);y2=L2'\(L2\x2);
    y=[y1;y2];
end