%FEM HW 1
%%-utt + ut = (pi^2+1)sin(pi*t), u(0)=u(1)=0, in (0,1)
N = 40; 
h = 1/N;
x = (0:h:1)';
%build KU = F
K = diag(repmat((2/h+2*h/3),1,N))+diag(repmat((-1/h+h/6),1,N-1),1)+diag(repmat((-1/h+h/6),1,N-1),-1);
K(N,N-1) = 0; K(N,N) = 1; %form B.C.: u_n = 0;
F = zeros(N,1);

for i = 1:N-1
    F(i) = (pi^2+1)*integral(@(t)sin(pi.*t).*(t-(i-1)/N)/h,(i-1)/N,i/N)+(pi^2+1)*integral(@(t)sin(pi.*t).*((i+1)/N-t)/h,i/N,(i+1)/N);
end
F(N) = 0;
U = K\F;

%use B.C. to fullfill the solution
Ufem = zeros(N+1,1);
Ufem(1) = 0;
Ufem(2:N+1) = U;

%compute error
e = Ufem - sin(pi*x);
err1 = 0;
for i = 1:N
   err1 = err1 + (e(i+1)-e(i))^2/h;  
end

plot(x, Ufem, 'r-',x,sin(pi*x),'g-');
ymin = min(U) - 0.3;    ymax = max(U) + 0.3;
axis([0, 1, ymin, ymax]);
xlabel('t');ylabel('u');
legend('FEM method','Exact Solution' ,'Location', 'SouthWest');
title("FEM method");



%%for order of convergence, take another N
N = 80; 
h = 1/N;
x = (0:h:1)';
%build KU = F
K = diag(repmat((2/h+2*h/3),1,N))+diag(repmat((-1/h+h/6),1,N-1),1)+diag(repmat((-1/h+h/6),1,N-1),-1);
K(N,N-1) = 0; K(N,N) = 1; %form B.C.: u_n = 0;
F = zeros(N,1);

for i = 1:N-1
    F(i) = (pi^2+1)*integral(@(t)sin(pi.*t).*(t-(i-1)/N)/h,(i-1)/N,i/N)+(pi^2+1)*integral(@(t)sin(pi.*t).*((i+1)/N-t)/h,i/N,(i+1)/N);
end
F(N) = 0;
U = K\F;

%use B.C. to fullfill the solution
Ufem = zeros(N+1,1);
Ufem(1) = 0;
Ufem(2:N+1) = U;

%compute error
e = Ufem - sin(pi*x);
err2 = 0;
for i = 1:N
   err2 = err2 + (e(i+1)-e(i))^2/h;  
end

%compute order of convergence
order =floor( log(err1/err2)*log(2) );