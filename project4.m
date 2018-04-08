%EE511 Project 4
%Author: Yong Wang <yongw@usc.edu>
n=100;
X = rand([n 1]);
Y = rand([n 1]);
figure;
plot(X,Y,'.');
hold on;
A = X.^2 + Y.^2;
plot(X(A<=1),Y(A<=1),'r.');
x=linspace(0,1);
y = sqrt(1 - x.^2);
plot(x,y);
legend(sprintf('%d (out)',length(find(A>1))),sprintf('%d (in)',length(find(A<=1))),'Edge','Location','NE');
xlabel("x");
ylabel('y');
epi = 4 * length(find(A<=1))/n;
title(sprintf('Estimated pi=%f',epi));
k=50;
P=zeros(k,1);
for i=1:k
	X = rand([n 1]);
	Y = rand([n 1]); 
    A = X.^2 + Y.^2;
    P(i,1)= 4 * length(find(A<1))/n;
end
figure;
histogram(P);
title(sprintf('Estimated pi, mean=%f',sum(P)/length(P)));
xlabel("pi value");
ylabel('Frequency');
%% Variance
N = [100 300 500 1000 2000 3000 5000];
V = zeros(length(N),1);
for j = 1:length(N)
    for i=1:k
        X = rand([N(j) 1]);
        Y = rand([N(j) 1]); 
        A = X.^2 + Y.^2;
        P(i,1)= 4 * length(find(A<1))/N(j);
    end
    V(j) = var(P);
    fprintf("n=%d mean=%f, var=%f\n",N(j), sum(P)/length(P),V(j));
end
figure;
bar(N,V);
title('Sample Variance');
xlabel("n (Sample number)");
ylabel('Variance of Estimates');

%% Monte Carlo integral by generating uniform points
n=1000;
x1=0.8;
x2=3.0;
y1 = 1/(1+sinh(x1*2).*log(x1));
y2 = 1/(1+sinh(x2*2).*log(x2));
X=(x2-x1).*rand(n,1)+x1;
Y=(y2-y1).*rand(n,1)+y1;
A=1./(1+sinh(X*2).*log(X));
figure;
plot(X,Y,'.');
hold on;
plot(X(A>Y),Y(A>Y),'r.');
x=linspace(0.8,3);
y=1./(1+sinh(x*2).*log(x));
plot(x,y);
legend(sprintf('%d (out)',length(find(A<Y))),sprintf('%d (in)',length(find(A>Y))),'f(x)','Location','NE');
xlabel('x');
ylabel('y');
title(sprintf('Integral = %d',abs(x2-x1)*abs(y2-y1)*length(find(A>Y))/n));
N=[1000];
V=zeros(length(N),1);
for j=1:length(N)
    int = zeros(k,1);
    for i=1:k
        n=N(j);
        x1=0.8;
        x2=3.0;
        y1 = 1/(1+sinh(x1*2).*log(x1));
        y2 = 1/(1+sinh(x2*2).*log(x2));
        X=(x2-x1).*rand(n,1)+x1;
        Y=(y2-y1).*rand(n,1)+y1;
        A=1./(1+sinh(X*2).*log(X));
        int(i,1) = abs(x2-x1)*abs(y2-y1)*length(find(A>Y))/n;
    end
    V(j) = var(int);
    fprintf("1: Mean=%f, Variance=%f\n",sum(int)/k, V(j));
end

n=1000;
x1=-pi;
x2=pi;
y1 = -pi;
y2 = pi;
z1 = exp(-1*0^4 - 0^4);
z2 = exp(-1*pi^4 - pi^4);
X=(x2-x1).*rand(n,1)+x1;
Y=(y2-y1).*rand(n,1)+y1;
Z=(z2-z1).*rand(n,1)+z1;
A=exp(-1*(X.^4) - Y.^4);
figure;
scatter3(X,Y,Z,'.');
hold on;
scatter3(X(A>Z),Y(A>Z),Z(A>Z),'r.');
x=linspace(x1,x2);
y=linspace(y1,y2);
surf(x,y,exp(-1*(x.^4) - y'.^4));
legend(sprintf('%d (out)',length(find(A<Z))),sprintf('%d (in)',length(find(A>Z))),'f(x,y)','Location','NE');
xlabel('x');
ylabel('y');
int = abs(x2-x1)*abs(y2-y1)*abs(z2-z1)*length(find(A>Z))/n;
title(sprintf('Integral = %f',int));
figure;
contour(x,y,exp(-1*(x.^4) - y'.^4));
hold on;
plot(-pi/2,y,'b.',pi/2,y,'b.');
title('Stratified region');
N=[1000];
V=zeros(length(N),1);
for j=1:length(N)
    int = zeros(k,1);
    for i=1:k
        n=N(j);
        X=(x2-x1).*rand(n,1)+x1;
        Y=(y2-y1).*rand(n,1)+y1;
        Z=(z2-z1).*rand(n,1)+z1;
        A=exp(-1*X.^4 - Y.^4);
        int(i,1) = abs(x2-x1)*abs(y2-y1)*abs(z2-z1)*length(find(A>Z))/n;
    end
    V(j) = var(int);
    fprintf("2: Mean=%f, Variance=%f\n",sum(int)/k, V(j));
end
%% Stratification
v = [0.8 1.2;1.2 1.6;1.6 2.4;2.4 3.0];
s = [500 300 150 50];
k=50;
int = zeros(k,1);
for i=1:k
    r = 0;
    for j=1:length(s)
        n = s(j);
        x1=v(j,1); x2=v(j,2);
        x=(x2-x1).*rand(n,1)+x1;
        y=1./(1+sinh(x*2).*log(x));
        r = r + sum(y.*(x2-x1))/n;
        if j == 1
            tt=y;
        end
    end
    int(i) = r;
end
fprintf("Stratification 1: integral=%f, var=%f\n", sum(int)/k, var(int));
v = [-pi -pi/2;-pi/2 pi/2;pi/2 pi];
s = [100 800 100];
k=50;
int = zeros(k,1);
for i=1:k
    r = 0;
    for j=1:length(s)
        n = s(j);
        x1=v(j,1); x2=v(j,2);
        y1=-pi;y2=pi;
        x=(x2-x1).*rand(n,1)+x1;
        y=(y2-y1).*rand(n,1)+y1;
        z=exp(-1*x.^4 - y.^4);
        r = r + sum(z.*(x2-x1)*(y2-y1))/n;
        if j == 1
            tt=z;
        end
    end
    int(i) = r;
end
fprintf("Stratification 2: integral=%f, var=%f\n", sum(int)/k, var(int));
%% Importance sampling
n=1000;
for i=1:k
    x=exprnd(1/3,n,1)+0.8;
    f=1./(1+sinh(x*2).*log(x));
    g=exppdf(x-0.8,1/3);
    int(i)=sum(f./g)/n;
end
fprintf("Importance sampling 1: integral=%f, var=%f\n", sum(int)/k, var(int));
n=1000;
mu=[0 0];
sigma=[0.5 0;0 0.5];
for i=1:k
    x=mvnrnd(mu,sigma,n);
    f=exp(-1*x(:,1).^4 - x(:,2).^4);
    g=mvnpdf(x,mu,sigma);
    int(i)=sum(f./g)/n;
end
fprintf("Importance sampling 2: integral=%f, var=%f\n", sum(int)/k, var(int));
%% last integral
x=linspace(-5,5);
y=linspace(-5,5);
[X,Y] = meshgrid(x,y);
Z=20+X.^2+Y.^2-10*(cos(2*pi*X)+cos(2*pi*Y));
k=50;
n=10;
int=zeros(k,1);
for t=1:k
    r=0;
    for i=2:size(X,1)
        for j=2:size(Y,2)
            dx = X(i,j)-X(i,j-1);
            dy = Y(i,j)-Y(i-1,j);
            x=dx*rand(n,1)+X(i,j-1);
            y=dy*rand(n,1)+Y(i-1,j);
            z=20+x.^2+y.^2-10*(cos(2*pi*x)+cos(2*pi*y));
            r = r + sum(z*dx*dy)/n;
            %fprintf("dx=%f dy=%f, %f-%f %f-%f %d %d\n",dx,dy,X(i,j),X(i-1,j),Y(i,j),Y(i,j-1),i,j);
        end
    end
    fprintf("Integral=%f\n",r);
    int(t)=r;
end
fprintf("Integral=%f Var=%f n=%d\n",sum(int)/k,var(int),n);
%% draw some plot for report document
figure;
x=linspace(0.8,3);
y=1./(1+sinh(x*2).*log(x));
plot(x,y);
hold on;
y=exppdf(x-0.8,1/3);
plot(x,y,'--');
title('Evaluate g(x) for importance sampling');
x=linspace(-pi,pi);
y=linspace(-pi,pi);
figure;
contour(x,y,exp(-1*(x.^4) - y'.^4));
title('Evaluate g(x,y) for importance sampling');
hold on;
mu = [0 0]; 
SIGMA = [0.5 0; 0 0.5]; 
[X,Y] = meshgrid(x,y);
z=mvnpdf([X(:) Y(:)], mu, SIGMA);
z=reshape(z,length(x),length(y));
contour(x,y,z,'r--');
figure;
x=linspace(-5,5);
y=linspace(-5,5);
[X,Y] = meshgrid(x,y);
Z=20+X.^2+Y.^2-10*(cos(2*pi*X)+cos(2*pi*Y));
surf(x,y,Z);
