function GFLS262(lx,ly,Vr,volReq,maxedge,minedge,E,nu,tau,sigma,alpha)
%% CREATE DESIGN DOMAIN
BDY = [-0.5*(lx), -0.5*(ly); 0.5*(lx), 0.5*(ly)];
[xn,yn] = meshgrid(BDY(1,1):BDY(2,1), BDY(1,2):BDY(2,2));
%% INITIAL MESHING
phi = -(sin(xn/BDY(2,1)*4.8*pi).*cos(yn/BDY(2,1)*4.8*pi)+0.5); phi(1:4,:)= -1;
d1 = 0.4; d2 = 1;
[c,Curv] = ContourPoints(contour(xn,yn,phi,[0 0]),d1,d2);
[p,t,t1,t2,Ve,pmid] = GenerateMesh(xn,yn,c,phi,maxedge,minedge,1.2,BDY,80);
%% MAIN LOOP
for iterNum = 1:200
    %% FINITE ELEMENT ANALYSIS AND SENSITIVITY ANALYSIS
    [Ce,J,vol] = FEA(t,t1,t2,p,Ve,BDY,E,nu);
    Curvfull = zeros(length(p),1);
    Curvfull(1:length(c)) = Curv(1:length(c));
    Obj(iterNum) = J; % OBJECTIVE FUNCTION VALUE
    volt(iterNum) = vol; % VOLUME OF SOLID MATERIAL
    Cn = E2N(t,p,-Ce,Ve); % INTERPOLATE SHAPE DERIVATIVE
    Vc = -(Cn(1:length(p))+tau*Curvfull); % VELOCITY CALCULATION
    dt = alpha/max(abs(Cn(:))); % STEP LENGTH
    V = max(volReq ,vol-Vr); % TARGET VOLUME IN CURRENT ITERATION
    %% PRINT RESULTS
    disp(['It.: ' num2str(iterNum) '  Obj.: ' sprintf('%6.4f',J) '  Vol.: ' sprintf('%6.4f' ,vol) ]); clf;
    patch('Faces',t1,'Vertices',p,'EdgeColor',[250 250 250]/255,'FaceColor',[92 158 173]/255);
    hold on; axis off equal tight
    patch('Faces',t2,'Vertices',p,'EdgeColor','k','FaceColor',[239 111 108]/255);
    contour(xn,yn,phi,[0 0],'linewidth',2,'color',[62 43 109]/255);
    saveas(gcf,['./fig', int2str(iterNum) '.png']);
    %% CONVERGENCE CHECK
    if iterNum>=20 && abs(vol-volReq)<0.005 && all(abs(Obj(end)-Obj(end-5:end-1))< 0.005*abs(Obj(end))) && Obj(end)<Obj(end-1)
        return;
    end
    %% LEVEL SET FUNCTION REINITIALIZATION & DISCRETE GRADIENT DERIVATION
    x = ones(length(t),1)';
    x(1:length(t1)) = 0;
    phiE = (~x).*(pdist(pmid,t2,t1))-x.*(pdist(pmid,t1,t1));
    phiN = E2N(t,p,phiE,Ve);
    [Gradx,Grady] = GradDeri(phiN,p,t,Ve,p);
    %% DESIGN UPDATE USING BISECTION METHOD
    l1 = 0; l2 = 1e9;
    while abs(l2-l1)/abs(l2) > 1e-9
        Vbs = Vc - (l1+l2)/2;
        phiNnew = phiN - dt * Vbs.* ((abs(Gradx)).^2+(abs(Grady)).^2).^0.5;
        phinew = imgaussfilt(griddata(p(:,1),p(:,2),phiNnew,xn,yn,'cubic'), sigma); % GAUSSIAN FILTERING
        [c,Curv] = ContourPoints(contour(xn,yn,phinew,[0 0]),d1,d2);
        if isempty(c)==0
            [~,~,t1,~,Ve1,~,] = GenerateMesh(xn,yn,c,phinew,maxedge,minedge, 1.2,BDY,1);
            volCurr = 1-sum(Ve1(1:length(t1)))/sum(Ve1);
        else
            volVoid = dot(max(0,sign(mean(phiNnew(t),2))),Ve');
            volCurr = (sum(Ve)-volVoid)/sum(Ve);
        end
        if volCurr > V
            l1 = (l1+l2)/2;
        else
            l2 = (l1+l2)/2;
        end
    end
    %% REMESHING THE DESIGN DOMAIN
    phi = phinew;
    [p,t,t1,t2,Ve,pmid] = GenerateMesh(xn,yn,c,phi,maxedge,minedge,1.2, BDY,80);
end
%% SUBFUNCTIONS:
%% FIND POINTS ON THE CONTOUR
function [c,Curv] = ContourPoints(c,d1,d2)
tol = 1e-12; num = 1; col = 1;
while col < size(c,2)
    idx = col+1:col+c(2,col);
    s(num).x = c(1,idx);
    s(num).y = c(2,idx);
    s(num).isopen = abs(diff(c(1,idx([1 end])))) >tol || abs(diff(c(2,idx([1 end])))) >tol;
    num = num+1;
    col = col+c(2,col)+1;
end
c = []; Curv = [];
for k = 1:num-1
    ct = [s(k).x; s(k).y];
    if length(ct)>4
        ct1 = ct;
        ndist = sqrt(sum(diff(ct,1,2).^2,1));
        for i = 1:size(ndist,2)
            if  ndist(i) < d1
                ct1(:,i) = [0;0];
                ct1(:,(i+1)) = 0.5*(ct(:,i)+ct(:,(i+1)));
            end
        end
        ct1(:,all(ct1==0,1)) = [];
        if  sqrt(sum((ct1(:,1)-ct1(:,end)).^2,1)) < d1
            ct1(:,end) = [];
        end
        if size(ct1,2)>2
            if s(k).isopen == 0
                ct1=[ct1 ct1(:,1)];
            end
            ndist = sqrt(sum(diff(ct1,1,2).^2,1));
            ct=ct1(:,1);
            for i = 1:size(ndist,2)
                if  ndist(i) > d2
                    ct = [ct 0.5*(ct1(:,i)+ct1(:,(i+1))) ct1(:,i+1)];
                else
                    ct = [ct ct1(:,i+1)];
                end
            end
            if s(k).isopen == 1
                c = [c; ct'];
                Curv = [Curv; CalculateCurv(ct')];
            else
                c = [c; ct(:,1:end-1)'];
                Curv = [Curv; CalculateCurv(ct(:,1:end-1)')];
            end
        end
    end
end
%% BODY FITTED MESH GENERATOR
function [p,t,t1,t2,Ve,pmid] = GenerateMesh(xn,yn,c,dN,maxedge,minedge, fscale,BDY,maxiter)
x = xn;
x(2:2:end,:) = x(2:2:end,:)+0.5;
pi = [x(:),yn(:)];
d = zeros(size(pi,1),1);
for i = 1:size(pi,1)
    d(i) = sqrt(min((pi(i,1)-c(:,1)).^2+(pi(i,2)-c(:,2)).^2));
end
r0 = 1./min(max(minedge,d),maxedge).^2;
pfix=[c; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),0];
rv = 0.5;
p = [pfix; pi(r0./max(r0)>rv,:)];
p = unique(p,'rows','stable');
p1 = 0;
beta = 0.2;
t = delaunayn(p);
for i = 1:maxiter
    if max(sum((p-p1).^2,2))>1e-6
        t = delaunayn(p);
        edges = unique(sort([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],2),'rows');
        p1 = p;
    end
    midpoint = (p(edges(:,1),:)+p(edges(:,2),:))/2;
    d = zeros(size(midpoint,1),1);
    for j = 1:size(midpoint,1)
        d(j) = sqrt(min((midpoint(j,1)-c(:,1)).^2+(midpoint(j,2)-c(:,2)).^2));
    end
    L = sqrt(sum((p(edges(:,1),:)-p(edges(:,2),:)).^2,2));
    L1 = min(max(minedge,d),maxedge);
    L0 = fscale*L1*sqrt(sum(L.^2)/sum(L1.^2));
    Fb = max(L0-L,0)./L *[1,1].*(p(edges(:,1),:)-p(edges(:,2),:));
    Fp = full(sparse(edges(:,[1,1,2,2]),ones(size(d))*[1,2,1,2],[Fb,-Fb],size(p,1),2));
    Fp(1:size(pfix,1),:) = 0;
    p = p+beta*Fp;
    p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1)));
    p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2)));
end
[p, t] = Uniquenode(p,t);
pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
dE = interp2(xn,yn,dN,pmid(:,1),pmid(:,2),'cubic');
t1 = t(dE>0,:);
t2 = t(dE<=0,:);
t = [t1;t2];
pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
for kk = 1:length(t)
    Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]);
end
%% REMOVE REPEATED NODES
function [P,T] = Uniquenode(p,t)
for kk = 1:length(t)
    Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]);
end
t((Ve==0),:) = [];
P = unique(p(unique(sort(t(:))),:),'rows','stable');
for i = 1:length(t)
    for j = 1:3
        T(i,j) = find(P(:,1)==p(t(i,j),1)& P(:,2)==p(t(i,j),2));
    end
end
%% CALCULATE NODE SENSITIVITY
function [dN] = E2N(t,p,x,Ve)
dN = zeros(length(p),1);
for i = 1:length(p)
    [row,~] = find(t==i);
    dN(i) = dot(Ve(row),x(row))/sum(Ve(row));
end
%% FINITE ELEMENT ANALYSIS
function [Ce,J,vol] = FEA(t,t1,t2,p,Ve,BDY,E,nu)
NT = length(t); KK = zeros(6,6*NT);
for i = length(t1)+1:NT
    KK(:,6*i-5:6*i) = GetmatrixKe(p(t(i,:),1),p(t(i,:),2),E,nu);
end
elemDof = zeros(NT,6);
elemDof(:,[1 3 5]) = 2*t-1;
elemDof(:,[2 4 6]) = 2*t;
iK = reshape(kron(elemDof,ones(6,1))',36*NT,1);
jK = reshape(kron(elemDof,ones(1,6))',36*NT,1);
sK = reshape(KK,36*NT,1);
NK = sparse(iK,jK,sK,2*length(p),2*length(p));
NK = (NK+NK')/2;
fixedNodes1 = find(p(:,1)==BDY(1,1));
fixedNodes2 = find(p(:,1)==BDY(2,1) & p(:,2)==BDY(1,2));
fixedDof = [2*fixedNodes1-1; 2*fixedNodes2];
forceNodes = find(p(:,1)==BDY(1,1) & p(:,2)==BDY(1,2));
SolidDof = [2*unique(sort(t2(:)))-1; 2*unique(sort(t2(:)))];
freeDofs = setdiff(SolidDof,fixedDof);
U = zeros(2*length(p),1);
F = sparse(2*forceNodes,1,-50,2*length(p),1);
U(freeDofs,:) = NK(freeDofs,freeDofs) \ F(freeDofs,1);
for i = length(t1)+1:NT
    Ce(i) = 0.5 .* sum((U(elemDof(i,:))'*KK(:,6*i-5:6*i)).*U(elemDof(i,:))',2);
end
J = sum(Ce);
vol = 1-sum(Ve(1:length(t1)))/sum(Ve);
%% ELEMENT STIFFNESS MATRIX
function [Ke] = GetmatrixKe(X,Y,E0,nu)
D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2];
J = [X(1)-X(3) Y(1)-Y(3);X(2)-X(3) Y(2)-Y(3)];
Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1); -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];
Ke = 1/2*det(J)*Be'*D*Be;
%% LEVEL SET FUNCTION REINITIALIZATION
function dist = pdist(pmid,input,t1)
if size(input) == size(t1)
    yinds = pmid(1:length(t1),1);
    xinds = pmid(1:length(t1),2);
else
    yinds = pmid(length(t1)+1:end,1);
    xinds = pmid(length(t1)+1:end,2);
end
for i = 1:size(pmid,1)
    dist(i) = sqrt(min((yinds-pmid(i,1)).^2+(xinds-pmid(i,2)).^2));
end
%% DISCRETE GRADIENT DERIVATION
function [Gradx,Grady] = GradDeri(f,p,t,Ve,c)
Gradx = zeros(length(p),1);
Grady = zeros(length(p),1);
for A0 = 1:length(c)
    [row , ~] = find(t==A0);
    ts = t(row,:);
    dx = zeros(length(ts),1);
    dy = zeros(length(ts),1);
    AT = Ve(row);
    for i = 1:size(ts,1)
        p1 = p(ts(i,1),:);
        p2 = p(ts(i,2),:);
        p3 = p(ts(i,3),:);
        phi1 = f(ts(i,1));
        phi2 = f(ts(i,2));
        phi3 = f(ts(i,3));
        D1 = [0 -1; 1 0]*(p1-p3)';
        D2 = [0 -1; 1 0]*(p2-p1)';
        Gradu = (phi2-phi1)/2/AT(i) .* D1+(phi3-phi1)/2/AT(i).* D2;
        dx(i) = Gradu(1)*(AT(i)/sum(AT));
        dy(i) = Gradu(2)*(AT(i)/sum(AT));
    end
    Gradx(A0) = sum(dx);
    Grady(A0) = sum(dy);
end
%% CALCULATE CURVATURE
function [C] = CalculateCurv(X)
N = size(X,1);
X = [X,zeros(N,1)];
R = NaN(N,1);
for i = 2:N-1
    D = cross(X(i-1,:)'-X(i,:)',X(i+1,:)'-X(i,:)');
    R(i) = norm(X(i-1,:)'-X(i+1,:)')*norm(X(i,:)'-X(i+1,:)')*norm(X(i,:)'-X(i-1,:)')/2/norm(D);
end
C = 1./R; C(isnan(C))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Zicheng Zhuang and Yiwei Weng            %
% Department of Building and Real Estate, Hong Kong Polytechnic University %
% Please send your comments to: zhuanginhongkong@outlook.com               %
%                                                                          %
% The program is proposed for educational purposes and is introduced in    %
% the paper - A 262-line Matlab code for the level set topology            %
% optimization based on the estimated gradient field in the body-fitted    %
% mesh, SMO, 2024                                                          %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guarantee that the code is     %
% free from errors. Furthermore, we shall not be liable in any event.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%