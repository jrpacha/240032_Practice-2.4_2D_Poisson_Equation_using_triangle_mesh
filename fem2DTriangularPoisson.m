clearvars
close all

%Define geometry
nodes=[
    0,0;
    0.5,0;
    0.5,0.5;
    1,0;
    1,0.5;
    1,1;
    ];
elem=[1,2,3;
      5,3,2;
      2,4,5;
      3,5,6;
      ];

numNod=size(nodes,1);
numElem=size(elem,1);
numbering=1;
plotElements(nodes, elem, numbering);

%Define Coefficients vector of the model equation:
%In this case we use the Poisson coefficients defined
%in the problem above

a11=1;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=1;
coeff=[a11,a12,a21,a22,a00,f];

%Compute the global stiff matrix
u=zeros(numNod,1);
K=zeros(numNod);
Q=zeros(numNod,1);
F=zeros(numNod,1);

for e=1:numElem
    [Ke, Fe] = linearTriangElem(coeff,nodes,elem,e);
    % 
    %Assemble the elements
    %
    rows=[elem(e,1); elem(e,2); elem(e,3)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke;
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end
fixedNods=[4,5,6];
freeNods=setdiff(1:numNod,fixedNods);
%Boundary conditions (B.C.)
%---------- Essential B.C.
u(fixedNods)=0.0;
%---------- Natural B.C.
u(freeNods)=0.0;

%Reduced system
Fm=F(freeNods)-K(freeNods,fixedNods)*u(fixedNods);
Fm=Fm+Q(freeNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Fm;
u(freeNods)=um;

%PostProcess: Compute secondary variables and plot
%results
Q=K*u-F;

fprintf('\n%6s%8s%12s%12s%11s\n',...
    '#Nod.','X','Y','U','Q')
fprintf('%4d%14.4e%12.4e%12.4e%12.4e\n',...
    [[1:numNod]',nodes,u,Q]')

title='Poisson Solution';
colorScale='jet';
plotContourSolution(nodes,elem,u,title,colorScale);
