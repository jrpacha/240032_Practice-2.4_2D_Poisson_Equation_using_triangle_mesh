function [Ke,Fe] = linearTriangElementShort(coeff,nodes,elem,e)
%For triangular elements, returns the element stiff matrix, 
%Ke, and the Fe vector for each element depending on the 
%coeff vector for the model equation.
%
% coeff: coefficient vector = [a11,a12,a21,a22,a00,f] for 
%        the model equation.
% nodes: matrix with the coordinates of the nodes.
%  elem: connectivity matrix defining the elements.
%     e: number of the present element

%Compute Ke, Fe for each element
%
v1=nodes(elem(e,1),:);
v2=nodes(elem(e,2),:);
v3=nodes(elem(e,3),:);

A = 0.5*det([1,v1;1,v2;1,v3]);

beta1=v2(2)-v3(2);
beta2=v3(2)-v1(2);
beta3=v1(2)-v2(2);
beta=[beta1,beta2,beta3];

gamma1 = v3(1)-v2(1);
gamma2 = v1(1)-v3(1);
gamma3 = v2(1)-v1(1);
gamma=[gamma1,gamma2,gamma3];
%
% Element stiff matrix & element forces
%
Ke = zeros(3);
Fe = zeros(3,1);

if coeff(1) ~= 0
    Ke = Ke + 0.25*coeff(1)*(beta'*beta)/A;
end

if coeff(2) ~= 0
    Ke = Ke + 0.25*coeff(2)*(beta'*gamma)/A;
end

if coeff(3) ~= 0
    Ke = Ke + 0.25*coeff(3)*(gamma'*beta)/A;
end

if coeff(4) ~= 0
    Ke = Ke + 0.25*coeff(4)*(gamma'*gamma)/A;
end

if coeff(5) ~= 0
    Ke = Ke + coeff(5)*A*[2,1,1;1,2,1;1,1,2]/12;
end
    
if coeff(6) ~= 0
    Fe = coeff(6)*A*[1; 1; 1]/3;
end


end