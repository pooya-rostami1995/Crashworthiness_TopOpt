function [variable2,N,vivi]=comp_ini
% Component geometry initialization
ini_val = [0.38 0.05 0.7];
x_int = 0.5;y_int = 0.5;
x0=x_int/2:x_int:2;    % x-coordinates of the centers of components
y0=y_int/2:y_int:1;    % y-coordinates of the centers of components
xn=length(x0);               % number of component groups in x direction               
yn=length(y0);               % number of component groups in y direction
x0=kron(x0,ones(1,2*yn));                                         
y0=repmat(kron(y0,ones(1,2)),1,xn);   
N=length(x0);                % total number of components in the design domain
L=repmat(ini_val(1),1,N);                    % vector of the half length of each component
t1=repmat(ini_val(2),1,N);                   % vector of the half width of component at point A
st=repmat([ini_val(3) -ini_val(3)],1,N/2);   % vector of the sine value of the inclined angle of each component 
variable=[x0;y0;L;t1;st];
variable2=variable(:,[1:N/2]);
Vi14=[2-x0(1);y0(1);L(1);t1(1);-st(1)];
Vi13=[2-x0(2);y0(2);L(2);t1(2);-st(2)];
Vi16=[2-x0(3);y0(3);L(3);t1(3);-st(3)];
Vi15=[2-x0(4);y0(4);L(4);t1(4);-st(4)];
Vi10=[2-x0(5);y0(5);L(5);t1(5);-st(5)];
Vi9=[2-x0(6);y0(6);L(6);t1(6);-st(6)];
Vi12=[2-x0(7);y0(7);L(7);t1(7);-st(7)];
Vi11=[2-x0(8);y0(8);L(8);t1(8);-st(8)];
Vi=[Vi9,Vi10,Vi11,Vi12,Vi13,Vi14,Vi15,Vi16];
vivi=[variable2,Vi];
vivi=vivi(:);
end