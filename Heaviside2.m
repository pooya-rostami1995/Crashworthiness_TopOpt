function H2=Heaviside2(Phi_max)
num1=find(Phi_max>=0);
H2(num1)=1;
num2=find(Phi_max<0);
H2(num2)=0;
end