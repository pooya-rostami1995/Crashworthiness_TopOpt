function [Phi,Phi_max]=form_comp(Var_num,xy00,p,N,x,y,nelx,nely)
    ccc=reshape(xy00,5,N/2);x0=ccc(1,:);y0=ccc(2,:);L=ccc(3,:);t1=ccc(4,:);
    st=ccc(5,:);Vi14=[2-x0(1);y0(1);L(1);t1(1);-st(1)];
    Vi13=[2-x0(2);y0(2);L(2);t1(2);-st(2)];Vi16=[2-x0(3);y0(3);L(3);t1(3);-st(3)];
    Vi15=[2-x0(4);y0(4);L(4);t1(4);-st(4)];Vi10=[2-x0(5);y0(5);L(5);t1(5);-st(5)];
    Vi9=[2-x0(6);y0(6);L(6);t1(6);-st(6)];Vi12=[2-x0(7);y0(7);L(7);t1(7);-st(7)];
    Vi11=[2-x0(8);y0(8);L(8);t1(8);-st(8)];Vi=[Vi9,Vi10,Vi11,Vi12,Vi13,Vi14,Vi15,Vi16];
    vivi=[ccc,Vi];vivi=vivi(:);LSgrid.x = x(:);LSgrid.y = y(:);     % coordinate of nodes
    %Forming Phi^s
    for i = 1:N  
        Phi{i} = tPhi(vivi(Var_num*i-Var_num+1:Var_num*i),LSgrid.x,LSgrid.y,p); 
    end
    %Union of components
    tempPhi_max = Phi{1};
    for i = 2:N 
        tempPhi_max = max(tempPhi_max,Phi{i});  
    end
    Phi_max = reshape(tempPhi_max,nely+1,nelx+1);
end