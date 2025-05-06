function [intrusion_max,A1]=LS_dyn_FEm(newsol,S1,S2,S)
    global NFE;
    if isempty(NFE)
        NFE=1;
    end
%% Finding elements that should be deleted
    Var_num = 5;nelx = 79;nely = 19;
    M = [ nely + 1 , nelx + 1 ];EW = 2 / nelx;  % length of element
    EH = 1 / nely;  % width of element
    [ x ,y ] = meshgrid( EW * [ 0 :  nelx] , EH * [ 0 : nely]);
    [~,N] = comp_ini;   % initialization of components
    p = 6;
    [~,Phi_max]=form_comp(Var_num,newsol,p,N,x,y,nelx,nely);
    Rho_x = Heaviside2(Phi_max);
    Rho_x=flip(Rho_x);A1=sum(sum(Rho_x))/(nelx*nely);
    Rho_x = reshape(Rho_x,M(1),M(2));
    %% preparing elements that should be deleted
   AA = reshape(1:(1+nelx)*(1+nely),(1+nelx),(1+nely));
   CC = flip(AA');final_matrix = [CC(:),Rho_x(:)];
   secondColumn = final_matrix(:,2);
   rows_zero = final_matrix(secondColumn==0, :);
   sh_del_elm = rows_zero(:,1);sh_del_elm=sort( sh_del_elm);
   %% Reading initial key file information and editing main matrix
   InpBas='00run2';
   fid = fopen([InpBas '.key']);
   Rows = textscan(fid, '%s', 'delimiter', '\n');
   Rowsl=Rows{1,1};
   gg=find(contains(S,'$#    n1'));
   gg2=find(contains(S,'*NODE'));
   Rows_change=Rowsl(gg+1:gg2-1);
   fid2 = fopen( 'text2.txt', 'wt' );
for i=1:size(Rows_change,1)
    fprintf( fid2, Rows_change{i});
    fprintf(fid2,'\n');
end
ed2=readtable('text2.txt','Delimiter',' ');
ed2=table2array(ed2);
for i=1:size(Rows_change,1)
    eddd=~isnan(ed2(i,:));
    eddf=ed2(i,:);
    ed3(i+1,:)={eddf(eddd)};
end
ed3{1,1}=[12,6];
for  i=1:size(Rows_change,1)/2
    paiirs{i,2}=ed3{2*i};
    paiirs{i,1}=ed3{2*i-1};
end
sh_del_elm_rev=sh_del_elm +1659;
for  i=1:size(paiirs,1)
    ABB=paiirs{i,1};
    if ismember(ABB(1),sh_del_elm_rev)==1
            paiirs{i,1}=[];
            paiirs{i,2}=[];
    end
end
paiirs_new=paiirs;empty = cellfun('isempty',paiirs_new);
paiirs_new(all(empty,2),:) = [];
LL1=cell2mat(paiirs_new(:,1));LL2=cell2mat(paiirs_new(:,2));LL1=LL1(2:end,:);
for i=1:size(LL1,1)
        LL3(2*i,1:2)=LL1(i,:);  
end
for i=1:size(LL1,1)+1
        LL3(2*i-1,1:10)=LL2(i,:);   
end
%% Preparing new input deck for LS-Dyna analysis
mm='00final_run.key';fid3 = fopen( mm, 'wt' );
% Writing stable section part 1
for i=1:size(S1,1)
    fprintf( fid3, S1(i));
    fprintf(fid3,'\n');
end
% write changed middle section
for i=1:size(LL3,1)
    for j=1:10
        fprintf(fid3,'    ');
        LL4=num2str(LL3(i,j));  
        if numel(LL4)==1
            fprintf(fid3,'   ');
        end
        if numel(LL4)==2
            fprintf(fid3,'  ');
        end
        if numel(LL4)==3
            fprintf(fid3,' ');
        end
        fprintf(fid3,LL4);
    end
   fprintf(fid3,'\n');
end
% Writing stable section part 2
for j=1:size(S2,1)
    fprintf( fid3, S2(j));
    fprintf(fid3,'\n');
end
fprintf(fid3,'\n');
fclose all;
%% Running LS_Dyna
inputFile='00final_run.key';
dynaExe = 'C:\LSDYNA\program\ls-dyna_smp_d_R11_0_winx64_ifort131.exe';
folder='C:\Users\Admin\Desktop\ACOR_crash_epsilon';
baseExecStr = sprintf('cd /d "%s" & %s I=%s O=d3hsp',...
                folder, dynaExe, inputFile);
disp(1)
system(baseExecStr);
disp(2)
%% Reading the nodout file for intrusion value
   OutBas='nodout';fid = fileread(OutBas);
   Rows_res = textscan(fid, '%s', 'delimiter', '\n');
   Rows_final=Rows_res{1,1};
   read_row=15:7:size(Rows_final,1);
   for i=1:size(read_row,2)
       Row_1=(Rows_final(read_row(i),:));
       rowss{i,1}=str2num(Row_1{1,1});
   end
   for i=1:size(read_row,2)
       LKK=rowss{i,1};
       disppp(i,1)=LKK(1,2);
   end
   %% Calculating intrusion
    if size(disppp,1)==21
        intrusion_max=abs(min(disppp));
    else
        intrusion_max=1000000000000000000;
    end
    if size(disppp,1)==21
        if min(disppp)==disppp(21,1)
            intrusion_max=1000000000000000000;
        end
    end
   fclose all;
   NFE=NFE+1;
end