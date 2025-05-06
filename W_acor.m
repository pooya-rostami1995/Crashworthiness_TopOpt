clc; clear; close all
global NFE; NFE = 1;Loop=1;
volfraca=[0.4,0.5,0.6,0.7,0.8];
[variable,N] = comp_ini;   % initialization of components
%% Creating stable lines .key files
S = readlines('00run2.key');
data=find(contains(S,'$#    n1'));
S1=S;S1(data+1:end)=[];
data2=find(contains(S,'*NODE'));
S2=S;S2(1:data2-1)=[];
for bb=1:size(volfraca,2)
volfrac=volfraca(bb);
%% Parameters for MMC
xy00 = variable(:);xval = xy00;
xmin = [0 ; 0  ; 0.01 ; 0.01 ; -1.0];
xmin = repmat(xmin,N/2,1);
xmax = [2 ; 1 ; 0.9 ; 0.07 ; 1.0];
xmax = repmat(xmax,N/2,1);  %number of constraint
Var_num=5;  % number of design variables for each component
%% ACOR Parameters
CostFunction=@(newsol,volfrac,S1,S2,S) Costfunf2(newsol,volfrac,S1,S2,S);   % Cost Function
nVar=Var_num*N/2;             % Number of Decision Variables
VarSize = [1 N*Var_num/2]; % Unknown Variables Matrix Size
nIndividual = 60;          % Individualulation Size
nSample=50;         % Sample Size
q=0.5;              % Selection pressure parameter
zet_param=0.8;        % Deviation-Distance Ratio
%% Initialization
% Create Empty Individual Structure
empty_individual.Vectt=[];
empty_individual.Cost=[];
% Create Individualulation Matrix
Individual=repmat(empty_individual,nIndividual,1);
best_res.Vectt = xval';
best_res.Cost = CostFunction(best_res.Vectt,volfrac,S1,S2,S);
%% Initialize Individuals
Vect=[20,20,20,20,20]; Vect2=repmat(Vect,1,N/2); 
for i = 1:nIndividual 
    Individual(i).Vectt = xy00' + ((xmax-xmin)'./Vect2)*randn;
    % Apply Variable Limits
    Individual(i).Vectt=max(Individual(i).Vectt,xmin'); % Clipping for lower bound
    Individual(i).Vectt=min(Individual(i).Vectt,xmax'); % Clipping for upper bound
    Individual(i).Cost = CostFunction(Individual(i).Vectt,volfrac,S1,S2,S);
    if Individual(i).Cost < best_res.Cost  
        best_res = Individual(i);
    end
end
[~, S_order]=sort([Individual.Cost]); % Sort Individuals
Individual=Individual(S_order);
% Solution Weights
w=1/(sqrt(2*pi)*q*nIndividual)*exp(-0.5*(((1:nIndividual)-1)/(q*nIndividual)).^2);
% Selection Probabilities
pp=w/sum(w);
%% ACOR Main Loop
while   NFE<10000
%% Evaluate
    [intrusion_max,A1]=LS_dyn_FEm(xy00,S1,S2,S); 
%% Main Loop   
    % Means
    s=zeros(nIndividual,nVar);
    for l=1:nIndividual
        s(l,:)=Individual(l).Vectt;
    end
    % Standard Deviations
    S_deviation=zeros(nIndividual,nVar);
    for l=1:nIndividual
        D=0;
        for r=1:nIndividual
            D=D+abs(s(l,:)-s(r,:));
        end
        S_deviation(l,:)=zet_param*D/(nIndividual-1);
    end
    %% Construction of statistical model
    newIndividual=repmat(empty_individual,nSample,1);
    for t=1:nSample
        newIndividual(t).Vectt=zeros(VarSize);
        for i=1:nVar
            rr=rand;  
            C=cumsum(pp);
            l=find(rr<=C,1,'first');
            newIndividual(t).Vectt(i)=s(l,i)+S_deviation(l,i)*randn; % Generate Gaussian Random Variable
        end
        newIndividual(t).Vectt=max(newIndividual(t).Vectt,xmin'); %Clipping to lower limit
        newIndividual(t).Vectt=min(newIndividual(t).Vectt,xmax'); % Clipping to upper limit
        newIndividual(t).Cost=CostFunction(newIndividual(t).Vectt,volfrac,S1,S2,S);  
    end
    Individual=[Individual
         newIndividual]; % Merge individuals
    %% Sort Individuals
    [~, S_order]=sort([Individual.Cost]);
    Individual=Individual(S_order);
    Individual=Individual(1:nIndividual);     % Delete other individuals
    %% Update Best Solution Ever Found
    best_res=Individual(1);
    xval = best_res.Vectt;
    xy00 = round(xval*1e4)/1e4;
    disp([' It.: ' sprintf('%4i\t',Loop) ' Obj.: ' sprintf('%6.3f\t',intrusion_max) ' Vol.: ' ...
        sprintf('%6.4f\t',A1) 'NFE.: ' sprintf('%4i\t',NFE)]);
    Loop = Loop+1;
    zet_param=zet_param*0.99;
end
f0vala(bb) = intrusion_max;fvala(bb) = A1;nfe(bb)=NFE;xvala(bb).Vectt=xy00;
end