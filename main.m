close all;
clear all;
clc;
format shortE;

NDOF    = 1                       % Number of degrees of freedom per node

CORD    = [0:10:100]'
NN      = size(CORD,1)            % Total number of nodes
NDIM    = size(CORD,2)            % Number of coordinates per node

NEN     = 2                       % Number of nodes per element

NRC     = 1                       % Number of sectional properties

ELECON  = [1 2 1 1000;
           2 3 1 1000;
           3 4 1 1000;
           4 5 1 1000;
           5 6 1 1000;
           6 7 1 1000;
           7 8 1 1000;
           8 9 1 1000;
           9 10 1 1000;
           10 11 1 1000];            % Connectivity matrix; Global node nos, Material nos, A
NE      = size(ELECON,1);         % Number of elements

BC      = [1 1 0];                 % Boundary Condition matrix; Global Node no, DoF no, BC value
ND      = size(BC,1);              % Number of applied BCs specified

PLOAD   = [11,1,10000e3];
NL      = size(PLOAD,1)           % Number of point loads specified
      
UDL     = []                      % UDL matrix; Element no, DoF no, q
NEUDL   = size(UDL,1)             % Number of elements over which UDL is applied 

MAT     = [2000]                     % Material matrix; E
NM      = size(MAT,1)             % Number of different materials
NP      = size(MAT,2)             % Number of properties of material need to be specified

[length, E, A] = precalcu(NE, ELECON, CORD, MAT)

[LM] = dofmat(NE, NEN, NDOF, ELECON)

syms l_var X

N1 = 1 - X/l_var;
N2 = X/l_var;
Nu(l_var) = [N1 N2];

N1_X = diff(N1, X);
N2_X = diff(N2, X);
Nu_X(l_var) = [N1_X N2_X];

BL(l_var) = Nu_X;

Fg_ext = globalforce_ext(NDOF, NL, NN, PLOAD, UDL, NEN, LM, length, NEUDL, Nu)
DoF = zeros(NN,1);

load_steps = 50;
iter = 0;
Strain=zeros(load_steps,NE);
stress=zeros(load_steps,NE,1);

for ls = 1:1:load_steps
    
    iter = iter + 1
    Fg_ext_iter = Fg_ext*ls/load_steps
        BCnew=BC;
    BCnew(:,3)=BC(:,3)*ls/load_steps;
    for i=1:ND
    DoF(BCnew(i,1),1)=BCnew(i,3);
    end
 
    delta_DoF = 1*ones(NN*NDOF,1);
    sub_iter = 0;
    
    while (norm(delta_DoF) > 1e-6)
        
        sub_iter = sub_iter + 1
        
        Kg_tan_iter = global_tan_stiffness(NDOF, NEN, NN, NE, length, E, A, LM, BL, DoF, Nu_X)
        
        Fg_int_iter = globalforce_int(NDOF, NEN, NN, NE, length, E, A, LM, BL, DoF, Nu_X)
        
        Fg_iter = Fg_int_iter + Fg_ext_iter
        
        [UBCdof, FBCdof, UE, KFE, KFF, FF] = bcapplied(ND, BC, NDOF, NN, Fg_iter, Kg_tan_iter)
        
        UF = inv(KFF)*(FF);
        
        delta_DoF(UBCdof) = 0;
        delta_DoF(FBCdof) = UF;
        delta_DoF
        norm_delta_DoF = norm(delta_DoF)
        
        DoF = DoF + delta_DoF

    end
% Postprocessing

disp_mat=[1,NE];
for i = 1:1:NE

lengthe = length(i);
    
    BLe = BL(lengthe);
    N_mat=[1-0.1*CORD(i) 1+0.1*CORD(i)];
    disp_mat(ls,i)=N_mat*[DoF(i);DoF(i+1)];
    
    DoFe = DoF(LM(i,:));
    Nu_Xe = Nu_X(lengthe);
    du_dXe = Nu_Xe*DoFe;
    
    BNLe = du_dXe*Nu_Xe;
    
    EXX = du_dXe + 0.5*du_dXe^2;
    Strain(ls,i)=Strain(ls,i)+EXX;
    
    SXX = E(i)*EXX;
    stress(ls,i)=stress(ls,i)+SXX;
    
    N = A(i)*SXX;
end
end
figure 

plot((Strain(:,10)),(stress(:,10)),'b')
 xlabel('strain')
ylabel('stress')
figure
plot(disp_mat(50,:),Strain(50,:),'r')
ylabel('strain')
xlabel('displacement')
