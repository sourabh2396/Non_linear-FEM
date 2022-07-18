function[UBCdof, FBCdof, UE, KFE, KFF, FF] = bcapplied(ND, BCnew, NDOF, NN, Fg, Kg)
UBCdof = zeros(ND,1);
UE = zeros(ND,1);
for i = 1:ND
    ig = (BCnew(i,1)-1)*NDOF + BCnew(i,2);
    UBCdof(i)=ig;
    UE(ig,1) = BCnew(i,3);
    
end
FBCdof=setdiff((1:NN),UBCdof);

KFE = Kg(FBCdof, UBCdof);
KFF = Kg(FBCdof, FBCdof);
FF = Fg(FBCdof, 1);
end