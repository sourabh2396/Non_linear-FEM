function[Kg_tan_iter] = global_tan_stiffness(NDOF, NEN, NN, NE, length, E, A, LM, BL, DoF, Nu_X)

Kg_tan_iter = zeros(NN*NDOF, NN*NDOF);

for i = 1:1:NE
    
    lengthe = length(i);
    
    BLe = BL(lengthe);
    
    DoFe = DoF(LM(i,:));
    Nu_Xe = Nu_X(lengthe);
    du_dXe = Nu_Xe*DoFe;
    
    BNLe = du_dXe*Nu_Xe;
    
    EXX = du_dXe + 0.5*du_dXe^2;
    
    SXX = E(i)*EXX;
    
    N = A(i)*SXX;
    
    KeLX = transpose(BLe + BNLe)*E(i)*A(i)*(BLe + BNLe);
    KeL = double(int( KeLX, 0 , lengthe ));
    
    KeNLX = transpose(Nu_Xe)*N*Nu_Xe;    
    KeNL = double(int( KeNLX, 0, lengthe ));
    
    Ke = KeL + KeNL;
    
    for j = 1:1:NEN*NDOF
        jg = LM(i,j);
        for m = 1:1:NEN*NDOF
            mg = LM(i,m);
            Kg_tan_iter(jg,mg) = Kg_tan_iter(jg,mg) + Ke(j,m);
        end
    end
end

end