function[Fg_int_iter] = globalforce_int(NDOF, NEN, NN, NE, length, E, A, LM, BL, DoF, Nu_X)

Fg_int_iter = zeros(NN*NDOF, 1);

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
    
    FeX = -transpose(BLe + BNLe)*N;
    
    Fe = double(int( FeX, 0, lengthe ));
    
    for j = 1:1:NEN*NDOF
        jg = LM(i,j);
        Fg_int_iter(jg,1) = Fg_int_iter(jg,1) + Fe(j,1);
    end
end

end