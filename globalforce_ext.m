function[Fg_ext] = globalforce_ext(NDOF, NL, NN, PLOAD, UDL, NEN, LM, length, NEUDL, Nu)

Fg_ext = zeros(NN*NDOF, 1);

for i = 1:1:NEUDL
    elem_no = UDL(i,1);
    qe = UDL(i,2);
    lengthe = length(elem_no);
    
    Nue = Nu(lengthe);
    FeX = transpose(Nue)*qe;
    
    Fe = double(int( FeX, 0, lengthe ));
    
    for j = 1:1:NEN*NDOF
        jg = LM(elem_no,j);
        Fg_ext(jg,1) = Fg_ext(jg,1) + Fe(j,1);
    end
end

for j = 1:1:NL
    jg = (PLOAD(j,1)-1)*NDOF + PLOAD(j,2);
    Fg_ext(jg,1) = Fg_ext(jg,1) + PLOAD(j,3);
end

end