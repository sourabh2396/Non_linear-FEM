function[LM] = dofmat(NE, NEN, NDOF, ELECON)

for i = 1:1:NE
    for j = 1:1:NEN
        node_no = ELECON(i,j);
        for k = 1:1:NDOF
            l = (j-1)*NDOF + k;
            LM(i,l) = (node_no - 1)*NDOF + k;
        end
    end
end

end