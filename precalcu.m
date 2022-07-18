function[length, E, A] = precalcu(NE, ELECON, CORD, MAT)
for i = 1:1:NE
    n1 = ELECON(i,1);
    n2 = ELECON(i,2);
    length(i,1) = CORD(n2,1) - CORD(n1,1);
    E(i,1) = MAT(ELECON(i, 3), 1);
    A(i,1) = ELECON(i, 4);
end
end