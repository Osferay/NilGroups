
g := Random(G);
h := Random(G);
for i in [1..100] do
    IsConjugateNilGroup(G, g, g^h);
od;