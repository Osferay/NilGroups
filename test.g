G := MinimalTauGroups(5,4,1);
g := Random(G);
efa := PcpsOfEfaSeries(G);
C := G;
L := G;
h := Identity(G);

for i in [2..Length(efa)] do
    #We have efa[i]=Gi/G(i+1)
    #Get Gi
    N := SubgroupByIgs(G , NumeratorOfPcp( efa[i] ));
    #Calculate G -> G/Gi
    nat := NaturalHomomorphismByNormalSubgroup( G, N );

    #Check a trivial case:
    if IsAbelian(G/N) then
        C := G/N;
        L := G;
        h := h*Igs(G)[i-1]^(Exponents(g)[i-1]);

    else
    
        #Get G(i-1)
        M := SubgroupByIgs(G, NumeratorOfPcp( efa[i-1] ));

        #Calculate gGi
        hGi := Image(nat, h);

        #Calculate G(i-1)/Gi under G/Gi
        fac := Image( nat, M );

        #Calculate L/Gi under G/Gi
        LN := Image(nat, L);

        #Calculate f: L/Gi -> G(i-1)/Gi
        gens := Igs( LN );
        imgs := List( gens , x -> Comm(hGi ,x) );
        f := GroupHomomorphismByImages( LN , fac, gens, imgs );

        #Calculate the conjugate, first solve the equation of exponents

        e := [ Exponents( nat( h^-1*g ) * Igs(f(LN))[1] )[i-1] ];
        for j in [1..Length(imgs)] do
            Add(e, Exponents(imgs[j])[i-1]);
        od;
        e := GcdRepresentation( e );

        #find v such that f(v)=w
        t := Identity( Image(nat) );
        for j in [2..Length(e)] do
            t := t*gens[j-1]^e[j];
        od;

        #find l such that nat(l)=v
        l := Identity(G);
        for j in [1..Length(Exponents(t))] do
            l := l*Cgs(G)[j]^(Exponents(t))[j];
        od;

        #Update the conjugate
        h := h^l;

        #Calculate C_(G/Gi)(gGi)
        C := Kernel(f);

        #Calculate the preimage of C in G
        L := PreImage(nat, C);
    fi;
od;

#At this point we have C_G/Gn(gGn) but we need C_G(g) then we do another step:

fac := Length(efa)
N := SubgroupByIgs(G, NumeratorOfPcp( efa[fac] ));
gens := Igs( L );
imgs := List( gens , x -> Comm(h ,x) );
f := GroupHomomorphismByImages( L, N, gens, imgs );
C := Kernel(f);

e := [ Exponents( h^-1*g*Igs(f(L))[1] )[fac] ];
for j in [1..Length(imgs)] do
    Add(e, Exponents(imgs[j])[fac] );
od;
e := GcdRepresentation( e );

#find v such that f(v)=w
t := Identity( G );
for j in [2..Length(e)] do
    t := t*gens[j-1]^e[j];
od;

#find l such that nat(l)=v
l := Identity(G);
for j in [1..Length(Exponents(t))] do
    l := l*Cgs(G)[j]^(Exponents(t))[j];
od;