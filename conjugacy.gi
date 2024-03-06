InstallGlobalFunction( "CentralizerFGNil", function(G,g)

    local   efa,    #Elementary/Free abelian series of G
            C,      #Centralizer of G/Gi
            L,      #L/Gi = Ci
            N,      #Gi
            nat,    #Natural homomorphism G -> G/Gi
            M,      #G(i-1)
            gGi,    #Image g->nat(g)=gGi
            fac,    #factor G(i-1)/Gi
            LN,     #factor L/Gi
            gens,   #Generators L/Gi
            imgs,   #gens -> [g,gen]
            i,      #Bucle variables
            f;      #homomorphism L/Gi -> G(i-1)/Gi 

    efa := PcpsOfEfaSeries(G);
    C := G;
    L := G;

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

        else
        
            #Get G(i-1)
            M := SubgroupByIgs(G, NumeratorOfPcp( efa[i-1] ));

            #Calculate gGi
            gGi := Image(nat, g);

            #Calculate G(i-1)/Gi under G/Gi
            fac := Image( nat, M );

            #Calculate L/Gi under G/Gi
            LN := Image(nat, L);

            #Calculate f: L/Gi -> G(i-1)/Gi
            gens := Igs( LN );
            imgs := List( gens , x -> Comm(gGi ,x) );
            f := GroupHomomorphismByImages( LN , fac, gens, imgs );

            #Calculate C_(G/Gi)(gGi)
            C := Kernel(f);

            #Calculate the preimage of C in G
            L := PreImage(nat, C);
        fi;
    od;

    #At this point we have C_G/Gn(gGn) but we need C_G(g) then we do another step:

    N := SubgroupByIgs(G, NumeratorOfPcp( efa[Length(efa)] ));
    gens := Igs( L );
    imgs := List( gens , x -> Comm(g ,x) );
    f := GroupHomomorphismByImages( L, N, gens, imgs );
    C := Kernel(f);

    return C;

end );


InstallGlobalFunction( "CanonicalConjugate", function(G, g)

    local   efa,    #Elementary/Free abelian series of G
            C,      #Centralizer of G/Gi
            L,      #L/Gi = Ci
            h,      #Conjugate element
            N,      #Gi
            nat,    #Natural homomorphism G -> G/Gi
            M,      #G(i-1)
            hGi,    #Image g->nat(g)=gGi
            fac,    #factor G(i-1)/Gi
            LN,     #factor L/Gi
            gens,   #Generators L/Gi
            imgs,   #gens -> [g,gen]
            f,      #homomorphism L/Gi -> G(i-1)/Gi 
            e,      #exponent vector of f(L)
            t,      #element in L/Gi such that f(t)=generator(f(L))
            i,j,  #Bucle variable
            l;      #element in L such that nat(l)=t

    efa := PcpsOfEfaSeries(G);
    C := G;
    L := G;
    h := g;

    for i in [2..Length(efa)] do
        #We have efa[i]=Gi/G(i+1)
        #Get Gi
        N := SubgroupByIgs(G , NumeratorOfPcp( efa[i] ));
        #Calculate G -> G/Gi
        nat := NaturalHomomorphismByNormalSubgroup( G, N );       
        
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
        e := [];
        for j in [1..Length(imgs)] do
            Add(e, Exponents(imgs[j])[i-1]);
        od;
        e := GcdRepresentation(e);

        #find v such that f(v)=w
        t := Identity( Image(nat) );
        for j in [1..Length(e)] do
            t := t*gens[j]^e[j];
        od;

        #find l such that nat(l)=v
        l := Identity(G);
        for j in [1..Length(Exponents(t))] do
            l := l*Cgs(G)[j]^(Exponents(t))[j];
        od;

        #Calculate C_(G/Gi)(hGi)
        C := Kernel(f);

        #Calculate the preimage of C in G
        L := PreImage(nat, C);
        
        #Update the conjugate
        h := h^l;
    od;

    #At this point we have hGn~gGn but we need h~g then we do another step:

    N := SubgroupByIgs(G, NumeratorOfPcp( efa[Length(efa)] ));
    gens := Igs( L );
    imgs := List( gens , x -> Comm(g ,x) );
    f := GroupHomomorphismByImages( L, N, gens, imgs );
    C := Kernel(f);

    #Solve the equation of exponents
    e := [];
    for j in [1..Length(imgs)] do
        Add(e, Exponents(imgs[j])[Length(efa)]);
    od;
    e := GcdRepresentation(e);

    #find v such that f(v)=w
    t := Identity( G );
    for j in [1..Length(e)] do
        t := t*gens[j]^e[j];
    od;

    #find l such that nat(l)=v
    l := Identity(G);
    for j in [1..Length(Exponents(t))] do
        l := l*Cgs(G)[j]^(Exponents(t))[j];
    od;
    l := l^-1;

    return rec( cent:=C, conj := h, exp:=l );

end );