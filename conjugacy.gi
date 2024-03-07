####################################################################
## Global function to calculate the centralizer of g in G         ##
####################################################################

InstallGlobalFunction( "CentralizerNilGroup", function(G,g)

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

####################################################################
## Local function to the preimage of the conjugating element      ##
####################################################################

PreImageConjugate := function(G, img, nat, gens, imgs, i)

    local   e,
            vAi,
            v;

    #We have U={g1,...,gn} and f(gi)=g^ei then we want
    #Find e=(e1,...,en) of f(U)={g^e1,...,g^en}
    e := [];
    for j in [1..Length(imgs)] do
        Add(e, Exponents(imgs[j])[i]);
    od;

    #As f(U)=g^(gcd(e1,...,en)) we have f(g-1*k)=g^n with gcd(e1,...,en)|n
    #we can find a vector x1,...,xn such that f(g1^x1...gn^xn)=g^n
    e := GcdRepresentation(e)*( Exponents( Image( nat, img )  )[i] / Gcd( e ) );

    #Using the vector x1,...,xn we find uAi=g1^x1...gn^xnAi
    vAi := Identity( Image(nat) );
    for j in [1..Length(e)] do
        vAi := vAi*gens[j]^e[j];
    od;

    #Now we find v such that nat(v)= vAi            
    v := Identity(G);
    for j in [1..Length(Exponents(vAi))] do
        v := v*Cgs(G)[j]^(Exponents(vAi))[j];
    od;
    
    return v^-1;

end;

####################################################################
## Global function to check if g and h are conjugate in G         ##
####################################################################

InstallGlobalFunction( "IsConjugateNilGroup", function(G,g,h)

    local   efa,    #efa series of G
            U,      #Inverse image of the centralizer
            v,      #Conjugating on each step
            vAi,    #Image of the conjugating element
            k,      #Conjugate on each step
            exps,   #Storage of the exponents
            nat,    #Natural homomorphism G->G/Gi
            fac,    #Factor G(i-1)/Gi
            gens,   #Generators of U/Ai
            imgs,   #Generators of f(U/Ai)
            f;      #Group homomorphism

    efa  := EfaSeries(G);
    U    := G;
    exps := [];
    k    := h;

    #Catch trivial case
    nat := NaturalHomomorphismByNormalSubgroup(G, efa[2]);
    if Image(nat, k) <> Image(nat, h) then
        return rec( conj := Identity(G), state := false );
    fi;

    for i in [2..Length(efa)-1] do

        #Take G/Gi
        nat := NaturalHomomorphismByNormalSubgroup(G, efa[i+1]);

        #Take G(i-1)/Gi
        fac := Image( nat, efa[i] );

        #Create the homomorpism f:U/Gi -> G(i-1)/Gi
        gens:= Igs( Image( nat, U ) );
        imgs := List( gens , x -> Comm( Image(nat, g) ,x) );
        f   := GroupHomomorphismByImages( Image( nat, U ), fac, gens, imgs);
        
        #Check if (g^-1*k)Ai in f(U/Ai)
        if Image( nat, g^-1*k )  = Identity( Image(nat) ) then
        
            #In this case we don't have to do anything but as gives problems with 0/0 we have to diferenciate

        elif Image( nat, g^-1*k ) in Image(f) then
            
            v := PreImageConjugate(G, g^-1*k, nat, gens, imgs, i);
            Add(exps, v);
            k := k^v;

        else
            return rec( conj := Identity(G), state := false );
        fi;

        U   := PreImage(nat, CentralizerNilGroup( Image(nat), Image(nat, g) ) );
    od;

    v := Identity(G);
    for i in exps do
        v := v*i;
    od;

    #Check the result is correct
    if not h^v=g then
        Error(" The function is not working properly. ");
    fi;

    #Return the result
    return rec( conj := v, state := true );

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