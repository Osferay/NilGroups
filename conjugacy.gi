###########################################################################
## Local function to calculate the centralizer of a set of elements in G ##
###########################################################################

CentralizerNilGroupSeries := function( G, elms, pcps )
    local   C,      #Centralizer of elms in G
            i,      #Bucle variable
            pcp,    #Factor on each step Gi/G(i+1)
            N,      #Subgroup Gi
            fac,    #Factor group C/Gi in each step
            gens,   #Generators of gen
            rels,   #Relation matrix of Gi/G(i+1)
            elm,    #Single elements on elms
            matrix, #Matrix representing the image of the homomorphism f
            null,   #Kernel of f
            stb;    #Elements corresponding to the kernel

    C := G;

    for i in [2..Length(pcps)] do

        pcp := pcps[i]; 
        N   := SubgroupByIgs( G, NumeratorOfPcp(pcp) );

        fac := Pcp(C, N); 
        gens:= AsList(fac);

        rels := ExponentRelationMatrix( pcp );
        stb  := gens;
        for elm in elms do
            if Length( stb ) <> 0 then 

                # set up matrix
                matrix := List( stb, h -> ExponentsByPcp( pcp, Comm(h,elm) ) );
                Append( matrix, rels );

                # get nullspace
                null := PcpNullspaceIntMat( matrix );
                null := null{[1..Length(null)]}{[1..Length(stb)]};
                # calculate elements corresponding to null
                stb  := List( null, x -> MappedVector( x, stb ) );
                stb  := Filtered( stb, x -> x <> x^0 );
            fi;
        od;

        stb := AddIgsToIgs( stb, Igs(N) );
        C   := SubgroupByIgs( G, stb );

    od;

    return(C);

end;

############################################################################
## Global function to calculate the centralizer of a set of elements in G ##
############################################################################

InstallGlobalFunction( "CentralizerNilGroup", function(G, elms) 

    if not IsList(elms) then
        elms := [elms];
    fi;

    return CentralizerNilGroupSeries(G, elms, PcpsOfEfaSeries(G));

end );

############################################################################
## Local function to check if two elements in G are conjugate             ##
############################################################################

IsConjugateNilGroupSeries := function(G, g, h, pcps )

    local   C,      #Centralizer of elms in G
            k,      #Conjugate element
            i,      #Bucle variable
            pcp,    #Factor on each step Gi/G(i+1)
            c,      #g^k in each step
            N,      #Subgroup Gi
            fac,    #Factor group C/Gi in each step
            gens,   #Generators of gen
            matrix, #Matrix representing the image of the homomorphism f
            exp,    #Exponents of c^-1*g in each pcp
            stb,    #Elements corresponding to the kernel and the preimages
            solv,   #Conjugating element in each step
            null;   #Kernel of f
            

    # the first layer
    if ExponentsByPcp(pcps[1], g) <> ExponentsByPcp(pcps[1], h) then 
        return false; 
    fi;

    C := G;
    k := One(G);

    for i in [2..Length(pcps)] do

        pcp := pcps[i];
        c   := g^k;

        N   := SubgroupByIgs( G, NumeratorOfPcp(pcp) );
        fac := Pcp(C, N);
        gens:= AsList(fac);
        exp := ExponentsByPcp( pcp, c^-1 * h );
        if Length(gens) = 0 then
            if exp = 0*exp then
                stb := rec( stab := gens, prei := c^0 );
            else
                return false;
            fi;

        else

            # set up matrix
            matrix := List( gens, x -> ExponentsByPcp( pcp, Comm(x,c) ) );
            Append( matrix, ExponentRelationMatrix( pcp ) );

            # get solution
            solv := PcpSolutionIntMat( matrix, -exp );
            if IsBool( solv ) then 
                return false; 
            else
                solv := solv{[1..Length(gens)]};

                # get nullspace
                null := PcpNullspaceIntMat( matrix );
                null := null{[1..Length(null)]}{[1..Length(gens)]};

                # calculate elements
                solv := MappedVector( solv, gens );
                gens := List( null, x -> MappedVector( x, gens ) );
                gens := Filtered( gens, x -> x <> x^0 );
                stb  := rec( stab := gens, prei := solv );
            fi;
        fi;

        # extract results
        k   := k * stb.prei;
        stb := AddIgsToIgs( stb.stab, Igs(N) );
        C   := SubgroupByIgs( G, stb );
    
    od;

    return k;

end;

####################################################################
## Global function to check if g and h are conjugate in G         ##
####################################################################

InstallGlobalFunction( "IsConjugateNilGroup", function(G,g,h)

    return IsConjugateNilGroupSeries(G,g,h,PcpsOfEfaSeries(G));

end );

####################################################################
## Helper functions of CanonicalConjugateNilGroups, solves the    ##
## the problem of finding m minimal with the well-defined order   ##
####################################################################

####################################################################
## Minimize the function a + x*coef                               ##
####################################################################

MinEq := function( a, coef)

    local   d,      #Gcd of coeficients
            m,      #Minimum of coeficients
            rep;    #Representation of the minimum

    d := Gcd(coef);
    m := a mod d;

    if m = a then 

        return rec( min := a, rep := coef*0);

    else

        rep := (m - a)/d;
        rep := rep * GcdRepresentation( coef );
        return rec( min := m, rep := rep);

    fi;
end;

####################################################################
## Minimize the number of relations used to speed up computations ##
####################################################################

MinimizeRelations := function(l, exp, o)

    local   min,    #Minimum number of relations
            i,      #Bucle variable
            d,      #Gcd of the minimum number of relations
            n,      #Gcd of d and the order
            z;      #Solution
            
    min := [];
    for i in [1..Length(l)] do 
        Add(min, l[i]);
        d := Gcd( min );
        n := Gcd( d , o);
        if exp mod n = 0 then
            z := GcdRepresentation( d/n, o/n );
            z := (exp/n * z[1]) mod (o/n);
            return ( z* GcdRepresentation( min ) mod o);
        fi;
    od;
end;

####################################################################
## Minimize the number of relations used to speed up computations ##
####################################################################

MinimalSolFFE := function(pcp, matrix, gens, h, c, o)

    local   gen,    #Generator of the pcp
            d,      #Depth of the pcp
            pos,    #Positions where the matrix is non zero
            a,      #List of commutators
            exp,    #Exponents
            b,      #Minimum of the equation
            rep,    #Minimum relations
            n,      #Integer variable
            i;      #Bucle variable
    
    gen := pcp[1]; 
    d   := Depth( gen ); 
    exp := Exponents( h^-1 * c )[d];

    #Order the generators to speed computations
    pos := PositionsProperty( matrix , x -> x <> 0*x ); 
    pos := Reversed(gens{pos});

    #Find the exponents of each ordered generator
    a   := List( pos, x -> Comm(c,x) );
    d   := List( a,   x -> Exponents( x )[d] );
    if ForAll( d, x -> x = d[1] ) then
        rep := d*0;
        rep[1] := 1;
    else
        rep := GcdRepresentation( d );
    fi;

    #Minimize the equation
    n   := d * rep;
    a   := Gcd( n, o ); 
    b   := exp mod a; 

    rep := MinimizeRelations( d, b - exp, o);
    pos := pos{[1..Length(rep)]};

    exp := pos[1]^0;
    for i in [1..Length(rep)] do
        rep[i] := rep[i] mod ( o/Gcd(d[i],o) );
        a := (-rep[i] mod o);
        if a < rep[i] then
            exp := exp * pos[i]^(-a);
        else
            exp := exp * pos[i]^( rep[i] );
        fi;
    od;
    
    return rec( solv := gen^( b ), exp := exp );

end;

####################################################################
## Helper function of CanonicalConjugateNilGroups, checks if      ##
## the an element is conjugated in a finite group                 ##
####################################################################

PcpSolutionFFEMat := function( matrix, exp, o)

    local   l,      #Flattened matrix
            i,      #Bucle variable
            n,      #Gcd of the flattened matrix and 
            pos,    #Gcd representation of the matrix
            solv;   #Solution

    #Check if the equation has solution
    l      := List( matrix, x -> x[1] );
    l      := Reversed(l);
    solv   := Gcd( l );
    pos    := GcdRepresentation( l );
    n      := Gcd( solv, o );
    if (exp mod n) <> 0 then return false; fi;

    #Now solve it , trying to use the minimum number of relations
    solv   := MinimizeRelations(l, exp, o);
    
    for i in [1..Length(solv)] do
        if solv[i] <> 0 then
            solv[i] := solv[i] mod (o/Gcd(l[i],o));
            pos := (-solv[i] mod o);
            if pos < solv[i] then
                solv[i] := - pos;
            fi;
        fi;
    od;
    
    return Reversed(solv);

end;

######################################################################
## Local function to calculate the canonical conjugate of elms in G ##
######################################################################

CanonicalConjugateNilGroupSeries := function(G, elms, pcps )

    local   igs,    #Igs of G.
            C,      #Centralizer of kanos in G
            h,      #CanonicalConjugate
            k,      #Conjugate element
            i,j,elm,#Bucle variable
            pcp,    #Factor on each step Gi/G(i+1)
            o,      #Factor order
            c,      #g^k in each step
            N,      #Subgroup Gi
            fac,    #Factor group C/Gi in each step
            gens,   #Generators of gen
            matrix, #Matrix representing the image of the homomorphism f
            rels,   #Relations of pcp
            exp,    #Exponents of c^-1*g in each pcp
            stb,    #Elements corresponding to the kernel and the preimages
            solv,   #Conjugating element in each step
            m,      #Minimal element to have h*m conjugate
            null;   #Kernel of f
            

    # the first layer 
    igs := Igs(G);
    C   := [];
    h   := [];
    k   := [];
    for elm in elms do
        Add( C, G);
        Add( h, G.1^( ExponentsByIgs( igs, elms[1])[1] ) );
        Add( k, One(G) );
    od;

    for i in [2..Length(pcps)] do

        pcp  := pcps[i];
        o    := FactorOrder( pcp[1] );
        N    := SubgroupByIgs( G, NumeratorOfPcp(pcp) );

        for j in [1..Length(elms)] do
        
            fac    := Pcp(C[j], N);
            gens   := AsList(fac);
            c      := elms[j]^k[j];
            matrix := List( gens, x -> [ ExponentsByIgs( igs, Comm(x,c) )[i] ] );
            exp := ExponentsByIgs( igs, h[j]^-1 * c )[i];

            if matrix = 0*matrix then 
                #This case is when f is the identity homomorphism.
                
                m    := pcp[1]^exp;
                h[j] := h[j] * m;
                Append(matrix, ExponentRelationMatrix( pcp ));

            elif o = 0 then
                # check if they are conjugated
                solv := PcpSolutionIntMat( matrix, [exp] );
                
                if IsBool( solv ) then 
                    #if not add the minimal element to be conjugated
                    m      := List( matrix, x -> x[1] );
                    m      := MinEq( exp , m ).min;
                    h[j]   := h[j] * pcp[1]^m;
                    
                    #Calculate the new exponent
                    exp    := exp - m;
                    solv   := PcpSolutionIntMat( matrix, [exp] );
                fi;

                #Calculate the conjugating element
                solv := MappedVector( solv, gens );
                k[j] := k[j] * solv;
            
            else
                solv := PcpSolutionFFEMat(matrix, exp, o);
                if IsBool(solv) then
                    m    := MinimalSolFFE(pcp, matrix, gens, h[j], c, o);
                    h[j] := h[j] * m.solv;
                    k[j] := k[j] * m.exp;
                else
                    solv := MappedVector( solv, Reversed( Reversed(gens){[1..Length(solv)]} ) );
                    k[j] := k[j]*solv;
                fi; 
                Append(matrix, ExponentRelationMatrix( pcp ));
            fi;
            # get the kernel
            null := PcpNullspaceIntMat( matrix );
            null := null{[1..Length(null)]}{[1..Length(gens)]};

            # calculate elements
            gens := List( null, x -> MappedVector( x, gens ) );
            gens := Filtered( gens, x -> x <> x^0 );   
            stb  := AddIgsToIgs( gens, NumeratorOfPcp(pcp) );
            C[j] := SubgroupByIgs( G, stb );
        od;
    od;

    return rec(conj := k, kano := h, cent := C);

end;

#######################################################################
## Global function to calculate the canonical conjugate of elms in G ##
#######################################################################

InstallGlobalFunction( "CanonicalConjugateNilGroup", function(G, elms)

    local   pcps;

    if not IsList(elms) then
        elms := [elms];
    fi;
    pcps := PcpsBySeries( PcpSeries(G) );
    return CanonicalConjugateNilGroupSeries(G, elms, pcps);
    # return CanonicalConjugateNilGroupSeries(G, elms, PcpsOfEfaSeries(G));

end );

######################################################################
## Local function to calculate solve the conjugacy problem in G     ##
## using canonical conjugate elements                               ##
######################################################################

IsCanonicalConjugateNilGroupSeries := function(G, elms, pcps )

    local   igs,    #Igs of G
            C,      #Centralizer of elms in G
            h,      #CanonicalConjugate
            k,      #Conjugate element
            i,j,elm,#Bucle variable
            pcp,    #Factor on each step Gi/G(i+1)
            c,      #g^k in each step
            N,      #Subgroup Gi
            o,      #Factor order
            fac,    #Factor group C/Gi in each step
            gens,   #Generators of gen
            matrix, #Matrix representing the image of the homomorphism f
            exp,    #Exponents of c^-1*g in each pcp
            stb,    #Elements corresponding to the kernel and the preimages
            solv,   #Conjugating element in each step
            m,      #Minimal element to have h*m conjugate
            exps,   #Exponent of the minimal solution
            ref,    #Reference value of m
            null;   #Kernel of f
            

    # the first layer 
    igs := Igs(G);
    C   := G;
    k   := [One(G)];
    ref := ExponentsByIgs( igs, elms[1])[1];
    h   := G.1^( ref );

    for i in [2..Length(elms)] do
        if ExponentsByIgs( igs, elms[i] )[1] <> ref then
            return false;
        fi;
        k[i] := One(G);
    od;

    for i in [2..Length(pcps)] do

        pcp    := pcps[i];
        N      := SubgroupByIgs( G, NumeratorOfPcp(pcp) );
        o      := FactorOrder( pcp[1] );

        fac    := Pcp(C, N);
        gens   := AsList(fac);

        #For the fist element
        c      := elms[1]^k[1];
        matrix := List( gens, x -> [ ExponentsByIgs( igs, Comm(x,c) )[i] ] );
        exp := ExponentsByIgs( igs, h^-1 * c )[i];

        if matrix = 0*matrix then 
            #This case is when f is the identity homomorphism.
            h   := h * pcp[1]^exp;
            Append( matrix, ExponentRelationMatrix( pcp ) );

        elif o = 0 then
            # get solution if necessary
            solv := PcpSolutionIntMat( matrix, [exp] );

            if IsBool( solv ) then 

                m      := List( matrix, x -> x[1] );
                m      := MinEq( exp , m ).min;
                h      := h * pcp[1]^m;

                #Calculate the new exponent
                exp    := exp - m;
                solv   := PcpSolutionIntMat( matrix, [exp] );
                
            fi;

            # calculate elements
            solv := MappedVector( solv, gens );
            # extract results
            k[1] := k[1] * solv;

        else

            solv := PcpSolutionFFEMat(matrix, exp, o);

            if IsBool(solv) then
            
                m := MinimalSolFFE(pcp, matrix, gens, h, c, o);
                h := h * m.solv;
                k[1] := k[1] * m.exp;
            else
            
                solv := MappedVector( solv, Reversed( Reversed(gens){[1..Length(solv)]} ) );
                k[1] := k[1]*solv;
            fi; 
            
            Append(matrix, ExponentRelationMatrix( pcp ));
        fi;
        
        
        # Calculate the centralizer
        null := PcpNullspaceIntMat( matrix );
        null := null{[1..Length(null)]}{[1..Length(gens)]};

        stb  := List( null, x -> MappedVector( x, gens ) );
        stb  := Filtered( stb, x -> x <> x^0 );

        stb  := AddIgsToIgs( stb, NumeratorOfPcp(pcp) );
        C    := SubgroupByIgs( G, stb );

        #For the other elements
        for j in [2..Length(elms)] do

            c      := elms[j]^k[j];
            exp    := ExponentsByIgs( igs, h^-1 * c )[i];
            matrix := List( gens, x -> [ ExponentsByIgs( igs, Comm(x,c) )[i] ] );

            if matrix = 0*matrix and exp <> 0 then 
                return false;

            elif o = 0 then
                # get solution if necessary
                solv := PcpSolutionIntMat( matrix, [exp] );
                if IsBool(solv) then return false; fi;
    
                # calculate elements
                solv := MappedVector( solv, gens );
                # extract results
                k[j] := k[j] * solv;

            else
                solv := PcpSolutionFFEMat( matrix, exp, o);
                if IsBool(solv) then return false; fi;

                gens := Reversed( Reversed(gens){[1..Length(solv)]} );
                solv := MappedVector( solv, gens );
                k[j] := k[j] * solv;
            fi;

        od;

    od;

    return rec(kano := h, conj := k, cent := C );

end;

######################################################################
## Global function to calculate solve the conjugacy problem in G    ##
## using canonical conjugate elements                               ##
######################################################################

InstallGlobalFunction( "IsCanonicalConjugateNilGroup", function(G, elms)

    local   pcps,
            can;

    if not IsList(elms) then
        elms := [elms];
    fi;

    pcps := PcpsBySeries( PcpSeries(G) );
    can  := IsCanonicalConjugateNilGroupSeries(G, elms, pcps);

    return can;
    
end );

#######################################################################
## Local function to calculate the normalizer of U in G              ##
#######################################################################

NormalizerNilGroupSeries := function(G, U, efa)

    local   genU,   #Generators of U
            N,      #Normalizer of U in G
            H,      #Intersection of the previous step
            Hi,     #Intersection on each step
            i,      #Bucle variable
            h,      #Generator of U/V
            nat;    #Natural homomorphism N->N/V
    
    genU := Reversed( Igs(U) );
    N    := G;
    H  := Subgroup(U, [ genU[1] ]);

    for i in [2..Length(genU)] do
    
        Hi := Subgroup(U, genU{[1..i]});
        h   := Pcp(Hi, H)[1];
        nat := NaturalHomomorphismByNormalSubgroup(N, H);
        N   := PreImage( nat, CentralizerNilGroup( Image(nat), nat(h) ) );
        H   := Hi;

    od;

    return N;

end;

#######################################################################
## Global function to calculate the normalizer of U in G             ##
#######################################################################

InstallGlobalFunction( "NormalizerNilGroup", function(G,U)

    if not IsSubgroup(G,U) then
        Error( "U has to be a subgroup of G.");
    fi;

    return NormalizerNilGroupSeries(G, U, EfaSeries(G));

end );


#######################################################################
## Local function to calculate the preimage of an element            ##
#######################################################################

PreImageByQuotient := function(G, hom, elm )

    local   pcp,
            gens,
            matrix,
            exp,
            solv;

    pcp    := Pcp( Image(hom) );
    gens   := Cgs(G);
    matrix := List( gens, x-> ExponentsByPcp( pcp, hom(x) ) );
    exp    := ExponentsByPcp( pcp, elm );

    solv   := PcpSolutionIntMat( matrix, exp );
    return MappedVector( solv, gens );

end;

#######################################################################
## Local function to calculate is two subgroups of G are conjugated  ##
#######################################################################

IsConjugateSubgroupsNilGroup := function(G, U, V)
    local   x,      #Conjugating element of U and V
            Ui,     #Conjuate of U in each step
            H,      #Intersection of previous step of U
            K,      #Intersection of previous step of V
            N,      #Normalizer of W
            i,      #Bucle variable
            Hi,     #Intersection in each step of U
            Ki,     #Intersection in each step of V
            h,      #Generator of U/W
            k,      #Generator of V/W
            nat,    #Natural homomorphism N->N/W
            conj,   #Conjugating of h and k
            xi,     #Conjugating element in each step
            genU,   #Generators of U
            genV;   #Generators of V

    #Catch trivial case
    genU := Reversed( Igs(U) );
    genV := Reversed( Igs(V) );
    
    if Length(genU) <> Length(genV) then
        return false;
    fi;

    #Start the algorithm
    x  := One(G);
    Ui := U;
    H  := Subgroup(U, [  ]);
    K  := Subgroup(V, [  ]);
    N  := G; 

    if H <> K then
        return false;
    fi;

    for i in [1..Length( genU )] do
        #Take the intersection
        Hi := Subgroup(Ui, genU{[1..i]});
        Ki := Subgroup(V, genV{[1..i]});

        #Get the generators of U/W and V/W
        h   := Pcp(Hi, H)[1];
        k   := Pcp(Ki, K)[1];

        #Define the homomorphism N-> N/W
        nat := NaturalHomomorphismByNormalSubgroupNC(N, H );
        conj:= IsConjugateNilGroup( Image(nat), nat(h), nat(k) );

        if IsBool(conj) then
            return false;
        else
            #They are conjugated, update values
            xi  := PreImageByQuotient( N, nat, conj );
            x   := x * xi;
            H   := Hi^ xi;
            Ui  := Ui^ xi;
            genU:= Reversed( Igs(Ui) );
            K   := Ki;

            #Update the normalizer
            N   := PreImage( nat, CentralizerNilGroup( Image(nat), nat(k) ) );
        fi;
        
    od;

    return x;

end;

#######################################################################
## Global function to calculate is two subgroups of G are conjugated ##
#######################################################################

InstallGlobalFunction( "IsConjugateSubgroups", function(G, U, V)

    if not (IsSubgroup(G,U) and IsSubgroup(G,V) ) then
        Error( "U and V have to be subgroups of G.");
    fi;

    if U = V then
        return One(G);
    else
        return IsConjugateSubgroupsNilGroup(G, U, V);
    fi;

end );    

#######################################################################
## Local function to calculate the reduced canonical form            ##
#######################################################################

ReducedCanonical := function(G, U, elms)

    local   nat,    #Natural homomrphism from G to G/U
            kan,    #Canonical conjugate of the image of g by nat
            k,      #Reduced preimage of the canonical conjugate
            v,      #Conjugating element
            N,      #Normalizer
            i;      #Bucle variable

    nat := NaturalHomomorphismByNormalSubgroup(G, U );
    kan := CanonicalConjugateNilGroup( nat(G), List( elms, nat ));
    k   := [];
    v   := [];

    for i in [1..Length(elms)] do
        k[i]:= PreImagesRepresentative(nat, kan.kano[i]);
        v[i]:= PreImagesRepresentative(nat, kan.conj[i]);
        k[i]:= ReducedPreimage(k[i] , Cgs(U));
    od;

    N   := PreImage( nat, kan.cent[1] );

    return rec( kan := k, v := v, N := N );

end;

#######################################################################
## Local function to calculate the canonical conjugate subgroup of  ##
## a subgroup in G                                                   ##
#######################################################################

CanonicalConjugateSubgroupNilGroup := function(G, U)

    local   gU,     #Generators of U 
            N,      #Normalizer of K
            x,      #Conjugating element
            gK,     #Generators of K
            Ui,     #Subgroup Ui in the current step
            K,      #Canonical subgroup
            i,      #Bucle variable
            u,      #Generator on each step
            kan;    #Canonical conjugates of h

    gU := Reversed( Cgs(U) );
    N  := G; 
    x  := One(G);
    Ui := Subgroup( U, [ ]);
    gK := [];

    for i in [1..Length(gU)] do

        u := NormedPcpElement( gU[i]^x );
        kan := ReducedCanonical(N, Ui, u);

        x   := x*kan.v[1];
        Ui  := SubgroupByIgs(U, gU{[1..i]} );
        Ui   := Ui^x;
        
        N   := kan.N;

        Add( gK, kan.kan[1] );

    od;
    
    K    := Subgroup( G, gK );

    return rec( kano := K, conj := x, norm := N);

end;

#######################################################################
## Global function to calculate the canonical conjugate subgroup of  ##
## a subgroup in G                                                   ##
#######################################################################

InstallGlobalFunction( "CanonicalConjugateSubgroup", function(G, U)

    if not IsSubgroup(G,U) then
        Error( "U has to be subgroups of G.");
    fi;

    return CanonicalConjugateSubgroupNilGroup(G, U); 

end );  

######################################################################
## Local function to calculate solve the conjugacy problem in G     ##
## using canonical conjugate elements                               ##
######################################################################

IsCanonicalConjugateSubgroupNilGroup := function(G, U, V)

    local   gU,     #Generators of U 
            gV,     #Generators of V
            N,      #Normalizer of K
            x,      #Conjugating element of U
            y,      #Conjugating element of V
            gK,     #Generators of K
            Ui,     #Subgroup Ui in the current step
            K,      #Canonical subgroup
            i,      #Bucle variable
            u,      #Generator on each step of U
            v,      #Generator on each step of V
            kU;     #Canonical conjugates of 

    gU := Reversed( Cgs(U) );
    gV := Reversed( Cgs(V) );
    N  := G; 
    x  := One(G);
    y  := One(G);
    Ui := Subgroup( U, [ ]);
    gK := [];

    if Length(gU) <> Length(gV) then
        return false;
    fi;

    for i in [1..Length(gU)] do

        u := NormedPcpElement( gU[i]^x );
        v := NormedPcpElement( gV[i]^y );

        kU := ReducedCanonical(N, Ui, [u,v]);

        if kU.kan[1] = kU.kan[2] then
            x   := x*kU.v[1];
            Ui  := SubgroupByIgs(U, gU{[1..i]} );
            Ui  := Ui^x;

            y   := y*kU.v[2];
        
            N   := kU.N;

            Add( gK, kU.kan[1] );
        else
            return false;
        fi;

    od;
    
    K    := Subgroup( G, gK );

    return rec( kano := K, conj := [x,y], norm := N);

end;

######################################################################
## Global function to calculate solve the conjugacy problem in G    ##
## using canonical conjugate elements                               ##
######################################################################

InstallGlobalFunction( "IsCanonicalConjugateSubgroups", function(G, U, V)

    if not (IsSubgroup(G,U) and IsSubgroup(G,V) ) then
        Error( "U and V have to be subgroups of G.");
    fi;

    if U = V then
        return One(G);
    else
        return IsCanonicalConjugateSubgroupNilGroup(G, U, V);
    fi;

end );    

ConjugateList := function(G, list1)

    local   x,
            Ui,
            gi,
            list2,
            xi,
            i;

    x := One(G);
    list2 := [];

    for i in [ 1..Length(list1) ] do
        
        Ui := CentralizerNilGroup( G, list1{[1..i]});
        xi := Random(Ui);
        x := x * xi;
        Add( list2 , list1[i]^x );
    
    od;

    return rec( new := list2, conj := x) ;

end;

IsMultipleConjugateNilGroup := function(G, list1, list2) 
    
    local   x,
            Ui,
            gi,
            hi,
            xi,
            i;

    if Length(list1) <> Length(list2) then
        return false;
    fi;

    x := One(G);

    for i in [ 1..Length(list1) ] do
        
        Ui := CentralizerNilGroup( G, list2{[1..i]});
        gi := list1[i]^x;
        hi := list2[i];
        xi := IsCanonicalConjugateNilGroup(Ui, [gi, hi]);
        if IsBool(xi) then
            return false;
        else
            xi := xi.conj[1]*xi.conj[2]^-1;
            x := x * xi;
        fi;
    
    od;

    return rec( U := Ui, x := x );

end;