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
## Helper function of CanonicalConjugateNilGroups, solves the     ##
## the problem of finding m minimal with the well-defined order   ##
####################################################################

MinimalSol := function(N, matrix, exp)

    local   l,
            a,
            b,
            n;

    #Flatten the matrix
    l := [];
    for i in [1..Length(matrix)] do
        Add(l,matrix[i][1]);
    od;

    #Coeficients of the problem
    a := Gcd(l);
    b := exp[1];
    n := AbsoluteValue(Int(b/a));

    #Minimize the solution under <<
    if a+b < 0 then
        a := a*( n + 1 );
    elif a+b > 0 then
        if 0 < -a*( n ) + b and -a*( n ) + b < a+b then 
            a := -a*( n );
        fi;
    fi;

    return rec( solv := Cgs(N)[1]^( a+b ), exp:= a+b );

end;

######################################################################
## Local function to calculate the canonical conjugate of elms in G ##
######################################################################

CanonicalConjugateNilGroupSeries := function(G, elms, pcps )

    local   C,      #Centralizer of elms in G
            h,      #CanonicalConjugate
            k,      #Conjugate element
            i,j,elm,#Bucle variable
            pcp,    #Factor on each step Gi/G(i+1)
            c,      #g^k in each step
            N,      #Subgroup Gi
            fac,    #Factor group C/Gi in each step
            gens,   #Generators of gen
            matrix, #Matrix representing the image of the homomorphism f
            exp,    #Exponents of c^-1*g in each pcp
            stb,    #Elements corresponding to the kernel and the preimages
            solv,   #Conjugating element in each step
            m,      #Minimal element to have h*m conjugate
            null;   #Kernel of f
            

    # the first layer 
    C := [];
    h := [];
    k := [];
    for elm in elms do
        Add( C, G);
        Add( h, G.1^( ExponentsByPcp(pcps[1], elm)[1] ) );
        Add( k, One(G) );
    od;

    for i in [2..Length(pcps)] do

        pcp  := pcps[i];
        N    := SubgroupByIgs( G, NumeratorOfPcp(pcp) );

        for j in [1..Length(elms)] do
        
            fac  := Pcp(C[j], N);
            gens := AsList(fac);
            c := h[j]^k[j];

            exp    := ExponentsByPcp( pcp, c^-1 * elms[j] );
            matrix := List( gens, x -> ExponentsByPcp( pcp, Comm(x,c) ) );
            Append( matrix, ExponentRelationMatrix( pcp ) );
        
            if matrix = 0*matrix then 
                #This case is when f is the identity homomorphism.
                m    := Cgs(N)[1]^exp[1];
                h[j] := h[j] * m;
                k[j] := k[j] * One(G);

            else
                # get solution if necessary
                solv := PcpSolutionIntMat( matrix, -exp );
                if IsBool( solv ) then 
                    
                    m    := MinimalSol( N, matrix, exp ).solv;
                    h[j] := h[j] * m;
                    exp  := ExponentsByPcp( pcp, (c*m)^-1 * elms[j] );
                    solv := PcpSolutionIntMat( matrix, -exp );

                fi;
                solv := solv{[1..Length(gens)]};

                # get nullspace
                null := PcpNullspaceIntMat( matrix );
                null := null{[1..Length(null)]}{[1..Length(gens)]};

                # calculate elements
                solv := MappedVector( solv, gens );
                gens := List( null, x -> MappedVector( x, gens ) );
                gens := Filtered( gens, x -> x <> x^0 );
                stb  := rec( stab := gens, prei := solv );

                # extract results
                k[j] := k[j] * stb.prei;
                stb  := AddIgsToIgs( stb.stab, Igs(N) );
                C[j] := SubgroupByIgs( G, stb );
            fi;
        od;
    od;

    return rec(conj := k, kano := h, cent := C);

end;

#######################################################################
## Global function to calculate the canonical conjugate of elms in G ##
#######################################################################

InstallGlobalFunction( "CanonicalConjugateNilGroup", function(G, elms)

    if not IsList(elms) then
        elms := [elms];
    fi;

    return CanonicalConjugateNilGroupSeries(G, elms, PcpsOfEfaSeries(G));

end );

######################################################################
## Local function to calculate solve the conjugacy problem in G     ##
## using canonical conjugate elements                               ##
######################################################################

IsCanonicalConjugateNilGroupSeries := function(G, elms, pcps )

    local   C,      #Centralizer of elms in G
            h,      #CanonicalConjugate
            k,      #Conjugate element
            i,j,elm,#Bucle variable
            pcp,    #Factor on each step Gi/G(i+1)
            c,      #g^k in each step
            N,      #Subgroup Gi
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
    C   := [];
    k   := [];
    ref := Exponents(elms[1])[1];
    h   := G.1^( ref );
    for elm in elms do
        if Exponents(elm)[1] <> ref then
            return false;
        fi;
        Add( C, G);
        Add( k, One(G) );
    od;

    for i in [2..Length(pcps)] do

        pcp  := pcps[i];
        N    := SubgroupByIgs( G, NumeratorOfPcp(pcp) );

        for j in [1..Length(elms)] do
            Print(j,"\n");
            fac  := Pcp(C[j], N);
            gens := AsList(fac);
            c := h^k[j];

            exp    := ExponentsByPcp( pcp, c^-1 * elms[j] );
            matrix := List( gens, x -> ExponentsByPcp( pcp, Comm(x,c) ) );
            Append( matrix, ExponentRelationMatrix( pcp ) );
        
            if matrix = 0*matrix then 
                #This case is when f is the identity homomorphism.
                exps := exp[1];
                m    := rec( solv := Cgs(N)[1]^exps);
                k[j] := k[j] * One(G);

            else
                # get solution if necessary
                solv := PcpSolutionIntMat( matrix, -exp );
                if IsBool( solv ) then 

                    m    := MinimalSol( N, matrix, exp );
                    exps := m.exp;
                    exp  := ExponentsByPcp( pcp, (c*m.solv)^-1 * elms[j] );
                    solv := PcpSolutionIntMat( matrix, -exp );

                else
                    exps := 0;
                fi;

                solv := solv{[1..Length(gens)]};

                # get nullspace
                null := PcpNullspaceIntMat( matrix );
                null := null{[1..Length(null)]}{[1..Length(gens)]};

                # calculate elements
                solv := MappedVector( solv, gens );
                gens := List( null, x -> MappedVector( x, gens ) );
                gens := Filtered( gens, x -> x <> x^0 );
                stb  := rec( stab := gens, prei := solv );

                # extract results
                k[j] := k[j] * stb.prei;
                stb  := AddIgsToIgs( stb.stab, Igs(N) );
                C[j] := SubgroupByIgs( G, stb );
            fi;

            if j = 1 then
                ref := exps;
                h   := h * m.solv;
            elif exps<>ref then
                return false;
            fi;
            
        od;
    od;

    c := k[1];
    for i in [1..Length(k)] do
        k[i] := k[i]^-1 * c;
    od;

    return rec(kano := h, conj :=k );

end;

######################################################################
## Global function to calculate solve the conjugacy problem in G    ##
## using canonical conjugate elements                               ##
######################################################################

InstallGlobalFunction( "IsCanonicalConjugateNilGroup", function(G, elms)

    if not IsList(elms) then
        elms := [elms];
    fi;

    return IsCanonicalConjugateNilGroupSeries(G, elms, PcpsOfEfaSeries(G));

end );

#######################################################################
## Local function to calculate the normalizer of U in G              ##
#######################################################################

NormalizerNilGroupSeries := function(G, U, efa)

    local   N,      #Normalizer of U in G
            H,      #Intersection of the previous step
            Hi,     #Intersection on each step
            i,      #Bucle variable
            h,      #Generator of U/V
            nat;    #Natural homomorphism N->N/V

    N    := G;
    H    := IntersectionEfaTerm( U, efa[Length(efa)-1] );

    for i in Reversed([1..Length(efa)-2]) do
    
        Hi  := IntersectionEfaTerm(U, efa[i]);

        if Hi <> H then
            h   := AsList( Pcp(Hi, H) );
            nat := NaturalHomomorphismByNormalSubgroup(N, H);
            N   := PreImage( nat, CentralizerNilGroup( Image(nat), nat(h) ) );
            H   := Hi;

        fi;

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
    Append( matrix, ExponentRelationMatrix( pcp ) );

    exp    := ExponentsByPcp( pcp, elm );

    solv   := PcpSolutionIntMat( matrix, exp );
    solv   := solv{[1..Length( gens )]};

    return MappedVector( solv, gens );

end;

#######################################################################
## Local function to calculate is two subgroups of G are conjugated  ##
#######################################################################

IsConjugateSubgroupsNilGroupSeries := function(G, U, V, efa)
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
            kan,    #Canonical conjugates of h and k
            xi;     #Conjugating element in each step

    x  := One(G);
    Ui := U;
    H  := IntersectionEfaTerm( Ui , efa[Length(efa)-1]);
    K  := IntersectionEfaTerm( V , efa[Length(efa)-1]);
    N  := G; 

    if H <> K then
        return false;
    fi;

    for i in Reversed([1..Length(efa)-2]) do
        Hi := IntersectionEfaTerm( Ui, efa[i]);
        Ki := IntersectionEfaTerm( V , efa[i]);
        
        if Hi <> H and Ki <> K then

            #Get the generators of U/W and V/W
            h   := AsList( Pcp(Hi, H) )[1];
            k   := AsList( Pcp(Ki, K) )[1];

            #Define the homomorphism N-> N/W
            nat := NaturalHomomorphismByNormalSubgroup(N, H );
            kan := IsCanonicalConjugateNilGroup( nat(N), [nat(k), nat(h)]);

            if IsBool(kan) then
                return false;
            else
                xi  := kan.conj[2];
                xi  := PreImageByQuotient( N, nat, xi );
                x   := x*xi;
                H   := Hi^xi;
                Ui  := Ui^xi;
                K   := Ki;
                N   := PreImage( nat, CentralizerNilGroup( Image(nat), nat(k) ) );
            fi;
        
        elif Hi = H and Ki = K then 
            #Process next
        else
            return false;
        fi;
    od;

    return x;

end;

#######################################################################
## Global function to calculate is two subgroups of G are conjugated ##
#######################################################################

InstallGlobalFunction( "IsConjugateSubgroupsNilGroup", function(G, U, V)

    if not (IsSubgroup(G,U) and IsSubgroup(G,V) ) then
        Error( "U and V have to be subgroups of G.");
    fi;

    if U = V then
        return One(G);
    else
        return IsConjugateSubgroupsNilGroupSeries(G, U, V, EfaSeries(G));
    fi;

end );    

CentralizerNilSubgroupGroupSeries := function(G, U, elms, pcps)

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

    pcp := pcps[1];
    N   := SubgroupByIgs( G, DenominatorOfPcp(pcp) );

    fac := Pcp(U, N); 
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

CentralizerNilSubgroupGroup := function(G, U, elms)

    if not IsList(elms) then
        elms := [elms];
    fi;

    if not IsNormal(G,U) then
        Error(" U has to be normal in G.");
    elif elms in U then
        return CentralizerNilGroup(U, elms);
    else
        return CentralizerNilSubgroupGroupSeries(G, U, elms, PcpsOfInducedEfaSeries(G, U) );
    fi;

    

end;