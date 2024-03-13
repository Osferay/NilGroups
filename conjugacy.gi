###########################################################################
## Local function to calculate the centralizer of a set of elements in G ##
###########################################################################

CentralizerNilGroupSeries := function( G, elms, pcps )
    local   C,
            i, 
            pcp, 
            N, 
            gen, 
            gens, 
            rels, 
            elm, 
            matrix, 
            null,
            stb;

    C := G;

    for i in [2..Length(pcps)] do

        pcp := pcps[i]; 
        N   := SubgroupByIgs( G, NumeratorOfPcp(pcp) );

        gen := Pcp(C, N); 
        gens:= AsList(gen);

        rels := ExponentRelationMatrix( pcp );
        for elm in elms do
            if Length( gens ) <> 0 then 

                # set up matrix
                matrix := List( gens, h -> ExponentsByPcp( pcp, Comm(h,elm) ) );
                Append( matrix, rels );

                # get nullspace
                null := PcpNullspaceIntMat( matrix );
                null := null{[1..Length(null)]}{[1..Length(gens)]};

                # calculate elements corresponding to null
                stb  := List( null, x -> MappedVector( x, gens ) );
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

    local C;

    if not IsList(elms) then
        elms := [elms];
    fi;

    C := CentralizerNilGroupSeries(G, elms, PcpsOfEfaSeries(G));

    return C;

end );

############################################################################
## Local function to check if two elements in G are conjugate             ##
############################################################################

IsConjugateNilGroupSeries := function(G, g, h, pcps )

    local   C,
            k,
            i,
            pcp,
            c,
            N,
            fac,
            gens,
            matrix,
            exp,
            stb,
            solv,
            null;
            

    # the first layer
    if ExponentsByPcp(pcps[1], g) <> ExponentsByPcp(pcps[1], h) then 
        return false; 
    fi;

    C := G;
    k := One(G);

    for i in [2..Length(pcps)] do

        pcp := pcps[i];
        c := g^k;

        N := SubgroupByIgs( G, NumeratorOfPcp(pcp) );
        fac := Pcp(C, N);
        gens := AsList(fac);
        exp := ExponentsByPcp( pcp, c^-1 * h );
        Print(exp, "\n");
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
        k := k * stb.prei;
        stb := AddIgsToIgs( stb.stab, Igs(N) );
        C := SubgroupByIgs( G, stb );
    
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
## Local function to calculate the canonical conjugate of g in G  ##
####################################################################

CanonicalConjugateNilGroupSeries := function(G, g, pcps )

    local   C,
            k,
            i,
            pcp,
            c,
            N,
            fac,
            gens,
            matrix,
            exp,
            stb,
            solv,
            l,
            a,
            b,
            al,
            null;
            

    # the first layer 
    C := G;
    h := G.1^( ExponentsByPcp(pcps[1], g)[1] );
    k := One(G);

    for i in [2..Length(pcps)] do

        pcp := pcps[i];
        c := h^k;

        N := SubgroupByIgs( G, NumeratorOfPcp(pcp) );
        fac := Pcp(C, N);
        gens := AsList(fac);

        exp := ExponentsByPcp( pcp, c^-1 * g );
        matrix := List( gens, x -> ExponentsByPcp( pcp, Comm(x,c) ) );
        Append( matrix, ExponentRelationMatrix( pcp ) );
        
        if matrix = 0*matrix then 

            #This case is when f is the identity homomorphism.
            h := h * Cgs(N)[1]^exp[1];
            k := k * One(G);

        else
            # get solution
            solv := PcpSolutionIntMat( matrix, -exp );
            if IsBool( solv ) then 
                l := [];
                # we have that they are non conjugate so we have to work
                for i in [1..Length(matrix)] do
                    Add(l,matrix[i][1]);
                od;
                a := Gcd(l);
                b := exp[1];

                if a+b < 0 then
                    a := a*( AbsoluteValue(Int(b/a)) + 1 );
                elif a+b > 0 then
                    if 0 < -a*( AbsoluteValue(Int(b/a)) ) + b and -a*( AbsoluteValue(Int(b/a)) ) < a+b then 
                        a := -a*( AbsoluteValue(Int(b/a)) );
                    fi;
                fi;

                al := Cgs(N)[1]^( a+b );

                h := h * al;
                exp := ExponentsByPcp( pcp, (c*al)^-1 * g );
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
            k := k * stb.prei;
            stb := AddIgsToIgs( stb.stab, Igs(N) );
            C := SubgroupByIgs( G, stb );
        fi;
    
    od;

    return rec(conj := k, kano := h);

end;

####################################################################
## Global function to calculate the canonical conjugate in G      ##
####################################################################

InstallGlobalFunction( "CanonicalConjugateNilGroup", function(G,g)

    return CanonicalConjugateNilGroupSeries(G,g,PcpsOfEfaSeries(G));

end );