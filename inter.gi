#######################################################################
## Helper function to calculate the intersection of U with a term    ##
## of the efa series of G                                            ##
#######################################################################

InstallGlobalFunction( "IntersectionSeriesTerm", function(U, term)

    local   gens,   #Generators of U
            pos;    #Position of the first generator that is in U and term

    #Trivial case
    if Size(term) = 1 or Size(U) = 1 then 
        return rec( cross := Subgroup(term!.ParentAttr, []), gens :=[] ); 
    fi;
    
    gens := Cgs(U);
    pos  := 1;
    while ( not gens[pos] in term ) do
        pos := pos+1;
        if pos > Length(gens) then 
            return rec( cross := Subgroup( term!.ParentAttr, [ ] ),
                        gens  := [ One( term!.ParentAttr ) ] );
        fi;
    od;

    return rec( cross := Subgroup( term!.ParentAttr, gens{[pos..Length(gens)]} ),
                gens  := gens{[pos..Length(gens)]} );

end );

#######################################################################
## Global function to calculate the intersection of U with a the     ##
## terms of the efa series of G                                      ##
#######################################################################

InstallGlobalFunction( "InducedIntersectionSeries", function(G, U) 
    
    local   series,     #Series of G
            iseries,    #Intersection with the series
            i,          #Term of the series
            cross;      #Intersection term

    series  := PcpSeries(G);
    iseries := [];
    i       := 1;
    cross   := IntersectionSeriesTerm(U, series[i]).cross;

    while Size(cross) <> 1 do 
        if not cross in iseries then Add(iseries, cross); fi;
        i     := i + 1;
        cross := IntersectionSeriesTerm(U, series[i]).cross;
    od;
    
    return iseries;

end );

#######################################################################
## Global function to calculate the pcps of the induced series       ##
##                                                                   ##
#######################################################################

InstallGlobalFunction( "PcpsOfInducedIntersectionSeries", function(G, U)

    local   series, #Series of G
            cross,  #Intersection series
            pcps,   #Pcps to return
            i;      #Bucle variable

    series := PcpSeries(G);
    cross  := InducedIntersectionSeries(U, series);
    pcps   := [];

    for i in [1..Length(cross)-1] do
        Add(pcps, Pcp( cross[i], cross[i+1] ) );
    od;

    return pcps;

end );

#######################################################################
## Local function to calculate the intersection of two cyclic        ##
## subgroups of G                                                    ##
#######################################################################

IntersectionCyclicSubgroupsNilGroups := function(G, U, V)

    local   gU,
            gV,
            d1,
            d2,
            e,
            gen;
    
    gU := Cgs(U)[1];
    #Trivial case
    if U = V then
        return [gU];
    fi;

    #Cyclic case
    gV := Cgs(V)[1];
    d1 := Depth( gU );
    d2 := Depth( gV );
    if d1 <> d2 then
        return [ One(G) ];
    else
        e := Lcm( Exponents( gU )[d1], Exponents( gV )[d2] );
        gU := gU^e;
        gV := gV^e;

        if gU in V and gV in U then
            if gU = gV then
                return [ gU ];
            else
                return [ Cgs(G)[d1]^e ];
            fi;
        else
            return [ One(G) ];
        fi;
    fi;

end;

#######################################################################
## Local function to calculate the intersection of two subgroups     ##
#######################################################################

IntersectionSubgroupsNilGroupsSeries := function(G, U, V, ser, Gn, gn, d)

    local   G2,     #Second term of ser
            U2,     #Intersection U and G2
            V2,     #Intersection V and G2
            I,      #Intersection U2 and V2
            gens,   #Generators of I
            p,      #Projection map
            pI,     #Intersection of p(U) and p(V)
            B,      #Generators of pI
            Un,     #Intersection U and Gn
            Vn,     #Intersection V and Gn
            an,     #Additional values for A
            H,      #Values hu for all generators
            A,      #Exponents of gn to equalize hu and hv
            i,j,    #Bucle variable
            b,      #Preimage of the generator of pI
            hu,     #Sifting of b in U
            hv,     #Sifting of b in V
            R,      #Representation of d2
            d2,     #Gcd of A{2,...n}
            d1,     #Gcd of A
            x;      #Value to add to the generators

    #Trivial cases
    if U = V then
        return Cgs(U);

    elif Size(U) = 1 then
        return [ One(G) ];

    elif Size(V) = 1 then
        return [ One(G) ];

    fi;

    #Cyclic case
    if Length( Cgs(U) ) = 1 and Length( Cgs(V) ) = 1 then
        return IntersectionCyclicSubgroupsNilGroups(G, U, V);
    fi;

    #General Case
    
    i := 2;
    G2 := ser[i];
    U2  := IntersectionSeriesTerm(U, G2).cross;
    V2  := IntersectionSeriesTerm(V, G2).cross;

    while U2 = U and V2 = V do
        i  := i + 1;
        G2 := ser[i];
        U2 := IntersectionSeriesTerm(U, G2).cross;
        V2 := IntersectionSeriesTerm(V, G2).cross;
    od;

    I   := IntersectionSubgroupsNilGroupsSeries(G, U2, V2, ser, Gn, gn, d);
    gens:= ShallowCopy(I);
    p   := NaturalHomomorphismByNormalSubgroup(G, Gn);
    pI  := IntersectionSubgroupsNilGroups( Image(p), p(U), p(V) );
    B   := pI;
    
    if B[1] in p(G2) then
        return I;
    
    else

        Un := IntersectionSeriesTerm(U, Gn).gens[1];
        Vn := IntersectionSeriesTerm(V, Gn).gens[1];
        an := [ Exponents(Un)[d], Exponents(Vn)[d], FactorOrder(gn)];

        H := [];
        A := [];

        for j in [1..Length(B)] do

            b    := PreImagesRepresentative( p, B[j]);
            hu   := b * (Sifting(G, U, b).y)^-1;
            H[j] := hu;
            hv   := b * (Sifting(G, V, b).y)^-1;
            A[j] := Last(Exponents(hu))-Last(Exponents(hv));
            
        od;
        
        Append( A, an );
        R  := GcdRepresentation( A{[2..Length(A)]} );
        d2 := A{[2..Length(A)]} * R;
        d1 := Gcd( A );
        
        if d2 = 0 and d1 <> 0 then
            

        elif d1 = 0 then
            I := Subgroup(G, I);
            x := Sifting(G, I, H[1]).y;
            Add( gens, x, 1);

        else

            b := d2/d1;
            B := R{[1..Length(R)-2]}*(-A[1]/d1);
            B[Length(B)] := B[Length(B)] * an[1];
            Add(B, b, 1);
            Add(H, gn);
            x := MappedVector(B, H);

            if x in G2 then 

            else
                I := Subgroup(G, I);
                x := Sifting(G, I, x).y;
                Add(gens, x, 1);
            fi;
        fi;
    fi;

    return gens;
    

end;

#######################################################################
## Global function to calculate the intersection of two subgroups    ##
#######################################################################

InstallGlobalFunction( "IntersectionSubgroupsNilGroups", function(G, U, V) 

    local   ser,    #Pcp series of G
            Gn,     #Last term of ser
            gn,     #Generator of gn
            d;      #Depth of gn
    
    ser := PcpSeries(G);
    Gn  := ser[Length(ser)-1];
    gn  := Pcp(Gn)[1];
    d   := Depth(gn);

    return IntersectionSubgroupsNilGroupsSeries(G, U, V, ser, Gn, gn, d);


end );