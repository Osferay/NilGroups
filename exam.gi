########################################################################################
##  Examples of finitely generated nilpotent groups used to check performance         ##
########################################################################################

SomeNilpotentGroups := function( n )
    local   G,          #The group to return
            et,         #Expresion tree
            a,b,c,d,e,x,#Expresion tree variables
            rengel,     #Right engel relations
            H,          #Denominator of the nilpotent qoutient
            ftl;        #From the left collectors
    
    if n = 1 then 
        #Torsion-free of hirsch length 5
        ftl := FromTheLeftCollector(5);
        SetConjugate( ftl, 2, 1, [2, 1, 3, 11] );
        SetConjugate( ftl, 3, 1, [3, 1, 4, 8]  );
        SetConjugate( ftl, 3, 2, [3, 1, 5, 15] );
        SetConjugate( ftl, 4, 1, [4, 1, 5, 16] );
        G := PcpGroupByCollector( ftl );
        return G;
    
    elif n = 2 then
        #[0, 0, 0, 0, 0, 4, 2, 2]
        et := ExpressionTrees( "a", "b", "x" );
        a := et[1];; b := et[2];; x := et[3];;
        rengel := LeftNormedComm( [a,x,x,x] );
        H := rec( generators := et, relations := [rengel] );
        G := NilpotentQuotient( H, [x] );
        return G;
    
    elif n = 3 then
        #[0, 0, 0, 0, 3872]
        ftl := FromTheLeftCollector(5);
        SetRelativeOrder( ftl, 5, 3872);
        SetConjugate( ftl, 2 , 1, [2, 1, 3, 22, 4, 88]);
        SetConjugate( ftl, 3 , 1, [3, 1, 4, 16, 5, 128]);
        SetConjugate( ftl, 3 , 2, [3, 1, 5, 15]);
        SetConjugate( ftl, 4 , 1, [4, 1, 5, 352]);
        G := PcpGroupByCollector( ftl );
        return G;

    elif n = 4 then
        ftl := FromTheLeftCollector(5);
        SetRelativeOrder( ftl, 5, 11264);
        SetRelativeOrder( ftl, 4, 352);
        SetPower( ftl, 4, [5,5120]);
        SetConjugate( ftl, 2 , 1, [2, 1, 3, 22, 4, 88]);
        SetConjugate( ftl, 3 , 1, [3, 1, 4, 16, 5, 128]);
        SetConjugate( ftl, 3 , 2, [3, 1, 5, 15]);
        SetConjugate( ftl, 4 , 1, [4, 1, 5, 32]);
        G := PcpGroupByCollector( ftl );
        return G;
    
    elif n = 5 then

        et := ExpressionTrees( "a", "b", "c", "d", "x" );
        a := et[1];; b := et[2];; c := et[3];; d := et[4];; x := et[5];;
        rengel := LeftNormedComm( [a,x,x,x] );
        H := rec( generators := et, relations := [rengel^3, b^356] );
        G := NilpotentQuotient( H, [x], 4);
        return G;

    elif n = 6 then
        et := ExpressionTrees( "a", "b", "x" );
        a := et[1];; b := et[2];; x := et[3];;
        rengel := LeftNormedComm( [a,x,x,x,x,x] );
        H := rec( generators := et, relations := [rengel] );
        G := NilpotentQuotient( H, [x], 9 );
        return G;

    elif n = 7 then 
        et := ExpressionTrees( 7 );
        a := et[1];; b := et[2];; c := et[3];; d := et[4];; e := et[5];; x := et[6];;
        rengel := LeftNormedComm( [a,x,x,x] );
        H := rec( generators := et, relations := [rengel^2, b^625, d^80, e^512] );
        G := NilpotentQuotient( H, [x], 2);
        return G;

    elif n = 8 then
        et := ExpressionTrees( 7 );
        a := et[1];; b := et[2];; c := et[3];; x := et[4];;
        rengel := LeftNormedComm( [a,x,x,x] );
        H := rec( generators := et, relations := [rengel, b^585*c^3, a^255] );
        G := NilpotentQuotient( H, [x], 2 );
        return G;

    else

        return fail;

    fi;

end;