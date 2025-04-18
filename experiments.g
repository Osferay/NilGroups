########################################################################################
##  Examples of finitely generated nilpotent groups used to check performance         ##
########################################################################################
LoadPackage( "nq" );

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

dir      := DirectoryCurrent();
filename := Filename( dir, "GapLog.txt" );
log      := Filename( dir, "ExperimentsLog.txt" );
PrintTo( filename, "Starting to log gap events\n" );
PrintTo( log     , "Login events\n");

############################################################################################
##  Experiments for canonical conjugacy with small nilpotent groups                       ##
############################################################################################

AppendTo( filename, "Results for the first experiment\n");
Print("Starting experiment one \n");
for n in [1..4] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);
    AppendTo( log, "Group generated.\n");

    for i in [1..100] do
        g := RandomElementRangeGenerators(G, 1);
        h := g^Random(G);
        AppendTo(log, "Elements ready." );
        
        t := Runtime();
        IsConjugateNilGroup(G, g, h);
        t := Runtime() - t;
        Add(ts2, t);
        AppendTo( log, " Old algorithm done", t );

        elms := [g,h];
        t := Runtime();
        can := IsCanonicalConjugateElements(G, elms);
        t := Runtime() - t;
        Add(ts1, t);
        AppendTo( log, " New algorithm done", t, "\n");

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment one \n");

############################################################################################
##  Experiments for canonical conjugacy with large nilpotent groups                       ##
############################################################################################

AppendTo( filename, "Results for the second experiment\n");
Print("Starting experiment two \n");
for n in [5,6] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);
    AppendTo( log, "Group generated.\n");

    for i in [1..100] do
        g := Random(G);
        h := g^Random(G);
        AppendTo(log, "Elements ready." );
        
        t := Runtime();
        IsConjugateNilGroup(G, g, h);
        t := Runtime() - t;
        Add(ts2, t);
        AppendTo( log, " Old algorithm done", t );

        elms := [g,h];
        t := Runtime();
        can := IsCanonicalConjugateElements(G, elms);
        t := Runtime() - t;
        Add(ts1, t);
        AppendTo( log, " New algorithm done", t, "\n");

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment two \n");

############################################################################################
##  Experiments for canonical conjugacy of subgroups with small nilpotent groups          ##
############################################################################################

AppendTo( filename, "Results for the third experiment\n");
Print("Starting experiment three \n");
for n in [1..4] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);
    AppendTo( log, "Group generated.\n");

    for i in [1..100] do
        U := Subgroup( G, [Random(G), Random(G)]);
        V := U^Random(G);
        AppendTo(log, "Subgroup ready." );
        
        t := Runtime();
        IsConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts2, t);
        AppendTo( log, " Old algorithm done", t );

        t := Runtime();
        can := IsCanonicalConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts1, t);
        AppendTo( log, " New algorithm done", t, "\n");

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment three \n");

############################################################################################
##  Experiments for canonical conjugacy of subgroups with large nilpotent groups          ##
############################################################################################

AppendTo( filename, "Results for the fourth experiment\n");
Print("Starting experiment four \n");
for n in [7,8] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);
    AppendTo( log, "Group generated.\n");

    for i in [1..10] do
        U := Subgroup(G, [RandomElementRangeGenerators(G,10), RandomElementRangeGenerators(G, 10+n)]);
        V := U^Random(G);
        AppendTo(log, "Subgroup ready." );
        
        t := Runtime();
        IsConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts2, t);
        AppendTo( log, " Old algorithm done", t );

        t := Runtime();
        can := IsCanonicalConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts1, t);
        AppendTo( log, " New algorithm done", t, "\n");

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment four \n");

############################################################################################
##  Experiments for canonical conjugacy of list of elements with small nilpotent groups   ##
############################################################################################

AppendTo( filename, "Results for the fifth experiment\n");
Print("Starting experiment five \n");
for n in [1..4] do
    G := SomeNilpotentGroups(n);
    ts1 := [];
    ts2 := [];
    AppendTo( log, "Group generated.\n");
    
    for i in [1..100] do
        l := RandomListElements(G, 3, "no_id");
        list := ShallowCopy(l);

        for i in [4..20] do
            g := Random(l);
            Add(list, Random(l)^Random(G));
        od;
        AppendTo(log, "List ready." );

        t := Runtime();
        IsConjugateList(G, list);
        t := Runtime() - t;
        Add(ts2, t);
        AppendTo( log, " Old algorithm done", t );

        t := Runtime();
        can := CanonicalConjugateList(G, list);
        t := Runtime() - t;
        Add(ts1, t);
        AppendTo( log, " New algorithm done", t, "\n");

    od;
    
    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment five. \n");

Unbind(dir);
Unbind(log);
Unbind(filename);