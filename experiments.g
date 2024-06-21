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