dir      := DirectoryCurrent();
filename := Filename( dir, "GapLog.txt");
PrintTo( filename, "Starting to log gap events\n");

############################################################################################
##  Experiments for canonical conjugacy with small nilpotent groups                       ##
############################################################################################

AppendTo( filename, "Results for the first experiment\n");
Print("Starting experiment one \n");
for n in [2..5] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);

    for i in [1..100] do
        g := Random(G);
        h := g^Random(G);
        
        t := Runtime();
        IsConjugate(G, g, h);
        t := Runtime() - t;
        Add(ts2, t);

        elms := [g,h];
        t := Runtime();
        can := IsCanonicalConjugateElements(G, elms);
        t := Runtime() - t;
        Add(ts1, t);
        for k in [1..Length(elms)] do
            if elms[k]^can.conj[k] <> can.kano then Error("Incorrect conjugating element."); fi;
        od;

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment one \n");

############################################################################################
##  Experiments for canonical conjugacy with large nilpotent groups                       ##
############################################################################################

AppendTo( filename, "Results for the first experiment\n");
Print("Starting experiment two \n");
for n in [6,7] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);

    for i in [1..100] do
        g := Random(G);
        h := g^Random(G);
        
        t := Runtime();
        IsConjugate(G, g, h);
        t := Runtime() - t;
        Add(ts2, t);

        elms := [g,h];
        t := Runtime();
        can := IsCanonicalConjugateElements(G, elms);
        t := Runtime() - t;
        Add(ts1, t);
        for k in [1..Length(elms)] do
            if elms[k]^can.conj[k] <> can.kano then Error("Incorrect conjugating element."); fi;
        od;

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment two \n");

############################################################################################
##  Experiments for canonical conjugacy of subgroups with small nilpotent groups          ##
############################################################################################

AppendTo( filename, "Results for the first experiment\n");
Print("Starting experiment three \n");
for n in [2..5] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);

    for i in [1..100] do
        U := Subgroup( G, [Random(G), Random(G)]);
        V := U^Random(G);
        
        t := Runtime();
        IsConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts2, t);

        t := Runtime();
        can := IsCanonicalConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts1, t);
        if IsBool(can) then Error("No conjugating."); fi;

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment three \n");

############################################################################################
##  Experiments for canonical conjugacy of subgroups with large nilpotent groups          ##
############################################################################################

AppendTo( filename, "Results for the first experiment\n");
Print("Starting experiment four \n");
for n in [8..9] do
    ts1 := [];
    ts2 := [];
    G := SomeNilpotentGroups(n);

    for i in [1..10] do
        U := Subgroup(G, [RandomElementRangeGenerators(G,10), RandomElementRangeGenerators(G, 10+n)]);
        V := U^Random(G);
        
        t := Runtime();
        IsConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts2, t);

        t := Runtime();
        can := IsCanonicalConjugateSubgroups(G, U, V);
        t := Runtime() - t;
        Add(ts1, t);
        if IsBool(can) then Error("No conjugating."); fi;

    od;

    AppendTo( filename, Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n, "\n");
    AppendTo( filename, Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n, "\n");
od;
Print("Finished experiment four \n");
Unbind(dir);
Unbind(filename);