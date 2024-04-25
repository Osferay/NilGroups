############################################################################################
##  Experiments with small nilpotent groups                                               ##
############################################################################################
for n in [1..5] do
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
        can := IsCanonicalConjugateNilGroup(G, elms);
        t := Runtime() - t;
        Print(t, "\n");
        Add(ts1, t);
        for k in [1..Length(elms)] do
            if elms[k]^can.conj[k] <> can.kano then Error("Incorrect conjugating element."); fi;
        od;

    od;

    Print( Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n);
    Print( Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n);
od;

############################################################################################
##  Experiments with large nilpotent groups                                               ##
############################################################################################
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
        can := IsCanonicalConjugateNilGroup(G, elms);
        t := Runtime() - t;
        Print(t, "\n");
        Add(ts1, t);
        for k in [1..Length(elms)] do
            if elms[k]^can.conj[k] <> can.kano then Error("Incorrect conjugating element."); fi;
        od;

    od;

    Print( Float( Sum(ts1)/Length(ts1) ), "Time consumed by new algorithm for n equal to ", n);
    Print( Float( Sum(ts2)/Length(ts2) ), "Time consumed by old algorithm for n equal to ", n);
od;