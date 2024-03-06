########################################################################################
##  Series for finitely generated torsion-free nilpotent groups, T-Groups in short    ##
########################################################################################

InstallGlobalFunction( "IsolatorSeries", function(G)

    local   C,      #Torsion subgroup
            i,      #bucle variable
            nat,    #Natural homomorphism G-->G/N 
            S,      #Isolator central series of G
            lcs;    #Lower central series of G

    lcs := LowerCentralSeriesOfGroup(G);
    S := [];
    C := G;

    for i in [1..Length(lcs)] do
        nat := NaturalHomomorphismByNormalSubgroup( G, lcs[i] );
        C := PreImages( nat, TorsionSubgroup( Image( nat ) ) );
        Add( S, C );
    od;

    return S;
end );

InstallGlobalFunction( "PcpsOfIsolatorSeries", function(G)

    return PcpsBySeries(IsolatorSeries(G));

end );

InstallGlobalFunction( "TauType", function(G)

    local   ics,
            i,      #bucle variable
            t;

    ics := IsolatorSeries(G);
    t := [];

    for i in [1..( Length(ics)-1 )] do
        Add(t, HirschLength( ics[i]/ics[i+1] ));
    od;

    return t;

end );

InstallGlobalFunction( "TwoStepCentralizer", function(G,i)

    local   ics,
            nat,
            C;

    ics := IsolatorSeries(G);

    if i > (Length(ics)-2) then
        Error( " Second argument too large. " );
    fi;

    nat := NaturalHomomorphismByNormalSubgroup(G, ics[i+2]);
    C := PreImages( nat, CentralizerPcpGroup( Image( nat ), Image ( nat, ics[i] )));

    return C;

end );

InstallGlobalFunction( "ExtendedIsolatorSeries", function(G)

    local   eis,
            C;

    eis := IsolatorSeries(G);
    C := TwoStepCentralizer(G, Length(eis)-2);

    if C<>eis[2] then
        Add(eis, C, 2);
    fi;

    return eis;

end );

InstallGlobalFunction( "ExtendedTauType", function(G)

    local   eis,
            i,      #bucle variable
            t;

    eis := ExtendedIsolatorSeries(G);
    t := [];
    
    for i in [1..( Length(eis)-1 )] do
        Add(t, HirschLength( eis[i]/eis[i+1] ));
    od;

    return t;

end );