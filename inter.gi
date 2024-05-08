#######################################################################
## Helper function to calculate the intersection of U with a term    ##
## of the efa series of G                                            ##
#######################################################################

InstallGlobalFunction( "IntersectionSeriesTerm", function(U, term)

    local   gens,   #Generators of U
            filt,   #Filtered generators that are in U and term
            pos;    #Position of the first generator that is in U and term

    #Trivial case
    if Size(term) = 1 then 
        return term; 
    fi;
    
    gens    := Igs(U);
    filt    := Filtered( gens, g -> g in term );

    if IsEmpty(filt) then
        #We have that the intersection is empty
        pos := Length(gens)+1; 
    else
        pos := Position(gens, filt[1]);
    fi;

    return Subgroup( term, gens{[pos..Length(gens)]} );

end );

#######################################################################
## Global function to calculate the intersection of U with a the     ##
## terms of the efa series of G                                      ##
#######################################################################

InstallGlobalFunction( "InducedIntersectionSeries", function(U, series) 
    
    local   iseries,    #Intersection with the series
            term,       #Term of the series
            iterm;      #Intersection term

    iseries := [];
    
    for term in series do
        iterm := IntersectionSeriesTerm(U, term);

        if Size(iterm) <> 1 and (not iterm in iseries) then
            Add(iseries, iterm);

        elif Size(iterm) = 1 then
            Add(iseries, iterm);
            break;
        fi;
    od;
    
    return iseries;

end );

#######################################################################
## Global function to calculate the pcps of the induced series       ##
##                                                                   ##
#######################################################################

InstallGlobalFunction( "PcpsOfInducedIntersectionSeries", function(U, series)

    local   iseries,    #Intersection series
            pcps,       #Pcps to return
            i;          #Bucle variable

    iseries := InducedIntersectionSeries(U, series);
    pcps := [];

    for i in [1..Length(iseries)-1] do
        Add(pcps, Pcp( iseries[i], iseries[i+1] ) );
    od;

    return pcps;

end );

IntersectionSubgroupsNilGroups := function(G, U, V)

    local   series, #Pcp series of G
            Gn,     #Last term of pcp series
            Un,
            Vn;

    series := Reversed( PcpSeries(G) );
    Gn     := series[1];

    Un     := IntersectionSeriesTerm(U, Gn);
    Vn     := IntersectionSeriesTerm(V, Gn);

end;