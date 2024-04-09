#######################################################################
## Helper function to calculate the intersection of U with a term    ##
## of the efa series of G                                            ##
#######################################################################

InstallGlobalFunction( "IntersectionEfaTerm", function(U, term)

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

InstallGlobalFunction( "InducedEfaSeries", function(G, U) 
    
    local   efa,
            iefa,
            term,
            iterm;

    efa  := EfaSeries(G);
    iefa := [];
    
    for term in efa do
        iterm := IntersectionEfaTerm(U, term);

        if Size(iterm) <> 1 and (not iterm in iefa) then
            Add(iefa, iterm);
        elif Size(iterm) = 1 then
            Add(iefa, iterm);
            break;
        fi;
    od;
    
    return iefa;

end );

InstallGlobalFunction( "PcpsOfInducedEfaSeries", function(G,U)

    local   iefa,
            pcps,
            i;

    iefa := InducedEfaSeries(G, U);
    pcps := [];

    for i in [1..Length(iefa)-1] do
        Add(pcps, Pcp( iefa[i], iefa[i+1] ) );
    od;

    return pcps;

end );