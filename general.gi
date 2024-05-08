####################################################################
## IsTGroup(G), checks if the group is nilpotent and torsion free ##
####################################################################

InstallMethod(  IsTauGroup, 
                "for pcp groups", [IsPcpGroup],
function( G )

    return ((IsNilpotentGroup(G) and IsTorsionFree(G) ));

end );

####################################################################
## Pref(x, y), generates the prefix of the reduction (+-x) mod y  ##
####################################################################

InstallGlobalFunction( Pref, function(x,y)

    local   s,
            pref;

    s := Minimum( x mod y, -x mod y);
    pref :=[];

    if IsInt((x-s)/y) then
        Add(pref, 1);
    fi;

    if IsInt((-x-s)/y) then
        Add(pref, -1);
    fi;

    return pref;

end );

####################################################################
## TauVector(G), generates the vector of relations of G           ##
####################################################################

InstallGlobalFunction( TauVector, function(G)

    local   pcp,
            i,j,k,  #Bucle variables
            v, 
            h;       #Hirsch number of G.

    if not (IsPcpGroup(G) and IsTauGroup(G)) then
        Error( "This method is non compatible with the group" );
    fi;

    h   := HirschLength(G);
    pcp := Pcp(G);
    v   := [];

    for i in [1..(h-2)] do
        for j in [(i+1)..(h-1)] do
            for k in [(j+1)..h] do
                Add( v,Exponents( pcp[j]^pcp[i] )[k] );
            od;
        od;
    od;

    return v;

end );

####################################################################
## Generates the matrix of relations of G using a vector v        ##
##                                      and Hirsch length h       ##
####################################################################

InstallGlobalFunction( MatrixRelationsByVector, function(v,h)

    local   t,      #Matrix of relations
            cont,   #count of the iteration
            N,      #Auxiliar matrix
            i,j,k;  #Bucle variables

    cont := 1;
    t := [];

    for i in [1..(h-2)] do
        N := NullMat(h,h);

        for j in [(i+1)..(h-1)] do

            for k in [(j+1)..h] do
                N[j][k] := v[cont];
                cont := cont + 1;
            od;

        od;

        Add(t,N);
    od;

    return t;

end );

####################################################################
## Generates the torsion free group G using a vector v and        ##
##                                          Hirsch length h       ##
####################################################################

InstallGlobalFunction( TauGroupByVector, function(v,h)

    local   cont,   #Count of iteration
            ftl,    #Collector of relations
            rel,    #Variable auxiliar
            G,      #Group to return
            i,j,k;  #Bucle variables

    #Create the relations
    ftl := FromTheLeftCollector(h);
    cont := 1;

    for i in [1..(h-2)] do

        for j in [(i+1)..(h-1)] do

            rel := [j,1];
            for k in [(j+1)..h] do

                if v[cont] <> 0 then
                    Add(rel,k);
                    Add(rel,v[cont]);
                fi;

                cont := cont + 1;
            od;

        SetConjugate( ftl, j, i, rel);
        od;
        
    od;

    #Create the group
    G := PcpGroupByCollector( ftl );

    return G;

end );

####################################################################
## Generates a random element of G using only generators between  ##
##                                                      m and n   ##
####################################################################

InstallGlobalFunction(  RandomElementRangeGenerators, 
                                function( arg )

    local   G,
            m,
            n,
            pcp,
            rel,
            i,
            g;

    G := arg[1];
    n := arg[2];
    pcp := Pcp(G);

    if Length( arg ) > 2 then
        m := arg[3];
    else
        m := Length( pcp );
    fi;

    pcp := Pcp(G);
    rel := RelativeOrdersOfPcp( pcp );
    g   := [];

    if m > Length( pcp ) then 
        Error("The number m is greater than the number of generators of G."); 
    fi;

        
    if Length( pcp ) = 0 then
        return One( G );
    fi;

    for i in [1..Length( pcp )] do
        if i in [n..m] then
            if rel[i] = 0 then
                g[i] := Random( Integers );
            else
                g[i] := Random( 0, rel[i]-1 );
            fi;
        else
            g[i] := 0;
        fi;
    od;

    return MappedVector( g, pcp );

end);

####################################################################
## Generates a random subgroup of G                               ##
## if the input n is given generates a subgroup with n generators ##
####################################################################

InstallGlobalFunction( RandomSubgroup, function( arg )

    local   G,
            igs,
            n,
            gens,
            nums,
            r,
            i;

    G    := arg[1];
    igs  := Igs(G);
    gens := [];
    nums := [];

    if Length(arg) > 1 then
        n := arg[2];
        if n > Length( igs ) then Error("The number n is greater than the number of generators of G."); fi;
    else
        n := Random( [1..Length( igs )] );
    fi;


    for i in [1..n] do
        r := Random([1..Length(igs)]);
        while r in nums do
            r := Random([1..Length(igs)]);
        od;
        Add( nums, r );
    od;

    Sort(nums);

    for i in [1..n] do
        Add( gens, RandomElementRangeGenerators(G, nums[i]) );
    od;

    return Subgroup(G, gens);

end );