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

#######################################################################
## Global function to calculate the a reduced element in a basis     ##
#######################################################################

InstallGlobalFunction( "ReducePcpElement", function( elm, basis )

    local   exp,    #Exponents of elm
            l,      #Leading exponents of the elements in the basis
            d,      #Depth of the elements in the basis
            v,      #Reduced preimage
            i,      #Bucle variable
            r,      #Residuo
            x;      #Expoenent of the element in the basis

    v     := elm ;
    basis := Reversed(basis);
    exp   := Exponents( v );
    l     := List( basis, LeadingExponent );
    d     := List( basis, Depth);
    

    for i in [1..Length(d)] do

        if exp[ d[i] ]<0 or exp[ d[i] ]>(l[i]-1) then
            r   :=  exp[ d[i] ] mod l[i];
            x   := (exp[ d[i] ] - r)/l[i];
            v   := v * basis[i]^-x;
            exp := Exponents( v );
        fi;
    od;

    return v;

end );

####################################################################
## Generates a random element of G using only generators between  ##
##                                                      m and n   ##
####################################################################

InstallGlobalFunction( "RandomElementRangeGenerators", 
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
                while g[i] = 0 do
                    g[i] := Random( Integers );
                od;
            else
                g[i] := Random( 1, rel[i]-1 );
            fi;
        else
            g[i] := 0;
        fi;
    od;

    return MappedVector( g, pcp );

end );

####################################################################
## Generates a random subgroup of G                               ##
## if the input n is given generates a subgroup with n generators ##
####################################################################

InstallGlobalFunction( "RandomSubgroup", function( arg )

    local   G,
            cgs,
            n,
            gens,
            nums,
            r,
            g,
            i,
            U;

    G    := arg[1];
    cgs  := Cgs(G);
    gens := [];
    nums := [];

    if Length(arg) > 1 then
        n := arg[2];
        if n > Length( cgs ) then Error("The number n is greater than the number of generators of G."); fi;
    else
        n := Random( [1..Length( cgs )] );
    fi;


    for i in [1..n] do
        r := Random([1..Length(cgs)]);
        while r in nums do
            r := Random([1..Length(cgs)]);
        od;
        Add( nums, r );
    od;

    Sort(nums);
    
    for i in [1..n] do
        g := RandomElementRangeGenerators(G, nums[i]);
        Add( gens, g );
    od;

    U := Subgroup(G, gens);
    return U;

end );

####################################################################
### Returns true if is less or equal in the order                ###
### 0 << 1 << ... << -1 << ...                                   ###
####################################################################

IntegerOrder := function(a,b)
    if a = 0 then
        return true;
    elif b = 0 then
        return false;
    elif a > 0 then
        if b < 0 then
            return true;
        else
            return a <= b;
        fi;
    else
        if b > 0 then
            return false;
        else
            return b <= a;
        fi;
    fi;
end;   

####################################################################
### Returns true if is less or equal in the order                ###
### extended to the exponents of g and h                         ###
####################################################################

ExponentOrder := function(g,h)
    local   e1,
            e2,
            i;

    e1  := Exponents(g);
    e2  := Exponents(h);

    for i in [1..Length(e1)] do
        if e1[i] <> e2[i] then
            return IntegerOrder( e1[i], e2[i] );
        fi;
    od;
    return true;
end;

####################################################################
### Returns true if g <<= h                                      ###
####################################################################

InstallGlobalFunction( "ConjugacyOrder" , function(g,h) 

    local   fam,ord;

    fam := FamilyObj( g );
    ord := OrderingByLessThanFunctionNC(fam, ExponentOrder);

    return IsLessThanUnder(ord, g, h);

end );

####################################################################
### Returns true if g <<= h                                      ###
####################################################################

InstallGlobalFunction( "Sifting", function(G, U, g)

    local   gen,
            n,
            y,
            B,
            d,
            exp,
            alp,
            l,
            i;

    gen := Cgs(U);
    n   := Length(gen);
    y   := g;
    B   := List( gen, x -> 0);
    d   := List( gen, Depth );
    exp := List( gen, Exponents );
    exp := List( [1..n], x -> exp[x][d[x]]);
    alp := List( [1..n], x -> Exponents(y)[d[x]]);
    l   := List( [1..n], x -> IntegerOrder( alp[x], exp[x]) );
    
    while not ForAll( l, x -> x = true ) do
        i    := PositionProperty(l, x -> x = false );
        B[i] := Int( Floor( Float( alp[i] / exp[i] ) ) );
        y    := (gen[i] ^ -B[i]) * y;
        alp := List( [1..n], x -> Exponents(y)[d[x]]);
        l   := List( [1..n], x -> IntegerOrder( alp[x], exp[x]) );
    od;

    return rec( y := y, B := B);

end );
