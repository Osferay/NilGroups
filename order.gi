IntegerOrder := function(a,b)
    if a = 0 then
        return true;
    elif b = 0 then
        return false;
    elif a > 0 then
        if b < 0 then
            return true;
        else
            return a<b;
        fi;
    else
        if b > 0 then
            return false;
        else
            return b<a;
        fi;
    fi;
end;   

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
    return false;
end;

InstallGlobalFunction( "ConjugacyOrder" , function(g,h) 

    local   fam,ord;

    fam := FamilyObj( g );
    ord := OrderingByLessThanFunctionNC(fam, ExponentOrder);

    return IsLessThanUnder(ord, g, h);

end );

PcpOrder := function(U,V)
    local   pcp1,
            pcp2,
            i;

    pcp1    := Pcp(U);
    pcp1    := Reversed( AsList(pcp1) );
    pcp2    := Pcp(V);
    pcp2    := Reversed( AsList(pcp2) );

    for i in [1..Length( pcp1 )] do
        if pcp1[i] <> pcp2[i] then
            return ExponentOrder( pcp1[i], pcp2[i] );
        fi;
    od;
    return false;

end;

InstallGlobalFunction( "SubgroupOrder", function(U,V) 

    local   fam, ord;

    fam := FamilyObj( U );
    ord := OrderingByLessThanFunctionNC( fam, PcpOrder);

    return IsLessThanUnder(ord, U, V);

end );