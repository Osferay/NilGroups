InstallGlobalFunction( IsomorphismTauGroup, function(G,H)

    local   TauT,   #Tau type of the groups
            v,      #vector of relations of G
            w,      #vector of relations of H
            g;      #Auxiliar variable

    if not ( ( HirschLength(G)=HirschLength(H) ) or ( TauType(G)=TauType(H) ) ) then
        return false;

    else
        if ( IsFreeAbelian(G) and IsFreeAbelian(H) ) then
            return true;
        else
            TauT    := TauType(G);
            v       := TauVector(G);
            w       := TauVector(H);

            if TauType = [2,1,1] then
                g := Gcd(v[1], v[3]);
                return ( w[1] = AbsoluteValue(v[1]) and w[3] = AbsoluteValue(w[3]) and ( w[4] = v[4] mod g or w[4] = -v[4] mod g) );

            elif TauType = [3,1,1] then
                g := Gcd(v[6], v[8]);
                return ( w[6] = AbsoluteValue(v[6]) and w[8] = AbsoluteValue(w[8]) );

            elif TauType = [2,1,1,1] then
                for i in [1..Length(v)] do
                    if not( AbsoluteValue(v[i]) = AbsoluteValue(w[i]) ) then
                        return false;
                    else 
                        for i in [1,4,6,8] do
                            if v[i]<0 then
                                return false;
                            else
                                return true;
                            fi;
                        od;
                    fi;
                od;
                
            elif TauType = [2,1,2] then

            elif TauType = [3,2] then

            elif TauType = [HirschLength(G)-1,1] then

            else
                Error( " The function only works for groups with small Hirsch number <= 5. ");
            fi;
        
        
        fi;
    fi;

end );

InstallGlobalFunction( SameGenus, function(G,H)

    local   TauT,   #TauType of the groups
            t,s,    #TauVectors of G and H
            d;      #Auxiliar variable


    #Catch trivial case of class less than 3 or diferent types
    if Length(TauType(G))<3 or ( TauType(G) <> TauType(H) )then
        return false
    
    else
        TauT := TauType(G);
        t := TauVector(G);
        s := TauVector(H);

        if TauT = [2,1,1] then

            d := Gcd(t[1],t[3]);

            return (t[1] = s[1]) and (t[3] = s[3]) and (Gcd(t[2],d) = Gcd(s[2],d));

        elif TauT = [3,1,1] then

            d := Gcd(t[6],t[8]);
            e := Gcd(t[2],t[5],d);

            return (t[2] = s[2]) and (t[6] = s[6]) and (t[8] = s[8]) and (Gcd(t[5],d) = Gcd(s[5],d)) and (Gcd(t[3],d) = Gcd(s[3],d));

        elif TauT = [2,1,1,1] then

            d := Gcd(t[1],t[6],t[8]);

            return (t[1] = s[1]) and (t[4] = s[4]) and (t[6] = s[6]) and (t[8] = s[8]); 


end;