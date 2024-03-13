####################################################################
## Deprecated function to calculate the centralizer of g in G     ##
####################################################################

# OLDCentralizerNilGroup := function(G,g)

#     local   efa,    #Elementary/Free abelian series of G
#             C,      #Centralizer of G/Gi
#             U,      #U/Gi = Ci
#             nat,    #Natural homomorphism G -> G/Gi
#             fac,    #factor G(i-1)/Gi
#             gens,   #Generators L/Gi
#             imgs,   #gens -> [g,gen]
#             i,      #Bucle variables
#             f;      #homomorphism L/Gi -> G(i-1)/Gi 

#     efa := EfaSeries(G);
#     C := G;
#     U := G;

#     for i in [2..Length(efa)-1] do

#         #Calculate G -> G/Gi
#         nat := NaturalHomomorphismByNormalSubgroup(G, efa[i+1]);

#         #Check a trivial case:
#         if IsAbelian( Image(nat) ) then
#             C := G/efa[i+1];
#             U := U;

#         else

#             #Take G(i-1)/Gi
#             fac := Image( nat, efa[i] );

#             #Create the homomorpism f:U/Gi -> G(i-1)/Gi
#             gens:= Igs( Image( nat, U ) );
#             imgs := List( gens , x -> Comm( Image(nat, g) ,x) );
#             f   := GroupHomomorphismByImages( Image( nat, U ), fac, gens, imgs);

#             #Calculate C_(G/Gi)(gGi)
#             C := Kernel(f);

#             #Calculate the preimage of C in G
#             U := PreImage(nat, C);

#         fi;
#     od;

#     #make a random check
#     if not Comm(g, Random(C)) = Identity(G) then
#         Error(" The function is not working properly. ");
#     fi;

#     return C;

# end;

####################################################################
## Deprecated function to calculate the centralizer of a set elms in G##
####################################################################

# InstallGlobalFunction( "MultiCentralizerNilGroup", function(G, elms)

#     local   C,i;

#     C := G;
#     for i in elms do
#         C := CentralizerNilGroup(C,i);
#     od;

#     return C;

# end );

####################################################################
## Deprecated function to the preimage of the conjugating element     ##
####################################################################

# PreImageConjugate := function(G, w, nat, gens, imgs)

#     local   e,      #Exponent vector of the images
#             k,      #Exponent of w
#             n,      #Exponent of the generator of U
#             j,
#             ndx,
#             vAi,
#             v;

#     #We have U={g1,...,gn} and f(gi)=g^ei then f(U)=g^e1,...,g^en, with gcd(e1,...,en)=n
#     #such that f(U)=<a^n>. We store such exponents in the vector e.
#     ndx := Length(Exponents(imgs[1]));
#     e := [];
#     for j in [1..Length(imgs)] do
#         Add(e, Exponents(imgs[j])[ndx]);
#     od;

#     #We solve the equation e1s1+...+ensn=k where k is the exponent of w
#     n := Gcd ( e );
#     k := Exponents( Image( nat, w )  )[ndx];
#     e := GcdRepresentation(e)*( k / n );

#     #Calculate w=f(vAi) for some vAi
#     vAi := Identity( Image(nat) );
#     for j in [1..Length(e)] do
#         vAi := vAi*gens[j]^e[j];
#     od;

#     #Now we find v such that nat(v)= vAi            
#     v := Identity(G);
#     for j in [1..Length(Exponents(vAi))] do
#         v := v*Cgs(G)[j]^(Exponents(vAi))[j];
#     od;
    
#     return v^-1;

# end;

####################################################################
## Deprecated function to check if g and h are conjugate in G         ##
####################################################################

# InstallGlobalFunction( "IsConjugateNilGroup", function(G,g,h)

#     local   efa,    #efa series of G
#             U,      #Inverse image of the centralizer
#             v,      #Conjugating on each step
#             k,      #Conjugate on each step
#             exps,   #Storage of the exponents
#             nat,    #Natural homomorphism G->G/Gi
#             fac,    #Factor G(i-1)/Gi
#             gens,   #Generators of U/Ai
#             imgs,   #Generators of f(U/Ai)
#             f,      #Group homomorphism
#             j,i;

#     efa  := EfaSeries(G);
#     U    := G;
#     v    := Identity(G);
#     k    := h;

#     #Make a check on inputs
#     if not (IsNilpotentGroup(G) or g in G or h in G) then
#         Error( "Wrong inputs ");
#     fi;

#     #Catch trivial case
#     nat := NaturalHomomorphismByNormalSubgroup(G, efa[2]);
#     if Image(nat, k) <> Image(nat, h) then
#         return rec( conj := Identity(G), state := false );
#     fi;

#     for i in [2..Length(efa)-1] do

#         #Take G/Gi
#         nat := NaturalHomomorphismByNormalSubgroup(G, efa[i+1]);

#         #Take G(i-1)/Gi
#         fac := Image( nat, efa[i] );

#         #Create the homomorpism f:U/Gi -> G(i-1)/Gi
#         gens:= Igs( Image( nat, U ) );
#         imgs := List( gens , x -> Comm( Image(nat, g) ,x) );
#         f   := GroupHomomorphismByImages( Image( nat, U ), fac, gens, imgs);
        
#         #Check if (g^-1*k)Ai in f(U/Ai)
#         if Image( nat, g^-1*k )  = Identity( Image(nat) ) then
        
#             #In this case we don't have to do anything but as gives problems with 0/0 we have to diferenciate

#         elif Image( nat, g^-1*k ) in Image(f) then
#             v := v * PreImageConjugate(G, g^-1*k, nat, gens, imgs);
#             k := h^v;
#         else
#             return rec( conj := Identity(G), state := false );
#         fi;

#         U:= PreImage(nat, Kernel(f) );
#     od;

#     #Check the result is correct
#     if h^v <> g then
#         Error(" The function is not working properly. ");
#     fi;

#     #Return the result
#     return rec( conj := v, state := true );

# end );

####################################################################
## Deprecated function to calculate the canonical conjugate of g in G        ##
####################################################################


# InstallGlobalFunction( "CanonicalConjugate", function(G, g)

#     local   efa,    #efa series of G
#             U,      #Inverse image of the centralizer
#             v,      #Conjugating on each step
#             h,      #Canonical conjugate on each step
#             nat,    #Natural homomorphism G->G/Gi
#             fac,    #Factor G(i-1)/Gi
#             gens,   #Generators of U/Ai
#             imgs,   #Generators of f(U/Ai)
#             f,      #Group homomorphism
#             ndx,
#             k,
#             al,
#             a,b,i,j;

#     efa  := EfaSeries(G);
#     U    := G;
#     nat  := NaturalHomomorphismByNormalSubgroup(G, efa[2]);
#     k    := Igs(G)[1]^Exponents( nat( g ) )[1];
#     h    := k;
#     v    := Identity(G);

#     for i in [2..Length(efa)-1] do

#         #Take G/Gi
#         nat := NaturalHomomorphismByNormalSubgroup(G, efa[i+1]);
#         ndx := Length(Exponents( nat(g) ));
#         #Take G(i-1)/Gi
#         fac := Image( nat, efa[i] );

#         #Create the homomorpism f:U/Gi -> G(i-1)/Gi
#         gens:= Igs( Image( nat, U ) );
#         imgs := List( gens , x -> Comm( Image(nat, h) ,x) );
#         f   := GroupHomomorphismByImages( Image( nat, U ), fac, gens, imgs);
#         if Size( Image(f) )=1 then

#             al := Cgs(G)[ndx]^Exponents( nat( h^-1*g ) )[ndx];

#         elif nat(h^-1*g) in Image(f) then

#             al := Identity(G);
#             v := v * PreImageConjugate(G, g^-1*(h*al), nat, gens, imgs);

#         else 
#             a := Exponents(( Igs( Image(f) )[1] ))[ndx];
#             b := Exponents(nat(h^-1*g))[ndx];
#             if a+b < 0 then
#                 a := a*( AbsoluteValue(Int(b/a)) + 1 );
#             elif a+b > 0 then
#                 if 0 < -a*( AbsoluteValue(Int(b/a)) ) + b and -a*( AbsoluteValue(Int(b/a)) ) < a+b then 
#                     a := -a*( AbsoluteValue(Int(b/a)) );
#                 fi;
#             fi;
#             al := Cgs(G)[ndx]^( a+b );
#             v := v * PreImageConjugate(G, g^-1*(h*al), nat, gens, imgs);

#         fi;
#         U := PreImage(nat, Kernel(f) );
#         k := k*al;
#         h := k^v;
#     od;

#     #Check the result is correct
#     if h <> g then
#         #Error(" The function is not working properly. ");
#     fi;

#     #Return the result
#     return rec( kano := k, conj := v );

# end );