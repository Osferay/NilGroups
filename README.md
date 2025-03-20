# Introduction

Functions to work with finitely presented nilpotent groups. This code is complementary to the paper [4].
To install al the functions save the files inside the GAP directory/pkg/NilGroups and then run in GAP.
```
Read( "*Gap Directory*/pkg/NilGroups/read.g");
```
**Some functions in this repository rely on the NQ package, so having it installed is mandatory.**

# Contents

### read.g:       File to read in gap.

### nilgrp.gd:    Definition of the functions.

### general.gi:   Offers some functions that are general for the other files.

### exam.gi:      Functions to generate the groups used in experiments.
  - **SomeNilpotentGroups( n )** : Retruns the group assigned to the number n, this groups are:
    - **n = 1** : Group Torsion-free of Hirsch length 5.
    - **n = 2** : Group with relative orders [0, 0, 0, 0, 0, 4, 2, 2].
    - **n = 3** : Group with relative orders [0, 0, 0, 0, 3872].
    - **n = 4** : Group with relative orders [0, 0, 0, 352, 11264].
    - **n = 5** : Group with 85 generators and Hirsch length 23.
    - **n = 6** : Group with 72 generators and Hirsch length 37.
    - **n = 7** : Group with 20 generators and Hirsch length 6.
    - **n = 8** : Group with 21 generators and Hirsch length 10.

### elmcon.gi: Functions to solve the element conjugacy problems:
  - *Centralizer problem*
    - **CentralizerNilGroup( G, elms )** : Returns the centralizer of the elements in the list elms.
  - *Conjugacy problem* 
    - **IsConjugateNilGroup( G, g, h )** : If the elements g and h are conjugate returns the conjugating element, otherwise returns false.
    - **IsCanonicalConjugateElements( G, elms )** : If the elements of the list elms are conjugate returns the canonical element, kano, the conjugating element, conj, and the centralizer cent. Otherwise returns false.
  - *Canonical conjugate* 
    - **CanonicalConjugateElements( G, elms )** : Returns the canonical conjugate of each element in the list elms. Also returs the centralizers and the conjugating elements.
  - *Conjugacy problem for lists of elements*
    - **IsConjugateList( G, list )** : Given a list of elements in G, returns a conjugate representative for each conjugacy class found and a list of indexes which represents the elements that belong to the respective conjugacy class.
    - **CanonicalConjugateList( G, list )** : Given a list of elements in G, returns a conjugate representative for each conjugacy class found and a list of indexes which represents the elements that belong to the respective conjugacy class.

Setting the information level to 1 of InfoConjugacy will show the progres of the algoritms.

### subgrpcon.gi: Functions to solve the subgroup conjugacy problems:
  - *Normalizer problem*
    - **NormalizerNilGroup( G, U )** : Returns the normalizer of the subgroup U in the group G.
  - *Subgroup Conjugacy problem*
    - **IsConjugateSubgroups( G, U, V )** : Returns the conjugating element if the two subgroups are conjugated. Otherwise returns false.
    - **IsCanonicalConjugateSubgroups( G, U, V )** : Returns the canonical conjugate of U and V, the conjugating elements and the normalizer, if both subgroups are conjugated. Otherwise reutrns false.
  - *Canonical conjugate subgroup* 
    - **CanonicalConjugateSubgroup( G, U )** : Returns the canonical conjugate of U, the conjugating element and the normalizer.

Setting the information level to 1 of InfoConjugacy will show the progres of the algoritms.
### inter.gi: Intersection functions.
  - **IntersectionSubgroupsNilGroups( G, U, V )** : Returns a list of generators for the intersection of U and V.
### prod.gi:  Product functions.
  - **SubgroupProductPair( G, U, V )** : Returns a product pair of U and V.
  - **ProductDecomposition( G, U, V, g)** : Returns, if posible, a descomposition of g as a multiplication of elements in U and V 
### experiments.g: 
Script to see the performance of the algorithms using the groups in exam.gi. The results are avaribale in [B. Eick and O. F. Ayala].

# Cites

1. Eddie H. Lo. Finding intersections and normalizers in finitely generated nilpotent groups. J. Symbolic Comput., 25(1):45–59, 1998.
2. Charles C. Sims. Computation with finitely presented groups. Volume 48 of Encyclopedia of Mathematics and its Applications. Cambridge University Press, Cambridge, 1994.
3. B. Eick and G. Ostheimer. On the orbit-stabilizer problem for integral matrix actions of polycyclic groups. Mathematics of Computation,72(243):1511–1529, 2003.
4. B. Eick and O. Fernández Ayala. The conjugacy problem and canonical representatives in finitely generated nilpotent groups. J. Symbolic Comput. 130:102422, 2025.
