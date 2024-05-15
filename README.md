# Introduction

Functions to work with infinite nilpotent groups and torsion-free nilpotent groups, Tau Groups in short.
To install al the functions save the files inside the GAP directory/pkg/NilGroups and then run in GAP.
'''
Read( "*Gap Directory*/pkg/NilGroups/read.g");
'''

# Contents

- read.g:       File to read in gap.
- nilgrp.gd:    Definition of the functions 
- general.gi:   Offers some functions that are general for the other files.
- exam.gi:      Functions to generate some minimal Tau groups and the groups used in experiments.
- series.gi:    Functions to generate the isolator series of a Tau group.
- iso.gi:       Functions to check if two Tau groups are isomorphic (**CURRENTLY NON WORKING**).
- order.gi:     Functions to order elements with the order 0<<1<<...<<-1...
- conjugacy.gi: Functions to solve the following problems:
  - Centralizer problem. CentralizerNilGroup(G, elms)
  - Conjugacy problem. IsConjugateNilGroup(G, g, h), IsCanonicalConjugateNilGroup(G, elms)
  - Canonical conjugate. CanonicalConjugateNilGroup(G, elms)
  - Normalizer problem. NormalizerNilGroup(G, U)
  - Subgroup Conjugacy problem. IsConjugateSubgroups(G, U, V), IsCanonicalConjugateSubgroups(G, U, V)
  - Canonical conjugate subgroup. CanonicalConjugateSubgroup(G, U)
- inter.gi: Intersection functions.
- deprecated.gi: Functions that are no longer used but can be recovered if necessary.
- experiments.g: Script to see the performance of the algorithms.

# Cites

> - Eddie H. Lo. Finding intersections and normalizers in finitely generated nilpotent groups. J. Symbolic Comput., 25(1):45–59, 1998.
> - Charles C. Sims. Computation with finitely presented groups, volume 48 of Encyclopedia of Mathematics and its Applications. Cambridge University Press, Cambridge, 1994.
> - B. Eick and G. Ostheimer. On the orbit-stabilizer problem for integral matrix actions of polycyclic groups. Mathematics of Computation,72(243):1511–1529, 2003.
> - B. Eick and A.-K. Engel. The isomorphism problem for torsion free nilpotent groups of Hirsch length at most 5, Groups Complex. Cryptol. 9 (1) (2017) 55--75.
