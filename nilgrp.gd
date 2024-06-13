LoadPackage( "nq" );

### General.gi
DeclareProperty( "IsTauGroup", IsGroup );
InstallTrueMethod( IsTauGroup, IsNilpotentGroup and IsTorsionFree);

DeclareGlobalFunction( "Pref" );
DeclareGlobalFunction( "TauVector" );
DeclareGlobalFunction( "ReducePcpElement" );
DeclareGlobalFunction( "RandomElementRangeGenerators" );    
DeclareGlobalFunction( "RandomSubgroup" );
DeclareGlobalFunction( "ConjugacyOrder" );
DeclareGlobalFunction( "Sifting" );

### Exam.gi
DeclareGlobalFunction( "MinimalTauGroups" );

### Series.gi
DeclareGlobalFunction( "IsolatorSeries" );
DeclareGlobalFunction( "PcpsOfIsolatorSeries" );
DeclareGlobalFunction( "TauType" );
DeclareGlobalFunction( "TwoStepCentralizer" );
DeclareGlobalFunction( "ExtendedIsolatorSeries" );
DeclareGlobalFunction( "ExtendedTauType" );

### Iso.gi
DeclareGlobalFunction( "IsomorphismTauGroup" );
DeclareGlobalFunction( "SameGenus" );

### Inter.gi

DeclareGlobalFunction( "IntersectionSeriesTerm" );
DeclareGlobalFunction( "InducedIntersectionSeries" );
DeclareGlobalFunction( "PcpsOfInducedIntersectionSeries" );
DeclareGlobalFunction( "IntersectionSubgroupsNilGroups" );

### prod.gi

DeclareGlobalFunction( "SubgroupProductPair" );
DeclareGlobalFunction( "ProductDecomposition" );

### conjugacy.gi
DeclareGlobalFunction( "CentralizerNilGroup" );
DeclareGlobalFunction( "IsConjugateNilGroup" );
DeclareGlobalFunction( "CanonicalConjugateElements" );
DeclareGlobalFunction( "IsCanonicalConjugateElements" );
DeclareGlobalFunction( "NormalizerNilGroup" );
DeclareGlobalFunction( "IsConjugateSubgroups" );
DeclareGlobalFunction( "CanonicalConjugateSubgroup");
DeclareGlobalFunction( "IsCanonicalConjugateSubgroups" );
DeclareGlobalFunction( "CanonicalConjugateList" );
DeclareGlobalFunction( "IsConjugateList" );