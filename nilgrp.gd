LoadPackage( "nq" );

### General.gi
DeclareProperty( "IsTauGroup", IsGroup );
InstallTrueMethod( IsTauGroup, IsNilpotentGroup and IsTorsionFree);

DeclareGlobalFunction( "Pref" );
DeclareGlobalFunction( "TauVector" );
DeclareGlobalFunction( "MatrixRelationsByVector" );
DeclareGlobalFunction( "TauGroupByVector" );
DeclareGlobalFunction( "RandomElementRangeGenerators" );    
DeclareGlobalFunction( "RandomSubgroup" );
DeclareGlobalFunction( "ConjugacyOrder" );
DeclareGlobalFunction( "Shifting" );

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

### conjugacy.gi
DeclareGlobalFunction( "CentralizerNilGroup" );
DeclareGlobalFunction( "IsConjugateNilGroup" );
DeclareGlobalFunction( "CanonicalConjugateNilGroup" );
DeclareGlobalFunction( "IsCanonicalConjugateNilGroup" );
DeclareGlobalFunction( "NormalizerNilGroup" );
DeclareGlobalFunction( "IsConjugateSubgroups" );
DeclareGlobalFunction( "CanonicalConjugateSubgroup");
DeclareGlobalFunction( "IsCanonicalConjugateSubgroups" );

Reread( "/home/oscar/NilGroups/general.gi" );
Reread( "/home/oscar/NilGroups/series.gi" );
Reread( "/home/oscar/NilGroups/exam.gi" );
#Reread( "/home/oscar/NilGroups/iso.gi" );
# Reread( "/home/oscar/NilGroups/inter.gi" );
Reread( "/home/oscar/NilGroups/conjugacy.gi" );