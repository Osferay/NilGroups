### General.gi
DeclareProperty( "IsTauGroup", IsGroup );
InstallTrueMethod( IsTauGroup, IsNilpotentGroup and IsTorsionFree);

DeclareGlobalFunction( "Pref" );
DeclareGlobalFunction( "TauVector" );
DeclareGlobalFunction( "MatrixRelationsByVector" );
DeclareGlobalFunction( "TauGroupByVector" );

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

### conjugacy.gi
DeclareGlobalFunction( "CentralizerNilGroup" );
DeclareGlobalFunction( "IsConjugateNilGroup" );
DeclareGlobalFunction( "CanonicalConjugateNilGroup" );
DeclareGlobalFunction( "NormalizerNilGroup" );
DeclareGlobalFunction( "IsConjugateSubgroupsNilGroup" );

Reread( "/home/oscar/T-Groups/general.gi" );
Reread( "/home/oscar/T-Groups/series.gi" );
Reread( "/home/oscar/T-Groups/exam.gi" );
#Reread( "/home/oscar/T-Groups/iso.gi" );
Reread( "/home/oscar/T-Groups/conjugacy.gi" );