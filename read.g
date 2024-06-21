dir   := DirectoryCurrent();
dir   := Filename(dir, "pkg/NilGroups");
dir   := Directory( dir );
if IsBool( DirectoryContents( dir ) ) then
    Error( "The current directory is empty." );
fi;

paths := [ "nilgrp.gd", "general.gi", "exam.gi", "inter.gi", "prod.gi", "elmcon.gi", "subgrpcon.gi"];

for path in paths do
    filename := Filename(dir, path);
    if not IsExistingFile(filename) then
        str := StringFormatted( "File {} is not in the current directory", filename);
        Error( str );
    else
        Read( filename );
    fi;
od;

Unbind( dir );
Unbind( paths ); 
Unbind( filename );
Unbind( str );