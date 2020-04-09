#!/usr/bin/perl
use File::Copy;
if ($#ARGV < 0) {die "No files to work on\n";}
$infile='';
$CopyrightHeaderLine0 = "/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */";
$CopyrightHeaderLine1 = "/*                                                                       */";
$CopyrightHeaderLine2 = "/*    This file is part of the HiGHS linear optimization suite           */";
$CopyrightHeaderLine3 = "/*                                                                       */";
$CopyrightHeaderLine4 = "/*    Written and engineered 2008-2020 at the University of Edinburgh    */";
$CopyrightHeaderLine5 = "/*                                                                       */";
$CopyrightHeaderLine6 = "/*    Available as open-source under the MIT License                     */";
$CopyrightHeaderLine7 = "/*                                                                       */";
$RemoveCopyrightHeader = 0;
while(<>) {
    if ($infile ne $ARGV) {
#
# New file, so move any previously processed file now
#
	if ($infile ne '') {
	    move($outfile, $infile);
	}
#
# New file - first line
#	
	$infile=$ARGV;
	$outfile="$infile.temp";
        open(outfile,">$outfile") || die "Can't open $outfile: $! \n";
#
# Print the copyright header
#
	print(outfile "$CopyrightHeaderLine0\n");
	print(outfile "$CopyrightHeaderLine1\n");
	print(outfile "$CopyrightHeaderLine2\n");
	print(outfile "$CopyrightHeaderLine3\n");
	print(outfile "$CopyrightHeaderLine4\n");
	print(outfile "$CopyrightHeaderLine5\n");
	print(outfile "$CopyrightHeaderLine6\n");
	print(outfile "$CopyrightHeaderLine7\n");
	print(outfile "$CopyrightHeaderLine0\n");
#
# Determine whether an old copyright header should be ignored
#
	$fpos=index($_,$CopyrightHeaderLine0);
#	print "Copyright header line 0 is at $fpos\n";
	if ($fpos<0) {
	    $RemoveCopyrightHeader = 0;
	    print "Insert  copyright header in $infile\n";
	    print (outfile);
	} else {
	    $RemoveCopyrightHeader = 1;
	    print "Replace copyright header in $infile\n";
	}
#	print "Line0: Rm=$RemoveCopyrightHeader|$_";
    } else {
#	print "Line+: Rm=$RemoveCopyrightHeader|$_";
#	
# Current file - second or subsequent line
#	
	if ($RemoveCopyrightHeader>0) {
#
# Removing an old copyright header 
# Look to see if we've reached the end
#
	    $fpos=index($_,$CopyrightHeaderLine0);
	    if ($fpos>=0) {
#
# Reached the end of the old copyright header so stop ignoring lines
#
		$RemoveCopyrightHeader = 0;
	    }
	} else {
#
# Not removing an old copyright header so write out the line
#
	    print (outfile);
	}
    }
    if (eof()) {
#
# End of last file so move now
#
	move($outfile, $infile);
    }
}
