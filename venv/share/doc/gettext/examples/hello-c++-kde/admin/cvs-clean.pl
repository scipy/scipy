#! /usr/bin/perl

#
# This script recursively (beginning with the current directory)
# wipes out everything not registered in CVS.
#
# written by Oswald Buddenhagen <ossi@kde.org>
#  inspired by the "old" cvs-clean target from Makefile.common
#
# This file is free software in terms of the BSD licence. That means
# that you can do anything with it except removing this license or
# the above copyright notice. There is NO WARRANTY of any kind.
#

sub rmrf()
{
  my $fn = shift;
  lstat ($fn);
  if (-d _) {
    if (opendir (DIR, $fn)) {
      for my $efn (grep (!/^\.\.?$/, readdir (DIR))) {
	&rmrf ($fn."/".$efn);
      }
      closedir (DIR);
      rmdir ($fn);
    }
  } else {
    unlink ($fn);
  }
}

sub newfiles()
{
  my ($indir, $incvs) = @_;
  for my $n (keys (%$incvs)) { delete $$indir{$n} }
  return sort (keys (%$indir));
}

sub cvsclean()
{
  my $dir = shift;
  my (%dirsdir, %filesdir, %dirscvs, %filescvs);
  my $dnam = $dir ? $dir : ".";
  if (!opendir (DIR, $dnam)) {
    print STDERR "Cannot enter \"".$dnam."\".\n";
    return;
  }
  for my $fn (grep (!/^\.\.?$/, readdir (DIR))) {
    if (-d $dir.$fn) {
      $fn eq "CVS" or $dirsdir{$fn} = 1;
    } else {
      $filesdir{$fn} = 1;
    }
  }
  closedir (DIR);
  if (!open (FILE, "<".$dir."CVS/Entries")) {
    print STDERR "No CVS information in \"".$dnam."\".\n";
    return;
  }
  while (<FILE>) {
    m%^D/([^/]+)/.*$% and $dirscvs{$1} = 1;
    m%^/([^/]+)/.*$% and $filescvs{$1} = 1;
  }
  close (FILE);
  if (open (FILE, "<".$dir."CVS/Entries.Log")) {
    while (<FILE>) {
      m%^A D/([^/]+)/.*$% and $dirscvs{$1} = 1;
      m%^A /([^/]+)/.*$% and $filescvs{$1} = 1;
      m%^R D/([^/]+)/.*$% and delete $dirscvs{$1};
      m%^R /([^/]+)/.*$% and delete $filescvs{$1};
    }
    close (FILE);
  }
  for my $fn (&newfiles (\%filesdir, \%filescvs)) {
    print ("F ".$dir.$fn."\n");
    &rmrf ($dir.$fn);
  }
  for my $fn (&newfiles (\%dirsdir, \%dirscvs)) {
    print ("D ".$dir.$fn."\n");
    &rmrf ($dir.$fn);
  }
  for my $fn (sort (keys (%dirscvs))) {
    &cvsclean ($dir.$fn."/");
  }
}

&cvsclean ("");
