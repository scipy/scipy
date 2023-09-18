#!/usr/bin/perl -w
# This script was originally based on the script of the same name from
# the KDE SDK (by dfaure@kde.org)
#
# This version is
#   Copyright (C) 2007, 2008 Adam D. Barratt
#   Copyright (C) 2012 Francesco Poli
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

=head1 NAME

licensecheck - simple license checker for source files

=head1 SYNOPSIS

B<licensecheck> B<--help>|B<--version>

B<licensecheck> [B<--no-conf>] [B<--verbose>] [B<--copyright>]
[B<-l>|B<--lines=>I<N>] [B<-i>|B<--ignore=>I<regex>] [B<-c>|B<--check=>I<regex>]
[B<-m>|B<--machine>] [B<-r>|B<--recursive>]
I<list of files and directories to check>

=head1 DESCRIPTION

B<licensecheck> attempts to determine the license that applies to each file
passed to it, by searching the start of the file for text belonging to
various licenses.

If any of the arguments passed are directories, B<licensecheck> will add
the files contained within to the list of files to process.

=head1 OPTIONS

=over 4

=item B<--verbose>, B<--no-verbose>

Specify whether to output the text being processed from each file before
the corresponding license information.

Default is to be quiet.

=item B<-l=>I<N>, B<--lines=>I<N>

Specify the number of lines of each file's header which should be parsed
for license information. (Default is 60).

=item B<-i=>I<regex>, B<--ignore=>I<regex>

When processing the list of files and directories, the regular
expression specified by this option will be used to indicate those which
should not be considered (e.g. backup files, VCS metadata).

=item B<-r>, B<--recursive>

Specify that the contents of directories should be added
recursively.

=item B<-c=>I<regex>, B<--check=>I<regex>

Specify a pattern against which filenames will be matched in order to
decide which files to check the license of.

The default includes common source files.

=item B<--copyright>

Also display copyright text found within the file

=item B<-m>, B<--machine>

Display the information in a machine readable way, i.e. in the form
<file><tab><license>[<tab><copyright>] so that it can be easily sorted
and/or filtered, e.g. with the B<awk> and B<sort> commands.
Note that using the B<--verbose> option will kill the readability.

=item B<--no-conf>, B<--noconf>

Do not read any configuration files. This can only be used as the first
option given on the command-line.

=back

=head1 CONFIGURATION VARIABLES

The two configuration files F</etc/devscripts.conf> and
F<~/.devscripts> are sourced by a shell in that order to set
configuration variables.  Command line options can be used to override
configuration file settings.  Environment variable settings are
ignored for this purpose.  The currently recognised variables are:

=over 4

=item B<LICENSECHECK_VERBOSE>

If this is set to I<yes>, then it is the same as the B<--verbose> command
line parameter being used. The default is I<no>.

=item B<LICENSECHECK_PARSELINES>

If this is set to a positive number then the specified number of lines
at the start of each file will be read whilst attempting to determine
the license(s) in use.  This is equivalent to the B<--lines> command line
option.

=back

=head1 LICENSE

This code is copyright by Adam D. Barratt <I<adam@adam-barratt.org.uk>>,
all rights reserved; based on a script of the same name from the KDE
SDK, which is copyright by <I<dfaure@kde.org>>.
This program comes with ABSOLUTELY NO WARRANTY.
You are free to redistribute this code under the terms of the GNU
General Public License, version 2 or later.

=head1 AUTHOR

Adam D. Barratt <adam@adam-barratt.org.uk>

=cut

use strict;
use warnings;
use Getopt::Long qw(:config gnu_getopt);
use File::Basename;

sub fatal($);
sub parse_copyright($);
sub parselicense($);

my $progname = basename($0);

# From dpkg-source
my $default_ignore_regex = '
# Ignore general backup files
(?:^|/).*~$|
# Ignore emacs recovery files
(?:^|/)\.#.*$|
# Ignore vi swap files
(?:^|/)\..*\.swp$|
# Ignore baz-style junk files or directories
(?:^|/),,.*(?:$|/.*$)|
# File-names that should be ignored (never directories)
(?:^|/)(?:DEADJOE|\.cvsignore|\.arch-inventory|\.bzrignore|\.gitignore)$|
# File or directory names that should be ignored
(?:^|/)(?:CVS|RCS|\.deps|\{arch\}|\.arch-ids|\.svn|\.hg|_darcs|\.git|
\.shelf|_MTN|\.bzr(?:\.backup|tags)?)(?:$|/.*$)
';

# Take out comments and newlines
$default_ignore_regex =~ s/^#.*$//mg;
$default_ignore_regex =~ s/\n//sg;

my $default_check_regex = '\.(c(c|pp|xx)?|h(h|pp|xx)?|f(77|90)?|p(l|m)|xs|sh|php|py(|x)|rb|java|vala|el|sc(i|e)|cs|pas|inc|dtd|xsl|mod|m|tex|mli?)$';

my $modified_conf_msg;

my ($opt_verbose, $opt_lines, $opt_noconf, $opt_ignore_regex, $opt_check_regex)
  = ('', '', '', '', '');
my $opt_recursive = 0;
my $opt_copyright = 0;
my $opt_machine = 0;
my ($opt_help, $opt_version);
my $def_lines = 60;

# Read configuration files and then command line
# This is boilerplate

if (@ARGV and $ARGV[0] =~ /^--no-?conf$/) {
    $modified_conf_msg = "  (no configuration files read)";
    shift;
} else {
    my @config_files = ('/etc/devscripts.conf', '~/.devscripts');
    my %config_vars = (
		       'LICENSECHECK_VERBOSE' => 'no',
		       'LICENSECHECK_PARSELINES' => $def_lines,
		      );
    my %config_default = %config_vars;

    my $shell_cmd;
    # Set defaults
    foreach my $var (keys %config_vars) {
	$shell_cmd .= qq[$var="$config_vars{$var}";\n];
    }
    $shell_cmd .= 'for file in ' . join(" ", @config_files) . "; do\n";
    $shell_cmd .= '[ -f $file ] && . $file; done;' . "\n";
    # Read back values
    foreach my $var (keys %config_vars) { $shell_cmd .= "echo \$$var;\n" }
    my $shell_out = `/bin/bash -c '$shell_cmd'`;
    @config_vars{keys %config_vars} = split /\n/, $shell_out, -1;

    # Check validity
    $config_vars{'LICENSECHECK_VERBOSE'} =~ /^(yes|no)$/
	or $config_vars{'LICENSECHECK_VERBOSE'} = 'no';
    $config_vars{'LICENSECHECK_PARSELINES'} =~ /^[1-9][0-9]*$/
	or $config_vars{'LICENSECHECK_PARSELINES'} = $def_lines;

    foreach my $var (sort keys %config_vars) {
	if ($config_vars{$var} ne $config_default{$var}) {
	    $modified_conf_msg .= "  $var=$config_vars{$var}\n";
	}
    }
    $modified_conf_msg ||= "  (none)\n";
    chomp $modified_conf_msg;

    $opt_verbose = $config_vars{'LICENSECHECK_VERBOSE'} eq 'yes' ? 1 : 0;
    $opt_lines = $config_vars{'LICENSECHECK_PARSELINES'};
}

GetOptions("help|h" => \$opt_help,
	   "version|v" => \$opt_version,
	   "verbose!" => \$opt_verbose,
	   "lines|l=i" => \$opt_lines,
	   "ignore|i=s" => \$opt_ignore_regex,
	   "recursive|r" => \$opt_recursive,
	   "check|c=s" => \$opt_check_regex,
	   "copyright" => \$opt_copyright,
	   "machine|m" => \$opt_machine,
	   "noconf" => \$opt_noconf,
	   "no-conf" => \$opt_noconf,
	   )
    or die "Usage: $progname [options] filelist\nRun $progname --help for more details\n";

$opt_lines = $def_lines if $opt_lines !~ /^[1-9][0-9]*$/;
$opt_ignore_regex = $default_ignore_regex if ! length $opt_ignore_regex;
$opt_check_regex = $default_check_regex if ! length $opt_check_regex;

if ($opt_noconf) {
    fatal "--no-conf is only acceptable as the first command-line option!";
}
if ($opt_help) { help(); exit 0; }
if ($opt_version) { version(); exit 0; }

die "Usage: $progname [options] filelist\nRun $progname --help for more details\n" unless @ARGV;

$opt_lines = $def_lines if not defined $opt_lines;

my @files = ();
my @find_args = ();
my $files_count = @ARGV;

push @find_args, qw(-maxdepth 1) unless $opt_recursive;
push @find_args, qw(-follow -type f -print);

while (@ARGV) {
    my $file = shift @ARGV;

    if (-d $file) {
	open FIND, '-|', 'find', $file, @find_args
	    or die "$progname: couldn't exec find: $!\n";

	while (<FIND>) {
	    chomp;
	    next unless m%$opt_check_regex%;
	    # Skip empty files
	    next if (-z $_);
	    push @files, $_ unless m%$opt_ignore_regex%;
	}
	close FIND;
    } else {
	next unless ($files_count == 1) or $file =~ m%$opt_check_regex%;
	push @files, $file unless $file =~ m%$opt_ignore_regex%;
    }
}

while (@files) {
    my $file = shift @files;
    my $content = '';
    my $copyright_match;
    my $copyright = '';
    my $license = '';
    my %copyrights;

    open (F, "<$file") or die "Unable to access $file\n";
    while (<F>) {
        last if ($. > $opt_lines);
        $content .= $_;
	$copyright_match = parse_copyright($_);
	if ($copyright_match) {
	    $copyrights{lc("$copyright_match")} = "$copyright_match";
	}
    }
    close(F);

    $copyright = join(" / ", values %copyrights);

    print qq(----- $file header -----\n$content----- end header -----\n\n)
	if $opt_verbose;

    # Remove Fortran comments
    $content =~ s/^[cC] //gm;
    $content =~ tr/\t\r\n/ /;
    # Remove C / C++ comments
    $content =~ s#(\*/|/[/*])##g;
    $content =~ tr% A-Za-z.,@;0-9\(\)/-%%cd;
    $content =~ tr/ //s;

    $license = parselicense($content);
    if ($opt_machine) {
	print "$file\t$license";
	print "\t" . ($copyright or "*No copyright*") if $opt_copyright;
	print "\n";
    } else {
	print "$file: ";
	print "*No copyright* " unless $copyright;
	print $license . "\n";
	print "  [Copyright: " . $copyright . "]\n"
	  if $copyright and $opt_copyright;
	print "\n" if $opt_copyright;
    }
}

sub parse_copyright($) {
    my $copyright = '';
    my $match;

    my $copyright_indicator_regex = '
	(?:copyright	# The full word
	|copr\.		# Legally-valid abbreviation
	|\x{00a9}	# Unicode character COPYRIGHT SIGN
	|\xc2\xa9	# Unicode copyright sign encoded in iso8859
	|\(c\)		# Legally-null representation of sign
	)';
    my $copyright_disindicator_regex = '
	\b(?:info(?:rmation)?	# Discussing copyright information
	|notice			# Discussing the notice
	|and|or                 # Part of a sentence
	)\b';

    if (m%$copyright_indicator_regex(?::\s*|\s+)(\S.*)$%ix) {
	$match = $1;

	# Ignore lines matching "see foo for copyright information" etc.
	if ($match !~ m%^\s*$copyright_disindicator_regex%ix) {
	    # De-cruft
	    $match =~ s/([,.])?\s*$//;
	    $match =~ s/$copyright_indicator_regex//igx;
	    $match =~ s/^\s+//;
	    $match =~ s/\s{2,}/ /g;
	    $match =~ s/\\@/@/g;
	    $copyright = $match;
	}
    }

    return $copyright;
}

sub help {
   print <<"EOF";
Usage: $progname [options] filename [filename ...]
Valid options are:
   --help, -h             Display this message
   --version, -v          Display version and copyright info
   --no-conf, --noconf    Don't read devscripts config files; must be
                          the first option given
   --verbose              Display the header of each file before its
                            license information
   --lines, -l            Specify how many lines of the file header
                            should be parsed for license information
                            (Default: $def_lines)
   --check, -c            Specify a pattern indicating which files should
                             be checked
                             (Default: '$default_check_regex')
   --machine, -m          Display in a machine readable way (good for awk)
   --recursive, -r        Add the contents of directories recursively
   --copyright            Also display the file's copyright
   --ignore, -i           Specify that files / directories matching the
                            regular expression should be ignored when
                            checking files
                            (Default: '$default_ignore_regex')

Default settings modified by devscripts configuration files:
$modified_conf_msg
EOF
}

sub version {
    print <<"EOF";
This is $progname, from the Debian devscripts package, version ###VERSION###
Copyright (C) 2007, 2008 by Adam D. Barratt <adam\@adam-barratt.org.uk>; based
on a script of the same name from the KDE SDK by <dfaure\@kde.org>.

This program comes with ABSOLUTELY NO WARRANTY.
You are free to redistribute this code under the terms of the
GNU General Public License, version 2, or (at your option) any
later version.
EOF
}

sub parselicense($) {
    my ($licensetext) = @_;

    my $gplver = "";
    my $extrainfo = "";
    my $license = "";

    if ($licensetext =~ /version ([^, ]+?)[.,]? (?:\(?only\)?.? )?(?:of the GNU (Affero )?(Lesser |Library )?General Public License )?(as )?published by the Free Software Foundation/i or
	$licensetext =~ /GNU (?:Affero )?(?:Lesser |Library )?General Public License (?:as )?published by the Free Software Foundation; version ([^, ]+?)[.,]? /i) {

	$gplver = " (v$1)";
    } elsif ($licensetext =~ /GNU (?:Affero )?(?:Lesser |Library )?General Public License, version (\d+(?:\.\d+)?)[ \.]/) {
	$gplver = " (v$1)";
    } elsif ($licensetext =~ /either version ([^ ]+)(?: of the License)?, or \(at your option\) any later version/) {
	$gplver = " (v$1 or later)";
    }

    if ($licensetext =~ /(?:675 Mass Ave|59 Temple Place|51 Franklin Steet|02139|02111-1307)/i) {
	$extrainfo = " (with incorrect FSF address)$extrainfo";
    }

    if ($licensetext =~ /permission (?:is (also granted|given))? to link (the code of )?this program with (any edition of )?(Qt|the Qt library)/i) {
	$extrainfo = " (with Qt exception)$extrainfo"
    }

    if ($licensetext =~ /(All changes made in this file will be lost|DO NOT (EDIT|delete this file)|Generated (automatically|by|from)|generated.*file)/i) {
	$license = "GENERATED FILE";
    }

    if ($licensetext =~ /is (free software.? you can redistribute it and\/or modify it|licensed) under the terms of (version [^ ]+ of )?the (GNU (Library |Lesser )General Public License|LGPL)/i) {
	$license = "LGPL$gplver$extrainfo $license";
    }

    if ($licensetext =~ /is free software.? you can redistribute it and\/or modify it under the terms of the (GNU Affero General Public License|AGPL)/i) {
	$license = "AGPL$gplver$extrainfo $license";
    }

    if ($licensetext =~ /is free software.? you (can|may) redistribute it and\/or modify it under the terms of (?:version [^ ]+ (?:\(?only\)? )?of )?the GNU General Public License/i) {
	$license = "GPL$gplver$extrainfo $license";
    }

    if ($licensetext =~ /is distributed under the terms of the GNU General Public License,/
	and length $gplver) {
	$license = "GPL$gplver$extrainfo $license";
    }

    if ($licensetext =~ /is distributed.*terms.*GPL/) {
	$license = "GPL (unversioned/unknown version) $license";
    }

    if ($licensetext =~ /This file is part of the .*Qt GUI Toolkit. This file may be distributed under the terms of the Q Public License as defined/) {
	$license = "QPL (part of Qt) $license";
    } elsif ($licensetext =~ /may be distributed under the terms of the Q Public License as defined/) {
	$license = "QPL $license";
    }

    if ($licensetext =~ /opensource\.org\/licenses\/mit-license\.php/) {
	$license = "MIT/X11 (BSD like) $license";
    } elsif ($licensetext =~ /Permission is hereby granted, free of charge, to any person obtaining a copy of this software and(\/or)? associated documentation files \(the (Software|Materials)\), to deal in the (Software|Materials)/) {
	$license = "MIT/X11 (BSD like) $license";
    } elsif ($licensetext =~ /Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, modify, and distribute this software and its documentation for any purpose/) {
	$license = "MIT/X11 (BSD like) $license";
    }

    if ($licensetext  =~ /Permission to use, copy, modify, and(\/or)? distribute this software for any purpose with or without fee is hereby granted, provided.*copyright notice.*permission notice.*all copies/) {
	$license = "ISC $license";
    }

    if ($licensetext =~ /THIS SOFTWARE IS PROVIDED .*AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY/) {
	if ($licensetext =~ /All advertising materials mentioning features or use of this software must display the following acknowledge?ment.*This product includes software developed by/i) {
	    $license = "BSD (4 clause) $license";
	} elsif ($licensetext =~ /(The name .*? may not|Neither the names? .*? nor the names of (its|their) contributors may) be used to endorse or promote products derived from this software/i) {
	    $license = "BSD (3 clause) $license";
	} elsif ($licensetext =~ /Redistributions of source code must retain the above copyright notice/i) {
	    $license = "BSD (2 clause) $license";
	} else {
	    $license = "BSD $license";
	}
    }

    if ($licensetext =~ /Mozilla Public License Version ([^ ]+)/) {
	$license = "MPL (v$1) $license";
    }

    if ($licensetext =~ /Released under the terms of the Artistic License ([^ ]+)/) {
	$license = "Artistic (v$1) $license";
    }

    if ($licensetext =~ /is free software under the Artistic [Ll]icense/) {
	$license = "Artistic $license";
    }

    if ($licensetext =~ /This program is free software; you can redistribute it and\/or modify it under the same terms as Perl itself/) {
	$license = "Perl $license";
    }

    if ($licensetext =~ /under the Apache License, Version ([^ ]+)/) {
	$license = "Apache (v$1) $license";
    }

    if ($licensetext =~ /(THE BEER-WARE LICENSE)/i) {
	$license = "Beerware $license";
    }

    if ($licensetext =~ /This source file is subject to version ([^ ]+) of the PHP license/) {
	$license = "PHP (v$1) $license";
    }

    if ($licensetext =~ /under the terms of the CeCILL /) {
	$license = "CeCILL $license";
    }

    if ($licensetext =~ /under the terms of the CeCILL-([^ ]+) /) {
	$license = "CeCILL-$1 $license";
    }

    if ($licensetext =~ /under the SGI Free Software License B/) {
	$license = "SGI Free Software License B $license";
    }

    if ($licensetext =~ /is in the public domain/i) {
	$license = "Public domain $license";
    }

    if ($licensetext =~ /terms of the Common Development and Distribution License(, Version ([^(]+))? \(the License\)/) {
	$license = "CDDL " . ($1 ? "(v$2) " : '') . $license;
    }

    if ($licensetext =~ /Microsoft Permissive License \(Ms-PL\)/) {
        $license = "Ms-PL $license";
    }

    if ($licensetext =~ /Permission is hereby granted, free of charge, to any person or organization obtaining a copy of the software and accompanying documentation covered by this license \(the \"Software\"\)/ or
	$licensetext =~ /Boost Software License([ ,-]+Version ([^ ]+)?(\.))/i) {
	$license = "BSL " . ($1 ? "(v$2) " : '') . $license;
    }

    if ($licensetext =~ /PYTHON SOFTWARE FOUNDATION LICENSE (VERSION ([^ ]+))/i) {
	$license = "PSF " . ($1 ? "(v$2) " : '') . $license;
    }

    if ($licensetext =~ /The origin of this software must not be misrepresented.*Altered source versions must be plainly marked as such.*This notice may not be removed or altered from any source distribution/ or
        $licensetext =~ /see copyright notice in zlib\.h/) {
	$license = "zlib/libpng $license";
    } elsif ($licensetext =~ /This code is released under the libpng license/) {
        $license = "libpng $license";
    }

    if ($licensetext =~ /Do What The Fuck You Want To Public License, Version ([^, ]+)/i) {
        $license = "WTFPL (v$1) $license";
    }

    if ($licensetext =~ /Do what The Fuck You Want To Public License/i) {
        $license = "WTFPL $license";
    }

    if ($licensetext =~ /(License WTFPL|Under (the|a) WTFPL)/i) {
        $license = "WTFPL $license";
    }

    $license = "UNKNOWN" if (!length($license));

    # Remove trailing spaces.
    $license =~ s/\s+$//;

    return $license;
}

sub fatal($) {
    my ($pack,$file,$line);
    ($pack,$file,$line) = caller();
    (my $msg = "$progname: fatal error at line $line:\n@_\n") =~ tr/\0//d;
    $msg =~ s/\n\n$/\n/;
    die $msg;
}
