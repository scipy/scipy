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
use Tie::File;
use Fcntl 'O_RDONLY';

sub fatal($);
sub parse_copyright($);
sub parselicense($);
sub remove_comments($);

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

my ($opt_verbose, $opt_lines, $opt_noconf) = ('', '', '');
my $opt_ignore_regex = $default_ignore_regex;
my $opt_check_regex = $default_check_regex;
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

push @find_args, qw(-not ( -path */LayoutTests/* -prune ) );
push @find_args, qw(-not ( -path */out/Debug/* -prune ) );
push @find_args, qw(-not ( -path */out/Release/* -prune ) );
push @find_args, qw(-not ( -path .git* -prune ) );
push @find_args, qw(-not ( -path .svn* -prune ) );

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
    my $header = '';
    my $copyright_match;
    my $copyright = '';
    my $license = '';
    my %copyrights;

    open (F, "<$file") or die "Unable to access $file\n";
    while (<F>) {
        last if ($. > $opt_lines);
        $header .= $_;
    }
    close(F);

    $copyright = join(" / ", values %copyrights);

    print qq(----- $file header -----\n$header----- end header -----\n\n)
	if $opt_verbose;

    remove_comments($header);
    $license = parselicense($header);

    # If no license in header, check footer (slow, because read file backwards)
    # Need for instance for Perl files, which often use the footer
    if ($license eq "UNKNOWN") {
        my $footer = '';
        tie(my @file_lines, "Tie::File", $file, autochomp => 0, mode => O_RDONLY) or die("Unable to access $file\n");
        # Avoid indexing error if header is entire file
        if ($#file_lines >= $opt_lines) {
            foreach (@file_lines[-$opt_lines .. -1]) {
                $footer .= $_;
            }
        }
        print qq(----- $file footer -----\n$header----- end footer -----\n\n)
            if $opt_verbose;
        remove_comments($footer);
        $license = parselicense($footer);
    }

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

sub remove_comments($) {
    $_ = $_[0];
    # Remove Fortran comments
    s/^[cC] //gm;
    # Remove .ASM comments
    s#^;\*?##gm;
    # Remove .S comments
    s#^@ ##gm;
    # Remove new lines
    tr/\t\r\n/ /;
    # Remove C / C++ comments
    s#(\*/|/[/*])##g;
    # Remove all characters not matching search
    tr% A-Za-z.,@;0-9\(\)/-%%cd;
    # Collapse multiple spaces into single space
    tr/ //s;
    $_[0] = $_;
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
    my $lgplver = "";
    my $extrainfo = "";
    my $license = "";

    if ($licensetext =~ /version ([^, ]+?)[.,]? (?:\(?only\)?.? )?(?:of the GNU (Affero )?General Public License )?(as )?published by the Free Software Foundation/i or
	$licensetext =~ /GNU (?:Affero )?General Public License (?:as )?published by the Free Software Foundation; version ([^, ]+?)[.,]? /i or
	$licensetext =~ /GNU (?:Affero )?General Public License,? [Vv]ersion (\d+(?:\.\d+)?)[ \.]/) {
	$gplver = " (v$1)";
    } elsif ($licensetext =~ /either version ([^ ]+)(?: of the License)?, or \(at your option\) any later version/) {
	$gplver = " (v$1 or later)";
    }

    if ($licensetext =~ /version ([^, ]+?)[.,]? (?:or later|or any later version) (?:of the GNU (?:Lesser |Library )General Public License )(as )?published by the Free Software Foundation/i or
	$licensetext =~ /(?:GNU (?:Lesser |Library )|(?:Lesser|Library) GNU )General Public License (?:(?:as )?published by the Free Software Foundation;)?,? (?:either )?[Vv]ersion ([^, ]+?)(?: of the license)?[.,]? (?:or later|or (?:\(at your option\) )?any later version)/i or
	$licensetext =~ /GNU (?:Lesser |Library )General Public License(?: \(LGPL\))?,? [Vv]ersion (\d+(?:\.\d+)?)[ \.]/) {
	$lgplver = " (v$1 or later)";
    }

    if ($licensetext =~ /permission (?:is (also granted|given))? to link (the code of )?this program with (any edition of )?(Qt|the Qt library)/i) {
	$extrainfo = " (with Qt exception)$extrainfo"
    }

    if ($licensetext =~ /(All changes made in this file will be lost|DO NOT (EDIT|delete this file)|Generated (automatically|by|from)|generated.*file|automatically generated)/i) {
	$license = "GENERATED FILE";
    }

    if ($licensetext =~ /is (free software.? you can redistribute it and\/or modify it|licensed) under the terms of (version [^ ]+ of )?the (GNU (Library |Lesser )General Public License|LGPL)/i or
        $licensetext =~ /(is distributed|may be used|can redistribute).*terms.*(LGPL|(Lesser|Library) GNU General Public License)/) {
        if ($lgplver) {
	    $license = "LGPL$lgplver$extrainfo $license";
        } else {
	    $license = "LGPL (unversioned/unknown version) $license";
        }
    }

    if ($licensetext =~ /is free software.? you (can|may) redistribute it and\/or modify it under the terms of (?:version [^ ]+ (?:\(?only\)? )?of )?the GNU General Public License/i) {
	$license = "GPL$gplver$extrainfo $license";
    } elsif ($licensetext =~ /is distributed under the terms of the GNU General Public License/
	and $gplver) {
	$license = "GPL$gplver$extrainfo $license";
    } elsif ($licensetext =~ /is distributed.*terms.*[^L]GPL/) {
        if ($gplver) {
	    $license = "GPL$gplver$extrainfo $license";
        } else {
	    $license = "GPL (unversioned/unknown version) $license";
        }
    }

    if ($licensetext =~ /This file is part of the .*Qt GUI Toolkit. This file may be distributed under the terms of the Q Public License as defined/) {
	$license = "QPL (part of Qt) $license";
    } elsif ($licensetext =~ /may be distributed under the terms of the Q Public License as defined/) {
	$license = "QPL $license";
    }

    if ($licensetext =~ /opensource\.org\/licenses\/mit/) {
	$license = "MIT/X11 (BSD like) $license";
    } elsif ($licensetext =~ /Permission is hereby granted, free of charge, to any person obtaining a copy of this software and(\/or)? associated documentation files \(the (Software|Materials)\), to deal in the (Software|Materials)/) {
	$license = "MIT/X11 (BSD like) $license";
    } elsif ($licensetext =~ /Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, modify, and distribute this software and its documentation for any purpose/) {
	$license = "MIT/X11 (BSD like) $license";
    } elsif ($licensetext =~ /Permission to use, copy, modify, distribute, and sell this software and its documentation for any purpose is hereby granted without fee/) {
	$license = "MIT/X11 (BSD like) $license";
    } elsif ($licensetext  =~ /MIT .* License/) {
	$license = "MIT/X11 (BSD like) $license";
    }

    if ($licensetext  =~ /This file is part of the Independent JPEG Group(')?s software.*For conditions of distribution and use, see the accompanying README file/i) {
	$license = "Independent JPEG Group License $license";
    }

    if ($licensetext  =~ /the University of Illinois Open Source License/){
	$license = "University of Illinois/NCSA Open Source License (BSD like) $license";
    }

    if ($licensetext  =~ /Permission to use, copy, modify, and(\/or)? distribute this software (and its documentation )?for any purpose (with or )?without fee is hereby granted, provided.*(copyright|entire) notice.*all copies/i) {
	$license = "ISC $license";
    }

    if ($licensetext =~ /THIS SOFTWARE IS PROVIDED .*AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY/ ||
        $licensetext =~ /THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABIL- ITY/) {
	if ($licensetext =~ /All advertising materials mentioning features or use of this software must display the following/) {
	    $license = "BSD (4 clause) $license";
	} elsif ($licensetext =~ /be used to endorse or promote products derived from this software/) {
	    $license = "BSD (3 clause) $license";
	} elsif ($licensetext =~ /Redistributions of source code must retain the above copyright notice/) {
	    $license = "BSD (2 clause) $license";
	} else {
	    $license = "BSD $license";
	}
    } elsif ($licensetext =~ /Use of this source code is governed by a BSD-style license/) {
        $license = "BSD-like $license";
    } elsif ($licensetext =~ /BSD terms apply/) {
        $license = "BSD-like $license";
    } elsif ($licensetext =~ /subject to the BSD License/) {
        # TODO(sbc): remove this case once we fix: http://crbug.com/177268
        $license = "BSD-like $license";
    } elsif ($licensetext =~ /license BSD/) {
        $license = "BSD-like $license";
    } elsif ($licensetext =~ /GOVERNED BY A BSD-STYLE SOURCE LICENSE/) {
        $license = "BSD-like $license";
    } elsif ($licensetext =~ /BSD 2 Clause License/) {
        $license = "BSD (2 clause) $license";
    } elsif ($licensetext =~ /BSD 3-Clause license/) {
        $license = "BSD (3 clause) $license";
    }

    if ($licensetext =~ /Mozilla Public License( Version|, v.) ([^ ]+[^., ]),?/) {
	$license = "MPL (v$2) $license";
    }

    if ($licensetext =~ /Released under the terms of the Artistic License ([^ ]+)/) {
	$license = "Artistic (v$1) $license";
    }

    if ($licensetext =~ /is free software under the Artistic [Ll]icense/) {
	$license = "Artistic $license";
    }

    if ($licensetext =~ /This (program|library) is free software; you can redistribute it and\/or modify it under the same terms as Perl itself/) {
	$license = "Perl $license";
    }

    if ($licensetext =~ /under the terms of the Apache ([^ ]+) License OR version 2 of the GNU/) {
	$license = "Apache (v$1) GPL (v2) $license";
    } elsif ($licensetext =~ /under the Apache License, Version ([^ ]+)/) {
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

    if ($licensetext =~ /under the SGI Free Software (B License|License B)/) {
	$license = "SGI Free Software License B $license";
    }

    if ($licensetext =~ /(in|into) the public domain/i) {
	$license = "Public domain $license";
    }

    if ($licensetext =~ /terms of the Common Development and Distribution License(, Version ([^(]+))? \(the License\)/) {
	$license = "CDDL " . ($1 ? "(v$2) " : '') . $license;
    }

    if ($licensetext =~ /Microsoft Permissive License \(Ms-PL\)/) {
        $license = "Ms-PL $license";
    }

    if ($licensetext =~ /as defined in and that are subject to the Apple Public Source License([ ,-]+Version ([^ ]+)?(\.))/) {
	$license = "APSL " . ($1 ? "(v$2) " : '') . $license;
    } elsif ($licensetext =~ /provided that if you redistribute the Apple Software in its entirety and without modifications, you must retain this notice and the following text and disclaimers in all such redistributions of the Apple Software/) {
	# https://fedoraproject.org/wiki/Licensing/Apple_MIT_License
	$license = "Apple MIT $license";
    }

    if ($licensetext =~ /Permission is hereby granted, free of charge, to any person or organization obtaining a copy of the software and accompanying documentation covered by this license \([\"]?the Software[\"]?\)/ or
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

    if ($licensetext =~ /under MIT license/) {
        $license = "MIT/X11 (BSD like) $license";
    }

    if ($licensetext =~ /License MIT(-| )License/) {
        $license = "MIT/X11 (BSD like) $license";
    }

    if ($licensetext =~ /As a special exception, you may create a larger work that contains part or all of the Bison parser skeleton and distribute that work under terms of your choice/) {
        $license = $license . "with Bison parser exception";
    }

    if ($licensetext =~ /As a special exception to the GNU General Public License, if you distribute this file as part of a program or library that is built using GNU Libtool, you may include this file under the same distribution terms that you use for the rest of that program/) {
        $license = $license . "with libtool exception";
    }

    if ($licensetext =~ /These materials are protected by copyright laws and contain material proprietary to the Khronos Group, Inc\. You may use these materials for implementing Khronos specifications, without altering or removing any trademark, copyright or other notice from the specification/) {
        $license = $license . "Khronos Group";
    }

    if ($licensetext =~ /This file is part of the FreeType project, and may only be used(,)? modified(,)? and distributed under the terms of the FreeType project license, LICENSE\.TXT\. By continuing to use, modify, or distribute this file you indicate that you have read the license and understand and accept it fully/) {
        $license = "FreeType (BSD like) $license";
    }
    if ($licensetext =~ /This software, and all works of authorship, whether in source or object code form as indicated by the copyright notice.*is made available, and may only be used, modified, and distributed under the FreeType Project License, LICENSE\.TXT\. Additionally, subject to the terms and conditions of the FreeType Project License, each contributor to the Work hereby grants to any individual or legal entity exercising permissions granted by the FreeType Project License and this section.*a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable.*patent license to make/) {
        $license = "FreeType (BSD like) with patent clause $license";
    }

    if ($licensetext =~ /Anti-Grain Geometry.*Permission to copy, use, modify, sell and distribute this software is granted provided this copyright notice appears in all copies. This software is provided as is without express or impl/) {
        $license = "Anti-Grain Geometry $license";
    }

    if ($licensetext =~ /Developed at SunSoft, a Sun Microsystems, Inc\. business\. Permission to use, copy, modify, and distribute this software is freely granted, provided that this notice is preserved\./) {
        $license = "SunSoft (BSD like) $license";
    }

    $license = "UNKNOWN" unless $license;

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
