#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $usage =
"$0 file_fasta

Parses RNA sequences to DNA by substituting 'U's with 'T'.
";

my $file_fasta=shift or die $usage;


parse_fasta($file_fasta);

exit;




sub parse_fasta{

    my ($file) = @_;
    my $is_q = 0;
    open (FASTA, "<$file") or die "can not open $file\n";
    while (<FASTA>){
      
	if (/^(>\S+)/){
	
	    print "$1\n";

    }elsif(/^@(\S+)/){

        print ">$1\n";
        $is_q = 1;

    }elsif(/^[AUTGCN]+$/){
		
	    s/U/T/g;
	    print uc($_);
        if ($is_q == 1){
            <FASTA>;
            <FASTA>;
        }
	}
    }
    
    close FASTA;
    return;
}
