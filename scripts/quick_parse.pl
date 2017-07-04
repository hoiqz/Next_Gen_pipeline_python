#!/usr/bin/perl

use strict;

my $file=$ARGV[0];
open FILE,$file;
my $accepted_file="$file.accepted.csv";
my $rejected_file="$file.rejected.csv";
open ACC,">$accepted_file";
open REJ,">$rejected_file";
my $header=0;
my %pid_arr;
while(<FILE>)
{
	chomp;
	my $full_line=$_;
#header
#	if($header==0)
#	{
#		my $header_line=$_;
#		print "HEADER line :$header_line\n";
#		print ACC "$header_line\n";
#		print REJ "$header_line\n";
#		$header++;
#		next;
#	}
	my @line_arr=split(',',$_);
	my $pid=$line_arr[1];
	my $specialty_grp=$line_arr[13];
	if ($specialty_grp =~/.ndocrinology/)
	{
		$pid_arr{$pid}=1;
#	print "$full_line\n";
	}

}
close FILE;
open FILE,$file;
while(<FILE>)
{
	chomp;
	my @line_arr=split(',',$_);
	my $pid=$line_arr[1];
	my $spec=$line_arr[13];
#now we check those patients that have went to endocrine lab before
#if they have attent only endocrine then ignore them
	if (exists $pid_arr{$pid} && $spec=~/.iabetes/){
#		print REJ "$_\n";
#		delete $pid_arr{$pid}; # delete this from endocrine list cuase he has been to diabetes clinic
		delete $pid_arr{$pid};

	}
}
close FILE;
open FILE, $file;
while(<FILE>)
{
	chomp;
	my @line_arr=split(',',$_);
	my $pid=$line_arr[1];
	my $spec=$line_arr[13];
	if (exists $pid_arr{$pid})
	{
		print REJ "$_\n";
	}
	else
	{
		print ACC "$_\n";
	}
}

close FILE;
close ACC;
close REJ;
