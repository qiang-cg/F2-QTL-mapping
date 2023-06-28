#!/usr/bin/perl -w

use strict;

my $in   = shift;
my $out  = shift;

open (R, $in ) || die "Reading $in failed\n";
#open (W1, ">".$out."\.p") || die "Open $out \.p for writing failed\n";
open (W2, ">".$out) || die "Open $out for writing failed\n";
 
my $flag=0;
while(<R>)
     {if($flag==1)
        { chomp $_;
          my @tm=split(/\t/, $_); 
          my $locus    = $tm[0]."_".$tm[1];
          my $chro     = $tm[0];
          $chro =~ s/Chr//;
          my $position = $tm[1];
          my $ml=1;
          if ($tm[4]=~ /\,/) {$ml=0;}         
 
          my $p1s=$tm[433];   #N
          my @tn1=split(/:/,$p1s);
          if($tn1[0] =~ /\./){next;}
          my @p1go=();
          if($tn1[0] =~ /\|/) { @p1go=split(/\|/,$tn1[0]); }
          else                { @p1go=split(/\//,$tn1[0]); }
                 
          my $p1ho=1;
          if($p1go[0] ne $p1go[1]){$p1ho=0;}
          my $p1dp=0;
          if($p1go[0] eq "."){$p1ho=0;}
          else               {$p1dp=$tn1[2];} 
          
          my $p2s=$tm[434];     #R
          my @tn2=split(/:/,$p2s);
          if($tn2[0] =~ /\./){next;}
          my @p2go=();
          if($tn2[0] =~ /\|/) { @p2go=split(/\|/,$tn2[0]); }
          else                { @p2go=split(/\//,$tn2[0]); }

          my $p2ho=1;
          if($p2go[0] ne $p2go[1]){$p2ho=0;}
          my $p2dp=0;
          if($p2go[0] eq "."){$p2ho=0;}
          else               {$p2dp=$tn2[2];}
          
          my @isg=split(/;/,$tm[7]);
          my @nf=split(/=/,$isg[2]);
          my $sn=$nf[1];
          if($sn<800) {next;}   #AN>800
          if (!$ml)   {next;}
          if($p1dp < 10 ) {next;}   #DP>=10
          if($p2dp < 10 ) {next;}
          if(!$p1ho || !$p2ho) {next;}
          if($p1go[0] eq $p2go[0]) {next;}
          my $a=$p1go[0];
          my $b=$p2go[0];
          my @geno_array=();
          my $tmn=scalar(@tm);
             $tmn--;
          my @snp=@tm[9..434];
          my $kn=0;
          foreach(@snp)
                 {my @tk=split(/:/,$_);
                 
                  my @tp=();
                  if($tk[0] =~ /\|/) {  @tp = split(/\|/,$tk[0]); }
                  else                { @tp = split(/\//,$tk[0]); }
                  
                  $geno_array[$kn]="-";
                  if($tp[0] eq $a && $tp[1] eq $b)  {$geno_array[$kn]="h";}           	
     	          if($tp[0] eq $b && $tp[1] eq $a)  {$geno_array[$kn]="h";} 
     	          if($tp[0] eq $a && $tp[1] eq $a)  {$geno_array[$kn]="a";}
     	          if($tp[0] eq $b && $tp[1] eq $b)  {$geno_array[$kn]="b";}
                  $kn++;
     	         }
     	#printf W1 "%s\t%s\t%s\t%s\t%s", $locus,$chro,$position,"0","2";
     	printf W2 "%s\t%s\t%s", $locus,$chro,$position;
     	foreach (@geno_array)
     	        {# printf W1 "\t%s",$_;
     	          printf W2 "\t%s",$_;
                }
        #print W1 "\n";              	  
     	print W2 "\n";
    }
    if( $_ =~ /\#CHROM/) {
		$flag=1;
	    my @tq=split(/\t/,$_);
		my @hd=@tq[9..434];
		printf W2 "%s\t%s\t%s", "locus","chro","position";
		foreach(@hd){printf W2 "\t%s",$_;}
		printf W2 "\n";
	}
}
  
close(R);
#close(W1);
close(W2);  
