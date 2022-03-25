#!/usr/bin/perl

use strict;

my($infile, $outfile, %pillar_data, $line, @data, $len, @anc_locs, @sort_anc_locs, $my_gene, $my_chrom, $my_loc, $subloc, $target_genome, $num_genomes, $loc, $parse, $val, @bits, %ortho_probs, $prob, $gene_set, $max_loc, @groups, @order,
$first_spp, $singles_only, $anc_chroms, $j, $i, @names, %taxa_present, %genes_present, $my_anc, $use, $k, $filename, $outgenome, $tag, @omit, $num_poss_track, @first_tracks, $track_num, $genome, $the_track, @genes, @best_tracks, $max, @probs, $site, $do_genome, $genome_steps, $use_degen, $stop, $other_track, @flip_tracks, $mytrack, @data2, $fliptrack, %track_ids, @flip_track_nums, $new_prob, $ploidy, %genome_names, $name_key, @powN, $combos, $loc);

if (@ARGV < 2) {
	print "Usage: POInT_extract_orthos.pl <PostProbsFile> <NewOrthologFile>\n";
	exit();
}
$infile=$ARGV[0];
$outfile=$ARGV[1];


$use_degen=0;



$num_genomes=0;


open(READDATA, "<" . $infile) or die;

$line=<READDATA>;
$line=~s/\n//;

@data=split(/\t/, $line);

$i=0;
while(!($data[$i] =~ /^#/)) {
    $genome_names{$data[$i]}=1;
    $i++;
}

$num_genomes=0;

foreach $name_key (keys %genome_names) {
    $num_genomes++;
    #print "Genome $name_key is $num_genomes\n";
}

$ploidy=$i/$num_genomes;
if ($ploidy == 2) {$combos=2;}
if ($ploidy == 3) {$combos=6;}
if ($ploidy == 4) {$combos=24;}


$powN[0]=1;
for ($i=1; $i<=$num_genomes; $i++) {
    $powN[$i]=$powN[$i-1]*$combos;
}

print "Ploidy of this event is $ploidy\n";
#exit();

#while(!($data[$num_genomes] =~ /^#/)) {$num_genomes++;}
$num_genomes=0;

#$num_genomes=$num_genomes/2;
$i=0;
while(!($data[$i] =~ /^#/)) {
    $names[$num_genomes]=$data[$i];
    print "Genome ",$num_genomes,  " is ", $names[$num_genomes], "\n";
    $num_genomes++;
    $i+=$ploidy;
}

$genome_steps=$ploidy*$num_genomes;

$loc=0;
for($j=$genome_steps; $j<$powN[$num_genomes]+$genome_steps; $j++) {
    $parse=$data[$j];
    $parse =~ s/^#//;
    #print "For ", $j-$genome_steps, " $parse ";
    #print "\n";
    for($i=0; $i<$num_genomes; $i++) {
        $parse =~ s/^_//;
        $parse =~ s/$names[$i]//;
        $parse =~ /_?(\d{2,3})/;
        #$data[$j] =~ /^#(\w{12})/;
        $val=$1;
        $parse =~ s/^$val//;
        @bits=split(//, $val);
        #print "| P", $parse, ": ", $i, " is ",$val;
        #$order[$j-10][0]=0;
        for($k=0; $k<$powN[1]; $k++) {
            $order[$j-$genome_steps][$i][$k]=$bits[$k];
            #print "NG $i, Pos $k= ", $order[$j-$genome_steps][$i][$k], "; ";
            
        }
        
    }
    #print "\n";
}



$site=0;

while ($line = <READDATA>) {
	$line =~ s/\n//;
    @data=split(/\t/, $line);
    
    %ortho_probs=();
    
    for($j=0; $j<$num_genomes; $j++) {
        for($i=0; $i<$ploidy; $i++) {
            $groups[$j][$i]=$data[$j*$ploidy+$i];
            #print "Taxa $j loc $i is ", $groups[$j][$i], "\n";
            #for($k=0; $k<3; $k++) {$subgenome_probs[$j][$i][$k]=0.0;}
        }
    }
    
    
    for($j=$genome_steps; $j<($powN[$num_genomes]+$genome_steps); $j++) {
        $prob = 1.0*$data[$j];
        
        $gene_set="";
        
        for($i=0; $i<$ploidy; $i++) {
            for($k=0; $k<$num_genomes; $k++) {
                $gene_set = $gene_set ."\t" . $groups[$k][$order[$j-$genome_steps][$k][$i]];
                
                #$subgenome_probs[$i][$order[$j-$offset][$i][$k]][$k]+=$prob;
            }
        }
        $gene_set =~ s/^\t//;
        
        
        if(exists $ortho_probs{$gene_set}) {
            $ortho_probs{$gene_set} = $ortho_probs{$gene_set} + $prob;
        }
        else {
            $ortho_probs{$gene_set} = $prob;
        }
    }
    $max=-1;
    $max_loc=-1;
    
    foreach $gene_set (keys %ortho_probs) {
        if ($max < $ortho_probs{$gene_set}) {
            $max = $ortho_probs{$gene_set};
            $max_loc=$gene_set;
        }
    }
    
    print $site, "\t", $max, "\t", $max_loc, "\n";
    
    @data2=split(/\t/, $max_loc);
    
    $j=0;
    for($i=0; $i<$ploidy; $i++) {
        for($k=0; $k<$num_genomes; $k++) {
            $best_tracks[$site][$k][$i]=$data2[$j];
            #print "$site: Genome: $k, pos: $i = ", $best_tracks[$site][$k][$i], "\n";
            $j++;
        }
    }
    $probs[$site]=$max;
    $site++;
}
close(READDATA);

open(WRITEDATA, ">" . $outfile) or die;
print WRITEDATA "SITE\tBESTPROB";
for ($genome=0; $genome<$num_genomes; $genome++) {
    for($j=0; $j<$ploidy; $j++) {
        print WRITEDATA  "\t", $names[$genome], $j+1;
    }
   
}
print WRITEDATA "\n";

for($i=0; $i<$site; $i++) {
    print WRITEDATA $i, "\t", $probs[$i];
    for ($genome=0; $genome<$num_genomes; $genome++) {
        for($j=0; $j<$ploidy; $j++) {
            print WRITEDATA  "\t", $best_tracks[$i][$genome][$j];
        }
    }
    print WRITEDATA "\n";
}
