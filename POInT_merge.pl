#!/usr/bin/perl

use strict;
use Bio::Tools::GFF;


my($i, $mixed, $used_genes, $infile, $outfile, $cnt_none, $line, @data, $genome, @genome_names, %ances_pillars, %merged_pillars, $num_genomes, @valid_pillars, $num_valid, $ances, $dupl1, $dupl2, @inferfiles, $file_of_lists, $ances_gff, $ances_gene, $cnt_present, %gene_chromos, $feature, $tag, $stream, @gene_starts, $gene_key, $ID, $full_name,
  $out_ances,  $nonsynt_cut, @nonsyns, %chromos, $chromo_counts, $chrom, $name, @orderfiles, $cnt_ns, @gene_chromos, @genes_in_order, @present_genes, $cnt, $new_order_file, @num_genes, $order_file, $dupl_level, $i, $j, @dupls, $dupls_used, @syn, $use, $id_file, @init_order, %id_gene_names, %fail_ids, @final_order, $valid_pillars, $size);

if (@ARGV < 3) {
        print "Usage: POInT_merge.pl <File of optimal homeologous genes> <GFF file for the non-polyploid_genome> <Outputfile>  (-s:#) (-d:#)\n"
        exit;
}

$file_of_lists=$ARGV[0];
$ances_gff=$ARGV[1];
$outfile=$ARGV[2];


$dupl_level=2;
$order_file="";
$id_file="";
$nonsynt_cut=0;

if (@ARGV > 3) {
    for($i=3; $i<@ARGV; $i++) {
        if ($ARGV[$i] =~ /^-s/i) {
            $ARGV[$i] =~ /^-s:(\d{1,3})/;
            $nonsynt_cut=$1;
        }
        if ($ARGV[$i] =~ /^-d/i) {
            $ARGV[$i] =~ /^-d:(\d{1,3})/;
            $dupl_level=$1;
        }
        if ($ARGV[$i] =~ /^-o/i) {
            $order_file=$ARGV[$i];
            $order_file =~ s/^-o://;
        }
        if ($ARGV[$i] =~ /^-i/i) {
            $id_file=$ARGV[$i];
            $id_file =~ s/^-i://;
        }
    }
}

if (($order_file ne "")  && ($id_file ne "")) {
    open(READDATA, "<" . $order_file) or die;
    $i=0;
    while($line =<READDATA>) {
        $line =~ s/\n//;
        @data=split(/\s/, $line);
        for($j=0; $j<@data; $j++) {
            $data[$j] =~ s/^-//;
            $init_order[$i+$j]=$data[$j];
        }
        $i= $i + @data;
    }
    close(READDATA);
    
    open(READDATA, "<" . $id_file) or die;
    while($line =<READDATA>) {
        $line =~ s/\n//;
        @data=split(/\t/, $line);
        if (!(exists $id_gene_names{$data[0]})) {
            $id_gene_names{$data[0]}=$data[1];
        }
        else {
            $fail_ids{$data[0]}=1;
        }
    }
    close(READDATA);
    $j=0;
    
    open(WRITEFAIL, ">Ancesorder_Fails" . $file_of_lists . ".txt") or die;
    for($i=0; $i<@init_order; $i++) {
        if (!(exists $id_gene_names{$init_order[$i]})) {
            print WRITEFAIL $init_order[$i], "\tNEVERFOUND\n";
        }
        else {
            if (!(exists $fail_ids{$init_order[$i]})) {
                $final_order[$j]=$id_gene_names{$init_order[$i]};
                $j++;
            }
            else {
                print WRITEFAIL $init_order[$i], "\tMULTIMAP\n";
            }
        }
    }
    close(WRITEFAIL);
    $valid_pillars=$j;
    print "Using provided ancestral order with $valid_pillars unique pillars\n;"
    
    
}

open(READDATA, "<" . $file_of_lists) or die;
$num_genomes=0;
while($line=<READDATA>) {
    $line=~ s/\n//;
    @data=split(/\t/, $line);
    $genome_names[$num_genomes]=$data[0];
    $inferfiles[$num_genomes]=$data[1];
    $orderfiles[$num_genomes]=$data[2];
    $num_genomes++;
}
close(READDATA);


for($genome=0; $genome<$num_genomes; $genome++) {
    $num_valid=0;
    $mixed=0;
    $used_genes=0;
    
    open(READDATA, "<" . $orderfiles[$genome]) or die;
    
    $line=<READDATA>;
    $cnt=0;
    while($line =<READDATA>) {
        $line =~ s/\n//;
        @data=split(/\t/, $line);
        $gene_chromos[$genome]{$data[1]}=$data[0];
        #print "Assigned chromosome ", $data[0], " to gene ", $data[1], " in the $genome th genome\n";
        $genes_in_order[$genome][$cnt]=$data[1];
        $cnt++;
    }
    $num_genes[$genome]=$cnt;
    close(READDATA);

    open(READDATA, "<" . $inferfiles[$genome]) or die;

    $line=<READDATA>;
    $line=<READDATA>;

    $cnt_none=0;

    while($line =<READDATA>) {
        $line =~ s/\n//;
        @data=split(/\t/, $line);
        $ances=$data[0];
        for($i=0; $i<$dupl_level; $i++) {
            $dupls[$i]=$data[2+(3*$i)];
            $syn[$i] =$data[3+(3*$i)];
        }
        
        $use=1;
        $cnt=0;
        for($i=0; $i<$dupl_level; $i++) {
            if ($syn[$i] eq "N") {$use=0;}
            if ($dupls[$i] ne "NONE") {$cnt++;}
        }

        if (($use == 1) && ($cnt>0)) {
            $valid_pillars[$num_valid] = $ances;
            for($i=0; $i<$dupl_level; $i++) {
                $valid_pillars[$num_valid] = $valid_pillars[$num_valid] . "\t" . $dupls[$i];
                $ances_pillars{$data[0]}{$genome_names[$genome]}[$i] = $dupls[$i];
                
                if ($dupls[$i] ne "NONE") {
                    $present_genes[$genome]{$dupls[$i]}=1;
                    $used_genes++;
                }
            }
            $num_valid++;
        }
        
        
    }
    close(READDATA);
    print $genome_names[$genome], ": Found $num_valid supported pillars\n";
    print genome_names[$genome], ": The valid pillars contain $used_genes total genes\n";

}


$stream =Bio::Tools::GFF->new(-file => $ances_gff, -gff_version => 3);
$chromo_counts=0;
while($feature = $stream->next_feature()) {
    # print $feature->primary_tag(), "\n"; # do something with feature
    if ($feature->primary_tag() eq "gene") {
        #print "Source tag: ", $feature->source_tag(), " at location : ", $feature->start(), " to ", $feature->end(), "\n";
        $name="";
        
        foreach $tag ($feature->get_all_tags()) {
            if ($tag eq "ID") {
                @data = $feature->get_tag_values($tag);
                $ID=$data[0];
            }
            if ($tag eq "Name") {
                @data = $feature->get_tag_values($tag);
                $name=$data[0];
            }
        }
 
	$full_name=$name;   
#            $full_name = $name . "_" . $feature->seq_id();
     
        if (!(exists $chromos{$feature->seq_id()})) {
            $chromos{$feature->seq_id()}=$chromo_counts;
            $chromo_counts++;
        }
        
        #print "Name: ", $name, " source tag: ", $feature->source_tag(), " at location : ", $feature->start(), " to ", $feature->end(), " ID: ", $feature->seq_id(), "\n";
        
        $gene_chromos{$full_name}=$chromos{$feature->seq_id()};
        
        if ($feature->start() < $feature->end()) {
            $gene_starts[$chromos{$feature->seq_id()}]{$full_name}=$feature->start();
        }
        else {
            $gene_starts[$chromos{$feature->seq_id()}]{$full_name}=$feature->end();
        }
        
    }
}

@present_genes=();

open(WRITEDATA, ">" . $outfile) or die;
$num_valid=0;

$out_ances=$outfile;
$out_ances =~ s/\.txt//;
$out_ances = $out_ances . "_ances_genes.txt";
open(WRITEANCES, ">" . $out_ances) or die;

print WRITEDATA  $genome_names[0];
for($genome=1; $genome<$num_genomes; $genome++) { print WRITEDATA "\t", $genome_names[$genome];}
print WRITEDATA "\n";

if (($order_file ne "")  && ($id_file ne "")) {
    for($i=0; $i<$valid_pillars; $i++) {
        $ances_gene=$final_order[$i];
        
        if (exists $ances_pillars{$ances_gene}) {
            $cnt_present=0;
            
            for($genome=0; $genome<$num_genomes; $genome++) {
                if (exists $ances_pillars{$ances_gene}{$genome_names[$genome]}) {$cnt_present++;}
            }
            print "Ancestral gene $ances_gene ", $gene_starts[$chrom]{$ances_gene}, " found in $cnt_present (syn) of the scaffolded genomes\n";
            
            if ($cnt_present == $num_genomes){
                $num_valid++;
                $line ="";
                for($j=0; $j<$dupl_level; $j++) {
                    for($genome=0; $genome<$num_genomes; $genome++) {
                        $line = $line . "\t" . $ances_pillars{$ances_gene}{$genome_names[$genome]}[$j];
                        $present_genes[$genome]{$ances_pillars{$ances_gene}{$genome_names[$genome]}[$j]}=1;
                    }
                }
                $line =~ s/^\t//;
                $line =~ s/\t$//;
                
		print WRITEANCES $num_valid-1, "\t", $ances_gene, "\t", $line,  "\n";
                print WRITEDATA $line, "\n";
                
            }
        }
    }
}
else {
    for ($chrom=0; $chrom<$chromo_counts; $chrom++) {
        foreach $ances_gene (sort {$gene_starts[$chrom]{$a} <=> $gene_starts[$chrom]{$b}} keys %{$gene_starts[$chrom]}) {
            if (exists $ances_pillars{$ances_gene}) {
                $cnt_present=0;
        
                for($genome=0; $genome<$num_genomes; $genome++) {
                    if (exists $ances_pillars{$ances_gene}{$genome_names[$genome]}) {$cnt_present++;}
                }
                 print "Ancestral gene $ances_gene ", $gene_starts[$chrom]{$ances_gene}, " found in $cnt_present (syn) of the scaffolded genomes\n";
                
                if ($cnt_present == $num_genomes){
                    $num_valid++;
                    $line ="";
                    for($i=0; $i<$dupl_level; $i++) {
                        for($genome=0; $genome<$num_genomes; $genome++) {
                            $line = $line . "\t" . $ances_pillars{$ances_gene}{$genome_names[$genome]}[$i];
                            $present_genes[$genome]{$ances_pillars{$ances_gene}{$genome_names[$genome]}[$i]}=1;
                        }
                    }
                    $line =~ s/^\t//;
                    $line =~ s/\t$//;
        	    print WRITEANCES $num_valid-1, "\t", $ances_gene, "\t", $line, "\n";
                    print WRITEDATA $line, "\n";
                
                }
            }
        }
    }
}
close(WRITEDATA);
close(WRITEANCES);
print "Wrote $num_valid total pillars to output file\n";

for($genome=0; $genome<$num_genomes; $genome++) {
    $new_order_file = $genome_names[$genome] . "_POInT_geneorders.txt";
    open(WRITEDATA, ">" . $new_order_file) or die;
    print WRITEDATA $genome_names[$genome], "\n";
    for($i=0; $i<$num_genes[$genome]; $i++) {
        if (exists $present_genes[$genome]{$genes_in_order[$genome][$i]}) {
            print WRITEDATA $gene_chromos[$genome]{$genes_in_order[$genome][$i]}, "\t", $genes_in_order[$genome][$i], "\n";
        }
    }
    close(WRITEDATA);
}

