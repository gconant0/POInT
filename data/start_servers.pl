#!/usr/bin/perl

use strict;
use Cwd;


my($infile, @events, $i, $line, @data, $command, $out, $arg, $restart, @pids, $j, $work_dir, $level, $model, $stop, $other_order, @data2, $root_model);


$infile=$ARGV[0];

$restart=0;


open(READDATA, "<" . $infile) or die;

if (@ARGV >1) {$restart=1;}

$line=<READDATA>;

$work_dir=Cwd::abs_path;

while($line=<READDATA>) {
    $line =~ s/\n//;
    @data=split(/\t/, $line);
    
    $model=$data[3];
    $root_model="";
    if ($model =~ /\+/) {
          @data2=split(/\+/, $model);
          $model=$data2[0];
          $root_model=$data2[1];
    }
 
    if ($restart==1) {
        $arg = "\"" . $data[1] . "\"";
        @pids=`ps aux |grep $arg`;
        
        for($j=0; $j<@pids; $j++) {
            $pids[$j] =~ s/\n//;
            if ($pids[$j] =~ /POInT_daemon/) {
                @data2=split(/\s{1,5}/, $pids[$j]);
                $arg=$data2[1];
                `kill -9 $arg`;
            }
        }
    }
    
     $arg = "\"" . $data[1] . "\"";
    $out=`ps aux |grep $arg`;
    
    
    print "Looking @ ", $data[1], ": $out\n";
    
    if(!($out =~ /POInT_daemon/)) {
        print "Event ", $data[0], " is not running\n";
        
        if (-e $data[1]) {`rm $data[1]`;}
        
        
        $arg=$data[0] ."/" . $model;
        $out = `cat $arg |grep Hierarchical`;
        $out =~ /Hierarchical\s{1,4}(\d)/;
        
        $level = $1;
        
        $command = "POInT_daemon -d:" . $level;
        
        $i=5;
        
        if ($data[5] =~ /ancesorder/) {
            $other_order="";

            for ($i=5; $i<(@data); $i++) {
		if (!($data[$i] =~ /\.png$/)) {
                	if ($data[$i] =~ /^-[kK]:/) {
                    		$other_order = $other_order . " " . $data[$i];
               	 	}
                	else {
                    	if ($data[$i] =~ /ancesorder/) {
                        	print "Adding $i ", $data[$i], " to $command\n";
                        	$command = $command . " -g:" . $data[$i];
                    	}
                	    else {
                        	$other_order = $other_order . " -a:" . $data[$i];
                    		}
                	}
		}
                
            }
            $command = $command . " -o:" . $data[4] . " -t:" . $data[2];
	   if ($root_model eq "") {$command = $command . " -m:" . $model . " -noopt -w:" . $data[1] . " " . $other_order;}
	   else {$command = $command . " -m:" . $model . " -r:" . $root_model . " -noopt -w:" . $data[1] . " " . $other_order;}
            
        }
        else {
            while (($data[$i] =~ /POInT_geneorder/) && ($i<@data)) {$i++;}
            $stop=$i;
            
            print "Found 5 to $stop files for point orders\n";
            
            for ($i=5; $i<$stop; $i++) {
                print "Adding $i ", $data[$i], " to $command\n";
                $command = $command . " -g:" . $data[$i];
                
            }
            $command = $command . " -o:" . $data[4] . " -t:" . $data[2];
	    if ($root_model eq "") {$command = $command . " -m:" . $model . " -noopt -w:" . $data[1];}
	    else { $command = $command . " -m:" . $model . " -r:" . $root_model . " -noopt -w:" . $data[1];}
            
            $other_order="";
            
            if ($stop <@data) {
                for($i=$stop; $i<@data; $i++) {
			if (!($data[$i] =~ /\.png$/)) {
                    		if ($data[$i] =~ /^-[kK]:/) {
                        		$other_order = $other_order . " " . $data[$i];
                    		}
                    		else 	{
                        		$other_order = $other_order . " -a:" . $data[$i];
                    		}
                    		print "Adding ", $data[$i], " to $other_order\n";
			}
                }
            }
            
            $command = $command . " " . $other_order;
        }
        my $pid = fork();
        if ($pid == 0)
        {
            # This is the child process.
            # exec() the external program.
            chdir $data[0];
            print "Running $command \n";
            exec($command)
            or die "could not exec $command: $!";
            
            
            
            #`$command &`;
            #$arg= $data[1];
            #$out=`chmod o+w $arg`;
            #my $mode = 0644;   chmod $mode, $arg;
            #print "Set permissions on $arg resulted in $out\n";
            
            chdir $work_dir;
        }
        elsif (!defined($pid))
        {
            die "could not fork";
        }
        
        
    }
    else {
        chdir $data[0];
        $arg= $data[1];
        my $mode = 0777;   chmod $mode, $arg;
        print "Set permissions on $arg resulted in $out\n";
        chdir $work_dir;
    }
    
    
    
}
