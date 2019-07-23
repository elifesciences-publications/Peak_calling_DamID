#!usr/bin/perl
#peak_finding_with_FDR
use warnings;

#print "\nFiles should be either tab delimited .gff files, or .bedgraph files, with values for each GATC fragment\n\n";
print "\nFiles should be tab delimited with values for each GATC fragment\n\n";

@infiles = glob("*.bedgraph");
@infiles = glob("*.gff") unless scalar(@infiles);

#$file_num = @infiles; if($file_num ==0){die "\nNo files detected! - make sure they are .gff or .bedgraph files\n"};
$file_num = @infiles; if($file_num ==0){die "\nNo files detected! - make sure they are .gff files\n"};
	
use Cwd;
my $dir = getcwd;
$dir = substr($dir, 2);    # gives current directory for calling up other perl programs

$rep_num = 0;
    
print "Files for processing are:\n\n";

foreach $a (@infiles){
	if($a =~ m/rep1/){$path1 = $a; chomp $path1; print "1st replicate is $a\n"; $rep_num = $rep_num + 1;}
	if($a =~ m/rep2/){$path2 = $a; chomp $path2; print "2nd replicate is $a\n"; $rep_num = $rep_num + 1;}
	if($a =~ m/rep3/){$path3 = $a; chomp $path3; print "3rd replicate is $a\n"; $rep_num = $rep_num + 1;}
	}

if($file_num != $rep_num){die "\n\nNumber of files does not match number of replicate data files\n- make sure all the replicate data files contain \"rep1\" or \"rep2\" etc\n";}


print "\nNumber of replicates is $rep_num\n";
if($rep_num == 1){$paths = $path1;}
if($rep_num == 2){$paths = "$path1\t$path2";}
if($rep_num == 3){$paths = "$path1\t$path2\t$path3";}

print "\nEnter name for analysis\n\n";
$exp_name = <STDIN>;
chomp $exp_name;

print "\nEnter FDR threshold \(usually 0.01 \(1 percent\)\)\n\n";
$FDR_thres = <STDIN>;
chomp $FDR_thres;


#print "\nEnter ratio threshold \(usually 1\)\n\n";
$data_thres = 0.2;  # threshold of all ratio values that make up a peak
#chomp $data_thres;


mkdir 'FDR_analysis_for_'."$exp_name", 0755 or die "\nCan't make analysis directory!\n";
mkdir 'FDR_analysis_for_'."$exp_name".'\logfiles', 0755 or die "\nCan't make analysis logfiles directory!\n";
mkdir 'FDR_analysis_for_'."$exp_name".'\logfiles\shuffled_data', 0755 or die "\nCan't make shuffled data analysis logfiles directory!\n";
mkdir 'FDR_analysis_for_'."$exp_name".'\logfiles\freq_data', 0755 or die "\nCan't make freq data analysis logfiles directory!\n";
mkdir 'FDR_analysis_for_'."$exp_name".'\FDR_results', 0755 or die "\nCan't make FDR results directory!\n";



#####################  make average and moving average files #######################

system("$dir/sub_programs/split_files_and_average.pl", "$exp_name", "$paths", "$rep_num");  #starts perl prog to make RPM normalised files and average file


##################### Now to start FDR analysis##################################################################

$replicate = 1;

while($replicate < ($rep_num + 1)){

chrom_analysis("2L");
chrom_analysis("2R");
chrom_analysis("3L");
chrom_analysis("3R");
chrom_analysis("4");
chrom_analysis("X");

sub chrom_analysis{  ####don't forget to close bracket!
	
	my $chrom = shift;
	my $chrom2 = "chr"."$chrom"; print "\n$chrom2\n";



################ now to randomise data set and calculate peak frequency ###############

system("$dir/sub_programs/random_permutation_analysis.pl", "$exp_name", "$chrom", "$replicate", "$data_thres", "$dir"); #starts perl prog to randomise data set and look for peaks


########################now to look at real data###################################################

system("$dir/sub_programs/real_data_analysis.pl", "$exp_name", "$chrom", "$replicate", "$data_thres", "$dir");  #starts perl prog to examine real data set and look for peaks


##############################now to merge FDR values and data########################################

system("$dir/sub_programs/merge_FDR_values_and_data.pl", "$exp_name", "$chrom", "$replicate", "$dir");  #starts perl prog to merge FDR values and data


##################### find all regions with FDR less than 0.0001 ####################################

system("$dir/sub_programs/assign_FDR_less_than_0.0001.pl", "$exp_name", "$chrom", "$replicate", "$data_thres", "$dir");  #Find regions with FDR < 0.0001

system("$dir/sub_programs/assign_each_GATC_lowest_fdr.pl", "$exp_name", "$chrom", "$replicate", "$dir");

system("$dir/sub_programs/make_peakfile.pl", "$exp_name", "$chrom", "$replicate", "$FDR_thres", "$dir");

}
$replicate = $replicate + 1;
}


system("$dir/sub_programs/merge_peak_files.pl", "$exp_name", "$rep_num", "$dir");

system("$dir/sub_programs/assign_final_peaks_FDR_and_height.pl", "$exp_name", "$rep_num", "$FDR_thres", "$dir");

system("$dir/sub_programs/compile_final_peak_file.pl", "$exp_name", "$FDR_thres", "$dir");

system("$dir/sub_programs/assign_peak_to_gene.pl", "$exp_name", "$FDR_thres", "$dir");



################################  make unique gene list #################################

system("$dir/sub_programs/make_unique_gene_list.pl", "$exp_name", "$dir");



