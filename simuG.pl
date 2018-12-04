#!/usr/bin/perl
use warnings FATAL => 'all';
use strict;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);
use List::Util qw(sum min max shuffle);
# use Data::Dumper;

##############################################################
#  script: simuG.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  version: 1.0.0
#  last edited: 2018.12.06
#  description: simuG.pl can simulate genome sequences with designated or random genomic variants of full spectrum (e.g. SNP, INDEL, CNV, inversions, and translocations).
##############################################################

################
# input parameters
################

# General options
my $help;
my $man;
my $version;
# Options for the reference genome
my $refseq;
my $excluded_chr_list;
# Options for SNP variants simulation
my $snp_vcf;
my $snp_count;
my $snp_model;
my $titv_ratio = 0.5;
# Options for INDEL variants simulation
my $indel_vcf;
my $indel_count;
my $indel_model;
my $ins_del_ratio = 1.0;
my $indel_size_powerlaw_alpha = 2.0;
my $indel_size_powerlaw_constant = 0.5;
# Options for copy-number variants simulation
my $cnv_vcf;
my $cnv_count;
my $cnv_gain_loss_ratio = 1.0;
my $cnv_max_copy_number = 10;
my $cnv_min_size = 100;
my $cnv_max_size = 10000;
# Options for inversions simulation
my $inversion_vcf;
my $inversion_count;
my $inversion_min_size = 1000;
my $inversion_max_size = 1000000;
my $inversion_breakpoint_gff;
# Options for translocation simulation
my $translocation_vcf;
my $translocation_count;
my $translocation_breakpoint_gff;
# Option for defining centromere for CNV/inversion/translocation simulation
my $centromere_gff;
# Options for setting random seed and output file prefix
my $seed;
my $prefix = "output_prefix";

GetOptions(
'help|h|?' => \$help,
'man|m' => \$man,
'version|ver' => \$version,
'refseq|r:s' => \$refseq,
'excluded_chr_list:s' => \$excluded_chr_list,
'snp_vcf:s' => \$snp_vcf,
'snp_count:i' => \$snp_count,
'snp_model:s' => \$snp_model,
'titv_ratio:f' => \$titv_ratio,
'indel_vcf:s' => \$indel_vcf,
'indel_count:i' => \$indel_count,
'indel_model:s' => \$indel_model,
'ins_del_ratio:f' => \$ins_del_ratio,
'indel_size_powerlaw_alpha:f' => \$indel_size_powerlaw_alpha,
'indel_size_powerlaw_constant:f' => \$indel_size_powerlaw_constant,
'cnv_vcf:s' => \$cnv_vcf,
'cnv_count:i' => \$cnv_count,
'cnv_gain_loss_ratio:f' => \$cnv_gain_loss_ratio,
'cnv_max_copy_number:i' => \$cnv_max_copy_number,
'cnv_min_size:i' => \$cnv_min_size,
'cnv_max_size:i' => \$cnv_max_size,
'inversion_vcf:s' => \$inversion_vcf,
'inversion_count:i' => \$inversion_count,
'inversion_min_size:i' => \$inversion_min_size,
'inversion_max_size:i' => \$inversion_max_size,
'inversion_breakpoint_gff:s' => \$inversion_breakpoint_gff,
'translocation_vcf:s' => \$translocation_vcf,
'translocation_count:i' => \$translocation_count,
'translocation_breakpoint_gff:s' => \$translocation_breakpoint_gff,
'centromere_gff:s' => \$centromere_gff,
'seed|s:i' => \$seed,
'prefix|p:s' => \$prefix,
);


## Option check and print out Usage if needed.
# If usage (or version number) was explicitly requested, print usage (or version number).

if (defined $help) {
  pod2usage(-verbose => 1);
}

if (defined $man) {
  pod2usage(-verbose => 2);
}

if (defined $version) {
  pod2usage(-verbose => 99, -sections => qw(NAME|DESCRIPTION|VERSION));
}

if (not defined $refseq) {
  pod2usage(-message => "Mandatory argument '-refseq' is missing! Exit!",
  -exitval => 1,
  -verbose => 1);
} elsif (not -e $refseq) {
  print "!!! Error! The defined input file $refseq does not exist!\n";
  print "!!! Exit!\n";
  die;
}

## Main program

print "\n";
my $local_time = localtime();
print "[$local_time]\n";
print "Starting simuG ..\n\n";
$local_time = localtime();
print "[$local_time]\n";
print "Check specified options ..\n";
# check specified options
if ((defined $snp_vcf) or (defined $snp_count) or (defined $indel_vcf) or (defined $indel_count)) {
  print "Running simuG for SNP/INDEL simulation >>\n";
  print "Ignore all options for CNV/inversion/translocation simulation.\n";
  undef $centromere_gff;
  undef $cnv_vcf;
  undef $cnv_count;
  undef $cnv_gain_loss_ratio;
  undef $cnv_max_copy_number;
  undef $cnv_min_size;
  undef $cnv_max_size;
  undef $inversion_vcf;
  undef $inversion_count;
  undef $inversion_min_size;
  undef $inversion_max_size;
  undef $inversion_breakpoint_gff;
  undef $translocation_vcf;
  undef $translocation_count;
  undef $translocation_breakpoint_gff;
} elsif ((defined $cnv_vcf) or (defined $cnv_count)) {
  print "Running simuG for CNV simulation >>\n\n";
  print "Ignore all options for inversion/translocation simulation.\n";
  undef $inversion_vcf;
  undef $inversion_count;
  undef $inversion_min_size;
  undef $inversion_max_size;
  undef $inversion_breakpoint_gff;
  undef $translocation_vcf;
  undef $translocation_count;
  undef $translocation_breakpoint_gff;
} elsif ((defined $inversion_vcf) or (defined $inversion_count)) {
  print "Running simuG for inversion simulation >>\n\n";
  print "Ignore all options for translocation simulation.\n";
  undef $translocation_vcf;
  undef $translocation_count;
  undef $translocation_breakpoint_gff;
} elsif ((defined $translocation_vcf) or (defined $translocation_count)) {
  print "Running simuG for translocation simulation >>\n\n";
} else {
  print "\n";
  print "!!! There seems no task for simuG to run !!!\n";
  print "!!! 1) If you want to run simuG for SNP/INDEL simulation, at least one of the following options need to be specified:\n";
  print "!!!    -snp_vcf, -snp_count, -indel_vcf, -indel_count\n";
  print "!!! 2) If you want to run simuG for CNV simulation, at least one of the following options need to be specified:\n";
  print "!!!    -cnv_vcf, -cnv_count\n";
  print "!!! 3) If you want to run simuG for inversion simulation, at least one of the following options need to be specified:\n";
  print "!!!    -inversion_vcf, -inversion_count\n";
  print "!!! 4) If you want to run simuG for translocation simulation, at least one of the following options need to be specified:\n";
  print "!!!    -translocation_vcf, translocation_count\n";
  print "!!! 5) If you want to check all available options for simuG, type: perl simuG.pl -h\n";
  print "!!! Exit !!!\n\n";
  $local_time = localtime();
  print "[$local_time]\n";
  exit;
}


# define excluded chromosome(s) if any
my %excluded_chr_list = ();
if (defined $excluded_chr_list) {
  my $excluded_chr_list_fh = read_file($excluded_chr_list);
  %excluded_chr_list = parse_list_file($excluded_chr_list_fh);
  $local_time = localtime();
  print "\n[$local_time]\n";
  print "Check for excluded chromosome(s) defined in $excluded_chr_list ..\n";

}

# set up reference genome as the template
my @refseq = ();
my %refseq = ();
my $refseq_fh = read_file($refseq);
parse_fasta_file($refseq_fh, \%refseq, \@refseq);
close $refseq_fh;

# remove excluded chromosomes if specified
foreach my $chr (@refseq) {
  if (exists $excluded_chr_list{$chr}) {
    print "Exclude chromosome: $chr\n";
    delete $refseq{$chr};
  }
}

# initialize the simulated genome
my %simseq = %refseq;
my %ref2sim_map = ();

# initialize the seed for random number generator
if (not defined $seed) {
    $seed = int(rand(2**31));
}
print "\nThis simulation use the random seed: $seed\n\n";
srand($seed);

if (defined $snp_vcf) {
  print "The option snp_vcf has been specified: snp_vcf = $snp_vcf\n";
  print "Ignore incompatible option: snp_count\n";
  undef $snp_count;
  print "Ignore incompatible option: snp_model\n";
  undef $snp_model;
  print "Ignore incompatible option: titv_ratio\n";
  undef $titv_ratio;
} elsif (defined $snp_count) {
  print "The option snp_count has been specified: snp_count = $snp_count\n";
  if (defined $snp_model) {
    print "The option snp_model has been specified: snp_model = $snp_model\n";
    print "Ignore incompatible option: titv_ratio\n";
    undef $titv_ratio;
  } else {
    print "The option titv_ratio has been specified: titv_ratio = $titv_ratio\n";
  }
}

if (defined $indel_vcf) {
  print "The option indel_vcf has been specified: indel_vcf = $indel_vcf\n";
  print "Ignore incompatible option: indel_count\n";
  undef $indel_count;
  print "Ignore incompatible option: indel_model\n";
  undef $indel_model;
  print "Ignore incompatible option: ins_del_ratio\n";
  undef $ins_del_ratio;
  print "Ignore incompatible option: indel_size_powerlaw_alpha\n";
  undef $indel_size_powerlaw_alpha;
  print "Ignore incompatible option: indel_size_powerlaw_constant\n";
  undef $indel_size_powerlaw_constant;
} elsif (defined $indel_count) {
  print "The option indel_count has been specified: indel_count = $indel_count\n";
  if (defined $indel_model) {
    print "The option indel_model has been specified: indel_model = $indel_model\n";
    print "Ignore incompatible option: ins_del_ratio\n";
    undef $ins_del_ratio;
    print "Ignore incompatible option: indel_size_powerlaw_alpha\n";
    undef $indel_size_powerlaw_alpha;
    print "Ignore incompatible option: indel_size_powerlaw_constant\n";
    undef $indel_size_powerlaw_constant;
  } else {
    print "The option ins_del_ratio has been specified: ins_del_ratio = $ins_del_ratio\n";
    print "The option indel_size_powerlaw_alpha has been specified: indel_size_powerlaw_alpha = $indel_size_powerlaw_alpha\n";
    print "The option indel_size_powerlaw_constant has been specified: indel_size_powerlaw_constant = $indel_size_powerlaw_constant\n";
  }
}

if (defined $cnv_vcf) {
  print "The option cnv_vcf has been specified: cnv_vcf = $cnv_vcf\n";
  print "Ignore incompatible option: cnv_count\n";
  undef $cnv_count;
  print "Ignore incompatible option: cnv_gain_loss_ratio\n";
  undef $cnv_gain_loss_ratio;
  print "Ignore incompatible option: cnv_max_copy_number\n";
  undef $cnv_max_copy_number;
  print "Ignore incompatible option: cnv_min_size\n";
  undef $cnv_min_size;
  print "Ignore incompatible option: cnv_max_size\n";
  undef $cnv_max_size;
  print "Ignore incompatible option: centromere_gff\n";
  undef $centromere_gff;
} elsif (defined $cnv_count) {
  print "The option cnv_count has been specified: cnv_count = $cnv_count\n";
  print "The option cnv_min_size has been specified: cnv_min_size = $cnv_min_size\n";
  print "The option cnv_min_size has been specified: cnv_max_size = $cnv_max_size\n";
  if (defined $centromere_gff) {
    print "The option centromere_gff has been specified: centromere_gff = $centromere_gff\n";
  }
}

if (defined $inversion_vcf) {
  print "The option inversion_vcf has been specified: inversion_vcf = $inversion_vcf\n";
  print "Ignore incompatible option: inversion_count\n";
  undef $inversion_count;
  print "Ignore incompatible option: inversion_min_size\n";
  undef $inversion_min_size;
  print "Ignore incompatible option: inversion_max_size\n";
  undef $inversion_max_size;
  print "Ignore incompatible option: inversion_breakpoint_gff\n";
  undef $inversion_breakpoint_gff;
  print "Ignore incompatible option: centromere_gff\n";
  undef $centromere_gff;
} elsif (defined $inversion_count) {
  print "The option inversion_count has been specified: inversion_count = $inversion_count\n";
  if (defined $inversion_breakpoint_gff) {
    print "The option inversion_breakpoint_gff has been specified: inversion_breakpoint_gff = $inversion_breakpoint_gff\n";
    print "Ignore incompatible option: inversion_min_size\n";
    undef $inversion_min_size;
    print "Ignore incompatible option: inversion_max_size\n";
    undef $inversion_max_size;
  } else {
    print "The option inversion_min_size has been specified: inversion_min_size = $inversion_min_size\n";
    print "The option inversion_max_size has been specified: inversion_max_size = $inversion_max_size\n";
  }
  if (defined $centromere_gff) {
    print "The option centromere_gff has been specified: centromere_gff = $centromere_gff\n";
  }
}

if (defined $translocation_vcf) {
  print "The option translocation_vcf has been specified: translocation_vcf = $translocation_vcf\n";
  print "Ignore incompatible option: translocation_count\n";
  undef $translocation_count;
  print "Ignore incompatible option: translocation_breakpoint_gff\n";
  undef $translocation_breakpoint_gff;
  print "Ignore incompatible option: centromere_gff\n";
  undef $centromere_gff;
} elsif (defined $translocation_count) {
  print "The option translocation_count has been specified: translocation_count = $translocation_count\n";
  if (defined $translocation_breakpoint_gff) {
    print "The option translocation_breakpoint_gff has been specified: translocation_breakpoint_gff = $translocation_breakpoint_gff\n";
  }
  if (defined $centromere_gff) {
    print "The option centromere_gff has been specified: centromere_gff = $centromere_gff\n";
  }
}

print "\n";

my %snp_vcf = ();
if (defined $snp_vcf) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Parsing the input vcf file: $snp_vcf\n\n";
  my $snp_vcf_fh = read_file($snp_vcf);
  %snp_vcf = parse_simple_vcf_file($snp_vcf_fh, 0, 'SNP');
  close $snp_vcf_fh;
}

my %indel_vcf = ();
if (defined $indel_vcf) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Parsing the input vcf file: $indel_vcf\n\n";
  my $indel_vcf_fh = read_file($indel_vcf);
  %indel_vcf = parse_simple_vcf_file($indel_vcf_fh, 0, 'INDEL');
  close $indel_vcf_fh;
}

my %vcf = ();
if ((defined $snp_vcf) and (defined $indel_vcf)) {
  %vcf = merge_vcf(\%snp_vcf, \%indel_vcf);
} elsif (defined $snp_vcf) {
  %vcf = %snp_vcf;
} elsif (defined $indel_vcf) {
  %vcf = %indel_vcf;
}

# introduce SNP/INDEL variants based on user-provided vcf(s)
if (%vcf) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing defined SNP/INDELs based on the input vcf file(s):\n";
  if (defined $snp_vcf) {
    print "> snp_vcf = $snp_vcf\n";
  }
  if (defined $indel_vcf) {
    print "> indel_vcf = $indel_vcf\n";
  }
  introduce_defined_snp_indel(\%vcf, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce random SNP variants
if (defined $snp_count) {
  $local_time = localtime();
  print "[$local_time] Introducing random SNPs based on the following parameters:\n";
  print "> snp_count = $snp_count\n";
  if (defined $snp_model) {
    print "> snp_model = $snp_model\n";

  } else {
    print "> titv_ratio = $titv_ratio\n";
  }
  introduce_random_snp($snp_count, $snp_model, $titv_ratio, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce random INDEL variants
if (defined $indel_count) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing random INDELs based on the following parameters:\n";
  print "> indel_count = $indel_count\n";
  if (defined $indel_model) {
    print "> indel_model = $indel_model\n";
  } else {
    print "> ins_del_ratio = $ins_del_ratio\n";
    print "> indel_size_powerlaw_alpha = $indel_size_powerlaw_alpha\n";
    print "> indel_size_powerlaw_constant = $indel_size_powerlaw_constant\n";
  }
  introduce_random_indel($indel_count, $indel_model, $ins_del_ratio, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce CNVs based on user-provided CNV VCF file.
if (defined $cnv_vcf) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing defined CNVs based on the input vcf file:\n";
  print "> cnv_vcf = $cnv_vcf\n";
  my $cnv_vcf_fh = read_file($cnv_vcf);
  my %sv = parse_sv_vcf_file($cnv_vcf_fh);
  my %cnv = extract_cnv_from_sv(\%sv);
  introduce_defined_cnv(\%cnv, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce random CNVs
if (defined $cnv_count) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing random CNVs with the following parameters:\n";
  print "> cnv_count = $cnv_count\n";
  print "> cnv_gain_loss_ratio = $cnv_gain_loss_ratio\n";
  print "> cnv_max_copy_number = $cnv_max_copy_number\n";
  print "> cnv_min_size = $cnv_min_size\n";
  print "> cnv_max_size = $cnv_max_size\n";
  if (defined $centromere_gff) {
    print "> centromere_gff = $centromere_gff\n";
  }
  introduce_random_cnv($cnv_count, $cnv_gain_loss_ratio, $cnv_min_size, $cnv_max_size, $cnv_max_copy_number, $centromere_gff, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce inversions based on the input inversion vcf file
if (defined $inversion_vcf) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing defined Inversions based on the input vcf file:\n";
  print "> inversion_vcf = $inversion_vcf\n";
  my $inversion_vcf_fh = read_file($inversion_vcf);
  my %sv = parse_sv_vcf_file($inversion_vcf_fh);
  my %inversion = extract_inversion_from_sv(\%sv);
  introduce_defined_inversion(\%inversion, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce random inversions
if (defined $inversion_count) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing random Inversions based on the following parameters:\n";
  print "> inversion_count = $inversion_count\n";
  if (defined $centromere_gff) {
    print "> centromere_gff = $centromere_gff\n";
  }
  if (defined $inversion_breakpoint_gff) {
    print "> inversion_breakpoint_gff = $inversion_breakpoint_gff\n";
  } else {
    print "> inversion_min_size = $inversion_min_size\n";
    print "> inversion_max_size = $inversion_max_size\n";
  }
  introduce_random_inversion($inversion_count, $inversion_min_size, $inversion_max_size, $centromere_gff, $inversion_breakpoint_gff, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce translocations based on the input translocation vcf file
if (defined $translocation_vcf) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing defined Translocations based on the input vcf file:\n";
  print "> translocation_vcf = $translocation_vcf\n";
  my $translocation_vcf_fh = read_file($translocation_vcf);
  my %sv = parse_sv_vcf_file($translocation_vcf_fh);
  my %translocation = extract_translocation_from_sv(\%sv);
  introduce_defined_translocation(\%translocation, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}

# introduce random translocations
if (defined $translocation_count) {
  $local_time = localtime();
  print "[$local_time]\n";
  print "Introducing random Translocations based on the following parameters:\n";
  print "> translocation_count = $translocation_count\n";
  if (defined $centromere_gff) {
    print "> centromere_gff = $centromere_gff\n";
  }
  if (defined $translocation_breakpoint_gff) {
    print "> translocation_breakpoint_gff = $translocation_breakpoint_gff\n";
  }
  introduce_random_translocation($translocation_count, $centromere_gff, $translocation_breakpoint_gff, \%refseq, \%simseq, \%ref2sim_map);
  print "\n";
}


# generate output files
$local_time = localtime();
print "[$local_time]\n";
print "Simulation completed! :) \n\n";
$local_time = localtime();
print "[$local_time]\n";
print "Generating output files .. \n\n";
generate_output_files($prefix, \@refseq, \%refseq, \%simseq, \%ref2sim_map);
$local_time = localtime();
print "[$local_time]\n";
print "Done! :) \n\n";

sub read_file {
  my $file = shift @_;
  my $fh;
  if ($file =~ /\.gz$/) {
    open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
  } else {
    open($fh, $file) or die "can't open $file";
  }
  return $fh;
}

sub write_file {
  my $file = shift @_;
  my $fh;
  if ($file =~ /\.gz$/) {
    open($fh, "| gzip -c >$file") or die "can't open pipe to $file\n";
  } else {
    open($fh, ">$file") or die "can't open $file\n";
  }
  return $fh;
}

sub parse_list_file {
  my $fh = shift @_;
  my %list = ();
  while (<$fh>) {
    chomp;
    /^\s*$/ and next;
    /^#/ and next;
    if (exists $list{$_}) {
      $list{$_}++;
    } else {
      $list{$_} = 1;
    }
  }
  return %list;
}


sub parse_fasta_file {
  my ($fh, $input_hashref, $input_arrayref) = @_;
  my $seq_name = "";
  while (<$fh>) {
    chomp;
    if (/^\s*$/) {
      next;
    } elsif (/^\s*#/) {
      next;
    } elsif (/^>(.*)/) {
      $seq_name = $1;
      push @$input_arrayref, $seq_name;
      $$input_hashref{$seq_name} = "";
    } else {
      $$input_hashref{$seq_name} .= $_;
    }
  }
}

sub profile_base_freq {
  my $genome_hashref = shift @_;
  my @base = qw(A T G C);
  my %base_count = ();
  my %base_freq = ();
  foreach my $chr (sort keys %$genome_hashref) {
    my $seq = uc $$genome_hashref{$chr};
    foreach my $base (@base) {
      my $bc  = () = $seq =~ /$base/g;
      if (exists $base_count{$base}) {
        $base_count{$base} += $bc;
      } else {
        $base_count{$base} = $bc;
      }
    }
  }
  my $total_base_count = sum(values %base_count);
  foreach my $base (@base) {
    $base_freq{$base} = $base_count{$base}/$total_base_count;
  }
  return %base_freq;
}

sub create_genome_space {
  my $genome_hashref = shift @_;
  my %genome_space = ();
  my $offset = 0;
  foreach my $chr (sort keys %$genome_hashref) {
    my $chr_length = length $$genome_hashref{$chr};
    my $start = $offset + 1;
    my $end = $start + $chr_length - 1;
    # print "chr=$chr, chr_length=$chr_length, genome_space_start=$start, genome_space_end=$end\n";
    $genome_space{'chr-wide'}{$chr}{"start"} = $start;
    $genome_space{'chr-wide'}{$chr}{"end"} = $end;
    $genome_space{'chr-wide'}{$chr}{"length"} = $chr_length;
    $offset = $end;
    if (not exists $genome_space{'genome-wide'}) {
      $genome_space{'genome-wide'}{"start"} = $start;
      $genome_space{'genome-wide'}{"end"} = $end;
      $genome_space{'genome-wide'}{"length"} = $chr_length;
    } else {
      $genome_space{'genome-wide'}{"length"} += $chr_length;
      $genome_space{'genome-wide'}{"end"} = $genome_space{'genome-wide'}{"length"};
    }
  }
  return %genome_space;
}


sub parse_simple_vcf_file {
  my ($fh, $qual_cutoff, $query_type) = @_;
  my %vcf = ();
  while (<$fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($ref_chr, $ref_start, $variant_id, $ref_allele, $alt_allele, $variant_qual, $variant_filter, $variant_info) = split /\t/, $_;
    if (($variant_qual eq ".") or ($variant_qual >= $qual_cutoff)) {
      my $variant_type;
      my $ref_allele_length = length $ref_allele;
      my $alt_allele_length = length $alt_allele;
      my $ref_end = $ref_start + $ref_allele_length - 1;
      if ($alt_allele =~ /,/) {
        print "!!! Warning! Multiple alternative variants found at the same site:\n";
        print "!!! $ref_chr:$ref_start $ref_allele=>$alt_allele QUAL=$variant_qual!\n";
        print "!!! Ignore all variants at this site.\n\n";
      } else {
        if ($ref_allele_length ne $alt_allele_length) {
          $variant_type = "INDEL";
        } else {
          $variant_type = "SNP";
        }
        if ((defined $query_type) and ($query_type ne $variant_type)) {
          next;
        } else {
          my $check_overlap_flag = 0;
          if (exists $vcf{$ref_chr}) {
            if (exists $vcf{$ref_chr}{$ref_start}) {
              $check_overlap_flag = 1;
              print "!!! Warning! Multiple variants were defined within the same region: $ref_chr:$ref_start-$ref_end in the input vcf file!\n";
              print "!!! Only keep the first instance: $ref_chr:$ref_start $vcf{$ref_chr}{$ref_start}{'ref_allele'}=>$vcf{$ref_chr}{$ref_start}{'alt_allele'} QUAL=$vcf{$ref_chr}{$ref_start}{'variant_qual'}.\n";
              print "!!! Ignore the variant: $ref_chr:$ref_start $ref_allele => $alt_allele QUAL=$variant_qual.\n\n";
            } else {
              foreach my $s (sort {$a <=> $b} keys %{$vcf{$ref_chr}}) {
                $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $vcf{$ref_chr}{$s}{'ref_start'}, $vcf{$ref_chr}{$s}{'ref_end'});
                if ($check_overlap_flag == 1) {
                  print "!!! Warning! Multiple variants were defined within the same region: $ref_chr:$ref_start-$ref_end in the input vcf file!\n";
                  print "!!! Only keep the first instance: $ref_chr:$s $vcf{$ref_chr}{$s}{'ref_allele'}=>$vcf{$ref_chr}{$s}{'alt_allele'} QUAL=$vcf{$ref_chr}{$s}{'variant_qual'}.\n";
                  print "!!! Ignore the variant: $ref_chr:$ref_start $ref_allele => $alt_allele QUAL=$variant_qual.\n\n";
                  last;
                }
              }
            }
          }
          if ($check_overlap_flag == 0) {
            $vcf{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
            $vcf{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
            $vcf{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
            $vcf{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
            $vcf{$ref_chr}{$ref_start}{'alt_allele'} = $alt_allele;
            $vcf{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
            $vcf{$ref_chr}{$ref_start}{'variant_id'} = $variant_id;
            $vcf{$ref_chr}{$ref_start}{'variant_qual'} = $variant_qual;
            $vcf{$ref_chr}{$ref_start}{'variant_info'} = $variant_info;
          }
        }
      }
    }
  }
  return %vcf;
}

sub merge_vcf {
  my ($snp_vcf_hashref, $indel_vcf_hashref) = @_;
  my %merged_vcf = %$snp_vcf_hashref;
  foreach my $ref_chr (sort keys %$indel_vcf_hashref) {
    foreach my $ref_start (sort {$a <=> $b} keys %{$$indel_vcf_hashref{$ref_chr}}) {
      my $ref_start = $$indel_vcf_hashref{$ref_chr}{$ref_start}{'ref_start'};
      my $ref_end = $$indel_vcf_hashref{$ref_chr}{$ref_start}{'ref_end'};
      my $check_overlap_flag = 0;
      if (exists $$snp_vcf_hashref{$ref_chr}) {
        if (exists $$snp_vcf_hashref{$ref_chr}{$ref_start}) {
          $check_overlap_flag = 1;
          print "!!! Warning! Both SNP and INDEL variants were defined within the same region: $ref_chr:$ref_start-$ref_end in the input vcf files!\n";
          print "!!! Only keep the SNP variant: $$snp_vcf_hashref{$ref_chr}{$ref_start}{'ref_allele'}=>$$snp_vcf_hashref{$ref_chr}{$ref_start}{'alt_allele'}.\n";
          print "!!! Ignore the INDEL variant $$indel_vcf_hashref{$ref_chr}{$ref_start}{'ref_allele'}=>$$indel_vcf_hashref{$ref_chr}{$ref_start}{'alt_allele'} within this region.\n\n";

        } else {
          foreach my $s (sort {$a <=> $b} keys %{$$snp_vcf_hashref{$ref_chr}}) {
            $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $$snp_vcf_hashref{$ref_chr}{$s}{'ref_start'}, $$snp_vcf_hashref{$ref_chr}{$s}{'ref_end'});
            if ($check_overlap_flag == 1) {
              print "!!! Warning! Both SNP and INDEL variants were defined within the same region: $ref_chr:$ref_start-$ref_end in the input vcf files!\n";
              print "!!! Only keep the SNP variant: $$snp_vcf_hashref{$ref_chr}{$s}{'ref_allele'}=>$$snp_vcf_hashref{$ref_chr}{$s}{'alt_allele'}.\n";
              print "!!! Ignore the INDEL variant $$indel_vcf_hashref{$ref_chr}{$ref_start}{'ref_allele'}=>$$indel_vcf_hashref{$ref_chr}{$ref_start}{'alt_allele'} within this region.\n\n";
              last;
            }
          }
        }
      }
      if ($check_overlap_flag == 0) {
        $merged_vcf{$ref_chr}{$ref_start} = \%{$$indel_vcf_hashref{$ref_chr}{$ref_start}};
      }
    }
  }
  return %merged_vcf;
}

sub parse_gff_file {
  my $fh = shift @_;
  my %gff = ();
  while (<$fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $_;
    my ($id) = ($attribute =~ /ID=([^;]+)/);
    $gff{$id}{'chr'} = $chr;
    $gff{$id}{'feature'} = $feature;
    $gff{$id}{'start'} = $start;
    $gff{$id}{'end'} = $end;
    $gff{$id}{'strand'} = $strand;
    $gff{$id}{'id'} = $id;
    $gff{$id}{'attribute'} = $attribute;
  }
  return %gff;
}

sub introduce_defined_snp_indel {
  my ($vcf_hashref, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  my $snp_count = 0;
  my $indel_count = 0;
  my %offset = ();
  foreach my $ref_chr (sort keys %$vcf_hashref) {
    $offset{$ref_chr} = 0;
  }
  foreach my $ref_chr (sort keys %$vcf_hashref) {
    foreach my $ref_start (sort {$a <=> $b} keys %{$$vcf_hashref{$ref_chr}}) {
      my $variant_id = $$vcf_hashref{$ref_chr}{$ref_start}{'variant_id'};
      my $variant_type = $$vcf_hashref{$ref_chr}{$ref_start}{'variant_type'};
      my $ref_end = $$vcf_hashref{$ref_chr}{$ref_start}{'ref_end'};
      my $ref_allele = $$vcf_hashref{$ref_chr}{$ref_start}{'ref_allele'};
      my $alt_allele = $$vcf_hashref{$ref_chr}{$ref_start}{'alt_allele'};
      my $ref_allele_length = length $ref_allele;
      my $alt_allele_length = length $alt_allele;
      # print "ref_chr=$ref_chr, ref_start=$ref_start, ref_end=$ref_end, offset=$offset{$ref_chr}\n";
      # check if there are pre-introduced SNP/INDEL variants overlapping at the same site already
      my $check_overlap_flag = 0;
      if (exists  $$ref2sim_map_hashref{$ref_chr}) {
        if (exists $$ref2sim_map_hashref{$ref_chr}{$ref_start}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} =~ /(SNP|INDEL)/) {
            $check_overlap_flag = 1;
          }
        } else {
          foreach my $s (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
            if ($$ref2sim_map_hashref{$ref_chr}{$s}{'variant_type'} =~ /(SNP|INDEL)/) {
              $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $$ref2sim_map_hashref{$ref_chr}{$s}{'ref_start'}, $$ref2sim_map_hashref{$ref_chr}{$s}{'ref_end'});
              if ($check_overlap_flag == 1) {
                last;
              }
            }
          }
        }
      }
      if ($check_overlap_flag == 0) {
        substr $$simseq_hashref{$ref_chr}, $ref_start - 1 + $offset{$ref_chr}, $ref_allele_length, $$vcf_hashref{$ref_chr}{$ref_start}{'alt_allele'};

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $variant_id;

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_start + $alt_allele_length - 1 + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "+";

        if ($variant_type eq "SNP") {
          $snp_count++;
        } else {
          # INDEL
          my $indel_size = $alt_allele_length - $ref_allele_length;
          $offset{$ref_chr} += $indel_size;
          $indel_count++;
          if ($indel_size > 0) {
            $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'indel_type'} = 'INSERTION';
            $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'indel_size'} = $indel_size;
          } else {
            $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'indel_type'} = 'DELETION';
            $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'indel_size'} = -$indel_size;
          }
        }
      } else {
        print "!!! Warning! Multiple variants were defined within the same region: $ref_chr:$ref_start-$ref_end in the input vcf file(s)\n";
        print "!!! Only keep the first instance: $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} => $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'alt_allele'}\n";
        print "!!! Ignore all the other variants defined within this region.\n";
      }
    }
  }
  # print "> Introduced $snp_count SNP variants and $indel_count INDEL variants based on the input vcf file(s).\n";
}

sub check_overlap_region {
  my ($region1_start, $region1_end, $region2_start, $region2_end) = @_;
  # print "region1_start=$region1_start, region1_end=$region1_end, region2_start=$region2_start, region2_end=$region2_end\n";
  my $flag = 0;
  if (($region1_start >= $region2_start) and ($region1_start <= $region2_end)) {
    $flag = 1;
  }
  if (($region1_end >= $region2_start) and ($region1_end <= $region2_end)) {
    $flag = 1;
  }
  if (($region1_start <= $region2_start) and ($region1_end >= $region2_end)) {
    $flag = 1;
  }
  if (($region1_start >= $region2_start) and ($region1_end <= $region2_end)) {
    $flag = 1;
  }
  return $flag;
}

sub adjust_variant_coordinates_in_simseq {
  my ($ref_chr, $ref_start, $offset, $ref2sim_map_hashref) = @_;
  if (exists $$ref2sim_map_hashref{$ref_chr}) {
    foreach my $s (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
      if ($s > $ref_start) {
        # adjustment needed
        $$ref2sim_map_hashref{$ref_chr}{$s}{'sim_start'} += $offset;
        $$ref2sim_map_hashref{$ref_chr}{$s}{'sim_end'} += $offset;
      }
    }
  }
}

sub cal_prob_interval {
  my $prob_hashref = shift @_;
  my %prob_interval = ();
  my $lower_bound = 0;
  my $upper_bound = 0;
  foreach my $key (sort {$$prob_hashref{$a} <=> $$prob_hashref{$b}} keys %$prob_hashref){
    $upper_bound = $lower_bound + $$prob_hashref{$key};
    $prob_interval{$key} = "$lower_bound-$upper_bound";
    $lower_bound = $upper_bound;
  }
  return %prob_interval;
}

sub sample_from_interval {
  my $prob_interval_hashref = shift @_;
  my $dice = rand();
  foreach my $key (keys %$prob_interval_hashref){
    my ($lower, $upper) = split /-/, $$prob_interval_hashref{$key};
    if(($dice >= $lower) and ($dice < $upper)){
      return $key;
    }
  }
}

sub sample_inserted_seq {
  my ($indel_size, $refseq_hashref) = @_;
  my $inserted_seq = "";
  # profile the base composition of the reference genome
  my %refseq_base_freq = profile_base_freq($refseq_hashref);
  my %refseq_base_prob_interval = cal_prob_interval(\%refseq_base_freq);
  # print Dumper(%refseq_base_prob_interval);
  for(my $i = 1; $i <= $indel_size; $i++) {
    $inserted_seq .= sample_from_interval(\%refseq_base_prob_interval);
  }
  return $inserted_seq;
}

sub parse_snp_model {
  my $fh = shift @_;
  my %model = ();
  while (<$fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    if (/titv_ratio=(\S+)/) {
      $model{'titv_ratio'} = $1;
      if ($model{'titv_ratio'} eq "NA") {
        print "\n!!! Error! The supplied SNP model is incomplete: titv_ratio = NA\n";
        print "!!! Exit!\n";
        die;
      }
    } else {
      my ($original_base, $new_base, $freq) = ($_ =~ /(\w)\-\>(\w)\t(\S+)/);
      # print "original_base=$original_base, new_base=$new_base, freq=$freq\n";
      $model{'substitution_freq'}{$original_base}{$new_base} = $freq;
      if ($model{'substitution_freq'}{$original_base}{$new_base} eq "NA") {
        print "\n!!! Error! The supplied SNP model is incomplete: ${original_base}\-\>${new_base} substitution frequency = NA\n";
        print "!!! Exit!\n";
        die;
      }
      if (exists $model{'total_substitution_freq'}{$original_base}) {
        $model{'total_substitution_freq'}{$original_base} += $freq;
      } else {
        $model{'total_substitution_freq'}{$original_base} = $freq;
      }
    }
  }

  my @base = qw(A T G C);
  foreach my $original_base (@base) {
    foreach my $new_base (@base) {
      if ($original_base ne $new_base) {
        $model{'substitution_prob'}{$original_base}{$new_base} = $model{'substitution_freq'}{$original_base}{$new_base}/$model{'total_substitution_freq'}{$original_base};
      }
    }
  }
  # print(Dumper(%model));
  return %model;
}

sub pupy {
  my $base = shift @_;
  my $result;
  if ($base =~ /(A|a)/) {
    $result = "purine";
  } elsif ($base =~ /(G|g)/) {
    $result = "purine";
  } elsif ($base =~ /(T|t)/) {
    $result = "pyrimidine";
  } elsif ($base =~ /(C|c)/) {
    $result = "pyrimidine";
  } else {
    $result = "unknown";
  }
  return $result;
}

sub determine_substitution_probability_for_snp {
  my ($titv_ratio, $snp_model) = @_;
  my %substitution_prob = ();
  my @base = qw(A T G C);
  # set up titv_ratio and substitution_probability to prepare for the derived base sampling
  if (defined $snp_model) {
    # based on model
    # print "> snp_model = $snp_model\n\n";
    my $snp_model_fh = read_file($snp_model);
    my %snp_model = parse_snp_model($snp_model_fh);
    $titv_ratio = $snp_model{'titv_ratio'};
    %substitution_prob = %{$snp_model{'substitution_prob'}};
  } else {
    # random
    # print "> titv_ratio = $titv_ratio\n\n";
    foreach my $original_base (@base) {
      foreach my $new_base (@base) {
        if ($original_base ne $new_base) {
          if (pupy($original_base) eq pupy($new_base)) {
            $substitution_prob{$original_base}{$new_base} = $titv_ratio/($titv_ratio + 1);
          } else {
            $substitution_prob{$original_base}{$new_base} = 1/(($titv_ratio + 1) * 2);
          }
        }
      }
    }
  }
  return %substitution_prob;
}

sub sample_genome_space {
  my $genome_space_hashref = shift @_;
  my $sample = 1 + int rand($$genome_space_hashref{'genome-wide'}{'length'});
  # print "sub sample_genome_space >> sample = $sample\n";
  my ($chr_sampled, $start_sampled) = genome_space_translator($sample, $genome_space_hashref);
  # print "sub sample_genome_space >> chr_sampled = $chr_sampled, start_sampled = $start_sampled\n";
  return ($chr_sampled, $start_sampled);
}

sub genome_space_translator {
  my ($sample, $genome_space_hashref) = @_;
  my ($chr_translated, $pos_translated);
  foreach my $chr (sort keys %{$$genome_space_hashref{'chr-wide'}}) {
    my $chr_start = $$genome_space_hashref{'chr-wide'}{$chr}{'start'};
    my $chr_end = $$genome_space_hashref{'chr-wide'}{$chr}{'end'};
    if (($sample >= $chr_start) and ($sample <= $chr_end)) {
      $chr_translated = $chr;
      $pos_translated = $sample - $chr_start + 1;
      last;
    }
  }
  return ($chr_translated, $pos_translated);
}

sub sample_alt_allele {
  my ($ref_allele, $substitution_prob_hashref) = @_;
  # print "ref_allele = $ref_allele\n";
  # print "substitution_prob = \n";
  # print Dumper(%$substitution_prob_hashref);
  my %conditional_substitution_prob = %{$$substitution_prob_hashref{$ref_allele}};
  # print "conditional_substitution_prob = \n";
  # print Dumper(%conditional_substitution_prob);
  my %prob_interval = cal_prob_interval(\%conditional_substitution_prob);
  # print "prob_interval = \n";
  # print Dumper(%prob_interval);
  my $alt_allele = sample_from_interval(\%prob_interval);
  # print "alt_allele = $alt_allele\n";
  return $alt_allele;
}


sub introduce_random_snp {
  my ($snp_count, $snp_model, $titv_ratio, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  my %refseq_genome_space = create_genome_space($refseq_hashref);
  my %substitution_prob = determine_substitution_probability_for_snp($titv_ratio, $snp_model);
  # print Dumper(%substitution_prob);
  # sample and introduce SNP in the same time
  for (my $i = 1; $i <= $snp_count; $i++) {
    # sample genome space
    SAMPLE_RANDOM_SNP:
    my ($ref_chr, $ref_start) = sample_genome_space(\%refseq_genome_space);
    # print "chr_sampled = $ref_chr, start_sampled = $ref_start\n";
    my $ref_end = $ref_start;
    my $ref_allele = substr $$refseq_hashref{$ref_chr}, $ref_start - 1, 1;
    my $alt_allele;
    # check if the sampled position contains any ambiguous bases
    if ($ref_allele =~ /(N|n)/) {
      goto SAMPLE_RANDOM_SNP;
    } else {
      # check if there are pre-introduced SNP/INDEL variants overlapping at the same site already
      my $check_overlap_flag = 0;
      if (exists $$ref2sim_map_hashref{$ref_chr}) {
        if (exists $$ref2sim_map_hashref{$ref_chr}{$ref_start}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} =~ /(SNP|INDEL)/) {
            $check_overlap_flag = 1;
            goto SAMPLE_RANDOM_SNP;
          }
        } else {
          foreach my $s (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
            if ($$ref2sim_map_hashref{$ref_chr}{$s}{'variant_type'} =~ /(SNP|INDEL)/) {
              $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $$ref2sim_map_hashref{$ref_chr}{$s}{'ref_start'}, $$ref2sim_map_hashref{$ref_chr}{$s}{'ref_end'});
              if ($check_overlap_flag == 1) {
                goto SAMPLE_RANDOM_SNP;
              }
            }
          }
        }
      }
      if ($check_overlap_flag == 0) {
        # register this SNP
        $alt_allele = sample_alt_allele($ref_allele, \%substitution_prob);
        # print "ref_allele = $ref_allele, alt_allele = $alt_allele\n";
        substr $$simseq_hashref{$ref_chr}, $ref_start - 1, 1, $alt_allele;

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_end;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "+";

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = "SNP_${i}";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = "SNP";
      }
    }
  }
}


sub parse_indel_model {
  my $fh = shift @_;
  my %model = ();
  while (<$fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    if (/ins_del_ratio=(\S+)/) {
      $model{'ins_del_ratio'} = $1;
      if ($model{'ins_del_ratio'} eq "NA") {
        print "!!! Error! The supplied INDEL model is incomplete: ins_del_ratio = NA\n";
        print "!!! Exit!\n";
        die;
      }
    } else {
      my ($indel_size, $freq) = ($_ =~ /(\d+)\t(\S+)/);
      if ($freq eq "NA") {
        print "!!! Error! The supplied INDEL model is incomplete: the frequency of ${indel_size}-bp INDEL is NA\n";
        print "!!! Exit!\n";
        die;
      }
      $model{'indel_size_prob'}{$indel_size} = $freq;
    }
  }
  return %model;
}

sub sample_indel_type {
  my $ins_del_ratio = shift @_;
  my $sample = rand(1);
  if ($sample < $ins_del_ratio/(1 + $ins_del_ratio)) {
    $sample = "INSERTION";
  } else {
    $sample = "DELETION";
  }
  return $sample;
}

sub sample_indel_size {
  my $indel_size_prob_hashref = shift @_;
  my %indel_size_prob_interval = cal_prob_interval($indel_size_prob_hashref);
  my $indel_size = sample_from_interval(\%indel_size_prob_interval);
  return $indel_size;
}

sub introduce_random_indel {
  my ($indel_count, $indel_model, $ins_del_ratio, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  my %refseq_genome_space = create_genome_space($refseq_hashref);
  my %indel_prob = ();
  if (defined $indel_model) {
    # print "> indel_model = $indel_model\n\n";
    my $indel_model_fh = read_file($indel_model);
    %indel_prob = parse_indel_model($indel_model_fh);
  } else {
    # print "> ins_del_ratio = $ins_del_ratio\n\n";
    $indel_prob{'ins_del_ratio'} = $ins_del_ratio;
    # initialize default INDEL size distribution
    my $cumulative_prob = 0;
    for(my $indel_size = 50; $indel_size > 0; $indel_size--) {
      if ($indel_size > 1) {
        $indel_prob{'indel_size_prob'}{$indel_size} = $indel_size_powerlaw_constant * $indel_size **(-$indel_size_powerlaw_alpha);
        $cumulative_prob += $indel_prob{'indel_size_prob'}{$indel_size};
      } else {
        $indel_prob{'indel_size_prob'}{$indel_size} = 1 - $cumulative_prob;
      }
    }
  }

  # sample INDEL first
  my %indel_samples = ();
  for (my $i = 1; $i <= $indel_count; $i++) {
    # sample indel type and size
    my $indel_type = sample_indel_type($indel_prob{'ins_del_ratio'});
    my $indel_size = sample_indel_size(\%{$indel_prob{'indel_size_prob'}});
    # print "indel_type = $indel_type, indel_size = $indel_size\n";
    my $ref_allele;
    my $ref_allele_length;
    my $alt_allele;
    my $alt_allele_length;
    # sample genome space
    SAMPLE_RANDOM_INDEL:
    my ($ref_chr, $ref_start) = sample_genome_space(\%refseq_genome_space);
    my $ref_end;
    # print "chr_sampled = $ref_chr, start_sampled = $ref_start\n";
    if ($indel_type eq "INSERTION") {
      $ref_end = $ref_start;
      $ref_allele = substr $$refseq_hashref{$ref_chr}, $ref_start - 1, 1;
      $ref_allele_length = 1;
    } else {
      # Deletion
      # actual deletion starts from ref_start + 1
      $ref_end = $ref_start + $indel_size;
      # check if the deletion will go beyond the chromosome end
      if ($ref_end > $refseq_genome_space{'chr-wide'}{$ref_chr}{'length'}) {
        goto SAMPLE_RANDOM_INDEL;
      }
      $ref_allele = substr $$refseq_hashref{$ref_chr}, $ref_start - 1, $indel_size + 1;
      $ref_allele_length = $indel_size + 1;
    }
    # check if the sampled INDEL contains any ambiguous bases
    if ($ref_allele =~ /(N|n)/) {
      goto SAMPLE_RANDOM_INDEL;
    }
    # check if there are pre-introduced SNP/INDEL variants overlapping at the same site already
    my $check_overlap_flag = 0;
    if (exists $$ref2sim_map_hashref{$ref_chr}) {
      if (exists $$ref2sim_map_hashref{$ref_chr}{$ref_start}) {
        if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} =~ /(SNP|INDEL)/) {
          $check_overlap_flag = 1;
          goto SAMPLE_RANDOM_INDEL;
        }
      } else {
        foreach my $s (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$s}{'variant_type'} =~ /(SNP|INDEL)/) {
            $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $$ref2sim_map_hashref{$ref_chr}{$s}{'ref_start'}, $$ref2sim_map_hashref{$ref_chr}{$s}{'ref_end'});
            if ($check_overlap_flag == 1) {
              goto SAMPLE_RANDOM_INDEL;
            }
          }
        }
      }
    }
    # check if there are pre-sampled INDEL variants overlapping at the same site already
    if (exists $indel_samples{$ref_chr}) {
      if (exists $indel_samples{$ref_chr}{$ref_start}) {
        $check_overlap_flag = 1;
        goto SAMPLE_RANDOM_INDEL;
      }
    } else {
      foreach my $s (sort {$a <=> $b} keys %{$indel_samples{$ref_chr}}) {
        $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $indel_samples{$ref_chr}{$s}{'ref_start'}, $indel_samples{$ref_chr}{$s}{'ref_end'});
        if ($check_overlap_flag == 1) {
          goto SAMPLE_RANDOM_INDEL;
        }
      }
    }
    # all check passed, register this INDEL
    if ($indel_type eq "INSERTION") {
      $alt_allele = sample_inserted_seq($indel_size, $refseq_hashref);
      $alt_allele = $ref_allele . $alt_allele;
    } else {
      $alt_allele = substr $$refseq_hashref{$ref_chr}, $ref_start - 1, 1;
    }
    # print "ref_allele = $ref_allele, alt_allele = $alt_allele\n";
    # register this INDEL for later introduction.
    $indel_samples{$ref_chr}{$ref_start}{'indel_id'} = "INDEL_${i}";
    $indel_samples{$ref_chr}{$ref_start}{'indel_type'} = $indel_type;
    $indel_samples{$ref_chr}{$ref_start}{'indel_size'} = $indel_size;
    $indel_samples{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
    $indel_samples{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
    $indel_samples{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
    $indel_samples{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
    $indel_samples{$ref_chr}{$ref_start}{'alt_allele'} = $alt_allele;
  }

  # introduce all sampled INDELs
  my %offset = ();
  foreach my $ref_chr (sort keys %$refseq_hashref) {
    $offset{$ref_chr} = 0;
  }
  for my $ref_chr (sort keys %indel_samples) {
    for my $ref_start (sort {$a <=> $b} keys %{$indel_samples{$ref_chr}}) {
      my $indel_id = $indel_samples{$ref_chr}{$ref_start}{'indel_id'};
      my $indel_type = $indel_samples{$ref_chr}{$ref_start}{'indel_type'};
      my $indel_size = $indel_samples{$ref_chr}{$ref_start}{'indel_size'};
      my $ref_end = $indel_samples{$ref_chr}{$ref_start}{'ref_end'};
      my $ref_allele = $indel_samples{$ref_chr}{$ref_start}{'ref_allele'};
      my $alt_allele = $indel_samples{$ref_chr}{$ref_start}{'alt_allele'};
      my $ref_allele_length = length $ref_allele;
      my $alt_allele_length = length $alt_allele;
      substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, $ref_allele_length, $alt_allele;

      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";

      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start + $offset{$ref_chr};
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_start + $alt_allele_length - 1 + $offset{$ref_chr};
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "+";

      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = "INDEL";
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $indel_id;
      if ($indel_type eq "INSERTION") {
        $offset{$ref_chr} += $indel_size;
        adjust_variant_coordinates_in_simseq($ref_chr, $ref_start, $indel_size, $ref2sim_map_hashref);
      } else {
        $offset{$ref_chr} -= $indel_size;
        adjust_variant_coordinates_in_simseq($ref_chr, $ref_start, -$indel_size, $ref2sim_map_hashref);
      }

    }
  }
}

sub parse_sv_vcf_file {
  my $fh = shift @_;
  my %sv = ();
  while (<$fh>) {
    chomp;
    /^\s*$/ and next;
    /^#/ and next;
    my ($ref_chr, $ref_start, $id, $ref_allele, $alt_allele, $variant_qual, $variant_filter, $variant_info) = split /\t/, $_;
    if (($variant_info !~ /SVTYPE=/) or ($variant_info !~ /EVENT=/)) {
      print "!!! Error! The tags SVTYPE and EVENT are mandatory in the input vcf for defined CNV, inversion, or translocation variants!\n";
      print "!!! Exit!\n";
      die;
    } else {
      my ($sv_type) = ($variant_info =~ /SVTYPE=([^;]+)/);
      my ($sv_event) = ($variant_info =~ /EVENT=([^;]+)/);
      if ($sv_type eq "INV") {
        my $ref_end;
        if ($variant_info =~ /END=([^;]+)/) {
          $ref_end = $1;
        } else {
          print "!!! Error! The mandatory tag 'END=' has not been specified in the input vcf for the defined INV:\n";
          print "!!! $_\n";
          print "!!! Exit!\n";
          die;
        }
        $sv{$sv_event}{'ref_chr'} = $ref_chr;
        $sv{$sv_event}{'ref_start'} = $ref_start;
        $sv{$sv_event}{'ref_end'} = $ref_end;
        $sv{$sv_event}{'sv_type'} = $sv_type;
        $sv{$sv_event}{'sv_event'} = $sv_event;
      } elsif ($sv_type eq "DEL") {
        my $ref_end;
        if ($variant_info =~ /END=([^;]+)/) {
          $ref_end = $1;
        } else {
          print "!!! Error! The mandatory tag 'END=' has not been specified in the input vcf for the defined DEL:\n";
          print "!!! $_\n";
          print "!!! Exit!\n";
          die;
        }
        $sv{$sv_event}{'ref_chr'} = $ref_chr;
        $sv{$sv_event}{'ref_start'} = $ref_start;
        $sv{$sv_event}{'ref_end'} = $ref_end;
        $sv{$sv_event}{'sv_type'} = $sv_type;
        $sv{$sv_event}{'sv_event'} = $sv_event;
      } elsif ($sv_type eq "BND") {
        $sv{$sv_event}{'sv_type'} = "to_be_classified";
        $sv{$sv_event}{'sv_event'} = $sv_event;
        # see the secion 5 of VCFv4.1 specification (https://samtools.github.io/hts-specs/VCFv4.1.pdf) for the detailed meaning of s, t, and p used below.
        my $s = $ref_allele;
        my $t;
        my $p;
        my $p_relative_strand; # the relative strand of p relative to its original sequence
        my $p_relative_position; # the relative positon of p relative to t: "before_t" or "after_t"
        # print "alt_allele = $alt_allele\n";
        if ($alt_allele =~ /\[$/) {
          $p_relative_strand = "+";
          $p_relative_position = "after_ref_allele";
          ($t, $p) = ($alt_allele =~ /(\S+)\[(\S+)\[/);
        } elsif ($alt_allele =~ /\]$/) {
          $p_relative_strand = "-";
          $p_relative_position = "after_ref_allele";
          ($t, $p) = ($alt_allele =~ /(\S+)\](\S+)\]/);
        } elsif ($alt_allele =~ /^\]/) {
          $p_relative_strand = "+";
          $p_relative_position = "before_ref_allele";
          ($p, $t) = ($alt_allele =~ /\](\S+)\](\S+)/);
        } elsif ($alt_allele =~ /^\[/) {
          $p_relative_strand = "-";
          $p_relative_position = "before_ref_allele";
          ($p, $t) = ($alt_allele =~ /\[(\S+)\[(\S+)/);
        } else {
          print "Unexpected ALT field in the input vcf file for defined CNVs/inversions/translocation:\n";
          print "$_\n";
          print "Exit!\n";
          die;
        }
        # print "s=$s, t=$t, p=$p\n";
        if (not exists $sv{$sv_event}{'BND'}) {
          @{$sv{$sv_event}{'BND'}} = ();
        }
        my %bnd = ();
        $bnd{'ref_chr'} = $ref_chr;
        $bnd{'ref_start'} = $ref_start;
        $bnd{'s'} = $s;
        $bnd{'t'} = $t;
        $bnd{'p'} = $p;
        $bnd{'p_relative_strand'} = $p_relative_strand;
        $bnd{'p_relative_position'} = $p_relative_position;
        push @{$sv{$sv_event}{'BND'}}, \%bnd;
      }
    }
  }
  return %sv;
}

sub check_cnv_overlap {
  my ($cnv_hashref, $chr, $start, $end, $check_donor) = @_;
  my $flag = 0;
  # print "chr=$chr, start=$start, end=$end, check_donor=$check_donor\n";
  foreach my $c (sort keys %$cnv_hashref) {
    foreach my $s (sort {$a <=> $b} keys %{$$cnv_hashref{$c}}) {
      if ($chr eq $c) {
        $flag = check_overlap_region($start, $end, $$cnv_hashref{$c}{$s}{'ref_start'}, $$cnv_hashref{$c}{$s}{'ref_end'});
        if ($flag == 1) {
          return $flag;
        }
      }
      if ($check_donor eq "yes") {
        if ($$cnv_hashref{$c}{$s}{'variant_type'} eq "DUP") {
          if ($$cnv_hashref{$c}{$s}{'donor_chr_in_ref'} eq $chr) {
            $flag = check_overlap_region($start, $end, $$cnv_hashref{$c}{$s}{'donor_start_in_ref'}, $$cnv_hashref{$c}{$s}{'donor_end_in_ref'});
            if ($flag == 1){
              return $flag;
            }
          }
        }
      }
    }
  }
  return $flag;
}

sub extract_cnv_from_sv {
  my $sv_hashref = shift @_;
  my %cnv = ();
  foreach my $sv_event (sort keys %$sv_hashref) {
    my $check_overlap_flag = 0;
    if ($$sv_hashref{$sv_event}{'sv_type'} eq "DEL") {
      my $ref_chr = $$sv_hashref{$sv_event}{'ref_chr'};
      my $ref_start = $$sv_hashref{$sv_event}{'ref_start'};
      my $ref_end = $$sv_hashref{$sv_event}{'ref_end'}; # ref_end - ref_start = deletion_size
      # check overlap with pre-defined CNVs
      my $check_overlap_flag = check_cnv_overlap(\%cnv, $ref_chr, $ref_start, $ref_end, "no");
      if ($check_overlap_flag == 1) {
        print "!!! Warning! Multiple overlapped CNVs found within the region $ref_chr:$ref_start-$ref_end!\n";
        print "!!! Only keep the first instance and ignore the others\n";
        next;
      } else {
        # register this CNV
        $cnv{$ref_chr}{$ref_start}{'cnv_id'} = $sv_event;
        $cnv{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $cnv{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $cnv{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $cnv{$ref_chr}{$ref_start}{'deletion_start'} = $ref_start + 1;
        $cnv{$ref_chr}{$ref_start}{'deletion_end'} = $ref_end;
        $cnv{$ref_chr}{$ref_start}{'deletion_size'} = $ref_end - $ref_start;
        $cnv{$ref_chr}{$ref_start}{'variant_type'} = "DEL";
        $cnv{$ref_chr}{$ref_start}{'cnv_type'} = "copy_number_loss";
      }
    } elsif ($$sv_hashref{$sv_event}{'sv_type'} eq "to_be_classified") {
      my $bnd_count = scalar @{$$sv_hashref{$sv_event}{'BND'}};
      if ($bnd_count == 2) {
        my %bnd = ();
        foreach my $b_hashref (sort @{$$sv_hashref{$sv_event}{'BND'}}) {
          my $s = $$b_hashref{'s'};
          my $p = $$b_hashref{'p'};
          my $t = $$b_hashref{'t'};
          my ($bnd_p_chr, $bnd_p_start) = split /:/, $p;
          my $p_relative_strand = $$b_hashref{'p_relative_strand'};
          my $p_relative_position = $$b_hashref{'p_relative_position'};
          if ($p_relative_position eq "after_ref_allele") {
            $bnd{'1'} = $b_hashref;
            $bnd{'1'}{'p_chr'} = $bnd_p_chr;
            $bnd{'1'}{'p_start'} = $bnd_p_start;
          } else {
            $bnd{'2'} = $b_hashref;
            $bnd{'2'}{'p_chr'} = $bnd_p_chr;
            $bnd{'2'}{'p_start'} = $bnd_p_start;
          }
        }
        # print Dumper(%bnd);
        # print "\n";
        # verify this is indeed a segmental duplication
        if ($bnd{'1'}{'ref_chr'} eq $bnd{'2'}{'ref_chr'}) {
          my $recipient_chr = $bnd{'1'}{'ref_chr'};
          my ($recipient_start, $recipient_end) = ($bnd{'1'}{'ref_start'}, $bnd{'2'}{'ref_start'});
          if ($recipient_end - $recipient_start == 1) {
            if ($bnd{'1'}{'p_chr'} eq $bnd{'2'}{'p_chr'}) {
              my $donor_chr_in_ref = $bnd{'1'}{'p_chr'};
              my ($donor_start_in_ref, $donor_end_in_ref) = sort {$a <=> $b} ($bnd{'1'}{'p_start'}, $bnd{'2'}{'p_start'});
              my $donor_strand_in_ref = $bnd{'1'}{'p_relative_strand'};
              # check overlap with pre-defined CNVs
              my $check_overlap_flag1 = check_cnv_overlap(\%cnv, $recipient_chr, $recipient_start, $recipient_end, "yes");
              my $check_overlap_flag2 = check_cnv_overlap(\%cnv, $donor_chr_in_ref, $donor_start_in_ref, $donor_end_in_ref, "no");
              if (($check_overlap_flag1 == 0) and ($check_overlap_flag2 == 0)) {
                $cnv{$recipient_chr}{$recipient_start}{'cnv_id'} = $$sv_hashref{$sv_event}{'sv_event'};
                $cnv{$recipient_chr}{$recipient_start}{'ref_chr'} = $recipient_chr;
                $cnv{$recipient_chr}{$recipient_start}{'ref_start'} = $recipient_start; #
                $cnv{$recipient_chr}{$recipient_start}{'ref_end'} = $recipient_start; # same as ref_start
                $cnv{$recipient_chr}{$recipient_start}{'donor_chr_in_ref'} = $donor_chr_in_ref;
                $cnv{$recipient_chr}{$recipient_start}{'donor_start_in_ref'} = $donor_start_in_ref;
                $cnv{$recipient_chr}{$recipient_start}{'donor_end_in_ref'} = $donor_end_in_ref;
                $cnv{$recipient_chr}{$recipient_start}{'donor_strand_in_ref'} = $donor_strand_in_ref;
                $cnv{$recipient_chr}{$recipient_start}{'duplication_size'} = $donor_end_in_ref - $donor_start_in_ref + 1;
                $cnv{$recipient_chr}{$recipient_start}{'variant_type'} = "DUP";
                $cnv{$recipient_chr}{$recipient_start}{'cnv_type'} = "copy_number_gain";
              }
            }
          }
        }
      }
    }
  }
  return %cnv;
}

sub introduce_defined_cnv {
  my ($cnv_hashref, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  my $cnv_count = 0;
  my %offset = ();
  foreach my $ref_chr (sort keys %$cnv_hashref) {
    $offset{$ref_chr} = 0;
  }
  foreach my $ref_chr (sort keys %$cnv_hashref) {
    foreach my $ref_start (sort {$a <=> $b} keys %{$$cnv_hashref{$ref_chr}}) {
      my $cnv_id = $$cnv_hashref{$ref_chr}{$ref_start}{'cnv_id'};
      my $variant_type = $$cnv_hashref{$ref_chr}{$ref_start}{'variant_type'};
      my $cnv_type = $$cnv_hashref{$ref_chr}{$ref_start}{'cnv_type'};
      if ($variant_type eq "DEL") {
        my $ref_start = $$cnv_hashref{$ref_chr}{$ref_start}{'ref_start'};
        my $ref_end = $$cnv_hashref{$ref_chr}{$ref_start}{'ref_end'};
        my $deletion_start = $$cnv_hashref{$ref_chr}{$ref_start}{'deletion_start'};
        my $deletion_end = $$cnv_hashref{$ref_chr}{$ref_start}{'deletion_end'};
        my $deletion_size = $$cnv_hashref{$ref_chr}{$ref_start}{'deletion_size'};
        # standard notation for the reference allele
        # my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, $deletion_size + 1;
        # short-handed notation for the reference allele
        my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # standard notation for the alternative allele
        # my $alt_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # short-handed notation for the alternative allele
        my $alt_allele = "<DEL>";
        substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr}, $deletion_size, "";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_start  + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "+";

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $cnv_id;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'cnv_type'} = "copy_number_loss";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'deletion_size'} = $deletion_size;
        $offset{$ref_chr} -= $deletion_size;
      } else {
        # cnv_type eq "DUP"
        my $donor_chr_in_ref = $$cnv_hashref{$ref_chr}{$ref_start}{'donor_chr_in_ref'};
        my $donor_start_in_ref = $$cnv_hashref{$ref_chr}{$ref_start}{'donor_start_in_ref'};
        my $donor_end_in_ref = $$cnv_hashref{$ref_chr}{$ref_start}{'donor_end_in_ref'};
        my $donor_strand_in_ref = $$cnv_hashref{$ref_chr}{$ref_start}{'donor_strand_in_ref'};
        my $duplication_size = $$cnv_hashref{$ref_chr}{$ref_start}{'duplication_size'};
        my $duplication_seq = substr $$refseq_hashref{$donor_chr_in_ref}, $donor_start_in_ref - 1, $duplication_size;
        # test the strand
        if ($donor_strand_in_ref eq "-") {
          $duplication_seq = revcom($duplication_seq);
        }
        # standard
        # standard notation for the reference allele
        # my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # short-handed notation for the reference allele
        my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # standard notation for the alternative allele
        # my $alt_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # $alt_allele = $alt_allele . $duplication_seq;
        # short-handed notation for the alternative allele
        my $alt_allele = "<DUP>";
        substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr}, 0, $duplication_seq;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_start + $duplication_size + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "+";

        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $cnv_id;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'cnv_type'} = "copy_number_gain";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'duplication_size'} = $duplication_size;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_chr_in_ref'} = $donor_chr_in_ref;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_start_in_ref'} = $donor_start_in_ref;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_end_in_ref'} = $donor_end_in_ref;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_strand_in_ref'} = $donor_strand_in_ref;
        $offset{$ref_chr} += $duplication_size;
      }
      $cnv_count++;
    }
  }
}

sub sample_cnv_type {
  my $cnv_gain_loss_ratio = shift @_;
  my $sample = rand(1);
  if ($sample < $cnv_gain_loss_ratio/(1 + $cnv_gain_loss_ratio)) {
    $sample = "copy_number_gain";
  } else {
    $sample = "copy_number_loss";
  }
  return $sample;
}

sub sample_cnv_size {
  my ($cnv_min_size, $cnv_max_size) = @_;
  my $cnv_size = $cnv_min_size + int rand($cnv_max_size - $cnv_min_size);
  return $cnv_size;
}

sub sample_cnv_copy_number {
  my $cnv_max_copy_number = shift @_;
  my $cnv_copy_number = 2 + int rand($cnv_max_copy_number - 1);
  return $cnv_copy_number;
}

sub sample_strand {
  my $forward2reverse_strand_ratio = shift @_;
  my $sample = rand(1);
  if ($sample < $forward2reverse_strand_ratio/(1 + $forward2reverse_strand_ratio)) {
    $sample = "+";
  } else {
    $sample = "-";
  }
  return $sample;
}

sub introduce_random_cnv {
  my ($cnv_count, $cnv_gain_loss_ratio, $cnv_min_size, $cnv_max_size, $cnv_max_copy_number, $centromere_gff, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  my %refseq_genome_space = create_genome_space($refseq_hashref);
  my $forward2reverse_strand_ratio = 1; # the ratio of forward vs. inverted inserted copies for DUP. Default = 1.0 (i.e. equal chance for the two possible inserted orientations).
  # define centromeres for SV simulation
  my %centromere = ();
  if (defined $centromere_gff) {
    my $centromere_gff_fh = read_file($centromere_gff);
    %centromere = parse_gff_file($centromere_gff_fh);
  }

  # sample CNV first
  my %cnv_samples = ();
  for (my $i = 1; $i <= $cnv_count; $i++) {
    # sample CNV type
    my $cnv_type = sample_cnv_type($cnv_gain_loss_ratio);
    # print "i = $i, cnv_type = $cnv_type\n";
    if ($cnv_type eq "copy_number_loss") {
      my $variant_type = "DEL";
      # sample deletion size
      my $deletion_size = sample_cnv_size($cnv_min_size, $cnv_max_size);
      # print "deletion_size = $deletion_size\n";
      # sample the deleted region from the genome space
      my $chr_end_margin = 1000;
      SAMPLE_RANDOM_DEL:
      # sample the deleted region (for copy number loss)
      my ($ref_chr, $ref_pos) = sample_genome_space(\%refseq_genome_space);
      my $ref_start = $ref_pos - 1;  # including the base immediately before the sampled region
      if ($ref_start < $chr_end_margin) {
	  goto SAMPLE_RANDOM_DEL;
      }
      my $ref_end = $ref_start + $deletion_size;
      # check if the sampled end position will go beyond the chromosome end
      if ($ref_end > ($refseq_genome_space{'chr-wide'}{$ref_chr}{'length'} - $chr_end_margin)) {
        goto SAMPLE_RANDOM_DEL;
      }
      # check if the sampled region overlapped with the defined centromeres
      if (defined $centromere_gff) {
        foreach my $centromere_id (sort keys %centromere) {
          if ($centromere{$centromere_id} eq $ref_chr) {
            my $centromere_check_flag = check_overlap_region($ref_start, $ref_end, $centromere{$centromere_id}{'start'}, $ref_end, $centromere{$centromere_id}{'end'});
            if ($centromere_check_flag == 1) {
              goto SAMPLE_RANDOM_DEL;
            }
          }
        }
      }
      my $ref_allele = substr $$refseq_hashref{$ref_chr}, $ref_start - 1, $deletion_size + 1;
      # check if the sampled region hit the uncertain part of the reference genome
      if ($ref_allele =~ /(N|n)/) {
        goto SAMPLE_RANDOM_DEL;
      }
      # check if there are pre-sampled CNV variants overlapping at the same site already
      my $check_overlap_flag = check_cnv_overlap(\%cnv_samples, $ref_chr, $ref_start, $ref_end, "yes");
      if ($check_overlap_flag == 1) {
        goto SAMPLE_RANDOM_DEL;
      } else {
        # all check passed, register this DEL
        # print "register sampled DEL at $ref_chr:$ref_start-$ref_end, deletion_size = $deletion_size\n";
        $cnv_samples{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $cnv_samples{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $cnv_samples{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $cnv_samples{$ref_chr}{$ref_start}{'cnv_id'} = "CNV_${i}.1";
        $cnv_samples{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
        $cnv_samples{$ref_chr}{$ref_start}{'cnv_type'} = $cnv_type;
        $cnv_samples{$ref_chr}{$ref_start}{'deletion_size'} = $deletion_size;
      }
    } else {
      my $variant_type = "DUP";
      my $duplication_size = sample_cnv_size($cnv_min_size, $cnv_max_size);
      # print "duplication_size = $duplication_size\n";
      # sample the donor region from the genome space
      SAMPLE_RANDOM_DONOR:
      my ($donor_chr_in_ref, $donor_start_in_ref) = sample_genome_space(\%refseq_genome_space);
      my $donor_end_in_ref = $donor_start_in_ref + $duplication_size - 1;
      # check if the sampled end position will go beyond the chromosome end
      if ($donor_end_in_ref > $refseq_genome_space{'chr-wide'}{$donor_chr_in_ref}{'length'}) {
        goto SAMPLE_RANDOM_DONOR;
      }
      # check if the sampled region overlapped with the defined centromeres
      if (defined $centromere_gff) {
        foreach my $centromere_id (sort keys %centromere) {
          if ($centromere{$centromere_id}{'chr'} eq $donor_chr_in_ref) {
            my $centromere_check_flag = check_overlap_region($donor_start_in_ref, $donor_end_in_ref, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
            if ($centromere_check_flag == 1) {
              goto SAMPLE_RANDOM_DONOR;
            }
          }
        }
      }
      my $duplication_seq = substr $$refseq_hashref{$donor_chr_in_ref}, $donor_start_in_ref - 1, $duplication_size + 1;
      # check if the sampled region hit the uncertain part of the reference genome
      if ($duplication_seq =~ /(N|n)/) {
        goto SAMPLE_RANDOM_DONOR;
      }
      # check if there are pre-sampled CNV variants overlapping at the same site already
      my $check_overlap_flag = check_cnv_overlap(\%cnv_samples, $donor_chr_in_ref, $donor_start_in_ref, $donor_end_in_ref, "no");
      if ($check_overlap_flag == 1) {
        goto SAMPLE_RANDOM_DONOR;
      } else {
        # donor sampling is completed
        # print "sampled donor: donor_chr_in_ref = $donor_chr_in_ref, donor_start_in_ref = $donor_start_in_ref, donor_end_in_ref = $donor_end_in_ref\n";
        # now sample inserted copy number
        my $cnv_extra_copy_number = sample_cnv_copy_number($cnv_max_copy_number) -  1;
        for (my $j = 1; $j <= $cnv_extra_copy_number; $j++) {
          # sample the inserted location for the donor sequence from the genome space
          SAMPLE_RANDOM_DUP_INSERT:
          my ($inserted_chr_in_ref, $inserted_start_in_ref) = sample_genome_space(\%refseq_genome_space);
          my $inserted_end_in_ref = $inserted_start_in_ref;
          # print "inserted_chr_in_ref = $inserted_chr_in_ref, inserted_start_in_ref = $inserted_start_in_ref\n";
          # check if the sampled region overlapped with the defined centromere
          if (defined $centromere_gff) {
            foreach my $centromere_id (sort keys %centromere) {
              if ($centromere{$centromere_id}{'chr'} eq $inserted_chr_in_ref) {
                my $centromere_check_flag = check_overlap_region($inserted_start_in_ref, $inserted_end_in_ref, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                if ($centromere_check_flag == 1) {
                  goto SAMPLE_RANDOM_DUP_INSERT;
                }
              }
            }
          }
          # check if there are pre-sampled CNV variants overlapping at the same site already
          my $check_overlap_flag2 = check_cnv_overlap(\%cnv_samples, $inserted_chr_in_ref, $inserted_start_in_ref, $inserted_end_in_ref, "yes");
          if ($check_overlap_flag2 == 1) {
            goto SAMPLE_RANDOM_DUP_INSERT;
          } else {
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'ref_chr'} = $inserted_chr_in_ref;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'ref_start'} = $inserted_start_in_ref;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'ref_end'} = $inserted_end_in_ref;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'donor_chr_in_ref'} = $donor_chr_in_ref;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'donor_start_in_ref'} = $donor_start_in_ref;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'donor_end_in_ref'} = $donor_end_in_ref;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'donor_strand_in_ref'} = sample_strand($forward2reverse_strand_ratio);
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'cnv_id'} = "CNV_"."$i.$j";
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'variant_type'} = $variant_type;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'cnv_type'} = $cnv_type;
            $cnv_samples{$inserted_chr_in_ref}{$inserted_start_in_ref}{'duplication_size'} = $duplication_size;
          }
        }
      }
    }
  }

  # introduce all sampled CNVs
  my %offset = ();
  foreach my $ref_chr (sort keys %$refseq_hashref) {
    $offset{$ref_chr} = 0;
  }
  for my $ref_chr (sort keys %cnv_samples) {
    for my $ref_start (sort {$a <=> $b} keys %{$cnv_samples{$ref_chr}}) {
      my $cnv_id = $cnv_samples{$ref_chr}{$ref_start}{'cnv_id'};
      my $variant_type = $cnv_samples{$ref_chr}{$ref_start}{'variant_type'};
      my $cnv_type = $cnv_samples{$ref_chr}{$ref_start}{'cnv_type'};
      my $ref_end = $cnv_samples{$ref_chr}{$ref_start}{'ref_end'};
      if ($variant_type eq "DEL") {
        my $deletion_start = $cnv_samples{$ref_chr}{$ref_start}{'deletion_start'};
        my $deletion_end = $cnv_samples{$ref_chr}{$ref_start}{'deletion_end'};
        my $deletion_size = $cnv_samples{$ref_chr}{$ref_start}{'deletion_size'};
        # standard notation for the reference allele
        # my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, $deletion_size + 1;
        # short-handed notation for the reference allele
        my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # standard notation for the alternative allele
        # my $alt_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # short-handed notation for the alternative allele
        my $alt_allele = "<DEL>";
        substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr}, $deletion_size, "";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_start  + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "+";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $cnv_id;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'cnv_type'} = $cnv_type;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'deletion_size'} = $deletion_size;
        $offset{$ref_chr} -= $deletion_size;
      } else {
        # variant_type eq "DUP"
        my $donor_chr_in_ref = $cnv_samples{$ref_chr}{$ref_start}{'donor_chr_in_ref'};
        my $donor_start_in_ref = $cnv_samples{$ref_chr}{$ref_start}{'donor_start_in_ref'};
        my $donor_end_in_ref = $cnv_samples{$ref_chr}{$ref_start}{'donor_end_in_ref'};
        my $donor_strand_in_ref = $cnv_samples{$ref_chr}{$ref_start}{'donor_strand_in_ref'};
        my $duplication_size = $cnv_samples{$ref_chr}{$ref_start}{'duplication_size'};
        my $duplication_seq = substr $$refseq_hashref{$donor_chr_in_ref}, $donor_start_in_ref - 1, $duplication_size;
        # test the orientation
        if ($donor_strand_in_ref eq "-") {
          $duplication_seq = revcom($duplication_seq);
        }
        # standard
        # standard notation for the reference allele
        # my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # short-handed notation for the reference allele
        my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # standard notation for the alternative allele
        # my $alt_allele = substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr} - 1, 1;
        # $alt_allele = $alt_allele . $duplication_seq;
        # short-handed notation for the alternative allele
        my $alt_allele = "<DUP>";
        substr $$simseq_hashref{$ref_chr}, $ref_start + $offset{$ref_chr}, 0, $duplication_seq;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_start + $duplication_size + $offset{$ref_chr};
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "+";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $cnv_id;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'cnv_type'} = "copy_number_gain";
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'duplication_size'} = $duplication_size;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_chr_in_ref'} = $donor_chr_in_ref;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_start_in_ref'} = $donor_start_in_ref;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_end_in_ref'} = $donor_end_in_ref;
        $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_strand_in_ref'} = $donor_strand_in_ref;
        $offset{$ref_chr} += $duplication_size;
      }
    }
  }
}


sub sample_inversion_size {
  my ($inv_min_size, $inv_max_size) = @_;
  my $inv_size = $inv_min_size + int rand($inv_max_size - $inv_min_size);
  return $inv_size;
}

sub extract_inversion_from_sv {
  my $sv_hashref = shift @_;
  my %inversion = ();
  foreach my $sv_event (sort keys %$sv_hashref) {
    my $check_overlap_flag = 0;
    if ($$sv_hashref{$sv_event}{'sv_type'} eq "INV") {
      my $ref_chr = $$sv_hashref{$sv_event}{'ref_chr'};
      my $ref_start = $$sv_hashref{$sv_event}{'ref_start'};
      my $ref_end = $$sv_hashref{$sv_event}{'ref_end'};
      my $inversion_id = $$sv_hashref{$sv_event}{'sv_event'};
      # check if there are pre-sampled SV variants overlapping at the same site already
      my $check_overlap_flag = 0;
      if (exists $inversion{$ref_chr}) {
        if (exists $inversion{$ref_chr}{$ref_start}) {
          $check_overlap_flag = 1;
        }
      } else {
        foreach my $s (sort {$a <=> $b} keys %{$inversion{$ref_chr}}) {
          $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $inversion{$ref_chr}{$s}{'ref_start'}, $inversion{$ref_chr}{$s}{'ref_end'});
          if ($check_overlap_flag == 1) {
            last;
          }
        }
      }
      if ($check_overlap_flag == 1) {
        print "!!! Warning! Multiple overlapped inversion found within $ref_chr:$ref_start-$ref_end!\n";
        print "!!! Only keep the first one and ignore the others\n";
        next;
      } else {
        # register this inversion
        $inversion{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $inversion{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $inversion{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $inversion{$ref_chr}{$ref_start}{'inversion_id'} = $inversion_id;
      }
    }
  }
  return %inversion;
}

sub introduce_defined_inversion {
  my ($inversion_hashref, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  # print Dumper(%$inversion_hashref);
  foreach my $ref_chr (sort keys %$inversion_hashref) {
    foreach my $ref_start (sort {$a <=> $b} keys %{$$inversion_hashref{$ref_chr}}) {
      my $inversion_id = $$inversion_hashref{$ref_chr}{$ref_start}{'inversion_id'};
      my $ref_end = $$inversion_hashref{$ref_chr}{$ref_start}{'ref_end'};
      my $inversion_size = $ref_end - $ref_start + 1;
      if ($ref_end > (length $$refseq_hashref{$ref_chr})) {
	  print "!!! Error! The end of the defined inversion $inversion_id goes beyond the end of the corresponding chromosome!\n";
	  print "!!! Exit!\n";
	  die;
      }
      my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start - 1, 1;
      my $alt_allele = "<INV>";
      my $inversion_seq = substr $$simseq_hashref{$ref_chr}, $ref_start - 1, $inversion_size;
      $inversion_seq = revcom($inversion_seq);
      substr $$simseq_hashref{$ref_chr}, $ref_start - 1, $inversion_size, $inversion_seq;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_end;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "-";
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = "INV";
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $inversion_id;
      $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'inversion_size'} = $inversion_size;
    }
  }
}

sub examine_breakpoint_for_inversion {
  my ($breakpoint_hashref, $centromere_gff) = @_;
  my %breakpoint_by_chr = ();
  my %valid_breakpoint_pair = ();
  foreach my $id (sort keys %$breakpoint_hashref) {
    my $chr = $$breakpoint_hashref{$id}{'chr'};
    my $feature = $$breakpoint_hashref{$id}{'feature'};
    $breakpoint_by_chr{$chr}{$feature}{$id}{'start'} = $$breakpoint_hashref{$id}{'start'};
    $breakpoint_by_chr{$chr}{$feature}{$id}{'end'} = $$breakpoint_hashref{$id}{'end'};
    $breakpoint_by_chr{$chr}{$feature}{$id}{'strand'} = $$breakpoint_hashref{$id}{'strand'};
  }
  foreach my $chr (sort keys %breakpoint_by_chr) {
    foreach my $feature (sort keys %{$breakpoint_by_chr{$chr}}) {
      my @breakpoint_on_chr = sort keys %{$breakpoint_by_chr{$chr}{$feature}};
      if ((scalar @breakpoint_on_chr) >= 2) {
        my @breakpoint_on_chr_positive_strand = ();
        my @breakpoint_on_chr_negative_strand = ();
        foreach my $b (@breakpoint_on_chr) {
          if ($breakpoint_by_chr{$chr}{$feature}{$b}{'strand'} eq "+") {
            push @breakpoint_on_chr_positive_strand, $b;
          } elsif ($breakpoint_by_chr{$chr}{$feature}{$b}{'strand'} eq "-") {
            push @breakpoint_on_chr_negative_strand, $b;
          } else {
            print "\n!!! Error! Strand was not defined for the breakpoint $b in the input inversion_breakpoint_gff\n";
            print "!!! Exit!\n";
            die;
          }
        }
        if (((scalar @breakpoint_on_chr_positive_strand) > 0) and ((scalar @breakpoint_on_chr_negative_strand) > 0)) {
          foreach my $b1 (@breakpoint_on_chr_positive_strand) {
            foreach my $b2 (@breakpoint_on_chr_negative_strand){
              my $inv_start;
              my $inv_end;
              if ($$breakpoint_hashref{$b1}{'start'} > $$breakpoint_hashref{$b2}{'end'}) {
                ($inv_start, $inv_end) = ($$breakpoint_hashref{$b2}{'start'}, $$breakpoint_hashref{$b1}{'end'});
              } else {
                ($inv_start, $inv_end) = ($$breakpoint_hashref{$b1}{'start'}, $$breakpoint_hashref{$b2}{'end'});
              }
              my $centromere_check_flag = 0;
              if (defined $centromere_gff) {
                my $centromere_gff_fh = read_file($centromere_gff);
                my %centromere = parse_gff_file($centromere_gff_fh);
                close $centromere_gff_fh;
                foreach my $centromere_id (sort keys %centromere) {
                  if ($centromere{$centromere_id}{'chr'} eq $chr) {
                    my $centromere_check_flag = check_overlap_region($inv_start, $inv_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                    if ($centromere_check_flag == 1) {
                      last;
                    }
                  }
                }
              }
              if ($centromere_check_flag == 0) {
                $valid_breakpoint_pair{'+-'}{$b1}{$b2} = 1;
                $valid_breakpoint_pair{'-+'}{$b2}{$b1} = 1;
              }
            }
          }
        }
      }
    }
  }
  if (not exists $valid_breakpoint_pair{'+-'}) {
    print "\n!!! Error! None of the defined breakpoints is valid for triggering inversions!\n";
    print "!!! Valid breakpoints should satisfy the following creteria:\n";
    print "!!! 1) The two breakpoints should come from the same chromosome.\n";
    print "!!! 2) The two breakpoints should come from opposite strands.\n";
    print "!!! 3) If the centromere for the corresponding chromosome has been defined, the two breakpoints should not enclose this centromere.\n";
    print "!!! Exit!\n";
    die;
  } else {
    return %valid_breakpoint_pair;
  }
}


sub introduce_random_inversion {
  my ($inversion_count, $inversion_min_size, $inversion_max_size, $centromere_gff, $inversion_breakpoint_gff, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  my %centromere = ();
  if (defined $centromere_gff) {
    my $centromere_gff_fh = read_file($centromere_gff);
    %centromere = parse_gff_file($centromere_gff_fh);
    close $centromere_gff_fh;
  }
  my %inversion = ();
  # print "sample inversion: i=$i\n";
  if (defined $inversion_breakpoint_gff) {
    my $inversion_breakpoint_gff_fh = read_file($inversion_breakpoint_gff);
    my %breakpoint = parse_gff_file($inversion_breakpoint_gff_fh);
    my %valid_breakpoint_pair = examine_breakpoint_for_inversion(\%breakpoint, $centromere_gff);
    for (my $i = 1; $i <= $inversion_count; $i++) {
      SAMPLE_RANDOM_INV1:
      # sample inversion based on defined breakpoints
      my $valid_breakpoint_pair_count = scalar (keys %{$valid_breakpoint_pair{'+-'}});
      if ($valid_breakpoint_pair_count < 1) {
        my $j = $i - 1;
        print "\n!!! Warning! No more valid breakpoint pairs can be found in the current simulation based on the defined breakpoint file: $inversion_breakpoint_gff\n";
        print "!!! Only $j inversions were introduced\n";
        last;
      } else {
        my @breakpoint1 = shuffle(sort keys %{$valid_breakpoint_pair{'+-'}});
        my $breakpoint1 = shift @breakpoint1;
        my @breakpoint2 = shuffle(sort keys %{$valid_breakpoint_pair{'+-'}{$breakpoint1}});
        my $breakpoint2 = shift @breakpoint2;
        # print "breakpoint1 = $breakpoint1\n";
        # print "breakpoint2 = $breakpoint2\n";
        my $ref_chr = $breakpoint{$breakpoint1}{'chr'};
        my $ref_start;
        my $ref_end;
        if ($breakpoint{$breakpoint1}{'start'} > $breakpoint{$breakpoint2}{'end'}) {
          ($ref_start, $ref_end) = ($breakpoint{$breakpoint2}{'start'}, $breakpoint{$breakpoint1}{'end'});
        } else {
          ($ref_start, $ref_end) = ($breakpoint{$breakpoint1}{'start'}, $breakpoint{$breakpoint2}{'end'});
        }
        # check if there are pre-sampled inversions overlapping at the same site already
        foreach my $inversion_id (sort keys %inversion) {
          if ($inversion{$inversion_id}{'ref_chr'} eq $ref_chr) {
            my $check_overlap_flag1 = check_overlap_region($breakpoint{$breakpoint1}{'start'}, $breakpoint{$breakpoint1}{'end'}, $inversion{$inversion_id}{'ref_start'}, $inversion{$inversion_id}{'ref_end'});
            if ($check_overlap_flag1 == 1) {
              delete $valid_breakpoint_pair{'+-'}{$breakpoint1};
              foreach my $bp2 (sort keys %{$valid_breakpoint_pair{'-+'}}) {
                if (exists $valid_breakpoint_pair{'-+'}{$bp2}{$breakpoint1}) {
                  delete $valid_breakpoint_pair{'-+'}{$bp2}{$breakpoint1};
                }
                if (scalar (keys %{$valid_breakpoint_pair{'-+'}{$bp2}}) == 0) {
                  delete $valid_breakpoint_pair{'-+'}{$bp2};
                }
              }
            }
            my $check_overlap_flag2 = check_overlap_region($breakpoint{$breakpoint2}{'start'}, $breakpoint{$breakpoint2}{'end'}, $inversion{$inversion_id}{'ref_start'}, $inversion{$inversion_id}{'ref_end'});
            if ($check_overlap_flag2 == 1) {
              delete $valid_breakpoint_pair{'-+'}{$breakpoint2};
              foreach my $bp1 (sort keys %{$valid_breakpoint_pair{'+-'}}) {
                if (exists $valid_breakpoint_pair{'+-'}{$bp1}{$breakpoint2}) {
                  delete $valid_breakpoint_pair{'+-'}{$bp1}{$breakpoint2};
                }
                if (scalar (keys %{$valid_breakpoint_pair{'+-'}{$bp1}}) == 0) {
                  delete $valid_breakpoint_pair{'+-'}{$bp1};
                }
              }
            }
            if (($check_overlap_flag1 == 1) or ($check_overlap_flag2 == 1)) {
              goto SAMPLE_RANDOM_INV1;
            }
          }
        }

        my $inversion_id = "INV_${i}";
        $inversion{$inversion_id}{'ref_chr'} = $ref_chr;
        $inversion{$inversion_id}{'ref_start'} = $ref_start;
        $inversion{$inversion_id}{'ref_end'} = $ref_end;
        $inversion{$inversion_id}{'breakpoint1'} = $breakpoint1;
        $inversion{$inversion_id}{'breakpoint2'} = $breakpoint2;
        # print "sampled inversion: ref_chr=$ref_chr, ref_start=$ref_start, ref_end=$ref_end, breakpoint1=$breakpoint1, breakpoint2=$breakpoint2\n";

        # delete used breakpoints
        delete $valid_breakpoint_pair{'+-'}{$breakpoint1};
        foreach my $bp2 (sort keys %{$valid_breakpoint_pair{'-+'}}) {
          if (exists $valid_breakpoint_pair{'-+'}{$bp2}{$breakpoint1}) {
            delete $valid_breakpoint_pair{'-+'}{$bp2}{$breakpoint1};
          }
          if (scalar (keys %{$valid_breakpoint_pair{'-+'}{$bp2}}) == 0) {
            delete $valid_breakpoint_pair{'-+'}{$bp2};
          }
        }
        delete $valid_breakpoint_pair{'-+'}{$breakpoint2};
        foreach my $bp1 (sort keys %{$valid_breakpoint_pair{'+-'}}) {
          if (exists $valid_breakpoint_pair{'+-'}{$bp1}{$breakpoint2}) {
            delete $valid_breakpoint_pair{'+-'}{$bp1}{$breakpoint2};
          }
          if (scalar (keys %{$valid_breakpoint_pair{'+-'}{$bp1}}) == 0) {
            delete $valid_breakpoint_pair{'+-'}{$bp1};
          }
        }
      }
    }
  } else {
    for (my $i = 1; $i <= $inversion_count; $i++) {
      # random sampling across the genome
      my %refseq_genome_space = create_genome_space($refseq_hashref);
      # sample the inverted region from the genome space
      my $chr_end_margin = 1000;
      SAMPLE_RANDOM_INV2:
      # sample the inverted region (for copy number loss)
      my ($ref_chr, $ref_start) = sample_genome_space(\%refseq_genome_space);
      # print "ref_chr = $ref_chr, ref_start = $ref_start\n";
      # check if the chromosome-end has been involved in the simulated inversion
      if ($ref_start <= $chr_end_margin) {
	  goto SAMPLE_RANDOM_INV2;
      }
      # sample inversion size
      my $inversion_size = sample_inversion_size($inversion_min_size, $inversion_max_size);
      my $ref_end = $ref_start + $inversion_size - 1;
      # print "ref_chr=$ref_chr, ref_start = $ref_start, inversion_size = $inversion_size\n";
      # check if the sampled end position will go beyond the chromosome end
      if ($ref_end > ($refseq_genome_space{'chr-wide'}{$ref_chr}{'length'} - $chr_end_margin)) {
        goto SAMPLE_RANDOM_INV2;
      }
      # check if the sampled region overlaped with the defined centromeres
      if (defined $centromere_gff) {
        foreach my $centromere_id (sort keys %centromere) {
          if ($centromere{$centromere_id}{'chr'} eq $ref_chr) {
            my $centromere_check_flag = check_overlap_region($ref_start, $ref_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
            if ($centromere_check_flag == 1) {
              goto SAMPLE_RANDOM_INV2;
            }
          }
        }
      }

      # check if there are pre-sampled inversions overlapping at the same site already
      my $check_overlap_flag = 0;
      foreach my $id (sort keys %inversion) {
        if ($inversion{$id}{'ref_chr'} eq $ref_chr) {
          $check_overlap_flag = check_overlap_region($ref_start, $ref_end, $inversion{$id}{'ref_start'}, $inversion{$id}{'ref_end'});
          if ($check_overlap_flag == 1) {
            goto SAMPLE_RANDOM_INV2;
          }
        }
      }

      # all check passed, register this INV
      my $inversion_id = "INV_${i}";
      $inversion{$inversion_id}{'ref_chr'} = $ref_chr;
      $inversion{$inversion_id}{'ref_start'} = $ref_start;
      $inversion{$inversion_id}{'ref_end'} = $ref_end;
      $inversion{$inversion_id}{'breakpoint1'} = "$ref_chr:$ref_start";
      $inversion{$inversion_id}{'breakpoint2'} = "$ref_chr:$ref_end";
      # print "sampled inversion: ref_chr=$ref_chr, ref_start=$ref_start, ref_end=$ref_end\n";
    }
  }
  # print Dumper(%inversion);

  # introduce all sampled inversions
  foreach my $inversion_id (sort keys %inversion) {
    my $ref_chr = $inversion{$inversion_id}{'ref_chr'};
    my $ref_start = $inversion{$inversion_id}{'ref_start'};
    my $ref_end = $inversion{$inversion_id}{'ref_end'};
    my $inversion_size = $ref_end - $ref_start + 1;
    # print "inversion_id = $inversion_id, ref_chr = $ref_chr, ref_start = $ref_start\n\n";
    my $ref_allele = substr $$simseq_hashref{$ref_chr}, $ref_start - 1, 1;
    my $alt_allele = "<INV>";
    my $inversion_seq = substr $$simseq_hashref{$ref_chr}, $ref_start - 1, $inversion_size;
    $inversion_seq = revcom($inversion_seq);
    substr $$simseq_hashref{$ref_chr}, $ref_start - 1, $inversion_size, $inversion_seq;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'} = "+";
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'} = $ref_chr;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'} = $ref_start;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'} = $ref_end;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'} = $alt_allele;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'} = "-";
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} = "INV";
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'} = $inversion_id;
    $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'inversion_size'} = $inversion_size;
  }
}

sub extract_translocation_from_sv {
  my $sv_hashref = shift @_;
  my %translocation = ();
  foreach my $sv_event (sort keys %$sv_hashref) {
    my $check_overlap_flag = 0;
    if ($$sv_hashref{$sv_event}{'sv_type'} eq "to_be_classified") {
      my $bnd_count = scalar @{$$sv_hashref{$sv_event}{'BND'}};
      # print "sv_event = $sv_event, bnd_count = $bnd_count\n";
      if ($bnd_count == 4) {
        my %bnd = ();
        foreach my $b_hashref (sort @{$$sv_hashref{$sv_event}{'BND'}}) {
          my $ref_chr = $$b_hashref{'ref_chr'};
          my $ref_start = $$b_hashref{'ref_start'};
          my $s = $$b_hashref{'s'};
          my $p = $$b_hashref{'p'};
          my $t = $$b_hashref{'t'};
          my ($bnd_p_chr, $bnd_p_start) = split /:/, $p;
          # print "p=$p, bnd_p_chr=$bnd_p_chr, bnd_p_start=$bnd_p_start\n";
          my $p_relative_strand = $$b_hashref{'p_relative_strand'};
          my $p_relative_position = $$b_hashref{'p_relative_position'};
          $bnd{$ref_chr}{$ref_start}{'s'} = $s;
          $bnd{$ref_chr}{$ref_start}{'t'} = $t;
          $bnd{$ref_chr}{$ref_start}{'p_chr'} = $bnd_p_chr;
          $bnd{$ref_chr}{$ref_start}{'p_start'} = $bnd_p_start;
          $bnd{$ref_chr}{$ref_start}{'p_relative_strand'} = $p_relative_strand;
          $bnd{$ref_chr}{$ref_start}{'p_relative_position'} = $p_relative_position;
        }
        # verify this is indeed a translocation
        my $translocation_test_flag = 0;
        my $translocation_resolution = "NA";
        # chromosome check
        my @bnd_ref_chr = sort keys %bnd;
        if ((scalar @bnd_ref_chr) != 2) {
          $translocation_test_flag = 1;
        } else {
        # breakpoint check
          my $bnd_ref_chr1 = $bnd_ref_chr[0];
          my $bnd_ref_chr2 = $bnd_ref_chr[1];
          foreach my $ref_chr (@bnd_ref_chr) {
            my @bnd_position = sort {$a <=> $b} keys %{$bnd{$ref_chr}};
            if ((scalar @bnd_position) != 2) {
              print "\n!!! Error! There seems mistakes in the defined translocation(s) invovled with $ref_chr.\n";
              print "!!! Please check the input translocation_vcf file\n";
              print "!!! Exit!\n";
              $translocation_test_flag = 1;
              die;
            } elsif (($bnd_position[1] - $bnd_position[0]) != 1) {
              print "\n!!! Error! There seems mistakes in the defined translocation(s) involved with $ref_chr.\n";
              print "!!! Please check the input translocation_vcf file\n";
              print "!!! Exit!\n";
              $translocation_test_flag = 1;
              die;
            }
          }
          # breakpoint and strand check
          if ($translocation_test_flag == 0) {
            my ($bnd_ref_chr1_bnd_position1, $bnd_ref_chr1_bnd_position2) = sort {$a <=> $b} keys %{$bnd{$bnd_ref_chr1}};
            my ($bnd_ref_chr2_bnd_position1, $bnd_ref_chr2_bnd_position2) = sort {$a <=> $b} keys %{$bnd{$bnd_ref_chr2}};

            my $bnd_ref_chr1_bnd_position1_p_chr = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position1}{'p_chr'};
            my $bnd_ref_chr1_bnd_position1_p_start = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position1}{'p_start'};
            my $bnd_ref_chr1_bnd_position1_p_relative_strand = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position1}{'p_relative_strand'};
            my $bnd_ref_chr1_bnd_position1_p_relative_position = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position1}{'p_relative_position'};

            my $bnd_ref_chr1_bnd_position2_p_chr = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position2}{'p_chr'};
            my $bnd_ref_chr1_bnd_position2_p_start = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position2}{'p_start'};
            my $bnd_ref_chr1_bnd_position2_p_relative_strand = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position2}{'p_relative_strand'};
            my $bnd_ref_chr1_bnd_position2_p_relative_position = $bnd{$bnd_ref_chr1}{$bnd_ref_chr1_bnd_position2}{'p_relative_position'};

            my $bnd_ref_chr2_bnd_position1_p_chr = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position1}{'p_chr'};
            my $bnd_ref_chr2_bnd_position1_p_start = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position1}{'p_start'};
            my $bnd_ref_chr2_bnd_position1_p_relative_strand = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position1}{'p_relative_strand'};
            my $bnd_ref_chr2_bnd_position1_p_relative_position = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position1}{'p_relative_position'};

            my $bnd_ref_chr2_bnd_position2_p_chr = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position2}{'p_chr'};
            my $bnd_ref_chr2_bnd_position2_p_start = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position2}{'p_start'};
            my $bnd_ref_chr2_bnd_position2_p_relative_strand = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position2}{'p_relative_strand'};
            my $bnd_ref_chr2_bnd_position2_p_relative_position = $bnd{$bnd_ref_chr2}{$bnd_ref_chr2_bnd_position2}{'p_relative_position'};

            # print "bnd_ref_chr1 = $bnd_ref_chr1\n";
            # print "bnd_ref_chr2 = $bnd_ref_chr2\n";
            # print "bnd_ref_chr1_bnd_position1_p_chr = $bnd_ref_chr1_bnd_position1_p_chr\n";
            # print "bnd_ref_chr1_bnd_position2_p_chr = $bnd_ref_chr1_bnd_position2_p_chr\n";
            # print "bnd_ref_chr2_bnd_position1_p_chr = $bnd_ref_chr2_bnd_position1_p_chr\n";
            # print "bnd_ref_chr2_bnd_position2_p_chr = $bnd_ref_chr2_bnd_position2_p_chr\n";
            # print "bnd_ref_chr1_bnd_position1_p_relative_strand = $bnd_ref_chr1_bnd_position1_p_relative_strand\n";
            # print "bnd_ref_chr1_bnd_position2_p_relative_strand = $bnd_ref_chr1_bnd_position2_p_relative_strand\n";
            # print "bnd_ref_chr2_bnd_position1_p_relative_strand = $bnd_ref_chr2_bnd_position1_p_relative_strand\n";
            # print "bnd_ref_chr2_bnd_position2_p_relative_strand = $bnd_ref_chr2_bnd_position2_p_relative_strand\n";
            # print "bnd_ref_chr1_bnd_position1_p_relative_position = $bnd_ref_chr1_bnd_position1_p_relative_position\n";
            # print "bnd_ref_chr1_bnd_position2_p_relative_position = $bnd_ref_chr1_bnd_position2_p_relative_position\n";
            # print "bnd_ref_chr2_bnd_position1_p_relative_position = $bnd_ref_chr2_bnd_position1_p_relative_position\n";
            # print "bnd_ref_chr2_bnd_position2_p_relative_position = $bnd_ref_chr2_bnd_position2_p_relative_position\n";

            if (($bnd_ref_chr1 ne $bnd_ref_chr2_bnd_position1_p_chr) or ($bnd_ref_chr1 ne $bnd_ref_chr2_bnd_position2_p_chr) or ($bnd_ref_chr2 ne $bnd_ref_chr1_bnd_position1_p_chr) or ($bnd_ref_chr2 ne $bnd_ref_chr1_bnd_position2_p_chr)) {
              print "\n!!! Error! There seems mistakes in the defined translocation(s) involved with $bnd_ref_chr1 and $bnd_ref_chr2.\n";
              print "!!! Please check the input translocation_vcf file\n";
              print "!!! Exit!\n";
              $translocation_test_flag = 1;
              die;
            } else {
              my $positive_strand_count = 0;
              if ($bnd_ref_chr1_bnd_position1_p_relative_strand eq "+") {
                $positive_strand_count++;
              }
              if ($bnd_ref_chr1_bnd_position2_p_relative_strand eq "+") {
                $positive_strand_count++;
              }
              if ($bnd_ref_chr2_bnd_position1_p_relative_strand eq "+") {
                $positive_strand_count++;
              }
              if ($bnd_ref_chr2_bnd_position2_p_relative_strand eq "+") {
                $positive_strand_count++;
              }
              if ($positive_strand_count == 4) {
                $translocation_resolution = "++++"; # strand choice for chr1_chr2_part1, chr1_chr2_part2, chr2_chr1_part1, chr2_chr1_part2;
              } elsif ($positive_strand_count == 0) {
                $translocation_resolution = "+--+";
              } else {
                print "\n!!! Error! There seems mistakes in the defined translocation(s) involved with $bnd_ref_chr1 and $bnd_ref_chr2.\n";
                print "!!! Please check the input translocation_vcf file\n";
                print "!!! Exit!\n";
                $translocation_test_flag = 1;
                die;
              }
            }
          }
        }

        if ($translocation_test_flag == 0) {
          # check if overlapped with pre-registed translocations
          my ($bnd_ref_chr1, $bnd_ref_chr2) = sort @bnd_ref_chr;
          my ($bnd_ref_chr1_bnd_position1, $bnd_ref_chr1_bnd_position2) = sort {$a <=> $b} keys %{$bnd{$bnd_ref_chr1}};
          my ($bnd_ref_chr2_bnd_position1, $bnd_ref_chr2_bnd_position2) = sort {$a <=> $b} keys %{$bnd{$bnd_ref_chr2}};

          my $check_overlap_flag = check_translocation_overlap(\%translocation, $bnd_ref_chr1, $bnd_ref_chr2);
          if ($check_overlap_flag == 1) {
            print "\n!!! Warning! Multiple translocation defined for $bnd_ref_chr1 or $bnd_ref_chr2!\n";
            print "!!! Only keep the first one and ignore the others\n";
          } else {
            # register this translocation
            my $ref_chr1 = $bnd_ref_chr1;
            my $ref_chr2 = $bnd_ref_chr2;

            my $ref_chr1_chr2_part1_chr;
            my $ref_chr1_chr2_part1_start;
            my $ref_chr1_chr2_part1_end;
            my $ref_chr1_chr2_part1_strand;

            my $ref_chr1_chr2_part2_chr;
            my $ref_chr1_chr2_part2_start;
            my $ref_chr1_chr2_part2_end;
            my $ref_chr1_chr2_part2_strand;

            my $ref_chr2_chr1_part1_chr;
            my $ref_chr2_chr1_part1_start;
            my $ref_chr2_chr1_part1_end;
            my $ref_chr2_chr1_part1_strand;

            my $ref_chr2_chr1_part2_chr;
            my $ref_chr2_chr1_part2_start;
            my $ref_chr2_chr1_part2_end;
            my $ref_chr2_chr1_part2_strand;

            if ($translocation_resolution eq "++++") {
              $ref_chr1_chr2_part1_chr = $ref_chr1;
              $ref_chr1_chr2_part1_start = 1;
              $ref_chr1_chr2_part1_end = $bnd_ref_chr1_bnd_position1;
              $ref_chr1_chr2_part1_strand = "+";

              $ref_chr1_chr2_part2_chr = $ref_chr2;
              $ref_chr1_chr2_part2_start = $bnd_ref_chr2_bnd_position2;
              $ref_chr1_chr2_part2_end = "chr_end";
              $ref_chr1_chr2_part2_strand = "+";

              $ref_chr2_chr1_part1_chr = $ref_chr2;
              $ref_chr2_chr1_part1_start = 1;
              $ref_chr2_chr1_part1_end = $bnd_ref_chr2_bnd_position1;
              $ref_chr2_chr1_part1_strand = "+";

              $ref_chr2_chr1_part2_chr = $ref_chr1;
              $ref_chr2_chr1_part2_start = $bnd_ref_chr1_bnd_position2;
              $ref_chr2_chr1_part2_end = "chr_end";
              $ref_chr2_chr1_part2_strand = "+";
            } else {
              $ref_chr1_chr2_part1_chr = $ref_chr1;
              $ref_chr1_chr2_part1_start = 1;
              $ref_chr1_chr2_part1_end = $bnd_ref_chr1_bnd_position1;
              $ref_chr1_chr2_part1_strand = "+";

              $ref_chr1_chr2_part2_chr = $ref_chr2;
              $ref_chr1_chr2_part2_start = 1;
              $ref_chr1_chr2_part2_end = $bnd_ref_chr2_bnd_position1;
              $ref_chr1_chr2_part2_strand = "-";

              $ref_chr2_chr1_part1_chr = $ref_chr2;
              $ref_chr2_chr1_part1_start = $bnd_ref_chr2_bnd_position2;
              $ref_chr2_chr1_part1_end = "chr_end";
              $ref_chr2_chr1_part1_strand = "-";

              $ref_chr2_chr1_part2_chr = $ref_chr1;
              $ref_chr2_chr1_part2_start = $bnd_ref_chr1_bnd_position2;
              $ref_chr2_chr1_part2_end = "chr_end";
              $ref_chr2_chr1_part2_strand = "+";
            }

            $translocation{$sv_event}{'ref_chr1'} = $ref_chr1;
            $translocation{$sv_event}{'ref_chr2'} = $ref_chr2;

            $translocation{$sv_event}{'ref_chr1_chr2_part1_chr'} = $ref_chr1_chr2_part1_chr;
            $translocation{$sv_event}{'ref_chr1_chr2_part1_start'} = $ref_chr1_chr2_part1_start;
            $translocation{$sv_event}{'ref_chr1_chr2_part1_end'} = $ref_chr1_chr2_part1_end;
            $translocation{$sv_event}{'ref_chr1_chr2_part1_strand'} = $ref_chr1_chr2_part1_strand;

            $translocation{$sv_event}{'ref_chr1_chr2_part2_chr'} = $ref_chr1_chr2_part2_chr;
            $translocation{$sv_event}{'ref_chr1_chr2_part2_start'} = $ref_chr1_chr2_part2_start;
            $translocation{$sv_event}{'ref_chr1_chr2_part2_end'} = $ref_chr1_chr2_part2_end;
            $translocation{$sv_event}{'ref_chr1_chr2_part2_strand'} = $ref_chr1_chr2_part2_strand;

            $translocation{$sv_event}{'ref_chr2_chr1_part1_chr'} = $ref_chr2_chr1_part1_chr;
            $translocation{$sv_event}{'ref_chr2_chr1_part1_start'} = $ref_chr2_chr1_part1_start;
            $translocation{$sv_event}{'ref_chr2_chr1_part1_end'} = $ref_chr2_chr1_part1_end;
            $translocation{$sv_event}{'ref_chr2_chr1_part1_strand'} = $ref_chr2_chr1_part1_strand;

            $translocation{$sv_event}{'ref_chr2_chr1_part2_chr'} = $ref_chr2_chr1_part2_chr;
            $translocation{$sv_event}{'ref_chr2_chr1_part2_start'} = $ref_chr2_chr1_part2_start;
            $translocation{$sv_event}{'ref_chr2_chr1_part2_end'} = $ref_chr2_chr1_part2_end;
            $translocation{$sv_event}{'ref_chr2_chr1_part2_strand'} = $ref_chr2_chr1_part2_strand;
          }
        }
      }
    }
  }
  return %translocation;
}

sub check_translocation_overlap {
  my ($translocation_hashref, $chr1, $chr2) = @_;
  my $flag = 0;
  foreach my $translocation_id (sort keys %$translocation_hashref) {
    if ((exists $$translocation_hashref{$translocation_id}{$chr1}) or (exists $$translocation_hashref{$translocation_id}{$chr2})) {
      $flag = 1;
      last;
    }
  }
  return $flag;
}


sub introduce_defined_translocation {
  my ($translocation_hashref, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  foreach my $translocation_id (sort keys %$translocation_hashref) {
    my $ref_chr1 = $$translocation_hashref{$translocation_id}{'ref_chr1'};
    my $ref_chr2 = $$translocation_hashref{$translocation_id}{'ref_chr2'};

    my $ref_chr1_chr2_part1_chr = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part1_chr'};
    my $ref_chr1_chr2_part1_start = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part1_start'};
    my $ref_chr1_chr2_part1_end = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part1_end'};
    my $ref_chr1_chr2_part1_strand = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part1_strand'};

    my $ref_chr1_chr2_part2_chr = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part2_chr'};
    my $ref_chr1_chr2_part2_start = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part2_start'};
    my $ref_chr1_chr2_part2_end = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part2_end'};
    my $ref_chr1_chr2_part2_strand = $$translocation_hashref{$translocation_id}{'ref_chr1_chr2_part2_strand'};

    my $ref_chr2_chr1_part1_chr = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part1_chr'};
    my $ref_chr2_chr1_part1_start = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part1_start'};
    my $ref_chr2_chr1_part1_end = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part1_end'};
    my $ref_chr2_chr1_part1_strand = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part1_strand'};

    my $ref_chr2_chr1_part2_chr = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part2_chr'};
    my $ref_chr2_chr1_part2_start = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part2_start'};
    my $ref_chr2_chr1_part2_end = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part2_end'};
    my $ref_chr2_chr1_part2_strand = $$translocation_hashref{$translocation_id}{'ref_chr2_chr1_part2_strand'};

    my $ref_chr1_chr2 = "${ref_chr1}_${ref_chr2}";
    my $ref_chr2_chr1 = "${ref_chr2}_${ref_chr1}";

    if ($ref_chr1_chr2_part1_strand eq $ref_chr1_chr2_part2_strand) {
      $ref_chr1_chr2_part2_end = length $$refseq_hashref{$ref_chr2};
      $ref_chr2_chr1_part2_end = length $$refseq_hashref{$ref_chr1};
    } else {
      $ref_chr2_chr1_part1_end = length $$refseq_hashref{$ref_chr2};
      $ref_chr2_chr1_part2_end = length $$refseq_hashref{$ref_chr1};
    }

    my $ref_chr1_chr2_part1_seq = substr $$refseq_hashref{$ref_chr1}, $ref_chr1_chr2_part1_start - 1, $ref_chr1_chr2_part1_end - $ref_chr1_chr2_part1_start + 1;
    my $ref_chr1_chr2_part2_seq = substr $$refseq_hashref{$ref_chr2}, $ref_chr1_chr2_part2_start - 1, $ref_chr1_chr2_part2_end - $ref_chr1_chr2_part2_start + 1;
    my $ref_chr2_chr1_part1_seq = substr $$refseq_hashref{$ref_chr2}, $ref_chr2_chr1_part1_start - 1, $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1;
    my $ref_chr2_chr1_part2_seq = substr $$refseq_hashref{$ref_chr1}, $ref_chr2_chr1_part2_start - 1, $ref_chr2_chr1_part2_end - $ref_chr2_chr1_part2_start + 1;

    if ($ref_chr1_chr2_part1_strand ne $ref_chr1_chr2_part2_strand) {
      $ref_chr1_chr2_part2_seq = revcom($ref_chr1_chr2_part2_seq);
      $ref_chr2_chr1_part1_seq = revcom($ref_chr2_chr1_part1_seq);
    }

    $$simseq_hashref{$ref_chr1_chr2} = $ref_chr1_chr2_part1_seq . $ref_chr1_chr2_part2_seq;
    $$simseq_hashref{$ref_chr2_chr1} = $ref_chr2_chr1_part1_seq . $ref_chr2_chr1_part2_seq;
    delete $$simseq_hashref{$ref_chr1};
    delete $$simseq_hashref{$ref_chr2};

    if ($ref_chr1_chr2_part1_strand eq $ref_chr1_chr2_part2_strand) {
      # translocation resolution: ++++
      # ref_chr1_chr2_part1
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_start'} = $ref_chr1_chr2_part1_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_strand'} = $ref_chr1_chr2_part1_strand;
      my $ref_chr1_chr2_part1_end_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr1_chr2_part1_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_allele'} = $ref_chr1_chr2_part1_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_strand'} = $ref_chr1_chr2_part1_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_allele'} = "$ref_chr1_chr2_part1_end_ref_allele\[$ref_chr2:$ref_chr1_chr2_part2_start\[";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_id'} = $translocation_id;
      # ref_chr1_chr2_part2
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_start'} = $ref_chr1_chr2_part2_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_end'} = $ref_chr1_chr2_part2_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_strand'} = $ref_chr1_chr2_part2_strand;
      my $ref_chr1_chr2_part2_start_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr1_chr2_part2_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_allele'} = $ref_chr1_chr2_part2_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_start'} = $ref_chr1_chr2_part1_end + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_end'} = $ref_chr1_chr2_part1_end + $ref_chr1_chr2_part2_end - $ref_chr1_chr2_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_strand'} = $ref_chr1_chr2_part2_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_allele'} = "\]$ref_chr1:$ref_chr1_chr2_part1_end\]$ref_chr1_chr2_part2_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part1
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_start'} = $ref_chr2_chr1_part1_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_end'} = $ref_chr2_chr1_part1_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_strand'} = $ref_chr2_chr1_part1_strand;
      my $ref_chr2_chr1_part1_end_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr2_chr1_part1_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_allele'} = $ref_chr2_chr1_part1_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_end'} = $ref_chr2_chr1_part1_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_strand'} = $ref_chr2_chr1_part1_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_allele'} = "$ref_chr2_chr1_part1_end_ref_allele\[$ref_chr1:$ref_chr2_chr1_part2_start\[";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part2
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_start'} = $ref_chr2_chr1_part2_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_end'} = $ref_chr2_chr1_part2_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_strand'} = $ref_chr2_chr1_part2_strand;
      my $ref_chr2_chr1_part2_start_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr2_chr1_part2_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_allele'} = $ref_chr2_chr1_part2_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_start'} = $ref_chr2_chr1_part1_end + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_end'} = $ref_chr2_chr1_part1_end + $ref_chr2_chr1_part2_end - $ref_chr2_chr1_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_strand'} = $ref_chr2_chr1_part2_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_allele'} = "\]$ref_chr2:$ref_chr2_chr1_part1_end\]$ref_chr2_chr1_part2_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_id'} = $translocation_id;
    } else {
      # translocation resolution: +--+
      # ref_chr1_chr2_part1
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_start'} = $ref_chr1_chr2_part1_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_strand'} = $ref_chr1_chr2_part1_strand;
      my $ref_chr1_chr2_part1_end_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr1_chr2_part1_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_allele'} = $ref_chr1_chr2_part1_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_strand'} = $ref_chr1_chr2_part1_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_allele'} = "$ref_chr1_chr2_part1_end_ref_allele\]$ref_chr2:$ref_chr1_chr2_part2_end\]";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_id'} = $translocation_id;
      # ref_chr1_chr2_part2
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_start'} = $ref_chr1_chr2_part2_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_end'} = $ref_chr1_chr2_part2_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_strand'} = $ref_chr1_chr2_part2_strand;
      my $ref_chr1_chr2_part2_end_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr1_chr2_part2_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_allele'} = $ref_chr1_chr2_part2_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_start'} = $ref_chr1_chr2_part1_end + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_end'} = $ref_chr1_chr2_part1_end + $ref_chr1_chr2_part2_end - $ref_chr1_chr2_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_strand'} = $ref_chr1_chr2_part2_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_allele'} = "$ref_chr1_chr2_part2_end_ref_allele\]$ref_chr1:$ref_chr1_chr2_part1_end\]";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part1
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_start'} = $ref_chr2_chr1_part1_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_end'} = $ref_chr2_chr1_part1_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_strand'} = $ref_chr2_chr1_part1_strand;
      my $ref_chr2_chr1_part1_start_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr2_chr1_part1_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_allele'} = $ref_chr2_chr1_part1_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_end'} = $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_strand'} = $ref_chr2_chr1_part1_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_allele'} = "\[$ref_chr1:$ref_chr2_chr1_part2_start\[$ref_chr2_chr1_part1_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part2
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_start'} = $ref_chr2_chr1_part2_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_end'} = $ref_chr2_chr1_part2_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_strand'} = $ref_chr2_chr1_part2_strand;
      my $ref_chr2_chr1_part2_start_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr2_chr1_part2_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_allele'} = $ref_chr2_chr1_part2_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_start'} = $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1 + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_end'} = $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1 + $ref_chr2_chr1_part2_end - $ref_chr2_chr1_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_strand'} = $ref_chr2_chr1_part2_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_allele'} = "\[$ref_chr2:$ref_chr2_chr1_part1_start\[$ref_chr2_chr1_part2_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_id'} = $translocation_id;
    }
  }
}

sub examine_breakpoint_for_translocation {
  my ($breakpoint_hashref, $refseq_hashref, $centromere_gff) = @_;
  my @breakpoint = sort keys %$breakpoint_hashref;
  my %valid_breakpoint_pair = ();
  foreach my $b1 (@breakpoint) {
    foreach my $b2 (@breakpoint) {
      if ($b1 ne $b2) {
        my $b1_feature = $$breakpoint_hashref{$b1}{'feature'};
        my $b1_chr = $$breakpoint_hashref{$b1}{'chr'};
        my $b1_start = $$breakpoint_hashref{$b1}{'start'};
        my $b1_end = $$breakpoint_hashref{$b1}{'end'};
        my $b1_strand = $$breakpoint_hashref{$b1}{'strand'};

        my $b2_feature = $$breakpoint_hashref{$b2}{'feature'};
        my $b2_chr = $$breakpoint_hashref{$b2}{'chr'};
        my $b2_start = $$breakpoint_hashref{$b2}{'start'};
        my $b2_end = $$breakpoint_hashref{$b2}{'end'};
        my $b2_strand = $$breakpoint_hashref{$b2}{'strand'};


        if (($b1_feature eq $b2_feature) and ($b1_chr ne $b2_chr)) {
          if (defined $centromere_gff) {
            # check for the case in which the translocated chromosomes will have zero or two centromeres
            my $centromere_gff_fh = read_file($centromere_gff);
            my %centromere = parse_gff_file($centromere_gff_fh);
            close $centromere_gff_fh;

            if ((($b1_strand eq "+") and ($b2_strand eq "+")) or (($b1_strand eq "-") and ($b2_strand eq "-"))) {
              # translocation resolution: ++++
              my $ref_chr1_chr2_part1_chr = $b1_chr;
              my $ref_chr1_chr2_part1_start = 1;
              my $ref_chr1_chr2_part1_end = $b1_start - 1;
              my $ref_chr1_chr2_part1_strand = "+";

              my $ref_chr1_chr2_part2_chr = $b2_chr;
              my $ref_chr1_chr2_part2_start = $b2_start;
              my $ref_chr1_chr2_part2_end = length $$refseq_hashref{$b2_chr};
              my $ref_chr1_chr2_part2_strand = "+";

              my $ref_chr2_chr1_part1_chr = $b2_chr;
              my $ref_chr2_chr1_part1_start = 1;
              my $ref_chr2_chr1_part1_end = $b2_start - 1;
              my $ref_chr2_chr1_part1_strand = "+";

              my $ref_chr2_chr1_part2_chr = $b1_chr;
              my $ref_chr2_chr1_part2_start = $b1_start;
              my $ref_chr2_chr1_part2_end = length $$refseq_hashref{$b1_chr};
              my $ref_chr2_chr1_part2_strand = "+";

              my $ref_chr1_chr2_part1_centromere_count = 0;
              my $ref_chr1_chr2_part2_centromere_count = 0;
              my $ref_chr2_chr1_part1_centromere_count = 0;
              my $ref_chr2_chr1_part2_centromere_count = 0;

              foreach my $centromere_id (sort keys %centromere) {
                if ($centromere{$centromere_id}{'chr'} eq $b1_chr) {
                  $ref_chr1_chr2_part1_centromere_count = check_overlap_region($ref_chr1_chr2_part1_start, $ref_chr1_chr2_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                  $ref_chr2_chr1_part2_centromere_count = check_overlap_region($ref_chr2_chr1_part2_start, $ref_chr2_chr1_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                }
                if ($centromere{$centromere_id}{'chr'} eq $b2_chr) {
                  $ref_chr1_chr2_part2_centromere_count = check_overlap_region($ref_chr1_chr2_part2_start, $ref_chr1_chr2_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                  $ref_chr2_chr1_part1_centromere_count = check_overlap_region($ref_chr2_chr1_part1_start, $ref_chr2_chr1_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                }
              }
              my $ref_chr1_chr2_centromere_count = $ref_chr1_chr2_part1_centromere_count + $ref_chr1_chr2_part2_centromere_count;
              my $ref_chr2_chr1_centromere_count = $ref_chr2_chr1_part1_centromere_count + $ref_chr2_chr1_part2_centromere_count;

              if (($ref_chr1_chr2_centromere_count == 1) and ($ref_chr2_chr1_centromere_count == 1)) {
                $valid_breakpoint_pair{$b1}{$b2}{'translocation_resolution'} = "++++";
                $valid_breakpoint_pair{$b2}{$b1}{'translocation_resolution'} = "++++";
              }
            } elsif ((($b1_strand eq "+") and ($b2_strand eq "-")) or (($b1_strand eq "-") and ($b2_strand eq "+"))) {
              # translocation resolution: +--+
              my $ref_chr1_chr2_part1_chr = $b1_chr;
              my $ref_chr1_chr2_part1_start = 1;
              my $ref_chr1_chr2_part1_end = $b1_start - 1;
              my $ref_chr1_chr2_part1_strand = "+";

              my $ref_chr1_chr2_part2_chr = $b2_chr;
              my $ref_chr1_chr2_part2_start = 1;
              my $ref_chr1_chr2_part2_end = $b2_end;
              my $ref_chr1_chr2_part2_strand = "-";

              my $ref_chr2_chr1_part1_chr = $b2_chr;
              my $ref_chr2_chr1_part1_start = $b2_end + 1;
              my $ref_chr2_chr1_part1_end = length $$refseq_hashref{$b2_chr};
              my $ref_chr2_chr1_part1_strand = "-";

              my $ref_chr2_chr1_part2_chr = $b1_chr;
              my $ref_chr2_chr1_part2_start = $b1_start;
              my $ref_chr2_chr1_part2_end = length $$refseq_hashref{$b1_chr};
              my $ref_chr2_chr1_part2_strand = "+";

              my $ref_chr1_chr2_part1_centromere_count = 0;
              my $ref_chr1_chr2_part2_centromere_count = 0;
              my $ref_chr2_chr1_part1_centromere_count = 0;
              my $ref_chr2_chr1_part2_centromere_count = 0;

              foreach my $centromere_id (sort keys %centromere) {
                if ($centromere{$centromere_id}{'chr'} eq $b1_chr) {
                  $ref_chr1_chr2_part1_centromere_count = check_overlap_region($ref_chr1_chr2_part1_start, $ref_chr1_chr2_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                  $ref_chr2_chr1_part2_centromere_count = check_overlap_region($ref_chr2_chr1_part2_start, $ref_chr2_chr1_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                }
                if ($centromere{$centromere_id}{'chr'} eq $b2_chr) {
                  $ref_chr1_chr2_part2_centromere_count = check_overlap_region($ref_chr1_chr2_part2_start, $ref_chr1_chr2_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                  $ref_chr2_chr1_part1_centromere_count = check_overlap_region($ref_chr2_chr1_part1_start, $ref_chr2_chr1_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                }
              }
              my $ref_chr1_chr2_centromere_count = $ref_chr1_chr2_part1_centromere_count + $ref_chr1_chr2_part2_centromere_count;
              my $ref_chr2_chr1_centromere_count = $ref_chr2_chr1_part1_centromere_count + $ref_chr2_chr1_part2_centromere_count;

              if (($ref_chr1_chr2_centromere_count == 1) and ($ref_chr2_chr1_centromere_count == 1)) {
                $valid_breakpoint_pair{$b1}{$b2}{'translocation_resolution'} = "+--+";
                $valid_breakpoint_pair{$b2}{$b1}{'translocation_resolution'} = "+--+";
              }
            } else {
              # when the strand of the corresponding breakpoint feature has not been specified (e.g. strand = ".")
              $valid_breakpoint_pair{$b1}{$b2}{'translocation_resolution'} = "NA";
              $valid_breakpoint_pair{$b2}{$b1}{'translocation_resolution'} = "NA";
            }
          } else {
            $valid_breakpoint_pair{$b1}{$b2}{'translocation_resolution'} = "NA";
            $valid_breakpoint_pair{$b2}{$b1}{'translocation_resolution'} = "NA";
          }
        }
      }
    }
  }
  if ((scalar (keys %valid_breakpoint_pair)) == 0) {
    print "\n!!! Error! None of defined breakpoints are valid for triggering translocation.\n";
    print "!!! Valid breakpoints need to satisfy the following criteria:\n";
    print "!!! 1) Be annotated as the same feature.\n";
    print "!!! 2) Come from at least two different chromosomes.\n";
    print "!!! 3) If centromeres have been defined, the breakpoints should not overlap with the defined centromeres.\n";
    print "!!! 4) If centromeres have been defined, the resulting rearranged chromosomes should still possess only one centromere per chromosome.\n";
    print "!!! Exit!\n";
    die;
  }
  return %valid_breakpoint_pair;
}

sub introduce_random_translocation {
  my ($translocation_count, $centromere_gff, $translocation_breakpoint_gff, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  my %centromere = ();
  if (defined $centromere_gff) {
    my $centromere_gff_fh = read_file($centromere_gff);
    %centromere = parse_gff_file($centromere_gff_fh);
  }

  my %translocation = ();
  my %available_chr = %$refseq_hashref;

  if (defined $translocation_breakpoint_gff) {
    my $translocation_breakpoint_gff_fh = read_file($translocation_breakpoint_gff);
    my %breakpoint = parse_gff_file($translocation_breakpoint_gff_fh);
    my %valid_breakpoint_pair = examine_breakpoint_for_translocation(\%breakpoint, $refseq_hashref, $centromere_gff);
    for (my $i = 1; $i <= $translocation_count; $i++) {
      # print "sample translocation: i = $i, translocation_count = $translocation_count\n";
      # check if there are still available chromsomes
      if ((scalar (keys %available_chr)) < 2) {
        my $j = $i - 1;
        print "\n!!! Warning! No more available chromosomes in the current simulation.\n";
        print "!!! Only $j translocation were introduced.\n";
        last;
      }
      # check if there are still valid $breakpoints
      if ((scalar (keys %valid_breakpoint_pair)) < 2) {
        my $j = $i - 1;
        print "\n!!! Warning! No more valid breakpoint pairs can be found in the current simulation based on the defined breakpoint file: $translocation_breakpoint_gff\n";
        print "!!! Only $j translocation were introduced.\n";
        last;
      }

      SAMPLE_RANDOM_TRA1:
      # sample inversion based on defined breakpoints
      my @breakpoint1 = shuffle(sort keys %valid_breakpoint_pair);
      my $breakpoint1 = shift @breakpoint1;
      my @breakpoint2 = shuffle(sort keys %{$valid_breakpoint_pair{$breakpoint1}});
      my $breakpoint2 = shift @breakpoint2;
      # print "breakpoint1 = $breakpoint1\n";
      # print "breakpoint2 = $breakpoint2\n";
      my $breakpoint1_chr = $breakpoint{$breakpoint1}{'chr'};
      my $breakpoint2_chr = $breakpoint{$breakpoint2}{'chr'};
      # check if the same chromosomes have been involved in previously simulated translocations
      if ((not exists $available_chr{$breakpoint1_chr}) or (not exists $available_chr{$breakpoint2_chr})) {
        # print "failed available chr check\n";
        delete $valid_breakpoint_pair{$breakpoint1}{$breakpoint2};
        delete $valid_breakpoint_pair{$breakpoint2}{$breakpoint1};
        goto SAMPLE_RANDOM_TRA1;
      }

      if ($valid_breakpoint_pair{$breakpoint1}{$breakpoint2}{'translocation_resolution'} eq "NA") {
        # random sampling of translocation resolution if it is unspecified
        if (rand(1) < 0.5) {
            $valid_breakpoint_pair{$breakpoint1}{$breakpoint2}{'translocation_resolution'} = "++++";
            $valid_breakpoint_pair{$breakpoint2}{$breakpoint1}{'translocation_resolution'} = "++++";
        } else {
            $valid_breakpoint_pair{$breakpoint1}{$breakpoint2}{'translocation_resolution'} = "+--+";
            $valid_breakpoint_pair{$breakpoint2}{$breakpoint1}{'translocation_resolution'} = "+--+";
        }
      }

      my $translocation_id = "TRA_${i}";
      my $ref_chr1;
      my $ref_chr2;
      my $ref_chr1_breakpoint_start;
      my $ref_chr1_breakpoint_end;
      my $ref_chr2_breakpoint_start;
      my $ref_chr2_breakpoint_end;

      if ($breakpoint1_chr lt $breakpoint2_chr) {
        $ref_chr1 = $breakpoint1_chr;
        $ref_chr2 = $breakpoint2_chr;
        $ref_chr1_breakpoint_start = $breakpoint{$breakpoint1}{'start'};
        $ref_chr1_breakpoint_end = $breakpoint{$breakpoint1}{'end'};
        $ref_chr2_breakpoint_start = $breakpoint{$breakpoint2}{'start'};
        $ref_chr2_breakpoint_end = $breakpoint{$breakpoint2}{'end'};
      } else {
        $ref_chr1 = $breakpoint2_chr;
        $ref_chr2 = $breakpoint1_chr;
        $ref_chr1_breakpoint_start = $breakpoint{$breakpoint2}{'start'};
        $ref_chr1_breakpoint_end = $breakpoint{$breakpoint2}{'end'};
        $ref_chr2_breakpoint_start = $breakpoint{$breakpoint1}{'start'};
        $ref_chr2_breakpoint_end = $breakpoint{$breakpoint1}{'end'};
      }

      if ($valid_breakpoint_pair{$breakpoint1}{$breakpoint2}{'translocation_resolution'} eq "++++") {
        # translocation_resolution: ++++
        $translocation{$translocation_id}{'ref_chr1'} = $ref_chr1;
        $translocation{$translocation_id}{'ref_chr2'} = $ref_chr2;

        $translocation{$translocation_id}{'ref_chr1_chr2_part1_chr'} = $ref_chr1;
        $translocation{$translocation_id}{'ref_chr1_chr2_part1_start'} = 1;
        $translocation{$translocation_id}{'ref_chr1_chr2_part1_end'} = $ref_chr1_breakpoint_start - 1;
        $translocation{$translocation_id}{'ref_chr1_chr2_part1_strand'} = "+";

        $translocation{$translocation_id}{'ref_chr1_chr2_part2_chr'} = $ref_chr2;
        $translocation{$translocation_id}{'ref_chr1_chr2_part2_start'} = $ref_chr2_breakpoint_start;
        $translocation{$translocation_id}{'ref_chr1_chr2_part2_end'} = length $$refseq_hashref{$ref_chr2};
        $translocation{$translocation_id}{'ref_chr1_chr2_part2_strand'} = "+";

        $translocation{$translocation_id}{'ref_chr2_chr1_part1_chr'} = $ref_chr2;
        $translocation{$translocation_id}{'ref_chr2_chr1_part1_start'} = 1;
        $translocation{$translocation_id}{'ref_chr2_chr1_part1_end'} = $ref_chr2_breakpoint_start - 1;
        $translocation{$translocation_id}{'ref_chr2_chr1_part1_strand'} = "+";

        $translocation{$translocation_id}{'ref_chr2_chr1_part2_chr'} = $ref_chr1;
        $translocation{$translocation_id}{'ref_chr2_chr1_part2_start'} = $ref_chr1_breakpoint_start;
        $translocation{$translocation_id}{'ref_chr2_chr1_part2_end'} = length $$refseq_hashref{$ref_chr1};
        $translocation{$translocation_id}{'ref_chr2_chr1_part2_strand'} = "+";
      } else {
        # translocation_resolution: +--+
        $translocation{$translocation_id}{'ref_chr1'} = $ref_chr1;
        $translocation{$translocation_id}{'ref_chr2'} = $ref_chr2;

        $translocation{$translocation_id}{'ref_chr1_chr2_part1_chr'} = $ref_chr1;
        $translocation{$translocation_id}{'ref_chr1_chr2_part1_start'} = 1;
        $translocation{$translocation_id}{'ref_chr1_chr2_part1_end'} = $ref_chr1_breakpoint_start - 1;
        $translocation{$translocation_id}{'ref_chr1_chr2_part1_strand'} = "+";

        $translocation{$translocation_id}{'ref_chr1_chr2_part2_chr'} = $ref_chr2;
        $translocation{$translocation_id}{'ref_chr1_chr2_part2_start'} = 1;
        $translocation{$translocation_id}{'ref_chr1_chr2_part2_end'} = $ref_chr2_breakpoint_end;
        $translocation{$translocation_id}{'ref_chr1_chr2_part2_strand'} = "-";

        $translocation{$translocation_id}{'ref_chr2_chr1_part1_chr'} = $ref_chr2;
        $translocation{$translocation_id}{'ref_chr2_chr1_part1_start'} = $ref_chr2_breakpoint_end + 1;
        $translocation{$translocation_id}{'ref_chr2_chr1_part1_end'} = length $$refseq_hashref{$ref_chr2};
        $translocation{$translocation_id}{'ref_chr2_chr1_part1_strand'} = "-";

        $translocation{$translocation_id}{'ref_chr2_chr1_part2_chr'} = $ref_chr1;
        $translocation{$translocation_id}{'ref_chr2_chr1_part2_start'} = $ref_chr1_breakpoint_start;
        $translocation{$translocation_id}{'ref_chr2_chr1_part2_end'} = length $$refseq_hashref{$ref_chr1};
        $translocation{$translocation_id}{'ref_chr2_chr1_part2_strand'} = "+";
      }
      # delete used breakpoints
      # delete $valid_breakpoint_pair{$breakpoint1}{$breakpoint2};
      # delete $valid_breakpoint_pair{$breakpoint2}{$breakpoint1};
      delete $valid_breakpoint_pair{$breakpoint1};
      delete $valid_breakpoint_pair{$breakpoint2};
      # update avialble_chr
      delete $available_chr{$ref_chr1};
      delete $available_chr{$ref_chr2};
    }
  } else {
    for (my $i = 1; $i <= $translocation_count; $i++) {
      # print "sample translocation: i=$i\n";
      # check if there are still valid $breakpoints
      if ((scalar (keys %available_chr)) < 2) {
        print "\n!!! Warning! No more available chromosomes in the current simulation.\n";
        print "!!! Only $i translocation were introduced.\n";
        last;
      }
      # random sampling across the genome
      my %refseq_genome_space = create_genome_space($refseq_hashref);
      SAMPLE_RANDOM_TRA2:
      my ($sample1_chr, $sample1_breakpoint_start) = sample_genome_space(\%refseq_genome_space);
      my ($sample2_chr, $sample2_breakpoint_start) = sample_genome_space(\%refseq_genome_space);

      if ($sample1_chr eq $sample2_chr) {
        goto SAMPLE_RANDOM_TRA2;
      } elsif ((not exists $available_chr{$sample1_chr}) or (not exists $available_chr{$sample2_chr})) {
        goto SAMPLE_RANDOM_TRA2;
      } else {
        my $sample1_breakpoint_end = $sample1_breakpoint_start;
        my $sample2_breakpoint_end = $sample2_breakpoint_start;

        # check if the sampled region overlaped with the defined centromeres
        if (defined $centromere_gff) {
          my $sample1_centromere_check_flag = 0;
          my $sample2_centromere_check_flag = 0;
          foreach my $centromere_id (sort keys %centromere) {
            if ($centromere{$centromere_id}{'chr'} eq $sample1_chr) {
              $sample1_centromere_check_flag = check_overlap_region($sample1_breakpoint_start, $sample1_breakpoint_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
            }
            if ($centromere{$centromere_id}{'chr'} eq $sample2_chr) {
              $sample2_centromere_check_flag = check_overlap_region($sample2_breakpoint_start, $sample2_breakpoint_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
            }
            if (($sample1_centromere_check_flag == 1) or ($sample2_centromere_check_flag == 1)) {
              goto SAMPLE_RANDOM_TRA2;
            }
          }
        }

        my $ref_chr1;
        my $ref_chr2;
        my $ref_chr1_breakpoint_start;
        my $ref_chr1_breakpoint_end;
        my $ref_chr2_breakpoint_start;
        my $ref_chr2_breakpoint_end;
        if ($sample1_chr lt $sample2_chr) {
          $ref_chr1 = $sample1_chr;
          $ref_chr1_breakpoint_start = $sample1_breakpoint_start;
          $ref_chr1_breakpoint_end = $sample1_breakpoint_end;
          $ref_chr2 = $sample2_chr;
          $ref_chr2_breakpoint_start = $sample2_breakpoint_start;
          $ref_chr2_breakpoint_end = $sample2_breakpoint_end;
        } else {
          $ref_chr1 = $sample2_chr;
          $ref_chr1_breakpoint_start = $sample2_breakpoint_start;
          $ref_chr1_breakpoint_end = $sample2_breakpoint_end;
          $ref_chr2 = $sample1_chr;
          $ref_chr2_breakpoint_start = $sample1_breakpoint_start;
          $ref_chr2_breakpoint_end = $sample1_breakpoint_end;
        }


        if (rand(1) < 0.5) {
          # try translocation resolution: ++++
          my $ref_chr1_chr2_part1_chr = $ref_chr1;
          my $ref_chr1_chr2_part1_start = 1;
          my $ref_chr1_chr2_part1_end = $ref_chr1_breakpoint_start - 1;
          my $ref_chr1_chr2_part1_strand = "+";

          my $ref_chr1_chr2_part2_chr = $ref_chr2;
          my $ref_chr1_chr2_part2_start = $ref_chr2_breakpoint_start;
          my $ref_chr1_chr2_part2_end = length $$refseq_hashref{$ref_chr2};
          my $ref_chr1_chr2_part2_strand = "+";

          my $ref_chr2_chr1_part1_chr = $ref_chr2;
          my $ref_chr2_chr1_part1_start = 1;
          my $ref_chr2_chr1_part1_end = $ref_chr2_breakpoint_start - 1;
          my $ref_chr2_chr1_part1_strand = "+";

          my $ref_chr2_chr1_part2_chr = $ref_chr1;
          my $ref_chr2_chr1_part2_start = $ref_chr1_breakpoint_start;
          my $ref_chr2_chr1_part2_end = length $$refseq_hashref{$ref_chr1};
          my $ref_chr2_chr1_part2_strand = "+";

          # check for rearranged chromosome with zero or two centromeres
          if (defined $centromere_gff) {
            my $ref_chr1_chr2_part1_centromere_count = 0;
            my $ref_chr1_chr2_part2_centromere_count = 0;
            my $ref_chr2_chr1_part1_centromere_count = 0;
            my $ref_chr2_chr1_part2_centromere_count = 0;
            foreach my $centromere_id (sort keys %centromere) {
              if ($centromere{$centromere_id}{'chr'} eq $ref_chr1) {
                $ref_chr1_chr2_part1_centromere_count = check_overlap_region($ref_chr1_chr2_part1_start, $ref_chr1_chr2_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                $ref_chr2_chr1_part2_centromere_count = check_overlap_region($ref_chr2_chr1_part2_start, $ref_chr2_chr1_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
              }
              if ($centromere{$centromere_id}{'chr'} eq $ref_chr2) {
                $ref_chr1_chr2_part2_centromere_count = check_overlap_region($ref_chr1_chr2_part2_start, $ref_chr1_chr2_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                $ref_chr2_chr1_part1_centromere_count = check_overlap_region($ref_chr2_chr1_part1_start, $ref_chr2_chr1_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
              }
            }
            my $ref_chr1_chr2_centromere_count = $ref_chr1_chr2_part1_centromere_count + $ref_chr1_chr2_part2_centromere_count;
            my $ref_chr2_chr1_centromere_count = $ref_chr2_chr1_part1_centromere_count + $ref_chr2_chr1_part2_centromere_count;

            if (($ref_chr1_chr2_centromere_count != 1) or ($ref_chr2_chr1_centromere_count != 1)) {
              goto SAMPLE_RANDOM_TRA2;
            }
          }

          my $translocation_id = "TRA_${i}";
          $translocation{$translocation_id}{'ref_chr1'} = $ref_chr1;
          $translocation{$translocation_id}{'ref_chr2'} = $ref_chr2;

          $translocation{$translocation_id}{'ref_chr1_chr2_part1_chr'} = $ref_chr1_chr2_part1_chr;
          $translocation{$translocation_id}{'ref_chr1_chr2_part1_start'} = $ref_chr1_chr2_part1_start;
          $translocation{$translocation_id}{'ref_chr1_chr2_part1_end'} = $ref_chr1_chr2_part1_end;
          $translocation{$translocation_id}{'ref_chr1_chr2_part1_strand'} = $ref_chr1_chr2_part1_strand;

          $translocation{$translocation_id}{'ref_chr1_chr2_part2_chr'} = $ref_chr1_chr2_part2_chr;
          $translocation{$translocation_id}{'ref_chr1_chr2_part2_start'} = $ref_chr1_chr2_part2_start;
          $translocation{$translocation_id}{'ref_chr1_chr2_part2_end'} = $ref_chr1_chr2_part2_end;
          $translocation{$translocation_id}{'ref_chr1_chr2_part2_strand'} = $ref_chr1_chr2_part2_strand;

          $translocation{$translocation_id}{'ref_chr2_chr1_part1_chr'} = $ref_chr2_chr1_part1_chr;
          $translocation{$translocation_id}{'ref_chr2_chr1_part1_start'} = $ref_chr2_chr1_part1_start;
          $translocation{$translocation_id}{'ref_chr2_chr1_part1_end'} = $ref_chr2_chr1_part1_end;
          $translocation{$translocation_id}{'ref_chr2_chr1_part1_strand'} = $ref_chr2_chr1_part1_strand;

          $translocation{$translocation_id}{'ref_chr2_chr1_part2_chr'} = $ref_chr2_chr1_part2_chr;
          $translocation{$translocation_id}{'ref_chr2_chr1_part2_start'} = $ref_chr2_chr1_part2_start;
          $translocation{$translocation_id}{'ref_chr2_chr1_part2_end'} = $ref_chr2_chr1_part2_end;
          $translocation{$translocation_id}{'ref_chr2_chr1_part2_strand'} = $ref_chr2_chr1_part2_strand;
          # update avialble_chr
          delete $available_chr{$ref_chr1};
          delete $available_chr{$ref_chr2};

        } else {
          # translocation resolution: +--+
          my $ref_chr1_chr2_part1_chr = $ref_chr1;
          my $ref_chr1_chr2_part1_start = 1;
          my $ref_chr1_chr2_part1_end = $ref_chr1_breakpoint_start - 1;
          my $ref_chr1_chr2_part1_strand = "+";

          my $ref_chr1_chr2_part2_chr = $ref_chr2;
          my $ref_chr1_chr2_part2_start = 1;
          my $ref_chr1_chr2_part2_end = $ref_chr2_breakpoint_end;
          my $ref_chr1_chr2_part2_strand = "-";

          my $ref_chr2_chr1_part1_chr = $ref_chr2;
          my $ref_chr2_chr1_part1_start = $ref_chr2_breakpoint_end + 1;
          my $ref_chr2_chr1_part1_end = length $$refseq_hashref{$ref_chr2};
          my $ref_chr2_chr1_part1_strand = "-";

          my $ref_chr2_chr1_part2_chr = $ref_chr1;
          my $ref_chr2_chr1_part2_start = $ref_chr1_breakpoint_start;
          my $ref_chr2_chr1_part2_end = length $$refseq_hashref{$ref_chr1};
          my $ref_chr2_chr1_part2_strand = "+";

          # check for rearranged chromosome with zero or two centromeres
          if (defined $centromere_gff) {
            my $ref_chr1_chr2_part1_centromere_count = 0;
            my $ref_chr1_chr2_part2_centromere_count = 0;
            my $ref_chr2_chr1_part1_centromere_count = 0;
            my $ref_chr2_chr1_part2_centromere_count = 0;
            foreach my $centromere_id (sort keys %centromere) {
              if ($centromere{$centromere_id}{'chr'} eq $ref_chr1) {
                $ref_chr1_chr2_part1_centromere_count = check_overlap_region($ref_chr1_chr2_part1_start, $ref_chr1_chr2_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                $ref_chr2_chr1_part2_centromere_count = check_overlap_region($ref_chr2_chr1_part2_start, $ref_chr2_chr1_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
              }
              if ($centromere{$centromere_id}{'chr'} eq $ref_chr2) {
                $ref_chr1_chr2_part2_centromere_count = check_overlap_region($ref_chr1_chr2_part2_start, $ref_chr1_chr2_part2_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
                $ref_chr2_chr1_part1_centromere_count = check_overlap_region($ref_chr2_chr1_part1_start, $ref_chr2_chr1_part1_end, $centromere{$centromere_id}{'start'}, $centromere{$centromere_id}{'end'});
              }
            }
            my $ref_chr1_chr2_centromere_count = $ref_chr1_chr2_part1_centromere_count + $ref_chr1_chr2_part2_centromere_count;
            my $ref_chr2_chr1_centromere_count = $ref_chr2_chr1_part1_centromere_count + $ref_chr2_chr1_part2_centromere_count;

            if (($ref_chr1_chr2_centromere_count != 1) or ($ref_chr2_chr1_centromere_count != 1)) {
              goto SAMPLE_RANDOM_TRA2;
            }
          }

          my $translocation_id = "TRA_${i}";
          $translocation{$translocation_id}{'ref_chr1'} = $ref_chr1;
          $translocation{$translocation_id}{'ref_chr2'} = $ref_chr2;

          $translocation{$translocation_id}{'ref_chr1_chr2_part1_chr'} = $ref_chr1_chr2_part1_chr;
          $translocation{$translocation_id}{'ref_chr1_chr2_part1_start'} = $ref_chr1_chr2_part1_start;
          $translocation{$translocation_id}{'ref_chr1_chr2_part1_end'} = $ref_chr1_chr2_part1_end;
          $translocation{$translocation_id}{'ref_chr1_chr2_part1_strand'} = $ref_chr1_chr2_part1_strand;

          $translocation{$translocation_id}{'ref_chr1_chr2_part2_chr'} = $ref_chr1_chr2_part2_chr;
          $translocation{$translocation_id}{'ref_chr1_chr2_part2_start'} = $ref_chr1_chr2_part2_start;
          $translocation{$translocation_id}{'ref_chr1_chr2_part2_end'} = $ref_chr1_chr2_part2_end;
          $translocation{$translocation_id}{'ref_chr1_chr2_part2_strand'} = $ref_chr1_chr2_part2_strand;

          $translocation{$translocation_id}{'ref_chr2_chr1_part1_chr'} = $ref_chr2_chr1_part1_chr;
          $translocation{$translocation_id}{'ref_chr2_chr1_part1_start'} = $ref_chr2_chr1_part1_start;
          $translocation{$translocation_id}{'ref_chr2_chr1_part1_end'} = $ref_chr2_chr1_part1_end;
          $translocation{$translocation_id}{'ref_chr2_chr1_part1_strand'} = $ref_chr2_chr1_part1_strand;

          $translocation{$translocation_id}{'ref_chr2_chr1_part2_chr'} = $ref_chr2_chr1_part2_chr;
          $translocation{$translocation_id}{'ref_chr2_chr1_part2_start'} = $ref_chr2_chr1_part2_start;
          $translocation{$translocation_id}{'ref_chr2_chr1_part2_end'} = $ref_chr2_chr1_part2_end;
          $translocation{$translocation_id}{'ref_chr2_chr1_part2_strand'} = $ref_chr2_chr1_part2_strand;

          # update avialble_chr
          delete $available_chr{$ref_chr1};
          delete $available_chr{$ref_chr2};
        }
      }
    }
  }

  foreach my $translocation_id (sort keys %translocation) {
    my $ref_chr1 = $translocation{$translocation_id}{'ref_chr1'};
    my $ref_chr2 = $translocation{$translocation_id}{'ref_chr2'};

    my $ref_chr1_chr2_part1_chr = $translocation{$translocation_id}{'ref_chr1_chr2_part1_chr'};
    my $ref_chr1_chr2_part1_start = $translocation{$translocation_id}{'ref_chr1_chr2_part1_start'};
    my $ref_chr1_chr2_part1_end = $translocation{$translocation_id}{'ref_chr1_chr2_part1_end'};
    my $ref_chr1_chr2_part1_strand = $translocation{$translocation_id}{'ref_chr1_chr2_part1_strand'};

    my $ref_chr1_chr2_part2_chr = $translocation{$translocation_id}{'ref_chr1_chr2_part2_chr'};
    my $ref_chr1_chr2_part2_start = $translocation{$translocation_id}{'ref_chr1_chr2_part2_start'};
    my $ref_chr1_chr2_part2_end = $translocation{$translocation_id}{'ref_chr1_chr2_part2_end'};
    my $ref_chr1_chr2_part2_strand = $translocation{$translocation_id}{'ref_chr1_chr2_part2_strand'};

    my $ref_chr2_chr1_part1_chr = $translocation{$translocation_id}{'ref_chr2_chr1_part1_chr'};
    my $ref_chr2_chr1_part1_start = $translocation{$translocation_id}{'ref_chr2_chr1_part1_start'};
    my $ref_chr2_chr1_part1_end = $translocation{$translocation_id}{'ref_chr2_chr1_part1_end'};
    my $ref_chr2_chr1_part1_strand = $translocation{$translocation_id}{'ref_chr2_chr1_part1_strand'};

    my $ref_chr2_chr1_part2_chr = $translocation{$translocation_id}{'ref_chr2_chr1_part2_chr'};
    my $ref_chr2_chr1_part2_start = $translocation{$translocation_id}{'ref_chr2_chr1_part2_start'};
    my $ref_chr2_chr1_part2_end = $translocation{$translocation_id}{'ref_chr2_chr1_part2_end'};
    my $ref_chr2_chr1_part2_strand = $translocation{$translocation_id}{'ref_chr2_chr1_part2_strand'};

    my $ref_chr1_chr2 = "${ref_chr1}_${ref_chr2}";
    my $ref_chr2_chr1 = "${ref_chr2}_${ref_chr1}";

    if ($ref_chr1_chr2_part1_strand eq $ref_chr1_chr2_part2_strand) {
      $ref_chr1_chr2_part2_end = length $$refseq_hashref{$ref_chr2};
      $ref_chr2_chr1_part2_end = length $$refseq_hashref{$ref_chr1};
    } else {
      $ref_chr2_chr1_part1_end = length $$refseq_hashref{$ref_chr2};
      $ref_chr2_chr1_part2_end = length $$refseq_hashref{$ref_chr1};
    }

    my $ref_chr1_chr2_part1_seq = substr $$refseq_hashref{$ref_chr1}, $ref_chr1_chr2_part1_start - 1, $ref_chr1_chr2_part1_end - $ref_chr1_chr2_part1_start + 1;
    my $ref_chr1_chr2_part2_seq = substr $$refseq_hashref{$ref_chr2}, $ref_chr1_chr2_part2_start - 1, $ref_chr1_chr2_part2_end - $ref_chr1_chr2_part2_start + 1;
    my $ref_chr2_chr1_part1_seq = substr $$refseq_hashref{$ref_chr2}, $ref_chr2_chr1_part1_start - 1, $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1;
    my $ref_chr2_chr1_part2_seq = substr $$refseq_hashref{$ref_chr1}, $ref_chr2_chr1_part2_start - 1, $ref_chr2_chr1_part2_end - $ref_chr2_chr1_part2_start + 1;

    if ($ref_chr1_chr2_part1_strand ne $ref_chr1_chr2_part2_strand) {
      $ref_chr1_chr2_part2_seq = revcom($ref_chr1_chr2_part2_seq);
      $ref_chr2_chr1_part1_seq = revcom($ref_chr2_chr1_part1_seq);
    }

    $$simseq_hashref{$ref_chr1_chr2} = $ref_chr1_chr2_part1_seq . $ref_chr1_chr2_part2_seq;
    $$simseq_hashref{$ref_chr2_chr1} = $ref_chr2_chr1_part1_seq . $ref_chr2_chr1_part2_seq;
    delete $$simseq_hashref{$ref_chr1};
    delete $$simseq_hashref{$ref_chr2};

    if ($ref_chr1_chr2_part1_strand eq $ref_chr1_chr2_part2_strand) {
      # translocation resolution: ++++
      # ref_chr1_chr2_part1
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_start'} = $ref_chr1_chr2_part1_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_strand'} = $ref_chr1_chr2_part1_strand;
      my $ref_chr1_chr2_part1_end_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr1_chr2_part1_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_allele'} = $ref_chr1_chr2_part1_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_strand'} = $ref_chr1_chr2_part1_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_allele'} = "$ref_chr1_chr2_part1_end_ref_allele\[$ref_chr2:$ref_chr1_chr2_part2_start\[";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_id'} = $translocation_id;
      # ref_chr1_chr2_part2
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_start'} = $ref_chr1_chr2_part2_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_end'} = $ref_chr1_chr2_part2_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_strand'} = $ref_chr1_chr2_part2_strand;
      my $ref_chr1_chr2_part2_start_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr1_chr2_part2_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'ref_allele'} = $ref_chr1_chr2_part2_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_start'} = $ref_chr1_chr2_part1_end + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_end'} = $ref_chr1_chr2_part1_end + $ref_chr1_chr2_part2_end - $ref_chr1_chr2_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_strand'} = $ref_chr1_chr2_part2_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'sim_allele'} = "\]$ref_chr1:$ref_chr1_chr2_part1_end\]$ref_chr1_chr2_part2_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_start}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part1
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_start'} = $ref_chr2_chr1_part1_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_end'} = $ref_chr2_chr1_part1_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_strand'} = $ref_chr2_chr1_part1_strand;
      my $ref_chr2_chr1_part1_end_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr2_chr1_part1_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'ref_allele'} = $ref_chr2_chr1_part1_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_end'} = $ref_chr2_chr1_part1_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_strand'} = $ref_chr2_chr1_part1_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'sim_allele'} = "$ref_chr2_chr1_part1_end_ref_allele\[$ref_chr1:$ref_chr2_chr1_part2_start\[";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_end}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part2
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_start'} = $ref_chr2_chr1_part2_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_end'} = $ref_chr2_chr1_part2_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_strand'} = $ref_chr2_chr1_part2_strand;
      my $ref_chr2_chr1_part2_start_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr2_chr1_part2_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_allele'} = $ref_chr2_chr1_part2_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_start'} = $ref_chr2_chr1_part1_end + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_end'} = $ref_chr2_chr1_part1_end + $ref_chr2_chr1_part2_end - $ref_chr2_chr1_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_strand'} = $ref_chr2_chr1_part2_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_allele'} = "\]$ref_chr2:$ref_chr2_chr1_part1_end\]$ref_chr2_chr1_part2_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_id'} = $translocation_id;
    } else {
      # translocation resolution: +--+
      # ref_chr1_chr2_part1
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_start'} = $ref_chr1_chr2_part1_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_strand'} = $ref_chr1_chr2_part1_strand;
      my $ref_chr1_chr2_part1_end_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr1_chr2_part1_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'ref_allele'} = $ref_chr1_chr2_part1_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_end'} = $ref_chr1_chr2_part1_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_strand'} = $ref_chr1_chr2_part1_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'sim_allele'} = "$ref_chr1_chr2_part1_end_ref_allele\]$ref_chr2:$ref_chr1_chr2_part2_end\]";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr1_chr2_part1_end}{'variant_id'} = $translocation_id;
      # ref_chr1_chr2_part2
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_start'} = $ref_chr1_chr2_part2_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_end'} = $ref_chr1_chr2_part2_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_strand'} = $ref_chr1_chr2_part2_strand;
      my $ref_chr1_chr2_part2_end_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr1_chr2_part2_end - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'ref_allele'} = $ref_chr1_chr2_part2_end_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_chr'} = $ref_chr1_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_start'} = $ref_chr1_chr2_part1_end + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_end'} = $ref_chr1_chr2_part1_end + $ref_chr1_chr2_part2_end - $ref_chr1_chr2_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_strand'} = $ref_chr1_chr2_part2_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'sim_allele'} = "$ref_chr1_chr2_part2_end_ref_allele\]$ref_chr1:$ref_chr1_chr2_part1_end\]";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr1_chr2_part2_end}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part1
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_chr'} = $ref_chr2;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_start'} = $ref_chr2_chr1_part1_start;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_end'} = $ref_chr2_chr1_part1_end;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_strand'} = $ref_chr2_chr1_part1_strand;
      my $ref_chr2_chr1_part1_start_ref_allele = substr $$refseq_hashref{$ref_chr2}, $ref_chr2_chr1_part1_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'ref_allele'} = $ref_chr2_chr1_part1_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_start'} = 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_end'} = $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_strand'} = $ref_chr2_chr1_part1_strand;
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'sim_allele'} = "\[$ref_chr1:$ref_chr2_chr1_part2_start\[$ref_chr2_chr1_part1_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr2}{$ref_chr2_chr1_part1_start}{'variant_id'} = $translocation_id;
      # ref_chr2_chr1_part2
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_chr'} = $ref_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_start'} = $ref_chr2_chr1_part2_start;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_end'} = $ref_chr2_chr1_part2_end;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_strand'} = $ref_chr2_chr1_part2_strand;
      my $ref_chr2_chr1_part2_start_ref_allele = substr $$refseq_hashref{$ref_chr1}, $ref_chr2_chr1_part2_start - 1, 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'ref_allele'} = $ref_chr2_chr1_part2_start_ref_allele;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_chr'} = $ref_chr2_chr1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_start'} = $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1 + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_end'} = $ref_chr2_chr1_part1_end - $ref_chr2_chr1_part1_start + 1 + $ref_chr2_chr1_part2_end - $ref_chr2_chr1_part2_start + 1;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_strand'} = $ref_chr2_chr1_part2_strand;
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'sim_allele'} = "\[$ref_chr2:$ref_chr2_chr1_part1_start\[$ref_chr2_chr1_part2_start_ref_allele";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_type'} = "TRA";
      $$ref2sim_map_hashref{$ref_chr1}{$ref_chr2_chr1_part2_start}{'variant_id'} = $translocation_id;
    }
  }
}


sub generate_output_files {
  my ($prefix, $refseq_arrayref, $refseq_hashref, $simseq_hashref, $ref2sim_map_hashref) = @_;
  # output fasta.gz file for the simulated genome
  my $output_simseq_fasta = "$prefix.simseq.genome.fa";
  my $output_simseq_fasta_fh = write_file($output_simseq_fasta);
  if ((defined $translocation_vcf) or (defined $translocation_count)) {
    foreach my $chr (sort keys %$simseq_hashref) {
      if (exists $$simseq_hashref{$chr}) {
        print $output_simseq_fasta_fh ">$chr\n$$simseq_hashref{$chr}\n";
      }
    }
  } else {
    foreach my $chr (@$refseq_arrayref) {
      if (exists $$simseq_hashref{$chr}) {
        print $output_simseq_fasta_fh ">$chr\n$$simseq_hashref{$chr}\n";
      }
    }
  }
  # generate the correspondance map for genomic variants introduced during simulation
  print "Generating the correspondance map for genomic variants introduced during simulation:\n";
  print "$prefix.refseq2simseq.map.txt\n";
  my $output_ref2sim_map = "$prefix.refseq2simseq.map.txt";
  my $output_ref2sim_map_fh = write_file($output_ref2sim_map);
  print $output_ref2sim_map_fh "ref_chr\tref_start\tref_end\tref_strand\tref_allele\tsim_chr\tsim_start\tsim_end\tref_strand\tsim_allele\tvariant_type\tvariant_id\tdonor_chr_in_ref\tdonor_start_in_ref\tdonor_end_in_ref\tdonor_strand_in_ref\n";

  if ((defined $translocation_vcf) or (defined $translocation_count)) {
    foreach my $ref_chr (@$refseq_arrayref) {
      if (not exists $$simseq_hashref{$ref_chr}) {
        foreach my $ref_start (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          print $output_ref2sim_map_fh "$ref_chr\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'}\t";
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} eq "DUP") {
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_chr_in_ref'}\t";
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_start_in_ref'}\t";
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_end_in_ref'}\t";
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_strand_in_ref'}\n";
          } else {
            print $output_ref2sim_map_fh ".\t.\t.\t.\n";
          }
        }
      }
    }
  } else {
    foreach my $ref_chr (@$refseq_arrayref) {
      if (exists $$simseq_hashref{$ref_chr}) {
        foreach my $ref_start (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          print $output_ref2sim_map_fh "$ref_chr\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_start'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_strand'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_strand'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'}\t";
          print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'}\t";
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} eq "DUP") {
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_chr_in_ref'}\t";
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_start_in_ref'}\t";
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_end_in_ref'}\t";
            print $output_ref2sim_map_fh "$$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_strand_in_ref'}\n";
          } else {
            print $output_ref2sim_map_fh ".\t.\t.\t.\n";
          }
        }
      }
    }
  }

  # generate reference-based vcf file for genomic variants introduced during simulation
  print "\nGenerating reference-based vcf file for genomic variants introduced during simulation:\n";
  my $gmt_time = gmtime();
  if ((defined $snp_vcf) or (defined $snp_count)) {
    print "$prefix.refseq2simseq.SNP.vcf\n";
    my $output_ref2sim_snp_vcf = "$prefix.refseq2simseq.SNP.vcf";
    my $output_ref2sim_snp_vcf_fh = write_file($output_ref2sim_snp_vcf);
    print $output_ref2sim_snp_vcf_fh "##fileformat=VCFv4.1\n";
    print $output_ref2sim_snp_vcf_fh "##fileDate=$gmt_time (GMT time)\n";
    print $output_ref2sim_snp_vcf_fh "##source=simuG.pl\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=variant_type,Number=1,Type=String,Description=\"The type of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=variant_id,Number=1,Type=String,Description=\"The id of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=ref_chr,Number=1,Type=String,Description=\"The reference genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=ref_start,Number=1,Type=String,Description=\"The reference genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=ref_end,Number=1,Type=String,Description=\"The reference genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=sim_chr,Number=1,Type=String,Description=\"The simulated genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=sim_start,Number=1,Type=String,Description=\"The simulated genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "##INFO=<ID=sim_end,Number=1,Type=String,Description=\"The simulated genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_snp_vcf_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    foreach my $ref_chr (@$refseq_arrayref) {
      if (exists $$simseq_hashref{$ref_chr}) {
        foreach my $ref_start (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} eq "SNP") {
            my $ref_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'};
            my $ref_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'};
            my $sim_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'};
            my $sim_chr = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'};
            my $sim_start = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'};
            my $sim_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'};
            my $variant_type = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'};
            my $variant_id = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'};
            print $output_ref2sim_snp_vcf_fh "$ref_chr\t$ref_start\t.\t";
            print $output_ref2sim_snp_vcf_fh "$ref_allele\t$sim_allele\t.\t.\tvariant_type=$variant_type;ref_chr=$ref_chr;ref_start=$ref_start;ref_end=$ref_end;sim_chr=$sim_chr;sim_start=$sim_start;sim_end=$sim_end\n";
          }
        }
      }
    }
  }
  if ((defined $indel_vcf) or (defined $indel_count)) {
    print "$prefix.refseq2simseq.INDEL.vcf\n\n";
    my $output_ref2sim_indel_vcf = "$prefix.refseq2simseq.INDEL.vcf";
    my $output_ref2sim_indel_vcf_fh = write_file($output_ref2sim_indel_vcf);
    print $output_ref2sim_indel_vcf_fh "##fileformat=VCFv4.1\n";
    print $output_ref2sim_indel_vcf_fh "##fileDate=$gmt_time (GMT time)\n";
    print $output_ref2sim_indel_vcf_fh "##source=simuG.pl\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=variant_type,Number=1,Type=String,Description=\"The type of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=variant_id,Number=1,Type=String,Description=\"The id of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=ref_chr,Number=1,Type=String,Description=\"The reference genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=ref_start,Number=1,Type=String,Description=\"The reference genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=ref_end,Number=1,Type=String,Description=\"The reference genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=sim_chr,Number=1,Type=String,Description=\"The simulated genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=sim_start,Number=1,Type=String,Description=\"The simulated genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "##INFO=<ID=sim_end,Number=1,Type=String,Description=\"The simulated genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_indel_vcf_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    foreach my $ref_chr (@$refseq_arrayref) {
      if (exists $$simseq_hashref{$ref_chr}) {
        foreach my $ref_start (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} eq "INDEL") {
            my $ref_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'};
            my $ref_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'};
            my $sim_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'};
            my $sim_chr = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'};
            my $sim_start = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'};
            my $sim_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'};
            my $variant_type = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'};
            my $variant_id = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'};
            print $output_ref2sim_indel_vcf_fh "$ref_chr\t$ref_start\t.\t";
            print $output_ref2sim_indel_vcf_fh "$ref_allele\t$sim_allele\t.\t.\tvariant_type=$variant_type;ref_chr=$ref_chr;ref_start=$ref_start;ref_end=$ref_end;sim_chr=$sim_chr;sim_start=$sim_start;sim_end=$sim_end\n";
          }
        }
      }
    }
  }
  if ((defined $cnv_vcf) or (defined $cnv_count)) {
    print "$prefix.refseq2simseq.CNV.vcf\n\n";
    my $output_ref2sim_cnv_vcf = "$prefix.refseq2simseq.CNV.vcf";
    my $output_ref2sim_cnv_vcf_fh = write_file($output_ref2sim_cnv_vcf);
    print $output_ref2sim_cnv_vcf_fh "##fileformat=VCFv4.1\n";
    print $output_ref2sim_cnv_vcf_fh "##fileDate=$gmt_time (GMT time)\n";
    print $output_ref2sim_cnv_vcf_fh "##source=simuG.pl\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=variant_type,Number=1,Type=String,Description=\"The type of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=variant_id,Number=1,Type=String,Description=\"The id of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=ref_chr,Number=1,Type=String,Description=\"The reference genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=ref_start,Number=1,Type=String,Description=\"The reference genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=ref_end,Number=1,Type=String,Description=\"The reference genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=sim_chr,Number=1,Type=String,Description=\"The simulated genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=sim_start,Number=1,Type=String,Description=\"The simulated genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=sim_end,Number=1,Type=String,Description=\"The simulated genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=donor_chr_in_ref,Number=1,Type=String,Description=\"The reference genome based chromosome ID of the donor region in the case of a segmental duplication\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=donor_start_in_ref,Number=1,Type=String,Description=\"The reference genome based start coordinate of the donor region in the case of a segmental duplication\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=donor_end_in_ref,Number=1,Type=String,Description=\"The reference genome based end coordinate of the donor region in the case of a segmental duplication\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=donor_strand_in_ref,Number=1,Type=String,Description=\"The strand of the donor region in the case of a segmental duplication\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"The event id of the introduced genomic variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"The type of strucrtural variant\">\n";
    print $output_ref2sim_cnv_vcf_fh "##INFO=<ID=END,Number=1,Type=String,Description=\"The end position of the variant described in this record\">\n";

    print $output_ref2sim_cnv_vcf_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    foreach my $ref_chr (@$refseq_arrayref) {
      if (exists $$simseq_hashref{$ref_chr}) {
        foreach my $ref_start (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} =~ /(DEL|DUP)/) {
            my $ref_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'};
            my $ref_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'};
            my $sim_chr = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'};
            my $sim_start = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'};
            my $sim_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'};
            my $sim_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'};
            my $variant_type = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'};
            my $variant_id = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'};
            if ($variant_type eq "DEL") {
              print $output_ref2sim_cnv_vcf_fh "$ref_chr\t$ref_start\t.\t";
              print $output_ref2sim_cnv_vcf_fh "$ref_allele\t$sim_allele\t.\t\tSVTYPE=DEL;EVENT=$variant_id;END=$ref_end\n";
            } else {
              my $donor_chr_in_ref = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_chr_in_ref'};
              my $donor_start_in_ref = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_start_in_ref'};
              my $donor_end_in_ref = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_end_in_ref'};
              my $donor_strand_in_ref = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'donor_strand_in_ref'};
              my $ref_start_after_duplication = $ref_start + 1;
              my $ref_allele_after_duplication = substr $$refseq_hashref{$ref_chr}, $ref_end, 1;
              if ($donor_strand_in_ref eq "+") {
                print $output_ref2sim_cnv_vcf_fh "$ref_chr\t$ref_start\t.\t";
                print $output_ref2sim_cnv_vcf_fh "$ref_allele\t$ref_allele\[$donor_chr_in_ref:$donor_start_in_ref\[\t.\t.\tSVTYPE=BND;EVENT=$variant_id\n";

                print $output_ref2sim_cnv_vcf_fh "$ref_chr\t$ref_start_after_duplication\t.\t";
                print $output_ref2sim_cnv_vcf_fh "$ref_allele_after_duplication\t\]$donor_chr_in_ref:$donor_end_in_ref\]$ref_allele_after_duplication\t.\t.\tSVTYPE=BND;EVENT=$variant_id\n";
              } else {
                print $output_ref2sim_cnv_vcf_fh "$ref_chr\t$ref_start\t.\t";
                print $output_ref2sim_cnv_vcf_fh "$ref_allele\t$ref_allele\]$donor_chr_in_ref:$donor_start_in_ref\]\t.\t.\tSVTYPE=BND;EVENT=$variant_id\n";
                print $output_ref2sim_cnv_vcf_fh "$ref_chr\t$ref_start_after_duplication\t.\t";
                print $output_ref2sim_cnv_vcf_fh "$ref_allele_after_duplication\t\]$donor_chr_in_ref:$donor_end_in_ref\]$ref_allele_after_duplication\t.\t.\tSVTYPE=BND;EVENT=$variant_id\n";
              }
            }
          }
        }
      }
    }
  }
  if ((defined $inversion_vcf) or (defined $inversion_count)) {
    print "$prefix.refseq2simseq.inversion.vcf\n\n";
    my $output_ref2sim_inversion_vcf = "$prefix.refseq2simseq.inversion.vcf";
    my $output_ref2sim_inversion_vcf_fh = write_file($output_ref2sim_inversion_vcf);
    print $output_ref2sim_inversion_vcf_fh "##fileformat=VCFv4.1\n";
    print $output_ref2sim_inversion_vcf_fh "##fileDate=$gmt_time (GMT time)\n";
    print $output_ref2sim_inversion_vcf_fh "##source=simuG.pl\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=variant_type,Number=1,Type=String,Description=\"The type of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=variant_id,Number=1,Type=String,Description=\"The id of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=ref_chr,Number=1,Type=String,Description=\"The reference genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=ref_start,Number=1,Type=String,Description=\"The reference genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=ref_end,Number=1,Type=String,Description=\"The reference genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=sim_chr,Number=1,Type=String,Description=\"The simulated genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=sim_start,Number=1,Type=String,Description=\"The simulated genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=sim_end,Number=1,Type=String,Description=\"The simulated genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"The event id of the introduced genomic variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"The type of strucrtural variant\">\n";
    print $output_ref2sim_inversion_vcf_fh "##INFO=<ID=END,Number=1,Type=String,Description=\"The end position of the variant described in this record\">\n";
    print $output_ref2sim_inversion_vcf_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    foreach my $ref_chr (@$refseq_arrayref) {
      if (exists $$simseq_hashref{$ref_chr}) {
        foreach my $ref_start (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} eq "INV") {
            my $ref_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'};
            my $ref_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'};
            my $sim_chr = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'};
            my $sim_start = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'};
            my $sim_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'};
            my $sim_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'};
            my $variant_type = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'};
            my $variant_id = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'};
            print $output_ref2sim_inversion_vcf_fh "$ref_chr\t$ref_start\t.\t";
            print $output_ref2sim_inversion_vcf_fh "$ref_allele\t$sim_allele\t.\t.\tvariant_type=$variant_type;variant_id=$variant_id;SVTYPE=INV;EVENT=$variant_id;END=$ref_end\n";
          }
        }
      }
    }
  }
  if ((defined $translocation_vcf) or (defined $translocation_count)) {
    print "$prefix.refseq2simseq.translocation.vcf\n\n";
    my $output_ref2sim_translocation_vcf = "$prefix.refseq2simseq.translocation.vcf";
    my $output_ref2sim_translocation_vcf_fh = write_file($output_ref2sim_translocation_vcf);
    print $output_ref2sim_translocation_vcf_fh "##fileformat=VCFv4.1\n";
    print $output_ref2sim_translocation_vcf_fh "##fileDate=$gmt_time (GMT time)\n";
    print $output_ref2sim_translocation_vcf_fh "##source=simuG.pl\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=variant_type,Number=1,Type=String,Description=\"The type of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=variant_id,Number=1,Type=String,Description=\"The id of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=ref_chr,Number=1,Type=String,Description=\"The reference genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=ref_start,Number=1,Type=String,Description=\"The reference genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=ref_end,Number=1,Type=String,Description=\"The reference genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=sim_chr,Number=1,Type=String,Description=\"The simulated genome based chromosome ID of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=sim_start,Number=1,Type=String,Description=\"The simulated genome based start coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=sim_end,Number=1,Type=String,Description=\"The simulated genome based end coordinate of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"The event id of the introduced genomic variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"The type of strucrtural variant\">\n";
    print $output_ref2sim_translocation_vcf_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    foreach my $ref_chr (@$refseq_arrayref) {
      if (not exists $$simseq_hashref{$ref_chr}) {
        foreach my $ref_start (sort {$a <=> $b} keys %{$$ref2sim_map_hashref{$ref_chr}}) {
          if ($$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'} eq "TRA") {
            my $ref_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_end'};
            my $ref_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'ref_allele'};
            my $sim_chr = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_chr'};
            my $sim_start = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_start'};
            my $sim_end = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_end'};
            my $sim_allele = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'sim_allele'};
            my $variant_type = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_type'};
            my $variant_id = $$ref2sim_map_hashref{$ref_chr}{$ref_start}{'variant_id'};
            print $output_ref2sim_translocation_vcf_fh "$ref_chr\t$ref_start\t.\t";
            print $output_ref2sim_translocation_vcf_fh "$ref_allele\t$sim_allele\t.\t.\tvariant_type=$variant_type;variant_id=$variant_id;SVTYPE=BND;EVENT=$variant_id\n";
          }
        }
      }
    }
  }
}


sub revcom {
  my $seq = shift @_;
  my $seq_revcom = reverse $seq;
  $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
  return $seq_revcom;
}



#-----------------------------------------------------------------
#----------------  Documentation / Usage / Help ------------------

=head1 NAME

simuG.pl - simulate genome sequences with designated or random genomic variants of full spectrum (SNPs, INDELs, CNVs, inversions, translocations).

=head1 SYNOPSIS

perl simuG.pl [options] [file ...]

=head1 OPTIONS

=over 8

=item B<-help> or B<-h>

Print help message. Example: -h.

=item B<-man> or B<-m>

Print more detailed help  message. Example: -m.

=item B<-version> or B<-v>

Print version information. Example: -v.

=item B<-refseq> or B<-r>

Specify the reference genome to be used as the template (in fasta or fasta.gz format). This option is mandatory. Default = "". Example: -refseq  ref.genome.fa(.gz).

=item B<-snp_vcf>

Specify the list of exact SNP variants to be introduced (in vcf or vcf.gz format). When specified, the options '-snp_count', '-snp_model', and '-titv_ratio' will be ignored. If there are also INDEL variants in the vcf file, they will be automatically ignored. Default = "". Example: -snp_vcf snp.vcf(.gz).

=item B<-snp_count>

Specify the number of SNP variants to be introduced. Default = "". Example: -snp_count 5000.

=item B<-snp_model>

Specify the SNP model file generated by the ancillary script vcf2model.pl. When specified, the option '-titv_ratio' will be ignored. Default = "". Example: -snp_model snp_model.txt.

=item B<-titv_ratio>

Specify the Ti/Tv ratio (transition/transversion ratio) used for simulate SNP variants. Default = 0.5. Example: -titv_ratio 2.0.

=item B<-indel_vcf>

Specify the list of exact INDEL variants to be introduced (in vcf or vcf.gz format). When specified, the options '-indel_count', '-indel_model', '-ins_del_ratio', '-indel_size_powerlaw_alpha', and '-indel_size_powerlaw_constant' will be ignored. If there are also SNP variants in the vcf file, they will be automatically ignored. Default = "". Example: -indel_vcf indel.vcf(.gz).

=item B<-indel_count>

Specify the number of INDEL variants to be introduced. Default = "". Example: -indel_count 500.

=item B<-indel_model>

Specify the INDEL model file generated by the ancillary script vcf2model.pl. When specified, the options '-ins_del_ratio', '-indel_size_powerlaw_alpha', and '-indel_size_powerlaw_constant' will be ignored. Default = "". Example: -indel_model indel_model.txt.

=item B<-ins_del_ratio>

Specify the Insertion/Deletion ratio used for simulate INDEL variants. Default = 1.0. Example: -ins_del_ratio 1.0.

=item B<-indel_size_powerlaw_alpha>

Specify the exponent factor alpha for power-law-fitted indel size distribution: p = C * (size) ** (alpha) for size >= 1 & size <= 50. Default = 2.0. Example: -indel_size_powerlaw_alpha 2.0.

=item B<-indel_size_powerlaw_constant>

Specify the exponent factor alpha for power-law-fitted indel size distribution: p = C * (size) ** (alpha) for size >= 1 & size <= 50. Default = 0.5. Example: -indel_size_powerlaw_constant 0.5.

=item B<-cnv_vcf>

Specify the list of exact CNV variants to be introduced (in vcf or vcf.gz format). When specified, the options '-cnv_count', '-cnv_gain_loss_ratio', '-cnv_max_copy_number', '-cnv_min_size', and '-cnv_max_size' will be ignored. Default = "". Example: -cnv_vcf cnv.vcf.

=item B<-cnv_count>

Specify the number of CNV variants to be introduced. Default = "". Example: -cnv_count 50.

=item B<-cnv_gain_loss_ratio>

Specify the relative ratio of DNA again over DNA loss. Default = 1.0. Example: -cnv_gain_loss_ratio 1.0.

=item B<-cnv_max_copy_number>

Specify the maximal copy number for CNV. Default = 10. Example: -cnv_max_copy_number 10.

=item B<-cnv_min_size>

Specify the minimal size (in basepair) for CNV variants. Default = 100. Example: -cnv_min_size 100.

=item B<-cnv_max_size>

Specify the maximal size (in basepair) for CNV variants. Default = 100000. Example: -cnv_max_size 100.

=item B<-inversion_vcf>

Specify the list of exact inversions to be introduced (in vcf or vcf.gz format). When specified, the options '-inversion_count', '-inversion_min_size', '-inversion_max_size', and '-inversion_breakpoint_gff' will be ignored. Default = "". Example: -inversion_vcf inversion.vcf(.gz).

=item B<-inversion_count>

Specify the number of inversions to be introduced. Default = "". Example: -inversion_count 5.

=item B<-inversion_min_size>

Specify the minimal size (in basepair) for inversion. Default = 1000. Example: -inversion_min_size 1000.

=item B<-inversion_max_size>

Specify the maximal size (in basepair) for inversion. Default = 1000000. Example: -inversion_max_size 1000000.

=item B<-inversion_breakpoint_gff>

Specify the list of potential breakpoints for triggering inversions (in gff3 or gff3.gz format). Default = "". Example: -inversion_breakpoint_gff inversion_breakpoint.gff(.gz).

=item B<-translocation_vcf>

Specify the list of exact translocations to be introduced (in vcf or vcf.gz format). When specified, the options '-translocation_count' and '-sv_breakpoint_gff' will be ignored. Default = "". Example: -translocation_vcf transloaction.vcf(.gz).

=item B<-translocation_count>

Specify the number of translocations to be introduced. Default = "". Example: -translocation_count 1.

=item B<-translocation_breakpoint_gff>

Specify the list of potential breakpoints for triggering translocations (in gff3 or gff3.gz format). Default = "". Example: -translocation_breakpoint_gff translocation_breakpoint.gff(.gz).

=item B<-centromere_gff>

Specify centromeres for constraining the CNV, inversion, and translocation simulation (in gff3 or gff3.gz format). Default = "". Example: -centromere_gff centromere.gff(.gz).

=item B<-excluded_chr_list>

Specify the name of chromosome(s) to be excluded for introducing genomic variants (a single-column list file in txt format). Default = "". Example: -excluded_chr_list excluded_chr_list.txt.

=item B<-seed> or B<-s>

Specify an integer as the random seed for the simulation. Default = "". Example: -seed 201812.

=item B<-prefix> or B<-p>

Specify the prefix for output files. Default = "output_prefix". Example: -prefix sim.

=back

=head1 DESCRIPTION

B<simuG.pl> can simulate genome sequences with designated or random genomic variants of full spectrum (e.g. SNPs, INDELs, CNVs, inversions, and translocations).

=head1 AUTHOR

B<Jia-Xing Yue> (GitHub ID: yjx1217)

=head1 VERSION

B<version> v1.0.0

=cut
