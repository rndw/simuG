# simuG

**simuG: a general-purpose genome simulator**

A simple, flexible, and powerful tool to simulate genome sequences with pre-defined or random genomic variants.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
Simulated genomes with pre-defined or random genomic variants can be very useful for benchmarking genomic and bioinformatics analyses. Here we introduce simuG as a light-weighted tool for simulating the full spectrum of genomic variants (SNPs, INDELs, CNVs, inversions, and translocations). The simplicity and versatility of simuG make it a unique general purpose genome simulator for a wide-range of simulation-based applications. 

## Citation
Jia-Xing Yue & Gianni Liti. (2018) simuG: a general-purpose genome simulator. (submitted) [bioRxiv preprint](https://www.biorxiv.org/content/early/2018/12/09/491498) 

## License
simuG is distributed under the MIT license.

## Installation
simuG is implemented in Perl5 and does not have any extra dependency. So as long as Perl5 has been installed on your system, whether it is Linux, Mac OSX or Windows, you should be able to directly run simuG via the command-line interface on your system after downloading it from GitHub:
```sh
git clone https://github.com/yjx1217/simuG.git
cd simuG
perl simuG.pl -h
perl vcf2model.pl -h
```

Please note that GNU Gzip (https://www.gnu.org/software/gzip/) needs to be pre-installed in your system if you want to run simuG.pl and vcf2model.pl with compressed input files (*.gz). 


## What’s Inside

Inside the downloaded simuG directory, you should see the following file structure:
```
.
├── LICENSE.md
├── README.md
├── Manual.pdf
├── simuG.pl
├── Testing_Example
│   ├── excluded_chr_list.yeast.txt
│   ├── sample.input.CNV.vcf.gz
│   ├── sample.input.INDEL.vcf.gz
│   ├── sample.input.inversion.vcf.gz
│   ├── sample.input.SNP.vcf.gz
│   ├── sample.input.translocation.vcf.gz
│   ├── SGDref.R64-2-1.centromere.gff3
│   ├── SGDref.R64-2-1.fa.gz
│   ├── sample.input.SNP.model.txt
│   ├── sample.input.INDEL.model.txt
│   └── Ty1_Ty3.breakpoint.gff3
└── vcf2model.pl
```

Two Perl scripts (*.pl) are provided, among which simuG.pl is the main program while vcf2model.pl is an ancillary script that can extract real-data based SNP and INDEL parameters from user-supplied SNP/INDEL variant calling VCF file. The directory Testing_Example contains some sample input files for walking through the testing examples.  


## Quick Start

Check the full list of available options.
```sh
perl simuG.pl -h
```

Simulate genome with pre-defined SNPs specified in the input VCF file.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -snp_vcf ./Testing_Example/sample.input.SNP.vcf(.gz) \ 
     -prefix output_prefix # the file name prefix for the output files
```

Simulate genome with pre-defined INDELs specified in the input VCF file.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -indel_vcf ./Testing_Example/sample.input.INDEL.vcf(.gz) \ 
     -prefix output_prefix # the file name prefix for the output files
```

Simulate genome with pre-defined CNVs specified in the input VCF file.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -cnv_vcf ./Testing_Example/sample.input.CNV.vcf(.gz) \ 
     -prefix output_prefix # the file name prefix for the output files
```

Simulate genome with pre-defined inversions specified in the input VCF file.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -inversion_vcf ./Testing_Example/sample.input.inversion.vcf(.gz) \ 
     -prefix output_prefix # the file name prefix for the output files
```

Simulate genome with pre-defined translocations specified in the input VCF file.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -translocation_vcf ./Testing_Example/sample.input.translocation.vcf(.gz) \ 
     -prefix output_prefix 
```

Simulate genome with 1000 random SNPs.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \
     -snp_count 1000 \ 
     -prefix output_prefix 
```

Simulate genome with 100 random INDELs.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -indel_count 100 \ 
     -prefix output_prefix 
```

Simulate genome with 10 random CNVs.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -cnv_count 10 \ 
     -prefix output_prefix 
```

Simulate genome with 5 random inversions.
```sh
perl simuG.pl \ 
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -inversion_count 5 \ 
     -prefix output_prefix 
```

Simulate genome with 2 random translocations.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -translocation_count 2 \ 
     -prefix output_prefix 
```

## Additional Examples

Some more advance examples are provided as follows:

Check the full list of available options of the ancillary script vcf2model.pl.
```sh
perl vcf2model.pl -h
```

Use vcf2model.pl to generate SNP and INDEL models based on the input SNP/INDEL variant calling VCF file derived from real data.
```sh
perl vcf2model.pl \
     -vcf input.real_data.SNP_INDEL.vcf(.gz) \ 
     -prefix output_prefix 
```

Use vcf2model.pl to generate SNP and INDEL models based on the input SNP/INDEL variant calling VCF file derived from real data while excluding variants called on the mitochondrial genome “chrMT” (defined in the excluded_chr_list.yeast.txt file) as well as variants with quality scores lower than 30 (if calculated).
```sh
perl vcf2model.pl \
     -vcf input.real_data.SNP_INDEL.vcf(.gz) \ 
     -qual 30 \
     -excluded_chr_list ./Testing_Example/excluded_chr_list.yeast.txt
     -prefix output_prefix 
```

Simulate genome with 1000 random SNPs and 100 random INDEL with titv_ratio = 2.0 and chrMT excluded (defined in the excluded_chr_list.yeast.txt file) for SNP simulation.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -snp_count 1000 \
     -titv_ratio 2.0 \
     -indel_count 100 \
     -excluded_chr_list ./Testing_Example/excluded_chr_list.yeast.txt \ 
     -prefix output_prefix 
```

Simulate genome with 100 random CNVs with a biased copy number gain/loss ratio of 2.0 and chrMT excluded (defined in the excluded_chr_list.yeast.txt file). Also, centromeres of the reference genome has been specified, so the simulated CNVs will not disrupt these specified centromeres.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa(.gz) \ 
     -cnv_count 100 \
     -cnv_gain_loss_ratio 2.0 \
     -centromere_gff ./Testing_Example/SGDref.R64-2-1.centromere.gff3 \
     -prefix output_prefix 
```

Simulate genome with 5 random inversions using only specified genomic features (i.e. full-length Ty1 and Ty3 transposable elements in this example) as the potential breakpoints. The mitochondrial genome chrMT (defined in the excluded_chr_list.yeast.txt file) has been excluded for this simulation. Also, centromeres of the reference genome has been specified, so the simulated inversions will not disrupt these specified centromeres.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa.gz \ 
     -inversion_count 5 \ 
     -inversion_breakpoint_gff ./Testing_Example/Ty1_Ty3.breakpoint.gff3 \
     -centromere_gff ./Testing_Example/SGDref.R64-2-1.centromere.gff3 \
     -prefix output_prefix 
```

Simulate genome with 2 random translocations using only specified genomic features (i.e. full-length Ty1 and Ty3 transposable elements in this example) as the potential breakpoints. The mitochondrial genome chrMT (defined in the excluded_chr_list.yeast.txt file) has been excluded for this simulation. Also, centromeres of the reference genome has been specified, so the simulated translocations will not disrupt these specified centromeres.
```sh
perl simuG.pl \
     -refseq ./Testing_Example/SGDref.R64-2-1.fa.gz \ 
     -translocation_count 2 \ 
     -translocation_breakpoint_gff ./Testing_Example/Ty1_Ty3.breakpoint.gff3 \
     -centromere_gff ./Testing_Example/SGDref.R64-2-1.centromere.gff3 \
     -prefix output_prefix 
```


## Full Option List
```
Options:
    -help or -h
            Print help message. 
	    Example: -h.

    -man or -m
            Print more detailed help message. 
	    Example: -m.

    -version or -v
            Print version information. 
	    Example: -v.

    -refseq or -r
            Specify the reference genome to be used as the template (in
            fasta or fasta.gz format). This option is mandatory. 
	    Default = "".
	    Example: -refseq ref.genome.fa(.gz).

    -snp_vcf
            Specify the list of exact SNP variants to be introduced (in vcf
            or vcf.gz format). When specified, the options '-snp_count',
            '-snp_model', and '-titv_ratio' will be ignored. If there are also
            INDEL variants in the vcf file, they will be automatically
            ignored. 
	    Default = "".
	    Example: -snp_vcf snp.vcf(.gz).

    -snp_count
            Specify the number of SNP variants to be introduced. 
	    Default = "". 
	    Example: -snp_count 5000.

    -snp_model
            Specify the SNP model file generated by the ancillary script
            vcf2model.pl. When specified, the option '-titv_ratio' will be
            ignored. 
	    Default = "".
	    Example: -snp_model snp_model.txt.

    -titv_ratio
            Specify the Ti/Tv ratio (transition/transversion ratio) used for
            simulate SNP variants. 
	    Default = 0.5. 
	    Example: -titv_ratio 2.0.

    -indel_vcf
            Specify the list of exact INDEL variants to be introduced (in
            vcf or vcf.gz format). When specified, the options
            '-indel_count', '-indel_model', '-ins_del_ratio',
            '-indel_size_powerlaw_alpha', and '-indel_size_powerlaw_constant'
            will be ignored. If there are also SNP variants in the vcf file,
            they will be automatically ignored.
	    Default = "". 
	    Example: -indel_vcf indel.vcf(.gz).

    -indel_count
            Specify the number of INDEL variants to be introduced. 
	    Default = "". 
	    Example: -indel_count 500.

    -indel_model
            Specify the INDEL model file generated by the ancillary script
            vcf2model.pl. When specified, the options '-ins_del_ratio',
            '-indel_size_powerlaw_alpha', and '-indel_size_powerlaw_constant'
            will be ignored. 
	    Default = "".
	    Example: -indel_model indel_model.txt.

    -ins_del_ratio
            Specify the Insertion/Deletion ratio used for simulate INDEL
            variants. 
	    Default = 1.0. 
	    Example: -ins_del_ratio 1.0.

    -indel_size_powerlaw_alpha
            Specify the exponent factor alpha for power-law-fitted indel
            size distribution: p = C * (size) ** (alpha) for size >= 1 &
            size <= 50. 
	    Default = 2.0. 
	    Example: -indel_size_powerlaw_alpha 2.0.

    -indel_size_powerlaw_constant
            Specify the exponent factor alpha for power-law-fitted indel
            size distribution: p = C * (size) ** (alpha) for size >= 1 &
            size <= 50. 
	    Default = 0.5. 
	    Example: -indel_size_powerlaw_constant 0.5.

    -cnv_vcf
            Specify the list of exact CNV variants to be introduced (in vcf
            or vcf.gz format). When specified, the options '-cnv_count',
            '-cnv_gain_loss_ratio', '-cnv_max_copy_number', '-cnv_min_size', 
	    and '-cnv_max_size' will be ignored.
	    Default = "". 
	    Example: -cnv_vcf cnv.vcf(.gz).

    -cnv_count
            Specify the number of CNV variants to be introduced. 
	    Default = "".
	    Example: -cnv_count 50.

    -cnv_gain_loss_ratio
            Specify the relative ratio of DNA again over DNA loss. 
	    Default = 1.0. 
	    Example: -cnv_gain_loss_ratio 1.0.

    -cnv_max_copy_number
            Specify the maximal copy number for CNV. 
	    Default = 10. 
	    Example: -cnv_max_copy_number 10.

    -cnv_min_size
            Specify the minimal size (in basepair) for CNV variants. 
	    Default = 100. 
	    Example: -cnv_min_size 100.

    -cnv_max_size
            Specify the maximal size (in basepair) for CNV variants. 
	    Default = 100000. 
	    Example: -cnv_max_size 100.

    -inversion_vcf
            Specify the list of exact inversions to be introduced (in vcf
            or vcf.gz format). When specified, the options
            '-inversion_count', '-inversion_min_size',
            '-inversion_max_size', '-inversion_breakpoint_gff' will be ignored.
	    Default = "".
            Example: -inversion_vcf inversion.vcf(.gz).

    -inversion_count
            Specify the number of inversions to be introduced. 
	    Default = "".
            Example: -inversion_count 5.

    -inversion_min_size
            Specify the minimal size (in basepair) for inversion. 
	    Default = 1000. 
	    Example: -inversion_min_size 1000.

    -inversion_max_size
            Specify the maximal size (in basepair) for inversion. 
	    Default = 1000000. 
	    Example: -inversion_max_size 1000000.

    -inversion_breakpoint_gff
            Specify the list of potential breakpoints for triggering
            inversions (in gff3 or gff3.gz format).
	    Default = "". 
	    Example: -inversion_breakpoint_gff inversion_breakpoint.gff(.gz).

    -translocation_vcf
            Specify the list of exact translocations to be introduced (in
            vcf or vcf.gz format). When specified, the options
            '-translocation_count' and '-translocation_breakpoint_gff' will be ignored.
            Default = "".
	    Example: -translocation_vcf transloaction.vcf(.gz).

    -translocation_count
            Specify the number of translocations to be introduced. 
	    Default = "".
	    Example: -translocation_count 1.

    -translocation_breakpoint_gff
            Specify the list of potential breakpoints for triggering
            translocations (in gff3 or gff3.gz format). 
	    Default = "".
	    Example: -translocation_breakpoint_gff translocation_breakpoint.gff(.gz).

    -centromere_gff
            Specify centromeres for constraining the CNV, inversion, and
            translocation simulation (in gff3 or gff3.gz format). 
	    Default = "".
	    Example: -centromere_gff centromere.gff(.gz).

    -excluded_chr_list
            Specify the name of chromosome(s) to be excluded for introducing
            genomic variants (a single-column list file in txt format). 
	    Default = "".
            Example: -excluded_chr_list excluded_chr_list.txt.

    -seed or -s
            Specify an integer as the random seed for the simulation.
            Default = "". 
	    Example: -seed 201812.

    -prefix or -p
            Specify the prefix for output files. 
	    Default = "output_prefix".
            Example: -prefix sim.
```

## Outputs

Upon the completion of the simulation, three files will be produced: 1) a simulated genome bearing introduced variants in FASTA format, 2) a tabular file showing the genomic locations of all introduced variants relative to both reference genome and simulated genome, 3) a VCF file showing the genomic locations of all introduced variants relative to the reference genome. 

**!!! IMPORTANT** Please note that in order to keep records on the exact genomic locations of introduced variants, simuG does not perform variant normalization for the generated VCF files. Depending on the immediate neighboring bases of the introduced genomic variants, this might have an impact if you directly compare simuG’s VCF outputs with those from other tools. Therefore, it is highly recommended to first normalize simuG’s VCF outputs as well as the VCF outputs from other tools using a VCF normalization tool such as vt (https://github.com/atks/vt) before making such comparison. You can read more about variant normalization here (https://genome.sph.umich.edu/wiki/Variant_Normalization).


## Visualization
For visualizing simulated gross chromosomal rearrangements such as inversions and translocations, we recommend using D-Genies (http://dgenies.toulouse.inra.fr/) to produce interactive dotplots between your input reference genome and simuG's simulated genome. D-Genies comes with a native online web interface (http://dgenies.toulouse.inra.fr/run) to handle user-submitted jobs, which is very convenient for the end users. Alternatively, you can also install D-Genies locally for batch job processing. Alternative tools such as Mummer (https://github.com/mummer4) and Gepard (http://cube.univie.ac.at/gepard) can also make nice dotplots for such pairwise genome comparison, although a little bit more work will be needed. 


## Contact
Please report any problems to the github issue tracker (at http://github.com/yjx1217/simuG). Alternatively, you can also write directly to me at yuejiaxing[at]gmail[dot]com.



 
