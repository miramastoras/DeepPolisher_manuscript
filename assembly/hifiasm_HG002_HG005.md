## Generating HiFiasm assemblies for training and testing DeepPolisher

The workflow for creating yak files from short reads:

https://github.com/human-pangenomics/hpp_production_workflows/blob/master/assembly/wdl/tasks/yak_no_stats.wdl

The hifiasm assembly workflow:

https://github.com/human-pangenomics/hpp_production_workflows/blob/master/assembly/wdl/workflows/trio_hifiasm_assembly_cutadapt_multistep.wdl


Short read data links for creating k-mer databases for assembly

HG003
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003_HiSeq300x_fastq/
HG004
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/

HG006
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/

HG007
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG007/


HG002 HiFi DeepConsensus 1.2
(We only used these four files to obtain a 40x subset of the reads; m64011_190714_120746.dc.q20.fastq.gz, m64011_190830_220126.dc.q20.fastq.gz, m64011_190901_095311.dc.q20.fastq.gz, m64012_190920_173625.dc.q20.fastq.gz)

https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64011_190714_120746.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64011_190728_111204.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64011_190830_220126.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64011_190901_095311.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64012_190920_173625.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64012_190921_234837.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64015_190920_185703.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/polishing/deepconsensus_runs/run_deepconsensus_hg002_20230630/m64015_190922_010918.dc.q20.fastq.gz



HG005 HiFi DeepConsensus 1.2
We did not include m64017_200730_190124.dc.q20.fastq.gz to obtain a ~40x subset of the reads.

https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200723_190224.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200730_190124.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200801_011415.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200802_073944.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200304_195708.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200309_192110.dc.q20.fastq.gz    

HG002 R941_Guppy6 (reads shorter than 100kb were filtered prior to assembly)

https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/ont/03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/ont/03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/ont/03_08_22_R941_HG002_3_Guppy_6.0.6_prom_sup.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/ont/03_08_22_R941_HG002_4_Guppy_6.0.6_prom_sup.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/ont/03_08_22_R941_HG002_5_Guppy_6.0.6_prom_sup.fastq.gz

HG005 R941_Guppy5 (reads shorter than 100kb were filtered prior to assembly)

https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG005/nanopore/Guppy_5.0.7/05_25_21_R941_GM24631_3X_Guppy_5.0.7_prom_sup.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG005/nanopore/Guppy_5.0.7/05_25_21_R941_GM24631_9X_2_Guppy_5.0.7_prom_sup.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG005/nanopore/Guppy_5.0.7/05_25_21_R941_GM24631_9X_3_Guppy_5.0.7_prom_sup.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG005/nanopore/Guppy_5.0.7/05_25_21_R941_GM24631_9X_Guppy_5.0.7_prom_sup.fastq.gz


Generating Hifiasm assemblies:

HG005
```
{
  "trioHifiasmAssembly.minOntReadLength": 100000,
  "trioHifiasmAssembly.threadCount": 64,
  "trioHifiasmAssembly.paternalYak": "/private/groups/patenlab/masri/hprc/polishing/HG005/parents/pat.HG006.yak",
  "trioHifiasmAssembly.trioHifiasm.offsetMem": [40, 20, 20],
  "trioHifiasmAssembly.maternalID": "HG007",
  "trioHifiasmAssembly.trioHifiasm.memCovRatios": [4.7, 3.8, 3.6],
  "trioHifiasmAssembly.paternalID": "HG006",
  "trioHifiasmAssembly.maternalYak": "/private/groups/patenlab/masri/hprc/polishing/HG005/parents/mat.HG007.yak",
  "trioHifiasmAssembly.filterAdapters": true,
  "trioHifiasmAssembly.childReadsHiFi": ["/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200723_190224.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200801_011415.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200802_073944.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200304_195708.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200309_192110.dc.q20.fastq.gz"],
  "trioHifiasmAssembly.trioHifiasm.dockerImage": "quay.io/masri2019/hpp_hifiasm:0.19.5",
  "trioHifiasmAssembly.homCov": 40,
  "trioHifiasmAssembly.childReadsONT": ["/private/groups/patenlab/masri/hprc/polishing/HG005/reads/R941_Guppy5/gt_100k/05_25_21_R941_GM24631_3X_Guppy_5.0.7_prom_sup.gt_100kb.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/R941_Guppy5/gt_100k/05_25_21_R941_GM24631_9X_2_Guppy_5.0.7_prom_sup.gt_100kb.fastq.gz", "/private/groups/patenlab/masri/hprc/polishing/HG005/reads/R941_Guppy5/gt_100k/05_25_21_R941_GM24631_9X_3_Guppy_5.0.7_prom_sup.gt_100kb.fastq.gz", "/private/groups/patenlab/masri/hprc/polishing/HG005/reads/R941_Guppy5/gt_100k/05_25_21_R941_GM24631_9X_Guppy_5.0.7_prom_sup.gt_100kb.fastq.gz"],
  "trioHifiasmAssembly.childID": "HG005.trio_hifiasm_0.19.5.DC_1.2_40x"
}
```

HG002
```
{
  "trioHifiasmAssembly.minOntReadLength": 100000,
  "trioHifiasmAssembly.threadCount": 64,
  "trioHifiasmAssembly.paternalYak": "/private/groups/patenlab/masri/hprc/polishing/HG002/parents/pat.HG003.yak",
  "trioHifiasmAssembly.trioHifiasm.offsetMem": [40, 20, 20],
  "trioHifiasmAssembly.maternalID": "HG004",
  "trioHifiasmAssembly.trioHifiasm.memCovRatios": [4.7, 3.8, 3.6],
  "trioHifiasmAssembly.paternalID": "HG003",
  "trioHifiasmAssembly.maternalYak": "/private/groups/patenlab/masri/hprc/polishing/HG002/parents/mat.HG004.yak",
  "trioHifiasmAssembly.filterAdapters": true,
  "trioHifiasmAssembly.childReadsHiFi": ["/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190714_120746.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190830_220126.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190901_095311.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64012_190920_173625.dc.q20.fastq.gz"],
  "trioHifiasmAssembly.trioHifiasm.dockerImage": "quay.io/masri2019/hpp_hifiasm:0.19.5",
  "trioHifiasmAssembly.homCov": 40,
  "trioHifiasmAssembly.childReadsONT": ["/private/groups/patenlab/masri/hprc/polishing/HG002/reads/R941_Guppy6/gt_100k/downsampled_40x/03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.gt_100kb.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/R941_Guppy6/gt_100k/downsampled_40x/03_08_22_R941_HG002_4_Guppy_6.0.6_prom_sup.gt_100kb.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/R941_Guppy6/gt_100k/downsampled_40x/03_08_22_R941_HG002_3_Guppy_6.0.6_prom_sup.gt_100kb.fastq.gz"],
  "trioHifiasmAssembly.childID": "HG002.trio_hifiasm_0.19.5.DC_1.2_40x"
}
```

wdl workflow used:
https://github.com/human-pangenomics/hpp_production_workflows/blob/master/assembly/wdl/workflows/trio_hifiasm_assembly_cutadapt_multistep.wdl
