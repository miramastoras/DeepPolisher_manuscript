## Training a DeepPolisher ONT model

### ONT Q25 data, shasta v0.12.1

Phased assemblies by GFAse provided by Konstantinos:
```
/private/groups/migalab/kkyriaki/experiments/shasta_assemblies/Q25_400speed_ML05/gfase/homology_new_chainer/phase_0.fasta
/private/groups/migalab/kkyriaki/experiments/shasta_assemblies/Q25_400speed_ML05/gfase/homology_new_chainer/phase_1.fasta
/private/groups/migalab/kkyriaki/experiments/shasta_assemblies/Q25_400speed_ML05/gfase/homology_new_chainer/unphased.fasta
```

Combine assemblies to same fasta file:
```
cd /private/groups/patenlab/mira/ONT_DeepPolisher/assemblies

cat phase_0.fasta phase_1.fasta unphased.fasta > HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.fasta
```
### Prepare alignments for training

Get long read aligner wdl inputs
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

```
{
  "longReadAlignmentScattered.correctBamOptions": "--primaryOnly -m0 -a0",
  "longReadAlignmentScattered.secphaseDockerImage": "mobinasri/secphase:v0.4.3",
  "longReadAlignmentScattered.preset": "map-ont",
  "longReadAlignmentScattered.secphaseVersion": "v0.4.3",
  "longReadAlignmentScattered.alignment.dockerImage": "mobinasri/long_read_aligner:v0.4.0",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.fasta",
  "longReadAlignmentScattered.secphaseOptions": "--ont",
  "longReadAlignmentScattered.enableRunningSecphase": true,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -I8g -Y -L -p0.5",
  "longReadAlignmentScattered.kmerSize": 15,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "mm2v2.26.secphasev0.4.3",
  "longReadAlignmentScattered.readFiles": ["/private/nanopore/basecalled/q27/q25_400speed/PAW42666andPAW42495_gt10q10k_doradoTrimmed.fastq.gz"]
}
```
