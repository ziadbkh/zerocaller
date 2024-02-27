#!/usr/bin/env bash

USAGE="""Usage: $0 \\
    [-i|--input_vcf]=<input_vcf> \\
    [-c|--vep_cache]=<vep_cache> \\
    [-r|--ref_genome]=<ref_genome> \\
    [-p|--prefix]=<prefix> \\
    [-o|--options]=<options> \\
    [--vep_plugins_dir]=<vep_plugins_dir> \\
    [--plugin_cmd]=<plugin_cmd>

    Note: --vep_plugins_dir and --plugin_cmd must
    always be included together
"""

# Script that uploads the purple output to a database
if [[ $# < 3 || $# > 7 ]]
then
    echo "$USAGE"
    exit 1
fi

plugin_cmds=()
prefix="output"

for i in "$@"
do case $i in
  -i=*|--input_vcf=*)
  input_vcf=$(realpath "${i#*=}");;

  -c=*|--vep_cache=*)
  vep_cache=$(realpath "${i#*=}");;

  -r=*|--ref_genome=*)
  ref_genome=$(realpath "${i#*=}");;

  -p=*|--prefix=*)
  prefix="${i#*=}";;

  -o=*|--options=*)
  options="${i#*=}";;

  --vep_plugins_dir=*)
  vep_plugins_dir=$(realpath "${i#*=}");;

  --plugin_cmd=*)
  plugin_cmds+=("--plugin ${i#*=}");;

  *)
  echo -e "Invalid argument: $i\n$USAGE"
  exit 1;;
  esac
done

# Unzip vep_cache directory
[[ -d /opt/vep/.vep/ ]] || mkdir -p /opt/vep/.vep/
tar -xzf $vep_cache -C /opt/vep/.vep/
chmod -R 755 /opt/vep/.vep/

echo "Unzipped vep_cache into /opt/vep/.vep/"
ls /opt/vep/.vep/

species=$(ls /opt/vep/.vep/ | egrep -v '(genome|Plugins)')
build=$(ls /opt/vep/.vep/$species)
cache_version=$(echo $build | cut -d_ -f1)
cache_dir="/$species/$build/"
ls /opt/vep/.vep/${species}/

echo "VEP cached unzipped: /opt/vep/.vep/$species/$build/"

# Unzip the ref genome
# [[ -d /opt/vep/.vep/genome ]] || mkdir -p /opt/vep/.vep/{genome,Plugins}
# tar zxf $ref_genome --no-same-owner -C /opt/vep/.vep/genome/
# ref_genome_path=$(find /opt/vep/.vep/genome/ -name "*.fa")
# if [[ -z $ref_genome_path ]]
# then
#     ref_genome_path=$(find /opt/vep/.vep/genome/ -name "*.fasta")
# fi

echo "Reference genome: $ref_genome"

# Move plugins to vep_data
mkdir -p /opt/vep/.vep/Plugins
if [[ ! -z $vep_plugins_dir ]]
then
  mv $vep_plugins_dir/* /opt/vep/.vep/Plugins/
  echo "Plugins moved to /opt/vep/.vep/Plugins/"
  ls /opt/vep/.vep/Plugins/
fi

# Run VEP
echo "Running VEP"
echo """
vep \
  --cache \\
  --cache_version ${cache_version} \\
  --dir /opt/vep/.vep/ \\
  --input_file $input_vcf \\
  --output_file ${prefix}.vep.vcf \\
  --stats_file ${prefix}.vep.html \\
  --fasta ${ref_genome} \\
  --species homo_sapiens \\
  --vcf \\
  --offline \\
  --fork $(nproc) \\
  --no_progress \\
  --canonical \\
  --polyphen b \\
  --sift b \\
  --symbol \\
  --numbers \\
  --terms SO \\
  --biotype \\
  --total_length \\
  --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,INTRON,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,Feature_type,cDNA_position,CDS_position,Existing_variation,DISTANCE,STRAND,CLIN_SIG,LoFtool,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index,HGVSc,HGVSp,DOMAINS \\
  --hgvs \\
  --domains \\
  --shift_hgvs 1 ${plugin_cmds[@]} ${options}
"""
vep \
  --cache \
  --cache_version ${cache_version} \
  --dir /opt/vep/.vep/ \
  --input_file $input_vcf \
  --output_file ${prefix}.vep.vcf \
  --stats_file ${prefix}.vep.html \
  --fasta ${ref_genome} \
  --species homo_sapiens \
  --vcf \
  --offline \
  --fork $(nproc) \
  --no_progress \
  --canonical \
  --polyphen b \
  --sift b \
  --symbol \
  --numbers \
  --terms SO \
  --biotype \
  --total_length \
  --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,INTRON,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,Feature_type,cDNA_position,CDS_position,Existing_variation,DISTANCE,STRAND,CLIN_SIG,LoFtool,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index,HGVSc,HGVSp,DOMAINS \
  --hgvs \
  --domains \
  --shift_hgvs 1 \
  --plugin LoFtool ${plugin_cmds[@]} ${options} \
&& bgzip ${prefix}.vep.vcf \
&& tabix ${prefix}.vep.vcf.gz




vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html \\