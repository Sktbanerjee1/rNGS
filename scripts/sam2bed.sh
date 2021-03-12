#!/bin/bash
#BSUB -q debugq
#BSUB -o %samsort.out
#BSUB -e %samsort.err
#BSUB -n 16

helpFunction()
{
   echo ""
   echo "Usage: $0 -i target -o out_dir -c chr_sizes"
   echo -e "\t-i input sam file"
   echo -e "\t-s sample_name"
   echo -e "\t-o Path to the directory where the outputs will be written"
   echo -e "\t-c sizes of the chromosomes"
   exit 1 # Exit script after printing help
}

while getopts "i:s:o:c:" opt
do
   case "$opt" in
      i ) target="$OPTARG" ;;
      s ) sample_name="$OPTARG" ;;
      o ) out_dir="$OPTARG" ;;
      c ) chr_sizes="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$target" ] || [ -z "$sample_name" ] || [ -z "$out_dir" ] || [ -z "$chr_sizes" ]
then
   echo "All the parameters are required!";
   helpFunction
fi

# Begin script in case all parameters are correct
now=$(date)
echo "_________________________"
# print run info
echo "${now}"
target=$(realpath ${target})
echo "ref: ${target}"
mkdir -p ${out_dir}
out_dir=$(realpath ${out_dir})
echo "out_dir: ${out_dir}"
echo "sample_name: ${sample_name}"
chr_sizes=$(realpath ${chr_sizes})
echo "chr_sizes: ${chr_sizes}"
echo "_________________________"
cd ${out_dir}

seqDepthDouble=`samtools view -F 0x04 ${target} | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >${sample}.seqDepth

# sorting sam file
echo "Running picard SortSam for ${sample_name}..."
picard SortSam \
I=${target} \
O=${sample_name}_sorted.sam \
SORT_ORDER=coordinate

# markduplicates
echo "running picard mark duplicates for ${sample_name}..."
picard MarkDuplicates \
I=${sample_name}_sorted.sam \
O=${sample_name}_sorted_markdup.sam \
METRICS_FILE=${sample_name}_rmdup.txt


# convert sam to bam
echo "Converting sam to bam..."
samtools view -bS -F 0x04 ${sample_name}_sorted_markdup.sam > ${sample_name}_sorted_markdup.bam

# sort bam
echo "Sorting bam..."
sambamba sort -n ${sample_name}_sorted_markdup.bam -o ${sample_name}_sorted.bam

# convert bam to bed
echo "Creating bed file..."
bedtools bamtobed -i ${sample_name}_sorted.bam -bedpe > ${sample_name}.bed

echo "creating bedgraph"
bedtools genomecov -bg -ibam ${sample_name}_sorted.bam > ${sample_name}.bedgraph

echo "Done!"