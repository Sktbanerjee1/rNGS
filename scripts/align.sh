#!/bin/bash
#BSUB -q debugq
#BSUB -o %align.out
#BSUB -e %align.err
#BSUB -n 16

helpFunction()
{
   echo ""
   echo "Usage: $0 -i idx -f ref -r rev -o out_dir -s sample"
   echo -e "\t-i generated genome idx"
   echo -e "\t-f paired fwd fastq/fastq.gz"
   echo -e "\t-r paired rev fastq/fastq.gz"
   echo -e "\t-o Path to the directory where the outputs will be written"
   echo -e "\t-s Name of the sample"
   exit 1 # Exit script after printing help
}

while getopts "i:f:r:o:s:" opt
do
   case "$opt" in
      i ) idx="$OPTARG" ;;
      f ) fwd="$OPTARG" ;;
      r ) rev="$OPTARG" ;;
      o ) out_dir="$OPTARG" ;;
      s ) sample="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$idx" ] || [ -z "$fwd" ] || [ -z "$rev" ] || [ -z "$out_dir" ] || [ -z "$sample" ]
then
   echo "All the parameters are required!";
   helpFunction
fi

# Begin script in case all parameters are correct
now=$(date)

echo "_________________________"
echo "${now}"
fwd=$(realpath ${fwd})
rev=$(realpath ${rev})
idx=$(realpath ${idx})
out_dir=${out_dir}align/${sample}/
echo "fwd: ${fwd}"
echo "rev: ${rev}"
echo "idx: ${idx}"
echo "out_dir: ${out_dir}"
echo "_________________________"
mkdir -p ${out_dir}
cd ${out_dir}

bowtie2 -p 10  --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
-I 10 -X 700 \
-x ${idx} -1 ${fwd} -2 ${rev} -S ${sample}_unsorted.sam >> ${sample}_alnStat.txt

echo "Done!"