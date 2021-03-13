#!/bin/bash
#BSUB -q debugq
#BSUB -o %align.out
#BSUB -e %align.err
#BSUB -n 16

helpFunction()
{
   echo ""
   echo "Usage: $0 -i idx -m run_mode -f Fwd_paied -r Rev_paired -u UnmetFile -o out_dir -s sample_prefix -p threads"
   echo -e "\t-i generated genome idx"
   echo -e "\t-f paired fwd fastq/fastq.gz"
   echo -e "\t-r paired rev fastq/fastq.gz"
   echo -e "\t-u fastq/fastq.gz (for SingleEnd only)"
   echo -e "\t-o Path to the directory where the outputs will be written"
   echo -e "\t-s Sample Prefix"
   echo -e "\t-m Value must be SingleEnd/PairedEnd"
   echo -e "\t-p Threads"
   exit 1 # Exit script after printing help
}

while getopts "i:f:r:o:s:m:u:p:" opt
do
   case "$opt" in
      i ) idx="$OPTARG" ;;
      f ) fwd="$OPTARG" ;;
      r ) rev="$OPTARG" ;;
      o ) out_dir="$OPTARG" ;;
      s ) sample="$OPTARG" ;;
      m ) run_mode="$OPTARG" ;;
      u ) unmet="$OPTARG" ;;
      p ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done



# Print helpFunction in case parameters are empty
if [ -z "$idx" ] || [ -z "$out_dir" ] || [ -z "$sample" ] || [ -z "$run_mode" ] || [ -z "$threads" ]
then
   echo "ERROR! Check Parameters!";
   helpFunction
fi

# Begin script in case all parameters are correct
now=$(date)
echo "_________________________"
echo "${now}"
echo "_________________________"

# configs
out_dir=${out_dir}align/${sample}/
idx=$(realpath ${idx})
mkdir -p ${out_dir}
out_dir=$(realpath ${out_dir})
echo "_________________________"
echo "idx: ${idx}"
echo "out_dir: ${out_dir}"
echo "_________________________"

# find fastq files
if [[ $run_mode == "PairedEnd" ]]; then
  # PairedEnd
  echo "Aligning in PairedEnd Mode"
  if [ -z "$fwd" ] || [ -z "$rev" ]
   then
      echo "ERROR! Supply fwd and rev sequence.";
      helpFunction
   fi
   fwd=$(realpath ${fwd})
   rev=$(realpath ${rev})
   echo "_________________________"
   echo "fwd: ${fwd}"
   echo "rev: ${rev}"
   echo "_________________________"
   echo "running bowtie2 PairedEnd with ${threads} threads"
   cd ${out_dir}
   # run bowtie2
   bowtie2 -p ${threads}  --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
   -I 10 -X 700 \
   -x ${idx} -1 ${fwd} -2 ${rev} -S ${sample}_unsorted.sam > ${sample}_alnStat.txt
  # end Trimmed mode
elif [[ $run_mode == "SingleEnd" ]]; then
  # SingleEnd
  echo "Aligning in PairedEnd Mode"
  if [ -z "$unmet" ]
   then
      echo "ERROR! Supply Unmet Sequence!";
      helpFunction
   fi
   unmet=$(realpath ${unmet})
   echo "_________________________"
   echo "unmet: ${unmet}"
   echo "_________________________"
   echo "running bowtie2 SingleEnd with ${threads} threads"
   cd ${out_dir}
   bowtie2 -p ${threads} --very-sensitive --end-to-end --no-mixed --no-discordant --phred33 \
   -I 10 -X 700 \
   -x ${idx} -U ${unmet} -S ${sample}_unsorted.sam > ${sample}_alnStat.txt
else
  echo "Invalid Run Mode Specified"
fi

echo "Done!"