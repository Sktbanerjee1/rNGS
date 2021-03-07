#!/bin/bash
#BSUB -q debugq
#BSUB -o %trim_fastq.out
#BSUB -e %trim_fatsq.err
#BSUB -n 16

helpFunction()
{
   echo ""
   echo "Usage: $0 -i target_dir -o out_dir -s sample_name -a adapter_file -w window -m min_len"
   echo -e "\t-i path to directory with fastq files"
   echo -e "\t-o Path to the directory where the outputs will be written"
   echo -e "\t-s Name of the sample"
   echo -e "\t-a Path to the adapter fasta file"
   echo -e "\t-w Sliding Window"
   echo -e "\t-m Minimum read length"
   exit 1 # Exit script after printing help
}

while getopts "i:o:s:a:w:m:" opt
do
   case "$opt" in
      i ) target_dir="$OPTARG" ;;
      o ) out_dir="$OPTARG" ;;
      s ) sample_name="$OPTARG" ;;
      a ) adapter_file="$OPTARG" ;;
      w ) window="$OPTARG" ;;
      m ) min_len="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$target_dir" ] || [ -z "$out_dir" ] || [ -z "$sample_name" ] || [ -z "$adapter_file" ] || [ -z "$window" ] || [ -z "$min_len" ]
then
   echo "All the parameters are required!";
   helpFunction
fi

now=$(date)
dt=$(date +"%d%m%y")
echo "_________________________"
echo "Trim Illumina reads with trimmomatic..."
echo "${now}"
target_dir=$(realpath ${target_dir})
echo "target: ${target_dir}"
out_dir=${out_dir}Fastq_trimmed_${dt}/${sample_name}
mkdir -p ${out_dir}
out_dir=$(realpath ${out_dir})
echo "sample_name: ${sample_name}"
adapter_file=$(realpath ${adapter_file})
echo "adapter_file: ${adpter_file}"
echo "sliding_window: ${window}"
echo "min_length: ${min_len}"
echo "out_dir: ${out_dir}"

echo "_________________________"
# get fastq files
files=$(find ${target_dir} -name "*${sample_name}*.fastq.gz" | awk -F/ '{print $NF}')
files=( ${files} )
len_files=${#files[@]}
echo "found: ${len_files} files"


if [[ $len_files -eq 2 ]]; then
  echo "Running paired End"
  trimmomatic PE -threads 10 -phred33 -trimlog ${out_dir}/log.txt ${target_dir}/${files[0]} ${target_dir}/${files[1]} ${out_dir}/${sample_name}_R1_paired.fastq.gz ${out_dir}/${sample_name}_R1_unpaired.fastq.gz ${out_dir}/${sample_name}_R2_paired.fastq.gz ${out_dir}/${sample_name}_R2_unpaired.fastq.gz ILLUMINACLIP:${adapter_file}:2:30:10 SLIDINGWINDOW:${window} MINLEN:${min_len} > ${out_dir}/${sample_name}_trimstat.txt
elif [[ $len_files -eq 1 ]]; then
  echo "Running single End"
  trimmomatic SE -phred33 -trimlog ${out_dir}/log.txt ${target_dir}/${files[0]} ${out_dir}/${sample_name}_trimmed.fastq.gz ILLUMINACLIP:${adapter_file}:2:30:10 SLIDINGWINDOW:${window} MINLEN:${min_len} > ${out_dir}/${sample_name}_trimstat.txt
fi


echo "Trimming complete!"
echo "Done!"

