#!/bin/bash
#BSUB -q debugq
#BSUB -o %genome_generate.out
#BSUB -e %genome_generate.err
#BSUB -n 16

helpFunction()
{
   echo ""
   echo "Usage: $0 -i ref -o out_dir -p threads"
   echo -e "\t-i path to reference genome"
   echo -e "\t-o Path to the directory where the outputs will be written"
   echo -e "\t-p threads"
   exit 1 # Exit script after printing help
}

while getopts "i:o:p:" opt
do
   case "$opt" in
      i ) ref="$OPTARG" ;;
      o ) out_dir="$OPTARG" ;;
      p ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$ref" ] || [ -z "$out_dir" ] || [ -z "$threads" ]
then
   echo "All the parameters are required!";
   helpFunction
fi

# Begin script in case all parameters are correct
now=$(date)
echo "Creating Genome Index with Bowtie2"
echo "_________________________"
echo "${now}"
ref=$(realpath ${ref})
echo "ref: ${ref}"
out_dir=$(realpath ${out_dir}genome_idx/)
echo "out_dir: ${out_dir}"
echo "threads: ${threads}"
echo "_________________________"
mkdir ${out_dir}
cd ${out_dir}
bowtie2-build --threads ${threads} -f ${ref} hg38idx
echo "Done!"