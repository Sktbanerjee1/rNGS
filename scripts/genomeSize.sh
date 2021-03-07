#!/bin/bash
#BSUB -q debugq
#BSUB -o %genome_size.out
#BSUB -e %genome_size.err
#BSUB -n 16

helpFunction()
{
   echo ""
   echo "Usage: $0 -i input -o out_dir"
   echo -e "\t-i path to sam/bam file"
   echo -e "\t-o Path to the directory where the outputs will be written"
   exit 1 # Exit script after printing help
}

while getopts "i:o:" opt
do
   case "$opt" in
      i ) sam="$OPTARG" ;;
      o ) out_dir="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$sam" ] || [ -z "$out_dir" ]
then
   echo "All the parameters are required!";
   helpFunction
fi

# Begin script in case all parameters are correct
now=$(date)
dt=$(date +"%d%m%y")
echo "Calculating chromosome sizes"
echo "_________________________"
echo "${now}"
sam=$(realpath ${sam})
echo "sam file: ${sam}"
out_dir=${out_dir}genome_size_${dt}/
mkdir ${out_dir}
out_dir=$(realpath ${out_dir})
echo "out_dir: ${out_dir}"
echo "_________________________"
cd ${out_dir}

samtools view -H ${sam}|grep @SQ|sed 's/@SQ\tSN:\|LN://g' > ${out_dir}/ChromSize.genome 

echo "Done!"