#All the Directories
data_dir="/home/user/Desktop/Hemant_DH307/GSE53697/data"
fastqc_output_dir="/home/user/Desktop/Hemant_DH307/GSE53697/fastqc_output"
trimmed_output_dir="/home/user/Desktop/Hemant_DH307/GSE53697/trimmed_data"
hisat2_output_dir="/home/user/Desktop/Hemant_DH307/GSE53697/HISAT2_output"
stringtie_output_dir="/home/user/Desktop/Hemant_DH307/GSE53697/StringTie_output"
reference_annotation="/home/user/Desktop/Hemant_DH307/E-MTAB-11855/genome.gtf"

# Step 1: FASTQC
for sample_dir in "$data_dir"/*; do
    sample_name=$(basename "$sample_dir")
    
    mkdir -p "$fastqc_output_dir/$sample_name"
    
    fastqc --threads 16 "$sample_dir"/* -o "$fastqc_output_dir/$sample_name"
done

# Step 2: Trimming data using Trimmomatic
for sample_dir in "$data_dir"/*; do
    sample_name=$(basename "$sample_dir")
    
    mkdir -p "$trimmed_output_dir/$sample_name"

    for file1 in "$sample_dir"/*_1.fastq.gz; do
        file2="${file1%_1.fastq.gz}_2.fastq.gz"
        base_name=$(basename "$file1" _1.fastq.gz)

        java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 "$file1" "$file2" \
            "$trimmed_output_dir/$sample_name/${base_name}_output_1_paired.fastq.gz" "$trimmed_output_dir/$sample_name/${base_name}_output_1_unpaired.fastq.gz" \
            "$trimmed_output_dir/$sample_name/${base_name}_output_2_paired.fastq.gz" "$trimmed_output_dir/$sample_name/${base_name}_output_2_unpaired.fastq.gz" \
            ILLUMINACLIP:TruSeqPE.fa:2:30:10 LEADING:10 TRAILING:10 -phred33
        
        rm "$trimmed_output_dir/$sample_name/${base_name}_output_1_unpaired.fastq.gz"
        rm "$trimmed_output_dir/$sample_name/${base_name}_output_2_unpaired.fastq.gz"
    done
done


# Step 3: HISAT2
for sample_dir in "$trimmed_output_dir"/*; do
    sample_name=$(basename "$sample_dir")

    mkdir -p "$hisat2_output_dir/$sample_name"

    for file1 in "$sample_dir"/*_1_paired.fastq.gz; do
        file2="${file1%_1_paired.fastq.gz}_2_paired.fastq.gz"
        base_name=$(basename "$file1" _1_paired.fastq.gz)

        hisat2 -q --threads 16 -x grch38/genome -1 "$file1" -2 "$file2" -S "$hisat2_output_dir/$sample_name/${base_name}.sam"
        samtools view -@ 16 -bS "$hisat2_output_dir/$sample_name/${base_name}.sam" -o "$hisat2_output_dir/$sample_name/${base_name}.bam"
        
        rm "$trimmed_output_dir/$sample_name/${base_name}_1_paired.fastq.gz"
        rm "$trimmed_output_dir/$sample_name/${base_name}_2_paired.fastq.gz"
        rm "$hisat2_output_dir/$sample_name/${base_name}.sam"
        
        samtools sort --threads 16 "$hisat2_output_dir/$sample_name/${base_name}.bam" -o "$hisat2_output_dir/$sample_name/${base_name}.sorted.bam"
        mv "$hisat2_output_dir/$sample_name/${base_name}.sorted.bam" "$hisat2_output_dir/$sample_name/final_merged.bam"
        samtools index "$hisat2_output_dir/$sample_name/final_merged.bam"
    done
done

echo "HISAT2 finished processing all input files."

# Step 4: StringTie
for sample_dir in "$hisat2_output_dir"/*; do
    if [ -f "$sample_dir/final_merged.bam" ]; then
        sample_name=$(basename "$sample_dir")

        echo "Processing sample: $sample_name"

        stringtie "$sample_dir/final_merged.bam" -o "$stringtie_output_dir/${sample_name}.gtf" -G "$reference_annotation" -e
        
        rm "$sample_dir/final_merged.bam"

        echo "Command: stringtie "$sample_dir/final_merged.bam" -o "$stringtie_output_dir/$sample_name/${sample_name}_transcripts.gtf" -G "$reference_annotation" -e"
    else
        echo "Skipping directory $sample_dir as it does not contain the final merged BAM file."
    fi
done

