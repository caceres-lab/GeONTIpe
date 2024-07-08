import yaml

with open("config.yaml") as f:
    config = yaml.safe_load(f)

direc = config["directory"]
Ref = config["Ref"]
option = config["option"]

rule all:
    input:
        f"{direc}/Infor/Individuals.txt",
        "tracking_pipeline/directories_ready.txt",
        f"{direc}/tracking_pipeline/Download_fastq.txt",
        f"{direc}/tracking_pipeline/statistics.txt",
        f"{direc}/FileList.txt",
        f"{direc}/Infor/Reference/Std_ref.mmi",
        f"{direc}/tracking_pipeline/progress.txt",
        f"{direc}/tracking_pipeline/completed_samples.txt",
        f"{direc}/tracking_pipeline/failed_samples.txt",
        "tracking_pipeline/merging_complete.txt",
        f"{direc}/Infor/probecreation/regions_probes.txt",
        f"{direc}/tracking_pipeline/windows_ready.txt",
        f"{direc}/Infor/probecreation/allregions.txt",
        f"{direc}/Infor/probecreation/generated_probes.txt",
        f"{direc}/tracking_pipeline/created_probes.txt",
        f"{direc}/tracking_pipeline/preparing_genotyping_complete.txt",
        f"{direc}/Infor/ExpectedDistProvesInverted.txt",
        f"{direc}/tracking_pipeline/genotyping_complete.txt",
        f"{direc}/Genotypes.txt",
        f"{direc}/SV.txt",
        f"{direc}/Reads.txt"

def set1_rules():
    rule individuals:
        input:
            f"{direc}/Infor/ListaTodos.txt"
        output:
            f"{direc}/Infor/Individuals.txt"
        shell:
            "cut -f1 {input} | sort | uniq > {output}"

    rule download_fastq:
        input:
            individuals=f"{direc}/Infor/Individuals.txt",
            ready="tracking_pipeline/directories_ready.txt"
        output:
            f"{direc}/tracking_pipeline/Download_fastq.txt"
        threads: 4
        shell:
            """
            cat {input.individuals} | xargs -I {{}} -P {threads} {direc}/Scripts/Download.sh {{}}
            """

    rule statistics:
        input:
            f"{direc}/tracking_pipeline/Download_fastq.txt"
        output:
            f"{direc}/tracking_pipeline/statistics.txt"
        params:
            quality=config["Quality_reads"],
            length=config["Read_length"]
        threads: 8
        shell:
            """
            > {output}
            cat {direc}/Infor/Individuals.txt | xargs -P {threads} -I {{}} bash -c '
                NanoPlot -t {threads} --fastq {direc}/{{}}/Descargas/{{}}/*.fastq.gz --minlength {params.length} --minqual {params.quality} -o {direc}/{{}}/Resultados/Nanoplot/{{}}/sta_{{}}
                echo {{}} >> {output}
            ' bash
            """

    rule split_file:
        input:
            individuals=f"{direc}/Infor/Individuals.txt",
            download_complete=f"{direc}/tracking_pipeline/Download_fastq.txt",
            stat=f"{direc}/tracking_pipeline/statistics.txt"
        output:
            f"{direc}/FileList.txt"
        params:
            split=config["number_of_splits"]
        threads: 4
        shell:
            """
            cat {input.individuals} | xargs -I {{}} -P {threads} {direc}/Scripts/Split.sh {{}} {params.split}
            """

    rule indexing_reference:
        input:
            ref=Ref
        output:
            ref_mmi=f"{direc}/Infor/Reference/Std_ref.mmi"
        threads: 1
        shell:
            """
            minimap2 -d {output.ref_mmi} {input.ref}
            """

    rule filter_and_mapping_reads:
        input:
            filelist=f"{direc}/FileList.txt"
        output:
            f"{direc}/tracking_pipeline/progress.txt"
        params:
            dir=config["directory"],
            quality=config["Quality_reads"],
            length=config["Read_length"],
            Reads=config["Kind_reads"],
            Inversion=config["Inversion_detection"],
            SV=config["SV_detection"],
            Mapq=config["Mapq_filter"],
            ref_mmi_file=f"{direc}/Infor/Reference/Std_ref.mmi"
        resources:
            runtime=3600, #24 hrs
            mem_mb=128000
        threads: 64
        shell:
            """
            > {output}
            cat {input.filelist} | xargs -I {{}} -P {threads} {direc}/Scripts/Filt_and_mapping.sh {{}} {params.quality} {params.length} {params.Reads} {params.Inversion} {params.SV} {params.ref_mmi_file} {params.Mapq} {params.dir}
            """

    rule checking_sample_mapping:
        input:
            progress=f"{direc}/tracking_pipeline/progress.txt",
            ind=f"{direc}/Infor/Individuals.txt"
        output:
            completed=f"{direc}/tracking_pipeline/completed_samples.txt",
            failed=f"{direc}/tracking_pipeline/failed_samples.txt"
        shell:
            """
            : > {output.completed}
            : > {output.failed}

            while IFS= read -r s; do
                introduced=$(grep "${{s}}[[:space:]]" {input.progress} | cut -f3 | sort | wc -l)
                completed=$(grep "${{s}}[[:space:]]" {input.progress} | cut -f3 | sort | grep "completed" | wc -l)

                if [ "$introduced" -eq "$completed" ]; then
                    echo "${{s}}" >> {output.completed}
                else
                    echo "${{s}}" >> {output.failed}
                fi
            done < {input.ind}
            """

    rule merging_bams:
        input:
            individuals=f"{direc}/tracking_pipeline/completed_samples.txt",
        output:
            touch("tracking_pipeline/merging_complete.txt")
        threads: 4
        shell:
            """
            cat {input.individuals} | xargs -I {{}} -P {threads} bash -c '
                samtools merge {direc}/{{}}/Resultados/mappeo/{{}}.bam {direc}/{{}}/Resultados/mappeo/Sorted*.bam &&
                samtools index {direc}/{{}}/Resultados/mappeo/{{}}.bam &&
                rm {direc}/{{}}/Resultados/mappeo/Sorted*.bam &&

                samtools merge {direc}/{{}}/Resultados/mappeo/{{}}_Finalnomapped.bam {direc}/{{}}/Resultados/mappeo/Unmappedsorted*.bam &&
                samtools index {direc}/{{}}/Resultados/mappeo/{{}}_Finalnomapped.bam &&
                rm {direc}/{{}}/Resultados/mappeo/Unmappedsorted*.bam
            '
            """
    
def set2_rules():
    rule checking_sample_mapping:
        input:
            ind=f"{direc}/Infor/Individuals.txt"
        output:
            touch("tracking_pipeline/merging_complete.txt"),
            f"{direc}/tracking_pipeline/Download_fastq.txt",
            f"{direc}/tracking_pipeline/statistics.txt",
            f"{direc}/FileList.txt",
            f"{direc}/tracking_pipeline/progress.txt",
            f"{direc}/tracking_pipeline/failed_samples.txt",
            completed=f"{direc}/tracking_pipeline/completed_samples.txt"
        shell:
            """
            : > {output.completed}

            while IFS= read -r s; do
                echo "${{s}}" >> {output.completed}
            done < {input.ind}

            # Crear todos los archivos de salida como archivos vacíos
            touch {output[0]}
            touch {output[1]}
            touch {output[2]}
            touch {output[3]}
            touch {output[4]}
            touch {output[5]}
            touch {output[6]}
            """

    rule created_probes_regions:
        input:
            cords=f"{direc}/Infor/allcoords.txt"
        output:
            f"{direc}/Infor/probecreation/regions_probes.txt"
        params:
            flanking_region=config["flanking_region"],
            routes=config["extra_route"]
        shell:
            """
            {direc}/Scripts/step1_probes.sh {params.flanking_region} {params.routes}
            """

    rule creating_probes_creating_windows:
        input:
            regions=f"{direc}/Infor/probecreation/regions_probes.txt"
        output:
            f"{direc}/tracking_pipeline/windows_ready.txt"
        threads:8
        params:
            probe_size=config["size_probes"],
            overlap=config["overlap_probes"],
            routes=config["extra_route"],
            minimum_unique_sequence=config["min_uni_seq"]
        shell:
            """
            cat {input.regions} | xargs -I {{}} -P {threads} {direc}/Scripts/generating_windows.sh {{}} {params.probe_size} {params.overlap} {params.routes} {params.minimum_unique_sequence}
            """

    rule creating_probes_sorting_windows:
        input:
            f"{direc}/tracking_pipeline/windows_ready.txt",
            cords=f"{direc}/Infor/allcoords.txt"
        output:
            f"{direc}/Infor/probecreation/allregions.txt"
        shell:
            """
            {direc}/Scripts/sorting_windows.sh
            """

    rule creating_probes:
        input:
            f"{direc}/Infor/probecreation/allregions.txt"
        output:
            f"{direc}/Infor/probecreation/generated_probes.txt"
        threads:8
        params:
            extra_region=config["extra_region"],
            Main_reference=config["Ref"],
            Secundary_reference=config["seq_ref"],
            routes=config["extra_route"]
        shell:
            """
            makeblastdb -in {params.Main_reference} -dbtype nucl
            makeblastdb -in {params.Secundary_reference} -dbtype nucl
            cat {input} | xargs -I {{}} -P {threads} {direc}/Scripts/generating_probes.sh {{}} {params.extra_region} {params.Main_reference} {params.Secundary_reference} {params.routes}
            """

def set3_rules():
    rule checking_sample_and_probes:
        input:
            ind=f"{direc}/Infor/Individuals.txt"
        output:
            touch("tracking_pipeline/merging_complete.txt"),
            f"{direc}/tracking_pipeline/Download_fastq.txt",
            f"{direc}/tracking_pipeline/statistics.txt",
            f"{direc}/FileList.txt",
            f"{direc}/tracking_pipeline/progress.txt",
            f"{direc}/tracking_pipeline/failed_samples.txt",
            f"{direc}/Infor/probecreation/regions_probes.txt",
            f"{direc}/tracking_pipeline/windows_ready.txt",
            f"{direc}/Infor/probecreation/allregions.txt",
            f"{direc}/Infor/probecreation/generated_probes.txt",
            completed=f"{direc}/tracking_pipeline/completed_samples.txt"
        shell:
            """
            : > {output.completed}

            while IFS= read -r s; do
                echo "${{s}}" >> {output.completed}
            done < {input.ind}

            # Crear todos los archivos de salida como archivos vacíos
            touch {output[0]}
            touch {output[1]}
            touch {output[2]}
            touch {output[3]}
            touch {output[4]}
            touch {output[5]}
            touch {output[6]}
            touch {output[7]}
            touch {output[8]}
            touch {output[9]}
            touch {output[10]}
            """

if option == "complete":
    set1_rules()
elif option == "genotype":
    set2_rules()
elif option == "genotype_wt_probes":
    set3_rules()
else:
    raise ValueError(f"Invalid option: {option}")

rule preparing_directories:
    input:
        ind=f"{direc}/Infor/Individuals.txt"
    output:
        touch("tracking_pipeline/directories_ready.txt")
    shell:
        """
        while IFS= read -r i; do
            mkdir -p {direc}/$i/Descargas/$i
            mkdir -p {direc}/$i/Resultados/Nanoplot
            mkdir -p {direc}/$i/Resultados/Filt_fq
        done < {input.ind}
        """

rule preparing_genotyping_inputs:
    input:
        f"{direc}/Infor/coords.txt"
    output:
        allcoords=f"{direc}/Infor/allcoords.txt",
        listaref=f"{direc}/Infor/ListaRef.txt"
    params:
        flanking_region=config["flanking_region"]
    shell:
        """
        > {output.allcoords}
        > {output.listaref}
        for i in $(cut -f1 {input}); do
            region="$(grep ${{i}} {input})"
            bp1_1="$(echo ${{region}} | sed 's/ /\t/g' | cut -f3)"
            bp1_2="$(echo ${{region}} | sed 's/ /\t/g' | cut -f4)"
            bp2_1="$(echo ${{region}} | sed 's/ /\t/g' | cut -f5)"
            bp2_2="$(echo ${{region}} | sed 's/ /\t/g' | cut -f6)"
            sizebp1="$((bp1_2 - bp1_1 + 1))"
            sizebp2="$((bp2_2 - bp2_1 + 1))"
            sizeinv="$((bp2_1 - bp1_2 + 1))"
            echo -e "${{region}}\t${{sizeinv}}\t${{sizebp1}}\t${{sizebp2}}" >> {output.allcoords}

            ID="$(echo ${{region}} | sed 's/ /\t/g' | cut -f1)"
            chr="$(echo ${{region}} | sed 's/ /\t/g' | cut -f2)"
            start="$((bp1_1 - {params.flanking_region}))"
            end="$((bp2_2 + {params.flanking_region}))"
            echo -e "${{ID}}\tchr${{chr}}:${{start}}-${{end}}" >> {output.listaref}
        done
        """

rule final_probes:
    input:
        f"{direc}/Infor/allcoords.txt",
        f"{direc}/Infor/probecreation/generated_probes.txt"
    output:
        f"{direc}/tracking_pipeline/created_probes.txt"
    shell:
        """
        {direc}/Scripts/final_steps.sh
        """

rule Preparing_genotyping:
    input:
        individuals=f"{direc}/tracking_pipeline/completed_samples.txt",
        merging_complete="tracking_pipeline/merging_complete.txt",
        created_probes=f"{direc}/tracking_pipeline/created_probes.txt"
    output:
        f"{direc}/tracking_pipeline/preparing_genotyping_complete.txt"
    threads: 4
    shell:
        """
        > {direc}/tracking_pipeline/preparing_genotyping_complete.txt
        cat {input.individuals} | xargs -I {{}} -P {threads} bash -c '
            cp -r {direc}/Genotyping {direc}/{{}} &&
            cp -r {direc}/Scripts/Genotyping_scripts/* {direc}/{{}}/Genotyping'
            echo {{}} >> {direc}/tracking_pipeline/preparing_genotyping_complete.txt
        """

rule expected_distances:
    input:
        f"{direc}/tracking_pipeline/created_probes.txt",
        merging_complete="tracking_pipeline/merging_complete.txt"
    output:
        f"{direc}/Infor/ExpectedDistProvesInverted.txt"
    shell:
        """
        {direc}/Scripts/invdist.sh
        """

rule Genotyping:
    input:
        merging_complete="tracking_pipeline/merging_complete.txt",
        preparing_genotyping_complete=f"{direc}/tracking_pipeline/preparing_genotyping_complete.txt",
        expected_distances=f"{direc}/Infor/ExpectedDistProvesInverted.txt",
        individuals=f"{direc}/tracking_pipeline/completed_samples.txt"
    output:
        f"{direc}/tracking_pipeline/genotyping_complete.txt"
    resources:
        runtime=3600, #24 hrs
        mem_mb=16000
    threads: 10
    params:
        local_analysis=config["local"],
        ftp_location=config["ftp_route"],
        task=config["task"],
        cov_probes=config["coverage_of_probes"],
        ID_prob=config["identity_probe"],
        Minor_prop=config["Minor_allele_proportion"],
        Low_conf_prop=config["Low_confidence_proportion"],
        Per_size_diff=config["Difference_respect_expected"],
        Min_SV=config["Minimum_size_SVs"],
        Split_SV=config["Ratio_to_split_SV"],
        SNP_file=config["SNP_selection_route"],
        Ref_file=config["Ref"],
        Q_base=config["Quality_by_base"],
        Min_freq=config["Min_freq_snps"],
        Min_cov=config["Min_coverage"],
        Dist_tree=config["Distance_of_tree"]
    shell:
        """
        cat {input.individuals} | xargs -I {{}} -P {threads} bash -c '{direc}/{{}}/Genotyping/runall.sh {direc}/{{}}/Genotyping {{}} {params.local_analysis} {params.ftp_location} {params.task} {params.cov_probes} {params.ID_prob} {params.Minor_prop} {params.Low_conf_prop} {params.Per_size_diff} {params.Min_SV} {params.Split_SV} {params.SNP_file} {params.Ref_file} {params.Q_base} {params.Min_freq} {params.Min_cov} {params.Dist_tree}
        touch {output}'
        """

rule Genotypes_and_SVs:
    input:
        individuals=f"{direc}/tracking_pipeline/completed_samples.txt",
        genotypes=f"{direc}/tracking_pipeline/genotyping_complete.txt"
    output:
        f"{direc}/Genotypes.txt",
        f"{direc}/SV.txt",
        f"{direc}/Reads.txt"
    shell:
        """
        {direc}/Scripts/FinalGenotype.sh
        """