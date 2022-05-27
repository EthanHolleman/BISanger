rule align_to_reference:
    conda:
        '../envs/tracy.yml'
    input:
        reference=lambda wildcards: sample_df.loc[wildcards.sample]['template_path'],
        read=lambda wildcards: sample_df.loc[wildcards.sample]['trace_path']
    output:
        alignment='{output_dir}/tracyAlignments/{sample}.align.json',
        abif='{output_dir}/tracyAlignments/{sample}.align.abif'
    params:
        # tracy appends extensions so need to fake it to not make output
        # json.json
        passed_output_name='{output_dir}/tracyAlignments/{sample}.align'
    shell:'''
    tracy align -r {input.reference} -k 5 -l 20 {input.read} -o {params.passed_output_name}
    '''


rule make_TC_Table:
    conda:
        '../envs/py.yml'
    input:
        alignment='{output_dir}/tracyAlignments/{sample}.align.json'
    output:
        TC_table='{output_dir}/conversionTables/{sample}.conversion.table.tsv',
        template_table='{output_dir}/templateTables/{sample}.template.table.tsv'
    params:
        sample_name=lambda wildcards: wildcards.sample,
        reference_name=lambda wildcards: sample_df.loc[wildcards.sample]['template_name'],
        sample_treatment=lambda wildcards: sample_df.loc[wildcards.sample]['treatment'],
        bisulfite=lambda wildcards: sample_df.loc[wildcards.sample]['bisulfite'],
        topo=lambda wildcards: sample_df.loc[wildcards.sample]['state'],
        expected_read_length=lambda wildcards: sample_df.loc[wildcards.sample]['expected_read_length'],
    shell:'''
    python workflow/scripts/T2C.py --A {input.alignment} --S {params.sample_name} --R {params.reference_name} \
    --T {params.sample_treatment} --B {params.bisulfite} --O {output.TC_table} \
    --C {params.topo} --E {params.expected_read_length} --Z {output.template_table}
    '''


rule concat_TC_tables:
    conda:
        '../envs/py.yml'
    input:
        expand(
            '{output_dir}/conversionTables/{sample}.conversion.table.tsv',
            output_dir=config['output_dir'], sample=sample_df.index.values
        )
    output:
        '{output}/concatTables/concatConversionTable.tsv'
    script:'../scripts/concat.py'


rule plot_conversion_data:
    conda:
        '../envs/R.yml'
    input:
        plots=expand(
            '{output_dir}/plots/rawTraces/{sample}.rawTrace.png',
            output_dir=config['output_dir'], sample=sample_df.index.values
        ),
        concat_sample_table='{output}/concatTables/concatConversionTable.tsv'
    output:
        png='{output}/plots/{output}.conversionPlots.png'
    script:'../scripts/plotT2C.R'


rule plot_raw_traces:
    conda:
        '../envs/R.yml'
    input:
        abif='{output}/tracyAlignments/{sample}.align.abif',
        template_table='{output}/templateTables/{sample}.template.table.tsv'
    output:
        png='{output}/plots/rawTraces/{sample}.rawTrace.png'
    params:
        read_length=lambda wildcards: sample_df.loc[wildcards.sample]['expected_read_length']
    script:'../scripts/rawTrace.R'


    