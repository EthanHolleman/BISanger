rule align_to_reference:
    conda:
        '../envs/tracy.yml'
    input:
        reference=lambda wildcards: sample_df.loc[wildcards.sample]['template_path'],
    output:
        alignment='{output_dir}/tracyAlignments/{sample}.align.json'
    shell:'''
    '''


rule make_TC_Table:
    conda:
        '../envs/TC.yml'
    input:
        alignment='{output_dir}/tracyAlignments/{sample}.align.json'
    output:
        TC_table='{output_dir}/conversionTables/{sample}.conversion.table.tsv'
    params:
        sample_name=lambda wildcards: config[wildcards.sample],
        reference_name=lambda wildcards: sample_df.loc[wildcards.sample]['template_name'],
        sample_treatment=lambda wildcards: sample_df.loc[wildcards.sample]['treatment'],
        bisulfite=lambda wildcards: sample_df.loc[wildcards.sample]['bisulfite']
    shell:'''
    python scripts/T2C.py {input.alignment} {params.sample_name} {params.reference_name}
    {params.sample_treatment} {params.bisulfite} {output.TC_table} 
    '''


# rule concat_TC_tables:
#     input:
#         expand('alltables')
#     output:
#         'alltables.tsv'
#     shell:'''
#     cat {input} > {output}
#     '''

    