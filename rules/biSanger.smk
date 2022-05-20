rule align_to_reference:
    conda:
        '../envs/tracy.yml'
    input:
        reference='',
        sample=''
    output:
        alignment=''


rule make_TC_Table:
    conda:
        '../envs/TC.yml'
    input:
        alignment=''
    output:
        TC_table=''
    parameters:
        sample_name='',
        reference_name='',
        sample_treatment=''
    shell:'''
    python scripts/T2C.py {input.alignment} {params.sample_name} {params.reference_name}
    {params.sample_treatment} {output.TC_table} 
    '''


rule concat_TC_tables:
    input:
        expand('alltables')
    output:
        'alltables.tsv'
    shell:'''
    cat {input} > {output}
    '''

    