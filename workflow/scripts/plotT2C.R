library(ggplot2)
library(RColorBrewer)
library(ggpubr)


# read concat tsv file here

read.concat <- function(table.path){

    as.data.frame(read.csv(table.path, sep='\t', header=TRUE))

}

# violin plots of replicates over conditions

violin.plots <- function(df){

    df.boiled <- subset(df, treatment=='boiled')
    print(nrow(df.boiled))
    df.cool <- subset(df, treatment=='ds')

    # Boiled samples linear vs supercoiled replicate comparison
    boiled.linear.vs.SC <- ggplot(df.boiled, aes(x=sample_name, y=T_to_C, fill=sample_name)) + 
                           geom_violin()  + facet_wrap(~topo_state) +
                           geom_jitter(shape=16, position=position_jitter(0.2)) +
                           scale_fill_brewer(palette='Dark2') + theme(legend.position='None') +
                           labs(title='Boiled')
    
    # Plot controls in same manor seperate by bisulfite treated and
    # not treated
    controls <- ggplot(df.cool, aes(x=sample_name, y=T_to_C, fill=sample_name)) +
                           geom_violin() + facet_wrap(~bisulfite) +
                           geom_jitter(shape=16, position=position_jitter(0.2)) +
                           scale_fill_brewer(palette='Dark2') + theme(legend.position='None')
    

    ggarrange(boiled.linear.vs.SC, controls, nrow=1, ncol=2)

}


line.plot <- function(df){


    df.boiled.bisulfite <- subset(df, topo_state=='supercoiled' & bisulfite=='True')
    df.boiled.linear <- subset(df, topo_state=='linear' & treatment=='boiled')
    write.csv(df.boiled.linear, 'test.linear.csv') 
    df.boiled.no <- subset(df, treatment=='ds')

    a <- ggplot(df.boiled.bisulfite, aes(x=read_index, y=T_to_C, color=sample_name)) + 
                           geom_line(alpha=0.7) +
                           scale_fill_brewer(palette='Dark2') +
                           labs(title='Supercoiled boiled bisulfite treated', x='', fill='') + theme_pubr() +
                           geom_point(pch=21)
    
    b <- ggplot(df.boiled.linear, aes(x=read_index, y=T_to_C, color=sample_name)) + 
                           geom_line() +
                           scale_fill_brewer(palette='Dark2') +
                           labs(title='Linear boiled bisulfite treated', x='', fill='') + theme_pubr() +
                           geom_point(pch=21)
    
    c <- ggplot( df.boiled.no, aes(x=read_index, y=T_to_C, color=sample_name)) + 
                           geom_line() +
                           scale_fill_brewer(palette='Dark2') +
                           labs(title='Double-stranded +/- bisulfite', fill='') + theme_pubr() + geom_point(pch=21)
    
    ggarrange(a, b, c, nrow = 3, ncol=1)


}


message(snakemake@input[['concat_sample_table']])

df <- read.concat(snakemake@input[['concat_sample_table']])
violin.plot <- violin.plots(df)
line.plot <- line.plot(df)
ggsave(snakemake@output[['png']], line.plot, dpi=300, width=10, height=6) 
