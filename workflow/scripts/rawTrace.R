library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(reshape)


read.concat <- function(table.path){

    as.data.frame(read.csv(table.path, sep='\t', header=TRUE))

}


drop.uncalled.values <- function(df, read.length){
    # Ab1 files have many points for same base
    # drop values which do not coresspond to the actual
    # base call
    subset(df, basenum!='NA' & basenum <= read.length)

}

plot.raw.trace <- function(df.base, template.table){

    df.trace <- df.base[c("peakC", "peakA", "peakT", "peakG", "basenum")]
    df.trace$T_to_C <- df.base$peak T / (df.base$peakT + df.base$peakC)
    df.trace$sum <- df.base$peakC + df.base$peakA + df.base$peakT + df.base$peakG

    df.trace$peakC <- df.trace$peakC / df.trace$sum
    df.trace$peakT <- df.trace$peakT / df.trace$sum
    df.trace$peakG <- df.trace$peakG / df.trace$sum
    df.trace$peakC <- df.trace$peakC / df.trace$sum
    

    df.melt <- melt(df.trace, 'basenum')


    all.bases <- ggplot(df.melt, aes(x=basenum, y=value, fill=variable)) +
        geom_bar(color='black', stat='identity', width=0.75) + theme_pubr() +
        scale_fill_brewer(palette='Dark2') +
        labs(title='Raw Sanger trace values by base')

    df.TC <- df.base[c("peakC", "peakT", "basenum")]
    df.melt.TC <- melt(df.TC, 'basenum')

    # scale so C positions are visible 
    max.value <- max(df.melt$value)
    template.table$value = (template.table$value * max.value) + (max.value * 0.25 * (template.table$value/template.table$value))

    TC <- ggplot() +
        geom_bar(
            data=template.table,
            aes(x=pos, y=value), fill='dodgerblue', stat='identity', alpha=0.3
            ) +
        geom_bar(
            data=df.melt,
            aes(x=basenum, y=value, fill=variable), 
            color='black', stat='identity', width=0.5
            ) +
        geom_line(
            data=df.melt,
            aes(x=basenum, y=T_to_C)
        )
        theme_pubr()

    
    ggarrange(all.bases, TC, nrow=2, ncol=1)


}


df <- read.concat(snakemake@input[['abif']])
df.template <- read.concat(snakemake@input[['template_table']])

read.length <- snakemake@params[['read_length']]
df.template <- subset(df.template, pos <= read.length)


df.drop <- drop.uncalled.values(df, snakemake@params[['read_length']])
plt <- plot.raw.trace(df.drop, df.template)
ggsave(snakemake@output[['png']], plt, dpi=300, width=35)