#!/usr/bin/Rscript
# Super fast manhattan plots
# GNU General Public License v3.0
# 
# @Yanick Hagemeijer (yanickhagemeijer@gmail.com) & Niek verweij
#
# Usage: 
# e.g. Rscript ~/Analyses/generic_manhattan.R GWASdata.tsv chrom pos minus_log10_p TRUE outputfile 0.00000005 snps.to.color.tsv 1000000 black
# WHere TRUE = if Pval is logtransformed or not. 
# where snps.to.color.tsv includes position of top-snp to color:
# chromosome	position	colour
# 19	0	grey
# 
# build 37
#[^)\]}\S]\s*?(\n[ \t]*?)+print

rm(list=ls()) # reset environment

###########################################################
# library loading/installing code
###########################################################
#lib_path = "/usr/local/lib/R/site-library/"
lib_path = .libPaths()
list_of_packages = c(
		#c("parallel", "methods", "stats4"),
		#c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb"),
        c("naturalsort", "data.table" ), #"GenomicRanges"), #"data.table", "foreach", "parallel", "doParallel"), #"dtplyr"),
        #c("reshape2"),
        #c("moments", "xtermStyle", "dataview"),
        #c("bit", "bit64", "gmp", "Rmpfr"),
        c()
)

missing_packages = list_of_packages[
	! (list_of_packages %in% installed.packages(lib.loc=lib_path)[,"Package"])
]

if(length(missing_packages)) {
    print(missing_packages)
    install.packages(missing_packages, lib=lib_path, dependencies=TRUE);
}

missing_packages = list_of_packages[
	! (list_of_packages %in% installed.packages()[,"Package"])
]
if(length(missing_packages)) {
	source("http://bioconductor.org/biocLite.R")

	for( lib in missing_packages) {
		biocLite(lib, lib.loc=lib_path)
	    library(lib, lib.loc=lib_path, character.only=TRUE);
	}
}

for( lib in list_of_packages) {
    library(lib, lib.loc=lib_path, character.only=TRUE);
}

###########################################################
# http://stackoverflow.com/questions/27574775/is-it-possible-to-use-the-r-data-table-funcion-foverlaps-to-find-the-intersectio
###########################################################

intersectBedFiles.foverlaps <- function() {
	nomatch_setting=c(NA,0L)[1]
	print("finding variants within <window-size> bp of user provided known variants")
	gwas_data = foverlaps(
		gwas_data,
		known_data,
		nomatch = nomatch_setting,
		mult= c("all", "first", "last")[3]
	)
	

	print("initial clean up")
	gwas_data[,
		`:=`(
			knownSTART = NULL,
			knownSTOP = NULL
		)
	]
	setkeyv(gwas_data, c("chromosome", "position", "STOP"))


	print("looking up GWS variants and marking all variants within <window-size> bp")
	leadsnps_binned = gwas_data[PVAL <= signif_cutoff & is.na(colour), .("chromosome"=chromosome, "gwsSTART"=position-known_window, "gwsSTOP"=position+known_window) ]
	setkeyv(leadsnps_binned, c("chromosome", "gwsSTART", "gwsSTOP"))
	
	gwas_data = foverlaps(
		gwas_data,
		leadsnps_binned,
		mult = "first",
		nomatch = nomatch_setting
	)


	print("setting undefined colours")
	setkeyv(gwas_data, c("colour"))
	gwas_data[
		is.na(colour),
		colour := ifelse(
			is.na(gwsSTART),
			default_plot_color[(.GRP %% length(default_plot_color)) +1],
			signif_plot_color
		),
		by = chromosome
	]


	print("overlap check")
	overlapping = gwas_data[, .("overlapping"=STOP+1L >= shift(position, n=1, fill=+Inf, type="lead")), by=chromosome]
	overlapping[is.na(overlapping), overlapping := FALSE]
	if(any(overlapping[,overlapping], na.rm = TRUE)) {
		print("WARNING:")
	 	print("This means you have overlapping ranges!!!!!!!")
	 	overlapping = overlapping[overlapping==TRUE,]
	 	print(overlapping[, .N, by=chromosome])
	}


	print("2nd cleanup")
	gwas_data[,
		`:=`(
			gwsSTART = NULL,
			gwsSTOP  = NULL,
			STOP     = NULL
		)
	]
	setkeyv(gwas_data, c("chromosome", "position"))


	print("extra checks")
	all_cols = colnames(gwas_data)
	without_color = all_cols[ all_cols != "colour" ]
	color_distinct = uniqueN(gwas_data, by=all_cols)
	monochrome_distinct = uniqueN(gwas_data, by=without_color)
	if(color_distinct > monochrome_distinct) {
		print("This means you have OVERLAPPING REGIONS with CONFLICTING COLOURS!!!!!!!")
		print("row counts")
		print(nrow(gwas_data))
		print(color_distinct)
		print(monochrome_distinct)
	}
	

	print("duplication check")
	print(gwas_data[duplicated(gwas_data)])


	print("remove NA's & make unique")
	setkeyv(gwas_data, c("chromosome", "position"))
	gwas_data = unique(gwas_data, by=c("chromosome", "position"))
	setkeyv(gwas_data, c("chromosome", "position"))
	print(nrow(gwas_data))
	print(head(gwas_data))

	return(gwas_data)
}

###########################################################
# pop items from command line argument vector
###########################################################
pop_command_arg = function(default=NA, index=1L) {
	if(length(args) > 0L) {
		popped <<- args[ index ]
		args <<- args[ -index ]
		return( popped )
	} else {
		return( default )
	}
}

###########################################################
# calculate static info for chromosomes
###########################################################
get_chromosome_sizes = function(online_fetch=FALSE) {
	print("calculating offset for each chromosome")
	
	# thnx UCSC ^^
	chromosome_sizes_online_source = "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes"

	# hardcoded by copying online source on 27-04-2017
	chromosome_sizes_local_copy = "chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
chr5	180915260
chr6	171115067
chr7	159138663
chr8	146364022
chr9	141213431
chr10	135534747
chr11	135006516
chr12	133851895
chr13	115169878
chr14	107349540
chr15	102531392
chr16	90354753
chr17	81195210
chr18	78077248
chr19	59128983
chr20	63025520
chr21	48129895
chr22	51304566
chrX	155270560
chrY	59373566
chrM	16571
chr6_ssto_hap7	4928567
chr6_mcf_hap5	4833398
chr6_cox_hap2	4795371
chr6_mann_hap4	4683263
chr6_apd_hap1	4622290
chr6_qbl_hap6	4611984
chr6_dbb_hap3	4610396
chr17_ctg5_hap1	1680828
chr4_ctg9_hap1	590426
chr1_gl000192_random	547496
chrUn_gl000225	211173
chr4_gl000194_random	191469
chr4_gl000193_random	189789
chr9_gl000200_random	187035
chrUn_gl000222	186861
chrUn_gl000212	186858
chr7_gl000195_random	182896
chrUn_gl000223	180455
chrUn_gl000224	179693
chrUn_gl000219	179198
chr17_gl000205_random	174588
chrUn_gl000215	172545
chrUn_gl000216	172294
chrUn_gl000217	172149
chr9_gl000199_random	169874
chrUn_gl000211	166566
chrUn_gl000213	164239
chrUn_gl000220	161802
chrUn_gl000218	161147
chr19_gl000209_random	159169
chrUn_gl000221	155397
chrUn_gl000214	137718
chrUn_gl000228	129120
chrUn_gl000227	128374
chr1_gl000191_random	106433
chr19_gl000208_random	92689
chr9_gl000198_random	90085
chr17_gl000204_random	81310
chrUn_gl000233	45941
chrUn_gl000237	45867
chrUn_gl000230	43691
chrUn_gl000242	43523
chrUn_gl000243	43341
chrUn_gl000241	42152
chrUn_gl000236	41934
chrUn_gl000240	41933
chr17_gl000206_random	41001
chrUn_gl000232	40652
chrUn_gl000234	40531
chr11_gl000202_random	40103
chrUn_gl000238	39939
chrUn_gl000244	39929
chrUn_gl000248	39786
chr8_gl000196_random	38914
chrUn_gl000249	38502
chrUn_gl000246	38154
chr17_gl000203_random	37498
chr8_gl000197_random	37175
chrUn_gl000245	36651
chrUn_gl000247	36422
chr9_gl000201_random	36148
chrUn_gl000235	34474
chrUn_gl000239	33824
chr21_gl000210_random	27682
chrUn_gl000231	27386
chrUn_gl000229	19913
chrUn_gl000226	15008
chr18_gl000207_random	4262"

	chromosome_sizes_source = ifelse(
		online_fetch == TRUE,
		chromosome_sizes_online_source,
		chromosome_sizes_local_copy
	)
	
	chrom_sizes = fread(
		chromosome_sizes_source,
		sep="\t",
		col.names = c("chromosome", "size"),
		colClasses = c("character", "numeric")
	)
	#chrom_sizes[,chromosome := sub("chr", "", chromosome, ignore.case=TRUE)]
	chrom_sizes[,chromosome := naturalfactor(chromosome, ordered=TRUE ) ]

	setkeyv(chrom_sizes, c("chromosome"))
	setkeyv(gwas_data, c("chromosome"))

	# print(nrow(chrom_sizes))
	# chroms = chroms[! grepl("[_MY]", as.character(chromosome), perl=TRUE, fixed=FALSE, ignore.case = TRUE) ,]
	chrom_sizes = chrom_sizes[ chromosome %in% unique(gwas_data[["chromosome"]]), ]


	# print(nrow(chrom_sizes))
	# print(chrom_sizes)
	

	# chrom_sizes = gwas_data[,
	# 	.( size = max(position) ),
	# 	by=chromosome
	# ]
	# chrom_sizes[, chromosome := naturalfactor(chromosome) ]
	setkeyv(chrom_sizes, c("chromosome"))

	chrom_sizes[, size_cumsum := cumsum(size) ]
	setkeyv(chrom_sizes, c("chromosome"))
	chrom_sizes[, mid_chrom := size_cumsum-(size/2)]
	setkeyv(chrom_sizes, c("chromosome"))

	# print(nrow(chrom_sizes))
	# print(chrom_sizes)
	return(chrom_sizes)
}



###########################################################
print("analysing user parameters")
###########################################################
standard_plot_symbol = 015
special_plot_symbol  = 018

standard_symbol_size = 0.20
special_symbol_size  = 0.95 #1.0


known_cols = c("chromosome", "position", "colour")
known_col_types    = list(character = known_cols)
pseudo_known_file = paste(
	paste(
		known_cols,
		collapse="\t"
	),
	"\n",
	sep=""
)

default_plot_color=c("#9ecae1", "#3182bd")
signif_plot_color = "firebrick3"
# color = default_plot_color

default_signif_cutoff = 5e-8
default_known_color  = "green"
default_known_window = 1e6L #1000000L



args = commandArgs(TRUE)

input  = pop_command_arg() #1
input_chromo_col = pop_command_arg() #2
input_position_col = pop_command_arg() #3
input_pval_col = pop_command_arg() #4
input_pval_col_logtransformed = pop_command_arg() #5
output = pop_command_arg() #6

signif_cutoff_str = pop_command_arg(default=default_signif_cutoff) #6
signif_cutoff = as.numeric(signif_cutoff_str) 

output_file = paste(output, ".", signif_cutoff_str, ".tiff", sep="")

known_file   = pop_command_arg(default=pseudo_known_file) #7
known_window = as.numeric(pop_command_arg(default=default_known_window))
known_color  = pop_command_arg(default=default_known_color) #8

print("start reading known variants data")
known_data   = fread(known_file, header=TRUE, sep="\t", showProgress=TRUE) #, select=known_cols,colClasses=known_col_types
if( ! "colour" %in% colnames(known_data) ) {
	known_data[,colour := known_color]
} else {
	known_data[,colour := ifelse(is.na(colour), known_color, colour) ]
}

extra_user_cols = setdiff( colnames(known_data), known_cols )
cat("dropping unrequested columns:", paste(c("",extra_user_cols), collapse="\n\t"), "\n", sep="")
known_data[, (extra_user_cols) := NULL ]
known_data = unique(known_data, by=known_cols)
known_data = na.omit(known_data)


###########################################################
print("start reading gwas summary statistics data")
###########################################################
column_types  = list(
    character = c( input_chromo_col ),
    numeric   = c( input_position_col, input_pval_col )
)

gwas_data = fread(
	input,
	sep="\t",
	header=TRUE,
	select=c(input_chromo_col, input_position_col, input_pval_col),
	colClasses=column_types,
	showProgress=TRUE
)

setnames(gwas_data, c(input_chromo_col, input_position_col, input_pval_col), c("chromosome", "position", "PVAL"))
setkeyv(gwas_data, c("PVAL"))
gwas_data = unique(subset(gwas_data, ! is.na(PVAL) ))
if (input_pval_col_logtransformed){ 
gwas_data$PVAL <- 10^(-as.numeric(gwas_data$PVAL)) 
} # convert to normal P value if -log10 values are provided


if(nrow(gwas_data) == 0L) {
	print("you didn't provide any summarry statistics .....")
	abra = cada + bra ###this will obviously cause an error .....
}

gwas_data[, STOP := position ]
gwas_data[,
	chromosome	:= substr(
		chromosome,
		regexpr(
			"[^0]",
			chromosome
		),
		nchar(chromosome)
	)
]
gwas_data[, chromosome := paste("chr", chromosome, sep="") ]

chrom_sizes = get_chromosome_sizes(online_fetch=FALSE)

gwas_data[, chromosome := naturalfactor(chromosome, levels=levels(chrom_sizes[["chromosome"]]) ) ]
setkeyv(gwas_data, c("chromosome", "position", "STOP" ))


known_data[, position := as.numeric(position) ]
known_data[, chromosome := substr(chromosome,regexpr("[^0]",chromosome),nchar(chromosome)) ]
known_data[, chromosome := paste("chr", chromosome, sep="") ]
known_data[,
	chromosome := naturalfactor(
		chromosome,
		levels=levels(chrom_sizes[["chromosome"]])
	)
]
known_data[,
	`:=`(
		knownSTART = position-known_window,
		knownSTOP  = position+known_window
	)
]
known_data[,position := NULL]
setkeyv(known_data, c("chromosome", "knownSTART", "knownSTOP"))

gwas_data = intersectBedFiles.foverlaps() # set the colours per variant based on chromosome/proximity to <certain other> variants
setkeyv(gwas_data, c("chromosome", "position"))
print("overlaps")
print(nrow(gwas_data))
print(head(gwas_data))

gwas_data_chromos = gwas_data[,levels(chromosome)]
#print(gwas_data_chromos)
print(chrom_sizes)
reference_chromos = chrom_sizes[,levels(chromosome)]
#print(reference_chromos)
mid_chrom_offset = chrom_sizes[chromosome %in% gwas_data_chromos,mid_chrom]
print(mid_chrom_offset)

#chrom_sizes[2, "chromosome", with=TRUE]

gwas_data_chromos = subset(gwas_data_chromos, gwas_data_chromos %in% gwas_data[,unique(chromosome)])



###########################################################
print("opening file to store plot in")
print(output_file)
###########################################################
tiff(file=output_file, width=21.56, height=8, unit="cm", res=1200, pointsize=5)


###########################################################
print("creating plot")
###########################################################
setkeyv(gwas_data, c( "chromosome", "PVAL" ), physical = TRUE)
print(head(gwas_data))
setkeyv(chrom_sizes, c( "chromosome" ), physical = TRUE)

par(mar = c(7, 4, 2, 2) + 0.2)
end_point = 0.5 + length(gwas_data_chromos) + length(gwas_data_chromos)-1 #this is the line which does the trick (together with barplot "space = 1" parameter)

ignore_print_statement = gwas_data[
	1, #unique(chromosome)[1],
	plot(
		position, #+ 0, #chrom_sizes[unique(chromosome)[1],shift(size_cumsum, n=1, fill=0)],
		# position + ifelse(as.integer(gwas_data[["chromosome"]][1]) == 1, 0, sum(max_posses[1:(as.integer(gwas_data[["chromosome"]][1])-1),chromosome_max])),
		PVAL,
		pch=ifelse(PVAL <= signif_cutoff, special_plot_symbol, standard_plot_symbol),
		cex=ifelse(PVAL <= signif_cutoff, special_symbol_size, standard_symbol_size),
		col=colour,
		ylim=c(0,-log10(min(gwas_data[["PVAL"]],na.rm=TRUE))),
		xlim=c(0,chrom_sizes[.N,size_cumsum]),
		#xlim=c(0,sum(max_posses[,chromosome_max])),
		xlab="Chromosome",
		ylab="", #"-log10 of P value", #"",
		#xaxt="n",
		#yaxt="n"
		axes=FALSE
	)
]

###########################################################
print("adding layout")
###########################################################
lines(
	c(0, chrom_sizes[.N,size_cumsum]),
	# c(0, max_posses[,sum(chromosome_max)]),
	c(-log(signif_cutoff, 10),-log(signif_cutoff, 10)),
	lty="dotted",
	lwd=0.5,
	col="black"
)

axis(
	1,
	at=c(0,chrom_sizes[,size_cumsum]),
	labels = FALSE,
	# labels=rep(
	# 	"",
	# 	times=gwas_data[,uniqueN(chromosome)]+1
	# ),
	tick=TRUE,
	lwd=0.5
)

text(
	chrom_sizes[gwas_data_chromos, mid_chrom], #seq(1.5,end_point,by=2),
	par("usr")[3] -0.3,#-0.25, 
    srt = 60,
    adj= 1,
    xpd = TRUE,
    labels = gwas_data_chromos,
    cex=1.5 #0.65
)


# axis(
# 	1,
# 	at=chrom_sizes[gwas_data_chromos, mid_chrom],
# 	labels=gwas_data_chromos,
# 	cex.axis=1.5,
# 	#srt=45,
# 	tick=FALSE
# )


# #ignore_print_statement = 
# gwas_data[,
# 	axis(
# 		1,
# 		at=chrom_sizes[.GRP,mid_chrom],
# 		labels=chrom_sizes[.GRP, "chromosome", with=TRUE],
# 		cex.axis=1.5,
# 		tick=FALSE
# 	),
# 	by = chromosome
# ]


# # xEND = c()
# # for(chr in levels(chrom_sizes[,levels(chromosome)])) {
# # 	chrom_end = max_posses[chromosome==chr,chromosome_max]
# # 	xEND = c(xEND, sum(xEND)+chrom_end)
# # }
# # axis(1, at=xEND, labels=rep("", times=gwas_data[,nlevels(chromosome)]), tick=TRUE, lwd=0.5)
# axis(1, at=c(0,chrom_sizes[,size_cumsum]), labels=rep("", times=gwas_data[,uniqueN(chromosome)]+1), tick=TRUE, lwd=0.5)



min_log_of_lowest_p = gwas_data[,ceiling(-log(min(PVAL, na.rm=TRUE), base=10))] #-log(min(p$PVAL, na.rm=TRUE), base=10)
naive_stepsize = min_log_of_lowest_p / 25.0
step_size = min(max(2, naive_stepsize), 6)

range = seq(from = 0, to = min_log_of_lowest_p, by = step_size)


axis(2, at=range, labels=range, pos=c(0,0), las=1, lwd=0.5)
mtext("-log10 of P value",side=2, at=min_log_of_lowest_p/2, line=1,cex=1.5)


###########################################################
print("adding data to plot")
###########################################################

#setkeyv(gwas_data, c("chromosome"))

setkeyv(chrom_sizes, c("chromosome"))
chrom_sizes[, shifted_size_cumsum := shift(size_cumsum, n=1, fill=0)]

ignore_print_statement = gwas_data[
	,
	points(
		position + chrom_sizes[.BY, shifted_size_cumsum, on="chromosome"],
		# position + ifelse(.GRP == 1, 0, sum(max_posses[1:(.GRP-1),chromosome_max])),
		-log10(PVAL),
		pch=ifelse(PVAL <= signif_cutoff, special_plot_symbol, standard_plot_symbol),
		cex=ifelse(PVAL <= signif_cutoff, special_symbol_size, standard_symbol_size),
		col=colour,
		#xaxt="n",
		#yaxt="n"
		#axes=FALSE
	),
	by=chromosome
]
#print(warnings())




###########################################################
print("saving plot")
###########################################################
graphics.off()

#print(warnings())
