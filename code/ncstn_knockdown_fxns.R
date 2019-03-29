##########################
# some libraries we need #
##########################
library( ggplot2 )
library( here )
library( stringr )

##############################
# manual palette for figures #
##############################
cbPalette <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

########################
# relative directories #
########################
input_rel_dir <- here::here( "data" )
output_rel_dir <- here::here( "results" )
output_figure_rel_dir <- here::here( "results", "figures" )
gprof_dir <- here::here( "data", "gprofiler" )

###################
# important files #
###################
file_link_map <- list(
'growth_curve' = list( 'https://ndownloader.figshare.com/files/14072483', 'knockdown_growth_curve.csv' ),
'hek001_array_confirm' = list( 'https://ndownloader.figshare.com/files/14072486', "knockdown_hek001_rtqpcr_confirm.csv" ),
'hek293_array_confirm' = list( 'https://ndownloader.figshare.com/files/14072492', "knockdown_hek293_rtqpcr_confirm.csv" ),
'stable_knockdown_rtqpcr' = list( 'https://ndownloader.figshare.com/files/14075861', "knockdown_postpuro_rtqpcr.csv" ),
'transient_knockdown_rtqpcr' = list( 'https://ndownloader.figshare.com/files/14075870', "knockdown_prepuro_rtqpcr.csv" ),
'knockout_rtqpcr' = list( 'https://ndownloader.figshare.com/files/14075966', "knockout_rtqpcr.csv" ),
'knockdown_luciferase_assay' = list( 'https://ndownloader.figshare.com/files/14075936', "knockdown_stimulate_luciferase.csv" ),
'knockout_luciferase_assay' = list( 'https://ndownloader.figshare.com/files/14075981', "knockout_stimulate_luciferase.csv" ),
'microarray_data' = list( 'https://ndownloader.figshare.com/files/14076008', "raw_knockdown_array.csv.gz" ),
"wt_coverage" = list( "https://ndownloader.figshare.com/files/14083946", "sic75_coverage.txt" ),
"wt_flagstat" = list( "https://ndownloader.figshare.com/files/14083955", "sic75_flagstat.txt" ),
"hetdel_coverage" = list( "https://ndownloader.figshare.com/files/14083967", "sic76_coverage.txt" ),
"hetdel_flagstat" = list( "https://ndownloader.figshare.com/files/14083973", "sic76_flagstat.txt" )
)

# microarray
annotation <- here::here( "data", "annotation.csv" )
raw_knockdown_array <- here::here( "data", "raw_knockdown_array.csv" )

# RT-qPCR
knockdown_prepuro_rtqpcr <- here::here( "data", "knockdown_prepuro_rtqpcr.csv" )
knockdown_postpuro_rtqpcr <- here::here( "data", "knockdown_postpuro_rtqpcr.csv" )
knockdown_hek001_rtqpcr_confirm <- here::here( "data", "knockdown_hek001_rtqpcr_confirm.csv" )
knockdown_hek293_rtqpcr_confirm <- here::here( "data", "knockdown_hek293_rtqpcr_confirm.csv" )
knockout_rtqpcr <- here::here( "data", "knockout_rtqpcr.csv" )

# growth curve
knockdown_growth_curve <- here::here( "data", "knockdown_growth_curve.csv" )

# luciferase stimulation
knockdown_stimulate_luciferase <- here::here( "data", "knockdown_stimulate_luciferase.csv" )
knockout_stimulate_luciferase <- here::here( "data", "knockout_stimulate_luciferase.csv" )

sic75_coverage <- here::here( "data", "sic75_coverage.txt" )
sic75_flagstat <- here::here( "data", "sic75_flagstat.txt" )
sic76_coverage <- here::here( "data", "sic76_coverage.txt" )
sic76_flagstat <- here::here( "data", "sic76_flagstat.txt" )

#####################################
# some common ggplot2 modifications #
# used in most plots                #
#####################################
gg_gprofile_quadplot <- theme(
    axis.title = element_text( size=14 ),
    axis.text = element_text( size=9 ),
    plot.title = element_text( size=15 )
)

gg_bigger_texts <- theme(
    axis.title = element_text( size=22 ),
    axis.text = element_text( size=20 ),
    legend.text = element_text( size=14 ),
    legend.title = element_text( size=15 ),
    plot.title = element_text( size=22 ),
	strip.text = element_text( size=15 )
)

gg_no_legend <- theme(
	legend.position='none'
)

gg_no_grid = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

gg_no_x_grid <- theme( 
  panel.grid.major.x = element_blank() )

gg_no_y_grid <- theme( 
  panel.grid.major.y = element_blank() )

gg_center_title <- theme(
	plot.title = element_text( hjust = 0.5 )
)

# for analyzing the RT-qPCR data
######################
# calculate delta Cq #
######################
deltaCqs <- function( input_df, reference="18S", refCol="Amplicon", meanCol="mean", sdCol="sd", remove_cols = c( "n", "mean", "sd" ), include_ref=FALSE )
{
	########################################
	# split into reference and test values #
	########################################
	reference_dat <- input_df[ input_df[,refCol] == reference , ]
	
	if ( include_ref == FALSE ) {
		nonref_dat <- input_df[ input_df[,refCol] != reference , ]
	} else {
		nonref_dat <- input_df
	}
	
	return_dat <- nonref_dat[ , -c( which( colnames( nonref_dat ) %in% remove_cols ) ) ]
	
	# dCq
	return_dat$mean <- nonref_dat[,meanCol] - reference_dat[,meanCol]
	return_dat$sd <- sqrt( nonref_dat[,sdCol]^2 + reference_dat[,sdCol]^2 )
	return( return_dat )
}

# luciferase data
#####################################################
# Fxn to calculate fold-change from luciferase data #
#####################################################
# important to note this works with the log2 data
# *NOT* the regular relative intensities
luciferase_foldchange <- function( tukey_comparison_string, summary_df ) {
	string_vect <- unlist( strsplit( tukey_comparison_string, "-" ) )
	
	lhs_mean <- summary_df[ string_vect[1], "mean" ]
	lhs_sd <- summary_df[ string_vect[1], "sd" ]
	rhs_mean <- summary_df[ string_vect[2], "mean" ]
	rhs_sd <- summary_df[ string_vect[2], "sd" ]
	
	log2fc <- lhs_mean - rhs_mean
	log2sd <- sqrt( sum( lhs_sd^2, rhs_sd^2 ) )
	
	fc_mean <- 2^log2fc
	fc_sd <- log(2) * fc_mean * log2sd
	
	return( data.frame( tukey_id=tukey_comparison_string, Foldchange_mean=fc_mean, Foldchange_sd=fc_sd ) )
}
