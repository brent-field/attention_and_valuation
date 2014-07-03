#!/usr/bin/perl

### Under development (peretually) ###
# to do:
# * add head orientation (for subjects who are rotated)
# * add 3dTqual
# * create Package with subroutines, make parser a sub-routine    http://mathforum.org/~ken/perl_modules.html
# * experiment with stim-maxlag or 1 or 2 in regression (not deconv) http://afni.nimh.nih.gov/afni/afniboard/OLDmessages/372.html
#
######################################

use strict;
use warnings;
use diagnostics;

# invoke Modules/Packages
use Cwd;
use File::Copy;
use Getopt::Std;

### DEFINES ###
use constant OFF => 0;
use constant ON => 1;
use constant FALSE => 0;
use constant TRUE => 1;
use constant LOG => 0;
use constant MANUAL => 1;
use constant AUTO => 2;
use constant WARNING => 1;
use constant ERROR => 2;
use constant ANATNUDGE => 1;
use constant ALLINEATE => 2;

use constant DELETE_PREVIOUS_FILES => ON;
use constant PARRALLELIZE => ON;
use constant DEBUG_ONLY => FALSE;


        ### Preprocessing variables ###
use constant BASELINE_POLYNOMIAL_MODE => MANUAL;  # AUTO or MANUAL; AUTO recommended. MANUAL means that the baseline value on *info.txt is used.
use constant WITHIN_SUBJECT_SEM => OFF;         # calculate within-subject SEM (BETWEEN will still be calculated)
use constant USE_PREEXISTING_DESPIKED_DATA_IN_REALIGN => FALSE;
use constant USE_PREEXISTING_TALAIRACHING_MODE => OFF;  # OFF, MANUAL, AUTO
use constant NUM_MOTION_PARAMETERS => 6;
use constant JUST_REGRESS_OUT_BASELINE => OFF;    # run 3dDeconvolve without regressors of interest, so as to see effects of baseline model in residual
use constant TR_FOR_REALIGNMENT => 25;

use constant TEMP => OFF;



#****** Program Steps ********
### Preprocessing Steps ###
use constant DICOM_TO_BRIK => OFF;    #Create AFNI files    was set to on 7-16-2013
use constant DESPIKE_PRE_VOLREG => OFF;
use constant REALIGN => OFF;
use constant DESPIKE_POST_VOLREG => OFF;  # recommended in most cases
use constant SPATIAL_SMOOTHING => OFF;
use constant BINARY_MASK => OFF;
use constant CALC_SIGNAL_CHANGE => OFF;
use constant MPRAGE_EPI_ALLIGNMENT => OFF;  # optiones are OFF, ANATNUDGE, ALLINEATE
use constant ADD_EDGE => OFF;                                # to check lpc_align.py quality. Leave off, unless specifically checking
use constant SKULL_STRIP => OFF;          #if you already have skull-stripped t1, much faster if this is off!!!
use constant TALAIRACH_MODE_EARLY => OFF;   # possible values are OFF, MANUAL, AUTO
use constant WAVER => OFF;
use constant GENERAL_LINEAR_TEST => OFF;
use constant ESTIMATE_IMPULSE_RESPONSE => OFF;

use constant TALAIRACH_MODE_LATE => OFF;   # usually off - possible values are OFF, MANUAL, AUTO
use constant STICK_FUNCTION_GLT => OFF;                                                                ### these are special preprocessing steps - generally not used
use constant GENERAL_LINEAR_TEST_BASELINEMODEL_ONLY => OFF;
use constant GENERAL_LINEAR_TEST_3BackBlock_only => OFF;
use constant GENERAL_LINEAR_TEST_3BackBlock_only_no_baseline => ON;


        ### ROI Preprocessing steps ###
use constant USE_POST_VR_DESPIKE_IN_ROI_ANALYSIS => ON;
use constant BINARY_MASK_ROI => OFF;
use constant CALC_SIGNAL_CHANGE_ROI => OFF;
use constant GENERAL_LINEAR_TEST_ROI => OFF;
use constant ESTIMATE_IMPULSE_RESPONSE_ROI => OFF;

        ### Group Anslysis ###
use constant GROUP_REGRESSION_STATS => OFF;
use constant GROUP_DECONV_STICK_STATS => OFF;
use constant GROUP_IRESP_ANALYSIS => OFF;
use constant ThreeD_ANOVA_AE => OFF;
use constant ThreeD_ANOVA_AE_PLOS => ON;
use constant ThreeD_ANOVA_AE_STICK => OFF;
use constant ThreeD_ANOVA_AE_PLOS_STICK => OFF;
use constant AE_PROBE_CONTRIB_TO_SPRITZ_IRESP => OFF;
use constant ThreeD_ANOVA_MED => OFF;
use constant ThreeD_ANOVA_ERIKSEN => OFF;
use constant ThreeD_ANOVA_ERIKSEN_SEQADJ => OFF;
use constant ThreeD_ANOVA_NSTROOP => OFF;
use constant NSTROOP_GROUP => OFF;

        ### ROI Group Anslysis ###
use constant GROUP_ROI_CLUSTER_STATS => OFF;
use constant GROUP_ROI_IRESP_ANALYSIS => OFF;

        ### cluster ANALYSIS  ####
use constant WRITE_INDIVIDUAL_MASKFILES => OFF;   ### if these files already exist, turning this off can speed processing
use constant USE_ALPHASIM_FAMILYWISE => OFF; ### MUST BE HANDCODED FOR EACH EXPERIMENT !!! ####
use constant WRITE_STAT_FILE => OFF; # WRITE the cluster by subject stat file ; slows processing some;


    ### define GLOBAL variables ###
use vars qw(@subjs $afni_root_path $clusterdir_name $scripts_dir @num_in_group @censor_files);
use vars qw($subj $groupdir_name @regressor_name $max_nudge @EPI_count $regressor_num @subject_group $group_num @unique_groups @contrast_name $talrch_type);
use vars qw($TRs_pre_regressor_for_estimating_baseline_time_courses $regressor_indicator_char);
our $studyname = "t";
our $log_output_buffer = "";                       # our = globally defined
our $log_file_open_flag = FALSE;
$regressor_name[0] = "";
$subj = 0;

    ### define LOCAL variables  ###
my $data_read = OFF;
my $num_subjs = 0;
my $regressor_count = 0;
my $regressor_read = OFF;
my @num_scans; #[0][0] = 0;
my $TR_duration = 2;

my $censor_fileroot = "_censor.1d";
my ($cmd, $char, $cwd, $ROI_num);
my ($contrast, $deconv_bucket_name, $p, @pfxs, $pfx, $pre_vr_despike_flag, $post_vr_despike_flag, $input, $subject_return_value, @runs, $subject);
my ($delete_intermediate_files, @preprocess_flag, @debridge_mask_filename, $r, $stder_betweensubj_cmd, $stder_withinsubj_cmd);
my ($mean_cmd, $sem_cmd, $iresp_ave_cmd, $stder_betweensubj_script_filename, $ttest_cmd, $test_num, $groupstat_script_filename, $line );
my ($reg_count, $output_coorindate_system, $deconv_input_file, $cat_cmd, $i, $j, @volregfiles, $spatial_blur_diameter, $run_num, $allnorm);
my ($contast_num, $regressor, @baseline_model_order, $HRF_regressor, $subj_dir, $dicom_dir, $logdir_name, $oldlogs_dirname);
my ($program_filename, @info_files, $contrast_dirname, $totalTRs, $catcmd, $concat_string_ideal, $concat_string_timings );
my ($concat_string_timings_shifted_for_iresp, $thresh_num, $clust_num, $cmd_line_args, $print_help_text);
my ($program_directory, $contrast_num, @contrast, @regressor_root_fname, @thresh, @debridge_mask_cm_filename, $bridge_mask_num);
my (@MPRAGE_series, @EPI_series, $dicom_root_path, @dicom_directory, $HRF_duration_sec, $num_contrasts, %group_hash, $num_in_group_counter, $subject_num);
my ($cluster_source_num, @cluster_source_subbrik, $subbrick);
my (%cmd_line_args, $clust_flag_num, @cluster_source_name, @min_cluster_size, $min_cluster_index);
my (@clust_source_flag, $clust_regressor_count, @clust_regressor_flag, $HRF_duration_TR, $TR_duration_ms, $num_slices_z_axis, @despike, $input_reg_file, $output_reg_file);
my ($cmd2, $cmd3, $line_num, $deconv_mask_file, $CensoredVoxelsFileString, $reg_file, $num_input_timeseries, $shifted_for_iresp_reg_file );
my ($INPUT_REG_FILE, @impulse_count, @reg_name_deconved, @test_names, $num_tests, @test_betacoef_subbrik, $first_in_group);
my ($iresp_ave_script_filename, $stder_withinsubj_script_filename, $num_of_ROIs, $iresp_ave_cmd_fim, $concat_iresp_ave_cmd, @nudgefiles, $nudgefile);
my ($clean_censor_file, $quality_script_filename, $pre_reg_TR );
my ($Blevel_subject, $Clevel_subject );    #ANOVA
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst, $starttime_string); # time of launch
my ($spatial_smooth_flag);

my $log_filename;

# Create and initialize flags for main() program blocks;
my ($preprocess_group, $preprocess_subject, $group_analysis, $roi_analysis_group, $roi_analysis_subject, $ROI_group_analysis, $cluster_analysis_all);
my ($cluster_analysis_specific,$quality_analysis,$plot_quality_analysis);
$preprocess_group = FALSE;
$preprocess_subject = FALSE;
$group_analysis = FALSE;
$roi_analysis_group = FALSE;
$roi_analysis_subject = FALSE;
$ROI_group_analysis = FALSE;
$cluster_analysis_all = FALSE;
$cluster_analysis_specific = FALSE;
$quality_analysis = FALSE;
$plot_quality_analysis = FALSE;
$line = "";
$test_betacoef_subbrik[0][0] = 0;
$test_names[0][0] = 0;
$subject_return_value = 0;
$TRs_pre_regressor_for_estimating_baseline_time_courses = 2;
$regressor_indicator_char = "@";
$spatial_smooth_flag = "";

# Get time of launch
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$starttime_string = sprintf "%4d-%02d-%02d_%02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec;

#### PARSE command line
%cmd_line_args = ();    # Hash keys and values of command line arguments
getopts("pi:gdRr:OCcQL", \%cmd_line_args);
#$subj = $cmd_line_args{i};
if (defined $cmd_line_args{r}) {$subj = $cmd_line_args{r}; }
$print_help_text = TRUE;
if (defined $cmd_line_args{p}) {$preprocess_group = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{i}) {$subj = $cmd_line_args{i}; $preprocess_subject = TRUE; $print_help_text = FALSE; }
if (defined $cmd_line_args{g}) {
        if (defined $cmd_line_args{p} or defined $cmd_line_args{i}) {
                print "Group processing must be done seperately in preprocessing at this time. \n" ;
                exit;
        } else {
                $group_analysis = TRUE; $print_help_text = FALSE;
        }
}
if (defined $cmd_line_args{d}) {$delete_intermediate_files = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{R}) {$roi_analysis_group = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{r}) {$roi_analysis_subject = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{O}) {$ROI_group_analysis = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{c}) {$cluster_analysis_specific = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{C}) {$cluster_analysis_all = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{Q}) {$quality_analysis = TRUE; $print_help_text = FALSE;}
if (defined $cmd_line_args{L}) {$plot_quality_analysis = TRUE; $print_help_text = FALSE;}

if (keys(%cmd_line_args) == 0) {
        print "$0: fMRI Analysis script for AFNI\n";
        print "Assumes that there is a *_info.txt in same directory.\n";
        print "<optional name>info.txt should be a Windows Text (space) delimited Excel spreadsheet.\n";
        print "Also assumes that regressor (1D) files are in each subjects directory.\n";
        print "By Brent Field. Not waranteed to be bug free.\n";
        print "Version 1.3\n";
        print "Flags \n";
        print "-p Preprocess a group of subjects (potentaily parrallelized) \n";
        print "-i <subject #> Preprocess subject. subject # is nth subject listed in *info.txt, starting with 0 \n";
        print "-g Run group analysis \n";
        print "-d Delete intermediate preprocessing files (to save space)\n";
        print "-R Region of interest analysis preprocessing for a group of subjects \n";
        print "-r <subject #> Region of interest analysis preprocessing for subject. subject # is nth subject listed in *info.txt, starting with 0 \n";
        print "-O Region of interest group analysis\n";
        print "-C Cluster analysis: Analyze all clusters, as specified in info.txt\n";
        print "-c <source> <threshold> <clust_num>: analyze a specific cluster\n";
        print "-Q Analyze raw data for quality problems and write out quality files\n";
        print "-L Plot data from the -Q quality analysis\n";

        exit;
}

$scripts_dir = cwd() ."/";

$_ = $0;
/(^.*\W)(.+\.pl)$/;                                         # Regular expression that seperately parses directory and filename
$program_directory = $1;
$program_filename = $2;


if (USE_PREEXISTING_TALAIRACHING_MODE == MANUAL) {$talrch_type = "mt" ;} else {$talrch_type = "at" ;}

parse_info_file();


if(chdir($afni_root_path) == 0) {
        if(mkdir($afni_root_path) == 0) {
                log_file("Bad file path specificied for afni root folder: $afni_root_path",ERROR); }
}

$groupdir_name = $afni_root_path."GroupAnalysis_PLoS/";  # change back to GroupAnalysis   7-16-2013
if(chdir($groupdir_name) == 0) {
        if(mkdir($groupdir_name) == 0) { log_file("Bad file path specificied for log file directory: $groupdir_name",ERROR); }
}

$clusterdir_name = $groupdir_name."Clusters/";
if(chdir($clusterdir_name) == 0) {
        if(mkdir($clusterdir_name) == 0) { log_file("Bad file path specificied for log file directory: $clusterdir_name",ERROR); }
}

$logdir_name = $afni_root_path."AnalysisLogs/";
if(chdir($logdir_name) == 0) {
        if(mkdir($logdir_name) == 0) { log_file("Bad file path specificied for log file directory: $logdir_name",ERROR); }
}
$oldlogs_dirname = $logdir_name."old/";
if(chdir($oldlogs_dirname) == 0) {
        if(mkdir($oldlogs_dirname) == 0) { log_file("Bad path specificied for old log files: $oldlogs_dirname",ERROR); }
} else {
                $cmd = "mv -u $logdir_name*??:??:?? $oldlogs_dirname";
                system($cmd);
}

print $logdir_name;
print $program_filename;
print $subjs[$subj];
print $starttime_string;

$log_filename = $logdir_name.$program_filename."_subj_$subjs[$subj].log.$starttime_string";
if(open(LOGFILE, ">".$log_filename)) {
        $log_file_open_flag = TRUE;
#        if ($preprocess_group == TRUE) {
                print "$logdir_name \n";
                #copy("$0","$logdir_name$program_filename.$starttime_string") or log_file("Can't copy $program_filename to log dir:  $log_filename",ERROR);

                #copy("$scripts_dir$info_files[0]","$logdir_name$info_files[0].$starttime_string") or log_file("Can't copy $info_files[0] to log file directory",ERROR);
#        }
} else { log_file("Can't open log file $log_filename: $!", ERROR) }



                                        ############## Preprocess Group #############

if ($preprocess_group == TRUE) {
        chdir($scripts_dir);
        for ($subj=0; $subj< $num_subjs; $subj++) {                                        # Main loop (loops on subject)
                if ($preprocess_flag[$subj] ne "Y") { next; }                                                                        # skip to next subject if processing isnt explicitly turned on

                log_file("Starting $0 to process participant $subjs[$subj]",LOG);
                $cmd = $program_filename." -i ".$subj;
                if (PARRALLELIZE) {
                        $cmd = "submit " . $cmd;
                        print "$cmd\n";
                        log_file($cmd,LOG);
                        system($cmd);
                } else {
                        $subject_return_value = execute_afni_command($cmd);
                }
                if ($subject_return_value == 1) { log_file("Problems processing subject $subjs[$subj]\n",WARNING); }
        }                         ##End subject loop - main loop ##
}   ### End preprocess Group


                                        ############## Preprocess / Process Subject #############

if ($preprocess_subject == TRUE) {

                # exceptions - remove regressors within particular subjects that have colinearity, Must be at end of regressor list in info.txt because this is eliminated
                if ($subjs[$subj] eq "med07" && $studyname eq "Nstroop_med") {
                        $regressor_count = $regressor_count - 1;
                }
                if ($subjs[$subj] eq "med07" && $studyname eq "Nstroop_med_SeqAdj") {
                        $regressor_count = $regressor_count - 1;
                }


        log_file("Getting events for participant $subjs[$subj]",LOG);
        $subj_dir = $afni_root_path .  $subjs[$subj] . "/";
        if(chdir($subj_dir) == 0) {
                if(mkdir($subj_dir) == 0) { log_file("Bad file path specificied for subject $subjs[$subj]: $subj_dir",ERROR); }
        }

        if (TEMP) {  # fill in with some subject-specific code here, if necessary


        }


        if (ADD_EDGE) {  # used to check quality of results of lpc_align.py
                        $cmd = "rm -f _ae._lpc* _ae.ExamineList.log";
                        execute_afni_command($cmd);
                        $cmd = "\@AddEdge _lpc.all_norm+orig _lpc.t1+orig t1-aligned+orig.";
                        execute_afni_command($cmd);
        }

        $contrast_dirname = $subj_dir."Contrasts/";
        if(chdir($contrast_dirname) == 0) {
                if(mkdir($contrast_dirname) == 0) { log_file("Bad file path specificied for contrast folder: $contrast_dirname",ERROR); }
                chdir($contrast_dirname);
        }

        open(STARTPOINTFILE, ">"."startpoints.1D") || log_file("Unable to open startpoints.1D for writing",ERROR);
        print STARTPOINTFILE "0\n";
        $totalTRs = 0;
        for ($run_num = 1; $run_num < $EPI_count[$subj]; $run_num++) {
                $totalTRs = $totalTRs + $num_scans[$subj][$run_num];   # note: look at N-1 runs, because don't need to know end of last run
                print STARTPOINTFILE $totalTRs . "\n";
        }
        close (STARTPOINTFILE);

        if (BASELINE_POLYNOMIAL_MODE == AUTO) {
                $baseline_model_order[$subj] = 1+ int((($totalTRs + $num_scans[$subj][$run_num]) * $TR_duration / $EPI_count[$subj]) / 150);  # 3dDeconlve baseline is 1 + 1 for every 150s. Use average nmber of TRs per run per subject
        } else {         # manual
                $baseline_model_order[$subj] = $baseline_model_order[0];
               }

        for ($contrast_num = 1; $contrast_num <= $num_contrasts; $contrast_num++) {
                        ## REgression Contrasts ##
                open(CONTRAST_FILE, ">".$contrast_name[$contrast_num] . ".1D") || log_file("Can't open contrast file for writing: $contrast_name[$contrast].1D",ERROR);
                print CONTRAST_FILE $EPI_count[$subj]*(1+$baseline_model_order[$subj]) . "\@0 ";   # 1 leading zero (regressor) per run plus 1 each for each term of baseline model (linear = 1, quadratic = 2, etc)
                for ($regressor = 1; $regressor <= $regressor_count; $regressor++) {
                        print CONTRAST_FILE $contrast[$contrast_num][$regressor] . " ";
                }
                print CONTRAST_FILE NUM_MOTION_PARAMETERS ."\@0";
                close(CONTRAST_FILE);

                ## Stick Deconvolution Contrasts ##
                open(CONTRAST_FILE, ">".$contrast_name[$contrast_num] . "_stick.1D") || log_file("Can't open contrast file for writing: $contrast_name[$contrast]_stick.1D",ERROR);
                print CONTRAST_FILE $EPI_count[$subj]*(1+$baseline_model_order[$subj]) . "\@0 ";   # 1 leading zero (regressor) per run plus 1 each for each term of baseline model (linear = 1, quadratic = 2, etc)
                for ($regressor = 1; $regressor <= $regressor_count; $regressor++) {
                        for ($HRF_regressor = 0; $HRF_regressor <= $HRF_duration_TR; $HRF_regressor++) {
                                print CONTRAST_FILE $contrast[$contrast_num][$regressor] . " ";
                        }
                }
                print CONTRAST_FILE NUM_MOTION_PARAMETERS ."\@0";
                close(CONTRAST_FILE);

        }

        $dicom_dir = $dicom_root_path . $dicom_directory[$subj] ."/";
        chdir($subj_dir);
        if (DICOM_TO_BRIK == ON) {
                check_file_exists($subj_dir."t1+orig.HEAD");
                log_file ("Creating Brik/Head files for $subjs[$subj] ",LOG);
                execute_afni_command("to3d -prefix t1 -session ".$subj_dir. " " . $dicom_dir . "*_" . $MPRAGE_series[$subj] ."_*");
        }

        if (DICOM_TO_BRIK == ON) {
                for ($run_num = 1; $run_num <= $EPI_count[$subj]; $run_num++) {

                        opendir(DIR, $dicom_dir);
                        $cmd = "\@runs = grep(/_$EPI_series[$subj][$run_num]_/,readdir(DIR))";
                        eval $cmd;
                        closedir(DIR);

                        check_file_exists($subj_dir . "run".$run_num."+orig.HEAD");
                        chdir($dicom_dir);

                        if ($run_num == 1) {       # get dicom header info
                                $cmd = "dicom_hdr -sexinfo ".$runs[1]." > EPI_dicom_header.txt";
                                execute_afni_command($cmd);
                                unless(open(DICOM_HEADER_FILE, "EPI_dicom_header.txt")) {
                                        log_file("Can't get EPI dicom header info",ERROR);
                                }
                                while ($line = <DICOM_HEADER_FILE>) {     # parse file
                                        if ($line =~ /ACQ Repetition Time\/\/(\d+)/) {
                                                $TR_duration_ms = $1;
                                        }
                                        if ($line =~ /sSliceArray\.lSize\s+=\s+(\d+)/) {
                                        print "getting num slices\n";
                                                $num_slices_z_axis = $1;
                                                last;
                                        }
                                }
                                close (DICOM_HEADER_FILE);
                        }
                        $cmd = "to3d -fim -time:zt $num_slices_z_axis $num_scans[$subj][$run_num] $TR_duration_ms alt+z2 -session ". $subj_dir . " -prefix run".$run_num. " ". join(" ",@runs);
                        #### Create fim files for run x with 30 levels for the z direction, $scan_number in the time direction, reading alternating slices in z direction starting at slice 1
                        execute_afni_command($cmd);

                        if ($run_num == $EPI_count[$subj]) { log_file("Processed " .$run_num ." EPIs for subject $subj\n",LOG); }
                }
        }  # Dicom to BRIK fim file



        chdir ($subj_dir);
        opendir(DIR, $subj_dir);
        @runs = grep(/run.{1,2}\+orig\.HEAD/, readdir(DIR));
        closedir(DIR);

                ### DESPIKE PRIOR to Volume registration - for serious noise issues ###
        if (DESPIKE_PRE_VOLREG == ON) {
                if (DESPIKE_POST_VOLREG == ON) {
                        log_file("Despiking on for both pre & post volume reg. Can only have one or other. ", ERROR)
                }
                $pre_vr_despike_flag = "_ds";
                foreach $r (@runs) {
                        if ($r =~ /(run.*)\+orig.*HEAD/) { $pfx = $1; }
                        check_file_exists($pfx . $pre_vr_despike_flag . "+orig.HEAD");
                        $cmd = "3dDespike -prefix " . $pfx . $pre_vr_despike_flag . " $r" ;
                        execute_afni_command($cmd);
                }
        }

                ### REALIGN ###
        @pfxs = ();
        $catcmd = "cat ";
        opendir(DIR, $subj_dir);                                ### Check to see if we have already made a despiked 3D+t orig file ####
        @despike = grep(/run.{1,2}_ds\+orig.HEAD/, readdir(DIR));
        closedir(DIR);
        if ((scalar(@despike) == scalar(@runs)) and USE_PREEXISTING_DESPIKED_DATA_IN_REALIGN and not(DESPIKE_POST_VOLREG)) { $pre_vr_despike_flag = "_ds"; }   ### use despiked data - all despiked files present; otherise use non-despiked ###

        foreach $r (@runs) {
                if ($r =~ /(run.*)\+orig.*HEAD/) {
                        $pfx = $1;
                        push(@pfxs, $pfx);
                }

                if (REALIGN == ON) {

                        ### if we ever develop a sagital EPI, can add 2dImReg about here   see afni.nimh.nih.gov/pub/dist/ASTON/afni10_volreg.pdf

                        check_file_exists($pfx . $pre_vr_despike_flag . "_vr+orig.HEAD");
                        $cmd = "3dvolreg -heptic -twopass -zpad 4 -tshift 4 -prefix " . $pfx . $pre_vr_despike_flag . "_vr -dfile " . $pfx . $pre_vr_despike_flag . "_vr.1D -maxdisp1D " . $pfx . $pre_vr_despike_flag ."_maxdisplacement ";
                        $cmd = $cmd . "-base " . $pfxs[0] . "+orig\'[".TR_FOR_REALIGNMENT."]\' -verbose " . $pfx . $pre_vr_despike_flag ."+orig";

                        execute_afni_command($cmd);

                        $catcmd = $catcmd . " " . $pfx .  $pre_vr_despike_flag . "_vr.1D";
                }
        }

        if (REALIGN == ON) {
                ### CONCATENATE MOTION ESTIMATE FILES ###
                $catcmd = $catcmd . " > vr_params.1D";
                print $catcmd . "\n";
                `$catcmd` unless DEBUG_ONLY;
        }



                ### DESPIKE after Volume registration - for minor spiking issues ###
        if (DESPIKE_POST_VOLREG == ON) {
                if (DESPIKE_PRE_VOLREG == ON) {
                        log_file("Despiking on for both pre & post volume reg. Can only have one or other. ", ERROR)
                }
                opendir(DIR, $subj_dir);                                ### Check to see if we have already made a despiked 3D+t orig file ####
                @volregfiles = grep(/run.{1,2}_vr\+orig.HEAD/, readdir(DIR));  #### this will despike _ds files if they are in the direction. Fix!!!
                closedir(DIR);
                if (scalar(@volregfiles) != scalar(@runs))  { log_file ("Number of volume registered runs dirrent than number of pre volumen registered runs", ERROR) };
                $post_vr_despike_flag = "_ds2";
                foreach $r (@volregfiles) {
                        if ($r =~ /(run.{1,2}vr)\+orig.*HEAD/) { $pfx = $1; }
                        check_file_exists($pfx . $post_vr_despike_flag . "+orig.HEAD");
                        $cmd = "3dDespike -prefix " . $pfx . $post_vr_despike_flag . " $r" ;
                        execute_afni_command($cmd);
                }
        }


        ### SPATIAL SMOOTHING ###
        if (SPATIAL_SMOOTHING == ON) {
        		$spatial_smooth_flag = "_sm".$spatial_blur_diameter;
                foreach $p (@pfxs) {
                        check_file_exists($p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . $spatial_smooth_flag."+orig.HEAD");
                        $cmd = "3dmerge -1blur_fwhm ".$spatial_blur_diameter." -doall ";
                        $cmd = $cmd . "-prefix " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "_".$spatial_smooth_flag ." ";
                        $cmd = $cmd . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "+orig";
                        execute_afni_command($cmd);
                        
                }
        }
        else {
         # need to add code to detect for spatial smoothing files
        	$spatial_smooth_flag = "";
        }

        ### CREATE BINARY MASK FROM LAST RUN ###
        ### Stuff to improve mask -clfrac at levels 0.2 0.25, or 0.33; -dilate
        if (BINARY_MASK == ON) {
                check_file_exists("mask+orig.HEAD");
                $cmd = "3dAutomask -prefix mask " . $pfxs[$#pfxs] .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . $spatial_smooth_flag."+orig";
                execute_afni_command($cmd);
        }

        ### CALCULATE PERCENT SIGNAL CHANGE and MEAN INTENSITIES ###
        if (CALC_SIGNAL_CHANGE == ON) {
                foreach $p (@pfxs) {
                        check_file_exists($p."_mean+orig.HEAD");
                        $cmd = "3dTstat -prefix " . $p . "_mean " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . $spatial_smooth_flag."+orig";
                        execute_afni_command($cmd);
                        check_file_written($p."_mean+orig.HEAD");
                }

                ### CALCULATE PCT SIGNAL CHANGE ###
                $allnorm = "";
                foreach $p (@pfxs) {
                        check_file_exists($p. $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . $spatial_smooth_flag."_norm+orig.HEAD");
                        $cmd = "3dcalc -a " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag .$spatial_smooth_flag."+orig ";
                        $cmd = $cmd . "-b " . $p . "_mean+orig ";
                        $cmd = $cmd . "-expr \"(a/b*100)\" ";                                # normalize to 100% signal
                        #$cmd = $cmd . "-expr \"((a-b)/b*100)\" ";                        #(used in NStroop)
                        $cmd = $cmd . "-prefix " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . $spatial_smooth_flag."_norm";
                        execute_afni_command($cmd);
                        check_file_written($p. $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . $spatial_smooth_flag."_norm+orig.HEAD");

                        $allnorm = $allnorm . " " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . $spatial_smooth_flag."_norm+orig";
                }

                ### CONCATENATE NORMALIZED DATA
                check_file_exists("all_norm+orig.HEAD");
                $cmd = "3dTcat -prefix all_norm " . $allnorm;
                execute_afni_command($cmd);
                check_file_written("all_norm+orig.HEAD");
        }  ### end of calculating percent signal change


        ### Align/warp anat data to epi ###   need to add iteration on run
        if (MPRAGE_EPI_ALLIGNMENT == ANATNUDGE) {   ########  We should  move this next too volume registation
                if (SKULL_STRIP == ON) { skull_strip() }
                check_file_exists("t1-nudge0+orig.HEAD");
                $cmd =  "3dAnatNudge -anat t1-ss+orig -epi run1+orig\'[".TR_FOR_REALIGNMENT."]\' -x 5 -y 5 -z 5 -prefix t1-nudge0 "; #-verb
                        execute_afni_command($cmd);

                        for ($i=0; $i<=40; $i++) {
                                        $j = $i + 1;
                                check_file_exists("t1-nudge".$j."+orig.HEAD");
                                        $cmd = "3dAnatNudge -anat t1-nudge$i+orig -epi run1+orig\'[".TR_FOR_REALIGNMENT."]\' -x 5 -y 5 -z 5 -prefix t1-nudge$j "; #-verb
                                        execute_afni_command($cmd);
                        }
                        get_max_nudge();
                        check_file_exists("t1-aligned+orig.HEAD");
                        check_file_exists("t1-aligned+orig.BRIK");
                        $cmd = "3drename t1-nudge".$max_nudge."+orig t1-aligned";
                        execute_afni_command($cmd);
        }
        if (MPRAGE_EPI_ALLIGNMENT == ALLINEATE) {   # better than nudge
                #Ziad's script
                        check_file_exists("t1_alepi+orig.HEAD");
                        $cmd = "lpc_align.py -anat t1+orig -epi run1+orig'[".TR_FOR_REALIGNMENT."]'";
                        execute_afni_command($cmd);

                        check_file_exists("t1-aligned+orig.HEAD");
                        check_file_exists("t1-aligned+orig.BRIK");
                        $cmd = "3drename t1_alepi+orig t1-aligned";
                        execute_afni_command($cmd);
        }

        ### AUTO-TALAIRACH EARLY - Pre 3dDeconvolve ###
        if (TALAIRACH_MODE_EARLY != OFF) {
                #get_max_nudge();

                if (TALAIRACH_MODE_EARLY == AUTO) {
                        autotalairach_structural_and_mask();

                        #Talairach functional/EPI
                        check_file_exists("all_norm_at+tlrc.HEAD");
                        check_file_exists("all_norm_at+tlrc.BRIK");

                        $cmd = "\@auto_tlrc -apar t1-aligned_at+tlrc -input all_norm+orig -dxyz 3";
                        execute_afni_command($cmd);
                        check_file_written("all_norm_at+tlrc.HEAD");

                }
                if (TALAIRACH_MODE_EARLY == MANUAL) {
                        manualtalairach_structural_and_mask();

                        #Talairach functional/EPI
                        check_file_exists("all_norm_mt+tlrc.HEAD");
                        check_file_exists("all_norm_mt+tlrc.BRIK");
                        $cmd = "adwarp -apar t1+tlrc. -dpar all_norm+orig -dxyz 3 -prefix all_norm_mt";
                        execute_afni_command($cmd);
                        check_file_written("all_norm_mt+tlrc.HEAD");
                }
        }     # End TALAIRACH Processing

        if (WAVER == ON) {
                for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                        $concat_string_ideal = "cat ";
                        $concat_string_timings = "cat ";
                        $concat_string_timings_shifted_for_iresp = "cat ";

                        for ($run_num = 1; $run_num <= $EPI_count[$subj]; $run_num++) {
                                $input_reg_file = $subjs[$subj]."_Rgrsr_Run".$run_num.$regressor_root_fname[$regressor_num];
                                $shifted_for_iresp_reg_file = $subjs[$subj]."_Rgrsr_Shifted_Run".$run_num.$regressor_root_fname[$regressor_num];

                                ### Create shifted files that that the impulse response shape can be calculated starting n TRs early (n specified in info.txt file)
                                if(open(SHIFTED_IRESP_REG_FILE, ">".$shifted_for_iresp_reg_file)) {
                                                #
                                } else {log_file("Can't open $shifted_for_iresp_reg_file for writing",ERROR);}
                                open(INPUT_REG_FILE, $input_reg_file) or log_file("Can't open $input_reg_file for reading",ERROR);
                                $line_num = 0;
                                foreach $line (<INPUT_REG_FILE>) {
                                        if (++$line_num > $TRs_pre_regressor_for_estimating_baseline_time_courses) {

                                                if ($line =~ /^\W*([01])/) {                 # look for 0 or 1 - exclude special character at end of file
                                                        print SHIFTED_IRESP_REG_FILE "$1\n";
                                                }
                                        }
                                }
                                for ($pre_reg_TR = 0; $pre_reg_TR < $TRs_pre_regressor_for_estimating_baseline_time_courses; $TRs_pre_regressor_for_estimating_baseline_time_courses++) {
                                      print SHIFTED_IRESP_REG_FILE "0\n"  ;
                                }

                                close(SHIFTED_IRESP_REG_FILE);
                                close(INPUT_REG_FILE);


                                $output_reg_file = $subjs[$subj]."_IdealWaveform_Run".$run_num.$regressor_root_fname[$regressor_num];


                                $cmd = "waver -GAM -dt 2.0 -numout $num_scans[$subj][$run_num] -input $input_reg_file > $output_reg_file";

                                execute_afni_command($cmd);

                                check_file_written($output_reg_file);

                                $concat_string_ideal = $concat_string_ideal.$output_reg_file." " ;
                                $concat_string_timings = $concat_string_timings.$input_reg_file." ";
                                $concat_string_timings_shifted_for_iresp = $concat_string_timings_shifted_for_iresp.$shifted_for_iresp_reg_file." ";

                        }

                        $cmd = $concat_string_ideal . "> ".$subjs[$subj]."_IdealWaveform_AllRuns".$regressor_root_fname[$regressor_num];
                        $cmd2 = $concat_string_timings . "> ".$subjs[$subj]."_StimulusTiming_AllRuns".$regressor_root_fname[$regressor_num];
                        $cmd3 = $concat_string_timings_shifted_for_iresp. "> ".$subjs[$subj]."_StimulusTiming_Shifted_AllRuns".$regressor_root_fname[$regressor_num];

                execute_afni_command($cmd);
                execute_afni_command($cmd2);
                execute_afni_command($cmd3);
                }   # end of loop on regressor_num
        }    # End of WAVER



        if (GENERAL_LINEAR_TEST == ON or STICK_FUNCTION_GLT == ON or ESTIMATE_IMPULSE_RESPONSE == ON or GENERAL_LINEAR_TEST_3BackBlock_only == ON or GENERAL_LINEAR_TEST_3BackBlock_only_no_baseline == ON) {   # run any of the 3dDeconvolve functions
                initialize_3dDeconv_variables();

                if (TALAIRACH_MODE_LATE != OFF) {
                        if (TALAIRACH_MODE_EARLY != OFF) { log_file("Can't have both early and late Talairaching on",ERROR); }
                        $deconv_input_file = "all_norm+orig";
                        $deconv_mask_file = "mask+orig";
                        $output_coorindate_system = "+orig"
                } else {
                        $deconv_input_file = "all_norm_".$talrch_type."+tlrc";
                        $deconv_mask_file = "mask_".$talrch_type."+tlrc";
                        $output_coorindate_system = "+tlrc"
                }


                if (GENERAL_LINEAR_TEST == ON) {
                        #$num_input_timeseries = $regressor_count + NUM_MOTION_PARAMETERS;

                        $deconv_bucket_name = $studyname."_deconv";
#                        $cmd = "3dDeconvolve -jobs 2 -nobout -input $deconv_input_file -mask $deconv_mask_file -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries $CensoredVoxelsFileString \\\n"; #-censor $CensoredVoxelsFileString
                        $cmd = "3dDeconvolve -jobs 2 -nobout -input $deconv_input_file -mask $deconv_mask_file -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries $CensoredVoxelsFileString \\\n"; #-censor $CensoredVoxelsFileString

                        if (JUST_REGRESS_OUT_BASELINE == OFF) {
                                for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                        $reg_file = $subjs[$subj]."_IdealWaveform_AllRuns".$regressor_root_fname[$regressor_num];
#                                        $cmd = $cmd . " -stim_file $regressor_num $reg_file -stim_label $regressor_num $regressor_root_fname[$regressor_num] \\\n";
                                        $cmd = $cmd . " -stim_file $regressor_num $reg_file -stim_label $regressor_num $regressor_indicator_char$regressor_name[$regressor_num] \\\n";
                                }
                                for ($regressor_num = $regressor_count +1 ; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                                        $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - $regressor_count)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                                }

                                $cmd = $cmd . " -num_glt $num_contrasts \\\n";

                                for ($contrast = 1; $contrast <= $num_contrasts; $contrast++) {
                                        $cmd = $cmd . " -glt 1 Contrasts/$contrast_name[$contrast].1D -glt_label $contrast $contrast_name[$contrast] \\\n";
                                }
                        } else {
                                $cmd = $cmd . "-fitts " . $studyname."_model";  # write the actual baseline model
                        }

                        check_file_exists($deconv_bucket_name.$output_coorindate_system.".HEAD");
                        check_file_exists($studyname."_residual".$output_coorindate_system.".HEAD");

                        $cmd = $cmd . " -errts " . $studyname."_residual -bucket $deconv_bucket_name -fout -tout -rout > ". $studyname . ".out";
                        execute_afni_command($cmd);
                        check_file_written($deconv_bucket_name.$output_coorindate_system.".HEAD");

                        late_talairach($deconv_bucket_name);

                }                # End of GENERAL_LINEAR _TEST

                if (GENERAL_LINEAR_TEST_3BackBlock_only == ON) {  # just compute the 3Back
                        $num_input_timeseries = 1 + NUM_MOTION_PARAMETERS;
                        $deconv_bucket_name = $studyname."_3BackBlock_deconv";
                        $cmd = "3dDeconvolve -jobs 2 -nobout -input $deconv_input_file -mask $deconv_mask_file -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries $CensoredVoxelsFileString \\\n"; #-censor $CensoredVoxelsFileString
                        $cmd = $cmd . " -stim_file 1 ".$subjs[$subj]."_IdealWaveform_AllRuns_3-Back.1D -stim_label 1 ".$regressor_indicator_char."3BackBlock";
                        # add instructions here
                        for ($regressor_num = 2; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                                $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - 1)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                        }
                        $cmd = $cmd . " -num_glt $num_contrasts \\\n";
                        check_file_exists($deconv_bucket_name.$output_coorindate_system.".HEAD");
                        check_file_exists($studyname."_residual".$output_coorindate_system.".HEAD");
                        $cmd = $cmd . " -errts " . $studyname."_residual -bucket $deconv_bucket_name -fout -tout -rout > ". $studyname . ".out";
                        execute_afni_command($cmd);
                        check_file_written($deconv_bucket_name.$output_coorindate_system.".HEAD");
                        late_talairach($deconv_bucket_name);
                }                # End of GENERAL_LINEAR _TEST


                if (GENERAL_LINEAR_TEST_3BackBlock_only_no_baseline == ON) {  # just compute the 3Back
                        if (INSTRUCTIONS == ON) {
                            $num_input_timeseries = 1 + NUM_MOTION_PARAMETERS;
                            $deconv_bucket_name = $studyname."_3BackBlock_NoBaseline_deconv";
                            $cmd = "3dDeconvolve -jobs 2 -nobout -input $deconv_input_file -mask $deconv_mask_file -polort 0 -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries $CensoredVoxelsFileString \\\n"; #-censor $CensoredVoxelsFileString
                            $cmd = $cmd . " -stim_file 1 ".$subjs[$subj]."_IdealWaveform_AllRuns_3-Back.1D -stim_label 1 ".$regressor_indicator_char."3BackBlock";
                            for ($regressor_num = 2; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                                    $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - 1)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                            }
                            $cmd = $cmd . " -num_glt $num_contrasts \\\n";
                            check_file_exists($deconv_bucket_name.$output_coorindate_system.".HEAD");
                            check_file_exists($studyname."_residual".$output_coorindate_system.".HEAD");
                            $cmd = $cmd . " -errts " . $studyname."_residual -bucket $deconv_bucket_name -fout -tout -rout > ". $studyname . ".out";
                            execute_afni_command($cmd);
                            check_file_written($deconv_bucket_name.$output_coorindate_system.".HEAD");
                            late_talairach($deconv_bucket_name);
                        } else {
                            $num_input_timeseries = 1 + NUM_MOTION_PARAMETERS;
                            $deconv_bucket_name = $studyname."_3BackBlock_NoBaseline_deconv";
                            $cmd = "3dDeconvolve -jobs 2 -nobout -input $deconv_input_file -mask $deconv_mask_file -polort 0 -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries $CensoredVoxelsFileString \\\n"; #-censor $CensoredVoxelsFileString
                            $cmd = $cmd . " -stim_file 1 ".$subjs[$subj]."_IdealWaveform_AllRuns_3-Back.1D -stim_label 1 ".$regressor_indicator_char."3BackBlock";
                            for ($regressor_num = 2; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                                    $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - 1)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                            }
                            $cmd = $cmd . " -num_glt $num_contrasts \\\n";
                            check_file_exists($deconv_bucket_name.$output_coorindate_system.".HEAD");
                            check_file_exists($studyname."_residual".$output_coorindate_system.".HEAD");
                            $cmd = $cmd . " -errts " . $studyname."_residual -bucket $deconv_bucket_name -fout -tout -rout > ". $studyname . ".out";
                            execute_afni_command($cmd);
                            check_file_written($deconv_bucket_name.$output_coorindate_system.".HEAD");
                            late_talairach($deconv_bucket_name);
                        }    
                }                # End of GENERAL_LINEAR _TEST



# the below can probably be removed at some point
                if (GENERAL_LINEAR_TEST_BASELINEMODEL_ONLY == ON) {
                        #$num_input_timeseries = $regressor_count + NUM_MOTION_PARAMETERS;
                        $deconv_bucket_name = $studyname."_deconv";

                        $cmd = "3dDeconvolve -jobs 4 -nobout -input $deconv_input_file -mask $deconv_mask_file -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries $CensoredVoxelsFileString \\\n"; #-censor $CensoredVoxelsFileString

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                $reg_file = $subjs[$subj]."_IdealWaveform_AllRuns".$regressor_root_fname[$regressor_num];
                                $cmd = $cmd . " -stim_file $regressor_num $reg_file -stim_label $regressor_num $regressor_indicator_char$regressor_name[$regressor_num] \\\n";
                        }
                        for ($regressor_num = $regressor_count +1 ; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                                $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - $regressor_count)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                        }

                        $cmd = $cmd . " -num_glt $num_contrasts \\\n";

                        for ($contrast = 1; $contrast <= $num_contrasts; $contrast++) {
                                $cmd = $cmd . " -glt 1 Contrasts/$contrast_name[$contrast].1D -glt_label $contrast $contrast_name[$contrast] \\\n";
                        }

                        check_file_exists($deconv_bucket_name.$output_coorindate_system.".HEAD");
                        check_file_exists($studyname."_residual".$output_coorindate_system.".HEAD");

                        $cmd = $cmd . " -errts " . $studyname."_residual -bucket $deconv_bucket_name -fout -tout -rout";
                        execute_afni_command($cmd);
                        check_file_written($deconv_bucket_name.$output_coorindate_system.".HEAD");

                        late_talairach($deconv_bucket_name);

                }                # End of GENERAL_LINEAR _TEST


                if (STICK_FUNCTION_GLT == ON) {

                        #$num_input_timeseries = $regressor_count + NUM_MOTION_PARAMETERS;

                        # -tshift added 2-10-08 (just prior to completeion of AE paper; -tshift not included in AE paper iresp data
                        $cmd = "3dDeconvolve -jobs 2 -input $deconv_input_file -mask $deconv_mask_file -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries -tshift \\\n";

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                $reg_file = $subjs[$subj]."_StimulusTiming_AllRuns".$regressor_root_fname[$regressor_num];
                                $cmd = $cmd . " -stim_file $regressor_num $reg_file -stim_label $regressor_num $regressor_indicator_char$regressor_name[$regressor_num] -stim_maxlag $regressor_num $HRF_duration_TR \\\n";
                        }

                        for ($regressor_num = $regressor_count +1 ; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                                $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - $regressor_count)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                        }

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                check_file_exists("stick_iresp_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("stick_iresp_stdev_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("stick_iresp_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("stick_iresp_stdev_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");

                                $cmd = $cmd . " -iresp $regressor_num stick_iresp_$regressor_name[$regressor_num] -sresp $regressor_num stick_iresp_stdev_$regressor_name[$regressor_num] \\\n";
                        }


                        $cmd = $cmd . " -num_glt $num_contrasts \\\n";

                        for ($contrast = 1; $contrast <= $num_contrasts; $contrast++) {
                                $cmd = $cmd . " -glt 1 Contrasts/$contrast_name[$contrast]_stick.1D -glt_label $contrast $contrast_name[$contrast]_stick \\\n";
                        }

                        $deconv_bucket_name = $studyname."_stick_deconv";
                        check_file_exists($deconv_bucket_name.$output_coorindate_system.".HEAD");
                        $cmd = $cmd . " -bucket $deconv_bucket_name -fout -tout -rout";

                        execute_afni_command($cmd);

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                check_file_written("stick_iresp_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_written("stick_iresp_stdev_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                        }

                        late_talairach($deconv_bucket_name);

                }                  # end of STICK_FUNCTION_GLT


                if (ESTIMATE_IMPULSE_RESPONSE == ON) {   ### This is currently outfitted to inlude several pre-event TRs

                        #$num_input_timeseries = $regressor_count + NUM_MOTION_PARAMETERS;

                        # -tshift added 2-10-08 (just prior to completeion of AE paper; -tshift not included in AE paper iresp data
                        $cmd = "3dDeconvolve -jobs 2 -input $deconv_input_file -mask $deconv_mask_file -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries -tshift \\\n";

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                $reg_file = $subjs[$subj]."_StimulusTiming_Shifted_AllRuns".$regressor_root_fname[$regressor_num];

                                #### Count number of events ###
                                open(REGFILE, $reg_file) or log_file("Can't open $reg_file. Must do so to count number of impulse", WARNING);
                                while($line = <REGFILE>) {
                                        if($line =~ /^1/) {$impulse_count[$regressor_num]++;}
                                }
                                close(REGFILE);

                                $cmd = $cmd . " -stim_file $regressor_num $reg_file -stim_label $regressor_num $regressor_indicator_char$regressor_name[$regressor_num] -stim_maxlag $regressor_num ". sprintf("%d",$HRF_duration_TR + $TRs_pre_regressor_for_estimating_baseline_time_courses) ."\\\n";
                        }

                        for ($regressor_num = $regressor_count +1 ; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                                $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - $regressor_count)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                        }

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                check_file_exists("iresp_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("iresp_stdev_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("iresp_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("iresp_stdev_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");

                                $cmd = $cmd . " -iresp $regressor_num iresp_$regressor_name[$regressor_num] -sresp $regressor_num iresp_stdev_$regressor_name[$regressor_num] \\\n";

                        }

                        $deconv_bucket_name = $studyname."_peri_stimulus_iresp_deconv";
                        check_file_exists($deconv_bucket_name.$output_coorindate_system.".HEAD");
                        $cmd = $cmd . " -bucket $deconv_bucket_name";

                        execute_afni_command($cmd);

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                check_file_written("iresp_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_written("iresp_stdev_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("iresp_stder_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                check_file_exists("iresp_stder_$regressor_name[$regressor_num]".$output_coorindate_system.".HEAD");
                                $cmd = "3dcalc -a iresp_stdev_$regressor_name[$regressor_num]+tlrc -expr \'a/sqrt($impulse_count[$regressor_num]) \' -prefix iresp_stder_$regressor_name[$regressor_num]";
                                execute_afni_command($cmd);
                        }

                        late_talairach($deconv_bucket_name);

                }                  # end of ESTIMSTE_IMPULSE_RESPONSE
        }                 # End of execute 3dDeconvolve processes


                ### AUTO-TALAIRACH Late - Post 3dDeconvolve ###
                ### Most of the time this will be run in the 3dDeconvolve section
                ### Here, limited ability to talairach just the _deconv (not iresp)
        if (TALAIRACH_MODE_LATE != OFF) {
                if (GENERAL_LINEAR_TEST == OFF and STICK_FUNCTION_GLT == OFF and ESTIMATE_IMPULSE_RESPONSE == OFF) {   #  no 3dDeconvolve functions were run
                        if (TALAIRACH_MODE_EARLY != OFF) { log_file("Can't have both early and late Talairaching on",ERROR); }

                        late_talairach($studyname."_deconv");

                }
        }     # End Late TALAIRACH Processing

} ## End of process_subject


                 #############################
                 ###### Group Analysis #######
                 #############################
if ($group_analysis == TRUE) {

        for ($r = 0; $r < $regressor_count; $r++) {
                $reg_name_deconved[$r] = $regressor_indicator_char . $regressor_name[$r+1] ;   # names of basic regressor in deconvolve file include underscore and ".1D"
        }

        @test_names = (@reg_name_deconved, @contrast_name[1..@contrast_name]);        # we want maps for both the regressors and the contrasts

        $num_tests = $#test_names; #  index number for last element in array test_names


                ### Create group Regression stats map: Loop through contrasts (and subjects) to create command strings **#
        if (GROUP_REGRESSION_STATS) {
                                        # First get subbriks
                $deconv_bucket_name = $studyname."_deconv+tlrc";
                for ($subj=0; $subj< $num_subjs; $subj++) {                                        # read
                        $cwd= $afni_root_path .  $subjs[$subj] . "/";
                        if(chdir($cwd) == 0) { log_file("Can't open directory $cwd",ERROR); }
                        execute_afni_command("3dinfo -verb $deconv_bucket_name > DeconvInfo.txt");
                        unless(open(DECONV_INFO_FILE, "DeconvInfo.txt")) { log_file("Can't open info file: $!" , ERROR); }
                        while ($line = <DECONV_INFO_FILE>) {
                                for ($test_num = 0; $test_num < $num_tests; $test_num++) {
                                                # this potential could be switched to the recently added 3dinfo -label2index <label> <dataset>

                                        if ($line =~ /(\d+) '$test_names[$test_num](_GLT)?#0_[cC]oef'/) { $test_betacoef_subbrik[$subj][$test_num] = $1;  }  # This Reg Expression still being tuned - may be bugs
                                }
                        }
                        close(DECONV_INFO_FILE);
                }                         ## End subject loop to learn sub briks ##


                chdir($groupdir_name);
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($test_num = 0; $test_num < $num_tests; $test_num++) {
#                                $first_in_group = TRUE;
                                check_file_exists("$test_names[$test_num]_$unique_groups[$group_num]+tlrc.HEAD");
                                $ttest_cmd = "3dttest -prefix $test_names[$test_num]_$unique_groups[$group_num] -base1 0 -set2 ";
#                                $residual_cmd = "3dcalc -prefix @test_names[$test_num]_@unique_groups[$group_num]_residual ";   # for calculate 3dFWHM and AlphaSim
#                                $expression_cmd = "-expr \'stdev(";
                                for ($subj=0; $subj< $num_subjs; $subj++) {
                                        if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                                $char = chr(97+$subj);
                                                $ttest_cmd = $ttest_cmd . "../$subjs[$subj]/$deconv_bucket_name\[" .$test_betacoef_subbrik[$subj][$test_num] ."] ";
#                                                $residual_cmd = $residual_cmd . "-$char \'../$subjs[$subj]/$deconv_bucket_name\[" .$test_betacoef_subbrik[$subj][$test_num] ."]\' ";
#                                                if ($first_in_group == FALSE) { $expression_cmd = $expression_cmd . ",$char"; }
#                                                else { $expression_cmd = $expression_cmd . "$char"; }
#                                                $first_in_group = FALSE;
                                        }
                                }

#                                $residual_cmd = $residual_cmd . $expression_cmd . ")**2\'";
                                if (PARRALLELIZE) {
                                        $groupstat_script_filename = "ttestscript_".$group_num."_".$test_num;
                                        if(open(GROUPSTAT_FILE, ">".$groupstat_script_filename)) {
                                                print GROUPSTAT_FILE $ttest_cmd . "\n";
#                                                print GROUPSTAT_FILE $residual_cmd ;
                                        } else {log_file("Can't open $groupstat_script_filename for writing",ERROR);}
                                        close(GROUPSTAT_FILE);
                                        `chmod 774 $groupstat_script_filename`;
                                        execute_afni_command("submit $groupdir_name\\$groupstat_script_filename");

#                                        $cmd = "echo $ttest_cmd > ttestscript_$group_num\_$test_num";
#                                        `$cmd`;
#                                        `chmod 774 ttestscript_$group_num\_$test_num`;
#                                        execute_afni_command("submit $groupdir_name\\ttestscript_$group_num\_$test_num");
                                } else {
                                        execute_afni_command($ttest_cmd);
#                                        execute_afni_command($residual_cmd);
                                }
                        }
                }
        }

                ### Create group stick deconcolution stats map: Loop through contrasts (and subjects) to create command strings **#
        if (GROUP_DECONV_STICK_STATS)        {
                                                        # First get subbriks
                $deconv_bucket_name = $studyname."_stick_deconv+tlrc";
                for ($subj=0; $subj< $num_subjs; $subj++) {                                        # read
                        $cwd= $afni_root_path .  $subjs[$subj] . "/";
                        if(chdir($cwd) == 0) { log_file("Can't open directory $cwd",ERROR); }

                        execute_afni_command("3dinfo -verb $deconv_bucket_name > DeconvInfo.txt");

                        unless(open(DECONV_INFO_FILE, "DeconvInfo.txt")) { log_file("Can't open info file: $!" , ERROR); }
                        while ($line = <DECONV_INFO_FILE>) {
                                $reg_count = 0;
                                for ($test_num = 0; $test_num < $num_tests; $test_num++) {
                                        $reg_count++;
                                        if ($reg_count <= $regressor_count) {   # if we are looking at the individual regressors, which have multiple time points (versus contrasts)
                                                if ($line =~ /(\d+) '$test_names[$test_num](_iresp)?(_stick)?(_GLT)?#2_[cC]oef'/) { $test_betacoef_subbrik[$subj][$test_num] = $1;  }

                                        } else {
                                                if ($line =~ /(\d+) '$test_names[$test_num](_iresp)?(_stick)?(_GLT)?#0_[cC]oef'/) { $test_betacoef_subbrik[$subj][$test_num] = $1;  }
                                        }
                                }
                        }
                        close(DECONV_INFO_FILE);
                }                         ## End subject loop to learn sub briks ##

                chdir($groupdir_name);
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($test_num = 0; $test_num < $num_tests; $test_num++) {
                                check_file_exists("$test_names[$test_num]_$unique_groups[$group_num]_stick+tlrc.HEAD");
                                $ttest_cmd = "3dttest -prefix $test_names[$test_num]_$unique_groups[$group_num]_stick -base1 0 -set2 ";

                                for ($subj=0; $subj< $num_subjs; $subj++) {
                                        if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                                $ttest_cmd = $ttest_cmd . "../$subjs[$subj]/$deconv_bucket_name\[" .$test_betacoef_subbrik[$subj][$test_num] ."] ";
                                        }
                                }
                                if (PARRALLELIZE) {
                                        $cmd = "echo $ttest_cmd > ttestscript_stick_$group_num\_$test_num";
                                        `$cmd`;
                                        `chmod 774 ttestscript_stick_$group_num\_$test_num`;
                                        execute_afni_command("submit $groupdir_name\\ttestscript_stick_$group_num\_$test_num");
                                } else { execute_afni_command($ttest_cmd); }
                        }
                }
        }


                ### Create group IRESP and WIthin-subject SEM map (loop through regressor and subject) ###
        if (GROUP_IRESP_ANALYSIS) {
                chdir($groupdir_name);
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                $first_in_group = TRUE;
                                $iresp_ave_script_filename = "iresp_ave_script_$group_num\_$regressor_num.sh";
                                $stder_withinsubj_script_filename = "stder_within_script_$group_num\_$regressor_num.sh";        ## within subject standard error - average of standard error within each subject
                                check_file_exists("iresp_$unique_groups[$group_num]$regressor_name[$regressor_num]+tlrc.HEAD");
                                check_file_exists("stder_withinsb_$unique_groups[$group_num]$regressor_name[$regressor_num]+tlrc.HEAD");
                                $iresp_ave_cmd = "3dcalc  -prefix iresp_$unique_groups[$group_num]_$regressor_name[$regressor_num] ";
                                $stder_withinsubj_cmd = "3dcalc  -prefix stder_withinsb_$unique_groups[$group_num]_$regressor_name[$regressor_num] ";
                                $mean_cmd = "-expr \'mean(";
                                for ($subj=0; $subj< $num_subjs; $subj++) {
                                        if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                                $char = chr(97+$subj);
                                                $iresp_ave_cmd = $iresp_ave_cmd . "-$char \'../$subjs[$subj]/iresp_$regressor_name[$regressor_num]+tlrc\' ";
                                                $stder_withinsubj_cmd = $stder_withinsubj_cmd . "-$char \'../$subjs[$subj]/iresp_stder_$regressor_name[$regressor_num]+tlrc\' ";
                                                if ($first_in_group == FALSE) { $mean_cmd = $mean_cmd . ",$char"; }
                                                else { $mean_cmd = $mean_cmd . "$char"; }
                                                $first_in_group = FALSE;
                                        }
                                }
                                $iresp_ave_cmd = $iresp_ave_cmd . $mean_cmd . ")\'";
                                $stder_withinsubj_cmd = $stder_withinsubj_cmd . $mean_cmd . ")\'";

                                if (PARRALLELIZE) {
                                        if(open(IRESP_AVE_FILE, ">".$iresp_ave_script_filename)) {
                                                print IRESP_AVE_FILE $iresp_ave_cmd;
                                        } else {log_file("Can't open $iresp_ave_script_filename for writing",ERROR);}
                                        close(IRESP_AVE_FILE);
                                        `chmod 774 $iresp_ave_script_filename`;
                                        execute_afni_command("submit $groupdir_name\\$iresp_ave_script_filename");

                                        if (WITHIN_SUBJECT_SEM == ON) {
                                                                                        if(open(STDER_WITHINSUBJ_FILE, ">".$stder_withinsubj_script_filename)) {
                                                                                                        print STDER_WITHINSUBJ_FILE $stder_withinsubj_cmd;
                                                                                        } else {log_file("Can't open $stder_withinsubj_script_filename for writing",ERROR);}
                                                                                        close(STDER_WITHINSUBJ_FILE);
                                                                                        `chmod 774 $stder_withinsubj_script_filename`;
                                                                                        execute_afni_command("submit $groupdir_name\\$stder_withinsubj_script_filename");
                                                                                }
                                } else {
                                        execute_afni_command($iresp_ave_cmd);
                                        if (WITHIN_SUBJECT_SEM == ON) {
                                                execute_afni_command($stder_withinsubj_cmd);
                                        }
                                }
                        }    # end for regressor_num
                }  # end for group_num
                                ### Now we iterate again just to execute between subject SEM calculation, which is dependant on IRESP
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                $first_in_group = TRUE;
                                $stder_betweensubj_script_filename = "stder_between_script_$group_num\_$regressor_num.pl";        ## between subject standard error
                                check_file_exists("stder_betweensb_$unique_groups[$group_num]_$regressor_name[$regressor_num]+tlrc.HEAD");
                                $stder_betweensubj_cmd = "3dcalc  -prefix stder_betweensb_$unique_groups[$group_num]_$regressor_name[$regressor_num] ";
                                $sem_cmd = "-expr \'sem(";
                                for ($subj=0; $subj< $num_subjs; $subj++) {
                                        if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                                $char = chr(97+$subj);
                                                $stder_betweensubj_cmd = $stder_betweensubj_cmd . "-$char \'../$subjs[$subj]/iresp_$regressor_name[$regressor_num]+tlrc\' ";
                                                if ($first_in_group == FALSE) {
                                                        $sem_cmd = $sem_cmd . ",$char";
                                                }
                                                else {
                                                        $sem_cmd = $sem_cmd . "$char";
                                                }
                                                $first_in_group = FALSE;
                                        }
                                }
                                $stder_betweensubj_cmd = $stder_betweensubj_cmd . $sem_cmd. ")\'";

                                if (PARRALLELIZE) {
                                        if(open(STDER_BETWEENSUBJ_FILE, ">".$stder_betweensubj_script_filename)) {
                                                print STDER_BETWEENSUBJ_FILE "#!/usr/bin/perl\n";
                                                print STDER_BETWEENSUBJ_FILE "\n\$i=0;\n";
                                                print STDER_BETWEENSUBJ_FILE "while(not(open(F, \"< \".\"iresp_$unique_groups[$group_num]_$regressor_name[$regressor_num]+tlrc.HEAD\"))) { \n";

                                                print STDER_BETWEENSUBJ_FILE "        sleep 10;\n";
                                                print STDER_BETWEENSUBJ_FILE "        if (\$i++ == 30) {last;}\n";
                                                print STDER_BETWEENSUBJ_FILE " }\n";
                                                print STDER_BETWEENSUBJ_FILE "`$stder_betweensubj_cmd`;\n";
                                        } else {log_file("Can't open $stder_betweensubj_script_filename for writing",ERROR);}
                                        close(STDER_BETWEENSUBJ_FILE);
                                        `chmod 774 $stder_betweensubj_script_filename`;
                                        execute_afni_command("submit $groupdir_name\\$stder_betweensubj_script_filename");
                                } else {
                                        execute_afni_command($stder_betweensubj_cmd);
                                }
                        }    # end for regressor_num
                }  # end for group_num

        } # END  if GROUP_IRESP_ANALYSIS

        if (ThreeD_ANOVA_AE == TRUE) {
                chdir($groupdir_name);
                $cmd = "3dANOVA3 -type 4 -alevels 2 -blevels 2 -clevels 16 \\\n";

                $Clevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Clevel_subject;
                                        $cmd = $cmd . "-dset 1 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[10] \\\n"; # Quinine 0
                                        $cmd = $cmd . "-dset 1 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[14] \\\n";        # Quinine 3
                                        $cmd = $cmd . "-dset 2 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[18] \\\n";  # Juice 0
                                        $cmd = $cmd . "-dset 2 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[22] \\\n";        # Juice 3
                                }
                        }
                }

                $cmd = $cmd . "-fa Solution \\\n";
                $cmd = $cmd . "-fb NBack \\\n";
                $cmd = $cmd . "-fab Solution_Nback_interact \\\n";
                $cmd = $cmd . "-amean 1 Quinine -amean 2 Juice \\\n";
                $cmd = $cmd . "-bmean 1 0back -bmean 2 3back \\\n";
                $cmd = $cmd . "-acont 1 -1 Spritz-Solution \\\n";   # Quinine-Juice
                $cmd = $cmd . "-bcont -1 1 Spritz-Nback \\\n";                # 0-3
                $cmd = $cmd . "-Abcont 1 : -1 1  Quinine-Nback \\\n";
                $cmd = $cmd . "-Abcont 2 : -1 1  Juice-Nback \\\n";

                $cmd = $cmd . "-bucket ANOVA-Spritz_undrgrd > ANOVA-Spritz_undrgrd.out \n";

                `rm -f ANOVA_Spritz_undrgrd*`;
                execute_afni_command($cmd);


                $cmd = "3dANOVA2 -type 3 -alevels 2 -blevels 16 \\\n";

                $Blevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Blevel_subject;
                                        $cmd = $cmd . "-dset 1 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[26] \\\n"; # Quinine 0
                                        $cmd = $cmd . "-dset 2 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[30] \\\n";  # Juice 0
                                }
                        }
                }

                $cmd = $cmd . "-fa Chckrbrd_Nback_F \\\n";
                $cmd = $cmd . "-amean 1 Flicker_0 \\\n";
                $cmd = $cmd . "-amean 2 Flicker_3 \\\n";
                $cmd = $cmd . "-adiff 2 1 Chckrbrd_Nback_cntrst \\\n";

                $cmd = $cmd . "-bucket ANOVA_Chkrbrd_undrgrd > ANOVA_Chkrbrd_undrgrd.out \n";

                `rm -f ANOVA_Chkrbrd_undrgrd*`;
                execute_afni_command($cmd);


                        ### NBack ###
                $cmd = "3dANOVA2 -type 3 -alevels 2 -blevels 16 \\\n";

                $Blevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Blevel_subject;
                                        $cmd = $cmd . "-dset 1 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[2] \\\n";
                                        $cmd = $cmd . "-dset 2 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[6] \\\n";
                                }
                        }
                }

                $cmd = $cmd . "-fa Probe_Nback_F \\\n";
                $cmd = $cmd . "-amean 1 Probe_0 \\\n";
                $cmd = $cmd . "-amean 2 Probe_3 \\\n";
                $cmd = $cmd . "-adiff 2 1 Probe_Nback_cntrst \\\n";

                $cmd = $cmd . "-bucket ANOVA_Probe_undrgrd > ANOVA_Probe_undrgrd.out \n";

                `rm -f ANOVA_Probe_results*`;
                execute_afni_command($cmd);

        } # END ThreeD_ANOVA_AE



        if (ThreeD_ANOVA_AE_PLOS == TRUE) {
                chdir($groupdir_name);
                $cmd = "3dANOVA3 -type 4 -alevels 4 -blevels 2 -clevels 16 \\\n";

                $Clevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Clevel_subject;
                                        $cmd = $cmd . "-dset 1 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[2] \\\n";  # Probe 0
                                        $cmd = $cmd . "-dset 1 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[6] \\\n";  # Probe 3
                                        $cmd = $cmd . "-dset 2 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[10] \\\n"; # Quinine 0
                                        $cmd = $cmd . "-dset 2 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[14] \\\n"; # Quinine 3
                                        $cmd = $cmd . "-dset 3 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[18] \\\n"; # Juice 0
                                        $cmd = $cmd . "-dset 3 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[22] \\\n"; # Juice 3
                                        $cmd = $cmd . "-dset 4 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[26] \\\n"; # Checkerboard 0
                                        $cmd = $cmd . "-dset 4 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[30] \\\n"; # Checkerboard 3     
                                }
                        }
                }

                $cmd = $cmd . "-fa Stimulus \\\n";
                $cmd = $cmd . "-fb Load \\\n";
                $cmd = $cmd . "-fab Stimulus_Load_interact \\\n";
                $cmd = $cmd . "-amean 1 Probe -amean 2 Quinine -amean 3 Juice -amean 4 Flicker \\\n";
                $cmd = $cmd . "-bmean 1 0back -bmean 2 3back \\\n";
                $cmd = $cmd . "-acontr 0 -1 1 0 Valence_ctrst \\\n";   # Quinine-Juice
                $cmd = $cmd . "-acontr -1 1 1 -1 liquidVSvisual_ctrst \\\n";            
                $cmd = $cmd . "-Abcontr 1 : -1 1 Load_Probe \\\n";                
                $cmd = $cmd . "-Abcontr 2 : -1 1 Load_Quinine \\\n";
                $cmd = $cmd . "-Abcontr 3 : -1 1 Load_Juice \\\n";
                $cmd = $cmd . "-Abcontr 4 : -1 1 Load_Flicker \\\n"; 



                $cmd = $cmd . "-mask FDR_mask+tlrc "    ;       

                $cmd = $cmd . "-bucket ANOVA_AllStimuli_PLOS > ANOVA_AllStimuli_PLOS.out \n";

                `rm -f ANOVA_AllStimuli_PLOS*`;
                execute_afni_command($cmd);

                $cmd = "3drefit -addFDR -FDRmask FDR_mask+tlrc ANOVA_AllStimuli_PLOS+tlrc";
                execute_afni_command($cmd);

        } # END ThreeD_ANOVA_AE_PLOS



        if (ThreeD_ANOVA_AE_STICK == TRUE) {
                chdir($groupdir_name);
                $cmd = "3dANOVA3 -type 4 -alevels 2 -blevels 2 -clevels 16 \\\n";

                $Clevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Clevel_subject;
                                        $cmd = $cmd . "-dset 1 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[36] \\\n"; # Quinine 0
                                        $cmd = $cmd . "-dset 1 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[52] \\\n";        # Quinine 3
                                        $cmd = $cmd . "-dset 2 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[68] \\\n";  # Juice 0
                                        $cmd = $cmd . "-dset 2 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[84] \\\n";        # Juice 3
                                }
                        }
                }

                $cmd = $cmd . "-fa Solution \\\n";
                $cmd = $cmd . "-fb NBack \\\n";
                $cmd = $cmd . "-fab Solution_Nback_interact \\\n";
                $cmd = $cmd . "-amean 1 Quinine -amean 2 Juice \\\n";
                $cmd = $cmd . "-bmean 1 0back -bmean 2 3back \\\n";
                $cmd = $cmd . "-acont 1 -1 Spritz-Solution \\\n";   # Quinine-Juice
                $cmd = $cmd . "-bcont 1 -1 Spritz-Nback \\\n";                # 0-3
                $cmd = $cmd . "-Abcont 1 : 1 -1  Quinine-Nback \\\n";
                $cmd = $cmd . "-Abcont 2 : 1 -1  Juice-Nback \\\n";

                $cmd = $cmd . "-bucket ANOVA-Spritz-baseline_undrgrd > ANOVA-Spritz-baseline_undrgrd.out \n";

                `rm -f ANOVA_Spritz_undrgrd*`;
                execute_afni_command($cmd);


                $cmd = "3dANOVA2 -type 3 -alevels 2 -blevels 16 \\\n";

                $Blevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Blevel_subject;
                                        $cmd = $cmd . "-dset 1 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[100] \\\n"; # Checkerborad 0
                                        $cmd = $cmd . "-dset 2 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[116] \\\n";  # Checkerboard 3
                                }
                        }
                }

                $cmd = $cmd . "-fa Chckrbrd_Nback_F \\\n";
                $cmd = $cmd . "-amean 1 Flicker_0 \\\n";
                $cmd = $cmd . "-amean 2 Flicker_3 \\\n";
                $cmd = $cmd . "-adiff 2 1 Chckrbrd_Nback_cntrst \\\n";

                $cmd = $cmd . "-bucket ANOVA_Chkrbrd-baseline_undrgrd > ANOVA_Chkrbrd-baseline_undrgrd.out \n";

                `rm -f ANOVA_Chkrbrd_undrgrd*`;
                execute_afni_command($cmd);


                        ### NBack ###
                $cmd = "3dANOVA2 -type 3 -alevels 2 -blevels 16 \\\n";

                $Blevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Blevel_subject;
                                        $cmd = $cmd . "-dset 1 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[4] \\\n"; # Probe 0
                                        $cmd = $cmd . "-dset 2 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[20] \\\n";  # Probe
                                }
                        }
                }

                $cmd = $cmd . "-fa Probe_Nback_F \\\n";
                $cmd = $cmd . "-amean 1 Probe_0 \\\n";
                $cmd = $cmd . "-amean 2 Probe_3 \\\n";
                $cmd = $cmd . "-adiff 2 1 Probe_Nback_cntrst \\\n";

                $cmd = $cmd . "-bucket ANOVA_Probe-baseline_undrgrd > ANOVA_Probe-baseline_undrgrd.out \n";

                `rm -f ANOVA_Probe_results*`;
                execute_afni_command($cmd);

        } # END ThreeD_ANOVA_AE_STICK


        if (ThreeD_ANOVA_AE_PLOS_STICK == TRUE) {
                chdir($groupdir_name);
                $cmd = "3dANOVA3 -type 4 -alevels 2 -blevels 2 -clevels 16 \\\n";

                $Clevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Clevel_subject;
                                        $cmd = $cmd . "-dset 1 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[36] \\\n"; # Quinine 0
                                        $cmd = $cmd . "-dset 1 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[52] \\\n";        # Quinine 3
                                        $cmd = $cmd . "-dset 2 1 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[68] \\\n";  # Juice 0
                                        $cmd = $cmd . "-dset 2 2 " . $Clevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[84] \\\n";        # Juice 3
                                }
                        }
                }

                $cmd = $cmd . "-fa Solution \\\n";
                $cmd = $cmd . "-fb NBack \\\n";
                $cmd = $cmd . "-fab Solution_Nback_interact \\\n";
                $cmd = $cmd . "-amean 1 Quinine -amean 2 Juice \\\n";
                $cmd = $cmd . "-bmean 1 0back -bmean 2 3back \\\n";
                $cmd = $cmd . "-acont 1 -1 Spritz-Solution \\\n";   # Quinine-Juice
                $cmd = $cmd . "-bcont 1 -1 Spritz-Nback \\\n";                # 0-3
                $cmd = $cmd . "-Abcont 1 : 1 -1  Quinine-Nback \\\n";
                $cmd = $cmd . "-Abcont 2 : 1 -1  Juice-Nback \\\n";

                $cmd = $cmd . "-bucket ANOVA_AllStimuli_stick_PLOS > ANOVA_AllStimuli_stick_PLOS.out \n";

                `rm -f ANOVA_AllStimuli_stick_PLOS*`;
                execute_afni_command($cmd);


                # $cmd = "3dANOVA2 -type 3 -alevels 2 -blevels 16 \\\n";

                # $Blevel_subject = 0;
                # for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                #         if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                #         for ($subj=0; $subj< $num_subjs; $subj++) {
                #                 if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                #                         ++$Blevel_subject;
                #                         $cmd = $cmd . "-dset 1 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[100] \\\n"; # Checkerborad 0
                #                         $cmd = $cmd . "-dset 2 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_stick_deconv+tlrc[116] \\\n";  # Checkerboard 3
                #                 }
                #         }
                # }

                # $cmd = $cmd . "-fa Chckrbrd_Nback_F \\\n";
                # $cmd = $cmd . "-amean 1 Flicker_0 \\\n";
                # $cmd = $cmd . "-amean 2 Flicker_3 \\\n";
                # $cmd = $cmd . "-adiff 2 1 Chckrbrd_Nback_cntrst \\\n";

                # $cmd = $cmd . "-bucket ANOVA_Chkrbrd-baseline_AllStimuli > ANOVA_Chkrbrd-baseline_AllStimuli.out \n";

                # `rm -f ANOVA_Chkrbrd_undrgrd*`;
                # execute_afni_command($cmd);


                        ### NBack ###
                # $cmd = "3dANOVA2 -type 3 -alevels 2 -blevels 16 \\\n";

                # $Blevel_subject = 0;
                # for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                #         if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                #         for ($subj=0; $subj< $num_subjs; $subj++) {
                #                 if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                #                         ++$Blevel_subject;
                #                         $cmd = $cmd . "-dset 1 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[4] \\\n"; # Probe 0
                #                         $cmd = $cmd . "-dset 2 " . $Blevel_subject . " ../$subjs[$subj]/".$studyname."_deconv+tlrc[20] \\\n";  # Probe
                #                 }
                #         }
                # }

                # $cmd = $cmd . "-fa Probe_Nback_F \\\n";
                # $cmd = $cmd . "-amean 1 Probe_0 \\\n";
                # $cmd = $cmd . "-amean 2 Probe_3 \\\n";
                # $cmd = $cmd . "-adiff 2 1 Probe_Nback_cntrst \\\n";

                # $cmd = $cmd . "-bucket ANOVA_Probe-baseline_undrgrd > ANOVA_Probe-baseline_undrgrd.out \n";

                # `rm -f ANOVA_Probe_results*`;
                # execute_afni_command($cmd);

        } # END ThreeD_ANOVA_AE_PLOS_STICK



        if (ThreeD_ANOVA_MED == TRUE) {
                chdir($groupdir_name);
                $cmd = "3dANOVA2 -type 3 -alevels 3 -blevels 12 \\\n";
                $Blevel_subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$Blevel_subject;
                                        $cmd = $cmd . "-dset 1 " . $Blevel_subject . " ../$subjs[$subj]/Medcue_deconv+tlrc[2] \\\n"; # Focused
                                        $cmd = $cmd . "-dset 2 " . $Blevel_subject . " ../$subjs[$subj]/Medcue_deconv+tlrc[6] \\\n";  # Open
                                        $cmd = $cmd . "-dset 3 " . $Blevel_subject . " ../$subjs[$subj]/Medcue_deconv+tlrc[10] \\\n";  # Compassion
                                }
                        }

                }
                $cmd = $cmd ."-fa MedState \\\n";
                $cmd = $cmd . "-amean 1 Focused \\\n";
                $cmd = $cmd . "-amean 2 Open \\\n";
                $cmd = $cmd . "-amean 3 Compassion \\\n";

                $cmd = $cmd . "-bucket ANOVA_MedState_meditator > ANOVA_MedState_meditator.out \n";
                `rm -f ANOVA_MedState_meditator*`;
                execute_afni_command($cmd);
        }

        if (ThreeD_ANOVA_ERIKSEN == TRUE) {
                chdir($groupdir_name);
#                $cmd = "3dANOVA2 -type 3 -alevels 3 -blevels 12 \\\n";
                $cmd = "3dANOVA3 -type 4 -alevels 3 -blevel 2 ";
                $subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        $cmd = $cmd . "-clevels " . $num_in_group[$group_num] . " \\\n";
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$subject;
                                        $cmd = $cmd . "-dset 1 1 " . $subject . " ../$subjs[$subj]/eriksen_med_deconv+tlrc[2] \\\n"; # Focused compat    the labels are probably wrong
                                        $cmd = $cmd . "-dset 2 1 " . $subject . " ../$subjs[$subj]/eriksen_med_deconv+tlrc[10] \\\n";  # Open compat
                                        $cmd = $cmd . "-dset 3 1 " . $subject . " ../$subjs[$subj]/eriksen_med_deconv+tlrc[18] \\\n";  # Compassion compat
                                        $cmd = $cmd . "-dset 1 2 " . $subject . " ../$subjs[$subj]/eriksen_med_deconv+tlrc[6] \\\n"; # Focused incompat
                                        $cmd = $cmd . "-dset 2 2 " . $subject . " ../$subjs[$subj]/eriksen_med_deconv+tlrc[14] \\\n";  # Open incompat
                                        $cmd = $cmd . "-dset 3 2 " . $subject . " ../$subjs[$subj]/eriksen_med_deconv+tlrc[22] \\\n";  # Compassion incompat
                                }
                        }
                }
                $cmd = $cmd ."-fab MedStateCompatInteract \\\n";
#                $cmd = $cmd . "-amean 1 Focused \\\n";
#                $cmd = $cmd . "-amean 2 Open \\\n";
#                $cmd = $cmd . "-amean 3 Compassion \\\n";

                $cmd = $cmd . "-bcontr -1 1 All_compatability \\\n";
                $cmd = $cmd . "-Abcontr 1 : -1 1 Foc_compatability \\\n";
                $cmd = $cmd . "-Abcontr 2 : -1 1 Neut_compatability \\\n";
                $cmd = $cmd . "-Abcontr 3 : -1 1 Open_compatability \\\n";


                $cmd = $cmd . "-bucket ANOVA_Eriksen_meditator > ANOVA_Eriksen_meditator.out & \n";
                `rm -f ANOVA_Eriksen_meditator*`;
                execute_afni_command($cmd);
        }        # end ThreeD_ANOVA_ERIKSEN == TRUE


                if (ThreeD_ANOVA_ERIKSEN_SEQADJ == TRUE) {
                chdir($groupdir_name);
                $cmd = "3dANOVA3 -type 4 -alevels 3 -blevel 4 ";
                $subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        $cmd = $cmd . "-clevels " . $num_in_group[$group_num] . " \\\n";
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$subject;
                                        $cmd = $cmd . "-dset 1 1 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[2] \\\n";   # F cC
                                        $cmd = $cmd . "-dset 1 2 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[6] \\\n";   # F cI
                                        $cmd = $cmd . "-dset 1 3 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[10] \\\n";  # F iC
                                        $cmd = $cmd . "-dset 1 4 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[14] \\\n";  # F iI
                                        $cmd = $cmd . "-dset 2 1 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[18] \\\n";  # N cC
                                        $cmd = $cmd . "-dset 2 2 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[22] \\\n";  # N cI
                                        $cmd = $cmd . "-dset 2 3 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[26] \\\n";  # N iC
                                        $cmd = $cmd . "-dset 2 4 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[30] \\\n";  # N iI
                                        $cmd = $cmd . "-dset 3 1 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[34] \\\n";  # O cC
                                        $cmd = $cmd . "-dset 3 2 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[38] \\\n";  # O cI
                                        $cmd = $cmd . "-dset 3 3 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[42] \\\n";  # O iC
                                        $cmd = $cmd . "-dset 3 4 " . $subject . " ../$subjs[$subj]/eriksen_PrevTrial_simple_deconv+tlrc[46] \\\n";  # O iI
                                }
                        }
                }

                $cmd = $cmd ."-fa Med \\\n";
                $cmd = $cmd ."-fb TrialType \\\n";
                $cmd = $cmd ."-fab Med-TrialTypeInteract \\\n";

#                $cmd = $cmd . "-amean 1 Focused \\\n";
#                $cmd = $cmd . "-amean 2 Neutral \\\n";
#                $cmd = $cmd . "-amean 3 Open \\\n";

#                $cmd = $cmd . "-bcontr -1 1 All_compatability \\\n";
                $cmd = $cmd . "-Abcontr 1 : 1 -3 1 1 Foc_iC_vs_others \\\n";
                $cmd = $cmd . "-Abcontr 2 : 1 -3 1 1 Neut_iC_vs_others \\\n";
                $cmd = $cmd . "-Abcontr 3 : 1 -3 1 1 Open_iC_vs_others \\\n";

                $cmd = $cmd . "-bucket ANOVA_Eriksen_SA_meditator > ANOVA_Eriksen_meditator.out \n";
                `rm -f ANOVA_Eriksen_SA_meditator*`;
                execute_afni_command($cmd);


                }   # end ThreeD_ANOVA_ERIKSEN_SEQADJ == TRUE




        if (ThreeD_ANOVA_NSTROOP == TRUE) {
                chdir($groupdir_name);
#                $cmd = "3dANOVA2 -type 3 -alevels 3 -blevels 12 \\\n";
                $cmd = "3dANOVA3 -type 4 -alevels 3 -blevel 3 ";
                $subject = 0;
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        $cmd = $cmd . "-clevels " . $num_in_group[$group_num] . " \\\n";
                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                        ++$subject;
                                        $cmd = $cmd . "-dset 1 1 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[2] \\\n"; # Focused compat
                                        $cmd = $cmd . "-dset 2 1 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[14] \\\n";  # Neut compat
                                        $cmd = $cmd . "-dset 3 1 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[26] \\\n";  # Open compat
                                        $cmd = $cmd . "-dset 1 2 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[6] \\\n"; # Focused incompat
                                        $cmd = $cmd . "-dset 2 2 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[18] \\\n";  # Neut incompat
                                        $cmd = $cmd . "-dset 3 2 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[30] \\\n";  # Open incompat
                                        $cmd = $cmd . "-dset 1 3 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[10] \\\n"; # Focused neut
                                        $cmd = $cmd . "-dset 2 3 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[22] \\\n";  # Neut neut
                                        $cmd = $cmd . "-dset 3 3 " . $subject . " ../$subjs[$subj]/Nstroop_med_deconv+tlrc[34] \\\n";  # Open neut
                                        ### should we add neutral stimuli as well??????
                                }
                        }
                }


                                 $cmd = $cmd ."-fa MedState \\\n";
                                 $cmd = $cmd ."-fb TrialType \\\n";
                                 $cmd = $cmd ."-fab MedTrialtypeInteract \\\n";

                $cmd = $cmd . "-bcontr -1 01 0 All_compatability \\\n";
                $cmd = $cmd . "-Abcontr 1 : -1 1 0 Foc_compatability \\\n";
                $cmd = $cmd . "-Abcontr 2 : -1 1 0 Neut_compatability \\\n";
                $cmd = $cmd . "-Abcontr 3 : -1 1 0 Open_compatability \\\n";

                $cmd = $cmd . "-amean 1 Focused \\\n";
                $cmd = $cmd . "-amean 2 Open \\\n";
                $cmd = $cmd . "-amean 3 Compassion \\\n";
                $cmd = $cmd . "-bmean 1 Congru \\\n";
                $cmd = $cmd . "-bmean 2 Incongru \\\n";
                $cmd = $cmd . "-bmean 3 Neutral \\\n";


                $cmd = $cmd . "-bucket ANOVA_NStroop_meditator > ANOVA_NStroop_meditator.out & \n";
                `rm -f ANOVA_NStroop_meditator*`;
                execute_afni_command($cmd);
        }        # end ThreeD_ANOVA_ERIKSEN == TRUE
        
        if (NSTROOP_GROUP == TRUE) {
                chdir($groupdir_name);
                $cmd = "3dttest -unpooled -prefix MedVsControl_LvsF_compat -set1 ";               
                	# load focused compatibility 
                $cmd = $cmd . "../med01/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med02/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med03/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med05/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med07/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med08/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med11/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med12/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med13/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med14/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med15/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../med16/Nstroop_med_deconv+tlrc'[78]' ";

				$cmd = $cmd . "-set2 ";
                
                $cmd = $cmd . "../cs01/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../cs02/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../cs07/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../cs08/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../cs10/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../cs14/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../cs17/Nstroop_med_deconv+tlrc'[78]' ";
                $cmd = $cmd . "../cs19/Nstroop_med_deconv+tlrc'[78]' ";
                
                execute_afni_command($cmd);


                $cmd = "3dttest -unpooled -prefix MedVsControl_LvsO_compat -set1 ";               
                	# load neutral-open compatibility 
                $cmd = $cmd . "../med01/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med02/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med03/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med05/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med07/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med08/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med11/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med12/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med13/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med14/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med15/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../med16/Nstroop_med_deconv+tlrc'[82]' ";

				$cmd = $cmd . "-set2 ";
                
                $cmd = $cmd . "../cs01/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../cs02/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../cs07/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../cs08/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../cs10/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../cs14/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../cs17/Nstroop_med_deconv+tlrc'[82]' ";
                $cmd = $cmd . "../cs19/Nstroop_med_deconv+tlrc'[82]' ";
                
  #              execute_afni_command($cmd);

                
                $cmd = "3dttest -unpooled -prefix MedVsControl_L_compat -set1 ";               
                	# load neutral compatability
                $cmd = $cmd . "../med01/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med02/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med03/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med05/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med07/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med08/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med11/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med12/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med13/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med14/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med15/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../med16/Nstroop_med_deconv+tlrc'[74]' ";

				$cmd = $cmd . "-set2 ";
                
                $cmd = $cmd . "../cs01/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../cs02/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../cs07/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../cs08/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../cs10/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../cs14/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../cs17/Nstroop_med_deconv+tlrc'[74]' ";
                $cmd = $cmd . "../cs19/Nstroop_med_deconv+tlrc'[74]' ";
                
     #           execute_afni_command($cmd);
                
                
                
        } # end execute_afni_command

                if (AE_PROBE_CONTRIB_TO_SPRITZ_IRESP == TRUE) {

                        chdir($afni_root_path."GroupAnalysis_main/");

                        if (FALSE) {
                                # Create waverd models for each regressor
                                $cmd = "waver -GAM -dt 2.0 -numout 225 -input 0-Probe_times_model.1D > 0-Probe_wavered_model.1D";
                                execute_afni_command($cmd);
                                $cmd = "waver -GAM -dt 2.0 -numout 225 -input 3-Probe_times_model.1D > 3-Probe_wavered_model.1D";
                                execute_afni_command($cmd);
                                $cmd = "waver -GAM -dt 2.0 -numout 225 -input 0_Juice_times_model.1D > 0-Juice_wavered_model.1D";
                                execute_afni_command($cmd);
                                $cmd = "waver -GAM -dt 2.0 -numout 225 -input 3_Juice_times_model.1D > 3-Juice_wavered_model.1D";
                                execute_afni_command($cmd);
                                $cmd = "waver -GAM -dt 2.0 -numout 225 -input 0_Quinine_times_model.1D > 0-Quinine_wavered_model.1D";
                                execute_afni_command($cmd);
                                $cmd = "waver -GAM -dt 2.0 -numout 225 -input 3_Quinine_times_model.1D > 3-Quinine_wavered_model.1D";
                                execute_afni_command($cmd);
                        }

                        if (FALSE) {
                                # Create3d_t dataset of model series times observed amplitude of HRS (Base to paek)
                                `rm -f *_iresp_peal+tlrc*`;
                                $cmd = "3dcalc -a iresp_undrgrd_0-Probe+tlrc[0] -b iresp_undrgrd_0-Probe+tlrc[1] -c iresp_undrgrd_0-Probe+tlrc[4] -d 0-Probe_wavered_model.1D -expr \"(c-((a+b)/2))*d\" -taxis 225:2 -prefix 0-Probe_iresp_peal"; # peak of Hemo response minus baseline average
                                execute_afni_command($cmd);
                                $cmd = "3dcalc -a iresp_undrgrd_3-Probe+tlrc[0] -b iresp_undrgrd_3-Probe+tlrc[1] -c iresp_undrgrd_3-Probe+tlrc[4] -d 3-Probe_wavered_model.1D -expr \"(c-((a+b)/2))*d\" -taxis 225:2 -prefix 3-Probe_iresp_peal"; # peak of Hemo response minus baseline average
                                execute_afni_command($cmd);

                                # convert these to the same size;
                                $cmd = "3dcalc -a iresp_undrgrd_Juice_0-Back+tlrc[4] -b 0-Juice_wavered_model.1D -expr \"a*b\" -taxis 225:2 -prefix 0-Juice_iresp_peal"; # peak of Hemo response minus baseline average
                                execute_afni_command($cmd);
                                $cmd = "3dcalc -a iresp_undrgrd_Juice_3-Back+tlrc[4] -b 3-Juice_wavered_model.1D -expr \"a*b\" -taxis 225:2 -prefix 3-Juice_iresp_peal"; # peak of Hemo response minus baseline average
                                execute_afni_command($cmd);
                                $cmd = "3dcalc -a iresp_undrgrd_Bitter_0-Back+tlrc[4] -b 0-Quinine_wavered_model.1D -expr \"a*b\" -taxis 225:2 -prefix 0-Quinine_iresp_peal"; # peak of Hemo response minus baseline average
                                execute_afni_command($cmd);
                                $cmd = "3dcalc -a iresp_undrgrd_Bitter_3-Back+tlrc[4] -b 3-Quinine_wavered_model.1D -expr \"a*b\" -taxis 225:2 -prefix 3-Quinine_iresp_peal"; # peak of Hemo response minus baseline average
                                execute_afni_command($cmd);
                        }

                        # Average, to create combined data series file

#                        `cp 0-Probe_iresp_peal 0-Probe_iresp_peal2`;
#                        `cp 3-Probe_iresp_peal 3-Probe_iresp_peal2`;

                        if (FALSE) {
                                `rm -f model_data*`;
                                $cmd = "3dMean -prefix model_data 0-Probe_iresp_peal+tlrc 0-Probe_iresp_peal+tlrc 0-Probe_iresp_peal+tlrc 0-Probe_iresp_peal+tlrc ";
                                $cmd = $cmd . "3-Probe_iresp_peal+tlrc 3-Probe_iresp_peal+tlrc 0-Juice_iresp_peal+tlrc 3-Juice_iresp_peal+tlrc 0-Quinine_iresp_peal+tlrc 3-Quinine_iresp_peal+tlrc  ";
                                execute_afni_command($cmd);
                        }


                        if (FALSE) {
                                #deonvolution

                                check_file_exists("probe_amplitude_model.HEAD");

                                $cmd = "3dDeconvolve -jobs 2 -input model_data+tlrc -mask group_mask+tlrc[0] -polort 0 -num_stimts 6 -tshift\\\n";
                                $cmd = $cmd . "-stim_file 1 0-Probe_times_shifted2_model.1D -stim_label 1 0-Probe -stim_maxlag 1 8 \\\n";
                                $cmd = $cmd . "-stim_file 2 3-Probe_times_shifted2_model.1D -stim_label 2 3-Probe -stim_maxlag 2 8 \\\n";
                                $cmd = $cmd . "-stim_file 3 0_Juice_times_shifted2_model.1D -stim_label 3 0-Juice -stim_maxlag 3 8 \\\n";
                                $cmd = $cmd . "-stim_file 4 3_Juice_times_shifted2_model.1D -stim_label 4 3-Juice -stim_maxlag 4 8 \\\n";
                                $cmd = $cmd . "-stim_file 5 0_Quinine_times_shifted2_model.1D -stim_label 5 0-Quinine -stim_maxlag 5 8 \\\n";
                                $cmd = $cmd . "-stim_file 6 3_Quinine_times_shifted2_model.1D -stim_label 6 3-Quinine -stim_maxlag 6 8 \\\n";
                                $cmd = $cmd . "-iresp 1 0-Probe \\\n";
                                $cmd = $cmd . "-iresp 2 3-Probe \\\n";
                                $cmd = $cmd . "-iresp 3 0-Juice \\\n";
                                $cmd = $cmd . "-iresp 4 3-Juice \\\n";
                                $cmd = $cmd . "-iresp 5 0-Quinine \\\n";
                                $cmd = $cmd . "-iresp 6 3-Quinine \\\n";
                                $cmd = $cmd . "-bucket probe_amplitude_model \\\n";


                                execute_afni_command($cmd);

                        }

                        if (TRUE) {
                                $cmd = "3dmaskave -mask Clusters/Spritz_Nback_FDR_26_2.17_clust_mask+tlrc -quiet -mrange 1 1 probe_amplitude_model+tlrc[19..27] > probe_amp_model_clust26_Juice-0 ";
                                execute_afni_command($cmd);
                                $cmd = "3dmaskave -mask Clusters/Spritz_Nback_FDR_26_2.17_clust_mask+tlrc -quiet -mrange 1 1 probe_amplitude_model+tlrc[28..36] > probe_amp_model_clust26_Juice-3 ";
                                execute_afni_command($cmd);
                                $cmd = "3dmaskave -mask Clusters/Spritz_Nback_FDR_26_2.17_clust_mask+tlrc -quiet -mrange 1 1 probe_amplitude_model+tlrc[37..45] > probe_amp_model_clust26_Quinine-0 ";
                                execute_afni_command($cmd);
                                $cmd = "3dmaskave -mask Clusters/Spritz_Nback_FDR_26_2.17_clust_mask+tlrc -quiet -mrange 1 1 probe_amplitude_model+tlrc[46..54] > probe_amp_model_clust_26_Quinine-3 ";
                                execute_afni_command($cmd);

                        }






                }  # AE_PROBE_CONTRIB_TO_SPRITZ_IRES


}  # End of Group Analysis



####################################################
###### Delete Preprocessing files , etc########
####################################################
if ($delete_intermediate_files) {
        print "Warning: deleting files used for intermediate processing. \n";
        print "If you are not done with preprocessing or group analysis, you may want these files.\n";
        print "Continue (Y/N)?\n";
        $input = <STDIN>;
        chomp $input;
        if ($input eq "Y") {
                for ($subj=0; $subj< $num_subjs; $subj++) {                                        # read
                        $cwd= $afni_root_path .  $subjs[$subj] . "/";
                        if(chdir($cwd) == 0) { log_file("Can't open directory $cwd",ERROR); }
                        print "Cleaning out preprocessing files for subject $subjs[$subj] \n";
                  #      $cmd = "rm -f all_norm*";
                  #      `$cmd`;
                        $cmd = "rm -f run?_vr*";
                        `$cmd`;
#                        $cmd = "rm -f t1-nudge?+*";
                        $cmd = "rm -f t1-nudge*";
                        `$cmd`;
                        $cmd = "rm -f t1-ss+*";
                        `$cmd`;
                        $cmd = "rm -f maskave_iresp*";
                        `$cmd`;
                        $cmd = "rm -f maskave_betas*";
                        `$cmd`;
                        $cmd = "rm -f run?_qual*";
                        `$cmd`;
                        $cmd = "rm -f qualscript_*";
                        `$cmd`;

                        $cmd = "rm -f *_residual+tlrc*";  # does not delete +orig residuals
                        `$cmd`;

                        # these can be removed
                        $cmd = "rm -f *t1-allin*";
                        `$cmd`;




                }
                if(chdir($groupdir_name) == 1) {
                        print "Cleaning out unnecessary files in GroupAnalysis directory \n";
                        $cmd = "rm -f *_within_*";
                        `$cmd`;

                }

                if(chdir($clusterdir_name) == 1) {
                        print "Cleaning out intermediate processing files for Cluster directory \n";
                        $cmd = "rm -f whereami_{*";
                        `$cmd`;
                        $cmd = "rm -f max_*";
                        `$cmd`;
                        $cmd = "rm -f GroupAverage_*";
                        `$cmd`;
                } else { log_file("Can't open directory $cwd",WARNING); }
        }
}  # End of Clean out preprocessing filed

                                        ############## Preprocess Group #############

if ($roi_analysis_group == TRUE) {
        chdir($scripts_dir);
        for ($subj=0; $subj< $num_subjs; $subj++) {                                        # Main loop (loops on subject)
                if ($preprocess_flag[$subj] ne "Y") { next; }                                                                        # skip to next subject if processing isnt explicitly turned on

                log_file("Starting $0 to calculate ROI info for $subjs[$subj]",LOG);
                $cmd = $program_filename." -r ".$subj;
                if (PARRALLELIZE) {
                        $cmd = "submit " . $cmd;
                        print "$cmd\n";
                        log_file($cmd,LOG);
                        system($cmd);
                } else {
                        $subject_return_value = execute_afni_command($cmd);
                }
                if ($subject_return_value == 1) { log_file("Problems processing subject $subjs[$subj]\n",WARNING); }
        }                         ##End subject loop - main loop ##
}   ### End preprocess Group


if ($roi_analysis_subject == TRUE) {

        ### RUN normal preprocessing processing, minus talairaching, through Waver. Then run this.
        ### Assumes ROIs have already been produced manualy in ROIs_epi+orig
        ### Bucket numbers of interest need to be manually specified (coef) for 3dROIstats command
        ##         add code if 3dfractionize is to be used right here
        ##  FIX SUPPORT FOR DESPIKING - currently will not processes despiked data, even though it looks like it

        $cwd= $afni_root_path .  $subjs[$subj] . "/";
        if(chdir($cwd) == 0) { log_file("Can't open directory $cwd",ERROR); }

        chdir ($cwd);
        opendir(DIR, $cwd);
        @runs = grep(/run.{1,2}\+orig\.HEAD/, readdir(DIR));
        closedir(DIR);

        foreach $r (@runs) {
                if ($r =~ /(run.*)\+orig.*HEAD/) {
                        $pfx = $1;
                        push(@pfxs, $pfx);
                }
        }

        if (USE_POST_VR_DESPIKE_IN_ROI_ANALYSIS) {
                $post_vr_despike_flag = "_ds2";
        }

        ### CREATE BINARY MASK FROM LAST RUN ###
        ### Stuff to improve mask -clfrac at levels 0.2 0.25, or 0.33; -dilate
        if (BINARY_MASK_ROI == ON) {
                check_file_exists("mask+orig.HEAD");
                $cmd = "3dAutomask -prefix mask " . $pfxs[$#pfxs] .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "+orig";
                execute_afni_command($cmd);
        }

        ### CALCULATE PERCENT SIGNAL CHANGE and MEAN INTENSITIES ###
        if (CALC_SIGNAL_CHANGE_ROI == ON) {
                foreach $p (@pfxs) {
                        check_file_exists($p."_mean+orig.HEAD");
                        $cmd = "3dTstat -prefix " . $p . "_mean " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "+orig";
                        execute_afni_command($cmd);
                        check_file_written($p."_mean+orig.HEAD");
                }

                        ### CALCULATE PCT SIGNAL CHANGE ###
                $allnorm = "";
                foreach $p (@pfxs) {
                        check_file_exists($p. $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "_norm+orig.HEAD");
                        $cmd = "3dcalc -a " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "+orig ";
                        $cmd = $cmd . "-b " . $p . "_mean+orig ";
                        $cmd = $cmd . "-expr \"(a/b*100)\" ";                                # normalize to 100% signal
                        #$cmd = $cmd . "-expr \"((a-b)/b*100)\" ";                        #(used in NStroop)
                        $cmd = $cmd . "-prefix " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "_norm";
                        execute_afni_command($cmd);
                        check_file_written($p. $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "_norm+orig.HEAD");

                        $allnorm = $allnorm . " " . $p .  $pre_vr_despike_flag . "_vr" . $post_vr_despike_flag . "_norm+orig";
                }

                ### CONCATENATE NORMALIZED DATA
                check_file_exists("all_norm_ROI+orig.HEAD");
                $cmd = "3dTcat -prefix all_norm_ROI " . $allnorm;
                execute_afni_command($cmd);
                check_file_written("all_norm_ROI+orig.HEAD");
        }  ### end of calculating percent signal change



        if (GENERAL_LINEAR_TEST_ROI == ON) {
                $num_input_timeseries = $regressor_count + NUM_MOTION_PARAMETERS;
                $deconv_bucket_name = $studyname."_deconv_ROI";

                $cmd = "3dDeconvolve -jobs 2 -nobout -input all_norm_ROI+orig -mask mask+orig -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries \\\n";

                for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                        $reg_file = $subjs[$subj]."_IdealWaveform_AllRuns".$regressor_root_fname[$regressor_num];
                        $cmd = $cmd . " -stim_file $regressor_num $reg_file -stim_label $regressor_num $regressor_indicator_char$regressor_name[$regressor_num] \\\n";
                }
                for ($regressor_num = $regressor_count +1 ; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                        $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - $regressor_count)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                }

                $cmd = $cmd . " -num_glt $num_contrasts \\\n";

                for ($contrast = 1; $contrast <= $num_contrasts; $contrast++) {
                        $cmd = $cmd . " -glt 1 Contrasts/$contrast_name[$contrast].1D -glt_label $contrast $contrast_name[$contrast] \\\n";
                }

                check_file_exists($deconv_bucket_name."+orig.HEAD");
                $cmd = $cmd . " -bucket $deconv_bucket_name -fout -tout -rout";
                execute_afni_command($cmd);
                check_file_written($deconv_bucket_name."+orig.HEAD");

        }                # End of GENERAL_LINEAR_TEST_ROI

        if (ESTIMATE_IMPULSE_RESPONSE_ROI == ON) {   ### This is currently outfitted to inlude several pre-event TRs

                $num_input_timeseries = $regressor_count + NUM_MOTION_PARAMETERS;
                $deconv_bucket_name = $studyname."_deconv_ROI";

                $cmd = "3dDeconvolve -jobs 2 -input all_norm_ROI+orig -mask mask+orig -polort $baseline_model_order[$subj] -concat Contrasts/startpoints.1D -num_stimts $num_input_timeseries -tsshift \\\n";

                for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                        $reg_file = $subjs[$subj]."_StimulusTiming_Shifted_AllRuns".$regressor_root_fname[$regressor_num];

                        #### Count number of events ###
                        open(REGFILE, $reg_file) or log_file("Can't open $reg_file. Must do so to count number of impulse", WARNING);
                        while($line = <REGFILE>) {
                                if($line =~ /^1/) {$impulse_count[$regressor_num]++;}
                        }
                        close(REGFILE);

                        $cmd = $cmd . " -stim_file $regressor_num $reg_file -stim_label $regressor_num $regressor_indicator_char$regressor_name[$regressor_num] -stim_maxlag $regressor_num ". sprintf("%d",$HRF_duration_TR + 2) ."\\\n";
                }

                for ($regressor_num = $regressor_count +1 ; $regressor_num <= $num_input_timeseries; $regressor_num++) {
                        $cmd = $cmd . " -stim_file $regressor_num 'vr_params.1D[".int($regressor_num - $regressor_count)."]' -stim_label $regressor_num movement".int($regressor_num - $regressor_count)."  -stim_base $regressor_num \\\n";
                }

                for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                        check_file_exists("iresp_ROI_$regressor_name[$regressor_num]+orig.HEAD");
                        check_file_exists("iresp_ROI_stdev_$regressor_name[$regressor_num]+orig.HEAD");
                        check_file_exists("iresp_ROI_$regressor_name[$regressor_num]+orig.BRIK");
                        check_file_exists("iresp_ROI_stdev_$regressor_name[$regressor_num]+orig.BRIK");

                        $cmd = $cmd . " -iresp $regressor_num iresp_ROI_$regressor_name[$regressor_num] -sresp $regressor_num iresp_ROI_stdev_$regressor_name[$regressor_num] \\\n";

                }

                $deconv_bucket_name = $studyname."_ROI_peri_stimulus_iresp_deconv";
                check_file_exists($deconv_bucket_name."+orig.HEAD");
                $cmd = $cmd . " -bucket $deconv_bucket_name";
                print $cmd;

                execute_afni_command($cmd);

                for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                        check_file_written("iresp_ROI_$regressor_name[$regressor_num]+orig.HEAD");
                        check_file_written("iresp_ROI_stdev_$regressor_name[$regressor_num]+orig.HEAD");
                        $cmd = "3dcalc -a iresp_ROI_stdev_$regressor_name[$regressor_num]+orig -expr \'a/sqrt($impulse_count[$regressor_num]) \' -prefix iresp_ROI_stder_$regressor_name[$regressor_num]";
                        execute_afni_command($cmd);

                        $num_of_ROIs = 2;   ## THIS needs to be placed in a more general location, or extracted directly from maximum value of ROI file
                        for ($ROI_num = 1; $ROI_num <= $num_of_ROIs; $ROI_num++) {
                                $cmd = "3dmaskave -mrange ".int($ROI_num-1)." $ROI_num -mask ROIs_epi+orig -q iresp_ROI_$regressor_name[$regressor_num]+orig > $subjs[$subj]_ROI".$ROI_num."_iresp_$regressor_name[$regressor_num].1D";
                                execute_afni_command($cmd);
                        }
                }

        }                  # end of ESTIMSTE_IMPULSE_RESPONSE_ROIR
}  # end of $roi_analysis_subject


                 #################################
                 ###### ROI group Analysis #######
                 #################################
if ($ROI_group_analysis == TRUE) {
        #@unique_groups = grep(!$group_hash{$_}++, @subject_group);

        if (GROUP_ROI_CLUSTER_STATS) {
                for ($subj=0; $subj< $num_subjs; $subj++) {                                        # read
                        $cwd= $afni_root_path .  $subjs[$subj] . "/";
                        if(chdir($cwd) == 0) { log_file("Can't open directory $cwd",ERROR); }

                        #        check_file_exists("Attention_deconv_ROI+orig.HEAD");
                        $cmd = "3dROIstats -mask ROIs_epi+orig -sigma -nzvoxels '".$studyname."_deconv_ROI+orig[0,4,8,12,16,20,24,28,32,36,40,44,48]' > $groupdir_name$subjs[$subj]_ROIstats.txt";
                        execute_afni_command($cmd);
                        $cmd = "3dclust -1clip 0.6 5.00 50 ROIs_epi+orig > $groupdir_name$subjs[$subj]_ROIclusterinfo.txt";
                        execute_afni_command($cmd);
                }
        } # end GROUP_ROI_CLUSTER_STATS


        if (GROUP_ROI_IRESP_ANALYSIS) {
                chdir($groupdir_name);
                for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                        if ($num_in_group[$group_num] == 1) {next;}                                                                                # Don't bother doing stats on group if just 1 subject
                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                $num_of_ROIs = 2;
                                for ($ROI_num = 1; $ROI_num <= $num_of_ROIs; $ROI_num++) {

                                        $first_in_group = TRUE;
                                        $iresp_ave_script_filename = "iresp_ave_script_$group_num\_$regressor_num\_$ROI_num.sh";
                                        check_file_exists("iresp_ROI_$unique_groups[$group_num]$regressor_name[$regressor_num]+tlrc.HEAD");
#                                        $f = "1deval ";
                                        $iresp_ave_cmd= "1deval ";
                                        $iresp_ave_cmd_fim = "3dcalc  -prefix iresp_ROI_$unique_groups[$group_num]_$regressor_name[$regressor_num] ";
                                        $concat_iresp_ave_cmd = "1dcat ";
                                        $mean_cmd = "-expr \'mean(";
                                        for ($subj=0; $subj< $num_subjs; $subj++) {
                                                if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group
                                                        $char = chr(97+$subj);
                                                        $iresp_ave_cmd = $iresp_ave_cmd . "-$char \'../$subjs[$subj]/$subjs[$subj]_ROI".$ROI_num."_iresp_$regressor_name[$regressor_num].1D\' ";
                                                        $iresp_ave_cmd_fim = $iresp_ave_cmd_fim . "-$char \'../$subjs[$subj]/iresp_ROI_$regressor_name[$regressor_num]+orig\' ";
                                                        if ($first_in_group == FALSE) { $mean_cmd = $mean_cmd . ",$char"; }
                                                        else { $mean_cmd = $mean_cmd . "$char"; }
                                                        $first_in_group = FALSE;
                                                        $concat_iresp_ave_cmd = $concat_iresp_ave_cmd . " ../$subjs[$subj]/$subjs[$subj]_ROI".$ROI_num."_iresp_$regressor_name[$regressor_num].1D ";
                                                }
                                        }
                                        $iresp_ave_cmd = $iresp_ave_cmd . $mean_cmd . ")\' > GroupAverage_$unique_groups[$group_num]_ROI".$ROI_num."_iresp_$regressor_name[$regressor_num].1D";
                                        $iresp_ave_cmd_fim = $iresp_ave_cmd_fim . $mean_cmd . ")\'";
                                        $concat_iresp_ave_cmd = $concat_iresp_ave_cmd . " > GroupConcat_$unique_groups[$group_num]_ROI".$ROI_num."_iresp_$regressor_name[$regressor_num].1D";

                                        if (PARRALLELIZE) {
                                                if(open(IRESP_AVE_FILE, ">".$iresp_ave_script_filename)) {
                                                        print IRESP_AVE_FILE $iresp_ave_cmd;
                                                        print IRESP_AVE_FILE $iresp_ave_cmd_fim;
                                                        print IRESP_AVE_FILE $concat_iresp_ave_cmd;
                                                } else {log_file("Can't open $iresp_ave_script_filename for writing",ERROR);}
                                                close(IRESP_AVE_FILE);
                                                `chmod 774 $iresp_ave_script_filename`;
                                                execute_afni_command("submit $groupdir_name\\$iresp_ave_script_filename");

                                        } else {

                                                execute_afni_command($iresp_ave_cmd);
                                                execute_afni_command($iresp_ave_cmd_fim);
                                                execute_afni_command($concat_iresp_ave_cmd);
                                        }
                                }        # end of loop on $ROI_num
                        }    # end for regressor_num
                }  # end for group_num
        } # end of GROUP_ROI_IRESP_ANALYSIS
}



if ($cluster_analysis_all == TRUE) {
        chdir($scripts_dir);
        for (my $source_num = 1; $source_num <= $clust_flag_num; $source_num++) {
                if ($clust_source_flag[$source_num] eq "Y") {
                        log_file("Starting $0 to process cluster ".$cluster_source_name[$source_num].$thresh[$source_num]."_".$source_num,LOG);
                        $cmd = $program_filename." -c ". $cluster_source_name[$source_num] ." ". $cluster_source_subbrik[$source_num] . " " . $thresh[$source_num]." 1 ".$source_num . " " .$min_cluster_size[$source_num];
                        if (PARRALLELIZE) {
                                my $file_name = "analyzecluster_".$cluster_source_name[$source_num].$thresh[$source_num]."_".$source_num.".sh";

                                if(open(ANAL_CLUST_FILE, ">".$file_name)) {
                                        print ANAL_CLUST_FILE $cmd . "\n";;
                                } else {log_file("Can't open $file_name for writing",ERROR);}
                                close(ANAL_CLUST_FILE);
                                `chmod 774 $file_name`;
                                execute_afni_command("submit $file_name");
                                sleep 5;   # need to space out launch - otherwise if same source data is used, creations of masks will interfere with each other

                        } else {
                                $subject_return_value = get_clusters($cluster_source_name[$source_num],$cluster_source_subbrik[$source_num],$thresh[$source_num],"1",$source_num,$min_cluster_size[$source_num]);
                        }
                        if ($subject_return_value == 1) { log_file("Problems processing subject $subjs[$subj]\n",WARNING); }

                }
        }
}  ## End of IF Cluster Analysis


if ($cluster_analysis_specific == TRUE) {
        if ($ARGV[0] eq "") {die("Must provide arguments in command line for specific cluster")}
        get_clusters($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5]);
}

if ($quality_analysis == TRUE) {
        write_data_quality_files();
}

if ($plot_quality_analysis == TRUE) {
        plot_data_quality_files();
}



##################################
########   Subfunctions   ########
##################################

sub skull_strip {
    check_file_exists("t1-ss+orig.HEAD");
    $cmd = "3dSkullStrip -input t1+orig -prefix t1-ss";
    execute_afni_command($cmd);
}

sub parse_info_file {

    my ($cluster_source_read,$num_in_group_counter);
    $cluster_source_read = FALSE;
    $num_in_group_counter = 0;

    opendir(DIR, $scripts_dir);
    @info_files = grep(/^[^\.].*info.txt/,readdir(DIR));        # find all *info.txt files, but not Mac ._ files


    if (scalar(@info_files) > 1) { log_file("More than one *info.txt file. Don't know which one to read",ERROR); }
    closedir(DIR);

    unless(open(INFOFILE, $scripts_dir . $info_files[0])) { log_file("Can't open info file: $!" , ERROR); }

    while ($line = <INFOFILE>) {
        if ($line =~ /^\s*Experimental_Design/) { $data_read = OFF; }
        if ($line =~ /^\s*Regressor_Filename_Root/g) {
                    $contrast_num = 0;
                    while ($line =~ m/(\S+)/g) {
                            $contrast_name[++$contrast_num] = $1;
                    }
                $regressor_read = ON;
         }
         else  {
                 if ($regressor_read == ON) {
                        if ($line =~ /^\s*$/) {   # blank line
                                $regressor_read = OFF;
                                next;
                        }
                        $line =~ m/(\S+)/g ;
                            $regressor_root_fname[++$regressor_count] = "_" .$1. ".1D";
                            $regressor_name[$regressor_count] = "$1";

                            $contrast_num = 0;
                            while ($line =~ m/(\S+)/g) {
                                    $contrast[++$contrast_num][$regressor_count] = $1;
                            }
                    }     # END if regressor read
        } # END else
        if ($line =~ /^\s*Source_StatMap_Filename_Root/g) {
                    $cluster_source_num = 0;
                    $clust_flag_num = 0;
                    while ($line =~ m/([^\s\+]+\+tlrc)\[?(\d+)?\]?/g) {

                            $cluster_source_name[++$cluster_source_num] = $1 ;

                                                        if (defined $2) {
                                                                        $cluster_source_subbrik[$cluster_source_num] = $2;
                                                        } else {
                                                                $cluster_source_subbrik[$cluster_source_num] = -1;
                                                                $subbrick = -1;
                                                        }
                    }
                $cluster_source_read = ON;
         }
         else  {
                 if ($cluster_source_read == ON) {
                        if ($line =~ /^\s*$/) {   # blank line
                                $cluster_source_read = OFF;
                                next;
                        }
                        if ($line =~ m/t-threshold/g) {
                                     while ($line =~ m/(\S+)/g) {
                                            $thresh[++$thresh_num] = $1;
                                    }
                        } else {
                                if ($line =~ m/Debridging_mask_root_filename\?/g) {
                                            while ($line =~ m/(\S+)/g) {
                                                    $debridge_mask_filename[++$bridge_mask_num] = $1 ;   # Period (.) is nothing
                                                    $debridge_mask_cm_filename[$bridge_mask_num] = $1 ;   # Period (.) is nothing
                                            }
                                    } else {
                                            if ($line =~ m/Min_cluster_size/g) {   ## minimum size in voxels, possibly already determined from Alphasim
                                                    while ($line =~ m/(\S+)/g) {
                                                            $min_cluster_size[++$min_cluster_index] = $1 ;
                                                    }
                                            } else {
                                                    if ($line =~ m/Analyze_clusters\?/g) {
                                                            while ($line =~ m/(\S+)/g) {
                                                                    $clust_source_flag[++$clust_flag_num] = uc($1);
                                                            }
                                                    } else {
                                                            $line =~ m/(\S+)/g ;
                                                            ++$clust_regressor_count;                  # may remove if above lines are re-added
                                                            $clust_num = 0;
                                                            while ($line =~ m/(\S+)/g) {
                                                                    $clust_regressor_flag[++$clust_num][$clust_regressor_count] = uc($1);
                                                            }
                                                    }
                                            }
                                    }
                            }
                    }     # END if regressor read
        }      # END else


        if ($line =~ /Notes/) { last; }
        if ($data_read == ON) {
                if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/g) {
                        $preprocess_flag[$num_subjs] = uc($1);
                        $subjs[$num_subjs] = $2;
                        $dicom_directory[$num_subjs] = $3;
                        $subject_group[$num_subjs] = $4;

                        $MPRAGE_series[$num_subjs] = "0" x ( 3 - length( $5 ) ) . $5;

                            $EPI_count[$num_subjs] = 0;
                            while ($line =~ m/(\S+)/g) {
                                    $EPI_series[$num_subjs][++$EPI_count[$num_subjs]] = "0" x ( 3 - length( $1 ) ) . $1;
                                $num_scans[$num_subjs][$EPI_count[$num_subjs]] = get_num_scans("\$EPI_series[".$num_subjs."][".$EPI_count[$num_subjs]."]",$dicom_root_path.$dicom_directory[$num_subjs]);
                            }
                        $num_subjs = $num_subjs + 1;
                }
        }
            if ($line =~ /Study_name\s+(\w*)/) {
                    if (JUST_REGRESS_OUT_BASELINE == ON) {
                        $studyname = $1 . "_JustBsln";  # for looking at effects of baseline
                } else {
                        $studyname = $1;
                }
                if ($studyname eq "") { log_file("No study name",WARNING); }
        }
        if ($line =~ /dicom_root_filelocation\s+(\S+)/) { $dicom_root_path = $1; }
        if ($line =~ /afnifiles_root_filelocation\s+(\S+)/) { $afni_root_path = $1; }
        if ($line =~ /spatial_smoothing_mm\s+(\S+)/) { $spatial_blur_diameter = $1; }
        if ($line =~ /baseline_model_order\s+(\S+)/) { $baseline_model_order[$num_subjs] = $1;}

        if ($line =~ /deconv_est_HRF_length_seconds\s+(\S+)/) {
                $HRF_duration_sec = $1;
                $HRF_duration_TR = $1 / $TR_duration;
        }
        if ($line =~ /deconv_est_pre-regressor-onset_HRF_length_seconds\s+(\S+)/) {$TRs_pre_regressor_for_estimating_baseline_time_courses = $1 / $TR_duration; }

        if ($line =~ /Preprocess_flag/) { $data_read = ON; }
    }

    close(INFOFILE);

    $num_contrasts = $contrast_num;
    @unique_groups = grep(!$group_hash{$_}++, @subject_group);
    for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
            $num_in_group_counter = 0;
            for ($subject_num = 0; $subject_num < $num_subjs; $subject_num++) {
                    if ($unique_groups[$group_num] eq $subject_group[$subject_num]) { $num_in_group_counter++; }
            }
            $num_in_group[$group_num] = $num_in_group_counter;
    }

    if ($num_subjs == 0) { log_file("No subjects read. Does info file exist (or is it currently being edited?",ERROR); }
    if ($dicom_root_path eq "") { log_file("No dicom root path specified (dicom_root_filelocation in col1, value in col2 of spreadsheet) ",ERROR); }
    if ($afni_root_path eq "") { log_file("No afni root path (afnifiles_root_filelocation in col1, value in col2 of spreadsheet)",ERROR); }
    if ($spatial_blur_diameter eq "") { log_file("No spatial blurring constant specified (spatial_smoothing_mm in col1, value in col2 of spreadsheet)",ERROR); }
    if ($baseline_model_order[0] eq "") {
            if (BASELINE_POLYNOMIAL_MODE == MANUAL) { log_file("No manual baseline model order specified, but code flag is set to manual baseline mode.\n(Spreadsheet:baseline_model_order in col1, value in col2 of spreadsheet)",ERROR); }
    }
}


sub get_max_nudge {  # used by tailaraching procedure
        opendir(DIR, ".");
        @nudgefiles = grep(/t1-nudge/, readdir(DIR));

        closedir(DIR);
        $max_nudge = 0;

        foreach $nudgefile (@nudgefiles) {
                $nudgefile =~ /t1\-nudge(\d+)\+orig/;
                if ($max_nudge < $1) { $max_nudge = $1 }
        }
}

sub autotalairach_structural_and_mask {
        #Talairach structural image
        check_file_exists("t1-aligned_at+tlrc.HEAD");
        check_file_exists("t1-aligned_at+tlrc.BRIK");

        #directory has moved
        $cmd = "\@auto_tlrc -no_ss -base /usr/csbmb/pkg/AFNI/afni/bindist/TT_icbm452+tlrc -input t1-aligned+orig";
        execute_afni_command($cmd);
        check_file_written("t1-aligned_at+tlrc.HEAD");

        #Talairach the mask
        check_file_exists("mask_at+tlrc.HEAD");
        check_file_exists("mask_at+tlrc.BRIK");

        $cmd = "\@auto_tlrc -apar t1-aligned_at+tlrc -input mask+orig -dxyz 3";
        execute_afni_command ($cmd);
        check_file_written("mask_at+tlrc.HEAD");

        $talrch_type = "at"
}

sub manualtalairach_structural_and_mask {
        #Talairach structural image
        check_file_exists("t1-aligned_mt+tlrc.HEAD");
        check_file_exists("t1-aligned_mt+tlrc.BRIK");

        $cmd = "adwarp -apar t1+tlrc. -dpar t1-aligned+orig -prefix t1-aligned_mt";
        execute_afni_command($cmd);
        check_file_written("t1-aligned_mt+tlrc.HEAD");

        #Talairach the mask
        check_file_exists("mask_mt+tlrc.HEAD");
        check_file_exists("mask_mt+tlrc.BRIK");
        $cmd = "adwarp -apar t1+tlrc. -dpar mask+orig -prefix mask_mt";
        execute_afni_command ($cmd);
        check_file_written("mask_mt+tlrc.HEAD");

        $talrch_type = "mt"
}

sub late_talairach { # executed after 3d_deconvolve
        # if TALAIRACH_MODE_LATE is set to OFF, noting is run
        my $deconv_output_filename = $_[0];

        if (TALAIRACH_MODE_LATE == AUTO) {
                get_max_nudge();

                autotalairach_structural_and_mask();

                #Talairach functional/EPI
                check_file_exists($deconv_output_filename."_at+tlrc.HEAD");
                check_file_exists($deconv_output_filename."_at+tlrc.BRIK");

                $cmd = "\@auto_tlrc -apar t1-aligned_at+tlrc -input ".$deconv_output_filename."+orig -dxyz 3";
                execute_afni_command($cmd);
                check_file_written($deconv_output_filename."_at+tlrc.HEAD");
        }

        if (TALAIRACH_MODE_LATE == MANUAL) {
                get_max_nudge();

                manualtalairach_structural_and_mask();

                #Talairach functional/EPI
                check_file_exists("all_norm_mt+tlrc.HEAD");
                check_file_exists("all_norm_mt+tlrc.BRIK");
                $cmd = "adwarp -apar t1+tlrc -dpar ".$deconv_output_filename."+orig -dxyz 3 -prefix all_norm_mt";
                execute_afni_command($cmd);
                check_file_written("all_norm_mt+tlrc.HEAD");
        }
}     # End late_talairachProcessing

sub plot_data_quality_files {
    my $run_num;
        for ($subj=0; $subj< $num_subjs; $subj++) {                                        # read
                $subj_dir = $afni_root_path .  $subjs[$subj] . "/";
                if(chdir($subj_dir) == 0) { log_file("Can't open directory $cwd",ERROR); }

                for ($run_num = 1; $run_num <= $EPI_count[$subj]; $run_num++) {
                        $cmd = "1dplot -xlabel $subjs[$subj] -ylabel qualityindex -sepscl -ynames run" . $run_num . "_stderr - run" . $run_num . "_qual.1d &" ;
                        print "$cmd\n";
                        `$cmd`;

                }
                close(DIR);
        } # END loop through subject
}  # End plot_data_quality_files()

sub write_data_quality_files {
    my ($cmd, $cmd2, $cmd3, $run_num, $subj, @censor_file);
    for ($subj=0; $subj< $num_subjs; $subj++) {                                        # read
    #        if ($preprocess_flag[$subj] ne "Y") { next; }                                                                        # skip to next subject if processing isnt explicitly turned on
        my $subj_dir= $afni_root_path .  $subjs[$subj] . "/";
        if(chdir($subj_dir) == 0) { log_file("Can't open directory $cwd",ERROR); }

        opendir(DIR, $subj_dir);

        for ($run_num = 1; $run_num <= $EPI_count[$subj]; $run_num++) {
                $cmd = "3dToutcount -automask -autoclip -polort 3 -range -automask run" . $run_num . "+orig > run".$run_num."_qual.1d";
                $cmd2 = "\@censor_file = grep(/^run".$run_num.$censor_fileroot."/,readdir(DIR))";
                eval $cmd2;
                if (@censor_file ne ""){
                        clean_files(@censor_file);
                        $clean_censor_file = "clean_run".$run_num.$censor_fileroot;
                } else {$clean_censor_file = ""};

                $cmd3 = "1dcat run".$run_num."_qual.1d $clean_censor_file > run".$run_num."_qual2.1d";

                if (PARRALLELIZE) {
                        $quality_script_filename = "qualscript_".$subj."_".$run_num.".sh";
                        if(open(QUALITY_FILE, ">".$quality_script_filename)) {
                                print QUALITY_FILE $cmd . "\n";
                                print QUALITY_FILE $cmd3 . "\n";
                        } else {log_file("Can't open $quality_script_filename for writing",ERROR);}
                        close(QUALITY_FILE);
                        `chmod 774 $quality_script_filename`;
                        execute_afni_command("submit $quality_script_filename");

                } else {
                        $subject_return_value = execute_afni_command($cmd);
                        $subject_return_value = execute_afni_command($cmd3);
                }

#                $cmd = "3dTqual -range -automask run" . $run_num . "_vr+orig > run".$run_num."_vr_qual";
                #print "$cmd\n";
                #`$cmd`;

        }
    }
}

sub initialize_3dDeconv_variables {
    my ($cmd, $clean_censor_files);

    if (JUST_REGRESS_OUT_BASELINE == OFF) {
    $num_input_timeseries = $regressor_count + NUM_MOTION_PARAMETERS;
    } else {
    $num_input_timeseries = 0;
    }
    # get censor files, if they exist. Must be named run<#>_censor.1d
    #Censor file is used by 3dDeconvolve to throw out bad sections of data
    opendir(DIR, $subj_dir);
    $cmd = "\@censor_files = grep(/^run\\d+".$censor_fileroot."/,readdir(DIR))";
    eval $cmd;
    if (scalar(@censor_files) == 0) {

    $CensoredVoxelsFileString = "";

    } else {
    if (scalar(@censor_files) != $EPI_count[$subj]) {
    log_file("Possibly incorect number (". scalar(@censor_files) .") of run censor files in subject $subjs[$subj]. Expecting 0 or $EPI_count[$subj].\nThere are $EPI_count[$subj] runs specified in *info.txt file \n",WARNING);
    }
    clean_files(@censor_files);
    $clean_censor_files = "clean_".join(' clean_' , @censor_files);

    $cmd = "cat $clean_censor_files > ". $subjs[$subj]."_censor_allruns.1d";                # make temporary file in row. May contain control/non-text characters is censor file was generated on Mac. 1dplot gets rid of that:
    execute_afni_command($cmd);
    $CensoredVoxelsFileString = "-censor " . $subjs[$subj]."_censor_allruns.1d";
    }
    closedir(DIR);
}

sub clean_files {
        # Clean AFNI 1d Files of spurious formatting characters as might occur with Mac files. Can have any number of rows or columns for input
        my @input_files = @_;
        my ($line, $file);

        if (scalar(@input_files) == 0) {log_file("No censor files passed to function clean_files for subject $subj",WARNING);}

        for ($file = 0; $file < scalar(@input_files); $file++) {
                my $input_file = $input_files[$file];

                unless(open(INPUTFILE, $input_file)) { die("Can't open input file: $!"); }
                unless(open(OUTPUTFILE, ">", "clean_".$input_file)) { die("Can't open output file: $!"); }

                while ($line = <INPUTFILE>) {
                        while ($line =~ m/(\d+)/g) {
                                print OUTPUTFILE $1 . "\n";
                        }
                }
        }
        close(INPUTFILE);
        close(OUTPUTFILE);
}

sub get_clusters {
        ### Input
        ### 1: source file name
        ### 2: subbrik
        ### 3: threshold of t-statistic
        ### 4: Nearest neighbor type (1= face, 2=edge, 3=corder)
        ### 5: Cluster source num (may want to remove or generalize this at some point)

        use Switch;

        # Global Variables
        use vars qw($Voxels $CM_RL $CM_AP $CM_IS $Mean_intens $Max_intens $Max_intens_RL $Max_intens_AP $Max_intens_IS);
        use vars qw($root_name $threshold $clust_name $loc $ave_beta_filename_clust $ave_iresp_filename_clust $clust_source_num $clust_num);

        # Local Variables
        my ($line, $cmd, $out_line, $cluster_masks_filename, $cluster_masks_CenterOfMass_filename,  $maskave_cmd);
        my ($iresp_ave_cmd, $iresp_SEM_btwnsbj_cmd, $beta_ave_cmd, $beta_SEM_btwnsbj_cmd, $first_regressor);
        my ($num_clusters, $subrick, $bucket_num, $parse_on, $char);
        my ($cluster_table_file_line, $whi_filename,  $clust_table_file_line);
        my $nn = "1.01";
        my $nearest_neighbor;
        my $subbrik = 0;

        $parse_on = FALSE;

        $Max_intens_RL = "";
        $Max_intens_AP = "";
        $Max_intens_IS = "";
        $Voxels = "";

        # parse function input
        my $source_filename = $_[0];
        $subbrik = $_[1];
        $threshold = $_[2];
        $nearest_neighbor = $_[3];
        $clust_source_num = $_[4];
        my $min_clust_size = $_[5];


        if ($source_filename =~ /(\S+)_(\S+)\+tlrc/) {  $root_name = $1 ."_"; } # root name, group and "+tlrc"

        if (int($subbrik) > -1) {
                $source_filename = "'../" . $source_filename ."[$subbrik]'";
        } else {
                $source_filename = "../" . $source_filename;
        }
        switch (int($nearest_neighbor)) {
                case 1 {$nn = "1.01"}
                case 2 {$nn = "1.42"}
                case 3 {$nn = "1.74"}
        }

        my $dilate_flag = "-1erode 50 -1dilate";
$dilate_flag = "";
        if (USE_ALPHASIM_FAMILYWISE == TRUE) {
                if ($min_clust_size < 6) { $dilate_flag = "" }  # can't see small clusters if dilation is turned on
        }
#        if (USE_ALPHASIM_FAMILYWISE == TRUE) {
#                if ($threshold eq "5.226") { $min_cluster_size = 6; }
#                if ($threshold eq "6.501") {
#                        $min_cluster_size = 3;
#                        $dilate_flag = "";  # can't see small clusters if dilation is turned on
#                }
#                if ($threshold eq "7.9") {
#                        $min_cluster_size = 1;
#                        $dilate_flag = "";
#                }
#        }

        chdir($clusterdir_name);

        open FOUT, ">Excel_summary_".$root_name.$threshold."_".$clust_source_num.".txt" or die "Can't write Excel summary output file\n";
        $out_line = "Group\tMask_orig_map\tRegressor\tMask_TVal\tCluster_Index\tCluster_Name\tMax_intens_RL\tMax_intens_AP\tMax_intens_IS\tVoxels\t";
        $out_line = $out_line."CM_RL\tCM_AP\tCM_IS\tLocation\t";
        $out_line = $out_line."beta\tiresp_TR1\t"."iresp_TR2\t"."iresp_TR3\t"."iresp_TR4\t"."iresp_TR5\t"."iresp_TR6\t"."iresp_TR7\t"."iresp_TR8\t"."iresp_TR9\t";
        $out_line = $out_line."betaSEM\tiresp_TR1_SEM\t"."iresp_TR2_SEM\t"."iresp_TR3_SEM\t"."iresp_TR4_SEM\t"."iresp_TR5_SEM\t"."iresp_TR6_SEM\t"."iresp_TR7_SEM\t"."iresp_TR8_SEM\t"."iresp_TR9_SEM\t";
        print FOUT $out_line."\n";


        if (WRITE_STAT_FILE == ON) {
                open STAT_FILE, ">Stats_".$root_name.$threshold."_".$clust_source_num.".txt" or die "Can't write stat file\n";
                $out_line = "Group\tSubject\tMask_orig_map\tRegressor\tMask_TVal\tCluster_Index\tCluster_Name\tMax_intens_RL\tMax_intens_AP\tMax_intens_IS\tVoxels\t";
                $out_line = $out_line."CM_RL\tCM_AP\tCM_IS\tLocation\t";
                $out_line = $out_line."beta\tiresp_TR1\t"."iresp_TR2\t"."iresp_TR3\t"."iresp_TR4\t"."iresp_TR5\t"."iresp_TR6\t"."iresp_TR7\t"."iresp_TR8\t"."iresp_TR9\t";
                print STAT_FILE $out_line."\n";
        }

        if ($debridge_mask_filename[$clust_source_num] eq ".") {   # If there are no bridging problems, if no mask created to break bridge

                $cluster_masks_filename = $root_name.$clust_source_num."_".$threshold."_clust_mask+tlrc";
                $cluster_masks_CenterOfMass_filename = $root_name.$clust_source_num."_".$threshold."_clust_mask_cm+tlrc";

                $cmd = "rm -f ". $cluster_masks_filename .".HEAD";
                execute_afni_command($cmd);
                $cmd = "rm -f ". $cluster_masks_CenterOfMass_filename .".HEAD";
                execute_afni_command($cmd);

                $cmd = "3dmerge $dilate_flag -dxyz=1 -1thresh $threshold -1clust_order $nn $min_clust_size -prefix $cluster_masks_filename $source_filename\n";
                execute_afni_command($cmd);

                $cmd = "3dclust -dxyz=1 $nn $min_clust_size $cluster_masks_filename > ".$root_name.$threshold."_clusters.txt";
                execute_afni_command($cmd);

                $cmd = "3dmerge $dilate_flag -dxyz=1 -1thresh $threshold -1clust $nn $min_clust_size -prefix $cluster_masks_CenterOfMass_filename $source_filename\n";
                execute_afni_command($cmd);

                $cmd = "3dclust -dxyz=1 $nn $min_clust_size $cluster_masks_CenterOfMass_filename > ".$root_name.$threshold."_cm_clusters.txt";
                execute_afni_command($cmd);

        } else { # if two areas are bridged, and a mask has been created to break the bridge
                $source_filename = $debridge_mask_filename[$clust_source_num]."_mask+tlrc";
                $cluster_masks_filename = $debridge_mask_filename[$clust_source_num]."_mask+tlrc";
                $cluster_masks_CenterOfMass_filename = $debridge_mask_cm_filename[$clust_source_num]."_mask_cm+tlrc";

                $cmd = "3dclust -dxyz=1 $nn $min_clust_size $cluster_masks_filename > ".$root_name.$threshold."_clusters.txt";
                execute_afni_command($cmd);

                $cmd = "3dclust -dxyz=1 $nn $min_clust_size $cluster_masks_CenterOfMass_filename > ".$root_name.$threshold."_cm_clusters.txt";
                execute_afni_command($cmd);
        }


                # Get number of clusters / masks
        $cmd = "3dBrickStat -max ".$cluster_masks_filename." > max_".$clust_source_num;
        execute_afni_command($cmd);

        unless(open(FIN, "max_".$clust_source_num)) { print STDERR "Can't open MAX file: $!\n"; }
        while ($line = <FIN>) {
                if ($line =~ /(\d*)/) { $num_clusters = $1; }
        }

        $maskave_cmd = "3dmaskave -mask ".$cluster_masks_filename;
        for ($group_num = 0; $group_num < scalar(@unique_groups); $group_num++) {
                if ($num_in_group[$group_num] == 1) {next;}
                for ($clust_num = 1; $clust_num <= $num_clusters; $clust_num ++) {

                        $first_regressor = TRUE;
                        $clust_name = $clust_num."-".$root_name.$subbrik."_".$threshold;

                        for ($regressor_num = 1; $regressor_num <= $regressor_count; $regressor_num++) {
                                if ($clust_regressor_flag[$clust_source_num][$regressor_num] eq "N") { next; }

                                $iresp_ave_cmd = "1deval ";
                                $iresp_SEM_btwnsbj_cmd = "1deval ";
                                $beta_ave_cmd = "1deval ";
                                $beta_SEM_btwnsbj_cmd = "1deval ";

                                $mean_cmd = "-expr \'mean(";
                                $sem_cmd = "-expr \'sem(";
                                $first_in_group = TRUE;

                                for ($subj = 0; $subj < $num_subjs; $subj ++) {

                                        if ($subject_group[$subj] eq $unique_groups[$group_num]) {                                                        # if member of current group

                                                $ave_iresp_filename_clust = "../../".$subjs[$subj]."/maskave_iresp_".$regressor_name[$regressor_num]."_".$root_name."_".$clust_source_num."_".$threshold."_clust".$clust_num.".1D";
                                                $ave_beta_filename_clust = "../../".$subjs[$subj]."/maskave_betas_".$regressor_name[$regressor_num]."_".$root_name."_".$clust_source_num."_".$threshold."_clust".$clust_num.".1D";

                                                $char = chr(97+$subj);

                                                $iresp_ave_cmd = $iresp_ave_cmd . "-$char \'".$ave_iresp_filename_clust."\' ";
                                                $iresp_SEM_btwnsbj_cmd = $iresp_SEM_btwnsbj_cmd . "-$char \'".$ave_iresp_filename_clust."\' ";
                                                $beta_ave_cmd = $beta_ave_cmd . "-$char \'".$ave_beta_filename_clust."\' ";
                                                $beta_SEM_btwnsbj_cmd = $beta_SEM_btwnsbj_cmd . "-$char \'".$ave_beta_filename_clust."\' ";

                                                if ($first_in_group == FALSE) {
                                                        $mean_cmd = $mean_cmd . ",$char";
                                                        $sem_cmd = $sem_cmd . ",$char";
                                                } else {
                                                        $mean_cmd = $mean_cmd . "$char";
                                                        $sem_cmd = $sem_cmd . "$char";

                                                        if ($first_regressor == TRUE) {  ## only want to do there these things for frist time we iterate within a cluster - not every time.

                                                                $bucket_num = get_AFNI_subbrik($regressor_indicator_char.$regressor_name[$regressor_num],"../../".$subjs[$subj]."/".$studyname."_deconv+tlrc");

                                                                        ##### Parse cluster table file #####
                                                                unless(open(F_CLUSTERS, $root_name.$threshold."_clusters.txt")) { print STDERR "Can't open cluster data file: $!\n"; }
                                                                while ($line = <F_CLUSTERS>) {
                                                                        if ($line =~ /-{2,}/) {                                                # Look for begining of cluster table
                                                                                $parse_on = TRUE;
                                                                        } else {
                                                                                if ($line =~ /(\d+)/) {          ###/[+-]?(\d+\.\d+|\d+\.|\.\d+|\d+)([eE][+-]?\d+)?/) # regular expression for any number
                                                                                        if ($parse_on == TRUE) {
                                                                                                $clust_table_file_line++;
                                                                                                if ($clust_table_file_line == $clust_num) {
                                                                                                        $line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/g;
                                                                                                        $Voxels = $1;
                                                                                                        $CM_RL = $2;
                                                                                                        $CM_AP = $3;
                                                                                                        $CM_IS = $4;
                                                                                                        $Mean_intens = $5;
                                                                                                        $Max_intens = $6;
                                                                                                        $Max_intens_RL = $7;
                                                                                                        $Max_intens_AP = $8;
                                                                                                        $Max_intens_IS = $9;
                                                                                                }
                                                                                        }  # END IF PARSE IS ON
                                                                                }# IF number
                                                                        } # Else
                                                                }        # End parse cluster table file
                                                                $parse_on = FALSE;
                                                                $clust_table_file_line = 0;

                                                                                                        ##### where am i ?  ####
                                                                $whi_filename = "whereami_{mask-".$root_name."-".$threshold."_clust".$clust_num."}.txt";
                                                                $cmd = "whereami ".$Max_intens_RL." ".$Max_intens_AP." ".$Max_intens_IS." > $whi_filename";
                                                                `$cmd`;
                                                                $loc = ".";
                                                                unless(open(F_WAI, $whi_filename)) { print STDERR "Can't open wherami file: $!\n"; }
                                                                while ($line = <F_WAI>) {
                                                                        if ($line =~ /Focus point: (.+)/) {
                                                                                        $loc= $1;
                                                                                        last;
                                                                        }
                                                                        if ($line =~ /Within (\d) mm: (.+)/) {
                                                                                        $loc= $2;
                                                                                        last;
                                                                        }
                                                                }        # End of Parse whereami file
                                                                close F_WAI;
                                                        }        # End of if first regressor        (only run for first subject and first regressor of each condition)
                                                }
                                                $first_in_group = FALSE;

                                                if (WRITE_INDIVIDUAL_MASKFILES == TRUE) {
                                                        #### Betas averaged #####
                                                        $cmd = $maskave_cmd . " -mrange ".$clust_num." ".$clust_num." -q ../../".$subjs[$subj]."/".$studyname."_deconv+tlrc[".$bucket_num."] >". $ave_beta_filename_clust;
                                                        execute_afni_command($cmd);

                                                        #### IRESP Means averaged ###
                                                        $cmd = $maskave_cmd . " -mrange ".$clust_num." ".$clust_num." -q ../../".$subjs[$subj]."/iresp_".$regressor_name[$regressor_num]."+tlrc > ". $ave_iresp_filename_clust;
                                                        execute_afni_command($cmd);
                                                } # End if write individual maskfile is turned on


                                                if (WRITE_STAT_FILE == ON) {
                                                        write_cluster_StatFile_line();
                                                        print "Subject $subj\n";
                                                }
                                        } # if member of current group
                                } # end loop through subjects


                                $iresp_ave_cmd = $iresp_ave_cmd . $mean_cmd . ")\' > GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_iresp_$regressor_name[$regressor_num]_mean.1D";
                                $iresp_SEM_btwnsbj_cmd = $iresp_SEM_btwnsbj_cmd . $sem_cmd. ")\' > GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_iresp_$regressor_name[$regressor_num]_SEM.1D";
                                $beta_ave_cmd = $beta_ave_cmd . $mean_cmd . ")\' > GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_betas_$regressor_name[$regressor_num]_mean.1D";
                                $beta_SEM_btwnsbj_cmd = $beta_SEM_btwnsbj_cmd . $sem_cmd. ")\' > GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_betas_$regressor_name[$regressor_num]_SEM.1D";

                                execute_afni_command($iresp_ave_cmd);
                                execute_afni_command($iresp_SEM_btwnsbj_cmd);  ## need to track down filenames
                                execute_afni_command($beta_ave_cmd);
                                execute_afni_command($beta_SEM_btwnsbj_cmd);  ## need to track down filenames

                                write_cluster_Excel_line();

                        } # END loop through primary regressors
                }    # loop through clusters
    } # end of loop through group num
} # End sub get_clusters

sub write_cluster_StatFile_line {
        # $0 = "mean" or "SEM"; i.e., calculation type
        #$computation_type = $_[0];

        my ($out_line, $line);

        $out_line = $unique_groups[$group_num]."\t".$subjs[$subj]."\t".$root_name."\t".$regressor_name[$regressor_num]."\t".$threshold."\t".$clust_num."\t".$clust_name."\t";
        $out_line = $out_line . $Max_intens_RL . "\t" . $Max_intens_AP. "\t" . $Max_intens_IS ."\t" . $Voxels . "\t" ;
        $out_line = $out_line . $CM_RL."\t".$CM_AP."\t".$CM_IS."\t".$loc."\t";

        #### Get BETA data ###
        unless(open(INFILE, $ave_beta_filename_clust)) { log_file("Can't open betas average file: $!" , ERROR); }
        while ($line = <INFILE>) {
                if ($line =~ /(\S+)/) {
                        $out_line = $out_line . $1 . "\t";
                }
        }

        #### Get IRESP data ###
        unless(open(INFILE, $ave_iresp_filename_clust)) { log_file("Can't open iresp average file: $!" , ERROR); }
        while ($line = <INFILE>) {
                if ($line =~ /(\S+)/) {
                        $out_line = $out_line . $1 . "\t";
                }
        }

        print STAT_FILE $out_line . "\n";  # FOUT is tab_delimted output file
        close INFILE;
}

sub write_cluster_Excel_line {
        # $0 = "mean" or "SEM"; i.e., calculation type
        my $computation_type = $_[0];
        my $line;

        my $out_line = $unique_groups[$group_num]."\t".$root_name."\t".$regressor_name[$regressor_num]."\t".$threshold."\t".$clust_num."\t".$clust_name."\t";
        $out_line = $out_line . $Max_intens_RL . "\t" . $Max_intens_AP. "\t" . $Max_intens_IS ."\t" . $Voxels . "\t" ;
        $out_line = $out_line . $CM_RL."\t".$CM_AP."\t".$CM_IS."\t". $loc."\t";

print "GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_betas_$regressor_name[$regressor_num]_mean.1D  \n";
        #### Get BETA data ###
        unless(open(IRESPFILE, "GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_betas_$regressor_name[$regressor_num]_mean.1D")) { log_file("Can't open betas average file: $!." , ERROR); }
        while ($line = <IRESPFILE>) {
                if ($line =~ /(\S+)/) {
                        $out_line = $out_line . $1 . "\t";
                }
        }

        #### Get IRESP data ###
        unless(open(IRESPFILE, "GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_iresp_$regressor_name[$regressor_num]_mean.1D")) { log_file("Can't open iresp average file: $!" , ERROR); }
        while ($line = <IRESPFILE>) {
                if ($line =~ /(\S+)/) {
                        $out_line = $out_line . $1 . "\t";
                }
        }

        #### Get BETA SEM data ###
        unless(open(IRESPFILE, "GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_betas_$regressor_name[$regressor_num]_SEM.1D")) { log_file("Can't open betas average file: $!" , ERROR); }
        while ($line = <IRESPFILE>) {
                if ($line =~ /(\S+)/) {
                        $out_line = $out_line . $1 . "\t";
                }
        }

        #### Get IRESP SEM data ###
        unless(open(IRESPFILE, "GroupAverage_$unique_groups[$group_num]_clust".$clust_num."_".$clust_source_num."_iresp_$regressor_name[$regressor_num]_SEM.1D")) { log_file("Can't open iresp average file: $!" , ERROR); }
        while ($line = <IRESPFILE>) {
                if ($line =~ /(\S+)/) {
                        $out_line = $out_line . $1 . "\t";
                }
        }

        print FOUT $out_line . "\n";  # FOUT is tab_delimted output file
        close IRESPFILE;
}

sub get_AFNI_subbrik {
        # $0 = regressor name;
        # $1 = bucket file (including path, if necesary)

        my $reg_name = $_[0];
        my $bucket_file = $_[1];

        my $brik_num = -1;
        my $line;

        execute_afni_command("3dinfo -verb $bucket_file > DeconvInfo.txt");

        unless(open(DECONV_INFO_FILE, "DeconvInfo.txt")) { log_file("Can't open info file: $!" , ERROR); }
        while ($line = <DECONV_INFO_FILE>) {
                if ($line =~ /(\d+) '$reg_name.*_[cC]oef'/) { $brik_num = $1;  }
        }
        close(DECONV_INFO_FILE);
        return $brik_num;
}

sub get_num_scans {   # get_num_scans(EPI_num, dicom_directory) returns the number of scans (PERL subs return last evaluated value if RETURN value not specified)
        opendir(DIR, $_[1]) || log_file("Cannot read directory $_[1]: $!",ERROR);
        my @runs;
        my $cmd = "\@runs = grep(/_". eval($_[0]) ."_/,readdir(DIR))";
        eval $cmd;
        closedir(DIR);
        scalar(@runs);
}

sub check_file_exists {
        my($filename) = @_;
        my $command;
        if(open(F, "< ".$filename)) {
                if (DELETE_PREVIOUS_FILES == ON) {
                        $command = "rm -f ". $filename;
                        print $command . "\n";
                        `$command` unless DEBUG_ONLY;
                }
                else {
                        log_file("File $filename already exists. May not be able to overwrite" , WARNING);
                }
        }
}

sub check_file_written {
        my($filename) = @_;
        open(F, "< ".$filename) or warn "Opening $filename: $!\n";
        unless (-C F < .0035) {                                                                # if file not written in the last .0035 days (5 minutes)
                log_file("$filename may not have been written (in past 5 minutes at least)." , ERROR);
                return -1;
        }
}

sub execute_afni_command {
        my $cmd = $_[0];
        print "$cmd\n";
        `$cmd` unless DEBUG_ONLY;
        log_file($cmd,LOG);
}

sub log_file {
    use Switch;
    my $msg_type = "";
    my $logstring = "";
    my $log_output_buffer = "";

    if (scalar(@_) > 1) {
            switch (int($_[1])) {
                    case LOG            {$msg_type = ""}
                    case WARNING        {$msg_type = "Warning: "}
                    case ERROR          {$msg_type = "Error: "}
            }
    }

    if ($log_file_open_flag == TRUE) {    # log file has succesfully been opened
            if ($log_output_buffer ne "") {
                    print LOGFILE $log_output_buffer;
                    $log_output_buffer = "";
            }
            $logstring = $msg_type . "$_[0]\n";
            print LOGFILE $logstring;
    }
    else {
            $log_output_buffer = $log_output_buffer .  $msg_type . $_[0] . "\n";
    }

    if ($_[1] == ERROR) { die $msg_type .$_[0]; }
}
