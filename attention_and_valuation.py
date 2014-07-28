__author__ = 'bfield'
import subprocess
import os.path
import sys
import csv
import itertools
import re
import numpy as np
from scipy.stats import ttest_rel

root_folder = "/jukebox/cohen/bfield/AttentionEnjoyment2anal4/AE_afni/"
group_folder = "/jukebox/cohen/bfield/AttentionEnjoyment2anal4/AE_afni/GroupAnalysis_PLoS/"
subjs = ["ae16","ae17","ae18","ae19","ae20","ae21","ae22","ae23","ae24","ae25","ae27","ae28","ae29","ae30","ae31","ae32"]
calc_flags = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p']
regressors = ['probe0','probe3','quinn0','quinn3','juice0','juice3','flick0','flick3']
regressor_type = ['beta','iresp']
marginal_means = ['spritz0_mean','spritz3_mean','probe_mean','juice_mean','quinn_mean','flick_mean','juice_diff','quin_diff']

iresp_briks = \
    {'probe0': '1,2,3,4,5',
     'probe3': '6,7,8,9,10',
     'quinn0': '11,12,13,14,15',
     'quinn3': '16,17,18,19,20',
     'juice0': '21,22,23,24,25',
     'juice3': '26,27,28,29,30',
     'flick0': '31,32,33,34,35',
     'flick3': '36,37,38,39,40'
    }
betas_brik =\
    {'probe0':'2',
     'probe3':'6',
     'quinn0':'10',
     'quinn3':'14',
     'juice0':'18',
     'juice3':'22',
     'flick0':'26',
     'flick3':'30',
     'load_Probe':'34',
     'load_Juice':'38',
     'load_Quin':'42',
     'load_Checkerboard':'54',
     'load_Salience':'46',
     'load_Valence':'50',
     'valence':'58'
    }    
masks = \
    {'main_stimulus': 'main_stimulus_q05_30vox_mask+tlrc',
     'main_load': 'main_load_q05_05vox_mask+tlrc',
     'interaction' : 'interact_q05_8vox_mask+tlrc',
     'load_probe' : 'cntrst_load_probe_mask+tlrc',
     'hedonic_v_visual' : 'cntrst_hedonic_vs_visual_mask+tlrc'
    }


def deconv_string(subj,brik_index_string,reg_type): 
    if reg_type == 'beta':  return root_folder + subj + '/' + 'AttentionEnjoy2_deconv+tlrc' + '\'[' + brik_index_string + ']\' '
    else:                   return root_folder + subj + '/' + 'AttentionEnjoy2_peri_stimulus_iresp_deconv+tlrc' + '\'[' + brik_index_string + ']\' ' 

def dataset_string(flag,brik_file,brik_index):
    briks = ""
    for subj in subjs:
        briks = briks + flag + ' ' + root_folder + subj + '/' + brik_file + '\'[' + str(brik_index) + ']\' '
    return briks    

def printr(string):         # print and return the same string
    print string
    return string 

def outfilepath_rm(outfile):
    if outfile != '':
        if outfile.find('/ae') == -1:   # -1 if group level, 0 if subject level
            return group_folder+outfile
        else:
            return outfile    
    return outfile 

def outfilepath(outfile):
    if outfile != '':
        if outfile.find('/ae') == -1:   # -1 if group level, 0 if subject level (which likely has the full path alread)
            return ' -prefix ' + group_folder + outfile
        else:
            return ' -prefix ' + outfile
    return outfile 

def afni(prog,infile,outfile,flags):
    if infile != '': infile = group_folder + infile + '+tlrc'
    if outfile !='':
        print('\nrm -f ' + outfilepath_rm(outfile) + '+tlrc.*')
        subprocess.call('rm -f ' + outfilepath_rm(outfile) + '+tlrc.*', shell=True)   
    print('\nafni_cmd: ' + prog + ' ' + outfilepath(outfile) + ' ' + flags + ' ' +infile + '\n')  
    subprocess.call(printr(prog + ' ' + outfilepath(outfile) + ' ' + flags + ' ' +infile),shell=True)

def add_3d(subj,brik1,brik2,outputname):  
	# beta
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_beta', ' -a '+deconv_string(subj, betas_brik[brik1],'beta') +' -b '+deconv_string(subj, betas_brik[brik2],'beta') +' -expr \'(a+b)\' ')
    #rint('debug: ' + '3drefit -sublabel 0 ' + outputname + '_beta ' + root_folder + subj + '/' + outputname + '+tlrc')
    subprocess.call('3drefit -sublabel 0 ' + outputname + '_beta ' + root_folder + subj + '/' + outputname + '_beta+tlrc', shell=True)    

    #iresp
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_iresp', ' -a '+deconv_string(subj, iresp_briks[brik1],'iresp') +' -b '+deconv_string(subj, iresp_briks[brik2],'iresp') +' -expr \'(a+b)\' ')    
    subprocess.call('3drefit -sublabel 0 TR1 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 1 TR2 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 2 TR3 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)          
    subprocess.call('3drefit -sublabel 3 TR4 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 4 TR5 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)        

def subtract_3d(subj,brik1,brik2,outputname):
	#beta
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_beta', ' -a '+deconv_string(subj, betas_brik[brik1],'beta') +' -b '+deconv_string(subj, betas_brik[brik2],'beta') +' -expr \'(a-b)\' ')    
    subprocess.call('3drefit -sublabel 0 ' + outputname + '_beta ' + root_folder + subj + '/' + outputname + '_beta+tlrc', shell=True)  

    #iresp
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_iresp', ' -a '+deconv_string(subj, iresp_briks[brik1],'iresp') +' -b '+deconv_string(subj, iresp_briks[brik2],'iresp') +' -expr \'(a-b)\' ')    
    subprocess.call('3drefit -sublabel 0 TR1 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 1 TR2 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 2 TR3 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)          
    subprocess.call('3drefit -sublabel 3 TR4 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 4 TR5 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)        

def interact_3d(subj,brik1,brik2,brik3,brik4,outputname): 
    # 2X2 iteraction: (brik1-brik2)-(brik3-brik4) = 0
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_beta', ' -a '+deconv_string(subj, betas_brik[brik1],'beta') +' -b '+ deconv_string(subj, betas_brik[brik2],'beta') + ' -c '+ deconv_string(subj, betas_brik[brik3],'beta') + ' -d '+ deconv_string(subj, betas_brik[brik4],'beta') +' -expr \'(a-b-c+d)\' ')    
    subprocess.call('3drefit -sublabel 0 ' + outputname + '_beta ' + root_folder + subj + '/' + outputname + '_beta+tlrc', shell=True)  
    
    #iresp
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_iresp', ' -a '+deconv_string(subj, betas_brik[brik1],'iresp') +' -b '+ deconv_string(subj, betas_brik[brik2],'iresp') + ' -c '+ deconv_string(subj, betas_brik[brik3],'iresp') + ' -d '+ deconv_string(subj, betas_brik[brik4],'iresp') +' -expr \'(a-b-c+d)\' ')    
    subprocess.call('3drefit -sublabel 0 TR1 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 1 TR2 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 2 TR3 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)          
    subprocess.call('3drefit -sublabel 3 TR4 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 4 TR5 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)   



def negative_3d(subj,inputname,outputname):
    #beta
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_beta', ' -a '+ root_folder + subj + '/' + inputname +'+tlrc -expr \'-a\' ')    
    subprocess.call('3drefit -sublabel 0 ' + outputname + '_beta ' + root_folder + subj + '/' + outputname + '_beta+tlrc', shell=True)   

    # #iresp
    # afni('3dcalc','',root_folder + subj + '/' + outputname + '_iresp', ' -a '+root_folder + subj + '/' + inputname_iresp +'+tlrc -expr \'-a\' ')    
    # subprocess.call('3drefit -sublabel 0 TR1 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    # subprocess.call('3drefit -sublabel 1 TR2 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    # subprocess.call('3drefit -sublabel 2 TR3 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)          
    # subprocess.call('3drefit -sublabel 3 TR4 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    # subprocess.call('3drefit -sublabel 4 TR5 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)          

def mean_3d(subj,brik1,brik2,outputname):
    # beta
    afni('3dcalc','',root_folder + subj + '/' + outputname + '_beta', ' -a '+deconv_string(subj, betas_brik[brik1],'beta') +' -b '+deconv_string(subj, betas_brik[brik2],'beta') +' -expr \'(a+b)/2\' ')    
    subprocess.call('3drefit -sublabel 0 ' + outputname + '_beta ' + root_folder + subj + '/' + outputname + '_beta+tlrc', shell=True)      

    # iresp
    afni('3dcalc','',root_folder + subj + '/' + outputname+'_iresp', ' -a '+deconv_string(subj, iresp_briks[brik1],'iresp') +' -b '+deconv_string(subj, iresp_briks[brik2],'iresp') +' -expr \'(a+b)/2\' ')    
    subprocess.call('3drefit -sublabel 0 TR1 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 1 TR2 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 2 TR3 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)          
    subprocess.call('3drefit -sublabel 3 TR4 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      
    subprocess.call('3drefit -sublabel 4 TR5 ' + root_folder + subj + '/' + outputname + '_iresp+tlrc', shell=True)      

def interaction_3d(subj,outputname):
    pass      

def make_cluster_mask(source_dataset,output_mask_dataset):
    subprocess.call('rm -f ' + outfilepath_rm(output_mask_dataset) + '+tlrc.*', shell=True)  
    afni('3dclust','','','-savemask '+group_folder+output_mask_dataset+' -quiet -1thresh 0.001 -dxyz=1 1.01 10 ' + group_folder+source_dataset + ' > '+group_folder+output_mask_dataset+'.1D') #'QuinineNegBold_masked_05+tlrc')  
  
def get_brik_indexes(regressor_type,regressor):
    if regressor_type == 'beta':     return betas_brik[regressor]
    elif regressor_type == 'iresp':  return iresp_briks[regressor]
    else: sys.stderr.write('unrecognized regressor type')

def subject_preprocessing_for_second_level_stats(): 
    for subj in subjs:
        #write marginal means
        mean_3d(subj,'juice0','quinn0','PLOS_spritz0_mean')
        mean_3d(subj,'juice3','quinn3','PLOS_spritz3_mean')
        mean_3d(subj,'probe0','probe3','PLOS_probe_mean')        
        mean_3d(subj,'juice0','juice3','PLOS_juice_mean')
        mean_3d(subj,'quinn0','quinn3','PLOS_quinn_mean')
        mean_3d(subj,'flick0','flick3','PLOS_flick_mean')

        afni('3dcalc','',root_folder + subj + '/PLOS_spritz_all_mean', ' -a '+ root_folder+subj+'/PLOS_spritz0_mean_beta+tlrc -b '+root_folder+subj+'/PLOS_spritz3_mean_beta+tlrc -expr \'(a+b)\'')

        # files for contrasts / ttests
        add_3d(subj,'quinn0','quinn3','PLOS_quinn_sum')      # valence
        add_3d(subj,'juice0','juice3','PLOS_juice_sum')      # valence
        add_3d(subj,'quinn0','juice0','PLOS_liquid_0')       # salience modulation
        add_3d(subj,'quinn3','juice3','PLOS_liquid_3')       # salience modulation       
        subtract_3d(subj,'quinn3','quinn0','PLOS_quin_diff') # Valence modulation
        subtract_3d(subj,'juice3','juice0','PLOS_juice_diff') # Valence modulation
        #subtract_3d(subj,'quinn0','quinn3','PLOS_inv_quin_diff') # Valence modulation       
        interact_3d(subj,'juice3','juice0','quinn3','quinn0','PLOS_liquid_load_interaction')    #iteraction
        #negative_3d(subj,'quinn_sum','PLOS_neg_quinn_sum')    


def count_file_lines(inputfile):
    with open(inputfile) as f:
        for i, l in enumerate(f):
            pass
        if 'i' in locals():
            return i+1
        else:
            return 0    

def write_mean_cluster(cluster_mask,subj_list = []): #,num_clusters):
    def cleanup():
        subprocess.call('rm -f ' +group_folder+'/Clusters/PLOS_' + cluster_mask + '*', shell=True)  
        for reg_type in regressor_type:   # beta versus iresp
            for subj in subjs:
                subprocess.call('rm -f ' +root_folder+subj+'/PLOS_*',shell-True) # + cluster_mask + '*', shell=True)   
    def filename_clustertable(): return group_folder+'Clusters/PLOS_'+cluster_mask+'_clustertable.1D'
    def filename_grandmean(reg,reg_type): return group_folder+'Clusters/PLOS_'+cluster_mask+'_'+reg+'_'+reg_type+'_grandmean.1D'
    def filename_grandstder(reg,reg_type): return group_folder+'Clusters/PLOS_'+cluster_mask+'_'+reg+'_'+reg_type+'_grandstder.1D'
    def write_afni_cluster_table():  # returns number of clusters
        subprocess.call(printr('rm -f ' + filename_clustertable()), shell=True) 
        afni('3dclust','','','-quiet -orient LPI -nosum -1thresh 0.001 -dxyz=1 1.01 10 ' + group_folder+cluster_mask + '+tlrc > ' + filename_clustertable() )
        return count_file_lines(filename_clustertable())
    def subj_cluster_mean_files(cluster_mask,regressor,cluster,reg_type,clust_num):
        files_string = ''
#        for subj_num in range(16): 
        for subj_num in range(len(subjs)):             
            files_string = files_string + '-' +calc_flags[subj_num]+' '+ root_folder + subjs[subj_num] + '/' +"PLOS_"+cluster_mask+'_'+regressor+'_mean_'+reg_type+'.1D[' +str(clust_num)+ '] '
        return files_string  
    def subj_cluster_marg_mean_files(cluster_mask,marg_mean,cluster,reg_type,clust_num):
        files_string = ''
        for subj_num in range(len(subjs)):        
            files_string = files_string + '-' +calc_flags[subj_num]+' '+ root_folder + subjs[subj_num] + '/' +"PLOS_"+cluster_mask+'_'+marg_mean+'_'+reg_type+'.1D[' +str(clust_num)+ '] '
        return files_string
       

    def write_means_by_regressor():
        def filename_clust_reg_stat(reg,cluster,reg_type,stat_type): return group_folder+'Clusters/PLOS_'+cluster_mask+'_'+reg+'_clust'+str(cluster)+'_'+reg_type+'_'+stat_type+'.1D '  
        for reg_type in regressor_type:   # ['beta','iresp']
            for subj in subjs:
                for reg in regressors:
                    subprocess.call(printr('3dROIstats -mask '+group_folder+cluster_mask+'+tlrc -quiet '+deconv_string(subj,get_brik_indexes(reg_type,reg),reg_type)+' > '+root_folder+subj+'/'+"PLOS_"+cluster_mask+'_'+reg+'_mean_'+reg_type+'.1D'),shell=True)   
                for marg_mean in marginal_means:
                    if reg_type == 'beta': 
                        subprocess.call(printr('3dROIstats -mask '+group_folder+cluster_mask+'+tlrc -quiet '+root_folder + subj + '/' + 'PLOS_'+marg_mean+'_'+reg_type+'+tlrc' + '\'[0]\' '+' > '+root_folder+subj+'/'+"PLOS_"+cluster_mask+'_'+marg_mean+'_'+reg_type+'.1D'),shell=True)                
                    if reg_type == 'iresp':     
                        subprocess.call(printr('3dROIstats -mask '+group_folder+cluster_mask+'+tlrc -quiet '+root_folder + subj + '/' + 'PLOS_'+marg_mean+'_'+reg_type+'+tlrc' + '\'[0..4]\' '+' > '+root_folder+subj+'/'+"PLOS_"+cluster_mask+'_'+marg_mean+'_'+reg_type+'.1D'),shell=True)                   
            for reg in regressors: 
                merge_mean_string = '1dcat ' 
                merge_stder_string = '1dcat '
                for cluster in range(num_clusters):  
                    subprocess.call(printr('1deval ' + subj_cluster_mean_files(cluster_mask,reg,cluster,reg_type,cluster) + '-expr \'mean(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)\' > ' + filename_clust_reg_stat(reg,cluster,reg_type,'mean') ), shell=True)
                    merge_mean_string = merge_mean_string + filename_clust_reg_stat(reg,cluster,reg_type,'mean')            
                    subprocess.call(printr('1deval ' + subj_cluster_mean_files(cluster_mask,reg,cluster,reg_type,cluster) + '-expr \'sem(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)\' > ' + filename_clust_reg_stat(reg,cluster,reg_type,'stder')), shell=True)  
                    merge_stder_string = merge_stder_string + filename_clust_reg_stat(reg,cluster,reg_type,'stder')                         
                subprocess.call(printr(merge_mean_string + ' > ' + filename_grandmean(reg,reg_type)), shell=True)        # pull togeter stats into one file to make life easier later                  
                subprocess.call(printr(merge_stder_string + ' > ' + filename_grandstder(reg,reg_type)), shell=True)        # pull togeter stats into one file to make life easier later    
            for marg_mean in marginal_means:    # we basically treat the marginal means as regressors from the convolution analys
                merge_mean_string = '1dcat '
                merge_stder_string = '1dcat '
                for cluster in range(num_clusters):
                	subprocess.call(printr('1deval ' + subj_cluster_marg_mean_files(cluster_mask,marg_mean,cluster,reg_type,cluster) + '-expr \'mean(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)\' > ' + filename_clust_reg_stat(marg_mean,cluster,reg_type,'mean') ), shell=True)
                	merge_mean_string = merge_mean_string + filename_clust_reg_stat(marg_mean,cluster,reg_type,'mean')            
                	subprocess.call(printr('1deval ' + subj_cluster_marg_mean_files(cluster_mask,marg_mean,cluster,reg_type,cluster) + '-expr \'sem(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)\' > ' + filename_clust_reg_stat(marg_mean,cluster,reg_type,'stder')), shell=True)  
                	merge_stder_string = merge_stder_string + filename_clust_reg_stat(marg_mean,cluster,reg_type,'stder')                         
                subprocess.call(printr(merge_mean_string + ' > ' + filename_grandmean(marg_mean,reg_type)), shell=True)        # pull togeter stats into one file to make life easier later                  
                subprocess.call(printr(merge_stder_string + ' > ' + filename_grandstder(marg_mean,reg_type)), shell=True)        # pull togeter stats into one file to make life easier later 
 

    def make_excel_file():
        def header_row_clustXreg(): return ['Cluster_mask']+['Regressor']+['Cluster_num']+['Volume_mm']+['Voxels']+['CM_RL']+['CM_AP']+['CM_IS']+['MI_RL']+['MI_AP']+['MI_IS']+['Beta mean']+['Beta SEM']+['IRESP Mean 1']+['IRESP Mean 2']+['IRESP Mean 3']+['IRESP Mean 4']+['IRESP Mean 5']+['IRESP stder 1']+['IRESP stder 2']+['IRESP stder 3']+['IRESP stder 4']+['IRESP stder 5']
        def header_row_clust(): return ['Cluster_mask']+['Cluster_num']+['Volume_mm']+['Voxels']+['MI_RL']+['MI_AP']+['MI_IS']+\
            ['Probe0 beta']+['Probe0 beta SEM']+['Probe0 IRESP1']+['Probe0 IRESP2']+['Probe0 IRESP3']+['Probe0 IRESP4']+['Probe0 IRESP5']+['Probe0 IRESP SEM1']+['Probe0 IRESP SEM2']+['Probe0 IRESP SEM3']+['Probe0 IRESP SEM4']+['Probe0 IRESP SEM5']+\
            ['Probe3 beta']+['Probe3 beta SEM']+['Probe3 IRESP1']+['Probe3 IRESP2']+['Probe3 IRESP3']+['Probe3 IRESP4']+['Probe3 IRESP5']+['Probe3 IRESP SEM1']+['Probe3 IRESP SEM2']+['Probe3 IRESP SEM3']+['Probe3 IRESP SEM4']+['Probe3 IRESP SEM5']+\
            ['Quinn0 beta']+['Quinn0 beta SEM']+['Quinn0 IRESP1']+['Quinn0 IRESP2']+['Quinn0 IRESP3']+['Quinn0 IRESP4']+['Quinn0 IRESP5']+['Quinn0 IRESP SEM1']+['Quinn0 IRESP SEM2']+['Quinn0 IRESP SEM3']+['Quinn0 IRESP SEM4']+['Quinn0 IRESP SEM5']+\
            ['Quinn3 beta']+['Quinn3 beta SEM']+['Quinn3 IRESP1']+['Quinn3 IRESP2']+['Quinn3 IRESP3']+['Quinn3 IRESP4']+['Quinn3 IRESP5']+['Quinn3 IRESP SEM1']+['Quinn3 IRESP SEM2']+['Quinn3 IRESP SEM3']+['Quinn3 IRESP SEM4']+['Quinn3 IRESP SEM5']+\
            ['Juice0 beta']+['Juice0 beta SEM']+['Juice0 IRESP1']+['Juice0 IRESP2']+['Juice0 IRESP3']+['Juice0 IRESP4']+['Juice0 IRESP5']+['Juice0 IRESP SEM1']+['Juice0 IRESP SEM2']+['Juice0 IRESP SEM3']+['Juice0 IRESP SEM4']+['Juice0 IRESP SEM5']+\
            ['Juice3 beta']+['Juice3 beta SEM']+['Juice3 IRESP1']+['Juice3 IRESP2']+['Juice3 IRESP3']+['Juice3 IRESP4']+['Juice3 IRESP5']+['Juice3 IRESP SEM1']+['Juice3 IRESP SEM2']+['Juice3 IRESP SEM3']+['Juice3 IRESP SEM4']+['Juice3 IRESP SEM5']+\
            ['Flick0 beta']+['Flick0 beta SEM']+['Flick0 IRESP1']+['Flick0 IRESP2']+['Flick0 IRESP3']+['Flick0 IRESP4']+['Flick0 IRESP5']+['Flick0 IRESP SEM1']+['Flick0 IRESP SEM2']+['Flick0 IRESP SEM3']+['Flick0 IRESP SEM4']+['Flick0 IRESP SEM5']+\
            ['Flick3 beta']+['Flick3 beta SEM']+['Flick3 IRESP1']+['Flick3 IRESP2']+['Flick3 IRESP3']+['Flick3 IRESP4']+['Flick3 IRESP5']+['Flick3 IRESP SEM1']+['Flick3 IRESP SEM2']+['Flick3 IRESP SEM3']+['Flick3 IRESP SEM4']+['Flick3 IRESP SEM5']+\
            ['spritz0_mean beta']+['spritz0_mean beta SEM']+['spritz0_mean IRESP1']+['spritz0_mean IRESP2']+['spritz0_mean IRESP3']+['spritz0_mean IRESP4']+['spritz0_mean IRESP5']+['spritz0_mean IRESP SEM1']+['spritz0_mean IRESP SEM2']+['spritz0_mean IRESP SEM3']+['spritz0_mean IRESP SEM4']+['spritz0_mean IRESP SEM5']+\
            ['spritz3_mean beta']+['spritz3_mean beta SEM']+['spritz3_mean IRESP1']+['spritz3_mean IRESP2']+['spritz3_mean IRESP3']+['spritz3_mean IRESP4']+['spritz3_mean IRESP5']+['spritz3_mean IRESP SEM1']+['spritz3_mean IRESP SEM2']+['spritz3_mean IRESP SEM3']+['spritz3_mean IRESP SEM4']+['spritz3_mean IRESP SEM5']+\
            ['probe_mean beta']+['probe_mean beta SEM']+['probe_mean IRESP1']+['probe_mean IRESP2']+['probe_mean IRESP3']+['probe_mean IRESP4']+['probe_mean IRESP5']+['probe_mean IRESP SEM1']+['probe_mean IRESP SEM2']+['probe_mean IRESP SEM3']+['probe_mean IRESP SEM4']+['probe_mean IRESP SEM5']+\
            ['juice_mean beta']+['juice_mean beta SEM']+['juice_mean IRESP1']+['juice_mean IRESP2']+['juice_mean IRESP3']+['juice_mean IRESP4']+['juice_mean IRESP5']+['juice_mean IRESP SEM1']+['juice_mean IRESP SEM2']+['juice_mean IRESP SEM3']+['juice_mean IRESP SEM4']+['juice_mean IRESP SEM5']+\
            ['quinn_mean beta']+['quinn_mean beta SEM']+['quinn_mean IRESP1']+['quinn_mean IRESP2']+['quinn_mean IRESP3']+['quinn_mean IRESP4']+['quinn_mean IRESP5']+['quinn_mean IRESP SEM1']+['quinn_mean IRESP SEM2']+['quinn_mean IRESP SEM3']+['quinn_mean IRESP SEM4']+['quinn_mean IRESP SEM5']+\
            ['flick_mean beta']+['flick_mean beta SEM']+['flick_mean IRESP1']+['flick_mean IRESP2']+['flick_mean IRESP3']+['flick_mean IRESP4']+['flick_mean IRESP5']+['flick_mean IRESP SEM1']+['flick_mean IRESP SEM2']+['flick_mean IRESP SEM3']+['flick_mean IRESP SEM4']+['flick_mean IRESP SEM5']+\
            ['juice_diff beta']+['juice_diff beta SEM']+['juice_diff IRESP1']+['juice_diff IRESP2']+['juice_diff IRESP3']+['juice_diff IRESP4']+['juice_diff IRESP5']+['juice_diff IRESP SEM1']+['juice_diff IRESP SEM2']+['juice_diff IRESP SEM3']+['juice_diff IRESP SEM4']+['juice_diff IRESP SEM5']+\
            ['quin_diff beta']+['quin_diff beta SEM']+['quin_diff IRESP1']+['quin_diff IRESP2']+['quin_diff IRESP3']+['quin_diff IRESP4']+['quin_diff IRESP5']+['quin_diff IRESP SEM1']+['quin_diff IRESP SEM2']+['quin_diff IRESP SEM3']+['quin_diff IRESP SEM4']+['quin_diff IRESP SEM5']+\
            ['juice ttest']+['quinn ttest']+['spritz ttest']+['probe ttest']+['flick ttest']+['valence ttest']+['juice_thresh']+['quinn_thresh']+['spritz_thresh']+['probe_thresh']+['flick_thresh']+['valence_thresh']+['CM_RL']+['CM_AP']+['CM_IS']+['ValenceMod_ttest']+['ValenceMod_thresh']
        def out_row(): return [cluster_mask] + [reg] + [str(cluster_num)] + [volume] + [cluster[0],cluster[1],cluster[2],cluster[3],cluster[13],cluster[14],cluster[15]] + list(mean_beta_row) + list(stder_beta_row) + list(mean_iresp_row) + list(stder_iresp_row)

        def subj_cluster_mean_file(subj,reg): return root_folder + subj + '/PLOS_'+cluster_mask+'_'+reg+'_beta.1D' 

        def betas(betafile): return betafile.readline().strip().split() 

        def calculate_ttests():
            def holm_threshold(prob_list):
                def thresh(p_val_num): return (0.05 / (p_val_num + 1))
                prob_list.sort()
                for p_val_num in range(len(prob_list)):
                    if prob_list[p_val_num] > thresh(p_val_num): return thresh(p_val_num)   # holm method says first time ttest pvals are larger than their associated thresh, then thresh is adjusted threshold for significance
                if 'p_val_num' in locals():
                    return thresh(p_val_num)  # return last threshold if no pvals never larger than threshold 
                else: return 0    

            juice_ttest = list(); quinn_ttest = list(); spritz_ttest = list(); valence_ttest = list(); probe_ttest = list(); flick_ttest = list(); valencemod_ttest = list();
            for cluster in range(num_clusters): 
                probe0 = list(); probe3 = list(); juice0 = list(); juice3 = list(); quinn0 = list(); quinn3 = list(); flick0 = list(); flick3 = list(); spritz0 = list(); spritz3 = list(); quinn = list(); juice = list(); juice_diff = list(); quin_diff = list();
                for subj in subjs:
                	with open(subj_cluster_mean_file(subj,'probe0_mean')) 	       as subj_probe0_data:   				probe0.append(float(betas(subj_probe0_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'probe3_mean')) 	       as subj_probe3_data:   				probe3.append(float(betas(subj_probe3_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'juice0_mean')) 	       as subj_juice0_data:   				juice0.append(float(betas(subj_juice0_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'juice3_mean')) 	       as subj_juice3_data:   				juice3.append(float(betas(subj_juice3_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'quinn0_mean')) 	       as subj_quinn0_data:   				quinn0.append(float(betas(subj_quinn0_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'quinn3_mean')) 	       as subj_quinn3_data:   				quinn3.append(float(betas(subj_quinn3_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'flick0_mean')) 	       as subj_flick0_data:   				flick0.append(float(betas(subj_flick0_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'flick3_mean')) 	       as subj_flick3_data:   				flick3.append(float(betas(subj_flick3_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'spritz0_mean'))         as subj_quinn0_data:  				spritz0.append(float(betas(subj_quinn0_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'spritz3_mean'))         as subj_quinn3_data:  				spritz3.append(float(betas(subj_quinn3_data)[cluster]))
                 	with open(subj_cluster_mean_file(subj,'quinn_mean'))           as subj_quinn_data:         			quinn.append(float(betas(subj_quinn_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'juice_mean')) 	       as subj_juice_data:       			juice.append(float(betas(subj_juice_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'juice_diff')) 	       as subj_juice_diff_data:     		juice_diff.append(float(betas(subj_juice_diff_data)[cluster]))
                	with open(subj_cluster_mean_file(subj,'quin_diff'))            as subj_quin_diff_data:   		    quin_diff.append(float(betas(subj_quin_diff_data)[cluster]))   

                juice_ttest.append(float(ttest_rel(juice0,juice3)[1]))
                quinn_ttest.append(float(ttest_rel(quinn0,quinn3)[1]))
                spritz_ttest.append(float(ttest_rel(spritz0,spritz3)[1]))
                probe_ttest.append(float(ttest_rel(probe0,probe3)[1]))
                flick_ttest.append(float(ttest_rel(flick0,flick3)[1]))
                valence_ttest.append(float(ttest_rel(quinn,juice)[1]))
                valencemod_ttest.append(float(ttest_rel(juice_diff,quin_diff)[1]))  # bug here

            return (juice_ttest, quinn_ttest, spritz_ttest, probe_ttest, flick_ttest, valence_ttest, valencemod_ttest, holm_threshold(juice_ttest), holm_threshold(quinn_ttest), holm_threshold(spritz_ttest), holm_threshold(probe_ttest), holm_threshold(flick_ttest), holm_threshold(valence_ttest), holm_threshold(valencemod_ttest))  



        excel_file = open(group_folder+'Clusters/'+cluster_mask+'table_by_cluster_and_regressor.csv','w')
        excel_csv = csv.writer(excel_file, delimiter=',') 
        excel_csv.writerow(header_row_clustXreg())  
        for reg in regressors + marginal_means: # regressors
            with open(filename_clustertable()) as clustertable_file:
                cluster_table = csv.reader(clustertable_file, delimiter=' ', skipinitialspace=True)
                with open(filename_grandmean(reg,'beta')) as beta_means_file:
                    with open(filename_grandstder(reg,'beta')) as beta_stder_file:
                        with open(filename_grandmean(reg,'iresp')) as iresp_means_file:
                            with open(filename_grandstder(reg,'iresp')) as iresp_stder_file:
                                rows_beta_mean = (row.strip().split() for row in beta_means_file)
                                transposed_beta_means = zip(*(row for row in rows_beta_mean if row))  
                                rows_beta_stder = (row.strip().split() for row in beta_stder_file)
                                transposed_beta_stder = zip(*(row for row in rows_beta_stder if row)) 
                                rows_iresp_mean = (row.strip().split() for row in iresp_means_file)
                                transposed_iresp_means = zip(*(row for row in rows_iresp_mean if row))  
                                rows_iresp_stder = (row.strip().split() for row in iresp_stder_file)
                                transposed_iresp_stder = zip(*(row for row in rows_iresp_stder if row)) 
                                cluster_num = 0
                                for mean_beta_row,stder_beta_row,mean_iresp_row,stder_iresp_row in itertools.izip(transposed_beta_means,transposed_beta_stder,transposed_iresp_means,transposed_iresp_stder): # go through each file, row by row, at the same time
                                    cluster_num += 1                
                                    cluster = clustertable_file.readline().strip().split()
                                    volume = int(cluster[0]) * 27
                                    excel_csv.writerow(out_row()) 
        excel_file.close()         


        juice_ttest, quinn_ttest, spritz_ttest, probe_ttest, flick_ttest, valence_ttest, valencemod_ttest, juice_thresh, quinn_thresh, spritz_thresh, probe_thresh, flick_thresh, valence_thresh, valencemod_thresh = calculate_ttests()        
        excel_file = open(group_folder+'Clusters/'+cluster_mask+'table_by_cluster' + reduced_filename_flag+ '.csv','w')
        excel_csv = csv.writer(excel_file, delimiter=',') 
        excel_csv.writerow(header_row_clust())
        cluster_num = 0      

        with open(filename_clustertable()) as clustertable_file:
            cluster_table = csv.reader(clustertable_file, delimiter=' ', skipinitialspace=True)
            for line in clustertable_file:
                cluster = line.strip().split()
                volume = int(cluster[0]) * 27
                cluster_num += 1                         
                outline = [cluster_mask] + [str(cluster_num)] + [volume] + [cluster[0],cluster[13],cluster[14],cluster[15]]  # add column for location                
                for reg in regressors + marginal_means: # regressors

                    with open(filename_grandmean(reg,'beta')) as beta_means_file:
                        with open(filename_grandstder(reg,'beta')) as beta_stder_file:
                            with open(filename_grandmean(reg,'iresp')) as iresp_means_file:
                                with open(filename_grandstder(reg,'iresp')) as iresp_stder_file:
                                    rows_beta_mean = (row.strip().split() for row in beta_means_file)
                                    transposed_beta_means = zip(*(row for row in rows_beta_mean if row))  
                                    rows_beta_stder = (row.strip().split() for row in beta_stder_file)
                                    transposed_beta_stder = zip(*(row for row in rows_beta_stder if row)) 
                                    rows_iresp_mean = (row.strip().split() for row in iresp_means_file)
                                    transposed_iresp_means = zip(*(row for row in rows_iresp_mean if row))  
                                    rows_iresp_stder = (row.strip().split() for row in iresp_stder_file)
                                    transposed_iresp_stder = zip(*(row for row in rows_iresp_stder if row)) 
                                    outline += list(transposed_beta_means[cluster_num-1]) + list(transposed_beta_stder[cluster_num-1]) + list(transposed_iresp_means[cluster_num-1]) + list(transposed_iresp_stder[cluster_num-1])
                excel_csv.writerow(outline+[juice_ttest[cluster_num-1]]+[quinn_ttest[cluster_num-1]]+[spritz_ttest[cluster_num-1]]+[probe_ttest[cluster_num-1]]+[flick_ttest[cluster_num-1]]+[valence_ttest[cluster_num-1]]+[juice_thresh]+[quinn_thresh]+[spritz_thresh]+[probe_thresh]+[flick_thresh]+[valence_thresh]+[cluster[1],cluster[2],cluster[3]]+[valencemod_ttest[cluster_num-1]]+[valencemod_thresh]) 
        excel_file.close()  

    reducde_filename_flag = ''
    if subj_list:
        subjs = subj_list
        reduced_filename_flag = 'reducde_subjects' + str(len(subjs))
    if False:   cleanup()
    num_clusters = write_afni_cluster_table()
    if True:   write_means_by_regressor()   
    if True:   make_excel_file()          
    

def write_contrast_table(input_dataset,cluster_mask_name,mask_dataset): # create a mask by FDR mask, but also another mask if provided
    def FDR_mask(): return group_folder + 'FDR_mask2+tlrc '
    def FDR_transformed(): return input_dataset+'_FDR'
    def FDR_masked_dataset(): return input_dataset+'_masked'

    afni('3dFDR','',FDR_transformed(),'-quiet -new -input '+group_folder+input_dataset+'+tlrc -mask ' + FDR_mask()) 
    if mask_dataset == '':
        afni('3dcalc','',FDR_masked_dataset(),'-a '+group_folder+FDR_transformed()+'+tlrc -b '+ FDR_mask()+' -expr a*\'step(b)\'')  # mask dataset      
    else:
        afni('3dcalc','',FDR_masked_dataset(),'-a '+group_folder+FDR_transformed()+'+tlrc -b '+ FDR_mask()+' -c ' +group_folder+mask_dataset+'+tlrc -expr a*\'step(b)*step(c)\'')  # mask dataset        
    print 'rm -f ' + outfilepath_rm(cluster_mask_name) + '+tlrc.*'
    subprocess.call('rm -f ' + outfilepath_rm(cluster_mask_name) + '+tlrc.*', shell=True)   
    afni('3dclust','','','-savemask ' + group_folder+cluster_mask_name + ' -quiet -1thresh 1.96 -dxyz=1 1.01 10 ' + group_folder+FDR_masked_dataset()+'+tlrc\'[0]\'')  
    if os.path.exists(group_folder+cluster_mask_name+'+tlrc.HEAD'):
        write_mean_cluster(cluster_mask_name)
    else:
        print 'Error: will not write clusters table because no clusters were found for ' + cluster_mask_name   




def write_quinine_mask():
    afni('3dbucket','ANOVA_AllStimuli_PLOS+tlrc\'[9]\'','Quinine_statmap','')
    afni('3dFDR','','Quinine_FDR','-quiet -input '+group_folder+'Quinine_statmap+tlrc -mask ' + group_folder + 'FDR_mask2+tlrc ')  
    #afni('3dclust','','Quinine_masked_05','-quiet -1thresh 1.96 -dxyz=1 1.01 10 ' + group_folder+'Quinine_FDR+tlrc\'[0]\'')
    afni('3dclust','','Quinine_masked_05','-quiet -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'Quinine_FDR+tlrc\'[0]\'')
def write_juice_mask():
    afni('3dbucket','ANOVA_AllStimuli_PLOS+tlrc\'[11]\'','Juice_statmap','')    
    afni('3dFDR','','Juice_FDR','-quiet -input '+group_folder+'Juice_statmap+tlrc -mask ' + group_folder + 'FDR_mask2+tlrc ') 
    #afni('3dclust','','Juice_masked_05','-quiet -1thresh 1.96 -dxyz=1 1.01 10 ' + group_folder+'Juice_FDR+tlrc\'[0]\'')   
    afni('3dclust','','Juice_masked_05','-quiet -1thresh 1.96 -dxyz=1 1.01 0 '+ group_folder+'Juice_FDR+tlrc\'[0]\'')   
def write_juice_or_quinn_mask():
    afni('3dcalc','','JuiceOrQuin_masked_05','-a '+group_folder+'Quinine_masked_05+tlrc -b '+group_folder+'Juice_masked_05+tlrc -expr \'or(a,b)\''); return 'JuiceOrQuin_masked_05+tlrc'
def write_juice_and_quinn_mask():
    afni('3dcalc','','JuiceAndQuin_masked_05','-a '+group_folder+'Quinine_masked_05+tlrc -b '+group_folder+'Juice_masked_05+tlrc -expr \'and(a,b)\''); return 'JuiceAndQuin_masked_05+tlrc'



def mask_by_interaction(infile,outfile):
    afni('3dcalc','',outfile,'-a ' + group_folder+infile +'+tlrc -b '+group_folder+'Interaction_masked_05+tlrc -expr \'a*step(b)\'') # mask by juice or quinine  
    return outfile + '+tlrc'    

def mask_by_juice_quinine(infile,outfile):
    afni('3dcalc','',outfile,'-a ' + group_folder+infile +'+tlrc -b '+group_folder+'JuiceOrQuin_masked_05+tlrc -expr \'a*step(b)\'') # mask by juice or quinine     
    return outfile + '+tlrc'

def write_stimulus_mask(): 
    afni('3dFDR','','Stimulus_FDR','-quiet -input '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[1]\' -mask ' + group_folder + 'FDR_mask2+tlrc ') 
    subprocess.call('rm -f ' + outfilepath_rm('Stimulus_masked_05') + '+tlrc.*', shell=True)   
    afni('3dclust','','',' -savemask '+ group_folder+'Stimulus_masked_05 -quiet -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'Stimulus_FDR+tlrc')   
    write_mean_cluster( 'Stimulus_masked_05')     

def make_stimulus_figure():
    subprocess.call('rm -f ' + outfilepath_rm('Stimulus_masked_05_figure') + '+tlrc.*', shell=True)
    afni('3dclust','','Stimulus_masked_01','-1thresh 2.573 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.01
    afni('3dclust','','Stimulus_masked_005','-1thresh 2.807 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.005
    afni('3dclust','','Stimulus_masked_001','-1thresh 3.227 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.001
    afni('3dclust','','Stimulus_masked_0005','-1thresh 3.482 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.0005
    afni('3dcalc','','Stimulus_masked_05_figure','-a '+group_folder+'Stimulus_masked_05+tlrc -b '+group_folder+'Stimulus_masked_01+tlrc -c '+group_folder+'Stimulus_masked_005+tlrc -expr \'step(a)+step(b)+step(c)\'')       


def write_valence_regions(): 
    afni('3dttest','','Valence', '-paired -set1 ' + dataset_string('','PLOS_juice_sum+tlrc','0') + ' -set2 ' + dataset_string('','PLOS_quinn_sum+tlrc','0'))
    write_contrast_table('Valence','Valence_stimulus-masked','Stimulus_masked_05')    

def write_quinine_regions(): write_contrast_table('Quinine_statmap','Quinine_stimulus-masked','Stimulus_masked_05')    
def write_juice_regions(): write_contrast_table('Juice_statmap','Juice_stimulus-masked','Stimulus_masked_05')    
def write_hedonic_versus_visual():   
    afni('3dFDR','','StimVsHedonic_FDR','-quiet -input '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[29]\' -mask ' + group_folder + 'FDR_mask2+tlrc ') 
    afni('3dcalc','','StimVHedonic_masked','-a '+group_folder+'StimVsHedonic_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Stimulus_masked_05+tlrc -expr \'a*step(b)*step(c)\'') 
    subprocess.call('rm -f ' + outfilepath_rm('StimVsHedonic_masked_05') + '+tlrc.*', shell=True)    
    afni('3dclust','','',' -savemask ' + group_folder+'StimVsHedonic_masked_05 -quiet -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'StimVHedonic_masked+tlrc')   
    write_mean_cluster( 'StimVsHedonic_masked_05')         

def write_quinine_negative_bold():
    def make_mask(cluster_mask_name):
        #write_quinine_mask()
        afni('3dcalc','','QuinineNegBold_masked_05','-a '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[9]\' -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Quinine_masked_05+tlrc -expr \'isnegative(a)*step(b)*c\'')  # mask dataset    
        make_cluster_mask('QuinineNegBold_masked_05+tlrc',cluster_mask_name)

    make_mask('QuinineNegBold_clusters')
    write_mean_cluster('QuinineNegBold_clusters') 

def write_juice_negative_bold():
    def make_mask(cluster_mask_name):
        write_juice_mask()
        afni('3dcalc','','JuiceNegBold_masked_05','-a '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[11]\' -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Juice_masked_05+tlrc -expr \'isnegative(a)*step(b)*c\'')  # mask dataset    
        make_cluster_mask('JuiceNegBold_masked_05+tlrc',cluster_mask_name)

    make_mask('JuiceNegBold_clusters')
    write_mean_cluster('JuiceNegBold_clusters')     

def write_liquid_vs_baseline():
    afni('3dttest','','Liquid_vs_baseline', '-base1 0 -set2 ' +dataset_string('','PLOS_spritz_all_mean+tlrc','0'))  
    afni('3dFDR','','Liquid_vs_baseline_FDR','-quiet -input '+group_folder+'Liquid_vs_baseline+tlrc -mask ' + group_folder + 'FDR_mask2+tlrc ') 
    afni('3dcalc','','Liquid_vs_baseline_masked','-a '+group_folder+'Liquid_vs_baseline_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Stimulus_masked_05+tlrc -expr \'a*step(b)*step(c)\'') 
    subprocess.call('rm -f ' + outfilepath_rm('Liquid_vs_baseline_masked_05') + '+tlrc.*', shell=True) 
    afni('3dclust','','',' -savemask ' + group_folder+'Liquid_vs_baseline_masked_05 -quiet -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'Liquid_vs_baseline_masked+tlrc\'[1]\'')      
    write_mean_cluster( 'Liquid_vs_baseline_masked_05')  


# def write_spritz_negative_bold():
#     def make_mask(cluster_mask_name):
#         #write_juice_mask()
#         afni('3dcalc','','SpritzNegBold_masked_05','-a '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[11]\' -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Spritz_masked_05+tlrc -expr \'isnegative(a)*step(b)*c\'')  # mask dataset    
#         make_cluster_mask('SpritzNegBold_masked_05+tlrc',cluster_mask_name)

#     make_mask('SpritzNegBold_clusters')
#     write_mean_cluster('SpritzNegBold_clusters')  

def write_load_mask():
    afni('3dFDR','','Load_FDR','-quiet -input '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[3]\' -mask ' + group_folder + 'FDR_mask2+tlrc ') 
    subprocess.call('rm -f ' + outfilepath_rm('Load_masked_05') + '+tlrc.*', shell=True)
    afni('3dclust','','',' -savemask '+group_folder+'Load_masked_05 -quiet -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'Load_FDR+tlrc\'[0]\'')
    write_mean_cluster( 'Load_masked_05')   

def write_salience_as_load_contrast(): 
    def make_mask(cluster_mask_name):
        afni('3dttest','','SalienceModulation_load', '-paired -set1 ' + dataset_string('','liquid_3+tlrc','0') + ' -set2 ' + dataset_string('','liquid_0+tlrc','0'))
            #afni('3dclust','','Load_mask','-1thresh 7.19999 -dxyz=1 1.01 0 ' + group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[3]\'')                    # make map of q=0.03

        afni('3dcalc','','SalienceMod_load_masked','-a '+group_folder+'SalienceModulation_load+tlrc -b '+group_folder+'FDR_mask2+tlrc -expr a*\'step(b)\'')  # mask dataset
            #subprocess.call('rm -f SalienceMod_FDR_temp+tlrc.*', shell=True)   
        afni('3dFDR','','SalienceMod_FDR_temp_load','-quiet -input '+group_folder+'SalienceMod_load_masked+tlrc -mask ' + group_folder + 'FDR_mask2+tlrc ')                #salience mask           
        afni('3dcalc','','SalienceMod_FDR_load','-a '+group_folder+'SalienceMod_FDR_temp_load+tlrc -b ' +group_folder+'Load_masked_05+tlrc -c ' +group_folder+'FDR_mask2+tlrc -expr a*\'step(b)*step(c)\'')  
                #afni('3dcalc','','SalienceMod_FDR','-a '+group_folder+'SalienceMod_FDR_temp+tlrc -b '+group_folder+'Load_mask+tlrc -c ' +group_folder+'Load_masked_05+tlrc -expr a*\'step(b)\'')   ######### c not included!!!!!!     

        subprocess.call('rm -f ' + outfilepath_rm('SalienceMod_FDR_masked_load_025') + '+tlrc.*', shell=True)
        afni('3dclust','','','-savemask '+group_folder+'SalienceMod_FDR_masked_load_025 -1thresh 2.242 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR_load+tlrc\'[1]\'')                    # make map of q=0.025
        subprocess.call('rm -f ' + outfilepath_rm('SalienceMod_FDR_masked_load_03') + '+tlrc.*', shell=True)
        afni('3dclust','','','-savemask '+group_folder+'SalienceMod_FDR_masked_load_03 -1thresh 2.170 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR_load+tlrc\'[1]\'')                    # make map of q=0.025

        subprocess.call('rm -f ' + outfilepath_rm('SalienceMod_FDR_masked_load_05') + '+tlrc.*', shell=True)               
        afni('3dclust','','','-savemask '+group_folder+'SalienceMod_FDR_masked_load_05 -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR_load+tlrc\'[1]\'')                    # make map of q=0.05

        # afni('3dclust','','SalienceMod_FDR_masked_03_load','-1thresh 2.170 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.03
        # afni('3dclust','','SalienceMod_FDR_masked_01_load','-1thresh 2.574 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.01
        # afni('3dclust','','SalienceMod_FDR_masked_005_load','-1thresh 2.807 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.005
        # afni('3dclust','','SalienceMod_FDR_masked_001_load','-1thresh 3.229 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.001
        # afni('3dclust','','SalienceMod_FDR_masked_0005_load','-1thresh 3.481 -dxyz=1 1.01 0 ' + group_folder+'SalienceMod_FDR+tlrc\'[1]\'')                    # make map of q=0.0005
        # afni('3dcalc','','SalienceMod_FDR_mask_combination_025_load','-a '+group_folder+'SalienceMod_FDR_masked_025_load+tlrc -b '+group_folder+'SalienceMod_FDR_masked_01_load+tlrc -c '+group_folder+'SalienceMod_FDR_masked_005_load+tlrc -expr \'step(a)+step(b)+step(c)\'')       
        
        #make_cluster_mask('SalienceMod_FDR_masked_load_025+tlrc',cluster_mask_name + '_025')
        make_cluster_mask('SalienceMod_FDR_masked_load_03+tlrc',cluster_mask_name + '_03')
        make_cluster_mask('SalienceMod_FDR_masked_load_05+tlrc',cluster_mask_name + '_05')
 
    #make_mask('SalienceMod_05_clusters') 
    #write_mean_cluster('SalienceMod_05_clusters')
    make_mask('SalienceMod_clusters') 
    #write_mean_cluster('SalienceMod_clusters_025')
    write_mean_cluster('SalienceMod_clusters_03')
    write_mean_cluster('SalienceMod_clusters_05')

    return  'Calculated salience modulation'

def write_probe_modulation():
	afni('3dttest','','ProbeMod', '-paired -set1 ' + dataset_string('','AttentionEnjoy2_deconv+tlrc','2') + ' -set2 ' + dataset_string('','AttentionEnjoy2_deconv+tlrc','6'))
	afni('3dbucket','ProbeMod+tlrc\'[1]\'','ProbeMod_statmap','') 
	write_contrast_table('ProbeMod_statmap','ProbeMod_masked_05','Load_masked_05')

def write_flickering_checkerboard_modulation():  # not writing
	afni('3dttest','','FlickMod', '-paired -set1 ' + dataset_string('','AttentionEnjoy2_deconv+tlrc','26') + ' -set2 ' + dataset_string('','AttentionEnjoy2_deconv+tlrc','30'))
	afni('3dbucket','FlickMod+tlrc\'[1]\'','FlickMod_statmap','') 
	write_contrast_table('FlickMod_statmap','FlickMod_masked_05','Load_masked_05')	


def write_interaction_mask():
    afni('3dFDR','','Interaction_FDR','-quiet -input '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[5]\' -mask ' + group_folder + 'FDR_mask2+tlrc ') 
    subprocess.call('rm -f ' + outfilepath_rm('Interaction_masked_05') + '+tlrc.*', shell=True)
    afni('3dclust','','', ' -savemask '+group_folder+'Interaction_masked_05 -quiet -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'Interaction_FDR+tlrc\'[0]\'') 
    write_mean_cluster( 'Interaction_masked_05') 

def write_valence_modulation():  # this is the interaction  
    def make_mask(cluster_mask_name):
#        afni('3dttest','','ValenceMod', '-paired -set1 ' + dataset_string('','juice_diff+tlrc','0') + ' -set2 ' + dataset_string('','inv_quin_diff+tlrc','0'))    
        afni('3dttest','','ValenceMod', '-base1 0 -set2 ' + dataset_string('','PLOS_liquid_load_interaction_beta+tlrc','0'))    

        afni('3dFDR','','ValenceMod_FDR','-quiet -new -input '+group_folder+'ValenceMod+tlrc\'[1]\' -mask ' + group_folder + 'FDR_mask2+tlrc ')               
        afni('3dcalc','','ValenceMod_masked','-a '+group_folder+'ValenceMod_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Interaction_masked_05+tlrc  -expr \'a*step(b)*step(c)\'')  # mask dataset            

        subprocess.call('rm -f ' + outfilepath_rm('ValenceMod_FDR_masked_05') + '+tlrc.*', shell=True) 
        afni('3dclust','','','-savemask '+group_folder+'ValenceMod_FDR_masked_05 -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'ValenceMod_masked+tlrc')                    #
        
        mask_by_juice_quinine('ValenceMod_FDR_masked_05','ValenceMod_MaskByJuiceQuinineFDR')
        mask_by_interaction('ValenceMod_MaskByJuiceQuinineFDR','ValenceMod_MaskByJuiceQuinineInteractionFDR')
        make_cluster_mask('ValenceMod_MaskByJuiceQuinineInteractionFDR+tlrc',cluster_mask_name)

    make_mask('ValenceMod_clusters')
    write_mean_cluster('ValenceMod_clusters')
    return  'Calculated valence modulation'    

def write_salience_as_interaction_contrast():  # probably shoudl depricate
    def make_mask(cluster_mask_name):
        afni('3dttest','','Salience', '-paired -set1 ' + dataset_string('','liquid_3+tlrc','0') + ' -set2 ' + dataset_string('','liquid_0+tlrc','0'))
        #afni('3dttest','','Salience', '-paired -set1 ' + dataset_string('','PLOS_juice_sum_beta+tlrc','0') + ' -set2 ' + dataset_string('','PLOS_neg_quinn_sum_beta+tlrc','0'))
        #afni('3dttest','','Salience', '-base1 0 -set2 ' + dataset_string('','PLOS_liquid_load_interaction+tlrc','0'))  # one sample ttest against zero

        afni('3dFDR','','Salience_FDR','-quiet -new -input '+group_folder+'Salience+tlrc\'[1]\' -mask ' + group_folder + 'FDR_mask2+tlrc ')   
        afni('3dcalc','','Salience_masked','-a '+group_folder+'Salience_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Interaction_masked_05+tlrc -expr \'a*step(b)*step(c)\'')  # mask dataset       
  
        # subprocess.call('rm -f ' + outfilepath_rm('Salience_FDR_masked_025') + '+tlrc.*', shell=True)
        # afni('3dclust','','','-savemask '+group_folder+'Salience_FDR_masked_025 -1thresh 2.242 -dxyz=1 1.01 0 ' + group_folder+'Salience_FDR+tlrc')                    # make map of q=0.03        
        # #afni('3dcalc','','Salience_FDR_mask_combination_025','-a '+group_folder+'Salience_FDR_masked_025+tlrc -b '+group_folder+'Salience_FDR_masked_01+tlrc -c '+group_folder+'Salience_FDR_masked_005+tlrc -expr \'step(a)+step(b)+step(c)\'')       
        # make_cluster_mask('Salience_FDR_masked_025+tlrc',cluster_mask_name)

        subprocess.call('rm -f ' + outfilepath_rm('Salience_FDR_masked_05') + '+tlrc.*', shell=True)
        afni('3dclust','','','-savemask '+group_folder+'Salience_FDR_masked_05 -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'Salience_FDR+tlrc')                    # make map of q=0.03        
        #afni('3dcalc','','Salience_FDR_mask_combination_025','-a '+group_folder+'Salience_FDR_masked_025+tlrc -b '+group_folder+'Salience_FDR_masked_01+tlrc -c '+group_folder+'Salience_FDR_masked_005+tlrc -expr \'step(a)+step(b)+step(c)\'')       
        make_cluster_mask('Salience_FDR_masked_05+tlrc',cluster_mask_name)

    make_mask('Salience_clusters') 
    write_mean_cluster('Salience_clusters') 
    return  'Calculated salience modulation'       

def write_interact_probe_modulation():
    afni('3dttest','','ProbeMod', '-paired -set1 ' + dataset_string('','AttentionEnjoy2_deconv+tlrc','2') + ' -set2 ' + dataset_string('','AttentionEnjoy2_deconv+tlrc','6'))
    afni('3dFDR','','InteractProbe_FDR','-quiet -new -input '+group_folder+'ProbeMod+tlrc\'[1]\' -mask ' + group_folder + 'FDR_mask2+tlrc ')   
    afni('3dcalc','','InteractProbe_masked','-a '+group_folder+'InteractProbe_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'Interaction_masked_05+tlrc -expr \'a*step(b)*step(c)\'')  # mask dataset 

    subprocess.call('rm -f ' + outfilepath_rm('InteractProbe_masked_05') + '+tlrc.*', shell=True)
    afni('3dclust','','','-savemask '+group_folder+'InteractProbe_masked_05 -1thresh 1.96 -dxyz=1 1.01 0 ' + group_folder+'InteractProbe_masked+tlrc')   
    write_mean_cluster('InteractProbe_masked_05')
    #write_contrast_table('ProbeMod','InteractProbeMod_masked_05','Interaction_masked_05')    


def write_visual_region():
    subprocess.call('rm -f ' + outfilepath_rm('Flicker_FDR_masked_0001') + '+tlrc.*', shell=True)
    subprocess.call('rm -f ' + outfilepath_rm('Flicker_FDR_masked_005') + '+tlrc.*', shell=True)    

    afni('3dFDR','','Flicker_FDR','-quiet -new -input '+group_folder+'ANOVA_AllStimuli_PLOS+tlrc\'[13]\' -mask ' +group_folder+'FDR_mask2+tlrc ')
    afni('3dclust','','','-savemask '+group_folder+'Flicker_FDR_masked_005 -1thresh 2.808 -dxyz=1 1.01 0 ' + group_folder+'Flicker_FDR+tlrc') 
    write_mean_cluster('Flicker_FDR_masked_005')  
    afni('3dclust','','','-savemask '+group_folder+'Flicker_FDR_masked_001 -1thresh 3.282 -dxyz=1 1.01 0 ' + group_folder+'Flicker_FDR+tlrc') 
    write_mean_cluster('Flicker_FDR_masked_001')      


def write_3BackBlock_regions():
    afni('3dttest','','3BackBlock_ttest', '-base1 0 -set2 ' + dataset_string('','AttentionEnjoy2_3BackBlock_deconv+tlrc','2') )
    afni('3dFDR','','3BackBlock_FDR','-quiet -new -input '+group_folder+'3BackBlock_ttest+tlrc\'[1]\' -mask ' + group_folder + 'FDR_mask2+tlrc')   
    afni('3dcalc','','3BackBlock_masked','-a '+group_folder+'3BackBlock_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -expr \'a*step(b)\'')  # mask dataset      

    # subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_masked_05') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'3BackBlock_masked_05 -1thresh 1.96 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_masked+tlrc')     
    # write_mean_cluster('3BackBlock_masked_05')

    # subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_masked_005') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'3BackBlock_masked_005 -1thresh 2.808 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_masked+tlrc')     
    # write_mean_cluster('3BackBlock_masked_005')

    subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_masked_001') + '+tlrc.*', shell=True)
    afni('3dclust','','','-savemask '+group_folder+'3BackBlock_masked_001 -1thresh 3.282 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_masked+tlrc')     
    write_mean_cluster('3BackBlock_masked_001')

def write_3BackBlock_NoBaseline_regions():
    afni('3dttest','','3BackBlock_NoBaseline_ttest', '-base1 0 -set2 ' + dataset_string('','AttentionEnjoy2_3BackBlock_NoBaseline_deconv+tlrc','2') )
    afni('3dFDR','','3BackBlock_NoBaseline_FDR','-quiet -new -input '+group_folder+'3BackBlock_NoBaseline_ttest+tlrc\'[1]\' -mask ' + group_folder + 'FDR_mask2+tlrc')   
    afni('3dcalc','','3BackBlock_NoBaseline_masked','-a '+group_folder+'3BackBlock_NoBaseline_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -expr \'a*step(b)\'')  # mask dataset    
    afni('3dcalc','','3BackBlock_NoBaseline_masked_positive','-a '+group_folder+'3BackBlock_NoBaseline_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'3BackBlock_NoBaseline_ttest+tlrc\'[0]\' -expr \'a*step(b)*ispositive(c)\'')  # mask dataset         
    afni('3dcalc','','3BackBlock_NoBaseline_masked_negative','-a '+group_folder+'3BackBlock_NoBaseline_FDR+tlrc -b '+group_folder+'FDR_mask2+tlrc -c '+group_folder+'3BackBlock_NoBaseline_ttest+tlrc\'[0]\' -expr \'a*step(b)*isnegative(c)\'')  # mask dataset    

    # subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_NoBaseline_masked_05') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'3BackBlock_NoBaseline_masked_05 -1thresh 1.96 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_NoBaseline_masked+tlrc')     
    # write_mean_cluster('3BackBlock_NoBaseline_masked_05')

    # subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_NoBaseline_masked_005') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'3BackBlock_NoBaseline_masked_005 -1thresh 2.808 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_NoBaseline_masked+tlrc')     
    # write_mean_cluster('3BackBlock_NoBaseline_masked_005')

    # subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_NoBaseline_masked_01') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'3BackBlock_NoBaseline_masked_01 -1thresh 1.645 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_NoBaseline_masked+tlrc')     
    # write_mean_cluster('3BackBlock_NoBaseline_masked_01')    

    subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_NoBaseline_masked_01_positive') + '+tlrc.*', shell=True)
    afni('3dclust','','','-savemask '+group_folder+'3BackBlock_NoBaseline_masked_01_positive -1thresh 1.645 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_NoBaseline_masked_positive+tlrc')     
    write_mean_cluster('3BackBlock_NoBaseline_masked_01_positive')    

    subprocess.call('rm -f ' + outfilepath_rm('3BackBlock_NoBaseline_masked_01_negative') + '+tlrc.*', shell=True)
    afni('3dclust','','','-savemask '+group_folder+'3BackBlock_NoBaseline_masked_01_negative -1thresh 1.645 -dxyz=1 1.01 10 ' + group_folder+'3BackBlock_NoBaseline_masked_negative+tlrc')     
    write_mean_cluster('3BackBlock_NoBaseline_masked_01_negative')        

def make_and_write_cluster(input,output,input_brik = 0, thresh = 1.96, voxels_in_thresh = 10):
    subprocess.call('rm -f ' + outfilepath_rm(output) + '+tlrc.*', shell=True)
    afni('3dclust','','','-savemask '+group_folder+output+' -1thresh ' + str(thresh) + ' -dxyz=1 1.01 ' + str(voxels_in_thresh) + ' ' + group_folder+input+'+tlrc\'['+str(input_brik)+']\'')     
    write_mean_cluster(output)    

def write_clusters_from_bartra_paper():

    # subprocess.call('rm -f ' + outfilepath_rm('bartra_rotated_Pos_not_greater_masked') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'bartra_rotated_Pos_not_greater_masked -1thresh 0.01 -dxyz=1 1.01 10 ' + group_folder+'bartra_rotation_Pos_not_greater+tlrc')     
    # write_mean_cluster('bartra_rotated_Pos_not_greater_masked')   

    # subprocess.call('rm -f ' + outfilepath_rm('bartra_rotated_POS_greater_only_masked') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'bartra_rotated_POS_greater_only_masked -1thresh 0.01 -dxyz=1 1.01 10 ' + group_folder+'bartra_rotation_Pos_greater_only+tlrc')     
    # write_mean_cluster('bartra_rotated_POS_greater_only_masked')      

    # subprocess.call('rm -f ' + outfilepath_rm('bartra_rotated_fig3C_conjunc_masked') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'bartra_rotated_fig3C_conjunc_masked -1thresh 0.01 -dxyz=1 1.01 10 ' + group_folder+'bartra_rotated_fig3C_conjunc+tlrc')     
    # write_mean_cluster('bartra_rotated_fig3C_conjunc_masked')      
   
    # subprocess.call('rm -f ' + outfilepath_rm('nucleus_accumbens_mask_EPIres_masked') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'nucleus_accumbens_mask_EPIres_masked -1thresh 0.01 -dxyz=1 1.01 10 ' + group_folder+'nucleus_accumbens_mask_EPIres+tlrc')     
    # write_mean_cluster('nucleus_accumbens_mask_EPIres_masked')     

    # subprocess.call('rm -f ' + outfilepath_rm('bartra_nucleus_accumbens_conj_pos_masked') + '+tlrc.*', shell=True)
    # afni('3dclust','','','-savemask '+group_folder+'bartra_nucleus_accumbens_conj_pos_masked -1thresh 0.01 -dxyz=1 1.01 10 ' + group_folder+'bartra_nucleus_accumbens_conj_pos+tlrc')     
    # write_mean_cluster('bartra_nucleus_accumbens_conj_pos_masked') 

   
    subprocess.call('rm -f ' + outfilepath_rm('bartra_rotated_POS_greater_only_debridged_masked') + '+tlrc.*', shell=True)
    afni('3dclust','','','-savemask '+group_folder+'bartra_rotated_POS_greater_only_debridged_masked -1thresh 0.01 -dxyz=1 1.01 10 ' + group_folder+'bartra_rotated_POS_greater_only_debridged+tlrc')     
    write_mean_cluster('bartra_rotated_POS_greater_only_debridged_masked') 

def write_weak_valence_clusters_for_trendlevel_analysis():
    make_and_write_cluster('ANOVA_AllStimuli_PLOS','Valence_veryweakthresh',27,2.130,15)   # 0.05

def write_ultraweak_valence_clusters_for_trendlevel_analysis():
    make_and_write_cluster('ANOVA_AllStimuli_PLOS','Valence_ultraweakthresh',27,1.753,20)   # 0.10

def write_results_for_atlas_CA_27_ML():
    #write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres') 
    #write_mean_cluster('CA_N27_ML_region27_LRectalGyrus_lowres')
    #write_mean_cluster('CA_N27_ML_region9_LMiddleOrbGyrus_lowres')
    #write_mean_cluster('CA_N27_ML_region25_LMidOrbitalGyrus_lowres')
    #write_mean_cluster('CA_N27_ML_region31_LAnteriorCingulate_lowres')
    #write_mean_cluster('CA_N27_ML_region21_LOlfactoryCortex_lowres')
    #write_mean_cluster('CA_N27_ML_region22_ROlfactoryCortex_lowres')
    #write_mean_cluster('CA_N27_ML_region26_RMidOrbitalGyrus_lowres')
    #write_mean_cluster('CA_N27_ML_region6_RSuperOrbitalGyrus_lowres') 
    #write_mean_cluster('CA_N27_ML_region28_RRectalGyrus_lowres')
    #write_mean_cluster('CA_N27_ML_region10_RMiddleOrbGyrus_lowres')
    #write_mean_cluster('CA_N27_ML_region32_RAnteriorCingulate_lowres')


    # # Drop subjects in order based on how litte the L SuperOrbital Gyrus survives the 3DAutomask cluster
    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae16","ae17","ae19","ae20","ae21","ae22","ae23","ae24","ae25","ae27","ae28","ae29","ae30","ae31","ae32"]) # drop ae18 
    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae16","ae17","ae19","ae20","ae21","ae22","ae23","ae25","ae27","ae28","ae29","ae30","ae31","ae32"]) # drop ae24
    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae16","ae17","ae19","ae20","ae21","ae22","ae23","ae25","ae27","ae28","ae29","ae31","ae32"]) # drop ae28 ae20
    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae16","ae17","ae19","ae20","ae21","ae22","ae23","ae25","ae27","ae29","ae31","ae32"]) # drop ae20
    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae16","ae17","ae19","ae21","ae22","ae23","ae25","ae27","ae29","ae31","ae32"]) # drop ae20

    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae16","ae17","ae19","ae22","ae23","ae27","ae29","ae31","ae32"]) # drop 21 & 25
    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae16","ae17","ae19","ae22","ae27","ae29","ae31","ae32"]) # drop 23
    # write_mean_cluster('CA_N27_ML_region5_LSuperOrbital_lowres',["ae17","ae19","ae22","ae27","ae29","ae31","ae32"]) # drop 16

    write_mean_cluster('CA_N27_ML_region6_RSuperOrbitalGyrus_lowres',["ae16","ae17","ae19","ae20","ae21","ae22","ae23","ae25","ae27","ae29","ae31","ae32"]) # drop ae20
    write_mean_cluster('CA_N27_ML_region6_RSuperOrbitalGyrus_lowres',["ae16","ae17","ae19","ae21","ae22","ae23","ae25","ae27","ae29","ae31","ae32"]) # drop ae20



def write_overlap_between_clusters_and_edge_of_brainsignal():

    clusters = [ \
                'CA_N27_ML_region5_LSuperOrbital',
                'CA_N27_ML_region6_RSuperOrbitalGyrus',
                'CA_N27_ML_region9_LMiddleOrbGyrus',
                'CA_N27_ML_region10_RMiddleOrbGyrus',
                'CA_N27_ML_region21_LOlfactoryCortex',
                'CA_N27_ML_region22_ROlfactoryCortex',
                'CA_N27_ML_region25_LMidOrbitalGyrus',
                'CA_N27_ML_region26_RMidOrbitalGyrus',
                'CA_N27_ML_region27_LRectalGyrus',
                'CA_N27_ML_region28_RRectalGyrus',
                'CA_N27_ML_region31_LAnteriorCingulate',
                'CA_N27_ML_region32_RAnteriorCingulate',
                ]

    def subj_mask_file(subj): return root_folder + subj + '/mask_at+tlrc '
    def cluster_file(cluster): return group_folder + cluster + '_lowres+tlrc.'
    def subject_overlap_file(subj,cluster): return root_folder + subj + '/overlap_brain_' + cluster + '.txt'
    def group_overlap_file(): return group_folder + 'overlap_summary.txt'

    subprocess.call('rm -f ' + group_overlap_file(), shell=True)   

    for subj in subjs:                                         # write each overlap file to one common file for later reorganization
        for cluster in clusters:
            afni('3dABoverlap','','', subj_mask_file(subj) + cluster_file(cluster) + ' > ' + subject_overlap_file(subj,cluster))
            subprocess.call('cat ' + subject_overlap_file(subj,cluster) + '>> ' + group_overlap_file(), shell=True)


def parse_overlap_file():
    # defines inherited from above, which can be removed when this code gets combined with the write_overlap_between_clusters_and_edge_of_brainsignal function
    def group_overlap_file(): return group_folder + 'overlap_summary.txt'


    # actual code
    def group_overlap_file_csv(): return group_folder + 'overlap_summary.csv'

    with open(group_overlap_file()) as raw_overlap_data:
        header_row = ['subj_num']+['ROI_num']+['ROI_name']+['Voxels_in_mask']+['Voxels_in_ROI']+['Voxles_ROI_or_mask']+['Voxels_ROI_&_mask']+['Mask_voxels_not_in_ROI']+['ROI_voxels_not_in_mask']+\
                        ['Mask_percent_not_in_ROI']+['ROI_percent_not_in_mask']+['Radius_X']+['Radius_Y']+['Radius_Z']
        excel_file = open(group_overlap_file_csv(),'w')
        excel_csv = csv.writer(excel_file, delimiter=',') 
        excel_csv.writerow(header_row)   

        while True:
            line1 = raw_overlap_data.readline()
            line2 = raw_overlap_data.readline()
            line3 = raw_overlap_data.readline()
            if not line1: break

            subj_num = re.search('/ae(\d+)/',line1).group(1)
            region_num = re.search('_region(\d+)_',line1).group(1)
            cluster_name = re.search('_region\d+_([a-zA-Z]+)_lowres',line1).group(1)
            clust_stats = re.search('^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+).(\d+)\s+(\d+).(\d+)\s+(\d+).(\d+)\s+(\d+).(\d+)\s+(\d+).(\d+)\s+',line3)

            excel_csv.writerow([subj_num]+[region_num]+[cluster_name]+[clust_stats.group(1)]+[clust_stats.group(2)]+[clust_stats.group(3)]+[clust_stats.group(4)]+[clust_stats.group(5)]+[clust_stats.group(6)]+\
                                [clust_stats.group(7)+'.'+clust_stats.group(8)]+[clust_stats.group(9)+'.'+clust_stats.group(10)]+[clust_stats.group(11)+'.'+clust_stats.group(12)]+[clust_stats.group(13)+'.'+clust_stats.group(14)]+\
                                [clust_stats.group(15)+'.'+clust_stats.group(16)])
                            





def make_brainpics_for_figures():
    def make_figure(FDR_dataset,mask_dataset,output_dataset):  # for two-tailed
        subprocess.call('rm -f ' + outfilepath_rm(output_dataset) + '+tlrc.*', shell=True)     
        afni('3dcalc','',output_dataset,'-a '+group_folder+FDR_dataset+'+tlrc -b '+ group_folder+mask_dataset+'+tlrc -expr \'-log10(fizt_t2p(a))*step(b)-1.30*step(b)\'') # 1.30 is just under prob of 0.05 - helps with intensity bar

    def make_figure_one_tailed(FDR_dataset,mask_dataset,output_dataset):
        subprocess.call('rm -f ' + outfilepath_rm(output_dataset) + '+tlrc.*', shell=True)     
        afni('3dcalc','',output_dataset,'-a '+group_folder+FDR_dataset+'+tlrc -b '+ group_folder+mask_dataset+'+tlrc -expr \'-log10(fizt_t2p(a)/2)*step(b)-1.30*step(b)\'') # 1.30 is just under prob of 0.05 - helps with intensity bar

    #make_figure('Stimulus_FDR','Stimulus_masked_05','Stimulus_masked_05_figure')  
    make_figure('Load_FDR','Load_masked_05','Load_masked_05_figure')  
    #make_figure('Interaction_FDR','Interaction_masked_05','Interaction_masked_05_figure')   
    #make_figure('StimVsHedonic_FDR','StimVsHedonic_masked_05','StimVsHedonic_masked_05_figure') 
    #make_figure('Stimulus_FDR','Stimulus_masked_05','Stimulus_masked_05_figure')
    ##make_figure('SalienceMod_FDR_single','SalienceMod_FDR_masked_025','SalienceMod_FDR_masked_025_figure')
    #make_figure('Flicker_FDR','Flicker_FDR_masked_005','Flicker_FDR_masked_005_fig')
    #make_figure('Flicker_FDR','Flicker_FDR_masked_001','Flicker_FDR_masked_001_fig')
    #make_figure('Liquid_vs_baseline_FDR','Liquid_vs_baseline_masked_05','Liquid_vs_baseline_masked_05_fig')

    #make_figure('ValenceMod_FDR','ValenceMod_FDR_masked_05','Valence_FDR_masked_05_fig')
    ##make_figure('Salience_FDR','Salience_FDR_masked_025','Salience_FDR_masked_025_fig')
    #make_figure('Salience_FDR','Salience_FDR_masked_05','Salience_FDR_masked_05_fig')
    #make_figure('ProbeMod_FDR','InteractProbe_masked_05','InteractProbe_masked_05_fig')

    # make_figure('3BackBlock_FDR','3BackBlock_masked_05','3BackBlock_masked_05_fig')
    # make_figure('3BackBlock_FDR','3BackBlock_masked_005','3BackBlock_masked_005_fig')
    # make_figure('3BackBlock_FDR','3BackBlock_masked_001','3BackBlock_masked_001_fig')

    #make_figure('3BackBlock_NoBaseline_FDR','3BackBlock_NoBaseline_masked_05','3BackBlock_NoBaseline_masked_05_fig')
    #make_figure('3BackBlock_NoBaseline_FDR','3BackBlock_NoBaseline_masked_005','3BackBlock_NoBaseline_masked_005_fig')
    #make_figure_one_tailed('3BackBlock_NoBaseline_FDR','3BackBlock_NoBaseline_masked_01','3BackBlock_NoBaseline_masked_01_fig')
    #make_figure_one_tailed('3BackBlock_NoBaseline_FDR','3BackBlock_NoBaseline_masked_01_positive','3BackBlock_NoBaseline_masked_01_fig_positive')
    #make_figure_one_tailed('3BackBlock_NoBaseline_FDR','3BackBlock_NoBaseline_masked_01_negative','3BackBlock_NoBaseline_masked_01_fig_negative')        




def main():
        print '\n<><><><><><><><>\n\n\n\n\n\n\n'

        # write all the support files: masks for regressors, etc
        if False: subject_preprocessing_for_second_level_stats()        
        if False: write_quinine_mask()
        if False: write_juice_mask()
        if False: write_juice_or_quinn_mask() 
        if False: write_juice_and_quinn_mask() 

        # GO THROUGH all= the contrasts
    # Stimulus ANOVA term
    	if False: write_stimulus_mask()  # stimulus-responsive regions
        if False: write_valence_regions()          # valence
        if False: write_juice_regions()          # juice        
        if False: write_quinine_regions()       # quinine 
        if False: write_hedonic_versus_visual()
        if False: write_quinine_negative_bold()        # quinine negative bold         # where did we get negative BOLD for quinine? hypothalamus -lateral part of ventromedial nuclues of hypothalams - food intake - we are real close to this
        if False: write_juice_negative_bold()        # juice negative bold   
        if False: write_liquid_vs_baseline()
 
    # Load ANOVA term
        if False: write_load_mask() # load-sensitive regions
        if False: write_salience_as_load_contrast()             # Salience modulation masked by juice or quinine   
        if False: write_probe_modulation()         # probe (n-back)
        if False: write_flickering_checkerboard_modulation()          # checkerboard

    # interaction
        if False: write_interaction_mask()  
        if False: write_valence_modulation()         # valence modulation
        #if False: write_salience_as_interaction_contrast()
        if False: write_interact_probe_modulation() 

    # Vision area 
        if False: write_visual_region()   

    # Chacing a valence effect
        if False: write_clusters_from_bartra_paper()
        if False: write_weak_valence_clusters_for_trendlevel_analysis()
        if False: write_ultraweak_valence_clusters_for_trendlevel_analysis()
        if True: write_results_for_atlas_CA_27_ML()
        if False: write_overlap_between_clusters_and_edge_of_brainsignal()
        if False: parse_overlap_file()        

    # Block effect
        if False: write_3BackBlock_regions()
        if False: write_3BackBlock_NoBaseline_regions()    

    # make afni files for figures 
        if False: make_brainpics_for_figures()


if __name__ == "__main__": main()   
   