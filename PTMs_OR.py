#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables

__author__ = ["Diego Mena Santos"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.3.0"
__maintainer__ = "Diego Mena Santos"
__email__ = "diego.mena@cnic.es"
__status__ = "Development"




import pandas as pd
import numpy as np
import statsmodels.api as sm
import configparser
import warnings
from operator import itemgetter
import argparse
import os
import sys
import logging
import openpyxl

warnings.filterwarnings("ignore")



""" Functions """


def check_folder(outfolder):
    if not os.path.exists(outfolder):
        os.makedirs(outfolder, exist_ok=True)


def read_file(path, raw_table):

    if raw_table==True:
        limma_table=pd.read_csv(path,sep='\t', header=[0,1])
        return limma_table
    else:
        exp_table=pd.read_csv(path,sep='\t', header=[0,1])
        return exp_table
    


def filter_table(raw_table, group_column_header, nm_label,nm_stat_header,
                  mod_stat_header,pgm_freq_header, stat_threshold, pgm_freq_threshold):

    nm=raw_table.loc[(raw_table[group_column_header]==nm_label)&
                 ((raw_table[nm_stat_header]<=stat_threshold))&
                 (raw_table[pgm_freq_header]>=pgm_freq_threshold)]
    
    mods=raw_table.loc[(raw_table[group_column_header]!=nm_label)&
                   ((raw_table[mod_stat_header]<=stat_threshold))&
                   (raw_table[pgm_freq_header]>=pgm_freq_threshold)]
    
    filtered_table=pd.concat([mods,nm],axis=0).reset_index(drop=True)

    return filtered_table



def prepare_table(filtered_table,samples,pgm_column_header, integrations_columns, columns_to_save,prot_column_header,
                  group_column_header):
    
    filtered_table_nodup=filtered_table.drop_duplicates(pgm_column_header)
    headers=pd.MultiIndex.from_product([integrations_columns,samples])
    new_headers=columns_to_save+headers.to_list()
    def_table=filtered_table_nodup.loc[:,new_headers].reset_index(drop=True)
    prot_list=def_table[prot_column_header]
    def_table[('pgm','BN')]=[f'A{i}_{prot_list[i]}' if 'NM' not in str(def_table[group_column_header][i]) else 
                             f'A{i}_NM_{prot_list[i]}' for i in range(len(def_table[group_column_header]))]

    return def_table


def norm_files(tabla,integrations_columns,prot_integration_label, nm_integration_label, mod_integration_label,
               group_column_header,nm_label,
               samples, prot_column_header):

    dfs=[]

    for i in integrations_columns:

        if i==prot_integration_label:
            headers=pd.MultiIndex.from_product([[prot_integration_label],samples])
            df=tabla.loc[:,[prot_column_header]+headers.to_list()].drop_duplicates(prot_column_header).set_index(prot_column_header).T
            df_norm=((df-df.mean())/df.std()).reset_index(names=['integration','Sample']).iloc[:,1:]
        
        elif i == nm_integration_label:
            headers=pd.MultiIndex.from_product([[nm_integration_label],samples])
            df=tabla.loc[tabla[group_column_header]==nm_label,headers.to_list()+[('pgm','BN')]].set_index(('pgm','BN')).T
            df_norm=((df-df.mean())/df.std()).reset_index(names=['integration','Sample']).iloc[:,1:]
        
        elif i== mod_integration_label:
            headers_mod=pd.MultiIndex.from_product([[mod_integration_label],samples])
            df=tabla.loc[tabla[group_column_header]!=nm_label,headers_mod.to_list()+[('pgm','BN')]].set_index(('pgm','BN')).T
            df_norm=((df-df.mean())/df.std()).reset_index(names=['integration','Sample']).iloc[:,1:]


        df.columns.name=None
        df_norm.set_index('Sample',inplace=True)

        dfs.append(df_norm)
    
    df_final=pd.concat(dfs,axis=1).reset_index()

    return df_final



def binary_regression(binary_table, formulas,variables, correct_param):

    coefs=[]
    pvalues=[]
    ci=[]
    nobs=[]

    for i in range(len(formulas)):

        try:

            model=sm.formula.glm(formulas[i].replace('-','_'), family=sm.families.Binomial(),data=binary_table).fit()
            ods_r=zip(np.exp(model.params))
            coeficiente=[i[0] for i in ods_r]

            pvalue=zip(model.pvalues)
            pvalue_params=[i[0] for i in pvalue]

            # exp_ci=np.exp(model.conf_int())
            cis=np.exp(model.conf_int()[0]).to_list()+np.exp(model.conf_int()[1]).to_list()

            coefs.append(coeficiente)
            pvalues.append(pvalue_params)
            ci.append(cis)
            nobs.append(model.nobs)

        except ValueError:
            coefs.append(['NO']*3)
            pvalues.append(['NO']*3)
            ci.append(['NO']*6)
            nobs.append(0)
            pass

    coefs_df=pd.DataFrame(coefs)
    pvalues_df=pd.DataFrame(pvalues)
    ci_df=pd.DataFrame(ci)
    # nobs_df=pd.DataFrame(nobs)

    df_results=pd.concat([coefs_df,pvalues_df,ci_df],axis=1)

    params=['OR','pvalues','CI_2.5','CI_97.5']

    if correct_param ==1:
        columns=['Intercept','Protein','Mod']
    else:
        columns=['Intercept','Mod']
    # columns=['Intercept','Protein']
    header=pd.MultiIndex.from_product([params,columns])

    df_results.columns=header
    # df_final=df_results.set_axis(mods+nms,axis='index')
    df_final=df_results.set_axis(variables,axis='index')
    df_final['Nº Obs']=nobs

    from operator import itemgetter
    df_final=df_final.reindex(sorted(df_final.columns, key=itemgetter(1)),axis=1)

    return df_final


def merge_report(report,data_table,columns_to_save):

    meta_data=data_table.set_index(('pgm','BN')).loc[:,columns_to_save]
    final_report=pd.merge(meta_data,report,left_index=True,right_index=True)

    return final_report


def plot_joiner(final_report, path_plots, path_plot_filtered, prot_column_header,
                filter_label,contrast, outfolder):


    if path_plots:


        ptmMapPath_nofilt = path_plots
        ptmMapPath_filt = path_plot_filtered
        plotted_q = [os.path.splitext(j)[0] for j in os.listdir(ptmMapPath_nofilt)]
        final_report[('PTMMap','NoFilt')]=[f"=HYPERLINK(\"{os.path.join(ptmMapPath_nofilt, j)}.html\", \"{j}\")" if j in plotted_q and
                                           contrast in os.path.dirname(path_plots) and contrast in os.path.dirname(path_plot_filtered)
                                           else '' for j in final_report[prot_column_header]]
        
        final_report[('PTMMap','Filt')]=[f"=HYPERLINK(\"{os.path.join(ptmMapPath_filt, j)}.html\", \"{j}\")" if j in plotted_q and 
                                         contrast in os.path.dirname(path_plots) and contrast in os.path.dirname(path_plot_filtered) 
                                         else '' for j in final_report[prot_column_header]]
        
        if contrast not in os.path.dirname(path_plots) or contrast not in os.path.dirname(path_plot_filtered):

            logging.warning('You have provide PTMMaps that do not belong to the contrast you are analyzing')


    else:

        ptmMapPath_nofilt = os.path.join(os.path.dirname(outfolder), f'PTMMaps/{contrast}/plots')
        ptmMapPath_filt = os.path.join(os.path.dirname(outfolder), f'PTMMaps/{contrast}/plots_{filter_label}')
        plotted_q = [os.path.splitext(j)[0] for j in os.listdir(ptmMapPath_nofilt)]
        ptmMapPathExcel_nofilt = f'../PTMMaps/{contrast}/plots'
        ptmMapPathExcel_filt = f'../PTMMaps/{contrast}/plots_{filter_label}'

        if os.path.exists(ptmMapPath_filt) and os.path.exists(ptmMapPath_nofilt):

            final_report[('PTMMap','NoFilt')]=[f"=HYPERLINK(\"{os.path.join(ptmMapPathExcel_nofilt, j)}.html\", \"{j}\")" if j in plotted_q else j for j in final_report[prot_column_header]]
            final_report[('PTMMap','Filt')]=[f"=HYPERLINK(\"{os.path.join(ptmMapPathExcel_filt, j)}.html\", \"{j}\")" if j in plotted_q else ''  for j in final_report[prot_column_header]]

        else:
            final_report[('PTMMap','NoFilt')]=''
            final_report[('PTMMap','Filt')]=''
            logging.warning('Your PTMMaps paths are wrong')

    return final_report



def report_format (report_plots, path_report, correct_param):

    header=list(zip(*report_plots.columns.to_list()))
    report_plots.columns = np.arange(0, report_plots.shape[1])
    report_plots=pd.concat([pd.DataFrame(header),report_plots])
    report_plots.to_excel(path_report,index=False, header=False)

    toFormat = [n+1 for n,i in enumerate(report_plots.iloc[:, -1]) if 'HYPERLINK' in i]
    toFormat2 = [n+1 for n,i in enumerate(report_plots.iloc[:, -2]) if 'HYPERLINK' in i]


    book = openpyxl.load_workbook(path_report)

    sheet = book['Sheet1']

    if correct_param==1:
        columns=['X','Y']
    else:
        columns=['T','U']

    for i in toFormat:
        sheet[f'{columns[0]}{i}'].font = openpyxl.styles.Font(color='0000FF', underline='single')

    for i in toFormat2:
        sheet[f'{columns[1]}{i}'].font = openpyxl.styles.Font(color='0000FF', underline='single')

    book.save(path_report)


def write_report(report,data_table,columns_to_save, path_report):

    meta_data=data_table.set_index(('pgm','BN')).loc[:,columns_to_save]
    final_report=pd.merge(meta_data,report,left_index=True,right_index=True)
    final_report.to_excel(path_report)

    return final_report



def main(path,path_exp,group_column_header, nm_label,nm_stat_header_contrast,
         mod_stat_header_contrast,pgm_freq_header, stat_threshold, pgm_freq_threshold,
         integrations_columns, columns_to_save, pgm_column_header,prot_column_header,
         prot_integration_label, nm_integration_label, mod_integration_label,path_report,
         binary_exp_table_label,groups, path_plots, path_plot_filtered, filter_label,outfolder,
         correct_param):

    logging.info('Reading files')

    limma_table=read_file(path,True)

    exp_table=read_file(path_exp,False)

    logging.info('Applying filters'+'\n')

    logging.info('pvalue threshold: '+str(stat_threshold))
    logging.info('pgmFrequency threshold: '+ str(pgm_freq_threshold)+'\n')

    samples=exp_table[groups]['Sample'].dropna().to_list()


    filtered_table=filter_table(limma_table,group_column_header, nm_label,nm_stat_header_contrast,
                                mod_stat_header_contrast, pgm_freq_header, stat_threshold, pgm_freq_threshold)
    

    logging.info('Nº pgms before filtering: '+ str(len(limma_table)))
    logging.info('Nº pgms after filtering: '+ str(len(filtered_table))+'\n')


    prepared_table=prepare_table(filtered_table,samples,pgm_column_header,integrations_columns,columns_to_save,
                            prot_column_header,group_column_header)
    


    logging.info('Standarizing data')
    
    
    norm_table=norm_files(prepared_table,integrations_columns,prot_integration_label, nm_integration_label, mod_integration_label,
               group_column_header,nm_label,
               samples, prot_column_header)
    

    binary_table=pd.merge(exp_table[groups],norm_table,left_on='Sample',right_on='Sample')

    binary_table.columns=[i.replace('-','_') for i in binary_table.columns]

    if correct_param ==1:

        formulas=[f'{binary_exp_table_label}~{i.split("_")[-1]}+{i}' for i in prepared_table[('pgm','BN')]]

    else:

        logging.info('You are not correcting by protein\'s Zq')
        
        formulas=[f'{binary_exp_table_label}~{i}' for i in prepared_table[('pgm','BN')]]

    variables=prepared_table[('pgm','BN')].to_list()


    logging.info('Applying Binary logistic regression')

    
    report=binary_regression(binary_table,formulas,variables, correct_param)

    logging.info('Writting file')


    merged_report=merge_report(report,prepared_table,columns_to_save)


    report_plots=plot_joiner(merged_report, path_plots,path_plot_filtered,prot_column_header,
                             filter_label,groups, outfolder)
    
    
    report_format(report_plots, path_report,correct_param)


    return report_plots






if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description='PTMs Odds ratio Calculator',
    epilog='''
    Example:
        python PTMs_OR.py -i infile -o outfile -e experiment table -c config

    ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "OR_calculator.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Zqs file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-e', '--experiment_table', required=True, help='Path to metadata.tsv file')
    parser.add_argument('-o', '--outfolder', required=True, help='Name of the output folder')

    args = parser.parse_args()

    path=args.infile
    path_exp=args.experiment_table
    path_report=args.outfolder
 


    # start main function



    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    params = dict(config.items('Binary Regression'))


    """Reading parameters"""

    group_column_header=tuple(eval(params['group_column_header']))
    nm_label=params['nm_label']
    nm_stat_header=eval(params['nm_statistic_header'])
    mod_stat_header=eval(params['mod_statistic_header'])
    pgm_freq_header=tuple(eval(params['pgm_freq_header']))
    stat_threshold=float(params['stats_threshold'])
    pgm_freq_threshold=float(params['pgm_freq_threshold'])
    pgm_column_header=tuple(eval(params['pgm_column_header']))
    prot_column_header=tuple(eval(params['protein_column_header']))

    prot_integration_label=params['prot_integration_label']
    nm_integration_label=params['nm_integration_label']
    mod_integration_label=params['mod_integration_label']
    binary_exp_table_label=params['binary_exp_table_label']
    correct_param=int(params['correct_by_protein'])

    integrations_columns= [prot_integration_label,nm_integration_label,mod_integration_label]
    columns_to_save=[pgm_column_header, group_column_header, pgm_freq_header, 
                     prot_column_header]+[tuple(eval(params[i])) for i in list(params.keys()) if '.save' in i]

    groups= eval(params['groups'])


    path_plots=params['plot_total']
    path_plot_filtered=params['plot_filtered']
    filter_label= params['filter_label']

    outfolder=os.path.join(path_report,'OR_results')

    check_folder(outfolder)
    log_file = os.path.join(outfolder,f'ORs_log.txt')

    for i in groups:
        
        nm_stat_header_contrast=nm_stat_header.copy()
        nm_stat_header_contrast[0]=nm_stat_header_contrast[0]+f'_{i}'
        nm_stat_header_contrast=tuple(nm_stat_header_contrast)

        mod_stat_header_contrast=mod_stat_header.copy()
        mod_stat_header_contrast[0]=mod_stat_header_contrast[0]+f'_{i}'
        mod_stat_header_contrast=tuple(mod_stat_header_contrast)


        path_report_contrast=os.path.join(outfolder,f'ORs_{i}.xlsx')

        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                        logging.StreamHandler()])
        
        print('\n')
        
        logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv]))+'\n')

        logging.info(f'Performing PTM OR calculation on {i} contrast'+'\n')

        main(path,path_exp,group_column_header, nm_label,nm_stat_header_contrast,
         mod_stat_header_contrast,pgm_freq_header, stat_threshold, pgm_freq_threshold,
         integrations_columns, columns_to_save, pgm_column_header,prot_column_header,
         prot_integration_label, nm_integration_label, mod_integration_label,path_report_contrast,
         binary_exp_table_label,i, path_plots, path_plot_filtered, filter_label, outfolder,
         correct_param)
        

        logging.info('End script'+'\n')
