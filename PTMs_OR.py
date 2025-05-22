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
warnings.filterwarnings("ignore")



""" Functions """


def read_file(path, raw_table):

    if raw_table==True:
        limma_table=pd.read_csv(path,sep='\t', header=[0,1])
        return limma_table
    else:
        exp_table=pd.read_csv(path,sep='\t')
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



def prepare_table(filtered_table,samples,pgm_column_header, integrations_columns, 
                  columns_to_save,prot_column_header, group_column_header):
    
    filtered_table_nodup=filtered_table.drop_duplicates(pgm_column_header)
    headers=pd.MultiIndex.from_product([integrations_columns,samples])
    new_headers=columns_to_save+headers.to_list()
    def_table=filtered_table_nodup.loc[:,new_headers].reset_index(drop=True)
    prot_list=def_table[prot_column_header]
    def_table[('pgm','BN')]=[f'A{i}_{prot_list[i]}' if 'NM' not in str(def_table[group_column_header][i]) else f'A{i}_NM_{prot_list[i]}' for i in range(len(def_table[group_column_header]))]

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



def binary_regression(binary_table, formulas,variables):

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
    columns=['Intercept','Protein','Mod']
    # columns=['Intercept','Protein']
    header=pd.MultiIndex.from_product([params,columns])

    df_results.columns=header
    # df_final=df_results.set_axis(mods+nms,axis='index')
    df_final=df_results.set_axis(variables,axis='index')
    df_final['Nº Obs']=nobs

    from operator import itemgetter
    df_final=df_final.reindex(sorted(df_final.columns, key=itemgetter(1)),axis=1)

    return df_final



def write_report(report,data_table,columns_to_save, path_report):

    meta_data=data_table.set_index(('pgm','BN')).loc[:,columns_to_save]
    final_report=pd.merge(meta_data,report,left_index=True,right_index=True)
    final_report.to_excel(path_report)

    return final_report



def main(path,path_exp,group_column_header, nm_label,nm_stat_header,
         mod_stat_header,pgm_freq_header, stat_threshold, pgm_freq_threshold,
         integrations_columns, columns_to_save, pgm_column_header,prot_column_header,
         prot_integration_label, nm_integration_label, mod_integration_label,path_report,
         binary_exp_table_label):

    print('\n')
    logging.info('Reading files')

    limma_table=read_file(path,True)

    exp_table=read_file(path_exp,False)
    samples=exp_table['Sample'].to_list()

    logging.info('Applying filters')

    print('\n'+'pvalue threshold:',stat_threshold)
    print('pgmFrequency threshold:', pgm_freq_threshold,'\n')

    filtered_table=filter_table(limma_table,group_column_header, nm_label,nm_stat_header,
                                mod_stat_header, pgm_freq_header, stat_threshold, pgm_freq_threshold)
    
    print('Nº pgms before filtering:', len(limma_table))
    print('Nº pgms after filtering:', len(filtered_table),'\n')


    

    prepared_table=prepare_table(filtered_table,samples,pgm_column_header,integrations_columns,columns_to_save,
                            prot_column_header,group_column_header)
    formulas=[f'{binary_exp_table_label}~{i.split("_")[-1]}+{i}' for i in prepared_table[('pgm','BN')]]
    variables=prepared_table[('pgm','BN')].to_list()
    
    logging.info('Standarizing data')


    norm_table=norm_files(prepared_table,integrations_columns,prot_integration_label, nm_integration_label, mod_integration_label,
               group_column_header,nm_label,
               samples, prot_column_header)
    

    binary_table=pd.merge(exp_table,norm_table,left_on='Sample',right_on='Sample')
    binary_table.columns=[i.replace('-','_') for i in binary_table.columns]

    # binary_table.to_excel(r'S:\U_Proteomica\LABS\LAB_ARR\LaCaixa\tejidos-secretomas\New_Comet_500\final_iSanxot\Intima\Q1\reports\WF\3_FDRoptimizer\file_to_check.xlsx')

    logging.info('Applying Binary loggistic regression')
    
    report=binary_regression(binary_table,formulas,variables)

    logging.info('Writting file')

    write_report(report,prepared_table,columns_to_save,path_report)


    return report


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
    parser.add_argument('-o', '--outfile', required=True, help='Name of the output file')

    args = parser.parse_args()

    path=args.infile
    path_exp=args.experiment_table
    path_report=args.outfile
    log_file = path_report[:-4] + '_binary_regression_log.txt'

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        handlers=[logging.FileHandler(log_file),
                                    logging.StreamHandler()])


    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))


    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    params = dict(config.items('Binary Regression'))


    """Reading parameters"""

    group_column_header=tuple(eval(params['group_column_header.save']))
    nm_label=params['nm_label']
    nm_stat_header=tuple(eval(params['nm_statistic_header']))
    mod_stat_header=tuple(eval(params['mod_statistic_header']))
    pgm_freq_header=tuple(eval(params['pgm_freq_header.save']))
    stat_threshold=float(params['stats_threshold'])
    pgm_freq_threshold=float(params['pgm_freq_threshold'])
    pgm_column_header=tuple(eval(params['pgm_column_header.save']))
    prot_column_header=tuple(eval(params['protein_column_header.save']))

    prot_integration_label=params['prot_integration_label']
    nm_integration_label=params['nm_integration_label']
    mod_integration_label=params['mod_integration_label']
    binary_exp_table_label=params['binary_exp_table_label']

    integrations_columns= [prot_integration_label,nm_integration_label,mod_integration_label]
    columns_to_save=[tuple(eval(params[i])) for i in list(params.keys()) if '.save' in i]


    main(path,path_exp,group_column_header, nm_label,nm_stat_header,
         mod_stat_header,pgm_freq_header, stat_threshold, pgm_freq_threshold,
         integrations_columns, columns_to_save, pgm_column_header,prot_column_header,
         prot_integration_label, nm_integration_label, mod_integration_label, path_report,
         binary_exp_table_label)

    logging.info('End script')

    # file_handler = logging.FileHandler(os.path.join(os.path.dirname(path_report),
    #                                                 os.path.basename(path).split('.tsv')[0]+'_ORs_calculator.log'))
    
    # logging.basicConfig(level=logging.INFO,
    #                 format='%(asctime)s - %(levelname)s - %(message)s',
    #                 datefmt='%m/%d/%Y %I:%M:%S %p',
    #                 handlers=[file_handler,
    #                             logging.StreamHandler()])