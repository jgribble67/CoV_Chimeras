import pandas as pd
import os
import fnmatch
import numpy as np

dir = "/home/denison-thelio/Current_projects/RNAseq/MERS-SARS2_recombination_control/Chimeric_Junctions/"
report = pd.DataFrame(columns=['sample',
                               'MERStoSARS2_junctions',
                               'MERStoSARS2_depth',
                               'SARS2toMERS_junctions',
                               'SARS2toMERS_depth',
                               'total_depth'])
sample_list = [line.rstrip('\n') for line in open("/home/denison-thelio/Current_projects/RNAseq/MERS-SARS2_recombination_control/Chimeric_Junctions/control_samples.txt")]
report['sample'] = sample_list
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, "*_MERS_to_SARS2_junctions.txt"):
        sample_name = str(file.split("_")[0] + "_" + file.split("_")[1])
        df_MERStoSARS2 = pd.read_csv(dir + file, skiprows=1, sep="\t")
        df_MERStoSARS2 = df_MERStoSARS2.T
        df_MERStoSARS2 = df_MERStoSARS2.reset_index()
        df_MERStoSARS2[['MERS_start', 'seq', 'SARS2_stop', 'unk', 'Depth']] = df_MERStoSARS2['index'].str.split("_", expand=True)
        df_MERStoSARS2 = df_MERStoSARS2[['MERS_start', 'SARS2_stop', 'Depth']]
        df_MERStoSARS2 = df_MERStoSARS2[df_MERStoSARS2['MERS_start'].apply(lambda x: x.isnumeric())]
        df_MERStoSARS2 = df_MERStoSARS2[df_MERStoSARS2['SARS2_stop'].apply(lambda x: x.isnumeric())]
        df_MERStoSARS2['MERS_start'] = pd.to_numeric(df_MERStoSARS2['MERS_start'])
        df_MERStoSARS2['SARS2_stop'] = pd.to_numeric(df_MERStoSARS2['SARS2_stop'])
        df_MERStoSARS2['Depth'] = pd.to_numeric(df_MERStoSARS2['Depth'])
        df_MERStoSARS2 = df_MERStoSARS2.loc[(df_MERStoSARS2['MERS_start'] < 30110) & (df_MERStoSARS2['SARS2_stop'] < 29870)]
        MERStoSARS2_junctions = df_MERStoSARS2.shape[0]
        depth = sum(df_MERStoSARS2['Depth'])
        report.loc[report['sample'] == sample_name, ['MERStoSARS2_junctions']] = MERStoSARS2_junctions
        report.loc[report['sample'] == sample_name, ['MERStoSARS2_depth']] = depth
        df_MERStoSARS2.to_csv(dir + sample_name + "_MERStoSARS2_junctions.txt", sep="\t", index=False)
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, "*_SARS2_to_MERS_junctions.txt"):
        sample_name = str(file.split("_")[0] + "_" + file.split("_")[1])
        df = pd.read_csv(dir + file, skiprows=1, sep="\t")
        df = df.T
        df = df.reset_index()
        df[['SARS2_start', 'seq', 'MERS_stop', 'unk', 'Depth']] = df['index'].str.split("_", expand=True)
        df = df[['SARS2_start', 'MERS_stop', 'Depth']]
        df = df[df['SARS2_start'].apply(lambda x: x.isnumeric())]
        df = df[df['MERS_stop'].apply(lambda x: x.isnumeric())]
        df['SARS2_start'] = pd.to_numeric(df['SARS2_start'])
        df['MERS_stop'] = pd.to_numeric(df['MERS_stop'])
        df['Depth'] = pd.to_numeric(df['Depth'])
        df = df.loc[(df['SARS2_start'] < 29870) & (df['MERS_stop'] < 30110)]
        SARS2toMERS_junctions = df.shape[0]
        report.loc[report['sample'] == sample_name, ['SARS2toMERS_junctions']] = SARS2toMERS_junctions
        report.loc[report['sample'] == sample_name, ['SARS2toMERS_depth']] = depth
        df.to_csv(dir + sample_name + "_SARS2toMERS_junctions.txt", sep="\t", index=False)
# for file in os.listdir(dir):
#     if fnmatch.fnmatch(file, "*_MERS_junctions.txt"):
#         sample_name = str(file.split("_")[0] + "_" + file.split("_")[1])
#         df = pd.read_csv(dir + file, skiprows=1, sep="\t")
#         df = df.T
#         df = df.reset_index()
#         df[['start', 'seq', 'stop', 'unk', 'Depth']] = df['index'].str.split("_", expand=True)
#         df = df[['start', 'stop', 'Depth']]
#         df = df[df['start'].apply(lambda x: x.isnumeric())]
#         df = df[df['stop'].apply(lambda x: x.isnumeric())]
#         df['start'] = pd.to_numeric(df['start'])
#         df['stop'] = pd.to_numeric(df['stop'])
#         df['Depth'] = pd.to_numeric(df['Depth'])
#         df = df.loc[(df['start'] < 30110) & (df['stop'] < 30110)]
#         MERS_junctions = df.shape[0]
#         report.loc[report['sample'] == sample_name, ['MERS_junctions']] = MERS_junctions
#         report.loc[report['sample'] == sample_name, ['MERS_depth']] = depth
#         df.to_csv(dir + sample_name + "_MERS_junctions_parsed.txt", sep="\t", index=False)
# for file in os.listdir(dir):
#     if fnmatch.fnmatch(file, "*_SARS2_junctions.txt"):
#         sample_name = str(file.split("_")[0] + "_" + file.split("_")[1])
#         df = pd.read_csv(dir + file, skiprows=1, sep="\t")
#         df = df.T
#         df = df.reset_index()
#         df[['start', 'seq', 'stop', 'unk', 'Depth']] = df['index'].str.split("_", expand=True)
#         df = df[['start', 'stop', 'Depth']]
#         df = df[df['start'].apply(lambda x: x.isnumeric())]
#         df = df[df['stop'].apply(lambda x: x.isnumeric())]
#         df['start'] = pd.to_numeric(df['start'])
#         df['stop'] = pd.to_numeric(df['stop'])
#         df['Depth'] = pd.to_numeric(df['Depth'])
#         df = df.loc[(df['start'] < 29870) & (df['stop'] < 29870)]
#         SARS2_junctions = df.shape[0]
#         report.loc[report['sample'] == sample_name, ['SARS2_junctions']] = SARS2_junctions
#         report.loc[report['sample'] == sample_name, ['SARS2_depth']] = depth
#         df.to_csv(dir + sample_name + "_SARS2_junctions_parsed.txt", sep="\t", index=False)
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, "*_virema_coverage.txt"):
        sample_name = str(file.split("_")[0] + "_" + file.split("_")[1])
        virus = str(file.split("_")[0])
        df = pd.read_csv(dir + file, sep="\t", header=0, index_col=False, names=['Genome', 'Position', 'Coverage'])
        if virus == "MERS":
            MERS_depth = df.loc[df['Genome'] == "JX869059.2", 'Coverage'].sum()
            report.loc[report['sample'] == sample_name, ['total_depth']] = MERS_depth
        if virus == "SARSCoV2":
            SARS2_depth = df.loc[df['Genome'] == "MT020881.1", 'Coverage'].sum()
            report.loc[report['sample'] == sample_name, ['total_depth']] = SARS2_depth
report.to_csv(dir + "MERS-SARS2_chimeric_junctions.txt", sep="\t", index=False)