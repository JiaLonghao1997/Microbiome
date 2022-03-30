import pandas as pd
import warnings
import os
Study = ['FR-CRC', 'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC']
root_dir = 'F://Zhaolab2020//BioNetwork//MicrobiomeNetwork//CRC//shaojun//'

if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    if not os.path.exists(root_dir + 'pre_data//species_loso'):
        os.makedirs(root_dir + 'pre_data//species_loso')
    root_dir = root_dir + 'pre_data//species_loso//'
    for study in Study:
        if not os.path.exists(root_dir + '{}'.format(study)):
            os.makedirs(root_dir + '{}'.format(study))

    tag = 'species'
    species_url = 'data//' + tag + '//feat_rel_crc.tsv'
    meta_url = 'data//meta//meta_crc.tsv'
    species = pd.read_csv(species_url, sep='\t', header=0)
    meta = pd.read_csv(meta_url, sep='\t', header=0)

    for study in Study:
        print(study)
        train_Sample_id = []
        train_label = []
        train_meta_row_id = []
        test_Sample_id = []
        test_label = []
        test_meta_row_id = []
        for i in range(len(meta[['Sample_ID', 'Study']])):
            if meta.loc[i, ['Study']].tolist()[0] == study:
                test_Sample_id.append(meta.loc[i, ['Sample_ID']].tolist()[0])
                test_meta_row_id.append(i)
                if meta.loc[i, ['Group']].tolist()[0] == 'CRC':
                    test_label.append(1)
                else:
                    test_label.append(-1)

            else:
                train_Sample_id.append(meta.loc[i, ['Sample_ID']].tolist()[0])
                train_meta_row_id.append(i)
                if meta.loc[i, ['Group']].tolist()[0] == 'CRC':
                    train_label.append(1)
                else:
                    train_label.append(-1)

        train_species_study = species[train_Sample_id].T
        train_species_study['Label'] = train_label
        train_species_study.to_csv(root_dir + '{}'.format(study) + '//train_data.csv')
        test_species_study = species[test_Sample_id].T
        test_species_study['Label'] = test_label
        test_species_study.to_csv(root_dir + '{}'.format(study) + '//test_data.csv')