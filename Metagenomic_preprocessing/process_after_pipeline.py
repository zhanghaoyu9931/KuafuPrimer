import os
import json
import pandas as pd
import argparse

def read_abun_file(directory, id):
    """
    读取每个ID对应的 abun 文件
    """
    abun_file_path = os.path.join(directory, f'clean_reads/{id}/{id}_step4_abun.txt')
    try:
        df = pd.read_csv(abun_file_path, sep='\t', usecols=['name', 'fraction_total_reads'])
        df.rename(columns={'fraction_total_reads': f'{id}_abun'}, inplace=True)
        
        # 合并 'uncultured' 行，对应的 'fraction_total_reads_{id}' 数值相加
        uncultured_rows = df[df['name'] == 'uncultured']
        if not uncultured_rows.empty:
            sum_uncultured = {'name': 'uncultured', f'{id}_abun': uncultured_rows[f'{id}_abun'].sum()}
            sum_uncultured = pd.DataFrame([sum_uncultured])
            df = df[df['name'] != 'uncultured']
            df = pd.concat([df, sum_uncultured], axis=0).reset_index(drop=True)
    
    except FileNotFoundError:
        print(f"Error: Abun file not found for ID '{id}' at '{abun_file_path}'.")
        df = pd.DataFrame({'name': ['nop'], f'{id}_abun': [1.0]})
    
    ### for ncbi kraken2 ###
    abun_file_path_ncbi = os.path.join(directory, f'clean_reads/{id}/{id}_step4_ncbi_abun.txt')
    try:
        df_ncbi = pd.read_csv(abun_file_path_ncbi, sep='\t', usecols=['name', 'fraction_total_reads'])
        df_ncbi.rename(columns={'fraction_total_reads': f'{id}_abun'}, inplace=True)
        
        # 合并 'uncultured' 行，对应的 'fraction_total_reads_{id}' 数值相加
        uncultured_rows = df_ncbi[df_ncbi['name'] == 'uncultured']
        if not uncultured_rows.empty:
            sum_uncultured = {'name': 'uncultured', f'{id}_abun': uncultured_rows[f'{id}_abun'].sum()}
            sum_uncultured = pd.DataFrame([sum_uncultured])
            df_ncbi = df_ncbi[df_ncbi['name'] != 'uncultured']
            df_ncbi = pd.concat([df_ncbi, sum_uncultured], axis=0).reset_index(drop=True)
    
    except FileNotFoundError:
        print(f"Error: Abun_ncbi file not found for ID '{id}' at '{abun_file_path_ncbi}'.")
        df_ncbi = pd.DataFrame({'name': ['nop'], f'{id}_abun': [1.0]})
    
    abun_file_path_ncbi_S = os.path.join(directory, f'clean_reads/{id}/{id}_step4_ncbi_abun_S.txt')
    try:
        df_ncbi_S = pd.read_csv(abun_file_path_ncbi_S, sep='\t', usecols=['name', 'fraction_total_reads'])
        df_ncbi_S.rename(columns={'fraction_total_reads': f'{id}_abun'}, inplace=True)
        
        # 合并 'uncultured' 行，对应的 'fraction_total_reads_{id}' 数值相加
        uncultured_rows = df_ncbi_S[df_ncbi_S['name'] == 'uncultured']
        if not uncultured_rows.empty:
            sum_uncultured = {'name': 'uncultured', f'{id}_abun': uncultured_rows[f'{id}_abun'].sum()}
            sum_uncultured = pd.DataFrame([sum_uncultured])
            df_ncbi_S = df_ncbi_S[df_ncbi_S['name'] != 'uncultured']
            df_ncbi_S = pd.concat([df_ncbi_S, sum_uncultured], axis=0).reset_index(drop=True)
    
    except FileNotFoundError:
        print(f"Error: Abun_ncbi file not found for ID '{id}' at '{abun_file_path_ncbi}'.")
        df_ncbi_S = pd.DataFrame({'name': ['nop'], f'{id}_abun': [1.0]})
    
    # 读取元信息
    meta_info_path = os.path.join(directory, f'tmpfiles/{id}_step1.json')
    try:
        with open(meta_info_path, 'r') as meta_file:
            meta_info = json.load(meta_file)
            sequencing_info = meta_info['summary']['sequencing']
            total_reads_before_filtering = meta_info['summary']['before_filtering']['total_reads']
            total_reads_after_filtering = meta_info['summary']['after_filtering']['total_reads']
    except FileNotFoundError:
        print(f"Error: Meta info file not found for ID '{id}' at '{meta_info_path}'.")
        return df, df_ncbi, df_ncbi_S, {'id': id,}

    ## 废弃silva的：读取 reads_after_bowtie2
    classified_reads_num = 0.0
    ## 废弃silva的：读取 reads of bacteria
    rp_file_path_step4 = os.path.join(directory, f'clean_reads/{id}/{id}_step4_rp.txt')
    try:
        rp_df = pd.read_csv(rp_file_path_step4, sep='\t', header=None)
        classified_bacterial_reads_num = rp_df.iloc[0, 1]
    except FileNotFoundError:
        print(f"Error: step4 RP file not found for ID '{id}' at '{rp_file_path_step4}'.")
        classified_bacterial_reads_num = 0.0
    
    # 读取 reads_after_bowtie2
    rp_file_path_ncbi = os.path.join(directory, f'clean_reads/{id}/{id}_step3_ncbi_rp.txt')
    try:
        rp_ncbi_df = pd.read_csv(rp_file_path_ncbi, sep='\t', header=None)
        classified_reads_num_ncbi = rp_ncbi_df.iloc[1, 1]
        unclassified_reads_num_ncbi = rp_ncbi_df.iloc[0, 1]
        reads_after_bowtie2 = classified_reads_num_ncbi + unclassified_reads_num_ncbi
    except:
        print(f"Error: step3_ncbi RP file not found for ID '{id}' at '{rp_file_path_ncbi}'.")
        reads_after_bowtie2, classified_reads_num_ncbi = 0.0, 0.0
    
    # 读取 reads of bacteria
    rp_file_path_ncbi_step4 = os.path.join(directory, f'clean_reads/{id}/{id}_step4_ncbi_rp.txt')
    try:
        rp_ncbi_df = pd.read_csv(rp_file_path_ncbi_step4, sep='\t', header=None)
        classified_bacterial_reads_num_ncbi = rp_ncbi_df.iloc[0, 1]
    except:
        print(f"Error: step4_ncbi RP file not found for ID '{id}' at '{rp_file_path_ncbi_step4}'.")
        classified_bacterial_reads_num_ncbi = 0.0
    
    
    return df, df_ncbi, df_ncbi_S, {
        'id': id,
        'sequencing_info': sequencing_info,
        'total_reads_before_filtering': total_reads_before_filtering,
        'total_reads_after_filtering': total_reads_after_filtering,
        'reads_after_bowtie2': reads_after_bowtie2,
        'reads_after_bowtie2_fraction': reads_after_bowtie2 / total_reads_before_filtering if total_reads_before_filtering > 0 else 0.0,    
        'reads_classified': classified_reads_num,
        'reads_classified_bacterial': classified_bacterial_reads_num,
        'reads_classified_ncbi': classified_reads_num_ncbi,
        'reads_classified_bacterial_ncbi': classified_bacterial_reads_num_ncbi,
    }

def merge_dataframes(dataframes):
    """
    合并多个数据框，以'name'列为键，取outer并以0.0填充nan
    """
        
    merged_df = dataframes[0]
    print(f'Parse {len(dataframes)} dataframes')
    for df in dataframes[1:]:
        print(f'Parse {df.columns}, {df.shape}')
        merged_df = pd.merge(merged_df, df, on='name', how='outer',).fillna(0.0)
    return merged_df

def main():
    """
    主要的处理函数：
    1. 读取并合并每个 ID 对应的 abun 文件，打印合并后的数据框
    2. 获得每个 ID 对应的中间reads数统计文件，打印metainfo
    """
    parser = argparse.ArgumentParser(description="Data processing script.")
    parser.add_argument("--metagenomic_dir", type=str, help="Path to the metagenomic directory.")
    parser.add_argument("--id_list", type=str, help="Path to the ID list file.")

    args = parser.parse_args()

    # 读取 ID 列表
    try:
        with open(args.id_list, 'r') as id_file:
            id_list = [line.strip() for line in id_file]
    except FileNotFoundError:
        print(f"Error: ID list file '{args.id_list}' not found.")
        exit(1)

    # 读取并合并每个 ID 对应的 abun 文件
    dataframes = []
    dataframes_ncbi = []
    dataframes_ncbi_S = []
    meta_info_dict = []
    for id in id_list:
        df, df_ncbi, df_ncbi_S, meta_info = read_abun_file(args.metagenomic_dir, id)
        dataframes.append(df)
        dataframes_ncbi.append(df_ncbi)
        dataframes_ncbi_S.append(df_ncbi_S)
        meta_info_dict.append(meta_info)

    merged_df = merge_dataframes(dataframes)
    merged_ncbi_df = merge_dataframes(dataframes_ncbi)
    merged_ncbi_df_S = merge_dataframes(dataframes_ncbi_S)
    meta_info_dict = pd.DataFrame(meta_info_dict)
    
    # 打印合并后的数据框
    print(f"Merged DataFrame: {merged_df.shape}")
    
    # 输出结果
    os.makedirs(f'{args.metagenomic_dir}/MetaAbun', exist_ok=True)
    out_abunTab = os.path.join(f'{args.metagenomic_dir}/MetaAbun', 'merged_abun_table.tsv')
    merged_df.rename(columns={'name': 'Genus'}, inplace=True)
    merged_ncbi_df.rename(columns={'name': 'Genus'}, inplace=True)
    merged_ncbi_df_S.rename(columns={'name': 'Species'}, inplace=True)
    
    # 合并重复的genus名称
    merged_df = merged_df.groupby(by='Genus').sum()
    merged_df = merged_df.loc[merged_df.mean(axis=1).sort_values(ascending=False).index, :] # 按照丰度排序
    merged_ncbi_df = merged_ncbi_df.groupby(by='Genus').sum()
    merged_ncbi_df = merged_ncbi_df.loc[merged_ncbi_df.mean(axis=1).sort_values(ascending=False).index, :] # 按照丰度排序
    merged_ncbi_df_S = merged_ncbi_df_S.groupby(by='Species').sum()
    merged_ncbi_df_S = merged_ncbi_df_S.loc[merged_ncbi_df_S.mean(axis=1).sort_values(ascending=False).index, :] # 按照丰度排序
    
    merged_df.to_csv(out_abunTab, sep='\t', index=True)
    merged_ncbi_df.to_csv(out_abunTab.replace('merged_abun_table', 'merged_abun_table_ncbi'), sep='\t', index=True)
    merged_ncbi_df_S.to_csv(out_abunTab.replace('merged_abun_table', 'merged_abun_table_ncbi_S'), sep='\t', index=True)
    meta_info_dict.to_csv(out_abunTab.replace('merged_abun_table', 'metainfo'), sep='\t', index=False)

if __name__ == "__main__":
    main()