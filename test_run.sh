# date: 20240831
# 写一个脚本，用于运行EcoPrimer
demo_input=input/demo_input_KuafuPrimer_ge_profile.csv
demo_design_output=output/demo_output_KuafuPrimer
# target_vs='v1v2;v1v3;v3v4;v4v4;v4v5;v6v8;v5v6;v5v7;v7v8'
target_vs='v3v4;v4v5;'

# 设计引物部分
python KuafuPrimer.py --input $demo_input \
    --out_root $demo_design_output \
    --target_vs $target_vs \
    --extend_bp_num 50 \
    --num_every_spe 10 \
    --step_search 5 \
    --NGS_mode Single_end \
    --input_type genera_profiling

# 评估引物部分
k_num=3
demo_PCR_output=output/demo_PCR_output_accurate

python Insilico_eva_primers.py --envi_forEva $demo_design_output'/samples_abundanceTab_clean.csv' \
    --primers_forEva $demo_design_output';' \
    --target_vs $target_vs \
    --K $k_num \
    --output $demo_PCR_output \
    --num_condidate_pris 3 \
    --very_fast

# screen 最优引物部分
demo_PCR_output=$demo_PCR_output'_K'$k_num

python Screen_best_PP.py --pcr_dir $demo_PCR_output \
    --rk_by Genus_accuracy