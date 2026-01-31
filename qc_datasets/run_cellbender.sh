


cd /mnt/sdb/scz_meta_analysis_processed/fernando_brennand/cellbender/

pixi run --manifest-path /home/deepak/pixi_envs/cellbender/pixi.toml cellbender remove-background --input /mnt/sdb/scz_meta_analysis_processed/fernando_brennand/cellbender/fernando.h5 --output /mnt/sdb/scz_meta_analysis_processed/fernando_brennand/cellbender/cellbendered_fernando.h5 --cuda

cd /mnt/sdb/scz_meta_analysis_processed/sawada_kato/cellbender/

pixi run --manifest-path /home/deepak/pixi_envs/cellbender/pixi.toml cellbender remove-background --input /mnt/sdb/scz_meta_analysis_processed/sawada_kato/cellbender/sawada.h5 --output /mnt/sdb/scz_meta_analysis_processed/sawada_kato/cellbender/cellbendered_sawada.h5 --cuda

