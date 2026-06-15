cargo build --release
target/release/gneiss-cli process --config data/UrbanNav-Shinjuku/20210519-google-pixel-4-rtk-ins.json --output rtk_ins.pos
target/release/gneiss-cli process --config data/UrbanNav-Shinjuku/20210519-google-pixel-4-spp-ins.json --output spp_ins.pos
python3 scripts/bench.py rtk_ins.pos data/UrbanNav-Shinjuku/ground_truth.csv
python3 scripts/bench.py spp_ins.pos data/UrbanNav-Shinjuku/ground_truth.csv
