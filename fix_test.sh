sed -i '' 's/let w = huber_weight(norm_res, tuning.huber_k);/let w = huber_weight(norm_res, tuning.huber_k); println!("norm_res={}, w={}, current_z={}", norm_res, w, current_z[i]);/g' crates/gneiss-rtk/src/engine/updater.rs
cargo test -p gneiss-rtk -- updater --nocapture
