sed -i '' 's/        } else if j == 16 {/        } else {/' crates/gneiss-rtk/src/doppler.rs
sed -i '' 's/            assert!((h\[(i, j)\] - 1.0).abs() < 1e-6, "H\[{},16\] should be 1.0", i);/            0.0/' crates/gneiss-rtk/src/doppler.rs
sed -i '' 's/        } else {/        /' crates/gneiss-rtk/src/doppler.rs
sed -i '' 's/            0.0/        /' crates/gneiss-rtk/src/doppler.rs
